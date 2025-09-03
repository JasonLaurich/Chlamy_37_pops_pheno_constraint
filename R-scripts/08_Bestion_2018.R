# Jason R Laurich
# April 3, 2025

# Revisited, checked, and cleaned September 3rd, 2025 (JRL)

# Going to work with the Bestion et al 2018 (https://onlinelibrary.wiley.com/doi/10.1111/ele.12932) data to fit
# Monod curves for phosphorous and TPCs

# Packages & functions ----------------------------------------------------

library(tidyr)
library(cowplot)
library(ggplot2)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)
library(rTPC)
library(nls.multstart)

# Upload & examine data ---------------------------------------------------

df <- read.csv("data-processed/11_Bestion_2018_raw_data.csv") # Raw growth data
head(df)
str(df)

df$Sp.fac <- as.factor(df$SpeciesName)

# So they co-varied phosphorous and temperature
# And they only provide their estimate of µ, so let's move right to fitting models.

# Going to split this data set into 2

# 1. Full phosphorous, variation in temperature. Actually let's go with P30 as well, so we can include Chlamydomonas

# 2. 20 C (consistent with our estimation of P limitation), variation in P

df.t <- df[df$Phosphate_c == 50,]
df.t.p30 <- df[df$Phosphate_c == 30,]

df.t <- droplevels(df.t)
df.t.p30 <-droplevels(df.t.p30)

df.t20 <- df[df$Temperature_c == 20,]
df.t20 <- droplevels(df.t20)

# TPC fitting -------------------------------------------------------------

bestion.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species ID
  DIC = numeric(),          # DIC
  P.conc = numeric(),       # Phosphorous concentration
  T.min.raw = numeric(),    # Minimum T (Jags raw)
  T.max.raw = numeric(),    # Maximum T (Jags raw)
  T.opt.raw = numeric(),    # Optimal T (Jags raw)
  r.max.raw = numeric(),    # Maximum growth rate (Jags raw)
  T.br.raw = numeric(),     # T breadth (Jags raw)
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  P.conc = numeric(),       # Phosphorous
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~3000)
  stringsAsFactors = FALSE            
)

# Let's do larger models for the final things (10 times larger)
ni.fit <- 330000    # iterations / chain
nb.fit <- 30000     # burn in periods for each chain
nt.fit <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 3         # number of chains, total of 3,000 estimates for each model. 

parameters.lactin2 <- c("cf.a", "cf.b", "cf.tmax", "cf.delta_t", "cf.sigma", "r.pred") # repeated here

inits.lactin.meta<- function() { # Pulling initial values centres from the start_vals function in rTPC
  list(
    cf.a = rnorm(1, mean = start.vals.lac[1], sd = 0.05),
    cf.tmax = rnorm(1, mean = start.vals.lac[3], sd = 1),
    cf.delta_t = rnorm(1, mean = start.vals.lac[4], sd = 1),
    cf.b = rnorm(1, mean = start.vals.lac[2], sd = 0.05),
    cf.sigma = runif(1, 0.1, 2)
  )
}

for (i in unique(df.t.p30$Sp.fac)){ # Most species don't have a rich enough dataset at P50, so we're starting with 30
  
  df.i <- df.t.p30 %>%
    filter(Sp.fac == i, !is.na(mu)) %>% 
    arrange(Temperature_c)
  
  max.temp <- max(df.i$Temperature_c, na.rm = TRUE)
  temp.at.max.growth <- df.i$Temperature_c[which.max(df.i$mu)]
  
  if (temp.at.max.growth == max.temp) {
    new_row <- df.i %>%
      filter(Temperature_c == max.temp) %>%
      mutate(
        Temperature_c = max.temp + 5,
        mu = 0
      )
    
    df.i <- bind_rows(df.i, new_row)
  }
  
  trait <- df.i$mu     # format the data for jags
  N.obs <- length(trait)
  temp <- df.i$Temperature_c
  
  Temp.xs <- seq(min(df.i$Temperature_c) - 5, max(df.i$Temperature_c) + 5, 0.1) # Temperature gradient we're interested in - upped the granularity here
  N.Temp.xs <-length(Temp.xs) # We'll reset this internally since the gradient varies substantially
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  start.vals.lac <- get_start_vals(df.i$Temperature_c, df.i$mu, model_name = 'lactin2_1995')
  
  lac_jag <- jags(
    data = jag.data, 
    inits = inits.lactin.meta, 
    parameters.to.save = parameters.lactin2, 
    model.file = "lactin_meta.txt",
    n.thin = nt.fit, 
    n.chains = nc.fit, 
    n.burnin = nb.fit, 
    n.iter = ni.fit, 
    DIC = TRUE, 
    working.directory = getwd()
  ) # ~ 10 min to run?
  
  print(paste("Done", i))
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)[-c(1:6),]   # generate the sequence of r.pred values
  df.jags$temp <- seq(min(df.i$Temperature_c) - 5, max(df.i$Temperature_c) + 5, 0.1)
  
  bestion.summ.df <- rbind(bestion.summ.df, data.frame(                        # Add summary data
    Sp.id = i,                                                                         # Species name 
    DIC = lac_jag$BUGSoutput$DIC,                                                      # DIC
    T.min.raw = df.jags$temp[min(which(df.jags$mean > 0))],                            # Minimum T
    T.max.raw = df.jags$temp[max(which(df.jags$mean > 0))],                            # Maximum T
    T.br.raw = df.jags$temp[max(which(df.jags$mean > (max(df.jags$mean) / 2)))] - 
      df.jags$temp[min(which(df.jags$mean > (max(df.jags$mean) / 2)))],                # T breadth
    T.opt.raw = df.jags$temp[which.max(df.jags$mean)],                                 # Optimal T
    r.max.raw = max(df.jags$mean)                                                      # Maximum growth rate   
  ))
  
  for (j in 1:6){
    fit.df <- rbind(fit.df, data.frame(                         # Model performance data
      Sp.id = i,                                                # Species id   
      Parameter = rownames(lac_jag$BUGSoutput$summary)[j],      # Model parameter (e.g. cf.a, cf.tmax, etc.)
      mean = lac_jag$BUGSoutput$summary[j,1],                   # Posterior mean
      Rhat = lac_jag$BUGSoutput$summary[j,8],                   # Rhat values
      n.eff = lac_jag$BUGSoutput$summary[j,9],                  # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
}

write.csv(bestion.summ.df, "data-processed/11a_Bestion_2018_TPCs.csv") # Save Bestion 2018 TPC summary table
write.csv(fit.df, "data-processed/11b_Bestion_2018_TPCs_fits.csv") # Save model fit summary table

# Monod curves (P) --------------------------------------------------------

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

bestion.summ.P.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species ID
  DIC = numeric(),          # DIC
  Temp = numeric(),         # Temperature
  K.s = numeric(),          # Half-saturation constant
  r.max = numeric(),        # Maximum population growth rate
  R.jag = numeric(),        # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),        # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.P.df <- data.frame(     # Save model fit estimates for examination
  Grad = character(),       # Specifiy abiotic gradient
  Sp.id = numeric(),        # Species id
  Temp = numeric(),         # Temperature
  Parameter = character(),  # Model parameter (e.g. K_s, r_max, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~3000)
  stringsAsFactors = FALSE            
)

for (i in unique(df.t20$Sp.fac)){ # For each species
  
  df.i <- df.t20 %>%
    filter(Sp.fac == i, !is.na(mu)) %>%
    arrange(Phosphate_c)
    
    trait <- df.i$mu      # format the data for jags
    N.obs <- length(trait)
    phos <- df.i$Phosphate_c
    
    S.pred <- seq(min(df.i$Phosphate_c) - 5, max(df.i$Phosphate_c) + 5, 0.1) # Phosphate gradient we're interested in - upped the granularity here
    S.pred <- S.pred[S.pred >= 0]  # Cut off values below 0
    N.S.pred <-length(S.pred) # We'll reset this internally since the gradient varies substantially
  
    jag.data <- list(trait = trait, N.obs = N.obs, S = phos, S.pred = S.pred, N.S.pred = N.S.pred)
    
    monod_jag <- jags( # Run the phosphorous Monod function. 
      data = jag.data,
      inits = inits.monod,
      parameters.to.save = parameters.monod,
      model.file = "monod.txt",
      n.thin = nt.fit,
      n.chains = nc.fit,
      n.burnin = nb.fit,
      n.iter = ni.fit,
      DIC = TRUE,
      working.directory = getwd()
    )
    
    print(paste("Done", i))
    
    df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3, (max(df.i$Phosphate_c) - S.pred[1])/0.5 + 9),]   # generate the sequence of r.pred values
    df.jags$phos <- seq(S.pred[1], max(df.i$Phosphate_c) + 5, 0.1)
    
    bestion.summ.P.df <- rbind(bestion.summ.P.df, data.frame(                                   # Add summary data
      Sp.id = i,                                                                                # Species name 
      DIC = monod_jag$BUGSoutput$DIC,                                                           # DIC
      Temp = 20,                                                                                # Temperature
      K.s = monod_jag$BUGSoutput$summary[1,1],                                                  # Half-saturation constant
      r.max = monod_jag$BUGSoutput$summary[3,1],                                                # Maximum population growth rate
      R.jag = df.jags$phos[which(df.jags$mean > 0.1)[1]],                                      # Minimum resource requirement for positive growth (from jags model)
      R.mth = 0.1*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.1)   # Minimum resource requirement for positive growth (from math)
    ))
    
    phos_sum <- monod_jag$BUGSoutput$summary[c(1:3, (max(df.i$Phosphate_c) - S.pred[1])/0.1 + 9),] # Have to create a new frame for summaries (not listed 1 to 6)
    
    for (j in 1:4){
      fit.P.df <- rbind(fit.P.df, data.frame(      # Model performance data
        Grad = "Phosphorous limitation",         # Abiotic gradient
        Sp.id = i,                               # Species       
        Temp = 20,                               # Temperature
        Parameter = rownames(phos_sum)[j],       # Model parameter (e.g. K_s, r_max, etc.)
        mean = phos_sum[j,1],                    # Posterior mean
        Rhat = phos_sum[j,8],                    # Rhat values
        n.eff = phos_sum[j,9]                  # Sample size estimates (should be ~3000)
      ))
      
    }
}

write.csv(bestion.summ.P.df, "data-processed/11c_Bestion_2018_Monod_phosphorous.csv") # Save Bestion 2018 TPC summary table
write.csv(fit.P.df, "data-processed/11d_Bestion_2018_Monod_phosphorous_fits.csv") # Save model fit summary table
