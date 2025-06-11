# Jason R Laurich

# June 9, 2025

# Going to upload data from the Edwards et al 2016 Limnol Oceanogr paper and run light Monods and TPCs on the raw data.


# Load packages -----------------------------------------------------------

library(tidyr)
library(cowplot)
library(tidyverse)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)
library(rTPC)
library(nls.multstart)

# Load and examine data ---------------------------------------------------

df <- read.csv("data-processed/27a_Edwards_2016_raw_data.csv") # Raw data file
head(df)
str(df)

df$Sp.fac <- as.factor(df$species)

# Lactin II TPCs ----------------------------------------------------------

edwards2016.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),               # Species ID
  DIC = numeric(),                 # DIC
  T.min.raw = numeric(),           # Minimum T (Jags raw)
  T.max.raw = numeric(),           # Maximum T (Jags raw)
  T.opt.raw = numeric(),           # Optimal T (Jags raw)
  r.max.raw = numeric(),           # Maximum growth rate (Jags raw)
  T.br.raw = numeric(),            # T breadth (Jags raw)
  stringsAsFactors = FALSE         # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
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

inits.lactin.cust<- function() { # Pulling initial values centres from the start_vals function in rTPC
  list(
    cf.a = rnorm(1, mean = start.vals.lac[1], sd = 0.05),
    cf.tmax = rnorm(1, mean = start.vals.lac[3], sd = 1),
    cf.delta_t = rnorm(1, mean = start.vals.lac[4], sd = 1),
    cf.b = rnorm(1, mean = start.vals.lac[2], sd = 0.05),
    cf.sigma = runif(1, 0.1, 2)
  )
}

for (i in unique(df$Sp.fac)){ # For each species
  
  df.i <- df %>%
    filter(Sp.fac == i, !is.na(growth.rate)) %>%
    filter(irradiance == max(irradiance, na.rm = TRUE)) %>%
    arrange(temperature)
  
  max.temp <- max(df.i$temp, na.rm = TRUE)
  temp.at.max.growth <- df.i$temp[which.max(df.i$growth.rate)]
  
  if (temp.at.max.growth == max.temp) {
    new_row <- df.i %>%
      filter(temperature == max.temp) %>%
      mutate(
        temperature = max.temp + 5,
        growth.rate = 0
      )
    
    df.i <- bind_rows(df.i, new_row)
  }
  
  trait <- df.i$growth.rate     # format the data for jags
  N.obs <- length(trait)
  
  Temp.xs <- seq(min(df.i$temperature) - 5, max(df.i$temperature) + 5, 0.1) # Temperature gradient we're interested in - upped the granularity here
  N.Temp.xs <-length(Temp.xs) # We'll reset this internally since the gradient varies substantially
  
  temp <- df.i$temperature
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  start.vals.lac <- get_start_vals(df.i$temperature, df.i$growth.rate, model_name = 'lactin2_1995')
  
  lac_jag <- jags(
    data = jag.data, 
    inits = inits.lactin.cust, 
    parameters.to.save = parameters.lactin2, 
    model.file = "lactin2_thomas.txt",
    n.thin = nt.fit, 
    n.chains = nc.fit, 
    n.burnin = nb.fit, 
    n.iter = ni.fit, 
    DIC = TRUE, 
    working.directory = getwd()
  ) # ~ 10 min to run?
  
  print(paste("Done", i))
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)[-c(1:6),]   # generate the sequence of r.pred values
  df.jags$temp <- seq(min(df.i$temperature) - 5, max(df.i$temperature) + 5, 0.1)
  
  edwards2016.summ.df <- rbind(edwards2016.summ.df, data.frame(                        # Add summary data
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

edwards2016.summ.df %>% 
  filter(rowSums(is.na(.)) > 0) %>% 
  print()

# These 4 populations didn't work, because their light values are not consistent.
# Go back and finish them manually? Maybe later

write.csv(edwards2016.summ.df, "data-processed/27b_Edwards2016_TPCs.csv") # let's save the file.
write.csv(fit.df, "data-processed/27c_Edwards2016_TPCs_fits.csv") # let's save the file.


# Light Monod curves ------------------------------------------------------

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

edwards2016.summ.l.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),                 # Species ID
  DIC = numeric(),                   # DIC
  K.s = numeric(),                   # Half-saturation constant
  r.max = numeric(),                 # Maximum population growth rate
  R.jag = numeric(),                 # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),                 # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
  stringsAsFactors = FALSE           # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~3000)
  stringsAsFactors = FALSE            
)

for (i in unique(df$Sp.fac)){ # For each species
  
  df.i <- df %>%
    filter(Sp.fac == i, !is.na(growth.rate)) %>%
    filter(temperature == temperature[which.min(abs(temperature - 20))]) %>%
    arrange(irradiance)
  
  
  trait <- df.i$growth.rate     # format the data for jags
  N.obs <- length(trait)
  
  light <- df.i$irradiance
  
  S.pred <- seq(min(df.i$irradiance) - 25, max(df.i$irradiance) + 25, 0.5) # Light gradient we're interested in - upped the granularity here
  S.pred <- S.pred[S.pred >= 0]  # Cut off values below 0
  N.S.pred <-length(S.pred) # We'll reset this internally since the gradient varies substantially
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = light, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the light Monod function. 
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
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3, (max(df.i$irradiance) - S.pred[1])/0.5 + 4),]   # generate the sequence of r.pred values
  df.jags$light <- seq(S.pred[1], max(df.i$irradiance) + 25, 0.5)
  
  edwards2016.summ.l.df <- rbind(edwards2016.summ.l.df, data.frame(                                 # Add summary data
    Sp.id = i,                                                                                  # Species name 
    DIC = monod_jag$BUGSoutput$DIC,                                                             # DIC
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                    # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                  # Maximum population growth rate
    R.jag = df.jags$light[which(df.jags$mean > 0.1)[1]],                                       # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.1*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.1)   # Minimum resource requirement for positive growth (from math)
  ))
  
  light_sum <- monod_jag$BUGSoutput$summary[c(1:3, (max(df.i$irradiance) - S.pred[1])/0.5 + 4),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = i,                                 # Species       
      Parameter = rownames(light_sum)[j],        # Model parameter (e.g. K_s, r_max, etc.)
      mean = light_sum[j,1],                     # Posterior mean
      Rhat = light_sum[j,8],                     # Rhat values
      n.eff = light_sum[j,9]                     # Sample size estimates (should be ~3000)
    ))
    
  }
}

write.csv(edwards2016.summ.l.df, "data-processed/27d_Edwards2016_light.csv") # let's save the file.
write.csv(fit.df, "data-processed/27e_Edwards2016_light_fits.csv") # let's save the file.
