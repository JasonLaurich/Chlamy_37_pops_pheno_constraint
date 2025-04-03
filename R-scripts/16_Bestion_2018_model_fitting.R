# Jason R Laurich
# April 3, 2025

# Going to work with the Bestion et al 2018 (https://onlinelibrary.wiley.com/doi/10.1111/ele.12932) data to fit
# Monod curves for phosphorous and TPCs

############# Packages ########################

library(tidyr)
library(cowplot)
library(ggplot2)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)
library(rTPC)
library(nls.multstart)

############# Upload and organize data #######################

df <- read.csv("data-processed/18_bestion_2018_growth_data.csv") # Summary file
head(df)
str(df)

df$Sp.fac <- as.factor(df$SpeciesName)

# So they co-varied phosphorous and temperature
# And they only provide their estimate of µ, so let's move right to fitting models.

# Going to split this data set into 2

# 1. Full phosphorous, variation in temperature. Actually let's go with P30 as well, so we can include Chlamydomonas

# 2. 20 C (consistent with Joey's estimation of P limitation), variation in P

df.t <- df[df$Phosphate_c == 50,]
df.t.p30 <- df[df$Phosphate_c == 30,]

df.t <- droplevels(df.t)
df.t.p30 <-droplevels(df.t.p30)

########### Let's run the TPCs for now #########################################

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

inits.lactin.cust<- function() { # Pulling initial values centres from the start_vals function in rTPC
  list(
    cf.a = rnorm(1, mean = start.vals.lac[1], sd = 0.05),
    cf.tmax = rnorm(1, mean = start.vals.lac[3], sd = 1),
    cf.delta_t = rnorm(1, mean = start.vals.lac[4], sd = 1),
    cf.b = rnorm(1, mean = start.vals.lac[2], sd = 0.05),
    cf.sigma = runif(1, 0.1, 2)
  )
}

Temp.xs <- seq(-5, 45, 0.1) # Temperature gradient we're interested in - upped the granularity here
N.Temp.xs <-length(Temp.xs)

plot.list.p30 <- list() # There's only 5-6 species, we'll plot them all. 

for (i in unique(df.t.p30$Sp.fac)){ # Most species don't have a rich enough dataset at P50, so we're starting with 30
  
  df.i <- df.t.p30[df.t.p30$Sp.fac == i, ]
  
  trait <- df.i$mu      # format the data for jags
  N.obs <- length(trait)
  temp <- df.i$Temperature_c
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  start.vals.lac <- get_start_vals(df.i$Temperature_c, df.i$mu, model_name = 'lactin2_1995')
  
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
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)[-c(1:6),]   # generate the sequence of r.pred values
  df.jags$temp <- seq(-5, 45, 0.1)

  p <- ggplot(data = df.jags, aes(x = temp)) +
    geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
    geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
    geom_point(data = df.i, aes(x = jitter(Temperature, 0.5), y = Growth.rate), color = "grey9", size = 2) + # Add observed data points with jitter for Temp
    scale_x_continuous(limits = c(-5, 45)) + 
    scale_y_continuous(limits = c(-2, 5)) + # Customize the axes and labels +
    labs(
      x = expression(paste("Temperature (", degree, "C)")),
      y = "Growth rate") +
    theme_classic() +
    geom_hline(yintercept = 0)
    
  plot.list.p50[[paste0("Species ", i)]] <- p
  
  thomas.summ.df <- rbind(thomas.summ.df, data.frame(                                  # Add summary data
    Sp.id = i,                                                                         # Species name 
    DIC = lac_jag$BUGSoutput$DIC,                                                      # DIC
    P.conc = 30,                                                                       # Phosphorous concentration
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
      P.conc = 30,                                              # Phosphorous
      Parameter = rownames(lac_jag$BUGSoutput$summary)[j],      # Model parameter (e.g. cf.a, cf.tmax, etc.)
      mean = lac_jag$BUGSoutput$summary[j,1],                   # Posterior mean
      Rhat = lac_jag$BUGSoutput$summary[j,8],                   # Rhat values
      n.eff = lac_jag$BUGSoutput$summary[j,9],                  # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
}

write.csv(thomas.summ.df, "data-processed/17_Thomas2012_TPCs.csv") # Save Thomas 2012 summary table
write.csv(fit.df, "data-processed/17a_Thomas2012_TPCs_fits.csv") # Save model fit summary table
