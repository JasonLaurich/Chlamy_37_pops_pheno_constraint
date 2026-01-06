# Jason R Laurich

# December 29th, 2025

# We are going to reestimate Edwards 2016 TPCs using nls.multstart

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(rTPC)
library(nls.multstart)
library(Deriv)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits

lactin2 <- function(temp, cf.a, cf.tmax, cf.delta_t, cf.b) { 
  exp(cf.a * temp) - exp(cf.a * cf.tmax - (cf.tmax - temp) / cf.delta_t) + cf.b
} # Define the Lactin II function

lactin2_deriv <- function(temp, cf.a, cf.b, cf.tmax, cf.delta_t) {
  rho <- cf.a
  T_max <- cf.tmax
  delta_T <- cf.delta_t
  
  term1 <- rho * exp(rho * temp)
  term2 <- (1 / delta_T) * exp(rho * T_max - (T_max - temp) / delta_T)
  
  return(term1 - term2)
} # Derivative of the Lactin II function

lactin2_halfmax <- function(temp, cf.a, cf.b, cf.tmax, cf.delta_t, r_half) {
  exp(cf.a * temp) - exp(cf.a * cf.tmax - (cf.tmax - temp) / cf.delta_t) + cf.b - r_half
} # OK we're going to modify the function to calculate T_breadth, based on a modified lactin.

# Load & examine the data ------------------------------------------------------------

df <- read.csv("data-processed/15_Edwards_2016_raw_data.csv") # Raw data file
head(df)
str(df)

df$Sp.fac <- as.factor(df$species)

df <- df %>% # The Edwards data has some species that show up in multiple references, so we need to create a unique sp.idx ref# combo
  mutate(unique.id = paste(species, reference, sep = "_")) # So that each entry is treated separately!

df$sp.num <- as.integer(factor(df$unique.id))

df <- df %>% 
  filter(irradiance<= 250) # After looking at the light models, insanely high lights are driving down growth rates. We'll cap this out at 250 to improve model fits.

# Fit TPCs ----------------------------------------------------------------

edwards.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species ID
  Sp.name = character(),    # Species name
  T.min = numeric(),        # Minimum T (calculus)
  T.max = numeric(),        # Maximum T (calculus)
  T.opt = numeric(),        # Optimal T (calculus)
  r.max = numeric(),        # Maximum growth rate (calculus)
  T.br = numeric(),         # T breadth (calculus)
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  est = numeric(),          # Estimate
  se = numeric(),           # Standard error
  p = numeric(),            # p-value
  stringsAsFactors = FALSE                  
)

for (i in c(1:63)){ # For each species, save 38 (no low temp growth data)
  
  df.i <- df %>%
    filter(sp.num == i, !is.na(growth.rate)) %>%
    group_by(temperature) %>%                     # group within each temperature
    slice_max(irradiance, n = 1, with_ties = FALSE) %>%  
    ungroup() %>%
    arrange(temperature)
  
  df.i <- droplevels(df.i)
  
  lac_nls <- nls_multstart(growth.rate ~ lactin2_1995(temp = temperature, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = c(0.01, -2.5, 5, 1),
                           start_upper = c(0.19, -0.5, max(df.i$temperature), 11),
                           lower = c(0, -3, min(df.i$temperature), 0.1),
                           upper = c(0.5, 0, max(df.i$temperature) + 5, 40),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  sum <- summary(lac_nls)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  cf.a <- df.nls[1,1] # Extract parameters
  cf.b <- df.nls[2,1]
  cf.tmax <- df.nls[3,1]
  cf.delta_t <- df.nls[4,1]
  
  # Find the T_opt: where the derivative crosses zero
  T_opt <- uniroot(
    function(temp) lactin2_deriv(temp, cf.a, cf.b, cf.tmax, cf.delta_t),
    interval = c(-10, 45)
  )$root
  
  r_max <- lactin2(temp=T_opt, cf.a=cf.a, cf.b=cf.b, cf.tmax=cf.tmax, cf.delta_t=cf.delta_t)
  
  Tmin <- uniroot(lactin2, interval = c(-1000000000, T_opt), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t)$root
  
  Tmax <- uniroot(lactin2, interval = c(T_opt,45), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t)$root
  
  r_half <- r_max/2 # calculate half of rmax and get the roots.
  
  Tlow <- uniroot(lactin2_halfmax, interval = c(Tmin, T_opt), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t, r_half = r_half)$root
  
  Thigh <- uniroot(lactin2_halfmax, interval = c(T_opt, Tmax), cf.a = cf.a, 
                   cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t, r_half = r_half)$root
  
  edwards.summ.df <- rbind(edwards.summ.df, data.frame( # Add summary data
    Sp.id = df.i$sp.num[1],                           # Species #
    Sp.name = df.i$species[1],                        # Species name
    T.min = Tmin,                                     # Minimum T (calculus)
    T.max = Tmax,                                     # Maximum T (calculus)
    T.opt = T_opt,                                    # Optimal T (calculus)
    r.max = r_max,                                    # Maximum growth rate (calculus)
    T.br = Thigh - Tlow                               # T breadth (calculus)
  ))
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(                         # Model performance data
      Sp.id = df.i$sp.num[1],                                   # Species #
      Sp.name = df.i$species[1],                                # Species name
      Parameter = df.nls$parameter[j],                          # Model parameter (e.g. a, b, tmax etc.)
      est = df.nls$Estimate[j],                                 # Estimate
      se = df.nls$`Std. Error`[j],                              # Error
      p = df.nls$`Pr(>|t|)`[j],                                 # p-values
      stringsAsFactors = FALSE            
    ))
  }
  
}

write.csv(edwards.summ.df, "data-processed/308a_Edwards_2016_TPCs_newest.csv") # let's save the file.
write.csv(fit.df, "data-processed/308b_Edwards_2016_TPCs_fits_newest.csv") # let's save the file.

# Light Monod curves ------------------------------------------------------

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

# Let's do larger models for the final things (10 times larger)
ni.fit <- 330000    # iterations / chain
nb.fit <- 30000     # burn in periods for each chain
nt.fit <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 3         # number of chains, total of 3,000 estimates for each model. 

edwards2016.summ.l.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),                 # Species ID
  Sp.name = character(),             # Species name
  K.s = numeric(),                   # Half-saturation constant
  r.max = numeric(),                 # Maximum population growth rate
  R.jag = numeric(),                 # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),                 # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
  stringsAsFactors = FALSE           # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~3000)
  stringsAsFactors = FALSE            
)

for (i in c(1:63)){ # For each species, save 38 (no low temp growth data)
  
  df.i <- df %>%
    filter(sp.num == i, !is.na(growth.rate)) %>%
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
    model.file = "monod.light.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  print(paste("Done", i))
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3, (max(df.i$irradiance) - S.pred[1])/0.5 + 25 + 4),]   # generate the sequence of r.pred values
  df.jags$light <- seq(S.pred[1], max(df.i$irradiance) + 25, 0.5)
  
  edwards2016.summ.l.df <- rbind(edwards2016.summ.l.df, data.frame(                             # Add summary data
    Sp.id = df.i$sp.num[1],                                                                     # Species #
    Sp.name = df.i$species[1],                                                                  # Species name
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                    # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                  # Maximum population growth rate
    R.jag = df.jags$light[which(df.jags$mean > 0.1)[1]],                                        # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.1*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.1)     # Minimum resource requirement for positive growth (from math)
  ))
  
  light_sum <- monod_jag$BUGSoutput$summary[c(1:3, (max(df.i$irradiance) - S.pred[1])/0.5 + 25 + 4),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = df.i$sp.num[1],                    # Species #
      Sp.name = df.i$species[1],                 # Species name      
      Parameter = rownames(light_sum)[j],        # Model parameter (e.g. K_s, r_max, etc.)
      mean = light_sum[j,1],                     # Posterior mean
      Rhat = light_sum[j,8],                     # Rhat values
      n.eff = light_sum[j,9]                     # Sample size estimates (should be ~3000)
    ))
    
  }
}

write.csv(edwards2016.summ.l.df, "data-processed/308c_Edwards_2016_Monod_light_new.csv") # save the metrics
write.csv(fit.df, "data-processed/308d_Edwards_2016_Monod_light_fits_new.csv") # save the fit data

