# Jason R Laurich
# April 14, 2025

# Going to work with the Kontopoulos et al 2020 (https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000894#abstract0) 
# data to fit TPCs for a bunch of species

# Load packages -----------------------------------------------------------

library(cowplot)
library(tidyverse)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)
library(rTPC)
library(nls.multstart)

# Load and examine data ---------------------------------------------------

df <- read.csv("data-processed/21_Kontopoulos2020_growth_data.csv") # Raw data file
head(df)
str(df)

df$Sp.fac <- as.factor(df$Species_standardised)

# TPC fitting -------------------------------------------------------------

# Set up model settings and parameters
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

kontopoulos.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),                 # Species ID (unique)
  Sp.name = character(),             # Species name
  DIC = numeric(),                   # DIC
  T.min.raw = numeric(),             # Minimum T (Jags raw)
  T.max.raw = numeric(),             # Maximum T (Jags raw)
  T.opt.raw = numeric(),             # Optimal T (Jags raw)
  r.max.raw = numeric(),             # Maximum growth rate (Jags raw)
  T.br.raw = numeric(),              # T breadth (Jags raw)
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

df <- df %>%
  mutate(unique.id = paste(Species_standardised, Reference, sep = "_")) # So that each entry is treated seperately!

for (i in unique(df$unique.id)[686:length(unique(df$unique.id))]) { # This allows me to start part-way through the list
  
  df.i <- df %>%
    filter(unique.id == i, !is.na(Trait_value)) %>%
    arrange(Temperature)
  
  max.temp <- max(df.i$Temperature, na.rm = TRUE)
  temp.at.max.growth <- df.i$Temperature[which.max(df.i$Trait_value)]
  
  if (temp.at.max.growth == max.temp) {
    new_row <- df.i %>%
      filter(Temperature == max.temp) %>%
      mutate(
        Temperature = max.temp + 5,
        Trait_value = 0
      )
    
    df.i <- bind_rows(df.i, new_row)
  }
  
  trait <- df.i$Trait_value     # format the data for jags
  N.obs <- length(trait)
  
  Temp.xs <- seq(min(df.i$Temperature) - 5, max(df.i$Temperature) + 5, 0.1) # Temperature gradient we're interested in - upped the granularity here
  N.Temp.xs <-length(Temp.xs) # We'll reset this internally since the gradient varies substantially
  
  temp <- df.i$Temperature
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  start.vals.lac <- get_start_vals(df.i$Temperature, df.i$Trait_value, model_name = 'lactin2_1995')
  
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
  df.jags$temp <- seq(min(df.i$Temperature) - 5, max(df.i$Temperature) + 5, 0.1)
  
  kontopoulos.summ.df <- rbind(kontopoulos.summ.df, data.frame(                        # Add summary data
    Sp.id = i,                                                                         # Species id (unique)
    Sp.name = df.i$Species_standardised[1],                                            # Species name
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
  
  print(nrow(kontopoulos.summ.df))
}

write.csv(kontopoulos.summ.df, "data-processed/21a_Kontopoulos2020_TPCs.csv") # Save Kontopoulos 2020 TPC summary table
write.csv(fit.df, "data-processed/21b_Kontopoulos2020_TPCs_fits.csv") # Save model fit summary table
