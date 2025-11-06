# Jason R Laurich
# November 6th, 2025

# Rerunning populations TPCs with the new growth data

# Will then compare my results to Joey's 2020 paper. 

# Packages & functions ----------------------------------------------------

library(tidyr)
library(dplyr)
library(R2jags)
library(mcmcplots)

# Load & examine the data ------------------------------------------------------------

df.r <- read.csv("data-processed/202_µ_estimates_light_new.csv")

head(df.r)
str(df.r)

light <- as.vector(as.numeric(as.character(unique(df.r$percentage)))) # for looping through light levels
ord.light<- sort(light)

df.r$pop.fac <- as.factor(df.r$population)
df.r$pop.num <- as.numeric(df.r$population.number)
df.r$well.ID <- as.factor(df.r$well.ID)
df.r$light <- as.numeric(df.r$percentage)

mat <- split(df.r, df.r$pop.num)  # Matrixify the data!

# Monod curves ------------------------------------------------------------

inits.monod.final <- function() { # In case I want to play with these in the future
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Repeated here

S.pred <- seq(0, 100, 0.05) # Light gradient we're interested in (percentages)
N.S.pred <-length(S.pred)

ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model.

# Summary dfs

summary.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Pop.fac = character(),    # Population name
  Pop.num = character(),    # Number assigned to population (not the same). This corresponds to the jags objects
  K.s = numeric(),          # Half-saturation constant
  r.max = numeric(),        # Maximum population growth rate
  R.jag = numeric(),        # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),        # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Grad = character(),       # Specifiy abiotic gradient
  Pop.fac = character(),    # Population name
  Pop.num = character(),    # Number assigned to population (not the same)       
  Parameter = character(),  # Model parameter (e.g. K_s, r_max, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~6000)
  stringsAsFactors = FALSE            
)

for (i in 1:length(mat)){ # for each population.
  
  df.i <- subset(mat[[i]])
  df.i <- droplevels(df.i)
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  light <- df.i$light
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = light, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the light Monod function. 
    data = jag.data,
    inits = inits.monod.final,
    parameters.to.save = parameters.monod,
    model.file = "monod.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  save(monod_jag, file = paste0("R2jags-objects/pop_", i, "_light_monod_new.RData")) # save the L limitation monod function
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3,2005),]   # generate the sequence of r.pred values
  df.jags$light <- seq(0, 100, 0.05)
  
  summary.df <- rbind(summary.df, data.frame(                                                   # Add summary data
    Pop.fac = df.i$pop.fac[1],                                                                  # Population name
    Pop.num = df.i$pop.num[1],                                                                  # Number assigned to population (not the same)
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                    # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                  # Maximum population growth rate
    R.jag = df.jags$light[which(df.jags$mean > 0.56)[1]],                                       # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.56*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.56)   # Minimum resource requirement for positive growth (from math)                                                   
  ))
  
  light_sum <- monod_jag$BUGSoutput$summary[c(1:3, 2005),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(        # Model performance data
      Grad = "Light limitation",               # Abiotic gradient
      Pop.fac = df.i$pop.fac[1],               # Population name
      Pop.num = df.i$pop.num[1],               # Number assigned to population (not the same)       
      Parameter = rownames(light_sum)[j],      # Model parameter (e.g. K_s, r_max, etc.)
      mean = light_sum[j,1],                   # Posterior mean
      Rhat = light_sum[j,8],                   # Rhat values
      n.eff = light_sum[j,9],                  # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
  
  message(sprintf("Done %d of 37", i))
  
}

write.csv(summary.df, "data-processed/214_Monod_light_pops_estimates.csv") # Save summary table
write.csv(fit.df, "data-processed/215_Monod_light_pops_fits.csv") # Save model fit summary table
