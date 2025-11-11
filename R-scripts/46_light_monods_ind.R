# Jason R Laurich

# October 6th, 2025

# I am going to fit Monod curve to each unique replicate of light data, so we have more power for statistical analyses. 

# Packages & functions ----------------------------------------------------

library(nls.multstart)
library(tidyverse)
library(gridExtra)
library(grid)
library(R2jags)
library(mcmcplots)

# Load the data -----------------------------------------------------------

df.r <- read.csv("data-processed/202_µ_estimates_light_new.csv")

head(df.r)
str(df.r)

df.r$pop.fac <- as.factor(df.r$population)
df.r$pop.num <- as.numeric(df.r$population.number)
df.r$well.ID <- as.factor(df.r$well.ID)

df <- read.csv("data-processed/06_light_rfus_time.csv") # For extracting light metrics.
head(df) # Percentage is also a measurement of light

df$percentage <- as.numeric(df$percentage) # From an examination of the csv, the light level 2 corresponds to "0.5-0.7"
df$percentage[is.na(df$percentage)] <- 0.6 # For now, let's set this to 0.6, but I need to talk with Joey about this. 

light <- as.vector(as.numeric(as.character(unique(df$percentage)))) # for looping through light levels
ord.light<- sort(light)

df.r$light.lvl <- as.factor(df.r$light.lvl)
df.r$light.numeric <- ord.light[as.numeric(df.r$light.lvl)] # Match up the numeric light level values.

df.r <- df.r %>%
  mutate(unique.id = str_c(population.number, str_sub(well.ID, 1, 3), sep = "."))

length(unique(df.r$unique.id)) # 148, and it is!

by_pop <- df.r %>%
  distinct(population.number, unique.id) %>%   # one row per (pop, unique.id)
  count(population.number, name = "n_unique")  # how many unique.id per pop

by_pop # all good

mat <- split(df.r, df.r$unique.id)  # Matrixify the data!

# Looping through all the replicates --------------------------------------

summary.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Pop.fac = character(),    # Population name
  Pop.num = character(),    # Number assigned to population (not the same). 
  unique.id = character(),  # Unique id (jag objects)
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
  unique.id = character(),  # Unique id (jag objects)
  Parameter = character(),  # Model parameter (e.g. K_s, r_max, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~6000)
  stringsAsFactors = FALSE            
)

inits.monod.final <- function() { # In case I want to play with these in the future
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Repeated here

S.pred <- seq(0, 100, 0.05) # Repeated here
N.S.pred <-length(S.pred)

# Repeated here
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model.

for (i in 1:length(mat)){ # for each replicated sample, can adjust if code crashes
  
  df.i <- subset(mat[[i]])
  df.i <- droplevels(df.i)
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  light <- df.i$percentage
  
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
  
  save(monod_jag, file = paste0("R2jags-objects/rep_", i, "_light_monod_new.RData")) # save the light limitation monod function
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3,2005),]   # generate the sequence of r.pred values
  df.jags$light <- seq(0, 100, 0.05)
  
  summary.df <- rbind(summary.df, data.frame(                                                   # Add summary data
    Pop.fac = df.i$pop.fac[1],                                                                  # Population name
    Pop.num = df.i$pop.num[1],                                                                  # Number assigned to population (not the same)
    unique.id = df.i$unique.id[1],                                                              # Unique id
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
      unique.id = df.i$unique.id[1],           # Unique id
      Parameter = rownames(light_sum)[j],      # Model parameter (e.g. K_s, r_max, etc.)
      mean = light_sum[j,1],                   # Posterior mean
      Rhat = light_sum[j,8],                   # Rhat values
      n.eff = light_sum[j,9],                  # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
  
  message(sprintf("Done %d of 148", i))
  
}

write.csv(summary.df, "data-processed/300a_Monod_light_estimates_new.csv") # Save summary table
write.csv(fit.df, "data-processed/300b_Monod_light_fits_new.csv") # Save model fit summary table
