# Jason R Laurich
# November 7th, 2025

# Rerunning populations salt tolerance curves with the new growth data

# Will then compare my results to Joey's 2020 paper. 

# Packages & functions ----------------------------------------------------

library(tidyr)
library(dplyr)
library(R2jags)
library(mcmcplots)

# Load & examine the data ------------------------------------------------------------

df.r <- read.csv("data-processed/203_µ_estimates_salt_new.csv")

head(df.r)
str(df.r)

df.r$pop.fac <- as.factor(df.r$population)
df.r$pop.num <- as.numeric(df.r$population.number)
df.r$well.ID <- as.factor(df.r$well.ID)
df.r$salt <- as.numeric(df.r$salt)

salt <- as.vector(as.numeric(as.character(unique(df.r$salt)))) # for looping through salt levels
ord.salt<- sort(salt)

mat <- split(df.r, df.r$pop.num)  # Matrixify the data!

# Fit tolerance curves ----------------------------------------------------

inits.salt <- function() { # Smaller a (prior 0.5-> 2)
  list(
    a = runif(1, 0.7, 1.8),  # Smaller window for a
    b = runif(1, 0.1, 5),  # Initial guess for b
    c = runif(1, 0.1, max(df.i$salt)),  # Initial guess for c
    sigma = runif(1, 0.1, 2)  # Initial guess for error
  )
}

parameters.salt <- c("a", "b", "c", "sigma", "r_pred_new") # Repeats

S.pred <- seq(0, 10, 0.005) # Repeats
N.S.pred <-length(S.pred)

# Repeated here
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model.

# Summary dfs

summary.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Pop.fac = character(),    # Population name
  Pop.num = character(),    # Number assigned to population (not the same). This corresponds to the jags objects
  r.max = numeric(),        # Maximum population growth rate (alpha)
  c.mod = numeric(),        # salt concentration at which r is half of alpha (extracted from model)
  c.pred = numeric(),       # salt concentration at which r is half of alpha (extracted from predicted values)
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
  salt <- df.i$salt
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the salt logistic growth curve function. 
    data = jag.data,
    inits = inits.salt,
    parameters.to.save = parameters.salt,
    model.file = "salt.tol.smalla.smallc.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  save(monod_jag, file = paste0("R2jags-objects/pop_", i, "_salt_final.RData")) # save the S tolerance function
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:4,2006),]   # generate the sequence of r.pred values
  df.jags$salt <- seq(0, 10, 0.005)
  
  summary.df <- rbind(summary.df, data.frame(                                                                 # Add summary data
    Pop.fac = df.i$pop.fac[1],                                                                                # Population name
    Pop.num = df.i$pop.num[1],                                                                                # Number assigned to population (not the same)
    r.max = monod_jag$BUGSoutput$summary[1,1],                                                                # Maximum population growth rate
    c.mod = monod_jag$BUGSoutput$summary[3,1],                                                                # salt concentration at which r is half of alpha (extracted from model)
    c.pred = df.jags$salt[which.min(abs(df.jags$mean - (monod_jag$BUGSoutput$summary[1,1] / 2)))]   # salt concentration at which r is half of alpha (extracted from predicted values)                                                   
  ))
  
  salt_sum <- monod_jag$BUGSoutput$summary[c(1:4, 2006),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(        # Model performance data
      Grad = "Salt stress",                    # Abiotic gradient
      Pop.fac = df.i$pop.fac[1],               # Population name
      Pop.num = df.i$pop.num[1],               # Number assigned to population (not the same)       
      Parameter = rownames(salt_sum)[j],       # Model parameter (e.g. K_s, r_max, etc.)
      mean = salt_sum[j,1],                    # Posterior mean
      Rhat = salt_sum[j,8],                    # Rhat values
      n.eff = salt_sum[j,9],                   # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
  
  message(sprintf("Done %d of 37", i))
  
}

write.csv(summary.df, "data-processed/216_salt_pops_estimates_final.csv") # Save summary table
write.csv(fit.df, "data-processed/217_salt_pops_fits_final.csv") # Save model fit summary table
