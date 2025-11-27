# Jason R Laurich
# November 26th, 2025

# OK we're chasing down the final few models that are still not performing.  
# Anc 4 - rep 1 - # 113/148
# 8 - rep 2 - # 98/148?

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(R2jags)
library(mcmcplots)
library(cowplot)
library(rTPC)

pred_lact <- function(temp, a, b, delta_t, tmax) {
  exp(a * temp) - exp(a * tmax - ((tmax - temp) / delta_t)) + b
}

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

df.r <- read.csv("data-processed/204_µ_estimates_temp_new.csv")

head(df.r)
str(df.r)

df.r$pop.fac <- as.factor(df.r$population)
df.r$pop.num <- as.numeric(df.r$population.number)
df.r$well.ID <- as.factor(df.r$well.ID)
df.r$temp <- as.numeric(df.r$temperature)

df.r <- df.r %>%
  mutate(unique.id = str_c(population.number, str_sub(well.ID, 1, 3), sep = "."))

length(unique(df.r$unique.id)) # 152, should be 1480/10 = 142

unique(df.r$pop.fac) # cc1629 is in there, that's why

df.r <- df.r %>%
  filter(!(pop.fac == "cc1629"))

df.r <- droplevels(df.r)

unique(df.r$pop.fac)

length(unique(df.r$unique.id)) # 148, should be 1480/10 = 148! Fixed for now. 

mat <- split(df.r, df.r$unique.id)  # Matrixify the data!

# OK let's start running models -------------------------------------------

inits.lactin <- function() { # Old starting values 
  list(
    cf.a = runif(1, 0.05, 0.15),  # More constrained initial values
    cf.tmax = runif(1, 37, 43),
    cf.delta_t = runif(1, 1, 5),
    cf.b = runif(1, -2.5, -1),
    cf.sigma = runif(1, 0.1, 2)
  )
}

# Let's set more generous MCMC settings still for our final models. 
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model. 

parameters.lactin2 <- c("cf.a", "cf.b", "cf.tmax", "cf.delta_t", "cf.sigma", "r.pred") # repeated here

Temp.xs <- seq(0, 45, 0.05) # Temperature gradient we're interested in - upped the granularity here
N.Temp.xs <-length(Temp.xs)

# Pop 8, rep 2 ------------------------------------------------------------

i <- 98

df.i <- subset(mat[[i]])
df.i <- droplevels(df.i)

trait <- df.i$r.exp     # format the data for jags
N.obs <- length(trait)
temp <- df.i$temp

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

###### Current model ######

lac_jag_9 <- jags(          # With the current priors that worked on 146/148 populations
  data = jag.data, 
  inits = inits.lactin, # choose one or the other
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin.txt",
  n.thin = nt.fit, 
  n.chains = nc.fit, 
  n.burnin = nb.fit, 
  n.iter = ni.fit, 
  DIC = TRUE, 
  working.directory = getwd()
)

head(lac_jag_9$BUGSoutput$summary)
mcmcplot(lac_jag_9)

save(lac_jag_9, file = "R2jags-objects/rep_98_lactin_newer.RData") # save the lactin2 model

###### Plotting ######

t.mu <- df.r %>% filter(unique.id == "31.D03")

curve.t <- tibble::tibble(
  res  = seq(0, 45, length.out = 200),
  rate = pred_lact(res, a = lac_jag_9$BUGSoutput$summary[1,1], b = lac_jag_9$BUGSoutput$summary[2,1], 
                   delta_t = lac_jag_9$BUGSoutput$summary[3,1], tmax = lac_jag_9$BUGSoutput$summary[5,1])
)

p.8.2 <- ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = t.mu,
             aes(x = temp, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "A: existing priors"
  ) +
  theme_classic() +
  ylim(-0.1, 6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.8.2

# OK that looks good! Let's move onto anc 4 - rep 1

# Anc 4, rep 1 ------------------------------------------------------------

i <- 113

df.i <- subset(mat[[i]])
df.i <- droplevels(df.i)

trait <- df.i$r.exp     # format the data for jags
N.obs <- length(trait)
temp <- df.i$temp

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

###### Current model ######

lac_jag_anc <- jags(          # With the current priors that worked on 146/148 populations
  data = jag.data, 
  inits = inits.lactin, # choose one or the other
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin.txt",
  n.thin = nt.fit, 
  n.chains = nc.fit, 
  n.burnin = nb.fit, 
  n.iter = ni.fit, 
  DIC = TRUE, 
  working.directory = getwd()
)

head(lac_jag_anc$BUGSoutput$summary)
mcmcplot(lac_jag_anc)

save(lac_jag_anc, file = "R2jags-objects/rep_113_lactin_newer.RData") # save the lactin2 model

###### Plotting ######

t.mu <- df.r %>% filter(unique.id == "35.B11")

curve.t <- tibble::tibble(
  res  = seq(0, 45, length.out = 200),
  rate = pred_lact(res, a = lac_jag_anc$BUGSoutput$summary[1,1], b = lac_jag_anc$BUGSoutput$summary[2,1], 
                   delta_t = lac_jag_anc$BUGSoutput$summary[3,1], tmax = lac_jag_anc$BUGSoutput$summary[5,1])
)

p.anc.1 <- ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = t.mu,
             aes(x = temp, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "B: existing priors"
  ) +
  theme_classic() +
  ylim(-0.1, 6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.anc.1

plot_grid(p.8.2, p.anc.1)

# Assembly into tables (only if models run in chunks) --------------------------------------

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

for (i in 1:length(mat)){ # for each replicated sample, can adjust if code crashes
  
  load(paste0("R2jags-objects/rep_", i, "_lactin_newer.RData"))
  
  df.i <- subset(mat[[i]])
  df.i <- droplevels(df.i)
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)[-c(1:6),]   # generate the sequence of r.pred values
  df.jags$temp <- seq(0, 45, 0.05)
  
  cf.a <- lac_jag$BUGSoutput$summary[1,1] # Extract parameters
  cf.b <- lac_jag$BUGSoutput$summary[2,1]
  cf.tmax <- lac_jag$BUGSoutput$summary[5,1]
  cf.delta_t <- lac_jag$BUGSoutput$summary[3,1]
  
  # Find the T_opt: where the derivative crosses zero
  T_opt <- uniroot(
    function(temp) lactin2_deriv(temp, cf.a, cf.b, cf.tmax, cf.delta_t),
    interval = c(10, 45)
  )$root
  
  r_max <- lactin2(temp=T_opt, cf.a=cf.a, cf.b=cf.b, cf.tmax=cf.tmax, cf.delta_t=cf.delta_t)
  
  Tmin <- uniroot(lactin2, interval = c(-5, T_opt), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t)$root
  
  Tmax <- uniroot(lactin2, interval = c(T_opt,45), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t)$root
  
  r_half <- r_max/2 # calculate half of rmax and get the roots.
  
  Tlow <- uniroot(lactin2_halfmax, interval = c(Tmin, T_opt), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t, r_half = r_half)$root
  
  Thigh <- uniroot(lactin2_halfmax, interval = c(T_opt, Tmax), cf.a = cf.a, 
                   cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t, r_half = r_half)$root
  
  summary.df <- rbind(summary.df, data.frame(         # Add summary data
    Pop.fac = df.i$pop.fac[1],                        # Population name
    Pop.num = df.i$pop.num[1],                        # Number assigned to population (not the same)
    unique.id = df.i$unique.id[1],                    # Unique id (jag objects)
    Model = "Lactin 2",                               # Model name
    T.min = Tmin,                                     # Minimum T (calculus)
    T.max = Tmax,                                     # Maximum T (calculus)
    T.br = Thigh - Tlow,                              # T breadth (calculus)
    T.opt = T_opt,                                    # Optimal T (calculus)
    r.max = r_max                                     # Maximum growth rate (calculus)
  ))
  
  for (j in 1:6){
    fit.df <- rbind(fit.df, data.frame(                         # Model performance data
      Model = "Lactin 2",                                       # Model name
      Pop.fac = df.i$pop.fac[1],                                # Population name
      Pop.num = df.i$pop.num[1],                                # Number assigned to population (not the same) 
      unique.id = df.i$unique.id[1],                            # Unique id (jag objects)
      Parameter = rownames(lac_jag$BUGSoutput$summary)[j],      # Model parameter (e.g. cf.a, cf.tmax, etc.)
      mean = lac_jag$BUGSoutput$summary[j,1],                   # Posterior mean
      Rhat = lac_jag$BUGSoutput$summary[j,8],                   # Rhat values
      n.eff = lac_jag$BUGSoutput$summary[j,9],                  # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
  
}

write.csv(summary.df, "data-processed/303a_TPC_estimates_final.csv") # Save summary table
write.csv(fit.df, "data-processed/303b_TPC_fits_final.csv") # Save model fit summary table

