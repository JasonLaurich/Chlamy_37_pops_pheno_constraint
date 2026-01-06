# Jason R Laurich

# December 18th, 2025

# Trying to figure out how to fix my interspecific TPC fits!

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(R2jags)
library(mcmcplots)
library(cowplot)
library(rTPC)
library(HDInterval)

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

df.t.raw <- read.csv('data-processed/10_Thomas_2012_raw_data.csv') # Thomas raw data
head(df.t.raw)

length(unique(df.t.raw$id.number))

mat <- split(df.t.raw, df.t.raw$id.number)  # Matrix

min(df.t.raw[df.t.raw$Growth.rate>0,]$Temperature) # Tmin => -1.8
max(df.t.raw[df.t.raw$Growth.rate>0,]$Temperature) # Tmax <= 37

# Examining initial values ------------------------------------------------

i <- 1

df.i <- df.t.raw %>%
  filter(id.number == i, !is.na(Growth.rate)) %>% 
  arrange(Temperature)

start.vals.lac <- get_start_vals(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995')

inits.lactin.meta<- function() { # Pulling initial values centres from the start_vals function in rTPC
  list(
    cf.a = rnorm(1, mean = start.vals.lac[1], sd = 0.05),
    cf.tmax = rnorm(1, mean = start.vals.lac[3], sd = 1),
    cf.delta_t = rnorm(1, mean = start.vals.lac[4], sd = 1),
    cf.b = rnorm(1, mean = start.vals.lac[2], sd = 0.05),
    cf.sigma = runif(1, 0.1, 2)
  )
}

n.iter <- 100 # I want to take a look at what these distributions look like!

init.vals <- replicate(
  n.iter,
  inits.lactin.meta(),
  simplify = TRUE
)

df.init <- as.data.frame(t(init.vals))
df.init <- df.init[, c("cf.a", "cf.b", "cf.tmax", "cf.delta_t")]

df.init <- as.data.frame(lapply(df.init, unlist))

str(df.init$cf.a)

p.a <-ggplot(df.init, aes(x=cf.a)) +
  geom_histogram(bins=20)

p.a

p.b <-ggplot(df.init, aes(x=cf.b)) +
  geom_histogram(bins=20)

p.b

p.tmax <-ggplot(df.init, aes(x=cf.tmax)) +
  geom_histogram(bins=20)

p.tmax

p.delta_t <-ggplot(df.init, aes(x=cf.delta_t)) +
  geom_histogram(bins=20)

p.delta_t

# OK let's try a different initialization

inits.lactin.meta2<- function() { # Pulling initial values centres from the start_vals function in rTPC
  list(
    cf.a = runif(1, start.vals.lac[1]*0.9, start.vals.lac[1]*1.1),
    cf.tmax = runif(1, start.vals.lac[3]*0.9, start.vals.lac[3]*1.1),
    cf.delta_t = runif(1, start.vals.lac[4]*0.9, start.vals.lac[4]*1.1),
    cf.b = runif(1, start.vals.lac[2]*1.1, start.vals.lac[2]*0.9),
    cf.sigma = runif(1, 0.1, 2)
  )
}

# OK we're going to try a constant range now

inits.lactin <- function() { # Old starting values 
  list(
    cf.a = runif(1, 0.05, 0.15),  # More constrained initial values
    cf.tmax = runif(1, 0, 43),
    cf.delta_t = runif(1, 0.1, 15),
    cf.b = runif(1, -2.5, -1),
    cf.sigma = runif(1, 0.1, 2)
  )
}

# Let's test out a couple populations? ------------------------------------

ni.fit <- 33000    # iterations / chain
nb.fit <- 3000     # burn in periods for each chain
nt.fit <- 30       # thinning interval : (33,000 - 3,000) / 30 = 1000 posterior estimates / chain
nc.fit <- 3        # number of chains, total of 3,000 estimates for each model. 

parameters.lactin2 <- c("cf.a", "cf.b", "cf.tmax", "cf.delta_t", "cf.sigma", "r.pred") # repeated here

Temp.xs <- seq(-5, 45, 0.1) # Temperature gradient we're interested in - upped the granularity here
N.Temp.xs <-length(Temp.xs)

df.i <- df.t.raw %>%
  filter(id.number == i, !is.na(Growth.rate)) %>% 
  arrange(Temperature)

trait <- df.i$Growth.rate     # format the data for jags
N.obs <- length(trait)
temp <- df.i$Temperature

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

lac_jag1 <- jags(
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
) # ~ 1 min to run?

lac_jag1$BUGSoutput$summary[1:6,] # terrible with the expanded cf.a and delta_t
mcmcplot(lac_jag1) # Looks terrible

df.jags <- data.frame(lac_jag1$BUGSoutput$summary)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.jags[1,1], b = df.jags[2,1], delta_t = df.jags[3,1], tmax = df.jags[5,1])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Thomas 2012: Lactin II TPC"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

### OK let's try the second set of inits?

lac_jag1.2 <- jags(
  data = jag.data, 
  inits = inits.lactin.meta2, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin_meta.txt",
  n.thin = nt.fit, 
  n.chains = nc.fit, 
  n.burnin = nb.fit, 
  n.iter = ni.fit, 
  DIC = TRUE, 
  working.directory = getwd()
) # ~ 1 min to run?

lac_jag1.2$BUGSoutput$summary[1:6,] # hmm still looks bad
mcmcplot(lac_jag1.2) # Looks terrible

df.jags <- data.frame(lac_jag1.2$BUGSoutput$summary)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.jags[1,1], b = df.jags[2,1], delta_t = df.jags[3,1], tmax = df.jags[5,1])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Thomas 2012: Lactin II TPC"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

### Okay I think we need to address the priors as well...

a.min <- start.vals.lac[1]*0.75
a.max <- start.vals.lac[1]*1.25

b.min <- start.vals.lac[2]*1.25
b.max <- start.vals.lac[2]*0.75

tmax.min <- start.vals.lac[3]*0.75
tmax.max <- start.vals.lac[3]*1.25

dt.min <- start.vals.lac[4]*0.75
dt.max <- start.vals.lac[4]*1.25

model_string <- sprintf("
  model {

    cf.a       ~ dunif(%f, %f)
    cf.tmax    ~ dunif(%f, %f)
    cf.delta_t ~ dunif(%f, %f)
    cf.b       ~ dunif(%f, %f)
    cf.sigma   ~ dunif(0, 5)
    cf.tau <- pow(cf.sigma, -2)

    for (i in 1:N.obs) {
      trait.mu[i] <- exp(cf.a * temp[i]) -
                     exp(cf.a * cf.tmax - (cf.tmax - temp[i]) / cf.delta_t) +
                     cf.b
      trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }

    for (j in 1:N.Temp.xs) {
      r.pred[j] <- exp(cf.a * Temp.xs[j]) -
                   exp(cf.a * cf.tmax - (cf.tmax - Temp.xs[j]) / cf.delta_t) +
                   cf.b
    }
  }
  ", a.min, a.max, tmax.min, tmax.max, dt.min, dt.max, b.min, b.max)

writeLines(model_string, con = "lactin_model.txt")

lac_jag1.3 <- jags(
  data = jag.data, 
  inits = inits.lactin.meta2, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin_model.txt",
  n.thin = nt.fit, 
  n.chains = nc.fit, 
  n.burnin = nb.fit, 
  n.iter = ni.fit, 
  DIC = TRUE, 
  working.directory = getwd()
) # ~ 1 min to run?

lac_jag1.3$BUGSoutput$summary[1:6,] # hmm that looks OK!
mcmcplot(lac_jag1.3) # Looks better?

df.jags <- data.frame(lac_jag1.3$BUGSoutput$summary)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = start.vals.lac[1], b = start.vals.lac[2], delta_t = start.vals.lac[4], tmax = start.vals.lac[3])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Thomas 2012: Lactin II TPC"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

lac_jag1.4 <- jags(
  data = jag.data, 
  inits = inits.lactin, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin.thomas.txt",
  n.thin = nt.fit, 
  n.chains = nc.fit, 
  n.burnin = nb.fit, 
  n.iter = ni.fit, 
  DIC = TRUE, 
  working.directory = getwd()
) # ~ 1 min to run?

lac_jag1.4$BUGSoutput$summary[1:6,] # hmm that looks OK!
mcmcplot(lac_jag1.4) # Looks better?

df.jags <- data.frame(lac_jag1.4$BUGSoutput$summary)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.jags[1,1], b = df.jags[2,1], delta_t = df.jags[3,1], tmax = df.jags[5,1])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Jags 1.4: mean values"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

# Let's take a closer look at delta_t:

d.t <- lac_jag1.4$BUGSoutput$sims.list$cf.delta_t  

d.t <- as.numeric(d.t)

d.t.df <- data.frame(stat = d.t)

ggplot(d.t.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.delta_t",
       x = "cf.delta_t",
       y = "Density") +
  geom_vline(xintercept = df.jags[3,1], linetype = "dashed", color = "red", size = 1)


# What about the median delta_t value? 
median(d.t)
mean(d.t)
hdi(d.t, credMass = 0.95)
hdi(d.t, credMass = 0.33)
mean(hdi(d.t, credMass = 0.33))
hdi(d.t, credMass = 0.1)
mean(hdi(d.t, credMass = 0.1))

ggplot(d.t.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.delta_t",
       x = "cf.delta_t",
       y = "Density") +
  geom_vline(xintercept = df.jags[3,1], linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = mean(hdi(d.t, credMass = 0.25)), linetype = "dashed", color = "black", size = 1)

# What about the mode of the density distribution directly?

dens <- density(d.t)
mode_delta <- dens$x[which.max(dens$y)]
mode_delta

ggplot(d.t.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.delta_t",
       x = "cf.delta_t",
       y = "Density") +
  geom_vline(xintercept = df.jags[3,1], linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens$x[which.max(dens$y)], linetype = "dashed", color = "black", size = 1)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.jags[1,1], b = df.jags[2,1], delta_t = dens$x[which.max(dens$y)], tmax = df.jags[5,1])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Jags 1.4: mean values, mode d.t"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

# What about all modes?

d.a <- lac_jag1.4$BUGSoutput$sims.list$cf.a  
d.a <- as.numeric(d.a)
d.a.df <- data.frame(stat = d.a)

dens.a <- density(d.a)

ggplot(d.a.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.a",
       x = "cf.a",
       y = "Density") +
  geom_vline(xintercept = mean(d.a), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.a$x[which.max(dens.a$y)], linetype = "dashed", color = "black", size = 1)

d.b <- lac_jag1.4$BUGSoutput$sims.list$cf.b  
d.b <- as.numeric(d.b)
d.b.df <- data.frame(stat = d.b)

dens.b <- density(d.b)

ggplot(d.b.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.b",
       x = "cf.b",
       y = "Density") +
  geom_vline(xintercept = mean(d.b), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.b$x[which.max(dens.b$y)], linetype = "dashed", color = "black", size = 1)

d.tmax <- lac_jag1.4$BUGSoutput$sims.list$cf.tmax  
d.tmax <- as.numeric(d.tmax)
d.tmax.df <- data.frame(stat = d.tmax)

dens.tmax <- density(d.tmax)

ggplot(d.tmax.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.tmax",
       x = "cf.tmax",
       y = "Density") +
  geom_vline(xintercept = mean(d.tmax), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.tmax$x[which.max(dens.tmax$y)], linetype = "dashed", color = "black", size = 1)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = dens.a$x[which.max(dens.a$y)], b = dens.b$x[which.max(dens.b$y)], delta_t = dens$x[which.max(dens$y)], tmax = dens.tmax$x[which.max(dens.tmax$y)])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Jags 1.4: all modes"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

#### OK let's try rerunning the model with these HDPIs as initials?
hd.a <- hdi(lac_jag1.4$BUGSoutput$sims.list$cf.a, credMass = 0.33)
hd.b <- hdi(lac_jag1.4$BUGSoutput$sims.list$cf.b, credMass = 0.33)
hd.tmax <- hdi(lac_jag1.4$BUGSoutput$sims.list$cf.tmax, credMass = 0.33)
hd.dt <- hdi(lac_jag1.4$BUGSoutput$sims.list$cf.delta_t, credMass = 0.33)

inits.lactin.2bayes <- function() { # Old starting values 
  list(
    cf.a = runif(1, hd.a[1], hd.a[2]),  # More constrained initial values
    cf.tmax = runif(1, hd.tmax[1], hd.tmax[2]),
    cf.delta_t = runif(1, hd.dt[1], hd.dt[2]),
    cf.b = runif(1, hd.b[1], hd.b[2]),
    cf.sigma = runif(1, 0.1, 2)
  )
}

lac_jag1.5 <- jags(
  data = jag.data, 
  inits = inits.lactin.2bayes, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin.thomas.txt",
  n.thin = nt.fit, 
  n.chains = nc.fit, 
  n.burnin = nb.fit, 
  n.iter = ni.fit, 
  DIC = TRUE, 
  working.directory = getwd()
) # ~ 1 min to run?

lac_jag1.5$BUGSoutput$summary[1:6,] # hmm that looks OK!
mcmcplot(lac_jag1.5) # Looks better? Nope

### What if I change the priors too?

hd.a2 <- hdi(lac_jag1.4$BUGSoutput$sims.list$cf.a, credMass = 0.36)
hd.b2 <- hdi(lac_jag1.4$BUGSoutput$sims.list$cf.b, credMass = 0.36)
hd.tmax2 <- hdi(lac_jag1.4$BUGSoutput$sims.list$cf.tmax, credMass = 0.36)
hd.dt2 <- hdi(lac_jag1.4$BUGSoutput$sims.list$cf.delta_t, credMass = 0.36)

model_string <- sprintf("
  model {

    cf.a       ~ dunif(%f, %f)
    cf.tmax    ~ dunif(%f, %f)
    cf.delta_t ~ dunif(%f, %f)
    cf.b       ~ dunif(%f, %f)
    cf.sigma   ~ dunif(0, 5)
    cf.tau <- pow(cf.sigma, -2)

    for (i in 1:N.obs) {
      trait.mu[i] <- exp(cf.a * temp[i]) -
                     exp(cf.a * cf.tmax - (cf.tmax - temp[i]) / cf.delta_t) +
                     cf.b
      trait[i] ~ dnorm(trait.mu[i], cf.tau)
    }

    for (j in 1:N.Temp.xs) {
      r.pred[j] <- exp(cf.a * Temp.xs[j]) -
                   exp(cf.a * cf.tmax - (cf.tmax - Temp.xs[j]) / cf.delta_t) +
                   cf.b
    }
  }
  ", hd.a2[1], hd.a2[2], hd.tmax2[1], hd.tmax2[2], hd.dt2[1], hd.dt2[2], hd.b2[1], hd.b2[2])

writeLines(model_string, con = "lactin_model.txt")

lac_jag1.6 <- jags(
  data = jag.data, 
  inits = inits.lactin.2bayes, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin_model.txt",
  n.thin = nt.fit, 
  n.chains = nc.fit, 
  n.burnin = nb.fit, 
  n.iter = ni.fit, 
  DIC = TRUE, 
  working.directory = getwd()
) # ~ 1 min to run?

lac_jag1.6$BUGSoutput$summary[1:6,] # hmm that looks OK!
mcmcplot(lac_jag1.6) # Looks better? Nope

df.jags <- data.frame(lac_jag1.6$BUGSoutput$summary)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.jags[1,1], b = df.jags[2,1], delta_t = df.jags[3,1], tmax = df.jags[5,1])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Thomas 2012: Lactin II TPC"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

# Let's try running the 1.4 model with larger specifications?

ni.fit2 <- 330000    # iterations / chain
nb.fit2 <- 30000     # burn in periods for each chain
nt.fit2 <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit2 <- 3         # number of chains, total of 3,000 estimates for each model. 

lac_jag1.7 <- jags(
  data = jag.data, 
  inits = inits.lactin, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin.thomas.txt",
  n.thin = nt.fit2, 
  n.chains = nc.fit2, 
  n.burnin = nb.fit2, 
  n.iter = ni.fit2, 
  DIC = TRUE, 
  working.directory = getwd()
) # ~ 10 min to run?

lac_jag1.7$BUGSoutput$summary[1:6,] # hmm that looks OK!
mcmcplot(lac_jag1.7) # Looks better?

df.jags <- data.frame(lac_jag1.7$BUGSoutput$summary)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.jags[1,1], b = df.jags[2,1], delta_t = df.jags[3,1], tmax = df.jags[5,1])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "1.7: bigger run"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

# What about all modes?

d.a <- lac_jag1.7$BUGSoutput$sims.list$cf.a  
d.a <- as.numeric(d.a)
d.a.df <- data.frame(stat = d.a)

dens.a <- density(d.a)

ggplot(d.a.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.a",
       x = "cf.a",
       y = "Density") +
  geom_vline(xintercept = mean(d.a), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.a$x[which.max(dens.a$y)], linetype = "dashed", color = "black", size = 1)

d.b <- lac_jag1.7$BUGSoutput$sims.list$cf.b  
d.b <- as.numeric(d.b)
d.b.df <- data.frame(stat = d.b)

dens.b <- density(d.b)

ggplot(d.b.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.b",
       x = "cf.b",
       y = "Density") +
  geom_vline(xintercept = mean(d.b), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.b$x[which.max(dens.b$y)], linetype = "dashed", color = "black", size = 1)

d.tmax <- lac_jag1.7$BUGSoutput$sims.list$cf.tmax  
d.tmax <- as.numeric(d.tmax)
d.tmax.df <- data.frame(stat = d.tmax)

dens.tmax <- density(d.tmax)

ggplot(d.tmax.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.tmax",
       x = "cf.tmax",
       y = "Density") +
  geom_vline(xintercept = mean(d.tmax), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.tmax$x[which.max(dens.tmax$y)], linetype = "dashed", color = "black", size = 1)

d.t <- lac_jag1.7$BUGSoutput$sims.list$cf.delta_t 
d.t <- as.numeric(d.t)
d.t.df <- data.frame(stat = d.t)

dens.t <- density(d.t)

ggplot(d.t.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.dt",
       x = "cf.dt",
       y = "Density") +
  geom_vline(xintercept = mean(d.t), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.t$x[which.max(dens.t$y)], linetype = "dashed", color = "black", size = 1)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = dens.a$x[which.max(dens.a$y)], b = dens.b$x[which.max(dens.b$y)], delta_t = dens.t$x[which.max(dens.t$y)], tmax = dens.tmax$x[which.max(dens.tmax$y)])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Jags 1.7: all modes, larger"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

# So it looks like size is not the issue here, but let's try ramping up 10x one more time just to confirm

ni.fit3 <- 3300000    # iterations / chain
nb.fit3 <- 300000     # burn in periods for each chain
nt.fit3 <- 3000       # thinning interval : (3,300,000 - 300,000) / 3000 = 1000 posterior estimates / chain
nc.fit3 <- 3          # number of chains, total of 3,000 estimates for each model. 

lac_jag1.8 <- jags(
  data = jag.data, 
  inits = inits.lactin, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin.thomas.txt",
  n.thin = nt.fit3, 
  n.chains = nc.fit3, 
  n.burnin = nb.fit3, 
  n.iter = ni.fit3, 
  DIC = TRUE, 
  working.directory = getwd()
) # ~ 10 min to run?

lac_jag1.8$BUGSoutput$summary[1:6,] # hmm that looks OK!
mcmcplot(lac_jag1.8) # Looks better?

df.jags <- data.frame(lac_jag1.8$BUGSoutput$summary)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.jags[1,1], b = df.jags[2,1], delta_t = df.jags[3,1], tmax = df.jags[5,1])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "1.8: bigger run x 10"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

# What about all modes?

d.a <- lac_jag1.8$BUGSoutput$sims.list$cf.a  
d.a <- as.numeric(d.a)
d.a.df <- data.frame(stat = d.a)

dens.a <- density(d.a)

ggplot(d.a.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.a",
       x = "cf.a",
       y = "Density") +
  geom_vline(xintercept = mean(d.a), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.a$x[which.max(dens.a$y)], linetype = "dashed", color = "black", size = 1)

d.b <- lac_jag1.8$BUGSoutput$sims.list$cf.b  
d.b <- as.numeric(d.b)
d.b.df <- data.frame(stat = d.b)

dens.b <- density(d.b)

ggplot(d.b.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.b",
       x = "cf.b",
       y = "Density") +
  geom_vline(xintercept = mean(d.b), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.b$x[which.max(dens.b$y)], linetype = "dashed", color = "black", size = 1)

d.tmax <- lac_jag1.8$BUGSoutput$sims.list$cf.tmax  
d.tmax <- as.numeric(d.tmax)
d.tmax.df <- data.frame(stat = d.tmax)

dens.tmax <- density(d.tmax)

ggplot(d.tmax.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.tmax",
       x = "cf.tmax",
       y = "Density") +
  geom_vline(xintercept = mean(d.tmax), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.tmax$x[which.max(dens.tmax$y)], linetype = "dashed", color = "black", size = 1)

d.t <- lac_jag1.8$BUGSoutput$sims.list$cf.delta_t 
d.t <- as.numeric(d.t)
d.t.df <- data.frame(stat = d.t)

dens.t <- density(d.t)

ggplot(d.t.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.dt",
       x = "cf.dt",
       y = "Density") +
  geom_vline(xintercept = mean(d.t), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.t$x[which.max(dens.t$y)], linetype = "dashed", color = "black", size = 1)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = dens.a$x[which.max(dens.a$y)], b = dens.b$x[which.max(dens.b$y)], delta_t = dens.t$x[which.max(dens.t$y)], tmax = dens.tmax$x[which.max(dens.tmax$y)])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Jags 1.7: all modes, larger"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

# OK my burn-in period seems like it maybe needs to be larger?

ni.fit4 <- 3000000    # iterations / chain
nb.fit4 <- 1000000    # burn in periods for each chain
nt.fit4 <- 2000       # thinning interval : (3,000,000 - 1,000,000) / 2000 = 1000 posterior estimates / chain
nc.fit4 <- 3          # number of chains, total of 3,000 estimates for each model. 

lac_jag1.9 <- jags(
  data = jag.data, 
  inits = inits.lactin, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin.thomas.txt",
  n.thin = nt.fit4, 
  n.chains = nc.fit4, 
  n.burnin = nb.fit4, 
  n.iter = ni.fit4, 
  DIC = TRUE, 
  working.directory = getwd()
) # ~ 10 min to run?

lac_jag1.8$BUGSoutput$summary[1:6,] # hmm that looks OK!
mcmcplot(lac_jag1.8) # Looks better?

df.jags <- data.frame(lac_jag1.8$BUGSoutput$summary)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.jags[1,1], b = df.jags[2,1], delta_t = df.jags[3,1], tmax = df.jags[5,1])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "1.8: bigger run x 10"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

# What about all modes?

d.a <- lac_jag1.8$BUGSoutput$sims.list$cf.a  
d.a <- as.numeric(d.a)
d.a.df <- data.frame(stat = d.a)

dens.a <- density(d.a)

ggplot(d.a.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.a",
       x = "cf.a",
       y = "Density") +
  geom_vline(xintercept = mean(d.a), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.a$x[which.max(dens.a$y)], linetype = "dashed", color = "black", size = 1)

d.b <- lac_jag1.8$BUGSoutput$sims.list$cf.b  
d.b <- as.numeric(d.b)
d.b.df <- data.frame(stat = d.b)

dens.b <- density(d.b)

ggplot(d.b.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.b",
       x = "cf.b",
       y = "Density") +
  geom_vline(xintercept = mean(d.b), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.b$x[which.max(dens.b$y)], linetype = "dashed", color = "black", size = 1)

d.tmax <- lac_jag1.8$BUGSoutput$sims.list$cf.tmax  
d.tmax <- as.numeric(d.tmax)
d.tmax.df <- data.frame(stat = d.tmax)

dens.tmax <- density(d.tmax)

ggplot(d.tmax.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.tmax",
       x = "cf.tmax",
       y = "Density") +
  geom_vline(xintercept = mean(d.tmax), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.tmax$x[which.max(dens.tmax$y)], linetype = "dashed", color = "black", size = 1)

d.t <- lac_jag1.8$BUGSoutput$sims.list$cf.delta_t 
d.t <- as.numeric(d.t)
d.t.df <- data.frame(stat = d.t)

dens.t <- density(d.t)

ggplot(d.t.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.dt",
       x = "cf.dt",
       y = "Density") +
  geom_vline(xintercept = mean(d.t), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.t$x[which.max(dens.t$y)], linetype = "dashed", color = "black", size = 1)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = dens.a$x[which.max(dens.a$y)], b = dens.b$x[which.max(dens.b$y)], delta_t = dens.t$x[which.max(dens.t$y)], tmax = dens.tmax$x[which.max(dens.tmax$y)])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Jags 1.7: all modes, larger"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)
