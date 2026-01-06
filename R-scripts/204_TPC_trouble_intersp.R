# Jason R Laurich

# December 28th, 2025

# A more detailed/ organized exploration of my TPC fit problems for interspecific datasets.

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(R2jags)
library(mcmcplots)
library(cowplot)
library(rTPC)
library(HDInterval)
library(nls.multstart)

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

df.t.fit <- read.csv('data-processed/10b_Thomas_2012_TPCs_fits.csv') # Thomas TPC fits
head(df.t.fit) # Don't look promising

df.t.fit <- df.t.fit %>%
  select(Sp.id, Parameter, mean) %>%
  pivot_wider(names_from = Parameter, values_from = mean)

# Examining initial values ------------------------------------------------

###### Existing model fit ######

i <- 1 # A problematic population, based on using the get_start_values() function to center inits etc.

tpc.t <- df.t.fit %>% filter(Sp.id == i)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = tpc.t$cf.a, b = tpc.t$cf.b, delta_t = tpc.t$cf.delta_t, tmax = tpc.t$cf.tmax)
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Thomas TPC : existing model"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

###### Checking initial values ######

df.i <- df.t.raw %>%
  filter(id.number == i, !is.na(Growth.rate)) %>% 
  arrange(Temperature)

start.vals.lac <- get_start_vals(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995')

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
    title = "Thomas TPC: init values (get_start...)"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

# So these do not look good....

###### Running models with universal priors and inits ######

inits.lactin <- function() { # Old starting values 
  list(
    cf.a = runif(1, 0.05, 0.15),  # More constrained initial values
    cf.tmax = runif(1, 0, 43),
    cf.delta_t = runif(1, 0.1, 15),
    cf.b = runif(1, -2.5, -1),
    cf.sigma = runif(1, 0.1, 2)
  )
}

ni.fit <- 33000    # iterations / chain
nb.fit <- 3000     # burn in periods for each chain
nt.fit <- 30       # thinning interval : (33,000 - 3,000) / 30 = 1000 posterior estimates / chain
nc.fit <- 3        # number of chains, total of 3,000 estimates for each model. 

parameters.lactin2 <- c("cf.a", "cf.b", "cf.tmax", "cf.delta_t", "cf.sigma", "r.pred") # repeated here

Temp.xs <- seq(-5, 45, 0.1) # Temperature gradient we're interested in - upped the granularity here
N.Temp.xs <-length(Temp.xs)

trait <- df.i$Growth.rate     # format the data for jags
N.obs <- length(trait)
temp <- df.i$Temperature

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

lac_jag1 <- jags(
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

lac_jag1$BUGSoutput$summary[1:6,] # we are having issues with chain convergence AND n.eff
mcmcplot(lac_jag1) # Awful, unsurprisingly

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
    title = "Thomas TPC: universal priors & inits"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

### 10 x

ni.fit2 <- 330000    # iterations / chain
nb.fit2 <- 30000     # burn in periods for each chain
nt.fit2 <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit2 <- 3         # number of chains, total of 3,000 estimates for each model. 

lac_jag2 <- jags(
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
) # ~ 1 min to run?

lac_jag2$BUGSoutput$summary[1:6,] # n.eff is still quite small.
mcmcplot(lac_jag2) # Looks better?

df.jags <- data.frame(lac_jag2$BUGSoutput$summary)

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
    title = "Thomas TPC: universal priors & inits, 10x"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

### 100 x

ni.fit3 <- 3300000    # iterations / chain
nb.fit3 <- 300000     # burn in periods for each chain
nt.fit3 <- 3000       # thinning interval : (3,300,000 - 300,000) / 3000 = 1000 posterior estimates / chain
nc.fit3 <- 3          # number of chains, total of 3,000 estimates for each model. 

lac_jag3 <- jags(
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
) # ~ longer to run

lac_jag3$BUGSoutput$summary[1:6,] # n.eff is still small, but better
mcmcplot(lac_jag3) # Looks better?

df.jags <- data.frame(lac_jag3$BUGSoutput$summary)

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
    title = "Thomas TPC: universal priors & inits, 100x"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

###### Look at parameter posteriors ######

d.a <- lac_jag3$BUGSoutput$sims.list$cf.a  
d.a <- as.numeric(d.a)
d.a.df <- data.frame(stat = d.a)

dens.a <- density(d.a)

p.a <- ggplot(d.a.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.a",
       x = "cf.a",
       y = "Density") +
  geom_vline(xintercept = mean(d.a), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.a$x[which.max(dens.a$y)], linetype = "dashed", color = "black", size = 1)

d.b <- lac_jag3$BUGSoutput$sims.list$cf.b  
d.b <- as.numeric(d.b)
d.b.df <- data.frame(stat = d.b)

dens.b <- density(d.b)

p.b <- ggplot(d.b.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.b",
       x = "cf.b",
       y = "Density") +
  geom_vline(xintercept = mean(d.b), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.b$x[which.max(dens.b$y)], linetype = "dashed", color = "black", size = 1)

d.tmax <- lac_jag3$BUGSoutput$sims.list$cf.tmax  
d.tmax <- as.numeric(d.tmax)
d.tmax.df <- data.frame(stat = d.tmax)

dens.tmax <- density(d.tmax)

p.tm <- ggplot(d.tmax.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.tmax",
       x = "cf.tmax",
       y = "Density") +
  geom_vline(xintercept = mean(d.tmax), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.tmax$x[which.max(dens.tmax$y)], linetype = "dashed", color = "black", size = 1)

d.t <- lac_jag3$BUGSoutput$sims.list$cf.delta_t 
d.t <- as.numeric(d.t)
d.t.df <- data.frame(stat = d.t)

dens.t <- density(d.t)

p.dt <- ggplot(d.t.df, aes(x = stat)) +
  geom_density(size = 1, color = "blue") +
  theme_minimal(base_size = 14) +
  labs(title = "Posterior density for cf.dt",
       x = "cf.dt",
       y = "Density") +
  geom_vline(xintercept = mean(d.t), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = dens.t$x[which.max(dens.t$y)], linetype = "dashed", color = "black", size = 1)

plot_grid(p.a, p.b, p.tm, p.dt,
          align= 'hv',
          nrow = 2)

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
    title = "Thomas TPC: universal priors & inits, 100x, modes"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

###### NLS? ######

start.vals.lac <- get_start_vals(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995')

lac_nls <- nls_multstart(Growth.rate ~ lactin2_1995(temp = Temperature, a, b, tmax, delta_t),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.lac - 10,
                         start_upper = start.vals.lac + 10,
                         lower = get_lower_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995'),
                         upper = get_upper_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

sum <- summary(lac_nls)

df.nls <- as.data.frame(sum$coefficients)
df.nls$parameter <- rownames(df.nls)
rownames(df.nls) <- NULL

df.nls

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.nls[1,1], b = df.nls[2,1], delta_t = df.nls[4,1], tmax = df.nls[3,1])
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Thomas TPC: nls, broad initials"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

# Validation: species 22 --------------------------------------------------

i <- 22

df.i <- df.t.raw %>%
  filter(id.number == i, !is.na(Growth.rate)) %>% 
  arrange(Temperature)

trait <- df.i$Growth.rate     # format the data for jags
N.obs <- length(trait)
temp <- df.i$Temperature

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

lac_jag22 <- jags(
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
) # ~ 1 min to run?

lac_jag22$BUGSoutput$summary[1:6,] # looks pretty good!
mcmcplot(lac_jag22) # Looks better?

df.jags <- data.frame(lac_jag22$BUGSoutput$summary)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.jags[1,1], b = df.jags[2,1], delta_t = df.jags[3,1], tmax = df.jags[5,1])
)

p.22.mean <- ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Population 22 - posterior means"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

d.a <- lac_jag22$BUGSoutput$sims.list$cf.a  
d.a <- as.numeric(d.a)
d.a.df <- data.frame(stat = d.a)

dens.a <- density(d.a)

d.b <- lac_jag22$BUGSoutput$sims.list$cf.b  
d.b <- as.numeric(d.b)
d.b.df <- data.frame(stat = d.b)

dens.b <- density(d.b)

d.tmax <- lac_jag22$BUGSoutput$sims.list$cf.tmax  
d.tmax <- as.numeric(d.tmax)
d.tmax.df <- data.frame(stat = d.tmax)

dens.tmax <- density(d.tmax)

d.t <- lac_jag22$BUGSoutput$sims.list$cf.delta_t 
d.t <- as.numeric(d.t)
d.t.df <- data.frame(stat = d.t)

dens.t <- density(d.t)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = dens.a$x[which.max(dens.a$y)], b = dens.b$x[which.max(dens.b$y)], delta_t = dens.t$x[which.max(dens.t$y)], tmax = dens.tmax$x[which.max(dens.tmax$y)])
)

p.22.mode <- ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Population 22 - posterior modes"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

start.vals.lac <- get_start_vals(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995')

lac_nls22 <- nls_multstart(Growth.rate ~ lactin2_1995(temp = Temperature, a, b, tmax, delta_t),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.lac - 10,
                         start_upper = start.vals.lac + 10,
                         lower = get_lower_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995'),
                         upper = get_upper_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

sum <- summary(lac_nls22)

df.nls <- as.data.frame(sum$coefficients)
df.nls$parameter <- rownames(df.nls)
rownames(df.nls) <- NULL

df.nls

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.nls[1,1], b = df.nls[2,1], delta_t = df.nls[4,1], tmax = df.nls[3,1])
)

p.22.nls <- ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Population 22 - nls"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

plot_grid(p.22.mean, p.22.mode, p.22.nls,
          align= 'hv',
          nrow = 1)

# Validation: species 40 --------------------------------------------------

i <- 40

df.i <- df.t.raw %>%
  filter(id.number == i, !is.na(Growth.rate)) %>% 
  arrange(Temperature)

trait <- df.i$Growth.rate     # format the data for jags
N.obs <- length(trait)
temp <- df.i$Temperature

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

lac_jag40 <- jags(
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
) # ~ 1 min to run?

lac_jag40$BUGSoutput$summary[1:6,] # Looks ok (except tmax)
mcmcplot(lac_jag40) # Looks better?

df.jags <- data.frame(lac_jag40$BUGSoutput$summary)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.jags[1,1], b = df.jags[2,1], delta_t = df.jags[3,1], tmax = df.jags[5,1])
)

p.40.mean <- ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Population 40 - posterior means"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

d.a <- lac_jag40$BUGSoutput$sims.list$cf.a  
d.a <- as.numeric(d.a)
d.a.df <- data.frame(stat = d.a)

dens.a <- density(d.a)

d.b <- lac_jag40$BUGSoutput$sims.list$cf.b  
d.b <- as.numeric(d.b)
d.b.df <- data.frame(stat = d.b)

dens.b <- density(d.b)

d.tmax <- lac_jag40$BUGSoutput$sims.list$cf.tmax  
d.tmax <- as.numeric(d.tmax)
d.tmax.df <- data.frame(stat = d.tmax)

dens.tmax <- density(d.tmax)

d.t <- lac_jag40$BUGSoutput$sims.list$cf.delta_t 
d.t <- as.numeric(d.t)
d.t.df <- data.frame(stat = d.t)

dens.t <- density(d.t)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = dens.a$x[which.max(dens.a$y)], b = dens.b$x[which.max(dens.b$y)], delta_t = dens.t$x[which.max(dens.t$y)], tmax = dens.tmax$x[which.max(dens.tmax$y)])
)

p.40.mode <- ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Population 40 - posterior modes"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

start.vals.lac <- get_start_vals(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995')

lac_nls40 <- nls_multstart(Growth.rate ~ lactin2_1995(temp = Temperature, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = start.vals.lac - 10,
                           start_upper = start.vals.lac + 10,
                           lower = get_lower_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995'),
                           upper = get_upper_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995'),
                           supp_errors = 'Y',
                           convergence_count = FALSE
)

sum <- summary(lac_nls40)

df.nls <- as.data.frame(sum$coefficients)
df.nls$parameter <- rownames(df.nls)
rownames(df.nls) <- NULL

df.nls

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.nls[1,1], b = df.nls[2,1], delta_t = df.nls[4,1], tmax = df.nls[3,1])
)

p.40.nls <- ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Population 40 - nls"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

plot_grid(p.40.mean, p.40.mode, p.40.nls,
          align= 'hv',
          nrow = 1)

# Validation: species 45 --------------------------------------------------

i <- 45

df.i <- df.t.raw %>%
  filter(id.number == i, !is.na(Growth.rate)) %>% 
  arrange(Temperature)

trait <- df.i$Growth.rate     # format the data for jags
N.obs <- length(trait)
temp <- df.i$Temperature

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

lac_jag45 <- jags(
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
) # ~ 1 min to run?

lac_jag45$BUGSoutput$summary[1:6,] # Looks great
mcmcplot(lac_jag45) # Looks better?

df.jags <- data.frame(lac_jag45$BUGSoutput$summary)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.jags[1,1], b = df.jags[2,1], delta_t = df.jags[3,1], tmax = df.jags[5,1])
)

p.45.mean <- ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Population 45 - posterior means"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

d.a <- lac_jag45$BUGSoutput$sims.list$cf.a  
d.a <- as.numeric(d.a)
d.a.df <- data.frame(stat = d.a)

dens.a <- density(d.a)

d.b <- lac_jag45$BUGSoutput$sims.list$cf.b  
d.b <- as.numeric(d.b)
d.b.df <- data.frame(stat = d.b)

dens.b <- density(d.b)

d.tmax <- lac_jag45$BUGSoutput$sims.list$cf.tmax  
d.tmax <- as.numeric(d.tmax)
d.tmax.df <- data.frame(stat = d.tmax)

dens.tmax <- density(d.tmax)

d.t <- lac_jag45$BUGSoutput$sims.list$cf.delta_t 
d.t <- as.numeric(d.t)
d.t.df <- data.frame(stat = d.t)

dens.t <- density(d.t)

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = dens.a$x[which.max(dens.a$y)], b = dens.b$x[which.max(dens.b$y)], delta_t = dens.t$x[which.max(dens.t$y)], tmax = dens.tmax$x[which.max(dens.tmax$y)])
)

p.45.mode <- ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Population 45 - posterior modes"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

start.vals.lac <- get_start_vals(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995')

lac_nls45 <- nls_multstart(Growth.rate ~ lactin2_1995(temp = Temperature, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = start.vals.lac - 10,
                           start_upper = start.vals.lac + 10,
                           lower = get_lower_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995'),
                           upper = get_upper_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995'),
                           supp_errors = 'Y',
                           convergence_count = FALSE
)

sum <- summary(lac_nls45)

df.nls <- as.data.frame(sum$coefficients)
df.nls$parameter <- rownames(df.nls)
rownames(df.nls) <- NULL

df.nls

curve.t <- tibble::tibble(
  res  = seq(-10, 50, length.out = 200),
  rate = pred_lact(res, a = df.nls[1,1], b = df.nls[2,1], delta_t = df.nls[4,1], tmax = df.nls[3,1])
)

p.45.nls <- ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = Temperature, y = Growth.rate),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Population 45 - nls"
  ) +
  theme_classic() +
  ylim(-0.1, 4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

plot_grid(p.45.mean, p.45.mode, p.45.nls,
          align= 'hv',
          nrow = 1)
