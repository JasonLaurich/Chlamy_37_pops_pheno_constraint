# Jason R Laurich

# September 30th, 2025

# OK I am going to chase down some potential reasons that the CIs in my models (TPCs, Monods, salt tolerance) are so wide...

# Load packages & specify functions -----------------------------------------------------------

library(nls.multstart)
library(tidyverse)
library(cowplot)
library(rTPC)
library(MuMIn)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)
library(investr)

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

# Load & examine the data -------------------------------------------------

df <- read.csv("data-processed/01_µ_estimates_temp.csv")

head(df)
str(df)

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$population.number)
df$well.ID <- as.factor(df$well.ID)
df$Temp <- df$temperature

mat <- split(df, df$pop.num)  # Matrixify the data!

df.med <- mat[[33]] # Going to pull out the median population in terms of thermal breadth at 0.56

head(df.med)

# TPC investigation -------------------------------------------------------

# We're going to try a couple of things first:

# 1. My initials, or using nls.multstart?
# 2. Expanding the priors? Combine that with 1?
# 3. Number of iterations?
# 4. How does this compare to using rTPC?

parameters.lactin2 <- c("cf.a", "cf.b", "cf.tmax", "cf.delta_t", "cf.sigma", "r.pred") # parameters

Temp.xs <- seq(0, 45, 0.05) # Temperature gradient we're interested in - upped the granularity here
N.Temp.xs <-length(Temp.xs)

# The MCMC settings I used:
ni.1 <- 330000   # iterations / chain
nb.1 <- 30000    # burn in periods for each chain
nt.1 <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.1 <- 6        # number of chains, total of 6,000 estimates for each model. 

# ~ 10x larger MCMC settings?
ni.2 <- 3300000   # iterations / chain
nb.2 <- 300000    # burn in periods for each chain
nt.2 <- 300      # thinning interval : (3,300,000 - 300,000) / 300 = 10,000 posterior estimates / chain
nc.2 <- 10        # number of chains, total of 6,000 estimates for each model. 

inits.lactin <- function() { # My initials
  list(
    cf.a = runif(1, 0.05, 0.15),  # More constrained initial values
    cf.tmax = runif(1, 37, 43),
    cf.delta_t = runif(1, 1, 5),
    cf.b = runif(1, -2.5, -1),
    cf.sigma = runif(1, 0.1, 2)
  )
}

start.vals.lac <- get_start_vals(df.med$temp, df.med$r.exp, model_name = 'lactin2_1995')

inits.lactin.meta<- function() { # Pulling initial values centres from the start_vals function in rTPC
  list(
    cf.a = rnorm(1, mean = start.vals.lac[1], sd = 0.05),
    cf.tmax = rnorm(1, mean = start.vals.lac[3], sd = 1),
    cf.delta_t = rnorm(1, mean = start.vals.lac[4], sd = 1),
    cf.b = rnorm(1, mean = start.vals.lac[2], sd = 0.05),
    cf.sigma = runif(1, 0.1, 2)
  )
}

###### Let's start running the models ######

# 1 : my original model

trait <- df.med$r.exp     # format the data for jags
N.obs <- length(trait)
temp <- df.med$Temp

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

lac_jag.1 <- jags(
  data = jag.data, 
  inits = inits.lactin, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin.txt",
  n.thin = nt.1, 
  n.chains = nc.1, 
  n.burnin = nb.1, 
  n.iter = ni.1, 
  DIC = TRUE, 
  working.directory = getwd()
)

save(lac_jag.1, file = "R2jags-objects/test_pop_33_lactin_normal.RData") # save the lactin2 model

# 2 : using the original settings, but using nls.mulstart to get my initial values.

lac_jag.2 <- jags(
  data = jag.data, 
  inits = inits.lactin.meta, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin_meta.txt",
  n.thin = nt.1, 
  n.chains = nc.1, 
  n.burnin = nb.1, 
  n.iter = ni.1, 
  DIC = TRUE, 
  working.directory = getwd()
)

save(lac_jag.2, file = "R2jags-objects/test_pop_33_lactin_small.nlsmult.RData") # save the lactin2 model

# 3 : my original model, but with more iteration

lac_jag.3 <- jags(
  data = jag.data, 
  inits = inits.lactin, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin.txt",
  n.thin = nt.2, 
  n.chains = nc.2, 
  n.burnin = nb.2, 
  n.iter = ni.2, 
  DIC = TRUE, 
  working.directory = getwd()
)

save(lac_jag.3, file = "R2jags-objects/test_pop_33_lactin_huge.RData") # save the lactin2 model

# 4 : using the larger settings, and using nls.mulstart to get my initial values.

lac_jag.4 <- jags(
  data = jag.data, 
  inits = inits.lactin.meta, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin_meta.txt",
  n.thin = nt.2, 
  n.chains = nc.2, 
  n.burnin = nb.2, 
  n.iter = ni.2, 
  DIC = TRUE, 
  working.directory = getwd()
)

save(lac_jag.4, file = "R2jags-objects/test_pop_33_lactin_huge.nlsmult.RData") # save the lactin2 model

# 5 : rTPC model as a comparison

lac_nls.5 <- nls_multstart(r.exp ~ lactin2_1995(temp = Temp, a, b, tmax, delta_t),
                         data = df.med,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.lac - 10,
                         start_upper = start.vals.lac + 10,
                         lower = get_lower_lims(df.med$Temp, df.med$r.exp, model_name = 'lactin2_1995'),
                         upper = get_upper_lims(df.med$Temp, df.med$r.exp, model_name = 'lactin2_1995'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(lac_nls.5)

nd <- data.frame(Temp = seq(0, 45, length.out = 901))  # Bootstrapping variation
pf <- predFit(lac_nls.5, newdata = nd,
              interval = "confidence", level = 0.95)

df.p.5 <- cbind(nd, as.data.frame(pf))   # columns: fit, lwr, upr

###### Load the objects and assign them unique names ######

load("R2jags-objects/test_pop_33_lactin_normal.RData")
df.lac.1 <- data.frame(lac_jag.1$BUGSoutput$summary)
df.p.1 <- df.lac.1[-c(1:6),]
df.p.1$temp <- seq(0, 45, 0.05)

load("R2jags-objects/test_pop_33_lactin_small.nlsmult.RData")
df.lac.2 <- data.frame(lac_jag.2$BUGSoutput$summary)
df.p.2 <- df.lac.2[-c(1:6),]
df.p.2$temp <- seq(0, 45, 0.05)

load("R2jags-objects/test_pop_33_lactin_huge.RData")
df.lac.3 <- data.frame(lac_jag.3$BUGSoutput$summary)
df.p.3 <- df.lac.3[-c(1:6),]
df.p.3$temp <- seq(0, 45, 0.05)

load("R2jags-objects/test_pop_33_lactin_huge.nlsmult.RData")
df.lac.4 <- data.frame(lac_jag.4$BUGSoutput$summary)
df.p.4 <- df.lac.4[-c(1:6),]
df.p.4$temp <- seq(0, 45, 0.05)

# Let's start by looking at plots without error bars — just the models and raw data

p.mod <- ggplot(df.med, aes(x = temperature, y = r.exp)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) +
  
  geom_line(data = df.p.1, aes(x = temp, y= mean), colour = "red3", linewidth = 1.5) +
  geom_line(data = df.p.2, aes(x = temp, y= mean), colour = "goldenrod2", linewidth = 1.5) +
  geom_line(data = df.p.3, aes(x = temp, y= mean), colour = "forestgreen", linewidth = 1.5) +
  geom_line(data = df.p.4, aes(x = temp, y= mean), colour = "dodgerblue3", linewidth = 1.5) +
  geom_line(data = df.p.5, aes(x = Temp, y= .fitted), colour = "mediumorchid", linewidth = 1.5) +
  
  ylim(-1,4) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "All models together, no error") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.mod # OK so it looks like the 3 models using the nls.multstart inits cluster together, and my inits are separate (but similar - those are the red and green ones)

# OK now let's look at plots with error bars

p.1 <- ggplot(df.med, aes(x = temperature, y = r.exp)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) +
  
  geom_line(data = df.p.1, aes(x = temp, y= mean), colour = "red3", linewidth = 1.5) +
  
  geom_ribbon(data = df.p.1,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              inherit.aes = FALSE, fill = "red3", alpha = 0.20) +
  
  ylim(-1,4) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "A - smaller model, fixed initial values") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.1

p.2 <- ggplot(df.med, aes(x = temperature, y = r.exp)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) +
  
  geom_line(data = df.p.2, aes(x = temp, y= mean), colour = "goldenrod2", linewidth = 1.5) +
  
  geom_ribbon(data = df.p.2,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              inherit.aes = FALSE, fill = "goldenrod2", alpha = 0.20) +
  
  ylim(-1,4) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "B - smaller model, flexible inits") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.2

p.3 <- ggplot(df.med, aes(x = temperature, y = r.exp)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) +
  
  geom_line(data = df.p.3, aes(x = temp, y= mean), colour = "forestgreen", linewidth = 1.5) +
  
  geom_ribbon(data = df.p.3,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              inherit.aes = FALSE, fill = "forestgreen", alpha = 0.20) +
  
  ylim(-1,4) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "C - huge model, fixed initial values") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.3

p.4 <- ggplot(df.med, aes(x = temperature, y = r.exp)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) +
  
  geom_line(data = df.p.4, aes(x = temp, y= mean), colour = "dodgerblue3", linewidth = 1.5) +
  
  geom_ribbon(data = df.p.4,
              aes(x = temp, ymin = X2.5., ymax = X97.5.),
              inherit.aes = FALSE, fill = "dodgerblue3", alpha = 0.20) +
  
  ylim(-1,4) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "D - huge model, flexible inits") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.4

p.5 <- ggplot(df.med, aes(x = temperature, y = r.exp)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) +
  
  geom_line(data = df.p.5, aes(x = Temp, y= fit), colour = "mediumorchid", linewidth = 1.5) +
  
  geom_ribbon(data = df.p.5,
              aes(x = Temp, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, fill = "mediumorchid", alpha = 0.20) +
  
  ylim(-1,4) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "E - nls_mulstart model fit") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.5

p.mod.error <- plot_grid(p.mod, p.1, p.2, p.3, p.4, p.5, align = 'hv', nrow = 2)
p.mod.error
