# Jason R Laurich

# January 6th, 2025

# OK, I am going to fit all of the TPCs to all of the relevant datasets: Thomas 2012, Bestion 2018, Lewington-Pearce 2019, Edwards 2016, and Levasseur 2025.
# In one file. We'll integrate data validation when identifying models to fit, and exclude those that fail to meet the following criteria:
# 1: n>= 5
# 2: there must be 1+ data point at higher temperature than that of the putative Topt
# 3: Tmin, as identified by the model must be > -2 C (the approximate freezing point of seawater)
# 4: (Bootstrapped) 95% CIs for Tmin and Tmax must be <10% of Tmax-Tmin

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(rTPC)
library(nls.multstart)
library(Deriv)
library(car)
library(boot)
library(minpack.lm)

library(R2jags)
library(mcmcplots)
library(bayestestR)

lactin2 <- function(temp, a, tmax, d.t, b) { 
  exp(a * temp) - exp(a * tmax - (tmax - temp) / d.t) + b
} # Define the Lactin II function

lactin2_deriv <- function(temp, a, b, tmax, d.t) {
  rho <- a
  T_max <- tmax
  delta_T <- d.t
  
  term1 <- rho * exp(rho * temp)
  term2 <- (1 / delta_T) * exp(rho * T_max - (T_max - temp) / delta_T)
  
  return(term1 - term2)
} # Derivative of the Lactin II function

lactin2_halfmax <- function(temp, a, b, tmax, d.t, r_half) {
  exp(a * temp) - exp(a * tmax - (tmax - temp) / d.t) + b - r_half
} # OK we're going to modify the function to calculate T_breadth, based on a modified lactin.

# Thomas 2012 ------------------------------------------------------------

###### Load the data ######

df.t.raw <- read.csv('data-processed/10_Thomas_2012_raw_data.csv') # Thomas raw data
head(df.t.raw)

length(unique(df.t.raw$id.number))

min(df.t.raw[df.t.raw$Growth.rate>0,]$Temperature) # Tmin => -1.8
max(df.t.raw[df.t.raw$Growth.rate>0,]$Temperature) # Tmax <= 37

df.t.raw %>% # Number of observations by species
  group_by(id.number) %>%
  summarise(n = n(), .groups = "drop") %>% 
  print(n = 200)

df.t.raw %>% # 24 species are lacking sufficient data. 
  group_by(id.number) %>%
  summarise(n = n(), .groups = "drop") %>% 
  filter(n <5) %>% 
  print(n = 200)

df.t <- df.t.raw %>% # Remove species with fewer than 5 data points
  group_by(id.number) %>%
  mutate(n = n()) %>% 
  filter(n >= 5) 

length(unique(df.t$id.number)) # Down to 170 out of 194 species. 

df.t <- df.t %>% # Remove species where there are no growth data above putative Topt. 
  group_by(id.number) %>%
  filter(!(max(Temperature) %in% Temperature[Growth.rate == max(Growth.rate)])) 

length(unique(df.t$id.number)) # Down to 169 out of 194 species.

df.t <-df.t %>% 
  rename(temp = Temperature,
         mu = Growth.rate)

###### Model fitting ######

thomas.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Study = character(),      # Study
  
  T.min = numeric(),        # Minimum T (calculus)
  T.min.min = numeric(),    # Minimum T (lower HDPI)
  T.min.max = numeric(),    # Minimum T (upper HDPI)
  T.min.na = numeric(),     # %% NA returns
  
  T.max = numeric(),        # Maximum T (calculus)
  T.max.min = numeric(),    # Maximum T (lower HDPI)
  T.max.max = numeric(),    # Maximum T (upper HDPI)
  T.max.na = numeric(),     # %% NA returns
  
  T.opt = numeric(),        # Optimal T (calculus)
  T.opt.min = numeric(),    # Optimal T (lower HDPI)
  T.opt.max = numeric(),    # Optimal T (upper HDPI)
  
  r.max = numeric(),        # Maximum growth rate (calculus)
  r.max.min = numeric(),    # Maximum growth rate (lower HDPI)
  r.max.max = numeric(),    # Maximum growth rate (upper HDPI)
  
  T.br = numeric(),         # T breadth (calculus)
  T.br.min = numeric(),     # T breadth (lower HDPI)
  T.br.max = numeric(),     # T breadth (upper HDPI)
  T.br.na = numeric(),      # T breadth na's

  a = numeric(),            # parameter: a
  b = numeric(),            # parameter: b
  tmax = numeric(),         # parameter: tmax
  d.t = numeric(),          # parameter: deltaT
  
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Study = character(),      # Study
  Parameter = character(),  # Model parameter (e.g. a, tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~6000)
  stringsAsFactors = FALSE             
)

# Specify Bayesian model structure
ni.fit <- 330000    # iterations / chain
nb.fit <- 30000     # burn in periods for each chain
nt.fit <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 3         # number of chains, total of 3,000 estimates for each model. 

parameters.lactin2 <- c("a", "b", "tmax", "d.t", "sigma", "r.pred") # repeated here

Temp.xs <- seq(-5, 42, 0.1) # Temperature gradient we're interested in - upped the granularity here
N.Temp.xs <-length(Temp.xs)

n <-0 # progression tracker

for (i in unique(df.t$id.number)){ # for each species ID
  
  n <- n + 1
  
  df.i <- df.t %>% 
    filter(id.number == i)
  df.i <- droplevels(df.i)
  
  lac_nls <- nls_multstart(mu ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = c(0.01, -2.5, 5, 1),
                           start_upper = c(0.19, -0.5, max(df.i$temp), 11),
                           lower = c(0, -3, min(df.i$temp), 0.1),
                           upper = c(0.5, 0, max(df.i$temp) + 5, 40),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  sum <- summary(lac_nls)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  lac.a <- df.nls[1,1] # Extract parameters
  lac.b <- df.nls[2,1]
  lac.tmax <- df.nls[3,1]
  lac.d.t <- df.nls[4,1]
  
  # Now we will use those terms to estimate priors for each species. 
  model_string <- sprintf("
  model {

    a       ~ dunif(%f, %f)
    tmax    ~ dunif(%f, %f)
    d.t     ~ dunif(%f, %f)
    b       ~ dunif(%f, %f)
    sigma   ~ dunif(0, 5)
    tau <- pow(sigma, -2)

    for (i in 1:N.obs) {
      trait.mu[i] <- exp(a * temp[i]) -
                     exp(a * tmax - (tmax - temp[i]) / d.t) +
                     b
      trait[i] ~ dnorm(trait.mu[i], tau)
    }

    for (j in 1:N.Temp.xs) {
      r.pred[j] <- exp(a * Temp.xs[j]) -
                   exp(a * tmax - (tmax - Temp.xs[j]) / d.t) +
                   b
    }
  }
  ", lac.a*0.9, lac.a*1.1, lac.tmax*0.9, lac.tmax*1.1, lac.d.t*0.9, lac.d.t*1.1, lac.b*1.1, lac.b*0.9)
  
  writeLines(model_string, con = "lactin_model.txt")
  
  # Custom initial values for the models. 
  inits.lactin <- function(lac.a, lac.tmax, lac.d.t, lac.b) {
    list(
      a = runif(1, lac.a*0.9, lac.a*1.1),  # More constrained initial values
      tmax = runif(1, lac.tmax*0.9, lac.tmax*1.1),
      d.t = runif(1, lac.d.t*0.9, lac.d.t*1.1),
      b = runif(1, lac.b*1.1, lac.b*0.9),
      sigma = runif(1, 0.1, 2)
    )
  }
  
  inits.fun <- function() inits.lactin(lac.a, lac.tmax, lac.d.t, lac.b)
  
  trait <- df.i$mu     # format the data for jags
  N.obs <- length(trait)
  temp <- df.i$temp
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  lac.jag <- jags(
    data = jag.data, 
    inits = inits.fun, 
    parameters.to.save = parameters.lactin2, 
    model.file = "lactin_model.txt",
    n.thin = nt.fit, 
    n.chains = nc.fit, 
    n.burnin = nb.fit, 
    n.iter = ni.fit, 
    DIC = TRUE, 
    working.directory = getwd()
  )
  
  jag.sum <- lac.jag$BUGSoutput$summary[c(1:4,476:477),]
  # mcmcplot(lac.jag)
  
  post <- as.data.frame(lac.jag$BUGSoutput$sims.list) # Get the posteriors
  
  jag.a <- median(post$a)    # Extract parameters
  jag.b <- median(post$b)
  jag.tmax <- median(post$tmax)
  jag.d.t <- median(post$d.t)
  
  # Topt
  calc_Topt <- function(a, b, tmax, d.t) {
    tryCatch(
      uniroot(
        function(temp) lactin2_deriv(temp, a, b, tmax, d.t),
        interval = c(-10, 45)
      )$root,
      error = function(e) NA
    )
  }
  
  T.opt <- mapply(
    calc_Topt,
    a = post$a,
    b = post$b,
    tmax = post$tmax,
    d.t = post$d.t
  )
  
  #rmax
  r.max <- mapply(
    lactin2,
    temp = T.opt,
    a = post$a,
    b = post$b,
    tmax = post$tmax,
    d.t = post$d.t
  )
  
  #Tmin
  Tmin.safe <- function(Topt, a, b, tmax, d.t, # Tmin equation
                        lower = -10) {
    
    f_low  <- lactin2(temp = lower, a = a, b = b, tmax = tmax, d.t = d.t)
    f_high <- lactin2(temp = Topt,  a = a, b = b, tmax = tmax, d.t = d.t)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2(temp, a, tmax, d.t, b),
      interval = c(lower, Topt)
    )$root
  }
  
  T.min <- mapply(
    Tmin.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t
  )
  
  #Tmax
  Tmax.safe <- function(Topt, a, b, tmax, d.t, # Tmax equation
                        upper = 50) {
    
    f_low  <- lactin2(temp = Topt, a = a, b = b, tmax = tmax, d.t = d.t)
    f_high <- lactin2(temp = upper,  a = a, b = b, tmax = tmax, d.t = d.t)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2(temp, a, tmax, d.t, b),
      interval = c(Topt, upper)
    )$root
  }
  
  T.max <- mapply(
    Tmax.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t
  )
  
  # Tbr
  
  r.half <- r.max/2
  
  Tmin.half.safe <- function(Topt, a, b, tmax, d.t, r_half, # Tmin equation
                        lower = -10) {
    
    f_low  <- lactin2_halfmax(temp = lower, a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    f_high <- lactin2_halfmax(temp = Topt,  a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2_halfmax(temp, a, b, tmax, d.t, r_half),
      interval = c(lower, Topt)
    )$root
  }
  
  T.min.half <- mapply(
    Tmin.half.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t,
    r_half = r.half
  )
  
  Tmax.half.safe <- function(Topt, a, b, tmax, d.t, r_half, # Tmax equation
                             upper = 50) {
    
    f_low  <- lactin2_halfmax(temp = Topt, a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    f_high <- lactin2_halfmax(temp = upper,  a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2_halfmax(temp, a, b, tmax, d.t, r_half),
      interval = c(Topt, upper)
    )$root
  }
  
  T.max.half <- mapply(
    Tmax.half.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t,
    r_half = r.half
  )
  
  T.br <- T.max.half - T.min.half
  
  thomas.summ.df <- rbind(thomas.summ.df, data.frame(                           # Add summary data
    Sp.id = df.i$id.number[1],                                                  # Species #
    Sp.name = df.i$Species.name[1],                                             # Species name
    Study = df.i$Study[1],                                                      # Study
    
    T.min = median(T.min, na.rm = T),                                           # Minimum T (calculus)
    T.min.min = hdi(T.min[!is.na(T.min)], credMass = 0.95)$CI_low,              # Minimum T (lower HDPI)
    T.min.max = hdi(T.min[!is.na(T.min)], credMass = 0.95)$CI_high,             # Minimum T (upper HDPI)
    T.min.na = mean(is.na(T.min)),                                              # % NA returns
    
    T.max = median(T.max, na.rm = T),                                           # Maximum T (calculus)
    T.max.min = hdi(T.max[!is.na(T.max)], credMass = 0.95)$CI_low,              # Maximum T (lower HDPI)
    T.max.max = hdi(T.max[!is.na(T.max)], credMass = 0.95)$CI_high,             # Maximum T (upper HDPI)
    T.max.na = mean(is.na(T.max)),                                              # % NA returns
    
    T.opt = median(T.opt, na.rm = T),                                           # Optimal T (calculus)
    T.opt.min = hdi(T.opt[!is.na(T.opt)], credMass = 0.95)$CI_low,              # Optimal T (lower HDPI)
    T.opt.max = hdi(T.opt[!is.na(T.opt)], credMass = 0.95)$CI_high,             # Optimal T (upper HDPI)
    
    r.max = median(r.max, na.rm = T),                                           # Maximum growth rate  (calculus)
    r.max.min = hdi(r.max[!is.na(r.max)], credMass = 0.95)$CI_low,              # Maximum growth rate  (lower HDPI)
    r.max.max = hdi(r.max[!is.na(r.max)], credMass = 0.95)$CI_high,             # Maximum growth rate  (upper HDPI)
    
    T.br = median(T.br, na.rm = T),                                             # T breadth (calculus)
    T.br.min = hdi(T.br[!is.na(T.br)], credMass = 0.95)$CI_low,                 # T breadth (lower HDPI)
    T.br.max = hdi(T.br[!is.na(T.br)], credMass = 0.95)$CI_high,                # T breadth (upper HDPI)
    T.br.na = mean(is.na(T.br)),                                                # % NA returns
    
    a = jag.a,                                                                  # parameter: a
    b = jag.b,                                                                  # parameter: b
    tmax = jag.tmax,                                                            # parameter: tmax
    d.t = jag.d.t                                                               # parameter: deltaT
  ))
  
  for (j in 1:6){
    fit.df <- rbind(fit.df, data.frame(                                         # Model performance data
      Sp.id = df.i$id.number[1],                                                # Species #
      Sp.name = df.i$Species.name[1],                                           # Species name
      Study = df.i$Study[1],                                                    # Study
      Parameter = rownames(jag.sum)[j],                                         # Model parameter (e.g. a, tmax, etc.)
      mean = jag.sum[j,1],                                                      # Posterior mean
      Rhat = jag.sum[j,8],                                                      # Rhat values
      n.eff = jag.sum[j,9],                                                     # Sample size estimates (should be ~3000)
      stringsAsFactors = FALSE            
    ))
  }
  
  print(paste("Done", n, "of ", length(unique(df.t$id.number))))
  
}



# Plotting & troubleshooting ----------------------------------------------


pred_lact <- function(temp, a, b, delta_t, tmax) {
  exp(a * temp) - exp(a * tmax - ((tmax - temp) / delta_t)) + b
}

curve.t <- tibble::tibble(
  res  = seq(0, 50, length.out = 200),
  rate = pred_lact(res, a = jag.a, b = jag.b, delta_t = jag.d.t, tmax = jag.tmax)
)

ggplot(curve.t, aes(x = res, y = rate)) +
  geom_line(size = 1.5, colour = "firebrick3") +
  geom_point(data = df.i,
             aes(x = temp, y = mu),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Temperature (°C)",
    y = "Exponential growth rate",
    title = "Thermal performance curve (Lactin II)"
  ) +
  theme_classic() +
  ylim(-0.1, 6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)
