# Jason R Laurich
# April 12, 2025

# Going to work with the Lewington-Pearce et al 2019 (https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/oik.06060) data to fit
# Monod curves for nitrogen and light and TPCs


# Load packages -----------------------------------------------------------

library(tidyr)
library(cowplot)
library(tidyverse)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)
library(rTPC)
library(nls.multstart)

# Load and examine the data -----------------------------------------------

df <- read.csv("data-processed/19_Lewington-Pearce2019_raw_data.csv") # Raw data file
head(df)
str(df)

df$Sp.fac <- as.factor(df$genus_species)

df$id <- paste0(df$block, ".", df$replicate)  # A unique identifier for each replicate.

df$log.fluorescence <- log(df$fluorescence + 0.001)

# What concentration of nitrogen should we be fitting TPCs at? Let's look at RFUs at D5

df %>% 
  group_by(Sp.fac, nitrate_level) %>% 
  filter(day == 5) %>% 
  summarize(mean.fluor = mean(fluorescence, na.rm = T)) %>% 
  print(n = 73)

# For most species 1000 is best, but not C. reinhardtii. However, that is the concentration 
# we have been using so we will stick with it.

# What about light?

df %>% 
  group_by(Sp.fac, light_level) %>% 
  filter(day == 5) %>% 
  summarize(mean.fluor = mean(fluorescence, na.rm = T)) %>% 
  print(n = 95)

# We're goint to go with 106 (105.7) — this is the highest for most sp.
# not including C. reinhardtii which has some fluorescences at light levels not 
# specified in the paper that are insanely high. 

# Calculate µ across Temperature and fit TPCs -------------------------------------

df.t <- df %>% 
  filter(light_level == 105.5, nitrate_level == 1000)

temp <- as.vector(unique(df.t$temperature))# for looping through nitrate levels
ord.temp<- sort(temp)

df.r.exp <- data.frame(                # Summary dataframe for r_exp
  Sp.id = character(),                 # Species
  temperature = numeric(),             # Temperature
  id = character(),                    # Unique ID
  r.exp = numeric()                    # Thresholded r.exp
)

for (i in unique(df.t$Sp.fac)){ # for every species
  
  for (t in ord.temp){ # at every temperature
    
    df.i <- df.t[df.t$Sp.fac == i, ]
    df.it <- df.i[df.i$temperature == t, ]
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in unique(df.it$id)){ # Doing this separately for each replicate
      
      df.it.wl <- subset(df.it, as.numeric(df.it$id) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$day), ]
      
      df.it.wl2 <- df.it.wl %>% 
        mutate(N0 = fluorescence[1])
      
      if (df.it.wl2$fluorescence[2] <  df.it.wl2$fluorescence[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.wl2 <- df.it.wl2[-1,]
        df.it.wl2$N0 <- df.it.wl2$fluorescence[1]
      }
      
      t.series <- unique(df.it.wl2$day) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl2[df.it.wl2$day <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(log.fluorescence~day, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each id x sp x T level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 10%
        
        if (ln.slopes[r] != 0 & round(ln.slopes[r], digits = 5) !=0) { # We also need to account for tiny values that are basically 0 (e.g. 5 e-16, but are messing up our loops)
          # percent.chg <- (ln.slopes[s] - ln.slopes[s + 1]) / ln.slopes[s] This was the reason we were getting weird negative results! If there was a stochastic drop.
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] # This should (along with the next line) fix it.
          
          if (percent.chg >= 1.10 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { # Because now, the drop-off ignores transiently negative slopes. 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
      
      df.it.wl.th <- df.it.wl2[df.it.wl2$day <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(fluorescence ~ N0 * exp(r*day),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp <- rbind(df.r.exp, data.frame(
          Sp.id = df.it.wl.th$Sp.fac[1],          
          temperature = df.it.wl2$temperature[1],        
          id = df.it.wl2$id[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          Sp.id = df.it.wl.th$Sp.fac[1],             # Species
          temperature = df.it.wl2$temperature[1],     # Temperature
          id = df.it.wl2$id[1],                       # ID
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
      }
      
    }
    
  }

}

write.csv(df.r.exp, "data-processed/19a_L-P_Temp_r_estimates.csv") # let's save the file.

#### TPCS ####

# df.r.exp <- read.csv("data-processed/19a_L-P_Temp_r_estimates.csv") # if needed

lewington.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),               # Species ID
  DIC = numeric(),                 # DIC
  T.min.raw = numeric(),           # Minimum T (Jags raw)
  T.max.raw = numeric(),           # Maximum T (Jags raw)
  T.opt.raw = numeric(),           # Optimal T (Jags raw)
  r.max.raw = numeric(),           # Maximum growth rate (Jags raw)
  T.br.raw = numeric(),            # T breadth (Jags raw)
  stringsAsFactors = FALSE         # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~3000)
  stringsAsFactors = FALSE            
)

# Let's do larger models for the final things (10 times larger)
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

Temp.xs <- seq(-5, 45, 0.1) # Temperature gradient we're interested in - upped the granularity here
N.Temp.xs <-length(Temp.xs)

for (i in unique(df.r.exp$Sp.id)){ # For each species
  
  df.i <- df.r.exp %>% 
    filter(Sp.id == i, !is.na(r.exp))
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  
  temp <- df.i$temperature
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  start.vals.lac <- get_start_vals(df.i$temperature, df.i$r.exp, model_name = 'lactin2_1995')
  
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
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)[-c(1:6),]   # generate the sequence of r.pred values
  df.jags$temp <- seq(-5, 45, 0.1)
  
  lewington.summ.df <- rbind(lewington.summ.df, data.frame(                            # Add summary data
    Sp.id = i,                                                                         # Species name 
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
}

# It does actually seem like there is only data for C. reinhardtii at 20 C! 

write.csv(lewington.summ.df, "data-processed/19b_Lewington2019_TPCs.csv") # Save Lewington-Pearce 2019 TPC summary table
write.csv(fit.df, "data-processed/19c_Lewington2019_TPCs_fits.csv") # Save model fit summary table

# Calculate µ across Nitrogen and fit Monod curves -------------------------------

# So we are going to do this at 20C to stay consistent with the rest of the dataset.
# As it turns out, light seems to have been held constant at 140.6 for the nitrate experiments

df.n <- df %>% 
  filter(light_level == 140.6, temperature == 20)

nit <- as.vector(unique(df.n$nitrate_level))# for looping through nitrate levels
ord.nit<- sort(nit)

df.r.exp.n <- data.frame(              # Summary dataframe for r_exp
  Sp.id = character(),                 # Species
  nit = numeric(),                     # Nitrate level
  id = character(),                    # Unique ID
  r.exp = numeric()                    # Thresholded r.exp
)

for (i in unique(df.n$Sp.fac)){ # for every species
  
  for (t in ord.nit){ # at every nitrogen level
    
    df.i <- df.n[df.n$Sp.fac == i, ]
    df.it <- df.i[df.i$nitrate_level == t, ]
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in unique(df.it$id)){ # Doing this separately for each replicate
      
      df.it.wl <- subset(df.it, as.numeric(df.it$id) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$day), ]
      
      df.it.wl2 <- df.it.wl %>% 
        mutate(N0 = fluorescence[1])
      
      if (df.it.wl2$fluorescence[2] <  df.it.wl2$fluorescence[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.wl2 <- df.it.wl2[-1,]
        df.it.wl2$N0 <- df.it.wl2$fluorescence[1]
      }
      
      t.series <- unique(df.it.wl2$day) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl2[df.it.wl2$day <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(log.fluorescence~day, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each id x sp x N level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 10%
        
        if (ln.slopes[r] != 0 & round(ln.slopes[r], digits = 5) !=0) { # We also need to account for tiny values that are basically 0 (e.g. 5 e-16, but are messing up our loops)
          # percent.chg <- (ln.slopes[s] - ln.slopes[s + 1]) / ln.slopes[s] This was the reason we were getting weird negative results! If there was a stochastic drop.
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] # This should (along with the next line) fix it.
          
          if (percent.chg >= 1.10 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { # Because now, the drop-off ignores transiently negative slopes. 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
      
      df.it.wl.th <- df.it.wl2[df.it.wl2$day <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(fluorescence ~ N0 * exp(r*day),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp.n <- rbind(df.r.exp.n, data.frame(
          Sp.id = df.it.wl.th$Sp.fac[1],          
          nit = df.it.wl2$nitrate_level[1],        
          id = df.it.wl2$id[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp.n <- rbind(df.r.exp.n, data.frame(
          Sp.id = df.it.wl.th$Sp.fac[1],             # Species
          nit = df.it.wl2$nitrate_level[1],          # Nitrogen
          id = df.it.wl2$id[1],                      # ID
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
      }
      
    }
    
  }
  
}

write.csv(df.r.exp.n, "data-processed/19d_Lewington2019_Nit_r_estimates.csv") # let's save the file.

#### Monod Curves ####

# df.r.exp.n <- read.csv("data-processed/19d_Lewington2019_Nit_r_estimates.csv") # if needed

lewington.summ.df.n <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),                 # Species ID
  DIC = numeric(),                   # DIC
  K.s = numeric(),                   # Half-saturation constant
  r.max = numeric(),                 # Maximum population growth rate
  R.jag = numeric(),                 # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),                 # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
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

# Put these here too, in case we want to run this section separately
ni.fit <- 330000    # iterations / chain
nb.fit <- 30000     # burn in periods for each chain
nt.fit <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 3         # number of chains, total of 3,000 estimates for each model. 

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

S.pred <- seq(0, 1000, 0.5) # Nitrogen gradient we are interested in here (concentration)
N.S.pred <-length(S.pred)   # We will keep this to also have 500 levels in the gradient?
# This won't really work - the problem with nitrogen is the combination of wide range of []s and
# fine-grained resolution at the lower end. We're going to estimate r every 0.5 so 2000 levels. 

for (i in unique(df.r.exp.n$Sp.id)){ # For each species
  
  df.i <- df.r.exp.n %>% 
    filter(Sp.id == i, !is.na(r.exp))
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  
  nit <- df.i$nit
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = nit, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the phosphorous Monod function. 
    data = jag.data,
    inits = inits.monod,
    parameters.to.save = parameters.monod,
    model.file = "monod.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3,2005),]   # generate the sequence of r.pred values
  df.jags$nit <- seq(0, 1000, 0.5)
  
  lewington.summ.df.n <- rbind(lewington.summ.df.n, data.frame(                                 # Add summary data
    Sp.id = i,                                                                                  # Species name 
    DIC = monod_jag$BUGSoutput$DIC,                                                             # DIC
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                    # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                  # Maximum population growth rate
    R.jag = df.jags$nit[which(df.jags$mean > 0.56)[1]],                                        # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.56*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.56)   # Minimum resource requirement for positive growth (from math)
  ))
  
  nit_sum <- monod_jag$BUGSoutput$summary[c(1:3, 2005),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(        # Model performance data
      Sp.id = i,                               # Species       
      Parameter = rownames(nit_sum)[j],        # Model parameter (e.g. K_s, r_max, etc.)
      mean = nit_sum[j,1],                     # Posterior mean
      Rhat = nit_sum[j,8],                     # Rhat values
      n.eff = nit_sum[j,9]                     # Sample size estimates (should be ~3000)
    ))
  
  }
}

# These numbers seem suspect. Need to investigate further
df.r.exp.n
lewington.summ.df.n

# There are no 0 N levels for most species, which means the growths are positive for many 
# even at the lowest N levels. 

write.csv(lewington.summ.df.n, "data-processed/19e_Lewington2019_Nit_Monods.csv") # Save Lewington-Pearce 2019 N Monod summary table
write.csv(fit.df, "data-processed/19f_Lewington2019_Nit_Monod_fits.csv") # Save model fit summary table

# Calculate µ across Light and fit Monod Curves ---------------------------

# So we are going to do this at 20C to stay consistent with the rest of the dataset.
# As it turns out, light seems to have been held constant at 140.6 for the nitrate experiments

df.l <- df %>% 
  filter(nitrate_level == 1000, temperature == 20)

light <- as.vector(unique(df.l$light_level))# for looping through light levels
ord.light<- sort(light)

df.r.exp.l <- data.frame(              # Summary dataframe for r_exp
  Sp.id = character(),                 # Species
  light = numeric(),                   # Light level
  id = character(),                    # Unique ID
  r.exp = numeric()                    # Thresholded r.exp
)

for (i in unique(df.l$Sp.fac)){ # for every species
  
  for (t in ord.light){ # at every nitrogen level
    
    df.i <- df.l[df.l$Sp.fac == i, ]
    df.it <- df.i[df.i$light_level == t, ]
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in unique(df.it$id)){ # Doing this separately for each replicate
      
      df.it.wl <- subset(df.it, as.numeric(df.it$id) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$day), ]
      
      df.it.wl2 <- df.it.wl %>% 
        mutate(N0 = fluorescence[1])
      
      if (df.it.wl2$fluorescence[2] <  df.it.wl2$fluorescence[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.wl2 <- df.it.wl2[-1,]
        df.it.wl2$N0 <- df.it.wl2$fluorescence[1]
      }
      
      t.series <- unique(df.it.wl2$day) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl2[df.it.wl2$day <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(log.fluorescence~day, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each id x sp x N level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 10%
        
        if (ln.slopes[r] != 0 & round(ln.slopes[r], digits = 5) !=0) { # We also need to account for tiny values that are basically 0 (e.g. 5 e-16, but are messing up our loops)
          # percent.chg <- (ln.slopes[s] - ln.slopes[s + 1]) / ln.slopes[s] This was the reason we were getting weird negative results! If there was a stochastic drop.
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] # This should (along with the next line) fix it.
          
          if (percent.chg >= 1.10 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { # Because now, the drop-off ignores transiently negative slopes. 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
      
      df.it.wl.th <- df.it.wl2[df.it.wl2$day <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(fluorescence ~ N0 * exp(r*day),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp.l <- rbind(df.r.exp.l, data.frame(
          Sp.id = df.it.wl.th$Sp.fac[1],          
          light = df.it.wl2$light_level[1],        
          id = df.it.wl2$id[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp.l <- rbind(df.r.exp.l, data.frame(
          Sp.id = df.it.wl.th$Sp.fac[1],             # Species
          light = df.it.wl2$light_level[1],          # Light
          id = df.it.wl2$id[1],                      # ID
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
      }
      
    }
    
  }
  
}

write.csv(df.r.exp.l, "data-processed/19g_Lewington2019_Light_r_estimates.csv") # let's save the file.

#### Monod Curves ####

# df.r.exp.l <- read.csv("data-processed/19g_Lewington2019_Light_r_estimates.csv") # if needed

lewington.summ.df.l <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),                 # Species ID
  DIC = numeric(),                   # DIC
  K.s = numeric(),                   # Half-saturation constant
  r.max = numeric(),                 # Maximum population growth rate
  R.jag = numeric(),                 # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),                 # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
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

# Put these here too, in case we want to run this section separately
ni.fit <- 330000    # iterations / chain
nb.fit <- 30000     # burn in periods for each chain
nt.fit <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 3         # number of chains, total of 3,000 estimates for each model. 

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

S.pred <- seq(0, 150, 0.3) # Light gradient we are interested in here
N.S.pred <-length(S.pred)   # We will keep this to also have 500 levels in the gradient?

for (i in unique(df.r.exp.l$Sp.id)){ # For each species
  
  df.i <- df.r.exp.l %>% 
    filter(Sp.id == i, !is.na(r.exp))
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  
  light <- df.i$light
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = light, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the light Monod function. 
    data = jag.data,
    inits = inits.monod,
    parameters.to.save = parameters.monod,
    model.file = "monod.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3,505),]   # generate the sequence of r.pred values
  df.jags$light <- seq(0, 150, 0.3)
  
  lewington.summ.df.l <- rbind(lewington.summ.df.l, data.frame(                                 # Add summary data
    Sp.id = i,                                                                                  # Species name 
    DIC = monod_jag$BUGSoutput$DIC,                                                             # DIC
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                    # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                  # Maximum population growth rate
    R.jag = df.jags$light[which(df.jags$mean > 0.56)[1]],                                       # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.56*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.56)   # Minimum resource requirement for positive growth (from math)
  ))
  
  light_sum <- monod_jag$BUGSoutput$summary[c(1:3, 505),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = i,                                 # Species       
      Parameter = rownames(light_sum)[j],        # Model parameter (e.g. K_s, r_max, etc.)
      mean = light_sum[j,1],                     # Posterior mean
      Rhat = light_sum[j,8],                     # Rhat values
      n.eff = light_sum[j,9]                     # Sample size estimates (should be ~3000)
    ))
    
  }
}

# These numbers seem suspect. Need to investigate further
df.r.exp.l
lewington.summ.df.l

# There are no 0 N levels for most species, which means the growths are positive for many 
# even at the lowest N levels. 

write.csv(lewington.summ.df.l, "data-processed/19h_Lewington2019_Light_Monods.csv") # Save Lewington-Pearce 2019 Light Monod summary table
write.csv(fit.df, "data-processed/19i_Lewington2019_Light_Monod_fits.csv") # Save model fit summary table

