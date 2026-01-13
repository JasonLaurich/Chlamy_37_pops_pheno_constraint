# Jason R Laurich

# December 29th, 2025

# We are going to reestimate Lewington-Pearce 2019 TPCs using nls.multstart

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(rTPC)
library(nls.multstart)
library(Deriv)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits

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

df <- read.csv("data-processed/12_Lewington-Pearce_2019_raw_data.csv") # Raw data file
head(df)
str(df)

df$Sp.fac <- as.factor(df$genus_species)

df$id <- paste0(df$block, ".", df$replicate)  # A unique identifier for each replicate.

df$log.fluorescence <- log(df$fluorescence + 0.001)

head(df)

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
      df.it.wl <- df.it.wl %>% 
        mutate(N0 = fluorescence[1])
      
      t.series <- unique(df.it.wl$day) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$day <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(log.fluorescence~day, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
      s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
      
      df.it.wl.th <- df.it.wl[df.it.wl$day <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
      # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
      
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
          temperature = df.it.wl$temperature[1],        
          id = df.it.wl$id[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          Sp.id = df.it.wl.th$Sp.fac[1],             # Species
          temperature = df.it.wl$temperature[1],     # Temperature
          id = df.it.wl$id[1],                       # ID
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
      }
      
    }
    
  }
  
}

write.csv(df.r.exp, "data-processed/307_Lewington-Pearce_2019_µ_estimates_temp_newer.csv") # let's save the file.

#### TPCS ####

# df.r.exp <- read.csv("data-processed/307_Lewington-Pearce_2019_µ_estimates_temp_newer.csv") # if needed

df.r.exp <- df.r.exp %>% 
  filter(Sp.id != "Chlamydomonas_reinhardtii")

lewington.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species ID
  Sp.name = character(),    # Species name
  T.min = numeric(),        # Minimum T (calculus)
  T.max = numeric(),        # Maximum T (calculus)
  T.opt = numeric(),        # Optimal T (calculus)
  r.max = numeric(),        # Maximum growth rate (calculus)
  T.br = numeric(),         # T breadth (calculus)
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  est = numeric(),          # Estimate
  se = numeric(),           # Standard error
  p = numeric(),            # p-value
  stringsAsFactors = FALSE                  
)

for (i in unique(df.r.exp$Sp.id)){ # For each species
  
  df.i <- df.r.exp %>% 
    filter(Sp.id == i, !is.na(r.exp)) %>% 
    arrange(temperature)
  df.i <- droplevels(df.i)
  
  lac_nls <- nls_multstart(r.exp ~ lactin2_1995(temp = temperature, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = c(0.01, -2.5, 5, 1),
                           start_upper = c(0.19, -0.5, max(df.i$temperature), 11),
                           lower = c(0, -3, min(df.i$temperature), 0.1),
                           upper = c(0.5, 0, max(df.i$temperature) + 5, 40),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  sum <- summary(lac_nls)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  cf.a <- df.nls[1,1] # Extract parameters
  cf.b <- df.nls[2,1]
  cf.tmax <- df.nls[3,1]
  cf.delta_t <- df.nls[4,1]
  
  # Find the T_opt: where the derivative crosses zero
  T_opt <- uniroot(
    function(temp) lactin2_deriv(temp, cf.a, cf.b, cf.tmax, cf.delta_t),
    interval = c(-10, 45)
  )$root
  
  r_max <- lactin2(temp=T_opt, cf.a=cf.a, cf.b=cf.b, cf.tmax=cf.tmax, cf.delta_t=cf.delta_t)
  
  Tmin <- uniroot(lactin2, interval = c(-1000000000, T_opt), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t)$root
  
  Tmax <- uniroot(lactin2, interval = c(T_opt,45), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t)$root
  
  r_half <- r_max/2 # calculate half of rmax and get the roots.
  
  Tlow <- uniroot(lactin2_halfmax, interval = c(Tmin, T_opt), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t, r_half = r_half)$root
  
  Thigh <- uniroot(lactin2_halfmax, interval = c(T_opt, Tmax), cf.a = cf.a, 
                   cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t, r_half = r_half)$root
  
  lewington.summ.df <- rbind(lewington.summ.df, data.frame( # Add summary data
    Sp.id = i,                                        # Species #
    Sp.name = df.i$Sp.id[1],                          # Species name
    T.min = Tmin,                                     # Minimum T (calculus)
    T.max = Tmax,                                     # Maximum T (calculus)
    T.opt = T_opt,                                    # Optimal T (calculus)
    r.max = r_max,                                    # Maximum growth rate (calculus)
    T.br = Thigh - Tlow                               # T breadth (calculus)
  ))
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(                         # Model performance data
      Sp.id = i,                                                # Species #
      Sp.name = df.i$Sp.id[1],                                  # Species name
      Parameter = df.nls$parameter[j],                          # Model parameter (e.g. a, b, tmax etc.)
      est = df.nls$Estimate[j],                                 # Estimate
      se = df.nls$`Std. Error`[j],                              # Error
      p = df.nls$`Pr(>|t|)`[j],                                 # p-values
      stringsAsFactors = FALSE            
    ))
  }
  
}

lewington.summ.df$Sp.id <- c(1:6)
fit.df$Sp.id <- c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4))

write.csv(lewington.summ.df, "data-processed/307a_Lewington-Pearce_2019_TPCs_newest.csv") # Save Lewington-Pearce 2019 TPC summary table
write.csv(fit.df, "data-processed/307b_Lewington-Pearce_2019_TPCs_fits_newest.csv") # Save model fit summary table

# Calculate µ across Light and fit Monods -------------------------------------

df.l <- df %>% 
  filter(nitrate_level == 1000)

light <- as.vector(unique(df.l$light_level))# for looping through light levels
ord.light<- sort(light)

temp <- as.vector(unique(df.l$temperature))
ord.temp <-sort(temp)

df.r.exp.l <- data.frame(              # Summary dataframe for r_exp
  Sp.id = character(),                 # Species
  light = numeric(),                   # Light level
  temp = numeric(),                    # Temperature
  id = character(),                    # Unique ID
  r.exp = numeric()                    # Thresholded r.exp
)

for (i in unique(df.l$Sp.fac)){ # for every species
  
  for (t in ord.light){ # at every light level
    
    for (x in ord.temp){ # at every temperature
      
      df.i <- df.l %>%
        filter(Sp.fac == i,
               light_level == t,
               temperature == x)
        
      df.it <- droplevels(df.i) # Drop unused levels to isolate well replicate IDs at given t 
      
      for (w in unique(df.it$id)){ # Doing this separately for each replicate
        
        df.it.wl <- subset(df.it, as.numeric(df.it$id) == w)
        
        df.it.wl <- df.it.wl[order(df.it.wl$day), ]
        df.it.wl <- df.it.wl %>% 
          mutate(N0 = fluorescence[1])
        
        t.series <- unique(df.it.wl$day) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
        t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
        
        ln.slopes <- c() # Re-initialize this too!
        
        for (z in t.series){
          
          df.it.wl.sl <- df.it.wl[df.it.wl$day <= z, ] # Subset the data to exclude time points above our window
          
          ln_slope <- lm(log.fluorescence~day, data = df.it.wl.sl)
          
          ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
          
        } # So now we have our slopes for each well.ID x Pop x P level
        
        s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
        
        df.it.wl.th <- df.it.wl[df.it.wl$day <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
        # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
        
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
            light = df.it.wl.th$light_level[1],  
            temp = df.it.wl.th$temperature[1], 
            id = df.it.wl$id[1],                
            r.exp = NA        
          ))
          
        }else{
          # Add data to our summary table
          df.r.exp.l <- rbind(df.r.exp.l, data.frame(
            Sp.id = df.it.wl.th$Sp.fac[1],             # Species
            light = df.it.wl.th$light_level[1],        # Light
            temp = df.it.wl.th$temperature[1],         # Temperature
            id = df.it.wl$id[1],                       # ID
            r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
          ))
        }
       }
      
    }
    
  }
  
}
    
write.csv(df.r.exp.l, "data-processed/307g_Lewington-Pearce_2019_µ_estimates_light_new.csv") # let's save the file.

# Light Monod curves ------------------------------------------------------

# df.r.exp.l <- read.csv("data-processed/307g_Lewington-Pearce_2019_µ_estimates_light_new.csv") # if needed

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

# Let's do larger models for the final things (10 times larger)
ni.fit <- 330000    # iterations / chain
nb.fit <- 30000     # burn in periods for each chain
nt.fit <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 3         # number of chains, total of 3,000 estimates for each model. 

lewington.summ.df.l <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),                 # Species ID
  Sp.name = character(),             # Species name
  K.s = numeric(),                   # Half-saturation constant
  r.max = numeric(),                 # Maximum population growth rate
  R.jag = numeric(),                 # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),                 # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
  stringsAsFactors = FALSE           # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~3000)
  stringsAsFactors = FALSE            
)

for (i in unique(df.r.exp.l$Sp.id)){ # For each species
  
  df.i <- df.r.exp.l %>% 
    filter(Sp.id == i, !is.na(r.exp)) %>% 
    arrange(light)
  df.i <- droplevels(df.i)
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  
  light <- df.i$light
  
  S.pred <- seq(0, 150, 0.5) # Light gradient we're interested in - upped the granularity here
  N.S.pred <-length(S.pred) # We'll reset this internally since the gradient varies substantially
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = light, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the light Monod function. 
    data = jag.data,
    inits = inits.monod,
    parameters.to.save = parameters.monod,
    model.file = "monod.light.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  print(paste("Done", i, "of 7"))
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3, (150/0.5) + 5),]   # generate the sequence of r.pred values
  df.jags$light <- seq(0, 150, 0.5)
  
  lewington.summ.df.l <- rbind(lewington.summ.df.l, data.frame( # Add summary data  
    Sp.id = i,                                                                                  # Species #
    Sp.name = df.i$Sp.id[1],                                                                    # Species name
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                    # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                  # Maximum population growth rate
    R.jag = df.jags$light[which(df.jags$mean > 0.1)[1]],                                        # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.1*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.1)     # Minimum resource requirement for positive growth (from math)
  ))
  
  light_sum <- monod_jag$BUGSoutput$summary[c(1:3, (150/0.5) + 5),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = i,                                 # Species #
      Sp.name = df.i$Sp.id[1],                   # Species name      
      Parameter = rownames(light_sum)[j],        # Model parameter (e.g. K_s, r_max, etc.)
      mean = light_sum[j,1],                     # Posterior mean
      Rhat = light_sum[j,8],                     # Rhat values
      n.eff = light_sum[j,9]                     # Sample size estimates (should be ~3000)
    ))
    
  }
  
}

lewington.summ.df.l$Sp.id <- c(1:7)
fit.df$Sp.id <- c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7,4))

write.csv(lewington.summ.df.l, "data-processed/307h_Lewington-Pearce_2019_Monod_light_new.csv") # Save Lewington-Pearce 2019 Light Monod summary table
write.csv(fit.df, "data-processed/307i_Lewington2019_Monod_light_fits_new.csv") # Save model fit summary table

# Calculate µ across nitrogen and fit Monods -------------------------------------

df.n <- df %>% 
  filter(light_level >= 105.5) # For some reason, light levels are not the same at different N levels...

nit <- as.vector(unique(df.n$nitrate_level))# for looping through nitrogen levels
ord.nit<- sort(nit)

temp <- as.vector(unique(df.n$temperature))
ord.temp <-sort(temp)

df.r.exp.n <- data.frame(              # Summary dataframe for r_exp
  Sp.id = character(),                 # Species
  nit = numeric(),                     # Nitrogen level
  temp = numeric(),                    # Temperature
  id = character(),                    # Unique ID
  r.exp = numeric()                    # Thresholded r.exp
)

for (i in unique(df.n$Sp.fac)){ # for every species
  
  for (t in ord.nit){ # at every nit level
    
    for (x in ord.temp){
      
      df.i <- df.n %>%
        filter(Sp.fac == i,
               nitrate_level == t,
               temperature == x)
      
      df.it <- droplevels(df.i) # Drop unused levels to isolate well replicate IDs at given t 
      
      for (w in unique(df.it$id)){ # Doing this separately for each replicate
        
        df.it.wl <- subset(df.it, as.numeric(df.it$id) == w)
        
        df.it.wl <- df.it.wl[order(df.it.wl$day), ]
        df.it.wl <- df.it.wl %>% 
          mutate(N0 = fluorescence[1])
        
        t.series <- unique(df.it.wl$day) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
        t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
        
        ln.slopes <- c() # Re-initialize this too!
        
        for (z in t.series){
          
          df.it.wl.sl <- df.it.wl[df.it.wl$day <= z, ] # Subset the data to exclude time points above our window
          
          ln_slope <- lm(log.fluorescence~day, data = df.it.wl.sl)
          
          ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
          
        } # So now we have our slopes for each well.ID x Pop x P level
        
        s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
        
        df.it.wl.th <- df.it.wl[df.it.wl$day <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
        # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
        
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
            nit = df.it.wl.th$nitrate_level[1],  
            temp = df.it.wl.th$temperature[1],    
            id = df.it.wl$id[1],                
            r.exp = NA        
          ))
          
        }else{
          # Add data to our summary table
          df.r.exp.n <- rbind(df.r.exp.n, data.frame(
            Sp.id = df.it.wl.th$Sp.fac[1],             # Species
            nit = df.it.wl.th$nitrate_level[1],        # Nitrogen
            temp = df.it.wl.th$temperature[1],         # Temperature
            id = df.it.wl$id[1],                       # ID
            r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
          ))
        }
      
      }
      
    }
    
  }
  
}

write.csv(df.r.exp.n, "data-processed/307c_Lewington-Pearce_2019_µ_estimates_nitrogen_new.csv") # let's save the file.

# Nitrogen Monod curves ------------------------------------------------------

# df.r.exp.n <- read.csv("data-processed/307c_Lewington-Pearce_2019_µ_estimates_nitrogen_new.csv") # if needed

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

# Let's do larger models for the final things (10 times larger)
ni.fit <- 330000    # iterations / chain
nb.fit <- 30000     # burn in periods for each chain
nt.fit <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 3         # number of chains, total of 3,000 estimates for each model. 

lewington.summ.df.n <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),                 # Species ID
  Sp.name = character(),             # Species name
  K.s = numeric(),                   # Half-saturation constant
  r.max = numeric(),                 # Maximum population growth rate
  R.jag = numeric(),                 # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),                 # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
  stringsAsFactors = FALSE           # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~3000)
  stringsAsFactors = FALSE            
)

for (i in unique(df.r.exp.n$Sp.id)){ # For each species
  
  df.i <- df.r.exp.n %>% 
    filter(Sp.id == i, !is.na(r.exp)) %>% 
    arrange(nit)
  df.i <- droplevels(df.i)
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  
  nit <- df.i$nit
  
  S.pred <- seq(0, 1000, 0.5) # Nitrogen gradient we're interested in - upped the granularity here
  N.S.pred <-length(S.pred) # We'll reset this internally since the gradient varies substantially
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = nit, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the nitrogen Monod function. 
    data = jag.data,
    inits = inits.monod,
    parameters.to.save = parameters.monod,
    model.file = "monod.nit.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  print(paste("Done", i, "of 7"))
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3, (1000/0.5) + 5),]   # generate the sequence of r.pred values
  df.jags$light <- seq(0, 1000, 0.5)
  
  lewington.summ.df.n <- rbind(lewington.summ.df.n, data.frame( # Add summary data  
    Sp.id = i,                                                                                  # Species #
    Sp.name = df.i$Sp.id[1],                                                                    # Species name
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                    # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                  # Maximum population growth rate
    R.jag = df.jags$light[which(df.jags$mean > 0.1)[1]],                                        # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.1*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.1)     # Minimum resource requirement for positive growth (from math)
  ))
  
  nit_sum <- monod_jag$BUGSoutput$summary[c(1:3, (1000/0.5) + 5),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = i,                                 # Species #
      Sp.name = df.i$Sp.id[1],                   # Species name      
      Parameter = rownames(nit_sum)[j],          # Model parameter (e.g. K_s, r_max, etc.)
      mean = nit_sum[j,1],                       # Posterior mean
      Rhat = nit_sum[j,8],                       # Rhat values
      n.eff = nit_sum[j,9]                       # Sample size estimates (should be ~3000)
    ))
    
  }
  
}

lewington.summ.df.n$Sp.id <- c(1:7)
fit.df$Sp.id <- c(rep(1,4), rep(2,4), rep(3,4), rep(4,4), rep(5,4), rep(6,4), rep(7,4))

write.csv(lewington.summ.df.n, "data-processed/307d_Lewington-Pearce_2019_Monod_nit_new.csv") # Save Lewington-Pearce 2019 Nitrogen Monod summary table
write.csv(fit.df, "data-processed/307e_Lewington2019_Monod_nit_fits_new.csv") # Save model fit summary table

# Calculate µ (P) ---------------------------------------------------------


