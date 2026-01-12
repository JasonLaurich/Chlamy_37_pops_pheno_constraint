# Jason R Laurich

# December 31st, 2025

# We are going to fit TPCs and Monod curves here. For temp, we'll compile data from all light, nitrogen, and phosphorous experiments.

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(rTPC)
library(nls.multstart)
library(Deriv)

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

df.l <- read.csv("data-processed/400a_Levasseur2025_l_rawdata.csv") # Raw data file
head(df.l)
str(df.l)

df.l$sp.num <- as.integer(factor(df.l$species_updated))
df.l$log.RFU <- log(df.l$RFU + 0.001)

df.lt <- df.l %>% 
  filter(value_resource == max(value_resource, na.rm = TRUE))

df.n <- read.csv("data-processed/400b_Levasseur2025_n_rawdata.csv") # Raw data file
head(df.n)
str(df.n)

df.n$sp.num <- as.integer(factor(df.n$species_updated))
df.n$log.RFU <- log(df.n$RFU + 0.001)

df.nt <- df.n %>% 
  filter(value_resource == max(value_resource, na.rm = TRUE))

df.p <- read.csv("data-processed/400c_Levasseur2025_p_rawdata.csv") # Raw data file
head(df.p)
str(df.p)

df.p$sp.num <- as.integer(factor(df.p$species_updated))
df.p$log.RFU <- log(df.p$RFU + 0.001)

df.pt <- df.p %>% 
  filter(value_resource == max(value_resource, na.rm = TRUE))

df.t <- rbind(df.lt, df.nt, df.pt)

df.t <- df.t %>% # We'll concatenate well id and resource level just in case the same well appears in different experiments
  mutate(unique.id = paste(well_id, level_of_resource, sep = "_")) # So that each entry is treated seperately!

df.l <- df.l %>% # We'll concatenate well id and resource level just in case the same well appears in different experiments
  mutate(unique.id = paste(well_id, level_of_resource, sep = "_")) # So that each entry is treated seperately!

df.n <- df.n %>% 
  filter(temperature == 20) # Filter out non 20C data. 

df.n <- df.n %>% # We'll concatenate well id and resource level just in case the same well appears in different experiments
  mutate(unique.id = paste(well_id, level_of_resource, sep = "_")) # So that each entry is treated seperately!

# Estimate µ --------------------------------------------------------------

temp <- as.vector(unique(df.t$temperature))# for looping through temps
ord.temp<- sort(temp)

df.r.exp <- data.frame(                # Summary dataframe for r_exp
  Sp.id = character(),                 # Species
  temperature = numeric(),             # Temperature
  id = character(),                    # Unique ID
  r.exp = numeric()                    # Thresholded r.exp
)

for (i in 1:20){ # for every species
  
  for (t in ord.temp){ # at every temperature
    
    df.i <- df.t[df.t$sp.num == i, ]
    df.it <- df.i[df.i$temperature == t, ]
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in unique(df.it$unique.id)){ # Doing this separately for each replicate
      
      df.it.wl <- subset(df.it, df.it$unique.id == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$Time_days), ]
      df.it.wl <- df.it.wl %>% 
        mutate(N0 = RFU[1])
      
      t.series <- unique(df.it.wl$Time_days) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$Time_days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(log.RFU~Time_days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
      s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
      
      df.it.wl.th <- df.it.wl[df.it.wl$Time_days <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
      # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*Time_days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp <- rbind(df.r.exp, data.frame(
          Sp.id = df.it.wl.th$species_updated[1],           
          temperature = df.it.wl$temperature[1],        
          id = df.it.wl$sp.num[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          Sp.id = df.it.wl.th$species_updated[1],    # Species
          temperature = df.it.wl$temperature[1],     # Temperature
          id = df.it.wl$sp.num[1],                   # ID
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
      }
      
    }
    
  }
  
}

write.csv(df.r.exp, "data-processed/401_Levasseur_2025_µ_estimates_temp.csv") # let's save the file.

# TPCs --------------------------------------------------------------------

# df.r.exp <- read.csv("data-processed/401_Levasseur_2025_µ_estimates_temp.csv") # if needed

levasseur.t.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
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

for (i in 1:20){ # For each species
  
  df.i <- df.r.exp %>% 
    filter(id == i, !is.na(r.exp)) %>% 
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
  
  levasseur.t.summ.df <- rbind(levasseur.t.summ.df, data.frame( # Add summary data
    Sp.id = df.i$id[1],                               # Species #
    Sp.name = df.i$Sp.id[1],                          # Species name
    T.min = Tmin,                                     # Minimum T (calculus)
    T.max = Tmax,                                     # Maximum T (calculus)
    T.opt = T_opt,                                    # Optimal T (calculus)
    r.max = r_max,                                    # Maximum growth rate (calculus)
    T.br = Thigh - Tlow                               # T breadth (calculus)
  ))
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(                         # Model performance data
      Sp.id = df.i$id[1],                                       # Species #
      Sp.name = df.i$Sp.id[1],                                  # Species name
      Parameter = df.nls$parameter[j],                          # Model parameter (e.g. a, b, tmax etc.)
      est = df.nls$Estimate[j],                                 # Estimate
      se = df.nls$`Std. Error`[j],                              # Error
      p = df.nls$`Pr(>|t|)`[j],                                 # p-values
      stringsAsFactors = FALSE            
    ))
  }
  
}

write.csv(levasseur.t.summ.df, "data-processed/402a_Levasseur_2025_TPCs_newer.csv") # Save Levasseur 2025 TPC summary table
write.csv(fit.df, "data-processed/402b_Levasseur_2025_TPCs_fits_newest.csv") # Save model fit summary table

# Estimate µ for light ----------------------------------------------------

light <- as.vector(unique(df.l$value_resource))# for looping through lights
ord.light<- sort(light)

temp <- as.vector(unique(df.l$temperature))
ord.temp <- sort(temp)

df.r.exp.l <- data.frame(              # Summary dataframe for r_exp
  Sp.id = character(),                 # Species
  light = numeric(),                   # Light
  temp = numeric(),                    # Temperature
  id = character(),                    # Unique ID
  r.exp = numeric()                    # Thresholded r.exp
)

for (i in 1:20){ # for every species
  
  for (t in ord.light){ # at every light level
    
    for (x in ord.temp){
      
      df.i <- df.l %>% 
        filter(sp.num == i,
               value_resource == t,
               temperature == x)
      
      df.it <- droplevels(df.i) # Drop unused levels to isolate well replicate IDs at given t
      
      for (w in unique(df.it$unique.id)){ # Doing this separately for each replicate
        
        df.it.wl <- subset(df.it, df.it$unique.id == w)
        
        df.it.wl <- df.it.wl[order(df.it.wl$Time_days), ]
        df.it.wl <- df.it.wl %>% 
          mutate(N0 = RFU[1])
        
        t.series <- unique(df.it.wl$Time_days) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
        t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
        
        ln.slopes <- c() # Re-initialize this too!
        
        for (z in t.series){
          
          df.it.wl.sl <- df.it.wl[df.it.wl$Time_days <= z, ] # Subset the data to exclude time points above our window
          
          ln_slope <- lm(log.RFU~Time_days, data = df.it.wl.sl)
          
          ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
          
        } # So now we have our slopes for each well.ID x Pop x P level
        
        s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
        
        df.it.wl.th <- df.it.wl[df.it.wl$Time_days <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
        # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
        
        r_exp <- nls_multstart(RFU ~ N0 * exp(r*Time_days),  # Exponential growth model (N0 is in our dataframe)
                               data = df.it.wl.th,
                               start_lower = c(r = -4.5), 
                               start_upper = c(r = 4.5),   
                               iter = 500,
                               supp_errors = 'Y',
                               control = nls.control(maxiter = 200))
        
        if (is.null(r_exp)){
          
          df.r.exp.l <- rbind(df.r.exp.l, data.frame(
            Sp.id = df.it.wl.th$species_updated[1],           
            light = df.it.wl$value_resource[1],  
            temp = df.it.wl$temperature[1],
            id = df.it.wl$sp.num[1],                
            r.exp = NA        
          ))
          
        }else{
          # Add data to our summary table
          df.r.exp.l <- rbind(df.r.exp.l, data.frame(
            Sp.id = df.it.wl.th$species_updated[1],    # Species
            light = df.it.wl$value_resource[1],        # Light
            temp = df.it.wl$temperature[1],            # Temperature
            id = df.it.wl$sp.num[1],                   # ID
            r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
          ))
        }
      
      }
      
    }
    
  }
  
}

write.csv(df.r.exp.l, "data-processed/403_Levasseur_2025_µ_estimates_light.csv") # let's save the file.

# Light Monods ------------------------------------------------------------

# df.r.exp.l <- read.csv("data-processed/403_Levasseur_2025_µ_estimates_light.csv") # if needed

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

levasseur.l.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species ID
  Sp.name = character(),    # Species name
  K.s = numeric(),          # Half-saturation constant
  r.max = numeric(),        # Maximum population growth rate
  R.jag = numeric(),        # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),        # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
  stringsAsFactors = FALSE  # Avoid factor conversion
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

for (i in 1:20){ # For each species
  
  df.i <- df.r.exp.l %>% 
    filter(id == i, !is.na(r.exp)) %>% 
    arrange(light)
  df.i <- droplevels(df.i)
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  
  light <- df.i$light
  
  S.pred <- seq(0, 250, 0.5) # Light gradient we're interested in - upped the granularity here
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
  
  print(paste("Done", i, "of 20"))
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3, (250/0.5) + 5),]   # generate the sequence of r.pred values
  df.jags$light <- seq(0, 250, 0.5)
  
  levasseur.l.summ.df <- rbind(levasseur.l.summ.df, data.frame( # Add summary data
    Sp.id = df.i$id[1],                                                                         # Species #
    Sp.name = df.i$Sp.id[1],                                                                    # Species name
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                    # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                  # Maximum population growth rate
    R.jag = df.jags$light[which(df.jags$mean > 0.1)[1]],                                        # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.1*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.1)     # Minimum resource requirement for positive growth (from math)
  ))
  
  light_sum <- monod_jag$BUGSoutput$summary[c(1:3, (250/0.5) + 5),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = df.i$id[1],                        # Species #
      Sp.name = df.i$Sp.id[1],                   # Species name     
      Parameter = rownames(light_sum)[j],        # Model parameter (e.g. K_s, r_max, etc.)
      mean = light_sum[j,1],                     # Posterior mean
      Rhat = light_sum[j,8],                     # Rhat values
      n.eff = light_sum[j,9]                     # Sample size estimates (should be ~3000)
    ))
    
  }
  
}

write.csv(levasseur.l.summ.df, "data-processed/404a_Levasseur_2025_light_monod.csv") # Save Levasseur 2025 light monod summary table
write.csv(fit.df, "data-processed/404b_Levasseur_2025_light_monod_fits.csv") # Save model fit summary table

# µ: nitrogen -------------------------------------------------------------

nit <- as.vector(unique(df.n$value_resource))# for looping through nitrogen levels
ord.nit<- sort(nit)

df.r.exp.n <- data.frame(              # Summary dataframe for r_exp
  Sp.id = character(),                 # Species
  nit = numeric(),                     # Nitrogen
  id = character(),                    # Unique ID
  r.exp = numeric()                    # Thresholded r.exp
)

for (i in 1:20){ # for every species
  
  for (t in ord.nit){ # at every nitrogen level
    
    df.i <- df.n[df.n$sp.num == i, ]
    df.it <- df.i[df.i$value_resource == t, ]
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in unique(df.it$unique.id)){ # Doing this separately for each replicate
      
      df.it.wl <- subset(df.it, df.it$unique.id == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$Time_days), ]
      df.it.wl <- df.it.wl %>% 
        mutate(N0 = RFU[1])
      
      t.series <- unique(df.it.wl$Time_days) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$Time_days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(log.RFU~Time_days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
      s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
      
      df.it.wl.th <- df.it.wl[df.it.wl$Time_days <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
      # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*Time_days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp.n <- rbind(df.r.exp.n, data.frame(
          Sp.id = df.it.wl.th$species_updated[1],           
          nit = df.it.wl$value_resource[1],        
          id = df.it.wl$sp.num[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp.n <- rbind(df.r.exp.n, data.frame(
          Sp.id = df.it.wl.th$species_updated[1],    # Species
          nit = df.it.wl$value_resource[1],          # Nitrogen
          id = df.it.wl$sp.num[1],                   # ID
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
      }
      
    }
    
  }
  
}

write.csv(df.r.exp.n, "data-processed/405_Levasseur_2025_µ_estimates_nitrogen.csv") # let's save the file.

# Nitrogen monods ---------------------------------------------------------

# df.r.exp.n <- read.csv("data-processed/405_Levasseur_2025_µ_estimates_nitrogen.csv") # if needed

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

levasseur.n.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species ID
  Sp.name = character(),    # Species name
  K.s = numeric(),          # Half-saturation constant
  r.max = numeric(),        # Maximum population growth rate
  R.jag = numeric(),        # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),        # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
  stringsAsFactors = FALSE  # Avoid factor conversion
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

for (i in 1:20){ # For each species
  
  df.i <- df.r.exp.n %>% 
    filter(id == i, !is.na(r.exp)) %>% 
    arrange(nit)
  df.i <- droplevels(df.i)
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  
  nit <- df.i$nit
  
  S.pred <- seq(0, 1000, 0.5) # Nitrogen gradient we're interested in - upped the granularity here
  N.S.pred <-length(S.pred) # We'll reset this internally since the gradient varies substantially
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = nit, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the nit Monod function. 
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
  
  print(paste("Done", i, "of 20"))
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3, (1000/0.5) + 5),]   # generate the sequence of r.pred values
  df.jags$light <- seq(0, 1000, 0.5)
  
  levasseur.l.summ.df <- rbind(levasseur.l.summ.df, data.frame( # Add summary data
    Sp.id = df.i$id[1],                                                                         # Species #
    Sp.name = df.i$Sp.id[1],                                                                    # Species name
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                    # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                  # Maximum population growth rate
    R.jag = df.jags$light[which(df.jags$mean > 0.1)[1]],                                        # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.1*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.1)     # Minimum resource requirement for positive growth (from math)
  ))
  
  nit_sum <- monod_jag$BUGSoutput$summary[c(1:3, (1000/0.5) + 5),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = df.i$id[1],                        # Species #
      Sp.name = df.i$Sp.id[1],                   # Species name     
      Parameter = rownames(nit_sum)[j],          # Model parameter (e.g. K_s, r_max, etc.)
      mean = nit_sum[j,1],                       # Posterior mean
      Rhat = nit_sum[j,8],                       # Rhat values
      n.eff = nit_sum[j,9]                       # Sample size estimates (should be ~3000)
    ))
    
  }
  
}

write.csv(levasseur.n.summ.df, "data-processed/406a_Levasseur_2025_nit_monod.csv") # Save Levasseur 2025 nit monod summary table
write.csv(fit.df, "data-processed/406b_Levasseur_2025_nit_monod_fits.csv") # Save model fit summary table
