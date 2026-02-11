# Jason R Laurich

# February 11th, 2026

# In this script I will fit light monod curves to the synthesis data set, working with estimates from a wide range of phytoplankton and studies
# In some cases, I will work directly with published estimates of µ, while in others will have to first estimate µ from 
# time-series data (as in scripts 01-05.R). 

# Where TPCs are available (fit within the same study), I will fit Monod curves at the temperature closest to Topt.
  # Else I will extract an approximate/putative Topt (temperature which has the highest growth rate) to use. 

# Inputs: 34_Lewington-Pearce_2019_raw_data.csv, 38_Edwards_2016_raw_data.csv, 39_Edwards_2016_TPCs.csv, 
  # 41_Levasseur2025_l_rawdata.csv, 45_Levasseur_2025_TPCs.csv
# Outputs: in processed-data : 47_Edwards_2016_light_monods.csv, 48_Edwards_2016_light_monods_fits.csv, 49_Lewington-Pearce_2019_µ_estimates_light.csv,
  # 50_Lewington_2019_light_monods.csv, 51_Lewington_2019_light_monods_fits.csv, 52_Levasseur_2025_µ_estimates_light.csv,
  # 53_Levasseur_2025_light_monods.csv, 54_Levasseur_2025_light_monods_fits.csv

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(Deriv)
library(car)
library(minpack.lm)
library(cowplot)

library(R2jags)
library(mcmcplots)
library(bayestestR)

# And R2jags settings
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

# Edwards 2016 ------------------------------------------------------------

###### Load the data ######

df.e.raw <- read.csv("processed-data/38_Edwards_2016_raw_data.csv") # Raw data file
head(df.e.raw)

df.e.raw <- df.e.raw %>% # The Edwards data has some species that show up in multiple references, so we need to create a unique sp.idx ref# combo
  mutate(unique.id = paste(species, reference, sep = "_")) # So that each entry is treated separately!

df.e <-df.e.raw %>% 
  rename(temp = temperature,
         light = irradiance,
         mu = growth.rate,
         Species.name = species) %>% 
  mutate(id.number = as.integer(factor(df.e.raw$unique.id)))

df.e <- df.e %>% 
  filter(light<= 250) # After looking at the light models, insanely high lights are driving down growth rates. We'll cap this out at 250 to improve model fits.

df.e.t <- read.csv("processed-data/39_Edwards_2016_TPCs.csv") # Raw data file

###### Model fitting ######

edwards.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Study = character(),      # Study
  
  K.s.mod = numeric(),      # Half saturation constant (model output)
  K.s.post = numeric(),     # Half saturation constant (posterior median)
  K.s.min = numeric(),      # Half saturation constant (lower HDPI)
  K.s.max = numeric(),      # Half saturation constant (upper HDPI)
  K.s.na = numeric(),       # % NA returns
  
  r.max.mod = numeric(),    # Maximum growth rate (model output)
  r.max.post = numeric(),   # Maximum growth rate (posterior median)
  r.max.min = numeric(),    # Maximum growth rate (lower HDPI)
  r.max.max = numeric(),    # Maximum growth rate (upper HDPI)
  r.max.na = numeric(),     # % NA returns
  
  R.post = numeric(),       # R* (m = 0.1) (posterior median)
  R.min = numeric(),        # R* (m = 0.1) (lower HDPI)
  R.max = numeric(),        # R* (m = 0.1) (upper HDPI)
  R.na = numeric(),         # % NA returns
  
  T.opt.TPC = numeric(),    # Topt (if there is a fitted TPC) — temp closest to this
  T.max.raw = numeric(),    # ~ Topt (temp in the raw data with highest growth rate)
  
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

n <-0 # progression tracker

for (i in unique(df.e$id.number[df.e$id.number >= 1])) { # for each species ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.e %>% 
    filter(id.number == i) %>% 
    filter(!is.na(mu))
  
  t.max <- df.i %>%                             # Find the temperature with the highest growth data
    filter(mu == max(mu, na.rm = TRUE)) %>% 
    pull(temp)
  
  if (i %in% df.e.t$Sp.id) {                    # pull T.opt from external dataset
    t.opt <- df.e.t %>%
      filter(Sp.id == i) %>%
      pull(T.opt) 
    fit <- "y"
  } else {                                      # fallback: Topt from raw data
    t.opt <- t.max
    fit <- "n"
  }
  
  t.sel <- df.i$temp[which.min(abs(df.i$temp - t.opt))] # Pick the temperature closest to Topt
  
  df.i <- df.i %>% 
    filter(temp == t.sel)
  
  trait <- df.i$mu    # format the data for jags
  N.obs <- length(trait)
  
  light <- df.i$light
  
  S.pred <- seq(0, max(df.i$light) + 25, 0.5) # Light gradient we're interested in - upped the granularity here
  N.S.pred <-length(S.pred) # We'll reset this internally since the gradient varies substantially
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = light, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod.jag <- jags( # Run the light Monod function. 
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
  
  post <- as.data.frame(monod.jag$BUGSoutput$sims.matrix) # The posteriors
  post$R <- 0.1*post$K_s/(post$r_max - 0.1)
  
  edwards.summ.df <- rbind(edwards.summ.df, data.frame(                         # Add data
    Sp.id = df.i$id.number[1],                                                  # Species id
    Sp.name = df.i$Species.name[1],                                             # Species name
    Study = df.i$reference[1],                                                  # Study
    
    K.s.mod = monod.jag$BUGSoutput$summary[1,1],                                # Half saturation constant (model output)
    K.s.post = median(post$K_s, na.rm = T),                                     # Half saturation constant (posterior median)
    K.s.min = hdi(post$K_s, ci = 0.95)$CI_low,                                  # Half saturation constant (lower HDPI)
    K.s.max = hdi(post$K_s, ci = 0.95)$CI_high,                                 # Half saturation constant (upper HDPI)
    K.s.na = mean(is.na(post$K_s)),                                             # % NA returns
    
    r.max.mod = monod.jag$BUGSoutput$summary[3,1],                              # Maximum growth rate (model output)
    r.max.post = median(post$r_max, na.rm = T),                                 # Maximum growth rate (posterior median)
    r.max.min = hdi(post$r_max, ci = 0.95)$CI_low,                              # Maximum growth rate (lower HDPI)
    r.max.max = hdi(post$r_max, ci = 0.95)$CI_high,                             # Maximum growth rate (upper HDPI)
    r.max.na = mean(is.na(post$r_max)),                                         # % NA returns
    
    R.post = median(post$R, na.rm = T),                                         # R* (m = 0.1) (posterior median)
    R.min = hdi(post$R, ci = 0.95)$CI_low,                                      # R* (m = 0.1) (lower HDPI)
    R.max = hdi(post$R, ci = 0.95)$CI_high,                                     # R* (m = 0.1) (upper HDPI)
    R.na = mean(is.na(post$R)),                                                 # % NA returns
    
    T.opt.TPC = ifelse(fit == "y", t.opt, NA_real_),                            # Topt (if there is a fitted TPC) — temp closest to this
    T.max.raw = t.max                                                           # ~ Topt (temp in the raw data with highest growth rate)
  ))
  
  light_sum <- monod.jag$BUGSoutput$summary[c(1:3, (max(df.i$light) - S.pred[1])/0.5 + 11),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = df.i$id.number[1],                 # Species #
      Sp.name = df.i$Species.name[1],            # Species name      
      Parameter = rownames(light_sum)[j],        # Model parameter (e.g. K_s, r_max, etc.)
      mean = light_sum[j,1],                     # Posterior mean
      Rhat = light_sum[j,8],                     # Rhat values
      n.eff = light_sum[j,9]                     # Sample size estimates (should be ~3000)
    ))
    
  }
  
  print(paste("Done", n, "of ", length(unique(df.e$id.number))))
  
}

write.csv(edwards.summ.df, "processed-data/47_Edwards_2016_light_monods.csv") # Edwards_2016 summary table
write.csv(fit.df, "processed-data/48_Edwards_2016_light_monods_fits.csv") # Save model fit summary table

# Lewington-Pearce 2019 ------------------------------------------------------------

###### Load the data ######

df <- read.csv("processed-data/34_Lewington-Pearce_2019_raw_data.csv") # Raw data file
head(df)

df$Sp.fac <- as.factor(df$genus_species)

df$id <- paste0(df$block, ".", df$replicate)  # A unique identifier for each replicate.

df$log.fluorescence <- log(df$fluorescence + 0.001)

head(df)

###### Calculate µ across light levels ######

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

write.csv(df.r.exp.l, "processed-data/49_Lewington-Pearce_2019_µ_estimates_light.csv") # let's save the file.

df.l.raw <- read.csv("processed-data/49_Lewington-Pearce_2019_µ_estimates_light.csv") # Raw data file

head(df.l.raw)

df.l <-df.l.raw %>% 
  rename(mu = r.exp,
         Species.name = Sp.id) %>% 
  mutate(id.number = as.integer(factor(df.l.raw$Sp.id)))

length(unique(df.l$id.number))

df.l.t <- read.csv("processed-data/36_Lewington_2019_TPCs.csv") # Raw data file

###### Model fitting ######

lewington.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  
  K.s.mod = numeric(),      # Half saturation constant (model output)
  K.s.post = numeric(),     # Half saturation constant (posterior median)
  K.s.min = numeric(),      # Half saturation constant (lower HDPI)
  K.s.max = numeric(),      # Half saturation constant (upper HDPI)
  K.s.na = numeric(),       # % NA returns
  
  r.max.mod = numeric(),    # Maximum growth rate (model output)
  r.max.post = numeric(),   # Maximum growth rate (posterior median)
  r.max.min = numeric(),    # Maximum growth rate (lower HDPI)
  r.max.max = numeric(),    # Maximum growth rate (upper HDPI)
  r.max.na = numeric(),     # % NA returns
  
  R.post = numeric(),       # R* (m = 0.1) (posterior median)
  R.min = numeric(),        # R* (m = 0.1) (lower HDPI)
  R.max = numeric(),        # R* (m = 0.1) (upper HDPI)
  R.na = numeric(),         # % NA returns
  
  T.opt.TPC = numeric(),    # Topt (if there is a fitted TPC) — temp closest to this
  T.max.raw = numeric(),    # ~ Topt (temp in the raw data with highest growth rate)
  
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

n <-0 # progression tracker

for (i in unique(df.l$id.number[df.l$id.number >= 1])) { # for each species ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.l %>% 
    filter(id.number == i) %>% 
    filter(!is.na(mu))
  
  t.max <- df.i %>%                             # Find the temperature with the highest growth data
    filter(mu == max(mu, na.rm = TRUE)) %>% 
    pull(temp)
  
  if (i %in% df.l.t$Sp.id) {                    # pull T.opt from external dataset
    t.opt <- df.l.t %>%
      filter(Sp.id == i) %>%
      pull(T.opt) 
    fit <- "y"
  } else {                                      # fallback: Topt from raw data
    t.opt <- t.max
    fit <- "n"
  }
  
  t.sel <- df.i$temp[which.min(abs(df.i$temp - t.opt))] # Pick the temperature closest to Topt
  
  df.i <- df.i %>% 
    filter(temp == t.sel)
  
  trait <- df.i$mu    # format the data for jags
  N.obs <- length(trait)
  
  light <- df.i$light
  
  S.pred <- seq(0, max(df.i$light) + 25, 0.5) # Light gradient we're interested in - upped the granularity here
  N.S.pred <-length(S.pred) # We'll reset this internally since the gradient varies substantially
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = light, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod.jag <- jags( # Run the light Monod function. 
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
  
  post <- as.data.frame(monod.jag$BUGSoutput$sims.matrix) # The posteriors
  post$R <- 0.1*post$K_s/(post$r_max - 0.1)
  
  lewington.summ.df <- rbind(lewington.summ.df, data.frame(                     # Add data
    Sp.id = df.i$id.number[1],                                                  # Species id
    Sp.name = df.i$Species.name[1],                                             # Species name
    
    K.s.mod = monod.jag$BUGSoutput$summary[1,1],                                # Half saturation constant (model output)
    K.s.post = median(post$K_s, na.rm = T),                                     # Half saturation constant (posterior median)
    K.s.min = hdi(post$K_s, ci = 0.95)$CI_low,                                  # Half saturation constant (lower HDPI)
    K.s.max = hdi(post$K_s, ci = 0.95)$CI_high,                                 # Half saturation constant (upper HDPI)
    K.s.na = mean(is.na(post$K_s)),                                             # % NA returns
    
    r.max.mod = monod.jag$BUGSoutput$summary[3,1],                              # Maximum growth rate (model output)
    r.max.post = median(post$r_max, na.rm = T),                                 # Maximum growth rate (posterior median)
    r.max.min = hdi(post$r_max, ci = 0.95)$CI_low,                              # Maximum growth rate (lower HDPI)
    r.max.max = hdi(post$r_max, ci = 0.95)$CI_high,                             # Maximum growth rate (upper HDPI)
    r.max.na = mean(is.na(post$r_max)),                                         # % NA returns
    
    R.post = median(post$R, na.rm = T),                                         # R* (m = 0.1) (posterior median)
    R.min = hdi(post$R, ci = 0.95)$CI_low,                                      # R* (m = 0.1) (lower HDPI)
    R.max = hdi(post$R, ci = 0.95)$CI_high,                                     # R* (m = 0.1) (upper HDPI)
    R.na = mean(is.na(post$R)),                                                 # % NA returns
    
    T.opt.TPC = ifelse(fit == "y", t.opt, NA_real_),                            # Topt (if there is a fitted TPC) — temp closest to this
    T.max.raw = t.max                                                           # ~ Topt (temp in the raw data with highest growth rate)
  ))
  
  light_sum <- monod.jag$BUGSoutput$summary[c(1:3, (max(df.i$light) - S.pred[1])/0.5 + 11),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = df.i$id.number[1],                 # Species #
      Sp.name = df.i$Species.name[1],            # Species name      
      Parameter = rownames(light_sum)[j],        # Model parameter (e.g. K_s, r_max, etc.)
      mean = light_sum[j,1],                     # Posterior mean
      Rhat = light_sum[j,8],                     # Rhat values
      n.eff = light_sum[j,9]                     # Sample size estimates (should be ~3000)
    ))
    
  }
  
  print(paste("Done", n, "of ", length(unique(df.l$id.number))))
  
}

write.csv(lewington.summ.df, "processed-data/50_Lewington_2019_light_monods.csv") # Lewington_2019 summary table
write.csv(fit.df, "processed-data/51_Lewington_2019_light_monods_fits.csv") # Save model fit summary table

# Levasseur 2025 ------------------------------------------------------------

###### Load the data ######

df.l <- read.csv("processed-data/41_Levasseur2025_l_rawdata.csv") # Raw data file
head(df.l)
str(df.l)

df.l$sp.num <- as.integer(factor(df.l$species_updated))
df.l$log.RFU <- log(df.l$RFU + 0.001)

df.l <- df.l %>% # We'll concatenate well id and resource level just in case the same well appears in different experiments
  mutate(unique.id = paste(well_id, level_of_resource, sep = "_")) # So that each entry is treated seperately!

light <- as.vector(unique(df.l$value_resource))# for looping through lights
ord.light<- sort(light)

temp <- as.vector(unique(df.l$temperature))
ord.temp <- sort(temp)

###### Estimate µ across light levels ######

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

write.csv(df.r.exp.l, "processed-data/52_Levasseur_2025_µ_estimates_light.csv") # let's save the file.

df.lv.raw <- read.csv("processed-data/52_Levasseur_2025_µ_estimates_light.csv") # Raw data file

head(df.lv.raw)

df.lv <-df.lv.raw %>% 
  rename(mu = r.exp,
         Species.name = Sp.id) %>% 
  mutate(id.number = as.integer(factor(df.lv.raw$Sp.id)))

length(unique(df.lv$id.number))

df.lv.t <- read.csv("processed-data/45_Levasseur_2025_TPCs.csv") # Raw data file

###### Model fitting ######

levasseur.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  
  K.s.mod = numeric(),      # Half saturation constant (model output)
  K.s.post = numeric(),     # Half saturation constant (posterior median)
  K.s.min = numeric(),      # Half saturation constant (lower HDPI)
  K.s.max = numeric(),      # Half saturation constant (upper HDPI)
  K.s.na = numeric(),       # % NA returns
  
  r.max.mod = numeric(),    # Maximum growth rate (model output)
  r.max.post = numeric(),   # Maximum growth rate (posterior median)
  r.max.min = numeric(),    # Maximum growth rate (lower HDPI)
  r.max.max = numeric(),    # Maximum growth rate (upper HDPI)
  r.max.na = numeric(),     # % NA returns
  
  R.post = numeric(),       # R* (m = 0.1) (posterior median)
  R.min = numeric(),        # R* (m = 0.1) (lower HDPI)
  R.max = numeric(),        # R* (m = 0.1) (upper HDPI)
  R.na = numeric(),         # % NA returns
  
  T.opt.TPC = numeric(),    # Topt (if there is a fitted TPC) — temp closest to this
  T.max.raw = numeric(),    # ~ Topt (temp in the raw data with highest growth rate)
  
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

n <-0 # progression tracker

for (i in unique(df.lv$id.number[df.lv$id.number >= 1])) { # for each species ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.lv %>% 
    filter(id.number == i) %>% 
    filter(!is.na(mu))
  
  t.max <- df.i %>%                             # Find the temperature with the highest growth data
    filter(mu == max(mu, na.rm = TRUE)) %>% 
    pull(temp)
  
  if (i %in% df.lv.t$Sp.id) {                    # pull T.opt from external dataset
    t.opt <- df.lv.t %>%
      filter(Sp.id == i) %>%
      pull(T.opt) 
    fit <- "y"
  } else {                                      # fallback: Topt from raw data
    t.opt <- t.max
    fit <- "n"
  }
  
  t.sel <- df.i$temp[which.min(abs(df.i$temp - t.opt))] # Pick the temperature closest to Topt
  
  df.i <- df.i %>% 
    filter(temp == t.sel)
  
  trait <- df.i$mu    # format the data for jags
  N.obs <- length(trait)
  
  light <- df.i$light
  
  S.pred <- seq(0, max(df.i$light) + 25, 0.5) # Light gradient we're interested in - upped the granularity here
  N.S.pred <-length(S.pred) # We'll reset this internally since the gradient varies substantially
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = light, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod.jag <- jags( # Run the light Monod function. 
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
  
  post <- as.data.frame(monod.jag$BUGSoutput$sims.matrix) # The posteriors
  post$R <- 0.1*post$K_s/(post$r_max - 0.1)
  
  levasseur.summ.df <- rbind(levasseur.summ.df, data.frame(                     # Add data
    Sp.id = df.i$id.number[1],                                                  # Species id
    Sp.name = df.i$Species.name[1],                                             # Species name
    
    K.s.mod = monod.jag$BUGSoutput$summary[1,1],                                # Half saturation constant (model output)
    K.s.post = median(post$K_s, na.rm = T),                                     # Half saturation constant (posterior median)
    K.s.min = hdi(post$K_s, ci = 0.95)$CI_low,                                  # Half saturation constant (lower HDPI)
    K.s.max = hdi(post$K_s, ci = 0.95)$CI_high,                                 # Half saturation constant (upper HDPI)
    K.s.na = mean(is.na(post$K_s)),                                             # % NA returns
    
    r.max.mod = monod.jag$BUGSoutput$summary[3,1],                              # Maximum growth rate (model output)
    r.max.post = median(post$r_max, na.rm = T),                                 # Maximum growth rate (posterior median)
    r.max.min = hdi(post$r_max, ci = 0.95)$CI_low,                              # Maximum growth rate (lower HDPI)
    r.max.max = hdi(post$r_max, ci = 0.95)$CI_high,                             # Maximum growth rate (upper HDPI)
    r.max.na = mean(is.na(post$r_max)),                                         # % NA returns
    
    R.post = median(post$R, na.rm = T),                                         # R* (m = 0.1) (posterior median)
    R.min = hdi(post$R, ci = 0.95)$CI_low,                                      # R* (m = 0.1) (lower HDPI)
    R.max = hdi(post$R, ci = 0.95)$CI_high,                                     # R* (m = 0.1) (upper HDPI)
    R.na = mean(is.na(post$R)),                                                 # % NA returns
    
    T.opt.TPC = ifelse(fit == "y", t.opt, NA_real_),                            # Topt (if there is a fitted TPC) — temp closest to this
    T.max.raw = t.max                                                           # ~ Topt (temp in the raw data with highest growth rate)
  ))
  
  light_sum <- monod.jag$BUGSoutput$summary[c(1:3, (max(df.i$light) - S.pred[1])/0.5 + 11),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = df.i$id.number[1],                 # Species #
      Sp.name = df.i$Species.name[1],            # Species name      
      Parameter = rownames(light_sum)[j],        # Model parameter (e.g. K_s, r_max, etc.)
      mean = light_sum[j,1],                     # Posterior mean
      Rhat = light_sum[j,8],                     # Rhat values
      n.eff = light_sum[j,9]                     # Sample size estimates (should be ~3000)
    ))
    
  }
  
  print(paste("Done", n, "of ", length(unique(df.lv$id.number))))
  
}

write.csv(levasseur.summ.df, "processed-data/53_Levasseur_2025_light_monods.csv") # Levasseur_2025 summary table
write.csv(fit.df, "processed-data/54_Levasseur_2025_light_monods_fits.csv") # Save model fit summary table