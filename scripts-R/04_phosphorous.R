# Jason R Laurich

# February 6th, 2026

# This file will fit µ to our time series growth data (RFUs) for C. reinhardtii populations for our phosphorous gradient data
# Then we will fit Monod curves to the data. 

# Inputs: in processed-data : 14_phosphorous_rfus_time.csv
# Outputs: in processed-data : 15_µ_estimates_phosphorous.csv, 16_phosphorous_monod_summary.csv, 17_phosphorous_monod_fits.csv

# Packages & functions ----------------------------------------------------

library(tidyverse)

library(R2jags)
library(mcmcplots)
library(bayestestR)
library(cowplot)

# Upload & examine the data -----------------------------------------------

df <- read_csv("processed-data/14_phosphorous_rfus_time.csv")
head(df) #RFU is density, days is time, phosphate_concentration does what it says on the tin

df <- df %>% 
  rename(phos = phosphate_concentration) %>% 
  select(well_plate, RFU, population, phos, days) %>% 
  filter(population != "COMBO") # A control treatment, not relevant to this experimental analysis

unique(df$population)
unique(df$phos)

df <- df %>% 
  
  mutate(logRFU = log(RFU + 0.001)) %>%   # take the log of RFUs for comparing log-linear slopes to identify the exponential growth phase
  
  rename(well.ID = well_plate)

N0.df <- df %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df <- df %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") 

# Estimate µ --------------------------------------------------------------

df.µ <- data.frame(                          # Initializing a data frame to store the results for each well, pop, and light
  
  population = character(),                  # population ID
  phos = numeric(),                          # phosphorous
  well.ID = character(),                     # well ID
  
  µ = numeric()                              # intrinsic growth rate estimate
)


n <- 0 # completion tracker

for (i in unique(df$well.ID[df$well.ID >= 1])) { # This allows code below to be run in chunks if necessary, just re-start at the next i. 
  
  n <- n + 1
  
  df.i <- df %>% 
    filter(well.ID == i)  # Filter data down to a single well time series
  
  df.i <- df.i[order(df.i$days), ]  # in case data are not ordered by time of observation
  
  t.series <- unique(df.i$days)           # Re-initialize this internally - we will only save summary data for each unique pop x light x well combo
  
  ln.slopes <- c() # Re-initialize this too!
  
  for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including time 0 
    
    df.i.sl <- df.i[df.i$days <= z, ]        # Subset the data to exclude time points above our window
    
    ln_slope <- lm(logRFU~days, data = df.i.sl) # calculate the log-linear slope
    
    ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
    
  }
  
  s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
  
  df.i.th <- df.i[df.i$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
  
  if(n_distinct(df.i.th$RFU) == 1) {
    
    µ.est <- 0    # If all of the thresholded values have the same RFU scores
    
  }else{
    
    µ.mod <- tryCatch(                          # Run the exponential growth model on the thresholded data
      nls_multstart(
        RFU ~ N0 * exp(r * days),
        data = df.i.th,
        start_lower = c(r = -4.5),
        start_upper = c(r = 4.5),
        iter = 500,
        supp_errors = "Y",
        control = nls.control(maxiter = 200)
      ),
      error = function(e) NULL
    )
    
    if (is.null(µ.mod)) {                   # if it fails, return NA
      µ.est <- NA_real_
    } else {
      µ.est <- coef(µ.mod)[["r"]]           # else return the model value for r
    } 
    
  } 
  
  df.µ <- rbind(df.µ, data.frame(                     # add the data to the summary data frame
    
    population = df.i$population[1],                  # population ID
    phos = df.i$phos[1],                              # phosphorous
    well.ID = df.i$well.ID[1],                        # well ID
    
    µ = µ.est                                         # intrinsic growth rate. 
  ))
  
  print(paste("Done", n, "of ", length(unique(df$well.ID))))
  
}

df.µ %>% 
  filter(is.na(µ))  # No null values

write.csv(df.µ, "processed-data/15_µ_estimates_phosphorous.csv",
          row.names = FALSE) # 1480 measurements

# Bayesian modelling ------------------------------------------------------

# Now we are going to fit Monod curves to each replicate using R2jags.

df.µ <- df.µ %>%
  mutate(rep.id = str_c(population, str_sub(well.ID, 1, 3), sep = ".")) # Need to create a unique replicate ID

length(unique(df.µ$rep.id)) # 156, should be 148! 

by_pop <- df.µ %>%
  distinct(population, rep.id) %>%   # one row per (pop, unique.id)
  count(population, name = "n_unique")  # how many unique.id per pop

by_pop # Pops 22 and anc2 are getting double the unique hits.

df.µ.22 <- df.µ[df.µ$population == "22",] # P20 is the issue

df.µ.anc2 <- df.µ[df.µ$population == "anc2",] # P20 is the issue
# These are perfectly swapped — but for now I'm going to remove them.  

df.µ <- df.µ %>%
  filter(!(population %in% c("22", "anc2") & phos == 20)) 

length(unique(df.µ$rep.id)) # 148, should be 1480/10 = 148! Fixed for now. 37 x 4 = 148 (pops, reps). Now missing some level replicates but that's OK, we deleted them. 

summary.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  
  population = character(), # population ID
  rep.ID = character(),     # rep ID
  
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
  
  R.post = numeric(),       # R* (m = 0.56) (posterior median)
  R.min = numeric(),        # R* (m = 0.56) (lower HDPI)
  R.max = numeric(),        # R* (m = 0.56) (upper HDPI)
  R.na = numeric(),         # % NA returns
  
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  
  population = character(), # population ID
  rep.ID = character(),     # rep ID
  
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~6000)
  stringsAsFactors = FALSE            
)

# Set generous MCMC settings still for our models. 
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model. 

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

S.pred <- seq(0, 50, 0.025) # Phosphorous gradient we are interested in here (concentration)
N.S.pred <-length(S.pred)

rep.ids <- unique(df.µ$rep.id) # Save the unique rep ids so we can run the for loop in chunks if needed

n <- 0 # for tracking progress

for (i in rep.ids[1:length(rep.ids)]) { # For all replicates, can account for multiple coding sessions to generate all of the objects
  
  n <- n + 1
  
  df.i <- df.µ %>% 
    filter(rep.id == i) %>% 
    filter(!is.na(µ))
  
  trait <- df.i$µ    # format the data for jags
  N.obs <- length(trait)
  
  phos <- df.i$phos
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = phos, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod.jag <- jags( # Run the phosphorous Monod function. 
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
  
  save(monod.jag, file = paste0("R2jags-models/rep_", i, "_phos_monod.RData")) # save the monod model
  # This folder is listed in gitignore, because the objects are too big to load
  
  post <- as.data.frame(monod.jag$BUGSoutput$sims.matrix) # The posteriors
  post$R <- 0.56*post$K_s/(post$r_max - 0.56)
  
  summary.df <- rbind(summary.df, data.frame(                                   # Add summary data
    
    population = df.i$population[1],                                            # population ID
    rep.ID = df.i$rep.id[1],                                                    # rep ID
    
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
    
    R.post = median(post$R, na.rm = T),                                         # R* (m = 0.56) (posterior median)
    R.min = hdi(post$R, ci = 0.95)$CI_low,                                      # R* (m = 0.56) (lower HDPI)
    R.max = hdi(post$R, ci = 0.95)$CI_high,                                     # R* (m = 0.56) (upper HDPI)
    R.na = mean(is.na(post$R))                                                  # % NA returns
    
  ))
  
  phos_sum <- monod.jag$BUGSoutput$summary[c(1:3),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:3){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      
      population = df.i$population[1],           # population ID
      rep.ID = df.i$rep.id[1],                   # rep ID    
      
      Parameter = rownames(phos_sum)[j],         # Model parameter (e.g. K_s, r_max, etc.)
      mean = phos_sum[j,1],                      # Posterior mean
      Rhat = phos_sum[j,8],                      # Rhat values
      n.eff = phos_sum[j,9]                      # Sample size estimates (should be ~3000)
    ))
    
  }
  
  print(paste("Done", n, "of ", length(unique(df.µ$rep.id))))
  
}

write.csv(summary.df, "processed-data/16_phos_monod_summary.csv",
          row.names = FALSE) # 148 measurements

write.csv(fit.df, "processed-data/17_phos_monod_fits.csv",
          row.names = FALSE) # 148 measurements
