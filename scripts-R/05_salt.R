# Jason R Laurich

# February 7th, 2026

# This file will fit µ to our time series growth data (RFUs) for C. reinhardtii populations for our salt gradient data
# Then we will fit salt tolerance (reverse logistic) curves to the data. 

# Inputs: in processed-data : 18_salt_rfus_time.csv
# Outputs: in processed-data : 19_µ_estimates_salt.csv, 20_salt_tol_summary.csv, 21_salt_tol_fits.csv
  # in R2jags-models : pop_i_salt_tol.RData (Salt tolerance R2jags objects, saved but not pushed)

# Packages & functions ----------------------------------------------------

library(tidyverse)

library(R2jags)
library(mcmcplots)
library(bayestestR)

# Upload & examine the data -----------------------------------------------

df <- read_csv("processed-data/18_salt_rfus_time.csv")
head(df) #RFU is density, days is time, treatment (S0 to S9) lists the concentrations in g/L

df$salt_level <- factor(df$treatment, levels = sort(unique(df$treatment)), ordered = TRUE) # Keep the numerical sorting.
df$salt <- as.numeric(df$salt_level) - 1# Subtract 1 to make it run from 0 to 9

df <- df %>% 
  select(RFU, population, salt, time) # For salt, we only have one estimate for each population. 

unique(df$population)
unique(df$salt)

df <- df %>% 
  
  mutate(logRFU = log(RFU + 0.001))  %>% # take the log of RFUs for comparing log-linear slopes to identify the exponential growth phase
  mutate(pop.ID = paste0(population, ".", salt))

length(unique(df$pop.ID)) # 370, which is correct. 37 populations x 10 salt levels

N0.df <- df %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(pop.ID) %>%
  summarize(N0 = RFU[which.min(time)]) %>%
  ungroup()

df <- df %>% # Recombine this with our dataframe
  left_join(N0.df, by = "pop.ID") 

# Estimate µ --------------------------------------------------------------

df.µ <- data.frame(                          # Initializing a data frame to store the results for each pop and salt
  
  population = character(),                  # population ID
  salt = numeric(),                          # salt
  pop.ID = character(),                      # pop ID
  
  µ = numeric()                              # intrinsic growth rate estimate
)


n <- 0 # completion tracker

for (i in unique(df$pop.ID[df$pop.ID >= 1])) { # This allows code below to be run in chunks if necessary, just re-start at the next i. 
  
  n <- n + 1
  
  df.i <- df %>% 
    filter(pop.ID == i)  # Filter data down to a single well time series
  
  df.i <- df.i[order(df.i$time), ]  # in case data are not ordered by time of observation
  
  t.series <- unique(df.i$time)           # Re-initialize this internally - we will only save summary data for each unique pop x salt 
  
  ln.slopes <- c() # Re-initialize this too!
  
  for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including time 0 
    
    df.i.sl <- df.i[df.i$time <= z, ]        # Subset the data to exclude time points above our window
    
    ln_slope <- lm(logRFU~time, data = df.i.sl) # calculate the log-linear slope
    
    ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
    
  }
  
  s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
  
  df.i.th <- df.i[df.i$time <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
  
  if(n_distinct(df.i.th$RFU) == 1) {
    
    µ.est <- 0    # If all of the thresholded values have the same RFU scores
    
  }else{
    
    µ.mod <- tryCatch(                          # Run the exponential growth model on the thresholded data
      nls_multstart(
        RFU ~ N0 * exp(r * time),
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
    salt = df.i$salt[1],                              # salt
    pop.ID = df.i$pop.ID[1],                          # pop ID
    
    µ = µ.est                                         # intrinsic growth rate. 
  ))
  
  print(paste("Done", n, "of ", length(unique(df$pop.ID))))
  
}

df.µ %>% 
  filter(is.na(µ))  # No null values

write.csv(df.µ, "processed-data/19_µ_estimates_salt.csv",
          row.names = FALSE) # 370 measurements

# Bayesian modelling ------------------------------------------------------

# Now we are going to fit salt tolerance curves to each population using R2jags.

inits.salt <- function() { # Smaller a (prior 0.5-> 2)
  list(
    a = runif(1, 0.7, 1.8),  # Smaller window for a
    b = runif(1, 0.1, 5),  # Initial guess for b
    c = runif(1, 0.1, 9),  # Initial guess for c
    sigma = runif(1, 0.1, 2)  # Initial guess for error
  )
}

parameters.salt <- c("a", "b", "c", "sigma", "r_pred_new") # Repeats

S.pred <- seq(0, 10, 0.005) # Repeats
N.S.pred <-length(S.pred)

# Repeated here
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model.

length(unique(df.µ$population)) # 37, should be 37!

summary.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  
  population = character(), # population ID
  
  r.max.mod = numeric(),    # Maximum population growth rate (alpha) (model output)
  r.max.post = numeric(),   # Maximum population growth rate (posterior median)
  r.max.min = numeric(),    # Maximum population growth rate (lower HDPI)
  r.max.max = numeric(),    # Maximum population growth rate (upper HDPI)
  r.max.na = numeric(),     # % NA returns
  
  c.mod = numeric(),        # salt concentration at which r is half of alpha (extracted from model)
  c.post = numeric(),       # salt concentration at which r is half of alpha (posterior median)
  c.post.min = numeric(),   # salt concentration at which r is half of alpha (lower HDPI)
  c.post.max = numeric(),   # salt concentration at which r is half of alpha (upper HDPI)
  c.post.na = numeric(),    # % NA returns
  
  b.mod = numeric(),        # decline rate (extracted from model)
  b.post = numeric(),       # decline rate (posterior median)
  b.post.min = numeric(),   # decline rate (lower HDPI)
  b.post.max = numeric(),   # decline rate (upper HDPI)
  b.post.na = numeric(),    # % NA returns
  
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  
  population = character(), # population ID
  
  Parameter = character(),  # Model parameter
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~6000)
  stringsAsFactors = FALSE            
)

rep.ids <- unique(df.µ$population) # Save the unique pop ids so we can run the for loop in chunks if needed

n <- 0 # for tracking progress

for (i in rep.ids[1:length(rep.ids)]) { # For all replicates, can account for multiple coding sessions to generate all of the objects
  
  n <- n + 1
  
  df.i <- df.µ %>% 
    filter(population == i) %>% 
    filter(!is.na(µ))
  
  trait <- df.i$µ    # format the data for jags
  N.obs <- length(trait)
  
  salt <- df.i$salt
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)
  
  salt.jag <- jags( # Run the salt logistic growth curve function. 
    data = jag.data,
    inits = inits.salt,
    parameters.to.save = parameters.salt,
    model.file = "salt.tolerance.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )

  save(salt.jag, file = paste0("R2jags-models/pop_", i, "_salt_tol.RData")) # save the monod model
  # This folder is listed in gitignore, because the objects are too big to load
  
  post <- as.data.frame(salt.jag$BUGSoutput$sims.matrix) # The posteriors
  
  summary.df <- rbind(summary.df, data.frame(                                   # Add summary data
    
    population = df.i$population[1],                                            # population ID
    
    r.max.mod = salt.jag$BUGSoutput$summary[1,1],                               # Maximum population growth rate (alpha) (model output)
    r.max.post = median(post$a, na.rm = T),                                     # Maximum population growth rate (posterior median)
    r.max.min = hdi(post$a, ci = 0.95)$CI_low,                                  # Maximum population growth rate (lower HDPI)
    r.max.max = hdi(post$a, ci = 0.95)$CI_high,                                 # Maximum population growth rate (upper HDPI)
    r.max.na = mean(is.na(post$a)),                                             # % NA returns
    
    c.mod = salt.jag$BUGSoutput$summary[3,1],                                   # salt concentration at which r is half of alpha (extracted from model)
    c.post = median(post$c, na.rm = T),                                         # salt concentration at which r is half of alpha (posterior median)
    c.post.min = hdi(post$c, ci = 0.95)$CI_low,                                 # salt concentration at which r is half of alpha (lower HDPI)
    c.post.max = hdi(post$c, ci = 0.95)$CI_high,                                # salt concentration at which r is half of alpha (upper HDPI)
    c.post.na = mean(is.na(post$c)),                                            # % NA returns
    
    b.mod = salt.jag$BUGSoutput$summary[2,1],                                   # decline rate (extracted from model)
    b.post = median(post$b, na.rm = T),                                         # decline rate (posterior median)
    b.post.min = hdi(post$b, ci = 0.95)$CI_low,                                 # decline rate (lower HDPI)
    b.post.max = hdi(post$b, ci = 0.95)$CI_high,                                # decline rate (upper HDPI)
    b.post.na = mean(is.na(post$b))                                             # % NA returns
    
  ))
  
  salt_sum <- salt.jag$BUGSoutput$summary[1:4,] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      
      population = df.i$population[1],           # population ID
      
      Parameter = rownames(salt_sum)[j],         # Model parameter (e.g. K_s, r_max, etc.)
      mean = salt_sum[j,1],                      # Posterior mean
      Rhat = salt_sum[j,8],                      # Rhat values
      n.eff = salt_sum[j,9]                      # Sample size estimates (should be ~6000)
    ))
    
  }
  
  print(paste("Done", n, "of ", length(unique(df.µ$population))))
  
}  

write.csv(summary.df, "processed-data/20_salt_tol_summary.csv",
          row.names = FALSE) # 37 measurements

write.csv(fit.df, "processed-data/21_salt_tol_fits.csv",
          row.names = FALSE) # 37 measurements
