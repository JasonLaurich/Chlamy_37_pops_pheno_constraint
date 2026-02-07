# Jason R Laurich

# February 7th, 2026

# This file will fit µ to our time series growth data (RFUs) for C. reinhardtii populations for our salt gradient data
# Then we will fit salt tolerance (reverse logistic) curves to the data. 

# Inputs: in processed-data : 18_salt_rfus_time.csv
# Outputs: in processed-data : 19_µ_estimates_salt.csv, 20_salt_tol_summary.csv, 21_salt_tol_fits.csv
  # in R2jags-models : rep_i_salt.RData (Salt tolerance R2jags objects, saved but not pushed)

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
