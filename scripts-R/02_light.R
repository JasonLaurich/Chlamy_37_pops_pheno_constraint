# Jason R Laurich

# February 5th, 2026

# This file will fit µ to our time series growth data (RFUs) for C. reinhardtii populations for our light gradient data
# Then we will fit Monod curves to the data. 

# Inputs: in processed-data : 06_light_rfus_time.csv
# Outputs: in processed-data : 07_µ_estimates_light.csv,

# Packages & functions ----------------------------------------------------

library(tidyverse)

library(R2jags)
library(mcmcplots)
library(bayestestR)
library(cowplot)

# Upload & examine the data -----------------------------------------------

df <- read_csv("processed-data/06_light_rfus_time.csv")
head(df) #RFU is density, days is time, percentage is light in percentages

df$percentage <- as.numeric(df$percentage) # From an examination of the csv, the light level 2 corresponds to "0.5-0.7"
df$percentage[is.na(df$percentage)] <- 0.6 # Set to 0.6

df <- df %>% 
  mutate(light = as.numeric(percentage) * 2.5) %>% 
  select(well_plate, RFU, population, light, days) %>% 
  filter(population != "COMBO") # A control treatment, not relevant to this experimental analysis

unique(df$population)
unique(df$light)

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
  light = numeric(),                         # light
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
  
  if(n_distinct(df.i.th$RFU == 1)) {
    
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
    light = df.i$light[1],                            # light
    well.ID = df.i$well.ID[1],                        # well ID
    
    µ = µ.est                                         # intrinsic growth rate. 
  ))
  
  print(paste("Done", n, "of ", length(unique(df$well.ID))))
  
}

df.µ %>% 
  filter(is.na(µ))  # No null values

write.csv(df.µ, "processed-data/07_µ_estimates_light.csv",
          row.names = FALSE) # 1480 measurements
