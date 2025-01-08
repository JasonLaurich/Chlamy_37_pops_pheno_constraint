# Jason R Laurich
# January 7, 2025

# This script will estimate r for each light condition in the same way that I estimated it for the TPC analysis (04_final_r_estimation)
# I'm going to use a sliding-window approach to identify the exponential phase of the logged linear data, then fit an exponential growth
# curve to un-logged data for that time period for each light level.

# I'll also play with a model that incorporates lag, as that seems to have been potentially important for calculating R* in Joey's experiment. 

# Then I'm going to fit Monod curves to those data (using R2jags) for each population to estimate R* (I*)

############# Packages ########################

library(tidyr)
library(cowplot)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)

############# Upload and examine data #######################

df <- read.csv("data-processed/10_light_rstar_rfus_time.csv")
head(df) # RFU is density, days is time, light_level is a factor (1 to 10). Percentage is also a measurement of light I think?
str(df)

df<-df[,-c(1,2,4,5,8,10,12,13,15,16,17)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df$percentage <- as.numeric(df$percentage) # From an examination of the csv, the light level 2 corresponds to "0.5-0.7"
df$percentage[is.na(df$percentage)] <- 0.6 # For now, let's set this to 0.6, but I need to talk with Joey about this. 
df$light_level <- factor(df$light_level, levels = sort(unique(df$light_level)), ordered = TRUE) # Keep the numerical sorting.
df$well_plate <- as.factor(df$well_plate)

# df$days <- df$days + 0.001 # Can't have 0s
df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac) # I don't recognize the COMBO group, I'm guessing this is a control?

df.exp <- subset(df, df$pop.fac != "COMBO") # We're going to work only with repeat data (this is most of the well_plates) which I'll use as replicates. 
df.exp$well.ID<-as.factor(df.exp$well_plate)

N0.df <- df.exp %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df.exp <- df.exp %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.exp <- split(df.exp, df.exp$pop.num)  # Each element is a data frame for one population in df.exp

############# Data exploration #####################



############# Loop through all populations ###################


############# Fit Monod curves to data ###################

