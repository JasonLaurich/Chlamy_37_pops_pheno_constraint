# Jason R Laurich
# January 16, 2025

# This script will estimate r for each phosphate concentration as for light (script 06)
# I'm going to use a sliding-window approach to identify the exponential phase of the logged linear data, then fit an exponential growth
# curve to un-logged data for that time period for each phosphorous level.

# Then I'm going to fit Monod curves to those data (using R2jags) for each population to estimate R* (P*)

############# Packages ########################

library(nls.multstart)
library(tidyr)
library(cowplot)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(R2jags)
library(mcmcplots)

############# Upload and examine data #######################

df <- read.csv("data-processed/12_phosphate_rstar_rfus_time.csv")
head(df) # RFU is density, days is time, nitrate_level is a factor (1 to 10). nitrate_concentration is what we want
str(df)

df<-df[,-c(1,2,4,5,8:10,12:14)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df$nitrate.conc <- as.numeric(df$nitrate_concentration)
df$nitrate_level <- factor(df$nitrate_level, levels = sort(unique(df$nitrate_level)), ordered = TRUE) # Keep the numerical sorting.
df$well_plate <- as.factor(df$well_plate)

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

############# Loop through all populations ###################

############# Exploration #######################

############# Fit Monod curves to data ###################