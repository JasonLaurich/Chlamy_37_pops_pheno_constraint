# Jason R Laurich
# November 20, 2024

# This script will estimate r one final way (which is the way we intend to calculate estimates used for the fitting of TPCs)
# For all sub-40 C populations, I will fit logged linear models to our data, within a window beginning at our first time point.

# As the width of that window increases, I will calculate the slope of the line, and interpret a decrease in the slope as the
# end of the exponential growth phase.
# I will then fit an exponential growth curve to that period of data for each population and well replicate of Chlamydomonas.

# For 40 C, our data shows a brief growth spike that represents anomalous data (a completion of cell division already initiated)
# that is not biologically meaningful. 
# Thus, for this treatment, I will threshold data to after maximum RFU count, and fit the same exponential growth curve to it.

############# Packages ########################

library(nls.multstart)
library(tidyr)
library(cowplot)
library(ggplot2)
library(dplyr)

############# Upload and examine data #######################

df <- read.csv("data-processed/chlamee-acute-rfu-time.csv")
head(df) #RFU is density, days is time, temperature is temperature
str(df)

df<-df[,-c(1,2,4,5,6,8,13)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df<-subset(df,df$temperature!=20)

# df$days <- df$days + 0.001 # Can't have 0s
df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac) # I don't recognize the cc1629 group, but we'll keep it for now

df.rep <- subset(df, df$plate_type == "repeat") # We're going to work only with repeat data (this is most of the well_plates) which I'll use as replicates. 
df.rep$well.ID<-as.factor(df.rep$well_plate)

N0.df <- df.rep %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df.rep <- df.rep %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.rep <- split(df.rep, df.rep$pop.num)  # Each element is a data frame for one population in df.rep

