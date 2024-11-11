# Jason R Laurich
# November 11, 2024

# This script will fit growth rates to Chlamydomonas reinhardtii at different temperatures,
# using the growthTools package (https://github.com/ctkremer/growthTools/tree/master)

############# Packages #####################

remotes::install_github("ctkremer/mleTools")
remotes::install_github("ctkremer/growthTools")

library(growthTools)
library(ggplot2)
library(dplyr)
library(tidyr)

############# Upload and examine data

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

mat.rep <- split(df.rep, df.rep$pop.num)  # Each element is a data frame for one population in df.rep

############# Data exploration #############

# We'll eventually loop this through all of my temperatures and populations, but for now:

df.it <- subset(mat.rep[[3]], temperature==28)
df.it <- droplevels(df.it)

df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == 1)

res<-get.growth.rate(df.it.wl$days, df.it.wl$logRFU, plot.best.Q = T, id = 'Population 11, 28C, Well B07_46')

# Extract some summary statistics
res$best.model # gr.sat, no lag period.
res$best.slope
res$best.model.rsqr
res$best.se

# OK, let's think about modelling this for all temperatures on populations.
# We'll want to record.
  # 1. best.model - are they all the same? Except for T = 40?
  # 2. best.slope and best.se, obviously!
  # 3. We are not going to want to plot anything yet, we'll likely do this for a couple of populations later.