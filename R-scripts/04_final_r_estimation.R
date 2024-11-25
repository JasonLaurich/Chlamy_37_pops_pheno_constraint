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

############## Data exploration #####################

df.it <- subset(mat.rep[[6]], temperature==34) # For one population and temperature, let's explore the model
df.it <- droplevels(df.it)

df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == 1)

head(df.it.wl)

t.series <- unique(df.it.wl$days) # First, we identify the string of unique days.
t.0 <- t.series[1] # Get the starting point
t.series <- t.series[-1] # Ok so this is what we will loop through! t.0 will always be the starting point, but the end point will change

ln.slopes <- c() # Store the logged linear slopes for each sliding window, from an lm
sl.direct <- c() # Store the directly calculated (rise/run)

for (t in t.series){ # This first loop will calculate the slopes for my first explorative population.
  
  df.it.wl.sl <- df.it.wl[df.it.wl$days <= t, ] # Subset the data to exclude time points above our window
  
  ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
  
  ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
  
  slope <- (df.it.wl.sl$logRFU[nrow(df.it.wl.sl)] - df.it.wl.sl$logRFU[1])/(df.it.wl.sl$days[nrow(df.it.wl.sl)] - df.it.wl.sl$days[1])
  sl.direct <- c(sl.direct, slope)
  
  for (s in 1:(length(ln.slopes) - 1)) { # Now we will loop through the ln.slopes vector to find when the slope decreases by more than 5%
    
    percent.chg <- (ln.slopes[s] - ln.slopes[s + 1]) / ln.slopes[s]
    if (percent.chg >= 0.05) {
      break  # Exit loop when condition is met. s will store the final time point data!
    }
  }
  
  # OK so we now have (in s) an indicator of when we should be thresholding growth data to limit it to the exponential phase.
  # This will work with ln.slopes AND t.series (because we excluded the t_0 data point from it)
}