# Jason Laurich
# Oct 24, 2024

################################################################################################################

# This script will upload Chlamydomonas growth (RFU) data, fit logarithmic growth curves to extract estimates of r
# for each of 37 populations at each of 6 temperatures, then fit TPCs to the data using nls and bayesian methods

######################### Upload packages ###################################

library(cowplot)
library(tidyverse)

######################### Upload and transform data #########################

df <- read.csv("data-processed/chlamee-acute-rfu-time.csv")
head(df) #RFU is density, days is time, temperature is temperature
str(df)

df$tfac <- as.factor(df$temperature)
df$pop <- as.factor(df$population)

str(df) # 38 populations?
levels(df$pop) # I don't recognize the cc1629 group, but we'll keep it for now

df$logRFU <- log(df$RFU) # Let's take the natural logarithm of RFU so we can fit an easy logarithmic growth curve.

######################### Estimate maximum growth rate, r ###################




######################### Fit TPCs using nls ################################



######################### Fit TPCs using JAGS ###############################