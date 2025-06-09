# Jason R Laurich

# June 9, 2025

# Going to upload data from the Edwards et al 2016 Limnol Oceanogr paper and run light Monods and TPCs on the raw data.


# Load packages -----------------------------------------------------------

library(tidyr)
library(cowplot)
library(tidyverse)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)
library(rTPC)
library(nls.multstart)

# Load and examine data ---------------------------------------------------

df <- read.csv("data-processed/27a_Edwards_2016_raw_data.csv") # Raw data file
head(df)
str(df)

df$Sp.fac <- as.factor(df$species)

# Lactin II TPCs ----------------------------------------------------------



# Light Monod curves ------------------------------------------------------


