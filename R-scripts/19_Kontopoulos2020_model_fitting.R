# Jason R Laurich
# April 14, 2025

# Going to work with the Kontopoulos et al 2020 (https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000894#abstract0) 
# data to fit TPCs for a bunch of species

# Load packages -----------------------------------------------------------

library(cowplot)
library(tidyverse)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)
library(rTPC)
library(nls.multstart)

# Load and examine data ---------------------------------------------------

df <- read.csv("data-processed/20_Narwani2015_summmary.csv") # Raw data file
head(df)
str(df)

df$Sp.fac <- as.factor(df$Species.name)

# TPC fitting -------------------------------------------------------------


