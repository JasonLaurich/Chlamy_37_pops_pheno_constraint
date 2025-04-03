# Jason R Laurich
# April 3, 2025

# Going to work with the Bestion et al 2018 (https://onlinelibrary.wiley.com/doi/10.1111/ele.12932) data to fit
# Monod curves for phosphorous and TPCs

############# Packages ########################

library(tidyr)
library(cowplot)
library(ggplot2)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)
library(rTPC)
library(nls.multstart)

############# Upload and organize data #######################