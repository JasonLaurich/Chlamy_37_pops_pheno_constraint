# Jason R Laurich
# March 5, 2025

# Modelling inter-gradient competitive/ thermal performance metrics here.

# Use same methods developed and tested in 12_gen-spec_glean-opp_toffs.R

# Start with examination of PCA, RDA data from 12
# Pick out some relationships to dive into

############# Packages ########################

library(dplyr)
library(ggplot2)
library(rPref)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(vegan)  # For PCA and RDA
library(ggrepel)
library(quantreg)
library(brms)

############# Upload and organize data #######################