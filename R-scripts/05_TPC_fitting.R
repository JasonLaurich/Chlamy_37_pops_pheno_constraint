# Jason R Laurich
# November 20, 2024

# Here, I will upload the final estimates of r for each population and well replicate at each temperature.
# And use them to fit thermal performance curves (TPCs).
# I will fit seven models in nls that reflect some of the best-performing phenomenological models, and that capture 
# much of the variance among available models (see Kontopoulos et al 2024, Nat Comm, https://doi.org/10.1038/s41467-024-53046-2)

# From the rTPC package: Ratkowsky (bounded), Rezende, Lactin 2, and modified Deutsch, and the Asrafi II, Mitchel-Angilletta and Atkin models (from scratch) 
# We'll calculate AICc values and plot the models, and select the best performing ones...
# Which we will then fit using bayesian methods in R2jags.

############# Packages ########################

library(nls.multstart)
library(tidyr)
library(cowplot)
library(ggplot2)

############# Upload and examine data #######################