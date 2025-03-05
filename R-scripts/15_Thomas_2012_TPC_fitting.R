# Jason R Laurich
# March 5, 2025

# In this script, I will work with the Thomas 2012 Science paper (doi: DOI: 10.1126/science.122483) supplement containing specific grwoth rates for 194 algae.

############# Packages ########################

library(tidyr)
library(cowplot)
library(ggplot2)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)

############# Upload and organize data #######################

df <- read.csv("data-processed/16b_thomas_2012_supp_raw.csv") # Summary file
head(df)
str(df)

length(unique(df$id.number))

mat <- split(df, df$id.number)  # Matrix

df2 <- read.csv("data-processed/16a_thomas_2012_supp_data_summ.csv") # Summary file

min(df2$Optimum)
max(df2$Optimum) # ~ Topt 1 to 39

min(df[df$Growth.rate>0,]$Temperature) # Tmin <= -1.8
max(df[df$Growth.rate>0,]$Temperature) # Tmax => 37

############# Test drive #######################

i <- sample(1:194,1) 

df.i <- subset(mat[[i]])

inits.lactin.thomas<- function() { # The final initial values set we landed on after experimenting. 
  list(
    cf.a = runif(1, 0.05, 0.15),  # More constrained initial values
    cf.tmax = runif(1, 1, 40),    # Much wider to accomodate interspecific variation
    cf.delta_t = runif(1, 1, 5),
    cf.b = runif(1, -2.5, -1),
    cf.sigma = runif(1, 0.1, 2)
  )
}

ni.fit <- 33000    # iterations / chain
nb.fit <- 3000     # burn in periods for each chain
nt.fit <- 30       # thinning interval : (33,000 - 3,000) / 30 = 1000 posterior estimates / chain
nc.fit <- 3        # number of chains, total of 3,000 estimates for each model. 


parameters.lactin2 <- c("cf.a", "cf.b", "cf.tmax", "cf.delta_t", "cf.sigma", "r.pred") # repeated here

Temp.xs <- seq(-5, 45, 0.1) # Temperature gradient we're interested in - upped the granularity here
N.Temp.xs <-length(Temp.xs)

trait <- df.i$Growth.rate     # format the data for jags
N.obs <- length(trait)
temp <- df.i$Temperature

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

lac_jag <- jags(
  data = jag.data, 
  inits = inits.lactin.thomas, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin2_thomas.txt",
  n.thin = nt.fit, 
  n.chains = nc.fit, 
  n.burnin = nb.fit, 
  n.iter = ni.fit, 
  DIC = TRUE, 
  working.directory = getwd()
)

############# Fit TPCs and save the summary statistics #######################