# Jason Laurich
# Oct 21, 2024

# This script will fits TPCs to data from Joey Bernhardt's postdoctoral research using the
# bayesTPC package.

# Uploading packages

library(nimble)
library(HDInterval)
library(MCMCvis)
library(coda) # makes diagnostic plots
library(IDPmisc) # makes nice colored pairs plots to look at joint posteriors
library(matrixStats)
library(truncnorm)
#devtools::install_github("johnwilliamsmithjr/bayesTPC") # If needed
library(bayesTPC)
library(cowplot)
library(tidyverse)

# Importing, interrogating, and re-arranging the data

data <- read.csv("data-processed/globe-chlamy-exponential-RFU-time.csv")
head(data)
str(data)

length(unique(data$population))
# 15 populations, we'll want to fit separate TPCs to each one. I imagine this will be easiest if we turn this dataframe into a matrix

#Remove any missing data
data <- subset(data, !is.na(data$RFU))

#OK, let's create a matrix to hold the data for each population (1 to 15 here)
pops<-as.data.frame(1:15)
mat <- lapply(1:nrow(pops), function(z) data[which(data$population==pops[z,1]),])

#Check it
for (i in c(5,10)){
  print(head(mat[[i]]))
}
#Looks good!


# Let's look at a few plots!

#This works with the data in matrix form
ggplot(aes(x = days, y = RFU, colour = factor(temperature)), data=mat[[3]]) + geom_point()
ggplot(aes(x = days, y = RFU, colour = factor(temperature)), data=mat[[9]]) + geom_point()

#So does this
for (i in 1:2){
  print(ggplot(aes(x = days, y = RFU, colour = factor(temperature)), data=mat[[i]]) + geom_point())
}

#OK, the plots look fine!

# We need to calculate the growth rate here first. This will be on our y-axis. For each population at each T,
# We need a measurement for replicate (ie. we are collapsing the temporal axis)
# I imagine we take a measurement of growth during the exponential phase here?

#OK, I want to try fitting growth rate a couple of ways

#1. nls

#2. growthrate package

#3. brms (bayesian)

#Let's consider which model of TPC to fit
get_models()

