# Jason R Laurich
# November 12, 2024

# This script will upload previously tabulated measurements of r for our 37 populations of Chlamydomonas at
# various temperatures, and fit TPCs to the data, using both Bayesian and nls methods.

############# Packages ########################

library(nls.multstart)
library(rTPC)
library(tidyr)
library(cowplot)
library(ggplot2)
library(dplyr)
library(MuMIn) #AICc

############# Upload and examine data #######################

# Both dataframes contain estimates of r for each well (within population and temperature)
df.exp <- read.csv("data-processed/01_pop_well_r_estimates.csv") # This is the exploratory dataset, where I calculated r using a variety of methods

df.gt <- read.csv("data-processed/03_r_growthtools.csv") # This dataset contains the r estimates fit with growthTools
# Note that data for 40C treatment here is worth considering further. We did not exclude the aberrant early time points
# at 40C which likely reflect a brief burst of reproduction in non-temperature acclimated Chlamydomonas cells.
# Using the r.exp.ln column of the df.exp for 40 C would give us r fit on post-burst data for 40 C.

str(df.gt)
df.gt$pop.fac <- as.factor(df.gt$population)
df.gt$pop.num <- as.numeric(df.gt$population.number)
df.gt$well.ID <- as.factor(df.gt$well.ID)
df.gt$Temp <- df.gt$temperature

mat.gt <- split(df.gt, df.gt$pop.num)  # Let's make this a matrix.

############# Explore fitting TPCs ###############################

# Need to think about which models to fit. I think I will start with (1) Briere - nice and simple, (2) Boatman, and 
# (3) Deutsch - these latter 2 seem like they might more realistically fit tmin and tmax, without requiring estimating
# the massive amount of parameters of say a Sharpe-Schoolfield model.

# Let's lay out the basic structure we will be using to loop later

for(i in 1:38){
  
  df.i <- subset(mat.gt[[i]])
  df.i <- droplevels(df.i) 
  
  # Fit model
  
  # Extract key model estimates

}

# Let's run some models on a random test population

i <- sample(1:38, 1)
df.i <- subset(mat.gt[[i]])
df.i <- droplevels(df.i) 

# We'll start with nls.multstart

# Briere

start.vals.briere <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'briere2_1999')

bri_nls <- nls_multstart(r.gt~briere2_1999(temp = Temp, Tmin, Tmax, a, b),
                        data = df.i,
                        iter = c(4,4,4,4),
                        start_lower = start.vals.briere - 10,
                        start_upper = start.vals.briere + 10,
                        lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'briere2_1999'),
                        upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'briere2_1999'),
                        supp_errors = 'Y',
                        convergence_count = F)

#Let's try plotting this
preds <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds <- broom::augment(bri_nls, newdata = preds)

ggplot(preds) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic()

summary(bri_nls)
AICc(bri_nls)

# Let's try calculating other parameters from this fitted curve. 
# For now, I'll just use this predicted data curve and extract estimates from there?

Topt <- preds$Temp[which.max(preds$.fitted)]
rmax <- max(preds$.fitted, na.rm = T)
Tbr <- diff(range(preds$Temp[preds$.fitted >= rmax / 2], na.rm = TRUE)) # Thermal breadth
# Ea seems harder to do, will look into later...

# Let's try the Boatman one as well.

start.vals.boatman <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'boatman_2017')

boat_nls <- nls_multstart(r.gt~boatman_2017(temp = Temp, rmax, tmin, tmax, a, b),
                                    data = df.i,
                                    iter = c(4,4,4,4,4),
                                    start_lower = start.vals.boatman - 10,
                                    start_upper = start.vals.boatman + 10,
                                    lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'boatman_2017'),
                                    upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'boatman_2017'),
                                    supp_errors = 'Y',
                                    convergence_count = FALSE)

summary(boat_nls)
AICc(boat_nls)

preds.boat <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.boat <- broom::augment(boat_nls, newdata = preds.boat)

ggplot(preds.boat) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic()

Topt.boat <- preds.boat$Temp[which.max(preds.boat$.fitted)]
Tbr.boat <- diff(range(preds.boat$Temp[preds.boat$.fitted >= rmax / 2], na.rm = TRUE)) # Thermal breadth

# Let's finish with the Deutsch model for now.

start.vals.deut <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'deutsch_2008')

deut_nls <- nls_multstart(r.gt ~ deutsch_2008(temp = Temp, rmax, topt, ctmax, a),
  data = df.i,
  iter = c(4, 4, 4, 4), 
  start_lower = start.vals.deut - 10,
  start_upper = start.vals.deut + 10,
  lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'deutsch_2008'),
  upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'deutsch_2008'),
  supp_errors = 'Y',
  convergence_count = FALSE
)

summary(deut_nls)
AICc(deut_nls)

preds.deut <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.deut <- broom::augment(deut_nls, newdata = preds.deut)

ggplot(preds.deut) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic()

############# Loop through entire dataset ################


############# Figure generation and file export ########################