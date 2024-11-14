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
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits

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

# OK so no we are going to try fitting these same 3 models using R2jags, which will fit Bayesian TPCs
# This model framework (see https://github.com/JoeyBernhardt/anopheles-rate-summation/blob/master/AnalysisDemo.R)
# is predicated on text files (e.g. briere.txt) that specify model structure and parameters. New ones must be
# made and saved in the directory before calling them when fitting the TPCs themselves. 

# Let's set up our MCMC model settings
ni <- 55000 # iterations / chain
nb <- 5000 # burn in periods for each chain
nt <- 50 # thinning interval : (55,000 - 5,000) / 50 = 1000 posterior estimates / chain
nc <- 5 # number of chains

# Temperature gradient across which to measure r
Temp.xs <- seq(0, 45, 0.1) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

# identify traits for jags
trait <- df.i$r.gt
N.obs <- length(trait)
temp <- df.i$Temp

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Briere settings - initial values tighter than priors, incorporate variation across chains.

# inits function
inits.bri<-function(){list(
  cf.q = runif( 1, 0.01, 0.1),
  cf.Tm = runif(1, 35, 45),
  cf.T0 = runif(1, 0, 10),
  cf.sigma = rlnorm(1, log(1), 0.5)
  )
}

parameters.bri <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "r.pred") # estimate these, will depend based on model

# jags MCMC, Briere function
bri_jag <- jags(data=jag.data, inits=inits.bri, parameters.to.save=parameters.bri, model.file="briere.txt",
                                 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

bri_jag$BUGSoutput$summary[1:5,] # Get estimates
mcmcplot(bri_jag) # Evaluate model performance
bri_jag$BUGSoutput$DIC # DIC

# plot!
plot(trait ~ jitter(Temp, 0.5), xlim = c(0, 45), data = df.i, 
     ylab = "Growth rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(bri_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bri_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bri_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Let's try the Boatman model

# inits function
inits.boat<-function(){list(
  rmax = runif(1, 1.5, 6.5),
  tmin = runif(1, 0, 10),
  tmax = runif(1, 35, 45),
  a = runif( 1, 0.01, 0.1),
  b = runif( 1, 0.01, 0.1),
  sigma = rlnorm(1, log(1), 0.5)
)
}

parameters.boat <- c("rmax", "tmin", "tmax", "a", "b", "sigma", "r.pred") # estimate these, will depend based on model

# jags MCMC, Boatman function
boat_jag <- jags(data=jag.data, inits=inits.boat, parameters.to.save=parameters.boat, model.file="boatman.txt",
                n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

boat_jag$BUGSoutput$summary[1:5,] # Get estimates
mcmcplot(boat_jag) # Evaluate model performance
boat_jag$BUGSoutput$DIC # DIC

# plot!
plot(trait ~ jitter(Temp, 0.5), xlim = c(0, 45), data = df.i, 
     ylab = "Growth rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(boat_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(boat_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(boat_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# That is not working great! We'll revisit this, but let's try the Deutsch for now:

# Initial values for the Deutsch model
inits.deut <- function() {
  list(
    rmax = runif(1, 1, 10),        # Random initial guess for rmax within prior range
    topt = runif(1, 20, 35),       # Initial guess for topt within a plausible range
    ctmax = runif(1, 36, 50),      # Initial guess for ctmax within a plausible range
    a = runif(1, 0.1, 5),          # Initial guess for a
    sigma = rlnorm(1, log(1), 0.5) # Initial guess for sigma using log-normal distribution
  )
}

# Parameters to save for Deutsch model
parameters.deut <- c("rmax", "topt", "ctmax", "a", "sigma", "r.pred")

# jags MCMC, Deutsch function
deut_jag <- jags(data=jag.data, inits=inits.deut, parameters.to.save=parameters.deut, model.file="deutsch.txt",
                 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

deut_jag$BUGSoutput$summary[1:5,] # Get estimates
mcmcplot(deut_jag) # Evaluate model performance
deut_jag$BUGSoutput$DIC # DIC

# plot!
plot(trait ~ jitter(Temp, 0.5), xlim = c(0, 45), data = df.i, 
     ylab = "Growth rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(deut_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(deut_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(deut_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# Also doesn't work well, possibly due to conditionality, and Topt, Tmax overlapping.

inits.sharpe <- function() {
  list(
    r_tref = runif(1, 1, 5),           # Initial rate at standard reference temperature
    e = runif(1, 0.5, 1.5),            # Activation energy (eV)
    el = runif(1, 1, 4),               # Low-temperature deactivation energy
    tl = runif(1, 10, 30),             # Temperature for low-temp 1/2 activity suppression
    eh = runif(1, 1, 4),               # High-temperature deactivation energy
    th = runif(1, 30, 50),             # Temperature for high-temp 1/2 activity suppression
    sigma = rlnorm(1, log(1), 0.5)     # Initial guess for sigma
  )
}

parameters.sharpe <- c("r_tref", "e", "el", "tl", "eh", "th", "sigma", "r.pred")

tref <- 20  # define tref

jag.data2 <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs, tref = tref)


sharpe_jag <- jags(
  data = jag.data2,
  inits = inits.sharpe,
  parameters.to.save = parameters.sharpe,
  model.file = "sharpe.txt",
  n.thin = nt,
  n.chains = nc,
  n.burnin = nb,
  n.iter = ni,
  DIC = TRUE,
  working.directory = getwd()
)

sharpe_jag$BUGSoutput$summary[1:5,] # Get estimates
mcmcplot(sharpe_jag) # Evaluate model performance
sharpe_jag$BUGSoutput$DIC # DIC

# plot!
plot(trait ~ jitter(Temp, 0.5), xlim = c(0, 45), data = df.i, 
     ylab = "Growth rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(sharpe_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(sharpe_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(sharpe_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

# OK, we're going to try the quadratic function

# Initial values for the Quadratic model
inits.quad <- function() {
  list(
    cf.q = runif(1, 0.01, 0.1),    # Coefficient for the quadratic term
    cf.T0 = runif(1, 5, 15),       # Lower temperature threshold
    cf.Tm = runif(1, 35, 45),      # Upper temperature threshold
    cf.sigma = rlnorm(1, log(1), 0.5) # Variability in data around the mean
  )
}

parameters.quad <- c("cf.q", "cf.T0", "cf.Tm", "cf.sigma", "r.pred")

quad_jag <- jags(
  data = jag.data,                # Data bundled for JAGS
  inits = inits.quad,         # Initial values function
  parameters.to.save = parameters.quad, # Parameters to monitor
  model.file = "quad.txt",    # JAGS model file
  n.thin = nt,
  n.chains = nc,
  n.burnin = nb,
  n.iter = ni,
  DIC = TRUE,
  working.directory = getwd()
)

quad_jag$BUGSoutput$summary[1:5,] # Get estimates
mcmcplot(quad_jag) # Evaluate model performance
quad_jag$BUGSoutput$DIC # DIC

# plot!
plot(trait ~ jitter(Temp, 0.5), xlim = c(0, 45), data = df.i, 
     ylab = "Growth rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(quad_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(quad_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(quad_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

############# Loop through entire dataset ################


############# Figure generation and file export ########################