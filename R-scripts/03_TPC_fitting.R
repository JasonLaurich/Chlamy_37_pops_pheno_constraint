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

# We'll start with nls.multstart, and do them all.

# I'm going to create a dataframe to store model AICc scores, and a list for the plots
model.AICc <- data.frame(
  model = character(),
  AICc = numeric(),
  stringsAsFactors = FALSE
)

nls.plot.list <- list()

# Beta

start.vals.beta <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'beta_2012')

beta_nls <- nls_multstart(r.gt~beta_2012(temp = Temp, a, b, c, d, e),
                         data = df.i,
                         iter = c(4,4,4,4,4),
                         start_lower = start.vals.beta - 10,
                         start_upper = start.vals.beta + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'beta_2012'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'beta_2012'),
                         supp_errors = 'Y',
                         convergence_count = F)

model.AICc <- rbind(model.AICc, data.frame(model = "Beta", AICc = AICc(beta_nls)))

preds.beta <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.beta <- broom::augment(beta_nls, newdata = preds.beta)

beta_plot <- ggplot(preds.beta) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Beta')

nls.plot.list[['beta']] <- beta_plot # store the plot

summary(beta_nls)

# Boatman

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

model.AICc <- rbind(model.AICc, data.frame(model = "Boatman", AICc = AICc(boat_nls)))

summary(boat_nls)

preds.boat <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.boat <- broom::augment(boat_nls, newdata = preds.boat)

boat_plot <- ggplot(preds.boat) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() + ggtitle('Boatman')

Topt.boat <- preds.boat$Temp[which.max(preds.boat$.fitted)]
Tbr.boat <- diff(range(preds.boat$Temp[preds.boat$.fitted >= rmax / 2], na.rm = TRUE)) # Thermal breadth

nls.plot.list[['boatman']] <- boat_plot # store the plot

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

model.AICc <- rbind(model.AICc, data.frame(model = "Briere", AICc = AICc(bri_nls)))

preds.bri <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.bri <- broom::augment(bri_nls, newdata = preds.bri)

bri_plot <- ggplot(preds.bri) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() + ggtitle('Briere')

nls.plot.list[['briere']] <- bri_plot # store the plot

summary(bri_nls)

# Let's try calculating other parameters from this fitted curve. 
# For now, I'll just use this predicted data curve and extract estimates from there?

Topt <- preds$Temp[which.max(preds$.fitted)]
rmax <- max(preds$.fitted, na.rm = T)
Tbr <- diff(range(preds$Temp[preds$.fitted >= rmax / 2], na.rm = TRUE)) # Thermal breadth
# Ea seems harder to do, will look into later...

# Delong

start.vals.delong <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'delong_2017')
start.vals.delong <- c(c = 0.1, eb = 1, ef = 1, tm = 30, ehc = 1)  # Above didn't work, trying this

del_nls <- nls_multstart(r.gt~delong_2017(temp = Temp, c, eb, ef, tm, ehc),
                         data = df.i,
                         iter = c(4,4,4,4,4),
                         start_lower = start.vals.delong - 10,
                         start_upper = start.vals.delong + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'delong_2017'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'delong_2017'),
                         supp_errors = 'Y',
                         convergence_count = F)

model.AICc <- rbind(model.AICc, data.frame(model = "Delong", AICc = AICc(del_nls)))

preds.del <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.del <- broom::augment(del_nls, newdata = preds.del)

del_plot <- ggplot(preds.del) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() + ggtitle('Delong')

nls.plot.list[['delong']] <- del_plot # store the plot

summary(del_nls)

# Deutsch

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

model.AICc <- rbind(model.AICc, data.frame(model = "Deutsch", AICc = AICc(deut_nls)))

preds.deut <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.deut <- broom::augment(deut_nls, newdata = preds.deut)

deut_plot <- ggplot(preds.deut) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Deutsch') + ylim(-3,5)

nls.plot.list[['deutsch']] <- deut_plot # store the plot

# Flinn

start.vals.flinn <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'flinn_1991')

flinn_nls <- nls_multstart(r.gt ~ flinn_1991(temp = Temp, a, b, c),
                          data = df.i,
                          iter = c(4, 4, 4), 
                          start_lower = start.vals.flinn - 10,
                          start_upper = start.vals.flinn + 10,
                          lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'flinn_1991'),
                          upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'flinn_1991'),
                          supp_errors = 'Y',
                          convergence_count = FALSE
)

summary(flinn_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Flinn", AICc = AICc(flinn_nls)))

preds.flinn <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.flinn <- broom::augment(flinn_nls, newdata = preds.flinn)

flinn_plot <- ggplot(preds.flinn) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Flinn')

nls.plot.list[['flinn']] <- flinn_plot # store the plot

# Gaussian

start.vals.gaus <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'gaussian_1987')

gaus_nls <- nls_multstart(r.gt ~ gaussian_1987(temp = Temp, rmax, topt, a),
                           data = df.i,
                           iter = c(4, 4, 4), 
                           start_lower = start.vals.gaus - 10,
                           start_upper = start.vals.gaus + 10,
                           lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'gaussian_1987'),
                           upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'gaussian_1987'),
                           supp_errors = 'Y',
                           convergence_count = FALSE
)

summary(gaus_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Gaussian", AICc = AICc(gaus_nls)))

preds.gaus <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.gaus <- broom::augment(gaus_nls, newdata = preds.gaus)

gaus_plot <- ggplot(preds.gaus) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Gaussian')

nls.plot.list[['gaussian']] <- gaus_plot # store the plot

# Hinshelwood

start.vals.hin <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'hinshelwood_1947')

hin_nls <- nls_multstart(r.gt ~ hinshelwood_1947(temp = Temp, a, e, b, eh),
                          data = df.i,
                          iter = c(4, 4, 4, 4), 
                          start_lower = start.vals.hin - 10,
                          start_upper = start.vals.hin + 10,
                          lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'hinshelwood_1947'),
                          upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'hinshelwood_1947'),
                          supp_errors = 'Y',
                          convergence_count = FALSE
)

summary(hin_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Hinshelwood", AICc = AICc(hin_nls)))

preds.hin <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.hin <- broom::augment(hin_nls, newdata = preds.hin)

hin_plot <- ggplot(preds.hin) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Hinshelwood')

nls.plot.list[['hinshelwood']] <- hin_plot # store the plot

# Joehnk

start.vals.joe <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'joehnk_2008')


joe_nls <- nls_multstart(r.gt ~ joehnk_2008(temp = Temp, rmax, topt, a, b, c),
                         data = df.i,
                         iter = c(4, 4, 4, 4, 4), 
                         start_lower = start.vals.joe - 10,
                         start_upper = start.vals.joe + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'joehnk_2008'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'joehnk_2008'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(joe_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Joehnk", AICc = AICc(joe_nls)))

preds.joe <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.joe <- broom::augment(joe_nls, newdata = preds.joe)

joe_plot <- ggplot(preds.joe) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Joenhk') +ylim(-3,5)

nls.plot.list[['Joenhk']] <- joe_plot # store the plot

# Johnson-Lewis

start.vals.john <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'johnsonlewin_1946')
start.vals.john[is.na(start.vals.john)] <- 0
upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'johnsonlewin_1946') # Need to do manually, r0 infinite
upper[upper == Inf] <- 10

john_nls <- nls_multstart(r.gt ~ johnsonlewin_1946(temp = Temp, r0, e, eh, topt),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.john - 10,
                         start_upper = start.vals.john + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'johnsonlewin_1946'),
                         upper = upper,
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(john_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Johnson-Lewis", AICc = AICc(john_nls)))

preds.john <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.john <- broom::augment(john_nls, newdata = preds.john)

john_plot <- ggplot(preds.john) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Johnson-Lewis')

nls.plot.list[['Johnson-Lewis']] <- john_plot # store the plot

# Kamykowski

start.vals.kam <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'kamykowski_1985')

kam_nls <- nls_multstart(r.gt ~ kamykowski_1985(temp = Temp, tmin, tmax, a, b, c),
                         data = df.i,
                         iter = c(4, 4, 4, 4, 4), 
                         start_lower = start.vals.kam - 10,
                         start_upper = start.vals.kam + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'kamykowski_1985'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'kamykowski_1985'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(kam_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Kamykowski", AICc = AICc(kam_nls)))

preds.kam <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.kam <- broom::augment(kam_nls, newdata = preds.kam)

kam_plot <- ggplot(preds.kam) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Kamykowski') +ylim (-3,5)

nls.plot.list[['Kamykowski']] <- kam_plot # store the plot

# Lactin 2

start.vals.lac <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'lactin2_1995')

lac_nls <- nls_multstart(r.gt ~ lactin2_1995(temp = Temp, a, b, tmax, delta_t),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.lac - 10,
                         start_upper = start.vals.lac + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'lactin2_1995'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'lactin2_1995'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(lac_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Lactin2", AICc = AICc(lac_nls)))

preds.lac <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.lac <- broom::augment(lac_nls, newdata = preds.lac)

lac_plot <- ggplot(preds.lac) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Lactin 2') +ylim(-3,5)

nls.plot.list[['Lactin 2']] <- lac_plot # store the plot

# LRF

start.vals.lrf <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'lrf_1991')

lrf_nls <- nls_multstart(r.gt ~ lrf_1991(temp = Temp, rmax, topt, tmin, tmax),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.lrf - 10,
                         start_upper = start.vals.lrf + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'lrf_1991'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'lrf_1991'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(lrf_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "LRF", AICc = AICc(lrf_nls)))

preds.lrf <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.lrf <- broom::augment(lrf_nls, newdata = preds.lrf)

lrf_plot <- ggplot(preds.lrf) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('LRF') +ylim(-3,5)

nls.plot.list[['LRF']] <- lrf_plot # store the plot

# Modified Gaussian

start.vals.mod <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'modifiedgaussian_2006')

mod_nls <- nls_multstart(r.gt ~ modifiedgaussian_2006(temp = Temp, rmax, topt, a, b),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.mod - 10,
                         start_upper = start.vals.mod + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'modifiedgaussian_2006'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'modifiedgaussian_2006'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(mod_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "ModifiedGaussian", AICc = AICc(mod_nls)))

preds.mod <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.mod <- broom::augment(mod_nls, newdata = preds.mod)

mod_plot <- ggplot(preds.mod) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('ModifiedGaussian') +ylim(-3,5)

nls.plot.list[['ModifiedGaussian']] <- mod_plot # store the plot

# Oneill

start.vals.onl <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'oneill_1972')

onl_nls <- nls_multstart(r.gt ~ oneill_1972(temp = Temp, rmax, ctmax, topt, q10),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.onl - 10,
                         start_upper = start.vals.onl + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'oneill_1972'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'oneill_1972'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(onl_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "ONeill", AICc = AICc(onl_nls)))

preds.onl <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.onl <- broom::augment(onl_nls, newdata = preds.onl)

onl_plot <- ggplot(preds.onl) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('ONeill') +ylim(-3,5)

nls.plot.list[['ONeill']] <- onl_plot # store the plot

# Pawar

start.vals.paw <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'pawar_2018')
start.vals.paw <- c(r_tref = 1.5, e = 0.7, eh = 1.5, topt = 34)  # Manually set, troublesome NA

tref <- 15

paw_nls <- nls_multstart(r.gt ~ pawar_2018(temp = Temp, r_tref, e, eh, topt, tref = tref),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = c(r_tref = 1.5, e = 0.6, eh = 1.0, topt = 30),  # Example bounds
                         start_upper = c(r_tref = 3.0, e = 1.0, eh = 2.5, topt = 40),
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'pawar_2018'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'pawar_2018'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

# I'm skipping Pawar for now, having trouble fitting it.

# Quadratic

start.vals.quad <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'quadratic_2008')

quad_nls <- nls_multstart(r.gt ~ quadratic_2008(temp = Temp, a, b, c),
                         data = df.i,
                         iter = c(4, 4, 4), 
                         start_lower = start.vals.quad - 10,
                         start_upper = start.vals.quad + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'quadratic_2008'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'quadratic_2008'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(quad_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Quadratic", AICc = AICc(quad_nls)))

preds.quad <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.quad <- broom::augment(quad_nls, newdata = preds.quad)

quad_plot <- ggplot(preds.quad) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Quadratic') +ylim(-3,5)

nls.plot.list[['Quadratic']] <- quad_plot # store the plot

# Ratkowsky

start.vals.rat <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'ratkowsky_1983')

rat_nls <- nls_multstart(r.gt ~ ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                          data = df.i,
                          iter = c(4, 4, 4, 4), 
                          start_lower = start.vals.rat - 10,
                          start_upper = start.vals.rat + 10,
                          lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'ratkowsky_1983'),
                          upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'ratkowsky_1983'),
                          supp_errors = 'Y',
                          convergence_count = FALSE
)

summary(rat_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Ratkowsky", AICc = AICc(rat_nls)))

preds.rat <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.rat <- broom::augment(rat_nls, newdata = preds.rat)

rat_plot <- ggplot(preds.rat) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Ratkowsky') +ylim(-3,5)

nls.plot.list[['Ratkowsky']] <- rat_plot # store the plot

# Rezende

start.vals.rez <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'rezende_2019')

rez_nls <- nls_multstart(r.gt ~ rezende_2019(temp = Temp, q10, a, b, c),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.rez - 10,
                         start_upper = start.vals.rez + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'rezende_2019'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'rezende_2019'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(rez_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Rezende", AICc = AICc(rez_nls)))

preds.rez <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.rez <- broom::augment(rez_nls, newdata = preds.rez)

rez_plot <- ggplot(preds.rez) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Rezende') +ylim(-3,5)

nls.plot.list[['Rezende']] <- rez_plot # store the plot

# Sharpe-Schoolfield Full

start.vals.ss <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'sharpeschoolfull_1981')
start.vals.ss <- c(r_tref = 1.5, e = 0.7, el = 0.1, tl =10, eh =0.1, th = 37)  # Manually set, troublesome NA


ss_nls <- nls_multstart(r.gt ~ sharpeschoolfull_1981(temp = Temp, r_tref, e, el, tl, eh, th, tref = 20),
                         data = df.i,
                         iter = c(4, 4, 4, 4, 4, 4), 
                         start_lower = start.vals.ss - 10,
                         start_upper = start.vals.ss + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'sharpeschoolfull_1981'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'sharpeschoolfull_1981'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(ss_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Sharpe-Schoolfield", AICc = AICc(ss_nls)))

preds.ss <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.ss <- broom::augment(ss_nls, newdata = preds.ss)

ss_plot <- ggplot(preds.ss) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Sharpe-Schoolfield') +ylim(-3,5)

nls.plot.list[['Sharpe-Schoolfield']] <- ss_plot # store the plot

# Spain

start.vals.spn <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'spain_1982')

spn_nls <- nls_multstart(r.gt ~ spain_1982(temp = Temp, a, b, c, r0),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.spn - 10,
                         start_upper = start.vals.spn + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'spain_1982'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'spain_1982'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(spn_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Spain", AICc = AICc(spn_nls)))

preds.spn <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.spn <- broom::augment(spn_nls, newdata = preds.spn)

spn_plot <- ggplot(preds.spn) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Spain') +ylim(-3,5)

nls.plot.list[['Spain']] <- spn_plot # store the plot

# Thomas 2012

start.vals.th1 <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'thomas_2012')

th1_nls <- nls_multstart(r.gt ~ thomas_2012(temp = Temp, a, b, c, topt),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.th1 - 10,
                         start_upper = start.vals.th1 + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'thomas_2012'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'thomas_2012'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(th1_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Thomas2012", AICc = AICc(th1_nls)))

preds.th1 <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.th1 <- broom::augment(th1_nls, newdata = preds.th1)

th1_plot <- ggplot(preds.th1) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Thomas2012') +ylim(-3,5)

nls.plot.list[['Thomas2012']] <- th1_plot # store the plot

# Thomas 2017

start.vals.th2 <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'thomas_2017')

th2_nls <- nls_multstart(r.gt ~ thomas_2017(temp = Temp, a, b, c, d, e),
                         data = df.i,
                         iter = c(4, 4, 4, 4, 4), 
                         start_lower = start.vals.th2 - 10,
                         start_upper = start.vals.th2 + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'thomas_2017'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'thomas_2017'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

# Not sure if this model is appropriate - I don't have mortality data!

# Weibull, last one

start.vals.wei <- get_start_vals(df.i$Temp, df.i$r.gt, model_name = 'weibull_1995')

wei_nls <- nls_multstart(r.gt ~ weibull_1995(temp = Temp, a, topt, b, c),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.wei - 10,
                         start_upper = start.vals.wei + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.gt, model_name = 'weibull_1995'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.gt, model_name = 'weibull_1995'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(wei_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Weibull", AICc = AICc(wei_nls)))

preds.wei <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.wei <- broom::augment(wei_nls, newdata = preds.wei)

wei_plot <- ggplot(preds.wei) + geom_point(aes(Temp, r.gt), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Weibull') +ylim(-3,5)

nls.plot.list[['Weibull']] <- wei_plot # store the plot

write.csv(model.AICc, "data-processed/04_rTPC_modelcomparison_population7.csv") # Save r estimates and other info

mod_grid<- plot_grid(plotlist = nls.plot.list, ncol = 6)
ggsave("figures/03_rTPC_modelcomparison_pop7.pdf", plot = mod_grid, width = 24, height = 16)

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


############# Figure generation and file export #######################
