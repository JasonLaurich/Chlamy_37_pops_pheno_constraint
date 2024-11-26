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
library(rTPC)
library(MuMIn)

##################### Upload and examine data #######################

df <- read.csv("data-processed/05_final_r_estimates.csv")

head(df)
str(df)

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$population.number)
df$well.ID <- as.factor(df$well.ID)
df$Temp <- df$temperature

mat <- split(df, df$pop.num)  # Matrixify the data!

############# Explore fitting TPCs - 1 population ###############################

i <- sample(1:38, 1) # We'll start with one population

df.i <- subset(mat[[i]])
df.i <- droplevels(df.i) 

model.AICc <- data.frame(
  model = character(),
  AICc = numeric(),
  stringsAsFactors = FALSE
)

nls.plot.list <- list()

# We'll start fitting models with nls_multstart(), selecting 4 from the rTPC package and 3 from the Kontopolous paper

# modified Deutsch

start.vals.deut <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'deutsch_2008')

deut_nls <- nls_multstart(r.exp ~ deutsch_2008(temp = Temp, rmax, topt, ctmax, a),
                          data = df.i,
                          iter = c(4, 4, 4, 4), 
                          start_lower = start.vals.deut - 10,
                          start_upper = start.vals.deut + 10,
                          lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'deutsch_2008'),
                          upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'deutsch_2008'),
                          supp_errors = 'Y',
                          convergence_count = FALSE
)

summary(deut_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Deutsch", AICc = AICc(deut_nls)))

preds.deut <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.deut <- broom::augment(deut_nls, newdata = preds.deut)

deut_plot <- ggplot(preds.deut) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Deutsch') + ylim(-3,5)

nls.plot.list[['deutsch']] <- deut_plot # store the plot

# Lactin 2

start.vals.lac <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'lactin2_1995')

lac_nls <- nls_multstart(r.exp ~ lactin2_1995(temp = Temp, a, b, tmax, delta_t),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.lac - 10,
                         start_upper = start.vals.lac + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'lactin2_1995'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'lactin2_1995'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(lac_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Lactin2", AICc = AICc(lac_nls)))

preds.lac <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.lac <- broom::augment(lac_nls, newdata = preds.lac)

lac_plot <- ggplot(preds.lac) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Lactin 2') +ylim(-3,5)

nls.plot.list[['Lactin 2']] <- lac_plot # store the plot

# Ratkowsky (bounded)

start.vals.rat <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'ratkowsky_1983')

rat_nls <- nls_multstart(r.exp ~ ratkowsky_1983(temp = Temp, tmin, tmax, a, b),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.rat - 10,
                         start_upper = start.vals.rat + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'ratkowsky_1983'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'ratkowsky_1983'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(rat_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Ratkowsky", AICc = AICc(rat_nls)))

preds.rat <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.rat <- broom::augment(rat_nls, newdata = preds.rat)

rat_plot <- ggplot(preds.rat) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Ratkowsky') +ylim(-3,5)

nls.plot.list[['Ratkowsky unbounded']] <- rat_plot # store the plot

bounded_ratkowsky <- function(temp, tmin, tmax, a, b) { # We're going to try to write a bounded function here.
  # Original Ratkowsky equation
  rate <- (a * (temp - tmin))^2 * (1 - exp(b * (temp - tmax)))^2
  
  # Bound the rate to 0 for temp > tmax
  rate[temp > tmax] <- 0
  
  # Ensure no negative rates
  rate[rate < 0] <- 0
  
  return(rate)
}

ratbd_nls <- nls_multstart(r.exp ~ bounded_ratkowsky(temp = Temp, tmin, tmax, a, b),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.rat - 10,
                         start_upper = start.vals.rat + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'ratkowsky_1983'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'ratkowsky_1983'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

model.AICc <- rbind(model.AICc, data.frame(model = "Ratkowsky bounded", AICc = AICc(ratbd_nls)))

preds.ratbd <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.ratbd <- broom::augment(ratbd_nls, newdata = preds.ratbd)

ratbd_plot <- ggplot(preds.ratbd) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Ratkowsky bounded') +ylim(-3,5)

nls.plot.list[['Ratkowsky bounded']] <- ratbd_plot # store the plot

# Rezende

start.vals.rez <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'rezende_2019')

rez_nls <- nls_multstart(r.exp ~ rezende_2019(temp = Temp, q10, a, b, c),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.rez - 10,
                         start_upper = start.vals.rez + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'rezende_2019'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'rezende_2019'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(rez_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Rezende", AICc = AICc(rez_nls)))

preds.rez <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.rez <- broom::augment(rez_nls, newdata = preds.rez)

rez_plot <- ggplot(preds.rez) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Rezende') +ylim(-3,5)

nls.plot.list[['Rezende']] <- rez_plot # store the plot

# Ashrafi II
# We are using the model formulation from the supplement to Kontopolous 2024
# 6. Ashrafi II (3 parameters)^4 : B(T) = a + b x T^3/2 + c x T^2

# This isn't in rTPC, so we'll have to write the function ourselves:

ashrafi_simple <- function(temp, a, b, c) {
  a + b * temp^(3/2) + c * temp^2
}

ash_nls <- nls_multstart(
  r.exp ~ ashrafi_simple(Temp, a, b, c),
  data = df.i,
  iter = c(4, 4, 4),  # Number of iterations for starting values
  start_lower = c(a = -10, b = -10, c = -1),  # Lower bounds for start values
  start_upper = c(a = 10, b = 10, c = 1),    # Upper bounds for start values
  lower = c(a = -30, b = -30, c = -15),   # Hard lower parameter bounds
  upper = c(a = 30, b = 30, c = 15),      # Hard upper parameter bounds
  supp_errors = 'Y',
  convergence_count = FALSE
)

summary(ash_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Ashrafi II", AICc = AICc(ash_nls)))

preds.ash <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.ash <- broom::augment(ash_nls, newdata = preds.ash)

ash_plot <- ggplot(preds.ash) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Ashrafi II simple') +ylim(-3,5)

nls.plot.list[['Ashrafi II simple']] <- ash_plot # store the plot

# Atkin

B(T) = B0 x (a - b x T)^T/10

# Mitchell-Angilletta



