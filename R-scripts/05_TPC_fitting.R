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
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)

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

# Ratkowsky, bounded so it doesn't jump back up after Tmax

start.vals.rat <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'ratkowsky_1983')

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
# We are using the model formulation from the supplement to Kontopolous 2024 :
# Atkin (3 parameters)^5: B(T) = B0 x (a - b x T)^T/10

atkin <- function(temp, B0, a, b) {
  B0 * (a - b * temp)^(temp / 10)
}

atk_nls <- nls_multstart(
  r.exp ~ atkin(Temp, B0, a, b),
  data = df.i,
  iter = c(4, 4, 4),
  start_lower = c(B0 = 0.1, a = 10, b = 0.01),
  start_upper = c(B0 = 5, a = 30, b = 0.1),
  lower = c(B0 = 0, a = 0, b = 0),
  upper = c(B0 = 10, a = 50, b = 1),
  supp_errors = 'Y',
  convergence_count = FALSE
)

summary(atk_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Atkin", AICc = AICc(atk_nls)))

preds.atk <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.atk <- broom::augment(atk_nls, newdata = preds.atk)

atk_plot <- ggplot(preds.atk) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Atkin') +ylim(-3,5)

nls.plot.list[['Atkin']] <- atk_plot # store the plot

# Mitchell-Angilletta
# # Mitchell-Angilletta (3 parameters)^6: B(T) = (a/(2 x b)) x [1 + cos (((T-T_pk)/b) x pi)]

mitchell <- function(temp, a, b, Tpk) {
  (a / (2 * b)) * (1 + cos(((temp - Tpk) / b) * pi))
}

mitch_nls <- nls_multstart(
  r.exp ~ mitchell(Temp, a, b, Tpk),
  data = df.i,
  iter = c(4, 4, 4),
  start_lower = c(a = 0.1, b = 1, Tpk = 20),
  start_upper = c(a = 10, b = 10, Tpk = 40),
  lower = c(a = 0, b = 0.1, Tpk = 0),
  upper = c(a = Inf, b = Inf, Tpk = 50),
  supp_errors = 'Y',
  convergence_count = FALSE
)

summary(mitch_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Mitchell-Angilletta", AICc = AICc(mitch_nls)))

preds.mitch <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.mitch <- broom::augment(mitch_nls, newdata = preds.mitch)

mitch_plot <- ggplot(preds.mitch) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Mitchell-Angilletta') +ylim(-3,5)

nls.plot.list[['Mitchell-Angilletta']] <- mitch_plot # store the plot

# I don't love how these are working. Let's try a few more

# Analytis-Kontomodimas

anal_kont <- function(temp, a, tmin, tmax) {
  a * (temp - tmin)^2 * (tmax - temp)
}

ankt_nls <- nls_multstart(
  r.exp ~ anal_kont(Temp, a, tmin, tmax),
  data = df.i,
  iter = c(4, 4, 4),
  start_lower = c(a = 0.1, tmin = 0, tmax = 30),
  start_upper = c(a = 10, tmin = 10, tmax = 50),
  lower = c(a = 0, tmin = -Inf, tmax = -Inf),
  upper = c(a = Inf, tmin = Inf, tmax = Inf),
  supp_errors = 'Y',
  convergence_count = FALSE
)

model.AICc <- rbind(model.AICc, data.frame(model = "Analytis-Kontodimas", AICc = AICc(ankt_nls)))

preds.ankt <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.ankt <- broom::augment(ankt_nls, newdata = preds.ankt)

ankt_plot <- ggplot(preds.ankt) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Analytis-Kontodimas') +ylim(-3,5)

nls.plot.list[['Analytis-Kontomodimas']] <- ankt_plot # store the plot

# Eubank

eubank <- function(temp, a, tpk, b) {
  a / ((temp - tpk)^2 + b)
}

eub_nls <- nls_multstart(
  r.exp ~ eubank(Temp, a, tpk, b),
  data = df.i,
  iter = c(4, 4, 4),
  start_lower = c(a = 0.1, tpk = 20, b = 0.1),
  start_upper = c(a = 10, tpk = 35, b = 5),
  lower = c(a = 0, tpk = 0, b = 0.01),
  upper = c(a = Inf, tpk = Inf, b = Inf),
  supp_errors = 'Y',
  convergence_count = FALSE
)

summary(eub_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Eubank", AICc = AICc(eub_nls)))

preds.eub <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.eub <- broom::augment(eub_nls, newdata = preds.eub)

eub_plot <- ggplot(preds.eub) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Eubank') +ylim(-3,6)

nls.plot.list[['Eubank']] <- eub_plot # store the plot

# Taylor-Sexton

tay_sex <- function(temp, bpk, tmin, tpk) {
  bpk * (-(temp - tmin)^4 + 2 * (temp - tmin)^2 * (tpk - tmin)^2) / (tpk - tmin)^4
}

tay_nls <- nls_multstart(
  r.exp ~ tay_sex(Temp, bpk, tmin, tpk),
  data = df.i,
  iter = c(4, 4, 4),
  start_lower = c(bpk = 0.1, tmin = 0, tpk = 20),
  start_upper = c(bpk = 10, tmin = 10, tpk = 40),
  lower = c(bpk = 0, tmin = -Inf, tpk = -Inf),
  upper = c(bpk = Inf, tmin = Inf, tpk = Inf),
  supp_errors = 'Y',
  convergence_count = FALSE
)

summary(tay_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Taylor-Sexton", AICc = AICc(tay_nls)))

preds.tay <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.tay <- broom::augment(tay_nls, newdata = preds.tay)

tay_plot <- ggplot(preds.tay) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Taylor-Sexton') +ylim(-3,5)

nls.plot.list[['Taylor-Sexton']] <- tay_plot # store the plot

# Briere II (rTPC)

start.vals.briere <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'briere2_1999')

bri_nls <- nls_multstart(r.exp~briere2_1999(temp = Temp, Tmin, Tmax, a, b),
                         data = df.i,
                         iter = c(4,4,4,4),
                         start_lower = c(Tmin = 5, Tmax = 30, a = 0, b = 0.1),
                         start_upper = c(Tmin = 15, Tmax = 50, a = 10, b = 5),
                         lower.lims <- c(a = 0, tmin = 0, tmax = 30, b = 0.1),
                         upper.lims <- c(a = 20, tmin = 20, tmax = 50, b = 10),
                         supp_errors = 'Y',
                         convergence_count = F)

summary(bri_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Briere", AICc = AICc(bri_nls)))

preds.bri <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.bri <- broom::augment(bri_nls, newdata = preds.bri)

bri_plot <- ggplot(preds.bri) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() + ggtitle('Briere')

nls.plot.list[['briere']] <- bri_plot # store the plot

# Thomas 1 (Norberg) - requested by Joey, this is in the rTPC package

start.vals.thom <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'thomas_2012')

thom_nls <- nls_multstart(r.exp~thomas_2012(temp = Temp, a, b, c, topt),
                         data = df.i,
                         iter = c(4,4,4,4),
                         start_lower = start.vals.thom - 10,
                         start_upper = start.vals.thom + 10,
                         lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'thomas_2012'),
                         upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'thomas_2012'),
                         supp_errors = 'Y',
                         convergence_count = FALSE)

summary(thom_nls)

model.AICc <- rbind(model.AICc, data.frame(model = "Thomas 1", AICc = AICc(thom_nls)))

preds.thom <- data.frame(Temp = seq(min(df.i$Temp - 2), max(df.i$Temp +2), length.out = 100))
preds.thom <- broom::augment(thom_nls, newdata = preds.thom)

thom_plot <- ggplot(preds.thom) + geom_point(aes(Temp, r.exp), df.i) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_classic() + ggtitle('Thomas 1')

nls.plot.list[['thomas']] <- thom_plot # store the plot

# Just for fun, let's fit the Deutsch model from scratch to make sure our approach to fitting our models is working equally well to rTPC

deutsch2 <- function(temp, rmax, topt, ctmax, a) {
  ifelse(
    temp < topt, 
    # Case when temp < topt
    rmax * exp(-((temp - topt) / (2 * a))^2),
    # Case when temp > topt
    rmax * (1 - ((temp - topt) / (topt - ctmax))^2)
  )
}

deut2_nls <- nls_multstart(
  r.exp ~ deutsch2(Temp, rmax, topt, ctmax, a),
  data = df.i,
  iter = c(4, 4, 4, 4),
  start_lower = c(rmax = 0.1, topt = 10, ctmax = 30, a = 0.1),
  start_upper = c(rmax = 10, topt = 30, ctmax = 50, a = 10),
  lower = c(rmax = 0, topt = -Inf, ctmax = -Inf, a = 0.01),
  upper = c(rmax = Inf, topt = Inf, ctmax = Inf, a = Inf),
  supp_errors = 'Y',
  convergence_count = FALSE
)

summary(deut2_nls)
AICc(deut2_nls) # So this is giving me the same numbers as the rTPC package which is great! This means that our custom models are working equally well.

write.csv(model.AICc, "data-processed/06_rTPC_modelcomparison_population2.csv") # Save model comparison data

mod_grid<- plot_grid(plotlist = nls.plot.list, ncol = 4) # plot it!
ggsave("figures/05_rTPC_modelcomparison_pop2.pdf", plot = mod_grid, width = 24, height = 16)

############### Loop through all populations, calculating AICc values ###################################

aicc_df <- data.frame( # list of models to consider
  model = c("Deutsch", "Lactin2", "Ratkowsky bounded", 
            "Rezende", "Ashrafi II", "Atkin", "Mitchell-Angilletta", 
            "Analytis-Kontodimas", "Eubank", "Taylor-Sexton", "Briere", "Thomas1"),
  stringsAsFactors = FALSE
)

for (i in 1:length(mat)){
  
  df.i <- subset(mat[[i]])
  df.i <- droplevels(df.i)
  
  aicc_pop <- c() # Temporary storage for this population
  
  # modified Deutsch
  
  start.vals <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'deutsch_2008')
  
  mod <- nls_multstart(r.exp ~ deutsch_2008(temp = Temp, rmax, topt, ctmax, a),
                            data = df.i,
                            iter = c(4, 4, 4, 4), 
                            start_lower = start.vals - 10,
                            start_upper = start.vals + 10,
                            lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'deutsch_2008'),
                            upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'deutsch_2008'),
                            supp_errors = 'Y',
                            convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Lactin 2
  
  start.vals <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'lactin2_1995')
  
  mod <- nls_multstart(r.exp ~ lactin2_1995(temp = Temp, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = start.vals - 10,
                           start_upper = start.vals + 10,
                           lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'lactin2_1995'),
                           upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'lactin2_1995'),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Bounded Ratkowsky
  
  start.vals <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'ratkowsky_1983')
  
  mod <- nls_multstart(r.exp ~ bounded_ratkowsky(temp = Temp, tmin, tmax, a, b),
                             data = df.i,
                             iter = c(4, 4, 4, 4), 
                             start_lower = start.vals - 10,
                             start_upper = start.vals + 10,
                             lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'ratkowsky_1983'),
                             upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'ratkowsky_1983'),
                             supp_errors = 'Y',
                             convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  #Rezende
  
  start.vals <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'rezende_2019')
  
  mod <- nls_multstart(r.exp ~ rezende_2019(temp = Temp, q10, a, b, c),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = start.vals - 10,
                           start_upper = start.vals + 10,
                           lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'rezende_2019'),
                           upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'rezende_2019'),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Ashrafi II
  
  mod <- nls_multstart(
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
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Atkin
  
  mod <- nls_multstart(
    r.exp ~ atkin(Temp, B0, a, b),
    data = df.i,
    iter = c(4, 4, 4),
    start_lower = c(B0 = 0.1, a = 10, b = 0.01),
    start_upper = c(B0 = 5, a = 30, b = 0.1),
    lower = c(B0 = 0, a = 0, b = 0),
    upper = c(B0 = 10, a = 50, b = 1),
    supp_errors = 'Y',
    convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Mitchell-Angilletta
  
  mod <- nls_multstart(
    r.exp ~ mitchell(Temp, a, b, Tpk),
    data = df.i,
    iter = c(4, 4, 4),
    start_lower = c(a = 0.1, b = 1, Tpk = 20),
    start_upper = c(a = 10, b = 10, Tpk = 40),
    lower = c(a = 0, b = 0.1, Tpk = 0),
    upper = c(a = Inf, b = Inf, Tpk = 50),
    supp_errors = 'Y',
    convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Analytis-Kontomodimas
  
  mod <- nls_multstart(
    r.exp ~ anal_kont(Temp, a, tmin, tmax),
    data = df.i,
    iter = c(4, 4, 4),
    start_lower = c(a = 0.1, tmin = 0, tmax = 30),
    start_upper = c(a = 10, tmin = 10, tmax = 50),
    lower = c(a = 0, tmin = -Inf, tmax = -Inf),
    upper = c(a = Inf, tmin = Inf, tmax = Inf),
    supp_errors = 'Y',
    convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Eubank
  
  mod <- nls_multstart(
    r.exp ~ eubank(Temp, a, tpk, b),
    data = df.i,
    iter = c(4, 4, 4),
    start_lower = c(a = 0.1, tpk = 20, b = 0.1),
    start_upper = c(a = 10, tpk = 35, b = 5),
    lower = c(a = 0, tpk = 0, b = 0.01),
    upper = c(a = Inf, tpk = Inf, b = Inf),
    supp_errors = 'Y',
    convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Taylor-Sexton
  
  mod <- nls_multstart(
    r.exp ~ tay_sex(Temp, bpk, tmin, tpk),
    data = df.i,
    iter = c(4, 4, 4),
    start_lower = c(bpk = 0.1, tmin = 0, tpk = 20),
    start_upper = c(bpk = 10, tmin = 10, tpk = 40),
    lower = c(bpk = 0, tmin = -Inf, tpk = -Inf),
    upper = c(bpk = Inf, tmin = Inf, tpk = Inf),
    supp_errors = 'Y',
    convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Briere II
  
  start.vals <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'briere2_1999')
  
  mod <- nls_multstart(r.exp~briere2_1999(temp = Temp, Tmin, Tmax, a, b),
                           data = df.i,
                           iter = c(4,4,4,4),
                           start_lower = start.vals - 10,
                           start_upper = start.vals + 10,
                           lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'briere2_1999'),
                           upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'briere2_1999'),
                           supp_errors = 'Y',
                           convergence_count = F)
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Thomas 1
  
  start.vals <- get_start_vals(df.i$Temp, df.i$r.exp, model_name = 'thomas_2012')
  
  mod <- nls_multstart(r.exp~thomas_2012(temp = Temp, a, b, c, topt),
                       data = df.i,
                       iter = c(4,4,4,4),
                       start_lower = start.vals - 10,
                       start_upper = start.vals + 10,
                       lower = get_lower_lims(df.i$Temp, df.i$r.exp, model_name = 'thomas_2012'),
                       upper = get_upper_lims(df.i$Temp, df.i$r.exp, model_name = 'thomas_2012'),
                       supp_errors = 'Y',
                       convergence_count = F)
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  aicc_df[[paste0("AICc_pop", i)]] <- aicc_pop #Add this population's AICc values as a new column
  
}

mean_aicc <- rowMeans(aicc_df[, -1], na.rm = TRUE)

aicc_df <- cbind(aicc_df[, 1, drop = FALSE], 
                 Mean_AICc = mean_aicc, 
                 aicc_df[, -1])

write.csv(aicc_df, "data-processed/07_nls_TPC_comparison_AICcs.csv") # Save model comparison data

############# Fit best models (Deutsch, Rezende, Lactin2) and extract shape parameters for all populations #####################

############# Bayesian modelling #########################

# OK, the next step is to try to fit our best performing models in R2jags - Deutsch, Rezende, and Lactin2
# This model framework (see https://github.com/JoeyBernhardt/anopheles-rate-summation/blob/master/AnalysisDemo.R)
# is predicated on text files (e.g. briere.txt) that specify model structure and parameters. New ones must be
# made and saved in the directory before calling them when fitting the TPCs themselves. 

# Let's set up our MCMC model settings
ni <- 55000 # iterations / chain
nb <- 5000 # burn in periods for each chain
nt <- 50 # thinning interval : (55,000 - 5,000) / 50 = 1000 posterior estimates / chain
nc <- 5 # number of chains

# Let's set up more generous MCMC settings
ni.2 <- 110000 # iterations / chain
nb.2 <- 10000 # burn in periods for each chain
nt.2 <- 100 # thinning interval : (110,000 - 10,000) / 100 = 1000 posterior estimates / chain

Temp.xs <- seq(0, 45, 0.1) # Temperature gradient we're interested in
N.Temp.xs <-length(Temp.xs)

i <- sample(1:38, 1) # OK we'll start with a random population again

df.i <- subset(mat[[i]])
df.i <- droplevels(df.i) 

head(df.i)

trait <- df.i$r.exp # Jags needs traits in a certain format
N.obs <- length(trait)
temp <- df.i$Temp

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

# Now we know the Briere function should work, because it's included as a txt file in the tutorial

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
bri_jag$BUGSoutput$summary[c(1:5, 110:120),]
mcmcplot(bri_jag) # Evaluate model performance
bri_jag$BUGSoutput$DIC # DIC

# plot!
plot(trait ~ jitter(Temp, 0.5), xlim = c(0, 45), data = df.i, 
     ylab = "Growth rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(bri_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bri_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bri_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

df.jags <- data.frame(bri_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:5),]
df.jags.plot$temp <- seq(0, 45, 0.1)

bri.jag.plot <- ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "orchid1", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "sienna", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temp, 0.5), y = trait), color = "gray12", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(1, 45)) + 
  scale_y_continuous(limits = c(-2, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Briere I model fit"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

bri.jag.plot

# OK, before moving on to my chosen models, I want to try a super simple model. Briere 2?

inits.bri2 <- function() {
  list(
    cf.a = runif(1, 0, 5),
    cf.tmin = runif(1, 0, 10),
    cf.tmax = runif(1, 30, 50),
    cf.b = runif(1, 1, 5),
    cf.sigma = runif(1, 0.5, 2)
  )
}

parameters.bri2 <- c("cf.a", "cf.tmin", "cf.tmax", "cf.b", "cf.sigma", "r.pred")

bri2_jag <- jags(data=jag.data, inits=inits.bri2, parameters.to.save=parameters.bri2, model.file="briere2stp.txt",
                n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())

bri2_jag$BUGSoutput$summary[1:6,] # Get estimates
mcmcplot(bri2_jag) # Evaluate model performance
bri2_jag$BUGSoutput$DIC # DIC

# plot!
plot(trait ~ jitter(Temp, 0.5), xlim = c(1, 45), data = df.i, 
     ylab = "Growth rate", xlab = expression(paste("Temperature (",degree,"C)")))
lines(bri2_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"] ~ Temp.xs, lty = 2)
lines(bri2_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"] ~ Temp.xs, lty = 2)
lines(bri2_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"] ~ Temp.xs)

preds <- data.frame( # Let's try this in ggplot - first generate some predictions
  Temp = Temp.xs,
  Mean = bri2_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"],
  LCL = bri2_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"],
  UCL = bri2_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"]
)

ggplot(data = preds, aes(x = Temp)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = Mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temp, 0.5), y = trait), color = "grey9", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(1, 45)) + 
  scale_y_continuous(limits = c(-0.5, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Briere 2 Model Fit"
  ) +
  theme_classic()

# OK, now let's try the Deutsch model!

# Deutsch

inits.deutsch <- function() {
  list(
    rmax = runif(1, 0, 5),  # More constrained initial values
    topt = runif(1, 30, 40),
    ctmax = runif(1, 37, 43),
    a = runif(1, 0.1, 5),
    sigma = runif(1, 0.5, 2)
  )
}

parameters.deutsch <- c("rmax", "topt", "ctmax", "a", "sigma", "r.pred")

deut_jag <- jags(
  data = jag.data, 
  inits = inits.deutsch, 
  parameters.to.save = parameters.deutsch, 
  model.file = "deutsch.txt",
  n.thin = nt.2, 
  n.chains = nc, 
  n.burnin = nb.2, 
  n.iter = ni.2, 
  DIC = TRUE, 
  working.directory = getwd()
)

deut_jag$BUGSoutput$summary[c(1:3, 455:457),] # Get estimates. For some reason rmax, sigma and topt are at the end!
mcmcplot(deut_jag) # Evaluate model performance
deut_jag$BUGSoutput$DIC # DIC

# This looks really nice!

df.jags <- data.frame(deut_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:3, 455:457),]
df.jags.plot$temp <- seq(0, 45, 0.1)

deut.jag.plot <- ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "azure3", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "firebrick", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temp, 0.5), y = trait), color = "darkgreen", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(0, 45)) + 
  scale_y_continuous(limits = c(-2, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Deutsch Model Fit"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)


# Lactin 2

inits.lactin2 <- function() {
  list(
    cf.a = runif(1, 0.01, 0.3),  # More constrained initial values
    cf.tmax = runif(1, 30, 40),
    cf.delta_t = runif(1, 1, 5),
    cf.b = runif(1, -5, 5),
    cf.sigma = runif(1, 0.5, 2)
  )
}

parameters.lactin2 <- c("cf.a", "cf.b", "cf.tmax", "cf.delta_t", "cf.sigma", "r.pred")

# This first one doesn't seem to work... not sure why? Has to do with the step function in the text file?

lac_jag <- jags(
  data = jag.data, 
  inits = inits.lactin2, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin.txt",
  n.thin = nt, 
  n.chains = nc, 
  n.burnin = nb, 
  n.iter = ni, 
  DIC = TRUE, 
  working.directory = getwd()
)

lac_jag$BUGSoutput$summary[1:6,] # Get estimates
mcmcplot(lac_jag) # Evaluate model performance
lac_jag$BUGSoutput$DIC # DIC

preds <- data.frame(
  Temp = Temp.xs,
  Mean = lac_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"],
  LCL = lac_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"],
  UCL = lac_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"]
)

ggplot(data = preds, aes(x = Temp)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = Mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temp, 0.5), y = trait), color = "grey9", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(1, 45)) + 
  scale_y_continuous(limits = c(-0.5, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Lactin 2 Model Fit"
  ) +
  theme_classic()

inits.lactin2 <- function() {
  list(
    cf.a = runif(1, 0.01, 0.2),  # More constrained initial values
    cf.tmax = runif(1, 35, 45),
    cf.delta_t = runif(1, 1, 5),
    cf.b = runif(1, -3, 3),
    cf.sigma = runif(1, 0.5, 2)
  )
}

# This next model doesn't really work that well.. The different chains are not converging.
# looks like our initial values and priors need to be tighter?


lac2_jag <- jags(
  data = jag.data, 
  inits = inits.lactin2, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin2.txt",
  n.thin = nt, 
  n.chains = nc, 
  n.burnin = nb, 
  n.iter = ni, 
  DIC = TRUE, 
  working.directory = getwd()
)

lac2_jag$BUGSoutput$summary[1:6,] # Get estimates
mcmcplot(lac2_jag) # Evaluate model performance
lac2_jag$BUGSoutput$DIC # DIC

preds <- data.frame(
  Temp = Temp.xs,
  Mean = lac2_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "mean"],
  LCL = lac2_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "2.5%"],
  UCL = lac2_jag$BUGSoutput$summary[6:(6 + N.Temp.xs - 1), "97.5%"]
)

ggplot(data = preds, aes(x = Temp)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = Mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temp, 0.5), y = trait), color = "grey9", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(1, 45)) + 
  scale_y_continuous(limits = c(-0.5, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Lactin 2 Model Fit"
  ) +
  theme_classic()

# Examine the results:

df.jags <- data.frame(lac2_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:6),]
df.jags.plot$temp <- seq(0, 45, 0.1)

ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temp, 0.5), y = trait), color = "grey9", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(0, 45)) + 
  scale_y_continuous(limits = c(-2, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Lactin 2 Model Fit"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

# OK, let's try tightening our priors and initial values...

inits.lactin2.1 <- function() {
  list(
    cf.a = runif(1, 0.01, 0.2),  # More constrained initial values
    cf.tmax = runif(1, 32, 38),
    cf.delta_t = runif(1, 0.1, 5),
    cf.b = runif(1, -1, 1),
    cf.sigma = runif(1, 0, 2)
  )
}

lac3_jag <- jags(
  data = jag.data, 
  inits = inits.lactin2.1, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin2_narrow.txt",
  n.thin = nt, 
  n.chains = nc, 
  n.burnin = nb, 
  n.iter = ni, 
  DIC = TRUE, 
  working.directory = getwd()
)

lac3_jag$BUGSoutput$summary[1:6,] # Get estimates
mcmcplot(lac3_jag) # Evaluate model performance
lac3_jag$BUGSoutput$DIC # DIC

# Examine the results:

df.jags <- data.frame(lac3_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:6),]
df.jags.plot$temp <- seq(0, 45, 0.1)

ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temp, 0.5), y = trait), color = "grey9", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(0, 45)) + 
  scale_y_continuous(limits = c(-2, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Lactin 2 Model Fit"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

inits.lactin2.2 <- function() {
  list(
    cf.a = runif(1, 0.05, 0.15),  # More constrained initial values
    cf.tmax = runif(1, 32, 38),
    cf.delta_t = runif(1, 1, 4),
    cf.b = runif(1, -0.5, 0.5),
    cf.sigma = runif(1, 0.1, 2)
  )
}

lac4_jag <- jags(
  data = jag.data, 
  inits = inits.lactin2.2, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin2_narrower.txt",
  n.thin = nt, 
  n.chains = nc, 
  n.burnin = nb, 
  n.iter = ni, 
  DIC = TRUE, 
  working.directory = getwd()
)

lac4_jag$BUGSoutput$summary[1:6,] # Get estimates
mcmcplot(lac4_jag) # Evaluate model performance
lac4_jag$BUGSoutput$DIC # DIC

df.jags <- data.frame(lac4_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:6),]
df.jags.plot$temp <- seq(0, 45, 0.1)

ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temp, 0.5), y = trait), color = "grey9", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(0, 45)) + 
  scale_y_continuous(limits = c(-2, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Lactin 2 Model Fit, Constrained"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

# Having a problem with independence, Tmax...

inits.lactin2.3 <- function() {
  list(
    cf.a = runif(1, 0.05, 0.15),  # More constrained initial values
    cf.tmax = runif(1, 37, 43),
    cf.delta_t = runif(1, 1, 4),
    cf.b = runif(1, -0.5, 0.5),
    cf.sigma = runif(1, 0.1, 2)
  )
}

lac5_jag <- jags(
  data = jag.data, 
  inits = inits.lactin2.3, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin2_narrower_lrgTmax.txt",
  n.thin = nt.2, 
  n.chains = nc, 
  n.burnin = nb.2, 
  n.iter = ni.2, 
  DIC = TRUE, 
  working.directory = getwd()
)

lac5_jag$BUGSoutput$summary[1:6,] # Get estimates
mcmcplot(lac5_jag) # Evaluate model performance
lac5_jag$BUGSoutput$DIC # DIC

df.jags <- data.frame(lac5_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:6),]
df.jags.plot$temp <- seq(0, 45, 0.1)

lact.jag.plot <- ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temp, 0.5), y = trait), color = "grey9", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(0, 45)) + 
  scale_y_continuous(limits = c(-2, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Lactin 2 Model Fit, Constrained, Doubled"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

# Rezende model?

inits.rezende <- function() {
  list(
    q10 = runif(1, 1.5, 5),    # Reasonable range for q10
    a = runif(1, 0.5, 5),      # Initial guess for a
    b = runif(1, 30, 40),      # Initial guess for threshold temperature
    c = runif(1, 0.1, 0.5),    # Initial guess for the decline parameter
    sigma = runif(1, 0.5, 2)   # Initial guess for error
  )
}

parameters.rezende <- c("q10", "a", "b", "c", "sigma", "r.pred")

rez_jag <- jags(
  data = jag.data, 
  inits = inits.rezende, 
  parameters.to.save = parameters.rezende, 
  model.file = "rezende.txt",
  n.thin = nt.2, 
  n.chains = nc, 
  n.burnin = nb.2, 
  n.iter = ni.2, 
  DIC = TRUE, 
  working.directory = getwd()
)

rez_jag$BUGSoutput$summary[c(1:5, 457),] # Get estimates. For some reason rmax, sigma and topt are at the end!
mcmcplot(rez_jag) # Evaluate model performance
rez_jag$BUGSoutput$DIC # DIC

df.jags <- data.frame(rez_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:3, 455:457),]
df.jags.plot$temp <- seq(0, 45, 0.1)

rez.jag.plot <- ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "azure3", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "blue4", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temp, 0.5), y = trait), color = "dodgerblue1", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(1, 45)) + 
  scale_y_continuous(limits = c(-2, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Rezende Model Fit"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

# Thomas model (Norberg)

inits.thomas <- function() {
  list(
    a = runif(1, 1, 5),      # Initial guess for a
    b = runif(1, -0.5, 0.5), # Initial guess for b
    c = runif(1, 5, 40),     # Initial guess for thermal niche width
    topt = runif(1, 30, 40), # Initial guess for topt
    sigma = runif(1, 0.5, 2) # Initial guess for error
  )
}

parameters.thomas <- c("a", "b", "c", "topt", "sigma", "r.pred")

thom_jag <- jags(
  data = jag.data, 
  inits = inits.thomas, 
  parameters.to.save = parameters.thomas, 
  model.file = "thomas.txt",
  n.thin = nt.2, 
  n.chains = nc, 
  n.burnin = nb.2, 
  n.iter = ni.2, 
  DIC = TRUE, 
  working.directory = getwd()
)

thom_jag$BUGSoutput$summary[c(1:4, 456:457),] # Get estimates. For some reason rmax, sigma and topt are at the end!
mcmcplot(thom_jag) # Evaluate model performance
thom_jag$BUGSoutput$DIC # DIC

df.jags <- data.frame(thom_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:4, 456:457),]
df.jags.plot$temp <- seq(0, 45, 0.1)

thom.jag.plot <- ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "orchid1", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "sienna", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temp, 0.5), y = trait), color = "gray12", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(1, 45)) + 
  scale_y_continuous(limits = c(-2, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Thomas 1 / Norberg Model Fit"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

# Let's make a summary plot with all 4 figures

plot_grid <- grid.arrange(deut.jag.plot, lact.jag.plot, rez.jag.plot, thom.jag.plot, ncol = 2)

ggsave("figures/06_R2jags_modelfits_pop26.pdf", plot_grid, width = 30, height = 30, units = "cm")

################## Model fitting ########################

# OK we are now going to loop through all of our populations, fitting and saving R2jags objects.
# We'll go with the Thomas and Lactin models for now.
# These are the two that allow negative grwowth rates at both ends, and Lactin 2
# is basically tied with Deutsch in our DIC comparisons for best model.

# Let's set more generous MCMC settings still for our final models. 
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model. 

summary.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Pop.fac = character(),    # Population name
  Pop.num = character(),    # Number assigned to population (not the same). This corresponds to the jags objects
  Model = character(),      # Model name
  T.min = numeric(),        # Minimum T
  T.max = numeric(),        # Maximum T
  T.br = numeric(),         # T breadth
  T.opt = numeric(),        # Optimal T
  r.max = numeric(),        # Maximum growth rate
  stringsAsFactors = FALSE  # Avoid factor conversion
)

dic.df <- data.frame(       # Save DIC scores for comparison
  Pop.fac = character(),    # Population name
  Pop.num = character(),    # Number assigned to population (not the same)       
  Model = character(),      # Model (lactin or thomas)
  DIC = numeric(),          # DIC           
  stringsAsFactors = FALSE            
)

inits.lactin.final <- function() { # The final initial values set we landed on after experimenting. 
  list(
    cf.a = runif(1, 0.05, 0.15),  # More constrained initial values
    cf.tmax = runif(1, 37, 43),
    cf.delta_t = runif(1, 1, 5),
    cf.b = runif(1, -2.5, -1),
    cf.sigma = runif(1, 0.1, 2)
  )
}

inits.thomas.final <- function() {
  list(
    a = runif(1, 1, 5),      # Initial guess for a
    b = runif(1, -0.5, 0.5), # Initial guess for b
    c = runif(1, 5, 40),     # Initial guess for thermal niche width
    topt = runif(1, 30, 40), # Initial guess for topt
    sigma = runif(1, 0.5, 2) # Initial guess for error
  )
}

parameters.lactin2 <- c("cf.a", "cf.b", "cf.tmax", "cf.delta_t", "cf.sigma", "r.pred") # repeated here
parameters.thomas <- c("a", "b", "c", "topt", "sigma", "r.pred") # repeated here

Temp.xs <- seq(0, 45, 0.05) # Temperature gradient we're interested in - upped the granularity here
N.Temp.xs <-length(Temp.xs)

for (i in 1:length(mat)){ # for each population
  
  df.i <- subset(mat[[i]])
  df.i <- droplevels(df.i)
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  temp <- df.i$Temp
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  lac_jag <- jags(
    data = jag.data, 
    inits = inits.lactin.final, 
    parameters.to.save = parameters.lactin2, 
    model.file = "lactin2_final.txt",
    n.thin = nt.fit, 
    n.chains = nc.fit, 
    n.burnin = nb.fit, 
    n.iter = ni.fit, 
    DIC = TRUE, 
    working.directory = getwd()
  )
  
  save(lac_jag, file = paste0("R2jags-objects/pop_", i, "_lactin2.RData")) # save the lactin2 model
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)[-c(1:6),]   # generate the sequence of r.pred values
  df.jags$temp <- seq(0, 45, 0.05)
  
  summary.df <- rbind(summary.df, data.frame(                             # Add summary data
    Pop.fac = df.i$pop.fac[1],                                            # Population name
    Pop.num = df.i$pop.num[1],                                            # Number assigned to population (not the same)
    Model = "Lactin 2",                                                   # Model name
    T.min = df.jags$temp[min(which(df.jags$mean > 0))],                   # Minimum T
    T.max = df.jags$temp[max(which(df.jags$mean > 0))],                   # Maximum T
    T.br = df.jags$temp[max(which(df.jags$mean > (max(df.jags$mean) / 2)))] - 
      df.jags$temp[min(which(df.jags$mean > (max(df.jags$mean) / 2)))],   # T breadth
    T.opt = df.jags$temp[which.max(df.jags$mean)],                        # Optimal T
    r.max = max(df.jags$mean)                                             # Maximum growth rate   
  ))
  
  dic.df <- rbind(dic.df, data.frame(         # Add DIC data
    Pop.fac = df.i$pop.fac[1],      # Population name
    Pop.num = df.i$pop.num[1],      # Number assigned to population (not the same)
    Model = "Lactin 2",             # Model name
    DIC = lac_jag$BUGSoutput$DIC    # DIC
  ))

}

for (i in 1:length(mat)){ # for each population
  
  df.i <- subset(mat[[i]])
  df.i <- droplevels(df.i)
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  temp <- df.i$Temp
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  lac_jag <- jags(
    data = jag.data, 
    inits = inits.lactin.final, 
    parameters.to.save = parameters.lactin2, 
    model.file = "lactin2_final.txt",
    n.thin = nt.fit, 
    n.chains = nc.fit, 
    n.burnin = nb.fit, 
    n.iter = ni.fit, 
    DIC = TRUE, 
    working.directory = getwd()
  )
  
  save(lac_jag, file = paste0("R2jags-objects/pop_", i, "_lactin2.RData")) # save the lactin2 model
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)[-c(1:6),]   # generate the sequence of r.pred values
  df.jags$temp <- seq(0, 45, 0.05)
  
  summary.df <- rbind(summary.df, data.frame(                             # Add summary data
    Pop.fac = df.i$pop.fac[1],                                            # Population name
    Pop.num = df.i$pop.num[1],                                            # Number assigned to population (not the same)
    Model = "Lactin 2",                                                   # Model name
    T.min = df.jags$temp[min(which(df.jags$mean > 0))],                   # Minimum T
    T.max = df.jags$temp[max(which(df.jags$mean > 0))],                   # Maximum T
    T.br = df.jags$temp[max(which(df.jags$mean > (max(df.jags$mean) / 2)))] - 
      df.jags$temp[min(which(df.jags$mean > (max(df.jags$mean) / 2)))],   # T breadth
    T.opt = df.jags$temp[which.max(df.jags$mean)],                        # Optimal T
    r.max = max(df.jags$mean)                                             # Maximum growth rate   
  ))
  
  dic.df <- rbind(dic.df, data.frame(         # Add DIC data
    Pop.fac = df.i$pop.fac[1],      # Population name
    Pop.num = df.i$pop.num[1],      # Number assigned to population (not the same)
    Model = "Lactin 2",             # Model name
    DIC = lac_jag$BUGSoutput$DIC    # DIC
  ))
  
}

write.csv(dic.df, "data-processed/08_DIC_values_BayesTPC.csv") # Save DIC table
write.csv(summary.df, "data-processed/09_TPC_shape_values_BayesTPC.csv") # Save DIC table
