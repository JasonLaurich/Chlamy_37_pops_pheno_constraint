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
library(rTPC)
library(nls.multstart)

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
    cf.a = runif(1, 0.05, 1),  # More constrained initial values
    cf.tmax = runif(1, 1, 40),    # Much wider to accomodate interspecific variation
    cf.delta_t = runif(1, 1, 10),
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
) # ~ 1 min to run...

lac_jag$BUGSoutput$summary[1:6,] # terrible with the expanded cf.a and delta_t
mcmcplot(lac_jag) # I can't live with that for sure

df.jags <- data.frame(lac_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:6),]
df.jags.plot$temp <- seq(-5, 45, 0.1)

lact.jag.plot <- ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temperature, 0.5), y = Growth.rate), color = "grey9", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(-5, 45)) + 
  scale_y_continuous(limits = c(-5, 5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Lactin 2 Model Fit, Constrained, Doubled"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

lact.jag.plot

# Let's try a frequentist approach

start.vals.lac <- get_start_vals(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995')

lac_nls <- nls_multstart(Growth.rate ~ lactin2_1995(temp = Temperature, a, b, tmax, delta_t),
                         data = df.i,
                         iter = c(4, 4, 4, 4), 
                         start_lower = start.vals.lac - 10,
                         start_upper = start.vals.lac + 10,
                         lower = get_lower_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995'),
                         upper = get_upper_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995'),
                         supp_errors = 'Y',
                         convergence_count = FALSE
)

summary(lac_nls)

preds.lac <- data.frame(Temperature = seq(min(df.i$Temperature - 2), max(df.i$Temperature +2), length.out = 100))
preds.lac <- broom::augment(lac_nls, newdata = preds.lac)

lac_plot <- ggplot(preds.lac) + geom_point(aes(Temperature, Growth.rate), df.i) +
  geom_line(aes(Temperature, .fitted), col = 'darkslateblue') + theme_classic() +ggtitle('Lactin 2') +ylim(-3,5) # This is substantially better. Why?

############# Fit TPCs and save the summary statistics #######################