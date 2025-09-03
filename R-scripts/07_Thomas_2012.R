# Jason R Laurich
# March 5, 2025

# Revisited, checked, and cleaned September 3rd, 2025 (JRL)

# In this script, I will work with the Thomas 2012 Science paper (doi: DOI: 10.1126/science.122483) supplement containing specific grwoth rates for 194 algae.

# Packages & functions ----------------------------------------------------

library(tidyr)
library(cowplot)
library(ggplot2)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(Deriv)
library(rTPC)
library(nls.multstart)

# Upload & examine data ---------------------------------------------------

df <- read.csv("data-processed/10_Thomas_2012_raw_data.csv") # Raw data containing growth rates at different temperatures for all species
head(df)
str(df)

length(unique(df$id.number))

mat <- split(df, df$id.number)  # Matrix

min(df[df$Growth.rate>0,]$Temperature) # Tmin => -1.8
max(df[df$Growth.rate>0,]$Temperature) # Tmax <= 37

# Testing -----------------------------------------------------------------

i <- sample(1:194,1) 

df.i <- subset(mat[[i]])

inits.lactin.meta<- function() { # The final initial values set we landed on after experimenting. 
  list(
    cf.a = runif(1, 0.05, 1),  # More constrained initial values
    cf.tmax = runif(1, 1, 40),    # Much wider to accommodate interspecific variation
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

lac_jag1 <- jags(
  data = jag.data, 
  inits = inits.lactin.meta, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin_meta.txt",
  n.thin = nt.fit, 
  n.chains = nc.fit, 
  n.burnin = nb.fit, 
  n.iter = ni.fit, 
  DIC = TRUE, 
  working.directory = getwd()
) # ~ 1 min to run...

lac_jag1$BUGSoutput$summary[1:6,] # terrible with the expanded cf.a and delta_t
mcmcplot(lac_jag1) # Looks terrible

df.jags <- data.frame(lac_jag1$BUGSoutput$summary)
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

lact.jag.plot # Hmm that fit looks OK! But there's some weird stuff going on...

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

lac_plot # Hmm this looks pretty good. 

# OK let's try customizing our initial values to fit something like the get_start_vals function does.

inits.lactin.thomas.cust<- function() { # Pulling initial values centres from the start_vals function in rTPC
  list(
    cf.a = rnorm(1, mean = start.vals.lac[1], sd = 0.05),
    cf.tmax = rnorm(1, mean = start.vals.lac[3], sd = 1),
    cf.delta_t = rnorm(1, mean = start.vals.lac[4], sd = 1),
    cf.b = rnorm(1, mean = start.vals.lac[2], sd = 0.05),
    cf.sigma = runif(1, 0.1, 2)
  )
}

lac_jag2 <- jags(
  data = jag.data, 
  inits = inits.lactin.thomas.cust, 
  parameters.to.save = parameters.lactin2, 
  model.file = "lactin2_thomas.txt",
  n.thin = nt.fit, 
  n.chains = nc.fit, 
  n.burnin = nb.fit, 
  n.iter = ni.fit, 
  DIC = TRUE, 
  working.directory = getwd()
) # ~ 1 min to run...

lac_jag2$BUGSoutput$summary[1:6,] # better, but cf.tmax and cf.b are still not looking great
mcmcplot(lac_jag2) # This might be workable, if I up the replication and stuff. 

df.jags <- data.frame(lac_jag2$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:6),]
df.jags.plot$temp <- seq(-5, 45, 0.1)

lact.jag.plot2 <- ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temperature, 0.5), y = Growth.rate), color = "grey9", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(-5, 45)) + 
  scale_y_continuous(limits = c(-1, 1)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Lactin 2 Model Fit, Constrained, Doubled"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

lact.jag.plot2

# OK let's run through extracting some summary metrics for practice here. 

df.jags.plot$temp[min(which(df.jags.plot$mean > 0))]                                     # Tmin
df.jags.plot$temp[max(which(df.jags.plot$mean > 0))]                                     # Tmax
df.jags.plot$temp[max(which(df.jags.plot$mean > (max(df.jags.plot$mean) / 2)))] - 
  df.jags.plot$temp[min(which(df.jags.plot$mean > (max(df.jags.plot$mean) / 2)))]        # T breadth
df.jags.plot$temp[which.max(df.jags.plot$mean)]                                          # Optimal T
max(df.jags.plot$mean)                                                                   # µmax

cf.a <- lac_jag2$BUGSoutput$summary[1,1] # Extract parameters
cf.b <- lac_jag2$BUGSoutput$summary[2,1]
cf.tmax <- lac_jag2$BUGSoutput$summary[5,1]
cf.delta_t <- lac_jag2$BUGSoutput$summary[3,1]

# TPC fitting -------------------------------------------------------------

thomas.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Pop.id = numeric(),       # Population id
  DIC = numeric(),          # DIC
  T.min.raw = numeric(),    # Minimum T (Jags raw)
  T.max.raw = numeric(),    # Maximum T (Jags raw)
  T.opt.raw = numeric(),    # Optimal T (Jags raw)
  r.max.raw = numeric(),    # Maximum growth rate (Jags raw)
  T.br.raw = numeric(),     # T breadth (Jags raw)
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Pop.id = numeric(),       # Population id
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~3000)
  stringsAsFactors = FALSE            
)

# Let's do larger models for the final things (10 times larger)
ni.fit <- 330000    # iterations / chain
nb.fit <- 30000     # burn in periods for each chain
nt.fit <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 3         # number of chains, total of 3,000 estimates for each model. 

parameters.lactin2 <- c("cf.a", "cf.b", "cf.tmax", "cf.delta_t", "cf.sigma", "r.pred") # repeated here

inits.lactin.meta<- function() { # Pulling initial values centres from the start_vals function in rTPC
  list(
    cf.a = rnorm(1, mean = start.vals.lac[1], sd = 0.05),
    cf.tmax = rnorm(1, mean = start.vals.lac[3], sd = 1),
    cf.delta_t = rnorm(1, mean = start.vals.lac[4], sd = 1),
    cf.b = rnorm(1, mean = start.vals.lac[2], sd = 0.05),
    cf.sigma = runif(1, 0.1, 2)
  )
}

for (i in 1:194){ # For each species
  
  df.i <- df %>%
    filter(id.number == i, !is.na(Growth.rate)) %>% 
    arrange(Temperature)
  
  max.temp <- max(df.i$Temperature, na.rm = TRUE)
  temp.at.max.growth <- df.i$Temperature[which.max(df.i$Growth.rate)]
  
  if (temp.at.max.growth == max.temp) {
    new_row <- df.i %>%
      filter(Temperature == max.temp) %>%
      mutate(
        Temperature = max.temp + 5,
        Growth.rate = 0
      )
    
    df.i <- bind_rows(df.i, new_row)
  }
  
  
  trait <- df.i$Growth.rate     # format the data for jags
  N.obs <- length(trait)
  temp <- df.i$Temperature
  
  Temp.xs <- seq(min(df.i$Temperature) - 5, max(df.i$Temperature) + 5, 0.1) # Temperature gradient we're interested in - upped the granularity here
  N.Temp.xs <-length(Temp.xs) # We'll reset this internally since the gradient varies substantially
  
  temp <- df.i$Temperature
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  start.vals.lac <- get_start_vals(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995')
  
  lac_jag <- jags(
    data = jag.data, 
    inits = inits.lactin.meta, 
    parameters.to.save = parameters.lactin2, 
    model.file = "lactin_meta.txt",
    n.thin = nt.fit, 
    n.chains = nc.fit, 
    n.burnin = nb.fit, 
    n.iter = ni.fit, 
    DIC = TRUE, 
    working.directory = getwd()
  ) # ~ 10 min to run?
  
  print(paste("Done", i, "of 194"))
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)[-c(1:6),]   # generate the sequence of r.pred values
  df.jags$temp <- seq(min(df.i$Temperature) - 5, max(df.i$Temperature) + 5, 0.1)
  
  thomas.summ.df <- rbind(thomas.summ.df, data.frame(                        # Add summary data
    Sp.id = i,                                                                         # Species name 
    DIC = lac_jag$BUGSoutput$DIC,                                                      # DIC
    T.min.raw = df.jags$temp[min(which(df.jags$mean > 0))],                            # Minimum T
    T.max.raw = df.jags$temp[max(which(df.jags$mean > 0))],                            # Maximum T
    T.br.raw = df.jags$temp[max(which(df.jags$mean > (max(df.jags$mean) / 2)))] - 
      df.jags$temp[min(which(df.jags$mean > (max(df.jags$mean) / 2)))],                # T breadth
    T.opt.raw = df.jags$temp[which.max(df.jags$mean)],                                 # Optimal T
    r.max.raw = max(df.jags$mean)                                                      # Maximum growth rate   
  ))
  
  for (j in 1:6){
    fit.df <- rbind(fit.df, data.frame(                         # Model performance data
      Sp.id = i,                                                # Species id   
      Parameter = rownames(lac_jag$BUGSoutput$summary)[j],      # Model parameter (e.g. cf.a, cf.tmax, etc.)
      mean = lac_jag$BUGSoutput$summary[j,1],                   # Posterior mean
      Rhat = lac_jag$BUGSoutput$summary[j,8],                   # Rhat values
      n.eff = lac_jag$BUGSoutput$summary[j,9],                  # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
}

write.csv(thomas.summ.df, "data-processed/10a_Thomas_2012_TPCs.csv") # Save Thomas 2012 summary table
write.csv(fit.df, "data-processed/10b_Thomas_2012_TPCs_fits.csv") # Save model fit summary table