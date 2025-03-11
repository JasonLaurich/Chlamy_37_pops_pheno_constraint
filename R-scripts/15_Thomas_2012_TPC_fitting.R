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

lac_jag1 <- jags(
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

# Let's try the Norberg function

start.vals.thom <- get_start_vals(df.i$Temperature, df.i$Growth.rate, model_name = 'thomas_2012')

jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

parameters.thomas <- c("a", "b", "c", "topt", "sigma", "r.pred") # repeated here

inits.thom.cust<- function() { # The final initial values set we landed on after experimenting. 
  list(
    a = rnorm(1, mean = start.vals.thom[1], sd = 1),
    b = rnorm(1, mean = start.vals.lac[2], sd = 1),
    c = rnorm(1, mean = start.vals.lac[3], sd = 1),
    topt = rnorm(1, mean = start.vals.lac[4], sd = 1),
    sigma = runif(1, 0.1, 2)
  )
}

thom_jag <- jags(
  data = jag.data, 
  inits = inits.thom.cust, 
  parameters.to.save = parameters.thomas, 
  model.file = "thomas_intersp.txt",
  n.thin = nt.fit, 
  n.chains = nc.fit, 
  n.burnin = nb.fit, 
  n.iter = ni.fit, 
  DIC = TRUE, 
  working.directory = getwd()
)

thom_jag$BUGSoutput$summary[c(1:4, 506:507),] # OK.. this looks pretty good!
mcmcplot(thom_jag) # Ehhh... there's some spread on the variable fits....

df.jags <- data.frame(thom_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:4, 506:507),]
df.jags.plot$temp <- seq(-5, 45, 0.1)

thom.jag.plot <- ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "orchid1", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "sienna", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(Temperature, 0.5), y = Growth.rate), color = "gray12", size = 2) + # Add observed data points with jitter for Temp
  scale_x_continuous(limits = c(-5, 45)) + 
  scale_y_continuous(limits = c(-5, 5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "Thomas 1 / Norberg Model Fit"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

thom.jag.plot # Yikes!!!!!!!!!!

thom_nls <- nls_multstart(Growth.rate~thomas_2012(temp = Temperature, a, b, c, topt),
                          data = df.i,
                          iter = c(4,4,4,4),
                          start_lower = start.vals.thom - 10,
                          start_upper = start.vals.thom + 10,
                          lower = get_lower_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'thomas_2012'),
                          upper = get_upper_lims(df.i$Temperature, df.i$Growth.rate, model_name = 'thomas_2012'),
                          supp_errors = 'Y',
                          convergence_count = FALSE)

preds.thom <- data.frame(Temperature = seq(min(df.i$Temperature - 2), max(df.i$Temperature +2), length.out = 100))
preds.thom <- broom::augment(thom_nls, newdata = preds.thom)

thom_plot <- ggplot(preds.thom) + geom_point(aes(Temperature, Growth.rate), df.i) +
  geom_line(aes(Temperature, .fitted), col = 'darkslateblue') + theme_classic() + ggtitle('Thomas 1')

thom_plot # OK so this also looks terrible — I think we'll stick with Lactin II.

# OK let's run through extracting some summary metrics for practice here. 

df.jags.plot$temp[min(which(df.jags.plot$mean > 0))]                                     # Tmin
df.jags.plot$temp[max(which(df.jags.plot$mean > 0))]                                     # Tmax
df.jags.plot$temp[max(which(df.jags.plot$mean > (max(df.jags.plot$mean) / 2)))] - 
  df.jags.plot$temp[min(which(df.jags.plot$mean > (max(df.jags.plot$mean) / 2)))]        # T breadth
df.jags.plot$temp[which.max(df.jags.plot$mean)]                                          # Optimal T
max(df.jags.plot$mean)                                                                   # µmax

lactin2 <- function(temp, cf.a, cf.tmax, cf.delta_t, cf.b) { 
  exp(cf.a * temp) - exp(cf.a * cf.tmax - (cf.tmax - temp) / cf.delta_t) + cf.b
} # Define the Lactin II function

lactin2_deriv <- function(temp, cf.a, cf.b, cf.tmax, cf.delta_t) {
  rho <- cf.a
  T_max <- cf.tmax
  delta_T <- cf.delta_t
  
  term1 <- rho * exp(rho * temp)
  term2 <- (1 / delta_T) * exp(rho * T_max - (T_max - temp) / delta_T)
  
  return(term1 - term2)
} # Derivative of the Lactin II function

cf.a <- lac_jag2$BUGSoutput$summary[1,1] # Extract parameters
cf.b <- lac_jag2$BUGSoutput$summary[2,1]
cf.tmax <- lac_jag2$BUGSoutput$summary[5,1]
cf.delta_t <- lac_jag2$BUGSoutput$summary[3,1]

# Find the T_opt: where the derivative crosses zero
T_opt <- uniroot(
  function(temp) lactin2_deriv(temp, cf.a, cf.b, cf.tmax, cf.delta_t),
  interval = c(-10, 45)
)$root

r_max <- lactin2(temp=T_opt, cf.a=cf.a, cf.b=cf.b, cf.tmax=cf.tmax, cf.delta_t=cf.delta_t)

Tmin <- uniroot(lactin2, interval = c(0, T_opt), cf.a = cf.a, 
                cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t)$root

Tmax <- uniroot(lactin2, interval = c(T_opt,45), cf.a = cf.a, 
                cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t)$root

# OK we're going to modify the function to calculate T_breadth, based on a modified lactin.
lactin2_halfmax <- function(temp, cf.a, cf.b, cf.tmax, cf.delta_t, r_half) {
  exp(cf.a * temp) - exp(cf.a * cf.tmax - (cf.tmax - temp) / cf.delta_t) + cf.b - r_half
}

r_half <- r_max/2 # calculate half of rmax and get the roots.

Tlow <- uniroot(lactin2_halfmax, interval = c(Tmin, T_opt), cf.a = cf.a, 
                cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t, r_half = r_half)$root

Thigh <- uniroot(lactin2_halfmax, interval = c(T_opt, Tmax), cf.a = cf.a, 
                 cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t, r_half = r_half)$root

# The calculus isn't working for some reason. We;ll just save the a, b, etc. values. 

############# Fit TPCs and save the summary statistics #######################

r <- sample(1:194, 12, rep = F) # OK so for these ID's, we're going to save them and plot them

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

plot.list <- list() # Initialize a plot list to save things to (for pops j)

# Let's do larger models for the final things (10 times larger)
ni.fit <- 330000    # iterations / chain
nb.fit <- 30000     # burn in periods for each chain
nt.fit <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 3         # number of chains, total of 3,000 estimates for each model. 

# I need to run this in chunks so that I can save the output periodically. Read the follow files to get the dfs we can add to.

thomas.summ.df <- read.csv("data-processed/17_Thomas2012_TPCs.csv") # TPC parameters
fit.df <- read.csv("data-processed/17a_Thomas2012_TPCs_fits.csv")   # Model fit details

for (i in 188:194){
  
  df.i <- subset(mat[[i]])
  
  trait <- df.i$Growth.rate     # format the data for jags
  N.obs <- length(trait)
  temp <- df.i$Temperature
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)
  
  start.vals.lac <- get_start_vals(df.i$Temperature, df.i$Growth.rate, model_name = 'lactin2_1995')
  
  lac_jag <- jags(
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
  ) # ~ 10 min to run?
  
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)[-c(1:6),]   # generate the sequence of r.pred values
  df.jags$temp <- seq(-5, 45, 0.1)
  
  if(i %in% r){
    
    save(lac_jag, file = paste0("R2jags-objects/thomas2012_pop_", i, "_lactin2.RData")) # save the lactin2 model
    
    p <- ggplot(data = df.jags, aes(x = temp)) +
      geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
      geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
      geom_point(data = df.i, aes(x = jitter(Temperature, 0.5), y = Growth.rate), color = "grey9", size = 2) + # Add observed data points with jitter for Temp
      scale_x_continuous(limits = c(-5, 45)) + 
      scale_y_continuous(limits = c(-2, 5)) + # Customize the axes and labels +
      labs(
        x = expression(paste("Temperature (", degree, "C)")),
        y = "Growth rate") +
      theme_classic() +
      geom_hline(yintercept = 0)
    
    plot.list[[paste0("Population_", j)]] <- p
    
  }
  
  thomas.summ.df <- rbind(thomas.summ.df, data.frame(                                  # Add summary data
    Pop.id = i,                                                                        # Population name 
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
      Pop.id = i,                                               # Population id       
      Parameter = rownames(lac_jag$BUGSoutput$summary)[j],      # Model parameter (e.g. cf.a, cf.tmax, etc.)
      mean = lac_jag$BUGSoutput$summary[j,1],                   # Posterior mean
      Rhat = lac_jag$BUGSoutput$summary[j,8],                   # Rhat values
      n.eff = lac_jag$BUGSoutput$summary[j,9],                  # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
}

write.csv(thomas.summ.df, "data-processed/17_Thomas2012_TPCs.csv") # Save Thomas 2012 summary table
write.csv(fit.df, "data-processed/17a_Thomas2012_TPCs_fits.csv") # Save model fit summary table
