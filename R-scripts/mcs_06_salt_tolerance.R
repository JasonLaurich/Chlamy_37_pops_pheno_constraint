# Jason R Laurich
# January 16, 2025

# Revisited, checked, and cleaned September 3rd, 2025 (JRL)

# This script will estimate µ for each salt concentration as for light (script 03)
# I'm going to use a sliding-window approach to identify the exponential phase of the logged linear data, then fit an exponential growth
# curve to un-logged data for that time period for each phosphorous level.

# Then I'm going to calculate salt tolerance
# This will be done using the method from Bernhardt et al 2020 (PNAS) — fitting a simplified reversed logistic growth function and estimating c

# Packages & functions ----------------------------------------------------

library(nls.multstart)
library(tidyr)
library(cowplot)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(R2jags)
library(mcmcplots)

# Upload & examine data ---------------------------------------------------

df <- read.csv("data-processed/09_salt_rfus_time.csv")
head(df) # RFU is density, time is time, treatment is the salt level s0 to S9.
str(df)

df<-df[,-c(2:4,8:9)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df$salt_level <- factor(df$treatment, levels = sort(unique(df$treatment)), ordered = TRUE) # Keep the numerical sorting.
df$salt <- as.numeric(df$salt_level) # These are just numbers for now reflecting the ordered factor.

df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac)

df.exp <- df
df.exp$well.ID<-as.factor(df.exp$unique_id)

N0.df <- df.exp %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(time)]) %>%
  ungroup()

df.exp <- df.exp %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.exp <- split(df.exp, df.exp$pop.num)  # Each element is a data frame for one population in df.exp

# Estimate µ --------------------------------------------------------------

df.r.exp <- data.frame( # Initializing a dataframe to store the results for each well, pop, and salt level
  population = character(),
  population.number = numeric(),
  salt.lvl = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

salt <- as.vector(as.numeric(as.character(unique(df.exp$salt)))) # for looping through nitrate levels
ord.salt<- sort(salt)

for (i in 1:length(mat.exp)){ # Looping through all of the populations
  
  for (t in ord.salt){ # and all of the salt levels
    
    df.it <- subset(mat.exp[[i]], salt==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$time), ] # Just in case the salt data is also disordered. 
      
      if (df.it.wl$RFU[2] <  df.it.wl$RFU[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.wl <- df.it.wl[-1,]
        df.it.wl$N0 <- df.it.wl$RFU[1]
      }
      
      t.series <- unique(df.it.wl$time) # Re-initialize this internally - we will only save summary data for each unique pop x salt x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$time <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~time, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x salt level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 5%
        
        if (ln.slopes[r] != 0 & round(ln.slopes[r], digits = 5) !=0) { # We also need to account for tiny values that are basically 0 (e.g. 5 e-16, but are messing up our loops)
          # percent.chg <- (ln.slopes[s] - ln.slopes[s + 1]) / ln.slopes[s] This was the reason we were getting weird negative results! If there was a stochastic drop.
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] # This should (along with the next line) fix it.
          
          if (percent.chg >= 1.10 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { # Because now, the drop-off ignores transiently negative slopes. 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
      
      df.it.wl.th <- df.it.wl[df.it.wl$time <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*time),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl.th$pop.fac[1],          
          population.number = df.it.wl$pop.num[1],      
          salt.lvl = df.it.wl$salt[1],        
          well.ID = df.it.wl$well.ID[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl$pop.fac[1],          # Population as factor
          population.number = df.it.wl$pop.num[1],   # Numeric population number
          salt.lvl = df.it.wl$salt[1],               # Salt level
          well.ID = df.it.wl$well.ID[1],             # Well ID (assuming one well ID per subset)
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
        
      }
      
    }
    
  }
  
}
# FIXME suggest adding percent done line. Got 19 warnings: "In summary.lm(ln_slope) : essentially perfect fit: summary may be unreliable", I suggest making a note in the script if these are ok.

write.csv(df.r.exp, "data-processed/09a_µ_estimates_salt.csv") # let's save the file.

# Visualization -----------------------------------------------------------

pops<-names(mat.exp) # We're going to plot out the growth curves for 5 populations - all 10 salt levels. 
ran <- sample(pops, 5, replace = F) # OK, let's select 5 random populations.

plot.list <- list() # make a list of the plots

for (i in ran){ # Looping through all of the populations
  
  for (t in ord.salt){ # and all of the salt levels
    
    df.it <- subset(mat.exp[[i]], salt==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    pred.df.mod <- data.frame( # Dataframe to store fitted data for salt levels growth curves. Will reset for each level, but that's OK, as we will save the plots first
      population.fac = character(),
      population.num = numeric(),
      phos.conc = numeric(),
      well.ID = character(),
      smt.days = numeric(),
      fit.RFU = numeric()
    )
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$time), ] # Not sure if this is needed for salt, but can't hurt.
      
      if (df.it.wl$RFU[2] <  df.it.wl$RFU[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.wl <- df.it.wl[-1,]
        df.it.wl$N0 <- df.it.wl$RFU[1]
      }
      
      t.series <- unique(df.it.wl$time) # Re-initialize this internally - we will only save summary data for each unique pop x salt x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$time <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~time, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x salt level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 10%
        
        if (ln.slopes[r] != 0 & round(ln.slopes[r], digits = 5) !=0) { # We also need to account for tiny values that are basically 0 (e.g. 5 e-16, but are messing up our loops)
          # percent.chg <- (ln.slopes[s] - ln.slopes[s + 1]) / ln.slopes[s] This was the reason we were getting weird negative results! If there was a stochastic drop.
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] # This should (along with the next line) fix it.
          
          if (percent.chg >= 1.10 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { # Because now, the drop-off ignores transiently negative slopes. 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
      
      df.it.wl.th <- df.it.wl[df.it.wl$time <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*time),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      smt.days <- seq(min(df.it.wl.th$time, na.rm = TRUE), max(df.it.wl.th$time, na.rm = TRUE), length.out = 500) # Get a smooth distribution of time points
      
      # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
      fit.RFU.mod <- predict(r_exp, newdata = data.frame(time = smt.days))
      pred.df.mod <- rbind(pred.df.mod, data.frame(
        population.fac = df.it.wl$pop.fac[1],
        population.num = i,
        salt.conc = t,
        well.ID = unique(df.it.wl$well.ID),
        smt.days = smt.days,
        fit.RFU = fit.RFU.mod
      ))
      
      # Now we plot it out
      
      p <- ggplot() + 
        geom_point(data=df.it, aes(x=time, y=RFU, colour= well.ID), size = 2) + # Observed data, not specific to individual well.IDs
        geom_line(data = pred.df.mod, aes(x = smt.days, y = fit.RFU, colour = well.ID), size = 1) + # Fitted line
        labs(title = paste("Population", i, "Salt concentration:", t), x = "Days", y = "RFU") +
        theme_minimal() +
        theme(legend.position = "none")
      
      plot.list[[paste0("Pop", df.it.wl.th$pop.fac[1], " Salt ", t)]] <- p
      
    }
  } # FIXME I set i to "9" and t to 1 for testing, and from the above loop got the warning "Warning message:In N0 * exp(r * time) : longer object length is not a multiple of shorter object length", which I want to make sure is not an issue. Otherwise it generate the plot just fine.
}

plot_grid <- arrangeGrob(grobs=plot.list, ncol = 10, nrow = 5)

grid.draw(plot_grid)   

ggsave("figures/07_µ_fits_salt.pdf", plot_grid, width = 40, height = 25, units = "cm")

# Salt tolerance curves ---------------------------------------------------

df.r <- read.csv("data-processed/09a_µ_estimates_salt.csv")

head(df.r)
str(df.r)

df.r$pop.fac <- as.factor(df.r$population)
df.r$pop.num <- as.numeric(df.r$population.number)
df.r$well.ID <- as.factor(df.r$well.ID)
df.r$salt.lvl <- as.numeric(df.r$salt.lvl)

mat <- split(df.r, df.r$pop.num)  # Matrixify the data!

i <- sample(1:37, 1) # We'll start by fitting a Monod curve to just one population

df.i <- subset(mat[[i]])
df.i <- droplevels(df.i) 

# Define the initial values for the parameters for the salt tolerance logistic curve
inits.salt <- function() {
  list(
    a = runif(1, 0.1, 5),  # Initial guess for a
    b = runif(1, 0.1, 5),  # Initial guess for b
    c = runif(1, 0.1, max(df.i$salt.lvl)),  # Initial guess for c
    sigma = runif(1, 0.1, 2)  # Initial guess for error
  )
}

parameters.salt <- c("a", "b", "c", "sigma", "r_pred_new") # Save these

trait <- df.i$r.exp     # format the data for jags
N.obs <- length(trait)
salt <- df.i$salt.lvl

S.pred <- seq(0, 10, 0.005) # Salt gradient we are interested in (here arbitrary 0 - 10, need to get actual []s) - we'll keep N.S.pred consistent across abiotic gradients for now
N.S.pred <-length(S.pred)

jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)

# We'll use the same larger number of chains and iterations we did for our final TPC models. 
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model. 

monod_jag <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)
# FIXME initialization worked, can't test below due to mcmcplot deprication.

mcmcplot(monod_jag) # Evaluate model performance
monod_jag$BUGSoutput$summary[c(1:4,2006),] # Get estimates

df.jags <- data.frame(monod_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:4,2006),]
df.jags.plot$salt <- seq(0, 10, 0.005)

salt.jag.plot <- ggplot(data = df.jags.plot, aes(x = salt)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) + # Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(salt.lvl, 0.5), y = trait), color = "grey9", size = 2) + # Add observed data points with jitter for salt
  scale_x_continuous(limits = c(0, 10)) + 
  scale_y_continuous(limits = c(-0.25, 2.25)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Salt (concentration)")),
    y = "Growth rate",
    title = "Logistic growth curve, salt stress"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

c <- df.jags.plot$salt[which.min(abs(df.jags.plot$mean - (monod_jag$BUGSoutput$summary[1,1] / 2)))]

c2 <- monod_jag$BUGSoutput$summary[3,1]

# OK let's loop through all of the populations: 

summary.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Pop.fac = character(),    # Population name
  Pop.num = character(),    # Number assigned to population (not the same). This corresponds to the jags objects
  r.max = numeric(),        # Maximum population growth rate (alpha)
  c.mod = numeric(),        # salt concentration at which r is half of alpha (extracted from model)
  c.pred = numeric(),       # salt concentration at which r is half of alpha (extracted from predicted values)
  stringsAsFactors = FALSE  # Avoid factor conversion
)

inits.salt <- function() { # In case I want to tweak these
  list(
    a = runif(1, 0.1, 5),  # Initial guess for a
    b = runif(1, 0.1, 5),  # Initial guess for b
    c = runif(1, 0.1, max(df.i$salt.lvl)),  # Initial guess for c
    sigma = runif(1, 0.1, 2)  # Initial guess for error
  )
}

parameters.salt <- c("a", "b", "c", "sigma", "r_pred_new") # Repeats

S.pred <- seq(0, 10, 0.005) # Repeats
N.S.pred <-length(S.pred)

# Repeated here
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model.

for (i in 1:length(mat)){ # for each population.
  
  df.i <- subset(mat[[i]])
  df.i <- droplevels(df.i)
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  salt <- df.i$salt.lvl
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the salt logistic growth curve function. 
    data = jag.data,
    inits = inits.salt,
    parameters.to.save = parameters.salt,
    model.file = "salt.tol.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  save(monod_jag, file = paste0("R2jags-objects/pop_", i, "_salt_tolerance.RData")) # save the salt tolerance reverse logistic model object
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:4,2006),]   # generate the sequence of r.pred values
  df.jags$salt <- seq(0, 10, 0.005)
  
  summary.df <- rbind(summary.df, data.frame(                                                                 # Add summary data
    Pop.fac = df.i$pop.fac[1],                                                                                # Population name
    Pop.num = df.i$pop.num[1],                                                                                # Number assigned to population (not the same)
    r.max = monod_jag$BUGSoutput$summary[1,1],                                                                # Maximum population growth rate
    c.mod = monod_jag$BUGSoutput$summary[3,1],                                                                # salt concentration at which r is half of alpha (extracted from model)
    c.pred = df.jags.plot$salt[which.min(abs(df.jags.plot$mean - (monod_jag$BUGSoutput$summary[1,1] / 2)))]   # salt concentration at which r is half of alpha (extracted from predicted values)                                                   
  ))
  
}

write.csv(summary.df, "data-processed/09b_salt_tolerance_estimates.csv") # Save summary table

# Analytical solutions & fit confirmation ---------------------------------

# OK, so we are going to upload all of the R2jags objects, and iteratively use calculus to estimate µmax and R*

i<-1 # for now, starting with one population.
# for (i in 1:37){
load(paste0("R2jags-objects/pop_", i, "_salt_tolerance.RData"))
#}

monod_jag$BUGSoutput$summary[c(1:4,2006),] # Get estimates# Look at the summary data

df.jags <- data.frame(monod_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:4,2006),]
df.jags.plot$salt <- seq(0, 10, 0.005)

salt.jag.plot <- ggplot(data = df.jags.plot, aes(x = salt)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) + # Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(salt.lvl, 0.5), y = r.exp), color = "grey9", size = 2) + # Add observed data points with jitter for Light
  scale_x_continuous(limits = c(0, 10)) + 
  scale_y_continuous(limits = c(-0.25, 2.25)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Salt (concentration)")),
    y = "Growth rate",
    title = "Reverse logistic curve, salt stress"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 1.4885590/2) +
  geom_vline(xintercept = 3.2465563)

salt.jag.plot # OK so µmax still seems lower (~1.3 compared to 3.8 ish for TPCs)

df.jags.plot$salt[which.min(abs(df.jags.plot$mean - (monod_jag$BUGSoutput$summary[1,1] / 2)))] #this is our solution from math.
# This it the calculation from Bernhardt et al 2020. Probably don't need to do anything more.
# The model also directly estimates µmax (which looks accurate given the raw data) so no need for any further math.

# Let's try solving this with math!
solve_salt <- function(a, b, c) {
  salt <- c + log((a / (a/2)) - 1) / b
  return(salt)
}

a <- monod_jag$BUGSoutput$summary[1,1]
b <- monod_jag$BUGSoutput$summary[2,1]
c <- monod_jag$BUGSoutput$summary[3,1]

c.mth <- solve_salt(a, b, c)
c.mth # Gives the same number. 


# Which leaves us with iteratively extracting model fit parameters. Let's make sure Rhat and n.eff values look good!

fit.df <- data.frame(       # Save model fit estimates for examination
  Grad = character(),       # Specifiy abiotic gradient
  Pop.fac = character(),    # Population name
  Pop.num = character(),    # Number assigned to population (not the same)       
  Parameter = character(),  # Model parameter (e.g. K_s, r_max, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~6000)
  stringsAsFactors = FALSE            
)

for (i in 1:37){
  
  df.i <- subset(mat[[i]])
  df.i <- droplevels(df.i)
  
  load(paste0("R2jags-objects/pop_", i, "_salt_tolerance.RData"))
  
  salt_sum <- monod_jag$BUGSoutput$summary[c(1:4, 2006),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(        # Model performance data
      Grad = "Salt stress",                    # Abiotic gradient
      Pop.fac = df.i$pop.fac[1],               # Population name
      Pop.num = df.i$pop.num[1],               # Number assigned to population (not the same)       
      Parameter = rownames(salt_sum)[j],       # Model parameter (e.g. K_s, r_max, etc.)
      mean = salt_sum[j,1],                    # Posterior mean
      Rhat = salt_sum[j,8],                    # Rhat values
      n.eff = salt_sum[j,9],                   # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
  
}

write.csv(fit.df, "data-processed/09c_salt_tolerance_fits.csv") # Save model fit summary table
# These are substantially worse at the moment (2 estimates < 1,000, 5 < 2,000, 10 < 3,000) but workable. 
