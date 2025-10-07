# Jason R Laurich
# January 15, 2025

# Revisited, checked, and cleaned September 3rd, 2025 (JRL)

# This script will estimate r for each nitrate concentration as for light (script 03)
# I'm going to use a sliding-window approach to identify the exponential phase of the logged linear data, then fit an exponential growth
# curve to un-logged data for that time period for each nitrogen level.

# Then I'm going to fit Monod curves to those data (using R2jags) for each population to estimate R* (N*)

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

df <- read.csv("data-processed/07_nitrogen_rfus_time.csv")
head(df) # RFU is density, days is time, nitrate_level is a factor (1 to 10). nitrate_concentration is what we want
str(df)

df<-df[,-c(1,2,4,5,8:10,12:14)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df$nitrate.conc <- as.numeric(df$nitrate_concentration)
df$nitrate_level <- factor(df$nitrate_level, levels = sort(unique(df$nitrate_level)), ordered = TRUE) # Keep the numerical sorting.
df$well_plate <- as.factor(df$well_plate)

# FIXME probably good to add reason (can't be zero) in this script too
df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac) # Removing the COMBO treatment, which is simply a control. 

df.exp <- subset(df, df$pop.fac != "COMBO")
df.exp$well.ID<-as.factor(df.exp$well_plate)

N0.df <- df.exp %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df.exp <- df.exp %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.exp <- split(df.exp, df.exp$pop.num)  # Each element is a data frame for one population in df.exp

# Estimate µ --------------------------------------------------------------

df.r.exp <- data.frame( # Initializing a dataframe to store the results for each well, pop, and nitrogen level
  population = character(),
  population.number = numeric(),
  nitrate.lvl = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

nitrate <- as.vector(as.numeric(as.character(unique(df.exp$nitrate.conc)))) # for looping through nitrate levels
ord.nit<- sort(nitrate)

# FIXME I suggest adding the reason for the loop here (fitting exponential) but to keep it short could refer to script where you did this with salt or temperature.
for (i in 1:length(mat.exp)){ # Looping through all of the populations
  
  for (t in ord.nit){ # and all of the nitrate levels
    
    df.it <- subset(mat.exp[[i]], nitrate.conc==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$days), ] # For whatever reason, the nitrogen data is not ordered properly by date in some cases.
      
      if (df.it.wl$RFU[2] <  df.it.wl$RFU[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.wl <- df.it.wl[-1,]
        df.it.wl$N0 <- df.it.wl$RFU[1]
      }
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x N level
      
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
      
      df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
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
          nitrate.lvl = df.it.wl$nitrate.conc[1],        
          well.ID = df.it.wl$well.ID[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl$pop.fac[1],          # Population as factor
          population.number = df.it.wl$pop.num[1],   # Numeric population number
          nitrate.lvl = df.it.wl$nitrate.conc[1],    # Nitrogen level
          well.ID = df.it.wl$well.ID[1],             # Well ID (assuming one well ID per subset)
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
        
      }
      
    }
    
  }
  
} # FIXME suggest adding percent done, also recieved the warning messages '1: In summary.lm(ln_slope) :essentially perfect fit: summary may be unreliable', would be good to flag if these are benign

write.csv(df.r.exp, "data-processed/07a_µ_estimates_nitrogen.csv") # let's save the file.

# Visualization -----------------------------------------------------------

pops<-names(mat.exp) # We're going to plot out the growth curves for 5 populations - all 10 nitrogen levels. 
ran <- sample(pops, 5, replace = F) # OK, let's select 5 random populations.

plot.list <- list() # make a list of the plots

for (i in ran){ # Looping through all of the populations
  
  for (t in ord.nit){ # and all of the N levels
    
    df.it <- subset(mat.exp[[i]], nitrate.conc==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    pred.df.mod <- data.frame( # Dataframe to store fitted data for N levels growth curves. Will reset for each level, but that's OK, as we will save the plots first
      population.fac = character(),
      population.num = numeric(),
      nitrate.conc = numeric(),
      well.ID = character(),
      smt.days = numeric(),
      fit.RFU = numeric()
    )
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$days), ] # For whatever reason, the nitrogen data is not ordered properly by date in some cases.
      
      if (df.it.wl$RFU[2] <  df.it.wl$RFU[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.wl <- df.it.wl[-1,]
        df.it.wl$N0 <- df.it.wl$RFU[1]
      }
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x N level
      
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
      
      df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      smt.days <- seq(min(df.it.wl.th$days, na.rm = TRUE), max(df.it.wl.th$days, na.rm = TRUE), length.out = 500) # Get a smooth distribution of time points
      
      # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
      fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
      pred.df.mod <- rbind(pred.df.mod, data.frame(
        population.fac = df.it.wl$pop.fac[1],
        population.num = i,
        nitrate.conc = t,
        well.ID = unique(df.it.wl$well.ID),
        smt.days = smt.days,
        fit.RFU = fit.RFU.mod
      ))
      
      # Now we plot it out
      
      p <- ggplot() + 
        geom_point(data=df.it, aes(x=days, y=RFU, colour= well.ID), size = 2) + # Observed data, not specific to individual well.IDs
        geom_line(data = pred.df.mod, aes(x = smt.days, y = fit.RFU, colour = well.ID), size = 1) + # Fitted line
        labs(title = paste("Population", i, "Nitrate concentration:", t), x = "Days", y = "RFU") +
        theme_minimal() +
        theme(legend.position = "none")
      
      plot.list[[paste0("Pop", df.it.wl.th$pop.fac[1], " Nitrate ", t)]] <- p
      
    }
  }
} # FIXME warnings on this loop: " In N0 * exp(r * days) :longer object length is not a multiple of shorter object length"

plot_grid <- arrangeGrob(grobs=plot.list, ncol = 10, nrow = 5)

grid.draw(plot_grid)   

ggsave("figures/05_µ_fits_nitrogen.pdf", plot_grid, width = 40, height = 25, units = "cm")

# Monod curves ------------------------------------------------------------

df.r <- read.csv("data-processed/07a_µ_estimates_nitrogen.csv")

head(df.r)
str(df.r)

df.r$pop.fac <- as.factor(df.r$population)
df.r$pop.num <- as.numeric(df.r$population.number)
df.r$well.ID <- as.factor(df.r$well.ID)
df.r$nitrate.conc <- as.numeric(df.r$nitrate.lvl)

mat <- split(df.r, df.r$pop.num)  # Matrixify the data!

i <- sample(1:37, 1) # We'll start by fitting a Monod curve to just one population

df.i <- subset(mat[[i]])
df.i <- droplevels(df.i) 

inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

trait <- df.i$r.exp     # format the data for jags
N.obs <- length(trait)
nit <- df.i$nitrate.conc

S.pred <- seq(0, 1000, 0.5) # Nitrogen gradient we are interested in here (concentration)
N.S.pred <-length(S.pred)

jag.data <- list(trait = trait, N.obs = N.obs, S = nit, S.pred = S.pred, N.S.pred = N.S.pred)

# We'll use the same larger number of chains and iterations we did for our final TPC models. 
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model. 

monod_jag <- jags( # Run the nitrogen Monod function. 
  data = jag.data,
  inits = inits.monod,
  parameters.to.save = parameters.monod,
  model.file = "monod.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)
# FIXME model started for me. Won't be able to run below due to mcmcplot:
mcmcplot(monod_jag) # Evaluate model performance
monod_jag$BUGSoutput$summary[c(1:3,2005),] # Get estimates

df.jags <- data.frame(monod_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:3,2005),]
df.jags.plot$nitrogen <- seq(0, 1000, 0.5)

nit.jag.plot <- ggplot(data = df.jags.plot, aes(x = nitrogen)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) + # Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(nitrate.conc, 0.5), y = trait), color = "grey9", size = 2) + # Add observed data points with jitter for N
  scale_x_continuous(limits = c(0, 1000)) + 
  scale_y_continuous(limits = c(-0.25, 2.25)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Nitrogen (concentration)")),
    y = "Growth rate",
    title = "Monod curve, nitrogen limitation"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

R <- df.jags.plot$nitrogen[which(df.jags.plot$mean > 0.56)[1]]
R

R2 <- 0.56*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.56)
R2

# OK let's loop through all of the populations: 

summary.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Pop.fac = character(),    # Population name
  Pop.num = character(),    # Number assigned to population (not the same). This corresponds to the jags objects
  K.s = numeric(),          # Half-saturation constant
  r.max = numeric(),        # Maximum population growth rate
  R.jag = numeric(),        # Minimum resource requirement for positive growth (from jags model)
  R.mth = numeric(),        # Minimum resource requirement for positive growth (analytical solution, R* = m*ks/(rmax-m))
  stringsAsFactors = FALSE  # Avoid factor conversion
)


inits.monod.final <- function() { # In case I want to play with these in the future
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Repeated here

S.pred <- seq(0, 1000, 0.5) # Repeated here
N.S.pred <-length(S.pred)

# Repeated here
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model.

for (i in 2:length(mat)){ # for each population. Fixed an error for population 1 at the data storage phase, so previously done.
  
  df.i <- subset(mat[[i]])
  df.i <- droplevels(df.i)
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  nit <- df.i$nitrate.conc
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = nit, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the nitrogen Monod function. 
    data = jag.data,
    inits = inits.monod,
    parameters.to.save = parameters.monod,
    model.file = "monod.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  save(monod_jag, file = paste0("R2jags-objects/pop_", i, "_nitrogen_monod.RData")) # save the nitrogen limitation monod function
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3,2005),]   # generate the sequence of r.pred values
  df.jags$nitrogen <- seq(0, 1000, 0.5)
  
  summary.df <- rbind(summary.df, data.frame(                                                   # Add summary data
    Pop.fac = df.i$pop.fac[1],                                                                  # Population name
    Pop.num = df.i$pop.num[1],                                                                  # Number assigned to population (not the same)
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                    # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                  # Maximum population growth rate
    R.jag = df.jags$nitrogen[which(df.jags$mean > 0.56)[1]],                                    # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.56*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.56)   # Minimum resource requirement for positive growth (from math)                                                   
  ))
  
} # FIXME same issue with package above

write.csv(summary.df, "data-processed/07b_Monod_nitrogen_estimates.csv") # Save summary table

# Analytical solutions and fit confirmation -------------------------------

# OK, so we are going to upload all of the R2jags objects, and iteratively use calculus to estimate µmax and R*

i<-1 # for now, starting with one population.
# for (i in 1:37){
load(paste0("R2jags-objects/pop_", i, "_nitrogen_monod.RData"))
#}

monod_jag$BUGSoutput$summary[c(1:3,2005),] # Get estimates# Look at the summary data

df.jags <- data.frame(monod_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:3,2005),]
df.jags.plot$nitrogen <- seq(0, 1000, 0.5)

nit.jag.plot <- ggplot(data = df.jags.plot, aes(x = nitrogen)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) + # Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(nitrate.conc, 0.5), y = r.exp), color = "grey9", size = 2) + # Add observed data points with jitter for Light
  scale_x_continuous(limits = c(0, 100)) + 
  scale_y_continuous(limits = c(-0.25, 2.25)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Nitrate (concentration)")),
    y = "Growth rate",
    title = "Monod curve, nitrogen limitation"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

nit.jag.plot # OK so µmax still seems lower (~1.3 compared to 3.8 ish for TPCs)

0.56*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.56) #this is our R* mathematical solution.
# This it the R* calculation from Bernhardt et al 2020. Probably don't need to do anything more.
# The model also directly estimates µmax (which looks accurate given the raw data) so no need for any further math.

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
  
  load(paste0("R2jags-objects/pop_", i, "_nitrogen_monod.RData"))
  
  nit_sum <- monod_jag$BUGSoutput$summary[c(1:3, 2005),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(        # Model performance data
      Grad = "Nitrogen limitation",            # Abiotic gradient
      Pop.fac = df.i$pop.fac[1],               # Population name
      Pop.num = df.i$pop.num[1],               # Number assigned to population (not the same)       
      Parameter = rownames(nit_sum)[j],        # Model parameter (e.g. K_s, r_max, etc.)
      mean = nit_sum[j,1],                     # Posterior mean
      Rhat = nit_sum[j,8],                     # Rhat values
      n.eff = nit_sum[j,9],                    # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
  
}

write.csv(fit.df, "data-processed/07c_Monod_nitrogen_fits.csv") # Save model fit summary table
# These look pretty good! A few models where n.eff is a touch low (4 estimates < 2,000, 10 < 3,000)
