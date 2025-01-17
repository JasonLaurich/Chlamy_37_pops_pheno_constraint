# Jason R Laurich
# January 16, 2025

# This script will estimate r for each salt concentration as for light (script 06)
# I'm going to use a sliding-window approach to identify the exponential phase of the logged linear data, then fit an exponential growth
# curve to un-logged data for that time period for each phosphorous level.

# Then I'm going to calculate salt tolerance
# This will be done using the method from Bernhardt et al 2020 (PNAS) — fitting a simplified logistic function and estimating c

############# Packages ########################

library(nls.multstart)
library(tidyr)
library(cowplot)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(R2jags)
library(mcmcplots)

############# Upload and examine data #######################

df <- read.csv("data-processed/13_chlamee_salt_rfus_time.csv")
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

############# Loop through all populations ###################

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

write.csv(df.r.exp, "data-processed/13a_salt_r_estimates.csv") # let's save the file.

############# Exploration #######################

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
  }
}

plot_grid <- arrangeGrob(grobs=plot.list, ncol = 10, nrow = 5)

grid.draw(plot_grid)   

ggsave("figures/11_salt_r_estimatation.pdf", plot_grid, width = 40, height = 25, units = "cm")

############# Fit Salt tolerance logistic curves to data ###################

df.r <- read.csv("data-processed/13a_salt_r_estimates.csv")

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

parameters.salt <- c("a", "b", "c", "sigma", "mu_pred") # Save these

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

mcmcplot(monod_jag) # Evaluate model performance
monod_jag$BUGSoutput$summary[c(1:3,2005),] # Get estimates

df.jags <- data.frame(monod_jag$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:3,2005),]
df.jags.plot$phos <- seq(0, 50, 0.025)

phos.jag.plot <- ggplot(data = df.jags.plot, aes(x = phos)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) + # Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  geom_point(data = df.i, aes(x = jitter(phos.conc, 0.5), y = trait), color = "grey9", size = 2) + # Add observed data points with jitter for P
  scale_x_continuous(limits = c(0, 50)) + 
  scale_y_continuous(limits = c(-0.25, 2.25)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Phosphorous (concentration)")),
    y = "Growth rate",
    title = "Monod curve, phosphorous limitation"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

R <- df.jags.plot$phos[which(df.jags.plot$mean > 0.56)[1]]
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

S.pred <- seq(0, 50, 0.025) # Repeated here
N.S.pred <-length(S.pred)

# Repeated here
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model.

for (i in 2:length(mat)){ # for each population. Starting at 2 again because I messed up the data entry phase, and did it manually for the first pop. Will work now properly in the loop
  
  df.i <- subset(mat[[i]])
  df.i <- droplevels(df.i)
  
  trait <- df.i$r.exp     # format the data for jags
  N.obs <- length(trait)
  phos <- df.i$phos.conc
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = phos, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod_jag <- jags( # Run the phosphorous Monod function. 
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
  
  save(monod_jag, file = paste0("R2jags-objects/pop_", i, "_phosphorous_monod.RData")) # save the phosphorous limitation monod function
  
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)[-c(1:3,2005),]   # generate the sequence of r.pred values
  df.jags$phos <- seq(0, 50, 0.025)
  
  summary.df <- rbind(summary.df, data.frame(                                                   # Add summary data
    Pop.fac = df.i$pop.fac[1],                                                                  # Population name
    Pop.num = df.i$pop.num[1],                                                                  # Number assigned to population (not the same)
    K.s = monod_jag$BUGSoutput$summary[1,1],                                                    # Half-saturation constant
    r.max = monod_jag$BUGSoutput$summary[3,1],                                                  # Maximum population growth rate
    R.jag = df.jags$phos[which(df.jags$mean > 0.56)[1]],                                        # Minimum resource requirement for positive growth (from jags model)
    R.mth = 0.56*monod_jag$BUGSoutput$summary[1,1]/(monod_jag$BUGSoutput$summary[3,1] - 0.56)   # Minimum resource requirement for positive growth (from math)                                                   
  ))
  
}

write.csv(summary.df, "data-processed/12b_phosphorous_Monod_estimates.csv") # Save summary table
