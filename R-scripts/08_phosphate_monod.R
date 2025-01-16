# Jason R Laurich
# January 16, 2025

# This script will estimate r for each phosphate concentration as for light (script 06)
# I'm going to use a sliding-window approach to identify the exponential phase of the logged linear data, then fit an exponential growth
# curve to un-logged data for that time period for each phosphorous level.

# Then I'm going to fit Monod curves to those data (using R2jags) for each population to estimate R* (P*)

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

df <- read.csv("data-processed/12_phosphate_rstar_rfus_time.csv")
head(df) # RFU is density, days is time, phosphate_level is a factor (1 to 10). phosphate_concentration is what we want
str(df)

df<-df[,-c(1,2,4:6,8,10,11,13,14,16,17)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df$phos.conc <- as.numeric(df$phosphate_concentration)
df$phos_level <- factor(df$phosphate_level, levels = sort(unique(df$phosphate_level)), ordered = TRUE) # Keep the numerical sorting.
df$well_plate <- as.factor(df$well_plate)

df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac) # I don't recognize the COMBO group, I'm guessing this is a control?

df.exp <- subset(df, df$pop.fac != "COMBO") 
df.exp$well.ID<-as.factor(df.exp$well_plate)

N0.df <- df.exp %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df.exp <- df.exp %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.exp <- split(df.exp, df.exp$pop.num)  # Each element is a data frame for one population in df.exp

############# Loop through all populations ###################

df.r.exp <- data.frame( # Initializing a dataframe to store the results for each well, pop, and phosphorous level
  population = character(),
  population.number = numeric(),
  phos.lvl = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

phos <- as.vector(as.numeric(as.character(unique(df.exp$phos.conc)))) # for looping through nitrate levels
ord.phos<- sort(phos)

for (i in 1:length(mat.exp)){ # Looping through all of the populations
  
  for (t in ord.phos){ # and all of the phosphate levels
    
    df.it <- subset(mat.exp[[i]], phos.conc==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$days), ] # Just in case the phosphate data is also disordered. 
      
      if (df.it.wl$RFU[2] <  df.it.wl$RFU[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.wl <- df.it.wl[-1,]
        df.it.wl$N0 <- df.it.wl$RFU[1]
      }
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x T x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
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
          phos.lvl = df.it.wl$phos.conc[1],          # Phosphorous level
          well.ID = df.it.wl$well.ID[1],             # Well ID (assuming one well ID per subset)
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
        
      }
      
    }
    
  }
  
}

write.csv(df.r.exp, "data-processed/12a_phosphorous_r_estimates.csv") # let's save the file.

############# Exploration #######################

pops<-names(mat.exp) # We're going to plot out the growth curves for 5 populations - all 10 phosphorous levels. 
ran <- sample(pops, 5, replace = F) # OK, let's select 5 random populations.

plot.list <- list() # make a list of the plots

for (i in ran){ # Looping through all of the populations
  
  for (t in ord.phos){ # and all of the P levels
    
    df.it <- subset(mat.exp[[i]], phos.conc==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    pred.df.mod <- data.frame( # Dataframe to store fitted data for light levels growth curves. Will reset for each level, but that's OK, as we will save the plots first
      population.fac = character(),
      population.num = numeric(),
      phos.conc = numeric(),
      well.ID = character(),
      smt.days = numeric(),
      fit.RFU = numeric()
    )
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$days), ] # Not sure if this is needed for P, but can't hurt.
      
      if (df.it.wl$RFU[2] <  df.it.wl$RFU[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.wl <- df.it.wl[-1,]
        df.it.wl$N0 <- df.it.wl$RFU[1]
      }
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x T x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
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
      
      smt.days <- seq(min(df.it.wl.th$days, na.rm = TRUE), max(df.it.wl.th$days, na.rm = TRUE), length.out = 500) # Get a smooth distribution of time points
      
      # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
      fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
      pred.df.mod <- rbind(pred.df.mod, data.frame(
        population.fac = df.it.wl$pop.fac[1],
        population.num = i,
        phos.conc = t,
        well.ID = unique(df.it.wl$well.ID),
        smt.days = smt.days,
        fit.RFU = fit.RFU.mod
      ))
      
      # Now we plot it out
      
      p <- ggplot() + 
        geom_point(data=df.it, aes(x=days, y=RFU, colour= well.ID), size = 2) + # Observed data, not specific to individual well.IDs
        geom_line(data = pred.df.mod, aes(x = smt.days, y = fit.RFU, colour = well.ID), size = 1) + # Fitted line
        labs(title = paste("Population", i, "Phosphate concentration:", t), x = "Days", y = "RFU") +
        theme_minimal() +
        theme(legend.position = "none")
      
      plot.list[[paste0("Pop", df.it.wl.th$pop.fac[1], " Phosphate ", t)]] <- p
      
    }
  }
}

plot_grid <- arrangeGrob(grobs=plot.list, ncol = 10, nrow = 5)

grid.draw(plot_grid)   

ggsave("figures/10_phosphorous_r_estimatation.pdf", plot_grid, width = 40, height = 25, units = "cm")

############# Fit Monod curves to data ###################