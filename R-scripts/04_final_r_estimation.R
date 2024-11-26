# Jason R Laurich
# November 20, 2024

# This script will estimate r one final way (which is the way we intend to calculate estimates used for the fitting of TPCs)
# For all sub-40 C populations, I will fit logged linear models to our data, within a window beginning at our first time point.

# As the width of that window increases, I will calculate the slope of the line, and interpret a decrease in the slope as the
# end of the exponential growth phase.
# I will then fit an exponential growth curve to that period of data for each population and well replicate of Chlamydomonas.

# For 40 C, our data shows a brief growth spike that represents anomalous data (a completion of cell division already initiated)
# that is not biologically meaningful. 
# Thus, for this treatment, I will threshold data to after maximum RFU count, and fit the same exponential growth curve to it.

############# Packages ########################

library(nls.multstart)
library(tidyr)
library(cowplot)
library(ggplot2)
library(dplyr)

############# Upload and examine data #######################

df <- read.csv("data-processed/chlamee-acute-rfu-time.csv")
head(df) #RFU is density, days is time, temperature is temperature
str(df)

df<-df[,-c(1,2,4,5,6,8,13)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df<-subset(df,df$temperature!=20)

# df$days <- df$days + 0.001 # Can't have 0s
df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac) # I don't recognize the cc1629 group, but we'll keep it for now

df.rep <- subset(df, df$plate_type == "repeat") # We're going to work only with repeat data (this is most of the well_plates) which I'll use as replicates. 
df.rep$well.ID<-as.factor(df.rep$well_plate)

N0.df <- df.rep %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df.rep <- df.rep %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.rep <- split(df.rep, df.rep$pop.num)  # Each element is a data frame for one population in df.rep

############## Data exploration #####################

df.it <- subset(mat.rep[[1]], temperature==10) # For one population and temperature, let's explore the model
df.it <- droplevels(df.it)

df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == 4)

head(df.it.wl)

t.series <- unique(df.it.wl$days) # First, we identify the string of unique days.
t.0 <- t.series[1] # Get the starting point
t.series <- t.series[-1] # Ok so this is what we will loop through! t.0 will always be the starting point, but the end point will change

ln.slopes <- c() # Store the logged linear slopes for each sliding window, from an lm
sl.direct <- c() # Store the directly calculated (rise/run)

for (t in t.series){ # This first loop will calculate the slopes for my first explorative population.
  
  df.it.wl.sl <- df.it.wl[df.it.wl$days <= t, ] # Subset the data to exclude time points above our window
  
  ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
  
  ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
  
  slope <- (df.it.wl.sl$logRFU[nrow(df.it.wl.sl)] - df.it.wl.sl$logRFU[1])/(df.it.wl.sl$days[nrow(df.it.wl.sl)] - df.it.wl.sl$days[1])
  sl.direct <- c(sl.direct, slope)
  
}

s <- length(ln.slopes) # Initialize the default to the final number, in case we don't see a drop off. 
  
for (r in 1:(length(ln.slopes) - 1)) { # Now we will loop through the ln.slopes vector to find when the slope decreases by more than 5%
    
  percent.chg <- ln.slopes[r] / ln.slopes[r + 1]
  if (percent.chg >= 1.05 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { # We don't want to include negative slopes that are driven by a brief decrease in RFUs at the start!
    s <- r # If the condition is met, reassign s to the corresponding drop-off point!
    break  # Exit loop when condition is met. s will store the final time point data!
  }
}

# OK so we now have (in s) an indicator of when we should be thresholding growth data to limit it to the exponential phase.
# This will work with ln.slopes AND t.series (because we excluded the t_0 data point from it)
  
df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
  
r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                          data = df.it.wl.th,
                          start_lower = c(r = 0.2), 
                          start_upper = c(r = 4.5),   
                          iter = 500,
                          supp_errors = 'Y',
                          control = nls.control(maxiter = 200))
  
summary(r_exp) # Looks good!
  
pred.time <- seq(min(df.it.wl.th$days), max(df.it.wl.th$days), length.out = 100) # Let's plot this
preds <- data.frame(days = pred.time)
preds$RFU.pred <- predict(r_exp, newdata = preds)
  
ggplot(df.it.wl.th, aes(x = days, y = RFU)) + # Plot observed and fitted data
  geom_point(color = "blue", size = 2) +
  geom_line(data = preds, aes(x = days, y = RFU.pred), color = "red", size = 1) +
  labs(title = "Exponential growth curve, Pop 16, Well B05_71, T 34C",
        x = "Days",
        y = "RFU")
  
# OK this is looking great! Let's try to loop this through all of my populations now.

# But first, for the 40 C data, we will have to use a different approach. 

df.it <- subset(mat.rep[[6]], temperature==40) # For one population and temperature, let's explore the model
df.it <- droplevels(df.it)

df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == 1)

# I think for now, we will just exclude data after that immediate growth spike that is not of biological interest

df.it.wl.late <- df.it.wl[df.it.wl$days > df.it.wl$days[which.max(df.it.wl$RFU)], ] # Do this based on the time point associated with max RFUs

r_ln <- lm(logRFU~days, data= df.it.wl.late)

summary(r_ln)

# OK now onto looping! We'll start with the non 40 C data. 

df.r.exp <- data.frame( # Initializing a dataframe to store the results for each well, pop, and temperature
  population = character(),
  population.number = numeric(),
  temperature = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

# Need to exclude 40 C data!

df.rep.34 <- subset(df.rep, df.rep$temperature!=40)
mat.rep.34 <- split(df.rep.34, df.rep.34$pop.num)  # Each element is a data frame for one population in df.rep

tmp <- as.vector(as.numeric(as.character(unique(df.rep.34$temperature)))) # for looping through temperatures
ord.tmp<- sort(tmp)

for (i in 1:length(mat.rep.34)){ # Looping through all of the populations
  
  for (t in ord.tmp){ # and all of the temperatures
    
    df.it <- subset(mat.rep.34[[i]], temperature==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x T x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (t in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= t, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x T level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 5%
        
        if (ln.slopes[r] != 0) { # Check for valid denominator to avoid division by zero
          # percent.chg <- (ln.slopes[s] - ln.slopes[s + 1]) / ln.slopes[s] This was the reason we were getting weird negative results! If there was a stochastic drop.
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] # This should (along with the next line) fix it.
          
          if (percent.chg >= 1.05 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { # Because now, the drop-off ignores transiently negative slopes. 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
        
        df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
        
        r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                               data = df.it.wl.th,
                               start_lower = c(r = 0.2), 
                               start_upper = c(r = 4.5),   
                               iter = 500,
                               supp_errors = 'Y',
                               control = nls.control(maxiter = 200))
        
        if (is.null(r_exp)){
          
          df.r.exp <- rbind(df.r.exp, data.frame(
            population = df.it.wl.th$pop.fac[1],          
            population.number = df.it.wl$pop.num[1],      
            temperature = df.it.wl$temperature[1],        
            well.ID = df.it.wl$well.ID[1],                
            r.exp = NA        
          ))
          
        }else{
          # Add data to our summary table
          df.r.exp <- rbind(df.r.exp, data.frame(
            population = df.it.wl.th$pop.fac[1],          # Population as factor
            population.number = df.it.wl$pop.num[1],      # Numeric population number
            temperature = df.it.wl$temperature[1],        # Temperature
            well.ID = df.it.wl$well.ID[1],                # Well ID (assuming one well ID per subset)
            r.exp = summary(r_exp)$parameters[1,1]        # The calculated r.exp value
          ))

        }
      
    }
    
  }
  
}

# Now we are going to work with the 40C data

df.rep.40 <- subset(df.rep, df.rep$temperature==40)
mat.rep.40 <- split(df.rep.40, df.rep.40$pop.num)  # Each element is a data frame for one population in df.rep

for (i in 1:length(mat.rep.40)){
  
  df.i <- subset(mat.rep.40[[i]])
  df.i <- droplevels(df.i) # Drop unused levels to isolate well replicate IDs
  
  for (w in 1:length(levels(df.i$well.ID))){
    
    df.i.wl <- subset(df.i, as.numeric(df.i$well.ID) == w)
    
    df.i.wl.late <- df.i.wl[df.i.wl$days > df.i.wl$days[which.max(df.i.wl$RFU)], ] # Do this based on the time point associated with max RFUs
    
    r_ln <- lm(logRFU~days, data= df.i.wl.late)
    
    df.r.exp <- rbind(df.r.exp, data.frame(
      population = df.i.wl.late$pop.fac[1],             # Population as factor
      population.number = df.i.wl.late$pop.num[1],      # Numeric population number
      temperature = df.i.wl.late$temperature[1],        # Temperature
      well.ID = df.i.wl.late$well.ID[1],                # Well ID (assuming one well ID per subset)
      r.exp = summary(r_ln)$coefficients[2,1]           # The calculated r.exp value
    ))
    
  }
  
}

# OK I've got my summary table! Let's save it then move on to plotting this out for a few randomly selected populations.
write.csv(df.r.exp, "data-processed/05_final_r_estimates.csv")

pops<-names(mat.rep)
ran <- sample(pops, 6, replace = F) # OK, let's select 6 random populations.

# Loop through these, generating first the trimmed data set for each population, T, and well, and then trimming based on slope thresholding
# Then fit the exponential growth model, generating a series of predicted data over a smoothed time period. Save these in a list
# Finally, plot these fitted models against the whole panel of raw data for each.

plot.list <- list() # make a list of the plots

for (i in ran){ # random populations
  
  for (t in ord.tmp){ # temperatures, now only between 10 and 34
    
    df.it <- subset(mat.rep.34[[i]], temperature==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    pred.df.mod34 <- data.frame( # Dataframe to store fitted data for the first 5 temperatures. Will reset for each T, but that's OK, as we will save the plots first
      population.fac = character(),
      population.num = numeric(),
      temperature = numeric(),
      well.ID = character(),
      smt.days = numeric(),
      fit.RFU = numeric()
    )
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      t.series <- unique(df.it.wl$days) 
      t.series <- t.series[-1] 
      
      ln.slopes <- c() 
      
      for (t in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= t, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x T level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 5%
        
        if (ln.slopes[r] != 0) { # Check for valid denominator to avoid division by zero
          # percent.chg <- (ln.slopes[s] - ln.slopes[s + 1]) / ln.slopes[s] This was the reason we were getting weird negative results! If there was a stochastic drop.
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] # This should (along with the next line) fix it.
          
          if (percent.chg >= 1.05 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { # Because now, the drop-off ignores transiently negative slopes. 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
      
      df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = 0.2), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      smt.days <- seq(min(df.it.wl.th$days, na.rm = TRUE), max(df.it.wl.th$days, na.rm = TRUE), length.out = 500) # Get a smooth distribution of time points
      
      # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
      fit.RFU.mod34 <- predict(r_exp, newdata = data.frame(days = smt.days))
      pred.df.mod34 <- rbind(pred.df.mod34, data.frame(
        population.fac = df.it.wl$pop.fac[1],
        population.num = i,
        temperature = t,
        well.ID = unique(df.it.wl$well.ID),
        smt.days = smt.days,
        fit.RFU = fit.RFU.mod34
      ))
      
      # Now we plot it out
      
      p <- ggplot() + 
        geom_point(data=df.it.wl, aes(x=days, y=RFU, colour= well.ID), size = 2) + # Observed data
        geom_line(data = pred.df.mod34, aes(x = days, y = fit.RFU, colour = well.ID), size = 1) + # Fitted line
        labs(title = paste("Population", i, "Temperature:", t), x = "Days", y = "RFU") +
        theme_minimal() +
        theme(legend.position = "none")
      
      plot.list[[paste0("Pop", df.it.wl.th$pop.fac[1], " Temp", t)]] <- p
      
    }
  }
}