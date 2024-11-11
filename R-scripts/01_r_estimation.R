# Jason Laurich
# Oct 24, 2024

################################################################################################################

# This script will upload Chlamydomonas growth (RFU) data, fit logistic growth curves to extract estimates of r
# for each of 37 populations at each of 6 temperatures.

######################### Upload packages ###################################

library(cowplot)
library(tidyverse)
library(nls.multstart)
library(MuMIn) #AICc
library(gridExtra)
library(ggplot2)

######################### Upload and transform data #########################

df <- read.csv("data-processed/chlamee-acute-rfu-time.csv")
head(df) #RFU is density, days is time, temperature is temperature
str(df)

df<-df[,-c(1,2,4,5,6,8,13)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df<-subset(df,df$temperature!=20)

df$days <- df$days + 0.001 # Can't have 0s
df$logRFU <- log(df$RFU + 0.001) # Let's take the natural logarithm of RFU so we can fit an easy logistic growth curve.

levels(df$pop.fac) # I don't recognize the cc1629 group, but we'll keep it for now

df.rep <- subset(df, df$plate_type == "repeat") # We're going to work only with repeat data (this is most of the well_plates) which I'll use as replicates. 
df.rep$well.ID<-as.factor(df.rep$well_plate)

# We're also going to load Joey's exponential dataset, wherein she has trimmed RFU data to the period of exponential growth.

df.exp <- read.csv("data-processed/chlamee-acute-exponential.csv")
str(df.exp)

df.exp$days <- df.exp$days + 0.001 # Can't have 0s
df.exp$logRFU <- log(df.exp$RFU + 0.001) # Let's take the natural logarithm of RFU so we can fit an easy logistic growth curve.

df.exp<-df.exp[,-c(1,2,4,5,6,8,13,17,18)]

df.exp$well.ID<-as.factor(df.exp$well_plate) # For well-level replication
df.exp$pop.fac <- as.factor(df.exp$population)
df.exp$pop.num <- as.numeric(df.exp$pop.fac)

# Split df.rep and df.exp by population
mat.rep <- split(df.rep, df.rep$pop.num)  # Each element is a data frame for one population in df.rep
mat.exp <- split(df.exp, df.exp$pop.num)  # Each element is a data frame for one population in df.exp


######################### Estimate maximum growth rate, r ###################

# We're going to start by fitting a logistic growth function to the data for one population as a proof of concept,
# then we will loop it through to calculate r.

#OK, I want to add in estimates for different wells. These should be treated like replicates, and it looks like there are ~ 4 for each pop x T.
# I think the column I want to sort by here is well_plate - for each 'well' there are repeat and single types? 
# What does repeat/single signify?

df.it <- subset(mat.rep[[3]], temperature==28)
df.it <- droplevels(df.it)

df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == 1)

log_tst <- nls_multstart(RFU ~ K / (1 + ((K - N0) / N0) * exp(-r * days)),  
                          data = df.it.wl,
                          start_lower = c(K = max(df.it.wl$RFU)*0.75, N0 = 1, r = 0.2),  # We'll base K on observed metrics so we can adjust it for each population x temperature
                          start_upper = c(K = max(df.it.wl$RFU)*1.25, N0 = 50, r = 2.5),   # Refined upper bounds for faster growth
                          iter = 500,
                          supp_errors = 'Y',
                          control = nls.control(maxiter = 200))

summary(log_tst)

# OK let's plot this, starting by generating a smoothed x axis
days.smth <- seq(min(df.it.wl$days, na.rm = TRUE),
                   max(df.it.wl$days, na.rm = TRUE), 
                   length.out = 500)

fit.smth <- predict(log_tst, newdata = data.frame(days = days.smth))

plt.df.smth <- data.frame(
  days = days.smth,
  fits = fit.smth
)

# Plot the raw data and fitted curve
ggplot() +
  geom_point(data = df.it.wl, aes(x = days, y = RFU), colour = "magenta4", size = 2, alpha = 0.6) +  # Raw data
  geom_line(data = plt.df.smth, aes(x = days.smth, y = fit.smth), colour = "darkorange2", linewidth = 1) +    # Smoothed logistic growth fit
  labs(title = "Representative logistic growth curve (Population 3, 28C)",
       x = "Days",
       y = "RFU") +
  theme_cowplot()

# OK, this looks great, but I want to check for wonky data points across the board.

# I think I'll tackle this by getting the residuals, and then filtering out extreme data before re-fitting the model. We'll keep both estimates for r when we loop.
df.it.wl$residuals <- residuals(log_tst)

# For now let's set the residual threshold at 2 sds
res.lim <- 2 * sd(df.it.wl$residuals)
df.cln <- subset(df.it.wl, abs(residuals) < res.lim)

# Now we'll refit the logistic model without the outliers
log_tst_cln <- nls_multstart(RFU ~ K / (1 + ((K - N0) / N0) * exp(-r * days)),  
                               data = df.cln,
                               start_lower = c(K = max(df.cln$RFU)*0.75, N0 = 1, r = 0.2),  
                               start_upper = c(K = max(df.cln$RFU)*1.25, N0 = 50, r = 2.5),  
                               iter = 500,
                               supp_errors = 'Y',
                               control = nls.control(maxiter = 200))

days.smth.cln <- seq(min(df.cln$days, na.rm = TRUE),
                 max(df.it$days, na.rm = TRUE), 
                 length.out = 500)

fit.smth.cln <- predict(log_tst_cln, newdata = data.frame(days = days.smth))

plt.df.smth.cln <- data.frame(
  days = days.smth.cln,
  fits = fit.smth.cln
)

# Plot the cleaned data and curve
ggplot() +
  geom_point(data = df.cln, aes(x = days, y = RFU), colour = "magenta4", size = 2, alpha = 0.6) +  # Raw data
  geom_line(data = plt.df.smth.cln, aes(x = days.smth.cln, y = fit.smth.cln), colour = "darkorange2", linewidth = 1) +    # Smoothed logistic growth fit
  labs(title = "Representative logistic growth curve without outliers (Population 3, 28C)",
       x = "Days",
       y = "RFU") +
  theme_cowplot()

AICc(log_tst,log_tst_cln) #The residual-cleaned model looks better and fits better.

summary(log_tst_cln) # K is roughly the same (1764.5 to 1758.6), as is r (3.1628 to 3.182), but the cleaned model fits better.

# OK, so now we want to write a looping function.
# Create a dataframe to store growth information for each population at each temperature.

tmp <- as.vector(as.numeric(as.character(unique(df$temperature)))) # for looping through temperatures
ord.tmp<- sort(tmp)

df.r.sum <- data.frame( # dataframe for storage
  population = character(),
  population.number = numeric(),
  temperature = numeric(),
  well.ID = character(),
  r.log = numeric()
)

# For extra fits we will add over time, we will instead initiate a vector and cbind them later.

r.log.cln <- vector()

#OK, now we are going to loop

for (i in 1:38){ #population
  for (t in ord.tmp){ # temperature
    
    df.it <- subset(mat.rep[[i]], temperature==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      log_r <- nls_multstart(RFU ~ K / (1 + ((K - N0) / N0) * exp(-r * days)),  #fit the model
                             data = df.it.wl,
                             start_lower = c(K = max(df.it.wl$RFU)*0.75, N0 = 1, r = 0.2), 
                             start_upper = c(K = max(df.it.wl$RFU)*1.25, N0 = 50, r = 3.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      # Add data to our summary table
      df.r.sum <- rbind(df.r.sum, data.frame(
        population = df.it.wl$pop.fac[1],             # Population as factor
        population.number = df.it.wl$pop.num[1],      # Numeric population number
        temperature = df.it.wl$temperature[1],        # Temperature
        well.ID = df.it.wl$well.ID[1],                # Well ID (assuming one well ID per subset)
        r.log = summary(log_r)$parameters[3,1]        # The calculated r.log value
      ))
      
      df.it.wl$residuals <- residuals(log_r) # Get the residuals so we can filter out wonky data. This is especially important here
      # Because high T growth curves spike super fast, then drops off (death?)
      
      # As before, let's threshold at 2 sds for now.
      res.lim <- 2 * sd(df.it.wl$residuals)
      df.cln <- subset(df.it.wl, abs(residuals) < res.lim)
      
      # Now we'll refit the logistic model without the outliers
      
      log_r_cln <- nls_multstart(RFU ~ K / (1 + ((K - N0) / N0) * exp(-r * days)),  
                                 data = df.cln,
                                 start_lower = c(K = max(df.cln$RFU)*0.75, N0 = 1, r = 0.2),  
                                 start_upper = c(K = max(df.cln$RFU)*1.25, N0 = 50, r = 3.5),  
                                 iter = 500,
                                 supp_errors = 'Y',
                                 control = nls.control(maxiter = 200))
      
      r.log.cln <- c(r.log.cln, summary(log_r_cln)$parameters[3,1])
    }
  }
}

df.r.sum$r.log.cln <- r.log.cln

# Without even graphing things, we can see that the problem is including data where populations are crashing after they have peaked.
# We would need to trim this data off, say by adding a 10% time buffer zone after the highest r has been recorded for each population.
# This is probably what Joey did to generate the exponential growth data. 
# The problem is especially acute at high temperatures, which fits the above hypothesis.

# OK I've looked at Joey's original code, it seems like she is treating the early spike in RFU under high T (40)
# as not real (ie. for her, days must be > 1 for exponential = T).

# Her code for sorting out what is exponential is as follows:
# mutate(exponential = case_when(temperature == 28 & days < 1 ~ "yes",
                  #             temperature == 34 & days < 1 ~ "yes",
                  #             temperature == 22 & days < 2 ~  "yes",
                  #             temperature == 10 & days < 7 ~  "yes",
                  #             temperature == 16 & days < 3 ~  "yes",
                  #             temperature == 40 & days < 7 & days > 1 ~  "yes",
                  #             TRUE  ~ "no")) %>% 

# After consulting with her, I'm going to (A) use her data, and (B), write my own script that iterates across populations and
# temperatures to identify where the growth starts to go down.
# For B, I want to estimate R both (1) at the 1st spike, and (2) during the declining phase as she did
# I actually think I want to threshold based on RFU, not days. So for example, exclude points after RFUs drop below ~25% of peak values.

# OK, working with the exponential growth data. I want to fit simple linear models to logged data and fit logistic growth curves
# to see how these approaches compare. 

# Storage vectors
r.exp.ln <- vector()

#OK, now we are going to loop

for (i in 1:38){ #population
  for (t in ord.tmp){ # temperature
    
    df.it <- subset(mat.exp[[i]], temperature==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      ln_lin <- lm(logRFU ~ days,  #fit simple linear model
                             data = df.it.wl)
      
      r.exp.ln <- c(r.exp.ln, summary(ln_lin)$coefficients[2,1])
    }
  }
}

df.r.sum$r.exp.ln <- r.exp.ln

# Logistic growth models will not work for the 40C group in the exp dataset. This is due to the fact that there
# are only ~3 data points for each time series, which are either decreasing or flat.

r.exp.log <- vector()

for (i in 1:38){ #population
  for (t in ord.tmp){ # temperature
    
    df.it <- subset(mat.exp[[i]], temperature==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      if (t!=40){
        
        log_exp <- nls_multstart(RFU ~ K / (1 + ((K - N.0) / N.0) * exp(-r * days)),  #changed to N.0 because N0 is in the dataframe
                                 data = df.it.wl,
                                 start_lower = c(K = max(df.it.wl$RFU)*0.75, N.0 = 1, r = 0.2), 
                                 start_upper = c(K = max(df.it.wl$RFU)*1.25, N.0 = 50, r = 3.5),   
                                 iter = 500,
                                 supp_errors = 'Y',
                                 control = nls.control(maxiter = 200))
        
        r.exp.log <- c(r.exp.log, summary(log_exp)$parameters[3,1])
        
        }else {
          r.exp.log <- c(r.exp.log, NA)
        }
    }
  }
}

df.r.sum$r.exp.log <- r.exp.log

# So I think for all of the non-40s, I like the fitting of logistic growth curves best.
# For the 40's, I guess we will now work with the ascending and descending phase?

# So for the 40C data, we'll fit logistic growth curves and linear regressions to the early data (days < 1), but we'll
# threshold based on max RFU count. 

# This will not give me enough datapoints to fit a logistic growth curve, so I think I'll just estimate r based on a linear 
# regression of logged RFU count

r.40.ln.early <- vector()
r.40.ln.full <- vector()

for (i in 1:38){ #population
  for (t in ord.tmp){ # temperature
    
    df.it <- subset(mat.rep[[i]], temperature==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      if (t==40){
      
        df.it.wl.early <- subset(df.it.wl, df.it.wl$days <= df.it.wl$days[which.max(df.it.wl$RFU)])
        
        ln_r_40_early <- lm(logRFU~days, data=df.it.wl.early) # Constrain data to early, increasing phase.
        
        r.40.ln.early <- c(r.40.ln.early, summary(ln_r_40_early)$coefficients[2,1])
        
        ln_r_40_full <- lm(logRFU~days, data=df.it.wl) # Fit to the whole dataset (linear regression)
        
        r.40.ln.full <- c(r.40.ln.full, summary(ln_r_40_full)$coefficients[2,1])
        
      }else {
        
        r.40.ln.early <- c(r.40.ln.early, NA)
        
        ln_r_40_full <- lm(logRFU~days, data=df.it.wl) # Fit to the whole dataset (linear regression)
        
        r.40.ln.full <- c(r.40.ln.full, summary(ln_r_40_full)$coefficients[2,1])
        
      }
    }
  }
}

df.r.sum$r.40.ln.early <- r.40.ln.early
df.r.sum$r.40.ln.full <- r.40.ln.full

write.csv(df.r.sum, "data-processed/01_pop_well_r_estimates.csv")

# OK, now for 3 populations, let's fit and plot models for all growth estimation approaches. Let's plot the raw data and model
# fits, and compare models using AIC.

pops<-names(mat.rep)
ran <- sample(pops, 3, replace = F) #Let's randomly select 3 populations to do a full modelling and plotting markup here. 

aic.results <- data.frame( # Build a dataframe to hold key information for each temperature, well, and model
  population.fac = character(),
  population.num = numeric(),
  temperature = numeric(),
  well.ID = character(),
  model = character(),
  AIC = numeric()
  )

AIC <- vector() # for storing additional aic scores which we will cbind to the original dataframe. 

pred.df.mod1 <- data.frame( # Dataframe to store fitted data for the first model
  population.fac = character(),
  population.num = numeric(),
  temperature = numeric(),
  well.ID = character(),
  smt.days = numeric(),
  fit.RFU = numeric()
)

pred.df.mod2 <- data.frame( # Dataframe to store fitted data for the second model
  population.fac = character(),
  population.num = numeric(),
  temperature = numeric(),
  well.ID = character(),
  smt.days = numeric(),
  fit.RFU = numeric()
)


for (i in ran){
  
  for (t in tmp){
    
    df.it <- subset(mat.rep[[i]], temperature==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      # OK, a sequence of days against which we will extract predicted model fits.
      smt.days <- seq(min(df.it.wl$days, na.rm = TRUE), max(df.it.wl$days, na.rm = TRUE), length.out = 500)
      
      # Fit model 1: logistic growth model
      log_r <- nls_multstart(RFU ~ K / (1 + ((K - N0) / N0) * exp(-r * days)),  #fit the model
                             data = df.it.wl,
                             start_lower = c(K = max(df.it.wl$RFU)*0.75, N0 = 1, r = 0.2), 
                             start_upper = c(K = max(df.it.wl$RFU)*1.25, N0 = 50, r = 3.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      # Store model 1 AIC
      aic.results <- rbind(aic.results, data.frame(
        population.fac = df.it.wl$pop.fac[1],
        population.num = i,
        temperature = t,
        well.ID = unique(df.it.wl$well.ID),
        model = "Model 1 (Logistic)",
        AIC = AICc(log_r)
      ))
      
      # Predict fitted RFU values for Model 1 and store
      fit.RFU.mod1 <- predict(log_r, newdata = data.frame(days = smt.days))
      pred.df.mod1 <- rbind(pred.df.mod1, data.frame(
        population.fac = df.it.wl$pop.fac[1],
        population.num = i,
        temperature = t,
        well.ID = unique(df.it.wl$well.ID),
        smt.days = smt.days,
        fit.RFU = fit.RFU.mod1
      ))
      
      # Threshold based on residuals
      df.it.wl$residuals <- residuals(log_r)
      res.lim <- 2 * sd(df.it.wl$residuals)
      df.cln <- subset(df.it.wl, abs(residuals) < res.lim)
      
      # Now we'll refit the logistic model without the outliers
      
      log_r_cln <- nls_multstart(RFU ~ K / (1 + ((K - N0) / N0) * exp(-r * days)),  
                                 data = df.cln,
                                 start_lower = c(K = max(df.cln$RFU)*0.75, N0 = 1, r = 0.2),  
                                 start_upper = c(K = max(df.cln$RFU)*1.25, N0 = 50, r = 3.5),  
                                 iter = 500,
                                 supp_errors = 'Y',
                                 control = nls.control(maxiter = 200))
      
      # Store model 2 AIC. We'll need a separate 2 column data frame which we can cbind after
      AIC <- c(AIC, AICc(log_r_cln))
      
      # Predict fitted RFU values for Model 2 and store
      fit.RFU.mod2 <- predict(log_r_cln, newdata = data.frame(days = smt.days))
      pred.df.mod2 <- rbind(pred.df.mod2, data.frame(
        population.fac = df.it.wl$pop.fac[1],
        population.num = i,
        temperature = t,
        well.ID = unique(df.it.wl$well.ID),
        smt.days = smt.days,
        fit.RFU = fit.RFU.mod2
      ))
      
    }
  }
}

aic.results$model2 <- rep("Model 2 (Trimmed logistic)", nrow(aic.results))
aic.results$AIC2 <- AIC

# OK, now we are going to fit linear models to the growth data from the exponential phase, storing fitted data and AIC values

pred.df.mod3 <- data.frame( # Dataframe to store fitted data for the third model
  population.fac = character(),
  population.num = numeric(),
  temperature = numeric(),
  well.ID = character(),
  smt.days = numeric(),
  fit.RFU = numeric()
)

AIC3 <- vector()

for (i in ran){
  
  for (t in tmp){
    
    df.it <- subset(mat.exp[[i]], temperature==t) # sort by t and pop
    df.it <- droplevels(df.it) # drop levels
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      # OK, a sequence of days against which we will extract predicted model fits.
      smt.days <- seq(min(df.it.wl$days, na.rm = TRUE), max(df.it.wl$days, na.rm = TRUE), length.out = 500)
      
      # Fit model 3: linear growth model on logged data
      ln_lin <- lm(logRFU ~ days,  #fit simple linear model
                   data = df.it.wl)
      
      # Store model 3 AIC
      AIC3 <- c(AIC3, AICc(ln_lin))
      
      # Predict fitted RFU values for Model 3 and store
      fit.RFU.mod3 <- predict(ln_lin, newdata = data.frame(days = smt.days))
      pred.df.mod3 <- rbind(pred.df.mod3, data.frame(
        population.fac = df.it.wl$pop.fac[1],
        population.num = i,
        temperature = t,
        well.ID = unique(df.it.wl$well.ID),
        smt.days = smt.days,
        fit.RFU = fit.RFU.mod3
      ))
      
    }
  }
}

aic.results$model3 <- rep("Model 3 (Linear logged)", nrow(aic.results))
aic.results$AIC3 <- AIC3

# Now I'm going to fit the 4th model (logistic growth curve on data limited to the exponential growth period) to all <40C data

pred.df.mod4 <- data.frame( # Dataframe to store fitted data for the fourth model
  population.fac = character(),
  population.num = numeric(),
  temperature = numeric(),
  well.ID = character(),
  smt.days = numeric(),
  fit.RFU = numeric()
)

AIC4 <- vector()

for (i in ran){
  
  for (t in tmp){
    
    df.it <- subset(mat.exp[[i]], temperature==t) # sort by t and pop
    df.it <- droplevels(df.it) # drop levels
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      if (t!=40){
        
        # OK, a sequence of days against which we will extract predicted model fits.
        smt.days <- seq(min(df.it.wl$days, na.rm = TRUE), max(df.it.wl$days, na.rm = TRUE), length.out = 500)
        
        log_exp <- nls_multstart(RFU ~ K / (1 + ((K - N.0) / N.0) * exp(-r * days)),  #changed to N.0 because N0 is in the dataframe
                                 data = df.it.wl,
                                 start_lower = c(K = max(df.it.wl$RFU)*0.75, N.0 = 1, r = 0.2), 
                                 start_upper = c(K = max(df.it.wl$RFU)*1.25, N.0 = 50, r = 3.5),   
                                 iter = 500,
                                 supp_errors = 'Y',
                                 control = nls.control(maxiter = 200))
        
        # Store model 4 AIC.
        AIC4 <- c(AIC4, AICc(log_exp))
        
        # Predict fitted RFU values for Model 4 and store
        fit.RFU.mod4 <- predict(log_exp, newdata = data.frame(days = smt.days))
        pred.df.mod4 <- rbind(pred.df.mod4, data.frame(
          population.fac = df.it.wl$pop.fac[1],
          population.num = i,
          temperature = t,
          well.ID = unique(df.it.wl$well.ID),
          smt.days = smt.days,
          fit.RFU = fit.RFU.mod4
        ))
        
        
      }else {
        
        AIC4 <- c(AIC4, NA) # Won't work for t == 40
        
        pred.df.mod4 <- rbind(pred.df.mod4, data.frame(
          population.fac = df.it.wl$pop.fac[1],
          population.num = i,
          temperature = t,
          well.ID = unique(df.it.wl$well.ID),
          smt.days = NA,
          fit.RFU = NA
        ))

      }
    }
  }
}

aic.results$model4 <- rep("Model 4 (Logistic growth, exponential phase)", nrow(aic.results))
aic.results$AIC4 <- AIC4

# OK now we'll run the 5th model, which fits lineral models to logged data for only the acsending portion of the 40C treatments.

pred.df.mod5 <- data.frame( # Dataframe to store fitted data for the fifth model
  population.fac = character(),
  population.num = numeric(),
  temperature = numeric(),
  well.ID = character(),
  smt.days = numeric(),
  fit.RFU = numeric()
)

AIC5 <- vector()

for (i in ran){
  
  for (t in tmp){
    
    df.it <- subset(mat.rep[[i]], temperature==t) # sort by t and pop
    df.it <- droplevels(df.it) # drop levels
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      if (t==40){
        
        df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
        
        df.it.wl.early <- subset(df.it.wl, df.it.wl$days <= df.it.wl$days[which.max(df.it.wl$RFU)])
        
        smt.days <- seq(min(df.it.wl.early$days, na.rm = TRUE), max(df.it.wl.early$days, na.rm = TRUE), length.out = 500)
        
        ln_r_40_early <- lm(logRFU~days, data=df.it.wl.early) # Constrain data to early, increasing phase.
        
        # Store model 5 AIC.
        AIC5 <- c(AIC5, AICc(ln_r_40_early))
        
        # Predict fitted RFU values for Model 5 and store
        fit.RFU.mod5 <- predict(ln_r_40_early, newdata = data.frame(days = smt.days))
        pred.df.mod5 <- rbind(pred.df.mod5, data.frame(
          population.fac = df.it.wl$pop.fac[1],
          population.num = i,
          temperature = t,
          well.ID = unique(df.it.wl$well.ID),
          smt.days = smt.days,
          fit.RFU = fit.RFU.mod5
        ))
        
      }else {
        
        AIC5 <- c(AIC5, NA) # Don't use for other temperatures
        
        pred.df.mod5 <- rbind(pred.df.mod5, data.frame(
          population.fac = df.it.wl$pop.fac[1],
          population.num = i,
          temperature = t,
          well.ID = unique(df.it.wl$well.ID),
          smt.days = NA,
          fit.RFU = NA
        ))
        
      }
    }
  }
}

aic.results$model5 <- rep("Model 5 (Log linear, early phase, 40C)", nrow(aic.results))
aic.results$AIC5 <- AIC5

# OK, onto the sixth and final model! Just fitting a logged linear regression to the full dataset for all temperatures.

pred.df.mod6 <- data.frame( # Dataframe to store fitted data for the sixth model
  population.fac = character(),
  population.num = numeric(),
  temperature = numeric(),
  well.ID = character(),
  smt.days = numeric(),
  fit.RFU = numeric()
)

AIC6 <- vector()

for (i in ran){
  
  for (t in tmp){
    
    df.it <- subset(mat.rep[[i]], temperature==t)
    df.it <- droplevels(df.it) # drop levels
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      # OK, a sequence of days against which we will extract predicted model fits.
      smt.days <- seq(min(df.it.wl$days, na.rm = TRUE), max(df.it.wl$days, na.rm = TRUE), length.out = 500)
      
      # Fit model 6: linear growth model on logged data across the full range
      ln_lin <- lm(logRFU ~ days,  #fit simple linear model
                   data = df.it.wl)
      
      # Store model 6 AIC
      AIC6 <- c(AIC6, AICc(ln_lin))
      
      # Predict fitted RFU values for Model 6 and store
      fit.RFU.mod6 <- predict(ln_lin, newdata = data.frame(days = smt.days))
      pred.df.mod6 <- rbind(pred.df.mod6, data.frame(
        population.fac = df.it.wl$pop.fac[1],
        population.num = i,
        temperature = t,
        well.ID = unique(df.it.wl$well.ID),
        smt.days = smt.days,
        fit.RFU = fit.RFU.mod6
      ))
      
    }
  }
}

aic.results$model6 <- rep("Model 6 (Linear logged, full data)", nrow(aic.results))
aic.results$AIC6 <- AIC6

write.csv(aic.results, "data-processed/02_sample_aic_results_r_estimation.csv") # AIC table for model comparison.

#############################################

# OK, now we want to generate some summary figures to see how these various models are fitting the data.

# Set the output directory
output_dir <- "figures" 

# Function to generate each model-temperature plot
generate_model_plot <- function(i, t, model_num) {
  # Select the appropriate raw data source (mat.rep or mat.exp) based on the model number
  if (model_num %in% c(1, 2, 5, 6)) {
    df.raw <- subset(mat.rep[[i]], temperature == t)  # Use mat.rep for models 1, 2, 5, and 6
  } else if (model_num %in% c(3, 4)) {
    df.raw <- subset(mat.exp[[i]], temperature == t)  # Use mat.exp for models 3 and 4
  }
  
  # Set y-axis variable based on whether the data should be logged
  y_var <- if (model_num %in% c(3, 5, 6)) "logRFU" else "RFU"
  
  # Retrieve the complete predictive data frame for the model
  df.pred <- get(paste0("pred.df.mod", model_num))
  
  # Subset predictions based on the current temperature `t` and population `i`
  df.pred <- subset(df.pred, temperature == t & population.num == i)
  df.pred$well.ID <- factor(df.pred$well.ID, levels = unique(df.pred$well.ID))
  
  # Skip specific model-temperature combinations by creating blank plots
  if ((model_num == 4 && t == 40) || (model_num == 5 && t != 40)) {
    return(ggplot() + theme_void())
  } else {
    # Generate the actual plot
    return(ggplot() +
             geom_point(data = df.raw, aes(x = days, y = !!sym(y_var), color = well.ID), alpha = 0.5) +
             geom_line(data = df.pred, aes(x = smt.days, y = fit.RFU, color = well.ID)) +
             labs(title = paste("Temperature:", t, "- Model", model_num), x = "Days", y = y_var) +
             theme_minimal() +
             theme(legend.position = "none"))
  }
}

# Function to arrange the 6x6 grid and save as PDF for each population
generate_and_save_plot <- function(i, population_label) {
  plot_list <- list()
  for (model_num in 1:6) {
    for (t in ord.tmp) {
      plot_list[[paste0("model", model_num, "_temp", t)]] <- generate_model_plot(i, t, model_num)
    }
  }
  
  # Arrange the plots in a 6x6 grid
  combined_plot <- arrangeGrob(
    grobs = plot_list,
    nrow = 6,
    ncol = 6,
    top = paste("Population:", population_label)
  )
  
  # Save the plot as PDF in the specified directory
  ggsave(
    filename = file.path(output_dir, paste0("01_pop", population_label, "_ex_r_modfits.pdf")),
    plot = combined_plot,
    width = 16,
    height = 12
  )
}

# Loop over the selected populations to create and save each PDF
for (i in ran) {
  population_label <- as.character(unique(mat.rep[[i]]$pop.fac)[1])  # Get character label for current population
  generate_and_save_plot(i, population_label)
}

