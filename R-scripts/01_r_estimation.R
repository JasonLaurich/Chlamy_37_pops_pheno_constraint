# Jason Laurich
# Oct 24, 2024

################################################################################################################

# This script will upload Chlamydomonas growth (RFU) data, fit logarithmic growth curves to extract estimates of r
# for each of 37 populations at each of 6 temperatures.

######################### Upload packages ###################################

library(cowplot)
library(tidyverse)
library(nls.multstart)

######################### Upload and transform data #########################

df <- read.csv("data-processed/chlamee-acute-rfu-time.csv")
head(df) #RFU is density, days is time, temperature is temperature
str(df)

df$pop <- as.factor(df$population)
df<-subset(df,df$temperature!=20)

df$days <- df$days + 0.00000001 # Can't have 0s

str(df) # 38 populations?
levels(df$pop) # I don't recognize the cc1629 group, but we'll keep it for now

df$logRFU <- log(df$RFU) # Let's take the natural logarithm of RFU so we can fit an easy logarithmic growth curve.

df.rep <- subset(df, df$plate_type == "repeat") # We're going to work only with repeat data (this is most of the well_plates) which I'll use as replicates. 
df.rep$well.ID<-as.factor(df.rep$well_plate)

# We're also going to load Joey's exponential dataset, wherein she has trimmed RFU data to the period of exponential growth.

df.exp <- read.csv("data-processed/chlamee-acute-exponential.csv")
str(df.exp)

df.exp$days <- df.exp$days + 0.00000001 # Can't have 0s

df.exp$logRFU <- log(df.exp$RFU) # Let's take the natural logarithm of RFU so we can fit an easy logarithmic growth curve.

df.exp$well.ID<-as.factor(df.exp$well_plate) # For well-level replication

######################### Estimate maximum growth rate, r ###################

# We're going to start by fitting a logarithmic growth function to the data for one population as a proof of concept,
# then we will loop it through to calculate r.

#OK, I want to add in estimates for different wells. These should be treated like replicates, and it looks like there are ~ 4 for each pop x T.
# I think the column I want to sort by here is well_plate - for each 'well' there are repeat and single types? 
# What does repeat/single signify?

df.it <- subset(df.rep, df.rep$temperature==28 & df.rep$pop==3)
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

AIC(log_tst,log_tst_cln) #The residual-cleaned model looks better and fits better.

summary(log_tst_cln) # K is roughly the same (1764.5 to 1758.6), as is r (3.1628 to 3.182), but the cleaned model fits better.

# OK, so now we want to write a looping function.
# Create a dataframe to store growth information for each population at each temperature.
# We'll want to record estimates for K and r (with SEs)

pop <- as.vector(unique(df$pop)) # for looping through population
tmp <- as.vector(as.numeric(as.character(unique(df$temperature)))) # for looping through temperatures

temp <- vector() # data storage
pops <- vector() # data storage
well.rep <- vector() # data storage, this ID vector will be specified inside the loop.

# Create vectors for model estimates
r.log <- vector()
SE.r.log <- vector()
K.r.log <- vector()
SE.K.r.log <- vector()

#for (i in pop){
#  pops <- c(pops, rep(i,6))
#}

df.r.sum<-as.data.frame(cbind(pops,temp,well.rep, r.log, SE.r.log, K.r.log, SE.K.r.log)) # empty table with no values.

# For extra fits we will add over time, we will instead initiate a vector and cbind them later.

r.log.cln <- vector()
SE.r.log.cln <- vector()
K.r.log.cln <- vector()
SE.K.r.log.cln <- vector()

#OK, now we are going to loop

for (i in pop){ #population
  for (t in tmp){ # temperature
    
    df.it <- subset(df.rep, df.rep$temperature==t & df.rep$pop==i) # get the dataset
    df.it <- droplevels(df.it) # drop superfluous levels (to isolate well replicate IDs within each population and temperature)
    
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
      df.r.sum[nrow(df.r.sum) + 1,] = list(df.it.wl[1,9], df.it.wl[1,11], levels(df.it.wl$well.ID)[w], summary(log_r)$parameters[3,1], summary(log_r)$parameters[3,2], summary(log_r)$parameters[1,1], summary(log_r)$parameters[1,2])
      
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
      SE.r.log.cln <- c(SE.r.log.cln, summary(log_r_cln)$parameters[3,2])
      K.r.log.cln <- c(K.r.log.cln, summary(log_r_cln)$parameters[1,1])
      SE.K.r.log.cln <- c(SE.K.r.log.cln, summary(log_r_cln)$parameters[1,2])
    }
  }
}

df.r.sum<-cbind(df.r.sum, r.log.cln, SE.r.log.cln, K.r.log.cln, SE.K.r.log.cln)
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

# OK, working with the exponential growth data. I want to fit simple linear models to logged data and fit logarithmic growth curves
# to see how these approaches compare. 

# Storage vectors
r.exp.ln <- vector()
SE.r.exp.ln <- vector()

r.exp.log <- vector()
SE.r.exp.log <- vector()
K.exp.log <- vector()
SE.K.exp.log <- vector()

#OK, now we are going to loop

for (i in pop){ #population
  for (t in tmp){ # temperature
    
    df.it <- subset(df.exp, df.exp$temperature==t & df.exp$pop==i) # get the dataset
    df.it <- droplevels(df.it) # drop superfluous levels (to isolate well replicate IDs within each population and temperature)
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      ln_lin <- lm(logRFU ~ days,  #fit the model
                             data = df.it.wl)
      
      r.exp.ln <- c(r.exp.ln, summary(ln_lin)$coefficients[2,1])
      SE.r.exp.ln <- c(SE.r.exp.ln, summary(ln_lin)$coefficients[2,2])
      
      log_exp <- nls_multstart(RFU ~ K / (1 + ((K - N.0) / N.0) * exp(-r * days)),  #changed to N.0 because N0 is in the dataframe
                             data = df.it.wl,
                             start_lower = c(K = max(df.it.wl$RFU)*0.75, N.0 = 1, r = 0.2), 
                             start_upper = c(K = max(df.it.wl$RFU)*1.25, N.0 = 50, r = 3.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      r.exp.log <- c(r.exp.log, summary(log_exp)$parameters[3,1])
      SE.r.exp.log <- c(SE.r.exp.log, summary(log_exp)$parameters[3,2])
      K.exp.log <- c(K.exp.log, summary(log_exp)$parameters[1,1])
      SE.K.exp.log <- c(SE.K.exp.log, summary(log_exp)$parameters[1,2])
    }
  }
}

df.r.sum<-cbind(df.r.sum, r.exp.ln, SE.r.exp.ln, r.exp.log, SE.r.exp.log, K.exp.log, SE.K.exp.log)

######################### Fit TPCs using nls ################################



######################### Fit TPCs using JAGS ###############################