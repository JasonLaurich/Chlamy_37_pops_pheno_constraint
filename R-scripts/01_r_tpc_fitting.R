# Jason Laurich
# Oct 24, 2024

################################################################################################################

# This script will upload Chlamydomonas growth (RFU) data, fit logarithmic growth curves to extract estimates of r
# for each of 37 populations at each of 6 temperatures, then fit TPCs to the data using nls and bayesian methods

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

# We're also going to load Joey's exponential dataset, wherein she has trimmed RFU data to the period of exponential growth.


######################### Estimate maximum growth rate, r ###################

# We're going to start by fitting a logarithmic growth function to the data for one population as a proof of concept,
# then we will loop it through to calculate r.

df.it <- subset(df, df$temperature==28 & df$pop==3)

log_tst <- nls_multstart(RFU ~ K / (1 + ((K - N0) / N0) * exp(-r * days)),  
                          data = df.it,
                          start_lower = c(K = max(df.it$RFU)*0.75, N0 = 1, r = 0.2),  # We'll base K on observed metrics so we can adjust it for each population x temperature
                          start_upper = c(K = max(df.it$RFU)*1.25, N0 = 50, r = 2.5),   # Refined upper bounds for faster growth
                          iter = 500,
                          supp_errors = 'Y',
                          control = nls.control(maxiter = 200))

summary(log_tst)

# First, generate a sequence of days for smooth curve plotting
days.smth <- seq(min(df.it$days, na.rm = TRUE),
                   max(df.it$days, na.rm = TRUE), 
                   length.out = 500)

fit.smth <- predict(log_tst, newdata = data.frame(days = days.smth))

plt.df.smth <- data.frame(
  days = days.smth,
  fits = fit.smth
)

# Plot the raw data and fitted curve
ggplot() +
  geom_point(data = df.it, aes(x = days, y = RFU), colour = "magenta4", size = 2, alpha = 0.6) +  # Raw data
  geom_line(data = plt.df.smth, aes(x = days.smth, y = fit.smth), colour = "darkorange2", size = 1) +    # Smoothed logistic growth fit
  labs(title = "Representative logistic growth curve (Population 3, 28C)",
       x = "Days",
       y = "RFU") +
  theme_cowplot()

# OK, this looks great, but I want to check for wonky data points across the board.

# I think I'll tackle this by getting the residuals, and then filtering out extreme data before re-fitting the model. We'll keep both estimates for r when we loop.
df.it$residuals <- residuals(log_tst)

# For now let's set the residual threshold at 2 sds
res.lim <- 2 * sd(df.it$residuals)
df.cln <- subset(df.it, abs(residuals) < res.lim)

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
  geom_line(data = plt.df.smth.cln, aes(x = days.smth.cln, y = fit.smth.cln), colour = "darkorange2", size = 1) +    # Smoothed logistic growth fit
  labs(title = "Representative logistic growth curve without outliers (Population 3, 28C)",
       x = "Days",
       y = "RFU") +
  theme_cowplot()

AIC(log_tst,log_tst_cln) #The residual-cleaned model looks better and fits better.

summary(log_tst_cln) # K has increased (from 1807 to 1845, and r has decreased from 2.465 to 2.372)
# Is that right?


# OK, so now we want to write a looping function.
# Create a dataframe to store growth information for each population at each temperature.
# We'll want to record estimates for K and r (with SEs)

tmp <- as.vector(as.numeric(as.character(unique(df$temperature)))) # for looping through temperatures

pop <- as.vector(unique(df$pop))
temp <- rep(tmp, length(pop)) 
pops <- vector()

for (i in pop){
  pops <- c(pops, rep(i,6))
}

df.r.sum<-as.data.frame(cbind(pops,temp))

# Create vectors for fits to record.
r.log <- vector()
SE.r.log <- vector()
K.r.log <- vector()
SE.K.r.log <- vector()

r.log.cln <- vector()
SE.r.log.cln <- vector()
K.r.log.cln <- vector()
SE.K.r.log.cln <- vector()


#OK, now we are going to loop

for (i in pop){ #population
  for (t in tmp){ # temperature
    df.it <- subset(df, df$temperature==t & df$pop==i) # get the dataset
    
    log_r <- nls_multstart(RFU ~ K / (1 + ((K - N0) / N0) * exp(-r * days)),  #fit the model
                             data = df.it,
                             start_lower = c(K = max(df.it$RFU)*0.75, N0 = 1, r = 0.2), 
                             start_upper = c(K = max(df.it$RFU)*1.25, N0 = 50, r = 3.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
    
    r.log <- c(r.log, summary(log_r)$parameters[3,1])
    SE.r.log <- c(SE.r.log, summary(log_r)$parameters[3,2])
    K.r.log <- c(K.r.log, summary(log_r)$parameters[1,1])
    SE.K.r.log <- c(SE.K.r.log, summary(log_r)$parameters[1,2])
    
    df.it$residuals <- residuals(log_r) # Get the residuals so we can filter out wonky data. This is especially important here
    # Because high T growth curves spike super fast, then drops off (death?)
    
    # For now let's set the residual threshold at 2 sds
    res.lim <- 2 * sd(df.it$residuals)
    df.cln <- subset(df.it, abs(residuals) < res.lim)
    
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

df.r.sum<-cbind(df.r.sum, r.log, r.log.cln)
# Without even graphing things, we can see that the problem is including data where populations are crashing after they have peaked.
# We would need to trim this data off, say by adding a 10% time buffer zone after the highest r has been recorded for each population.
# This is probably what Joey did to generate the exponential growth data. 
# The problem is especially acute at high temperatures, which fits the above hypothesis.

# OK I've looked at Joey's original code, it seems like she is treating the early spike in RFU under high T (40)
# as not real (ie. for her, days must be > 1 for exponential = T).

# I need to talk with her about that, for now I'm going to use her data. 

######################### Fit TPCs using nls ################################



######################### Fit TPCs using JAGS ###############################