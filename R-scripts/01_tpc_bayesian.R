# Jason Laurich
# Oct 21, 2024

# This script will fits TPCs to data from Joey Bernhardt's postdoctoral research using the
# bayesTPC package.

# Uploading packages

library(nimble)
library(HDInterval)
library(MCMCvis)
library(coda) # makes diagnostic plots
library(IDPmisc) # makes nice colored pairs plots to look at joint posteriors
library(matrixStats)
library(truncnorm)
#devtools::install_github("johnwilliamsmithjr/bayesTPC") # If needed
library(bayesTPC)
library(cowplot)
library(tidyverse)
library(growthrates) #one of the methods I'll use to estimate r
library(brms)
library(nls.multstart)

###########################################################################################

# Importing, interrogating, and re-arranging the data

data <- read.csv("data-processed/globe-chlamy-exponential-RFU-time.csv")
head(data)
str(data)

length(unique(data$population))
# 15 populations, we'll want to fit separate TPCs to each one. I imagine this will be easiest if we turn this dataframe into a matrix

#Remove any missing data
data <- subset(data, !is.na(data$RFU))

#OK, let's create a matrix to hold the data for each population (1 to 15 here)
pops<-as.data.frame(1:15)
mat <- lapply(1:nrow(pops), function(z) data[which(data$population==pops[z,1]),])

#Check it
for (i in c(5,10)){
  print(head(mat[[i]]))
}
#Looks good!

###########################################################################################

# Let's look at a few plots!

#This works with the data in matrix form
ggplot(aes(x = days, y = RFU, colour = factor(temperature)), data=mat[[3]]) + geom_point()
ggplot(aes(x = days, y = RFU, colour = factor(temperature)), data=mat[[9]]) + geom_point()

#So does this
for (i in 1:2){
  print(ggplot(aes(x = days, y = RFU, colour = factor(temperature)), data=mat[[i]]) + geom_point())
}

#OK, the plots look fine!

###########################################################################################

# We need to calculate the growth rate here first. This will be on our y-axis. For each population at each T,
# We need a measurement for replicate (ie. we are collapsing the temporal axis)
# I imagine we take a measurement of growth during the exponential phase here?

#OK, I want to try fitting growth rate a couple of ways

#Let's isolate data for a single Chlamy population at a single temperature. We'll try to 
#estimate growth and fit TPC's for that dataset before we scale up.

df3_34<-subset(mat[[3]],mat[[3]]$temperature==34)
ggplot(aes(x = days, y = RFU), data=df3_34) + geom_point()

# I believe for the spline fitting, all 0's must be removed, and there is an NA row at the end with days=0.
df3_34<-subset(df3_34, df3_34$days!=0)

#Let's do this for mat[[3]] as well

###########################################################################################

#1. nls

###########################################################################################

#2. growthrate package

# I am running the following code according to the following tutorial:
# https://tpetzoldt.github.io/growthrates/doc/Introduction.html#easy-linear-method

# We are going to try a logistic growth model

p     <- c(y0 = 0.01, mumax = 10, K = 600)
lower <- c(y0 = 1e-6, mumax = 0,   K = 0)
upper <- c(y0 = 0.05, mumax = 25,   K = 1200)

gr_test <- fit_growthmodel(FUN = grow_logistic, p = p, df3_34$days, df3_34$RFU,
                        lower = lower, upper = upper)

coef(gr_test)
plot(gr_test)
#mumax 8.14, fitting for all data I believe (e.g not factoring in replication)

#Let's try fitting nonparametric smoothing splines to the data

spl_test <- fit_spline(df3_34$RFU, df3_34$days)

par(mfrow = c(1, 2))
plot(spl_test, log = "y")
plot(spl_test)

coef(spl_test)
# this gives highly different results (mumax of 0.020261) — not sure why, need to think about this.

# For now, let's try to fit logarithmic growth curves to all of the temperatures for the 3rd population

p     <- c(y0 = 0.01, mumax = 5, K = 600)
lower <- c(y0 = 1e-6, mumax = 0,   K = 0)
upper <- c(y0 = 0.05, mumax = 25,   K = 1200)

## fit growth models to all data for Pop 3

grt_pop3 <- all_growthmodels(
  RFU ~ grow_logistic(days, parms) | temperature,
  data = mat[[3]], p = p, lower = lower, ncores = 2)

coef(grt_pop3)
par(mfrow = c(2, 3))
plot(grt_pop3)

# This is obviously not fitting the data for 10C or 40 C. I suspect the latter is caused by it being a negative slope.

df3_10<-subset(mat[[3]], mat[[3]]$temperature==10)

gr_test10 <- fit_growthmodel(FUN = grow_logistic, p = p, df3_10$days, df3_10$RFU,
                           lower = lower, upper = upper)

coef(gr_test10)
plot(gr_test10)

# That actually fits it fine...

str(mat[[3]])

#OK I think the issue is temperature is coded here as an integer, but we need it as a factor.
# Let's create another column

for (i in 1:15){
  mat[[i]]$temp<-as.factor(mat[[i]]$temperature)
}

grt_pop3.1 <- all_growthmodels(
  RFU ~ grow_logistic(days, parms) | temp, # temperature is a factor now
  data = mat[[3]], p = p, lower = lower, ncores = 2)

coef(grt_pop3.1)
par(mfrow = c(2, 3))
plot(grt_pop3.1)

#Still not working! I'm just going to do it manually, but first let's see what's happening with the 40 C group

df3_40<-subset(mat[[3]], mat[[3]]$temperature==40)

gr_test40 <- fit_growthmodel(FUN = grow_logistic, p = p, df3_40$days, df3_40$RFU,
                             lower = lower, upper = upper)

coef(gr_test40)
plot(gr_test40)

# I think the issue is our bounds on mumax, upper and lower.

p2     <- c(y0 = 0.01, mumax = -5, K = 50)
lower2 <- c(y0 = 1e-6, mumax = -25,   K = 0)
upper2 <- c(y0 = 0.05, mumax = 25,   K = 100)

gr_test40.1 <- fit_growthmodel(FUN = grow_logistic, p = p, df3_40$days, df3_40$RFU,
                             lower = lower2, upper = upper)

coef(gr_test40.1)
plot(gr_test40.1)

#OK, for loop incoming. We'll print the results and plot them all together

par(mfrow = c(2, 3))

pop3_mu<- matrix(vector(), 0, 4, dimnames=list(c(), c('temp', 'y0', 'mumax', 'K')))


for (t in c(10,16,22,28,34,40)){
  df<-subset(mat[[3]], mat[[3]]$temperature==t)
  log_gr <- fit_growthmodel(FUN = grow_logistic, p = p, df$days, df$RFU,
                               lower = lower, upper = upper)
  
  coef(log_gr)
  plot(log_gr)
  
  coef<-as.vector(c(t, coef(log_gr)[1:3]))
  pop3_mu<-rbind(pop3_mu, coef)
}

pop3_mu

# Aside from t = 40, for which the results are obviously bonkers, I think this looks OK.

# Let's practice fitting a TPC to this data (setting my for t40 to 0)?
pop3_mu[6,3]<-0

#OK, let's plot this relationship
par(mfrow = c(1, 1))
plot(mumax~temp, data=pop3_mu)

#list of models that can be implemented in bayes_TPC
get_models()

# Let's try a couple - ratkowsky, weibull, gaussian, and briere

get_formula("ratkowsky")
get_default_priors("ratkowsky")

# bayes_TPC can only work with Trait and Temp headers.

df.3.gr<-list(Trait = pop3_mu[,3], Temp=pop3_mu[,1])

Pop3_gr_rat <- b_TPC(data = df.3.gr, ## data
                    model = 'ratkowsky', ## model to fit
                    niter = 11000, ## total iterations
                    burn = 1000, ## number of burn in samples
                    samplerType = 'AF_slice', ## slice sampler
) 

summary(Pop3_gr_rat)

s1<-as.data.frame(Pop3_gr_rat$samples)
par(mfrow=c(2,2))
for(i in 1:4) acf(s1[,i], lag.max=50, main="", ylab = paste("ACF: ", names(s1)[i], sep=""))
plot(s1$T_max)

ppo_plot(Pop3_gr_rat, burn = 1000)

# I don't know what all this model structure is about (for example the slice sampler)
# Let's simplify!

Pop3_gr_rat.1 <- b_TPC(data = df.3.gr, ## data
                     model = 'ratkowsky', ## model to fit
                     niter = 11000, ## total iterations
                     burn = 1000, ## number of burn in samples
                     samplerType = 'RW', ## slice sampler,
                     thin = 100 #This should thin the model by recording data every 10 observations
) 

summary(Pop3_gr_rat.1)
s1.1<-as.data.frame(Pop3_gr_rat$samples)

Pop3_gr_rat$mcmc$getWAIC()
Pop3_gr_rat.1$mcmc$getWAIC()

head(s1)

plot(Pop3_gr_rat)


get_formula("gaussian")
get_default_priors("gaussian")

Pop3_gr_gau <- b_TPC(data = df.3.gr, ## data
                     model = 'gaussian', ## model to fit
                     niter = 11000, ## total iterations
                     burn = 1000, ## number of burn in samples
) 

summary(Pop3_gr_gau)
par(mfrow=c(1,1))

s2<-as.data.frame(Pop3_gr_gau$samples)
plot(s2$T_opt)
plot(s2$rmax)
# This is not working. The model is not fitting properly

# I'm not liking this bayes_TPC. It seems glitchy or maybe not fully fleshed out.
# Thinning isn't working, and the documentation isn't 100% of the way there. 

# I feel like it might be better to just fit this in nls and brms as classic non-linear models

# Do I need a package? Or can I do it manually?

#nls?

# OK so I am going to organize this as a testing/learning oppurtunity for me. If I understand nls, it basically
# enables fitting models (linear or not) to data using least squares regression. The twist is the explicit modelling
# of mathematical formulae, with terms that are then fit using the data. I'm going to try this gradually.

#1. Fit a linear regression to the Pop 3 data

plot(Trait~Temp, data=df.3.gr)

lm3 <- nls(Trait ~ m*Temp + b, #slope term is m, b is intercept
           start = list(m = 0.2, b = 0.5),
           data = df.3.gr) 

summary(lm3)
abline(a = 2.26774, b = 0.06799)

#OK, cool! That seems like the basic idea.
# Now the challenge is obviously expanding the complexity of that - account for inflection points (e.g. when T < x, do this...)
# and non-linear fits. 
# I don't know how to do that, but will figure it out.



elisa_mod2 <-
  nls(density ~ phi_1 / (1 + exp((phi_2 - log(pconc)) / phi_3)),
      data = elisa_assay,
      start = list(phi_1 = 2, phi_2 = 1, phi_3 = 1))


###########################################################################################


#3. brms (bayesian)

#Let's consider which model of TPC to fit
get_models()

