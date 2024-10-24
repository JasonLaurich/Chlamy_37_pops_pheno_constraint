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
library(rTPC)

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

#OK, because we know that the data were taken during the exponential growth phase, we should be able to just
# fit a slope to the logged data to estimate r, no?

df3_34$logRFU<-log(df3_34$RFU)
ggplot(aes(x=days, y=logRFU), data=df3_34) + geom_point()

loglin<-lm(logRFU~days, data=df3_34)
summary(loglin)
#gives me an r of 2.87134
#IS THIS PER CAPITA? No right?

grtnls<- nls(logRFU~ m*days + b,
             data=df3_34)

summary(grtnls) # Does the exact same thing! Not surprisingly.

#Let's loop this through my various temperatures here
mat[[3]]$logRFU<-log(mat[[3]]$RFU)

pop3_r<- matrix(vector(), 0, 2, dimnames=list(c(), c('temp', 'r')))

for (t in c(10,16,22,28,34,40)){
  df<-subset(mat[[3]], mat[[3]]$temperature==t)
  lm_gr <- lm(logRFU~days, data=df)
  
  coef<-as.vector(c(t, lm_gr$coefficients[2]))
  pop3_r<-rbind(pop3_r, coef)
}

#With nls should be the same
pop3_r2<- matrix(vector(), 0, 2, dimnames=list(c(), c('temp', 'r')))

for (t in c(10,16,22,28,34,40)){
  df<-subset(mat[[3]], mat[[3]]$temperature==t)
  lm_nls <- nls(logRFU~r*days + b, data=df)
  
  par<-as.data.frame(lm_nls$m$getPars())
  coef<-as.vector(c(t, par[1,1]))
  pop3_r2<-rbind(pop3_r2, coef)
}

# Yep, they are!

ggplot(aes(x=temp, y=r), data=pop3_r) + geom_point()

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
# I don't know how to do that, but I will figure it out!

#2. Let's try a simple model. The briere one seems like a good starting point!

bri3 <- nls (Trait ~ q*Temp*(Temp - Tmin)*sqrt((Tmax > Tmin)*abs(Tmax-Temp))*(Tmax>Temp)*(Temp>Tmin), 
             start = list(q = 0.5, Tmax = 32, Tmin = 5),
             data= df.3.gr)

bri3.1 <- nls (Trait ~ q*Temp*(Temp - Tmin)*sqrt((Tmax > Tmin)*abs(Tmax-Temp))*(Tmax>Temp)*(Temp>Tmin), 
               start = list(q = 0.2, Tmax = 35, Tmin = 10),
               data=df.3.gr)

# Neither is working. Try the nl_multstart option?

start_vals <- get_start_vals(df.3.gr$Temp, df.3.gr$Trait, model_name = 'briere2_1999')

bri3.2 <- nls_multstart(Trait~briere2_1999(temp = Temp, Tmin, Tmax, a, b),
                     data = df.3.gr,
                     iter = c(4,4,4,4),
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = get_lower_lims(df.3.gr$Temp, df.3.gr$Trait, model_name = 'briere2_1999'),
                     upper = get_upper_lims(df.3.gr$Temp, df.3.gr$Trait, model_name = 'briere2_1999'),
                     supp_errors = 'Y',
                     convergence_count = F)

summary(bri3.2)

#Let's get the predictions and store them for plotting!
preds <- data.frame(Temp = seq(min(df.3.gr$Temp), max(df.3.gr$Temp), length.out = 100))
preds <- broom::augment(bri3.2, newdata = preds)

df3gr<-as.data.frame(df.3.gr)

ggplot(preds) + geom_point(aes(Temp, Trait), df3gr) +
  geom_line(aes(Temp, .fitted), col = 'darkslateblue') + theme_bw()

# OK wow, that is a TPC! Look at me go

#So, could we have fitted this without using any of nls_multstart's under the hood mechanics?

#Back to the model specified by bayes_TPC? It's too different than what nls fit.
#Let's try replicating the nls_multstart model specification, and plug in what we know for the parameters...

bri4 <- nls (Trait ~ a*Temp*(Temp - Tmin)*(Tmax - Temp)^(1/b), 
             start = list(a = 4.007e-03, b = 2.500e+00, Tmax = 4.000e+01 , Tmin = 4.499e+00),
             data= df3gr)

# So it won't fit... why? 'Error in nlsModel(formula, mf, start, wts, scaleOffset = scOff, nDcentral = nDcntr) : 
# singular gradient matrix at initial parameter estimates''?

# But it will work in multstart, from which I pulled these parameter estimates...
# Let's try crossing nlsmultstart and this last nls model

# This should work then, if this equation adequately reflects the Briere function wrapper?

bri3.3 <- nls_multstart(Trait ~ a*Temp*(Temp - Tmin)*(Tmax - Temp)^(1/b),
                        data = df3gr,
                        iter = c(4,4,4,4),
                        start_lower = start_vals - 10,
                        start_upper = start_vals + 10,
                        lower = get_lower_lims(df3gr$Temp, df3gr$Trait, model_name = 'briere2_1999'), # try putting in computed values
                        upper = get_upper_lims(df3gr$Temp, df3gr$Trait, model_name = 'briere2_1999'),
                        supp_errors = 'Y',
                        convergence_count = F)




###########################################################################################

#3. brms (bayesian)

# Let's start by calculating growth rates using brms, like we did with nls

brm3.34 <- brm(logRFU~days,
               data=df3_34,
               family=gaussian())

summary(brm3.34) # Excellent, very close to the nls/lm data! Automatically fits an intercept
plot(brm3.34)

# Ok, let's fit this for all of the temperatures

pop3_br<- matrix(vector(), 0, 2, dimnames=list(c(), c('temp', 'r')))

for (t in c(10,16,22,28,34,40)){
  df<-subset(mat[[3]], mat[[3]]$temperature==t)
  brm_r <- brm(logRFU~days, data=df)
  
  fit<-as.data.frame(brm_r$fit)
  coef<-as.vector(c(t, mean(fit$b_days)))
  pop3_br<-rbind(pop3_br, coef)
}

# I feel like this takes forever, may as well just fit the growth curves with lm or nls
# Very similar results to nls, which is must faster.

#Once again, let's try fitting a simple linear model to the data, as we did with nls

brm_lin <- brm(r ~ temp, data = pop3_br, family = gaussian())

summary(brm_lin)
plot(brm_lin)
plot(conditional_effects(brm_lin), points=T)
#Ok obviousyl that is terrible, but... it worked!

priorbri <- prior(normal(5,10), nlpar="Tmin") +
  prior(normal(35,10), nlpar="Tmax") +
  prior(normal(0,1), nlpar ="a") +
  prior(normal(2.5, 2), nlpar= "b")

priorbri.1 <- prior(normal(10,2), nlpar="Tmin") +
  prior(normal(10,2), nlpar="Tmax") +
  prior(normal(10,2), nlpar ="a") +
  prior(normal(10,2), nlpar= "b")

brm_bri <- brm(bf(r ~ a*temp*(temp - Tmin)*(Tmax - temp)^(1/b), Tmax + Tmin + a + b ~ 1, nl = T),
               data = pop3_br, prior = priorbri.1)

summary(brm_bri)


prior1 <- prior(normal(1, 2), nlpar = "a") +
  prior(normal(0, 2), nlpar = "b")

fit1 <- brm(bf(r ~ a * exp(b * temp), a + b ~ 1, nl = TRUE),
            data = pop3_br, prior = prior1)

summary(fit1)

#Let's try in jags, using the code Joey worked on in the Anopheles-rate-summation project
# https://github.com/JoeyBernhardt/anopheles-rate-summation/blob/master/AnalysisDemo.R

##### MCMC Settings
# Number of posterior distribution elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 25000 # number of iterations in each chain
nb <- 5000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Derived Quantity Settings
Temp.xs <- seq(0, 45, 0.1) # temperature gradient to calculate derived quantities over
N.Temp.xs <-length(Temp.xs)

#####  inits Function
inits<-function(){list(
  cf.q = 0.01,
  cf.Tm = 35,
  cf.T0 = 5,
  cf.sigma = rlnorm(1))}

#####  Parameters to Estimate
parameters <- c("cf.q", "cf.T0", "cf.Tm","cf.sigma", "z.trait.mu.pred")

##### Organize Data for JAGS
pop3_br<-as.data.frame(pop3_br)

trait <- pop3_br$r
N.obs <- length(pop3_br$r)
temp <- pop3_br$temp

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs)

##### Run JAGS - **select correct model file** - Briere function, truncated normal distribution
model_bite_rate_constant <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere.txt",
                                 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=T, working.directory=getwd())
