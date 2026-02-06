# Jason R Laurich

# February 3rd, 2026

# This file will fit µ to our time series growth data (RFUs) for C. reinhardtii populations, then fit various TPCs to the µ data.
# We fit 12 different TPC curves to our populations using the rTPC package, before comparing the best-performing models in R2jags
# We then iteratively fit the most suitable model (Lactin II) to all of our replicates. 

# Inputs: in processed-data : 01_temp_rfus.time.csv
# Outputs: in processed-data : 02_µ_estimates_temp.csv, 03_nls_rTPC_comparison_AICcs.csv, 04_TPC_summary.csv, 05_TPC_fits.csv
  # in R2jags-models : rep_i_lactin2.RData (Lactin II TPCs R2jags objects, saved but not pushed)

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(rTPC)
library(nls.multstart)
library(Deriv)
library(MuMIn)

library(R2jags)
library(mcmcplots)
library(bayestestR)

lactin2 <- function(temp, a, b, tmax, d.t) { 
  exp(a * temp) - exp(a * tmax - (tmax - temp) / d.t) + b
} # Define the Lactin II function

lactin2_deriv <- function(temp, a, b, tmax, d.t) {
  rho <- a
  T_max <- tmax
  delta_T <- d.t
  
  term1 <- rho * exp(rho * temp)
  term2 <- (1 / delta_T) * exp(rho * T_max - (T_max - temp) / delta_T)
  
  return(term1 - term2)
} # Derivative of the Lactin II function

lactin2_halfmax <- function(temp, a, b, tmax, d.t, r_half) {
  exp(a * temp) - exp(a * tmax - (tmax - temp) / d.t) + b - r_half
} # Calculate thermal breadth at 1/2 µ_max. 

# Upload & examine the data -----------------------------------------------

df <- read_csv("processed-data/01_temp_rfus_time.csv")
head(df) #RFU is density, days is time, temperature is temperature

df <- df %>% 
  select(well_plate, RFU, population, plate_type, temperature, days) %>% 
  rename(temp = temperature)

unique(df$population)
unique(df$temp)

df <- df %>% 
  filter(temp != 20,                      # 20C only includes single measurement, which we exclude.
         population != "cc1629",          # cc1629 is not relevant to the larger data set and was not replicated across other gradients, so we will remove that here.
         plate_type == "repeat") %>%      # Only keep data from plates that were read multiple times
  
  mutate(logRFU = log(RFU + 0.001)) %>%   # take the log of RFUs for comparing log-linear slopes to identify the exponential growth phase

  rename(well.ID = well_plate)

N0.df <- df %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df <- df %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") 

# Estimate µ --------------------------------------------------------------

df.µ <- data.frame(                          # Initializing a data frame to store the results for each well, pop, and temperature
  
  population = character(),                  # population ID
  temp = numeric(),                          # temperature
  well.ID = character(),                     # well ID
  
  µ = numeric()                              # intrinsic growth rate estimate
)

n <- 0 # completion tracker

for (i in unique(df$well.ID[df$well.ID >= 1])) { # This allows code below to be run in chunks if necessary, just re-start at the next i. 
  
  n <- n + 1
  
  df.i <- df %>% 
    filter(well.ID == i)  # Filter data down to a single well time series
  
  df.i <- df.i[order(df.i$days), ]  # in case data are not ordered by time of observation
  
  t.series <- unique(df.i$days)           # Re-initialize this internally - we will only save summary data for each unique pop x T x well combo
  
  ln.slopes <- c() # Re-initialize this too!
  
  for (z in t.series[2:length(t.series)]){   # Can't consider the slope just including time 0 
    
    df.i.sl <- df.i[df.i$days <= z, ]        # Subset the data to exclude time points above our window
    
    ln_slope <- lm(logRFU~days, data = df.i.sl) # calculate the log-linear slope
    
    ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
    
  }
  
  s <- max(2, which.max(ln.slopes))  # We need at least 3 data points to fit growth curves
  
  df.i.th <- df.i[df.i$days <= t.series[s+1], ] # Get the thresholded data according to our sliding window approach
  
  if(n_distinct(df.i.th$RFU) == 1) {
    
    µ.est <- 0    # If all of the thresholded values have the same RFU scores
  }else{
    
  
  µ.mod <- tryCatch(                          # Run the exponential growth model on the thresholded data
    nls_multstart(
      RFU ~ N0 * exp(r * days),
      data = df.i.th,
      start_lower = c(r = -4.5),
      start_upper = c(r = 4.5),
      iter = 500,
      supp_errors = "Y",
      control = nls.control(maxiter = 200)
    ),
    error = function(e) NULL
  )
  
  if (is.null(µ.mod)) {                   # if it fails, return NA
    µ.est <- NA_real_
  } else {
    µ.est <- coef(µ.mod)[["r"]]           # else return the model value for r
  } 
  
  }
  
  df.µ <- rbind(df.µ, data.frame(                     # add the data to the summary data frame
    
    population = df.i$population[1],                  # population ID
    temp = df.i$temp[1],                              # temperature
    well.ID = df.i$well.ID[1],                        # well ID
    
    µ = µ.est                                         # intrinsic growth rate. 
  ))
  
  print(paste("Done", n, "of ", length(unique(df$well.ID))))
  
}

df.µ %>% 
  filter(is.na(µ))  # No null values

write.csv(df.µ, "processed-data/02_µ_estimates_temp.csv",
          row.names = FALSE) # 888 measurements

# Explore model fitting in rTPC -------------------------------------------

# We're going to fit a whole bunch of models using the rTPC package. We selected a wide range of differnt types of models that broadly perform well,
# according to a recent analysis (Kontopoulos et al 2024)

# For some of these, we have adjusted existing rTPC models or written our own

bounded_ratkowsky <- function(temp, tmin, tmax, a, b) {                         # Modified Ratkowsky that sets growth above tmax to 0
  
  rate <- (a * (temp - tmin))^2 * (1 - exp(b * (temp - tmax)))^2                # Original Ratkowsky equation
  
  rate[temp > tmax] <- 0                                                        # Bound the rate to 0 for temp > tmax
  
  rate[rate < 0] <- 0                                                           # Ensure no negative rates
  
  return(rate)
}

ashrafi_simple <- function(temp, a, b, c) {                                     # Ashrafi II, not present in rTPC
  
  a + b * temp^(3/2) + c * temp^2                                               # Based on Kontopoulos 2024 supplement (3 parameters) : B(T) = a + b x T^3/2 + c x T^2
  
}

atkin <- function(temp, B0, a, b) {                                             # Atkin, not present in rTPC
  
  B0 * (a - b * temp)^(temp / 10)                                               # Based on Kontopoulos 2024 supplement (3 parameters) : B(T) = B0 x (a - b x T)^T/10
  
}

mitchell <- function(temp, a, b, Tpk) {                                         # Mitchell-Angilletta, not present in rTPC
  
  (a / (2 * b)) * (1 + cos(((temp - Tpk) / b) * pi))                            # Based on Kontopoulos 2024 supplement (3 parameters) : B(T) = (a/(2 x b)) x [1 + cos (((T-T_pk)/b) x pi)]
  
}

anal_kont <- function(temp, a, tmin, tmax) {                                    # Analytis-Kontomodimas, not present in rTPC
  
  a * (temp - tmin)^2 * (tmax - temp)                                           # Based on Kontopoulos 2024 supplement
  
}

eubank <- function(temp, a, tpk, b) {                                           # Eubank, not present in rTPC
  
  a / ((temp - tpk)^2 + b)
  
}

tay_sex <- function(temp, bpk, tmin, tpk) {                                     # Taylor-Sexton, not present in rTPC
  
  bpk * (-(temp - tmin)^4 + 2 * (temp - tmin)^2 * (tpk - tmin)^2) / (tpk - tmin)^4
  
}


aicc_df <- data.frame( # list of models to consider
  model = c("Deutsch", "Lactin2", "Ratkowsky bounded", 
            "Rezende", "Ashrafi II", "Atkin", "Mitchell-Angilletta", 
            "Analytis-Kontodimas", "Eubank", "Taylor-Sexton", "Briere", "Thomas1"),
  stringsAsFactors = FALSE
)

for (i in unique(df$population)){ # We are going to run these comparisons at the population level to limit the number of models we have to fit.
  
  df.i <- df.µ %>% 
    filter(population == i)
  
  aicc_pop <- c() # Temporary storage for this population
  
  # modified Deutsch
  
  start.vals <- get_start_vals(df.i$temp, df.i$µ, model_name = 'deutsch_2008')
  
  mod <- nls_multstart(µ ~ deutsch_2008(temp = temp, rmax, topt, ctmax, a),
                       data = df.i,
                       iter = c(4, 4, 4, 4), 
                       start_lower = start.vals - 10,
                       start_upper = start.vals + 10,
                       lower = get_lower_lims(df.i$temp, df.i$µ, model_name = 'deutsch_2008'),
                       upper = get_upper_lims(df.i$temp, df.i$µ, model_name = 'deutsch_2008'),
                       supp_errors = 'Y',
                       convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Lactin 2
  
  start.vals <- get_start_vals(df.i$temp, df.i$µ, model_name = 'lactin2_1995')
  
  mod <- nls_multstart(µ ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                       data = df.i,
                       iter = c(4, 4, 4, 4), 
                       start_lower = start.vals - 10,
                       start_upper = start.vals + 10,
                       lower = get_lower_lims(df.i$temp, df.i$µ, model_name = 'lactin2_1995'),
                       upper = get_upper_lims(df.i$temp, df.i$µ, model_name = 'lactin2_1995'),
                       supp_errors = 'Y',
                       convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Bounded Ratkowsky
  
  start.vals <- get_start_vals(df.i$temp, df.i$µ, model_name = 'ratkowsky_1983')
  
  mod <- nls_multstart(µ ~ bounded_ratkowsky(temp = temp, tmin, tmax, a, b),
                       data = df.i,
                       iter = c(4, 4, 4, 4), 
                       start_lower = start.vals - 10,
                       start_upper = start.vals + 10,
                       lower = get_lower_lims(df.i$temp, df.i$µ, model_name = 'ratkowsky_1983'),
                       upper = get_upper_lims(df.i$temp, df.i$µ, model_name = 'ratkowsky_1983'),
                       supp_errors = 'Y',
                       convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  #Rezende
  
  start.vals <- get_start_vals(df.i$temp, df.i$µ, model_name = 'rezende_2019')
  
  mod <- nls_multstart(µ ~ rezende_2019(temp = temp, q10, a, b, c),
                       data = df.i,
                       iter = c(4, 4, 4, 4), 
                       start_lower = start.vals - 10,
                       start_upper = start.vals + 10,
                       lower = get_lower_lims(df.i$temp, df.i$µ, model_name = 'rezende_2019'),
                       upper = get_upper_lims(df.i$temp, df.i$µ, model_name = 'rezende_2019'),
                       supp_errors = 'Y',
                       convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Ashrafi II
  
  mod <- nls_multstart(
    µ ~ ashrafi_simple(temp, a, b, c),
    data = df.i,
    iter = c(4, 4, 4),  # Number of iterations for starting values
    start_lower = c(a = -10, b = -10, c = -1),  # Lower bounds for start values - have to set ourselves (not present in rTPC)
    start_upper = c(a = 10, b = 10, c = 1),    # Upper bounds for start values
    lower = c(a = -30, b = -30, c = -15),   # Hard lower parameter bounds
    upper = c(a = 30, b = 30, c = 15),      # Hard upper parameter bounds
    supp_errors = 'Y',
    convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Atkin
  
  mod <- nls_multstart(
    µ ~ atkin(temp, B0, a, b),
    data = df.i,
    iter = c(4, 4, 4),
    start_lower = c(B0 = 0.1, a = 10, b = 0.01),
    start_upper = c(B0 = 5, a = 30, b = 0.1),
    lower = c(B0 = 0, a = 0, b = 0),
    upper = c(B0 = 10, a = 50, b = 1),
    supp_errors = 'Y',
    convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Mitchell-Angilletta
  
  mod <- nls_multstart(
    µ ~ mitchell(temp, a, b, Tpk),
    data = df.i,
    iter = c(4, 4, 4),
    start_lower = c(a = 0.1, b = 1, Tpk = 20),
    start_upper = c(a = 10, b = 10, Tpk = 40),
    lower = c(a = 0, b = 0.1, Tpk = 0),
    upper = c(a = Inf, b = Inf, Tpk = 50),
    supp_errors = 'Y',
    convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Analytis-Kontomodimas
  
  mod <- nls_multstart(
    µ ~ anal_kont(temp, a, tmin, tmax),
    data = df.i,
    iter = c(4, 4, 4),
    start_lower = c(a = 0.1, tmin = 0, tmax = 30),
    start_upper = c(a = 10, tmin = 10, tmax = 50),
    lower = c(a = 0, tmin = -Inf, tmax = -Inf),
    upper = c(a = Inf, tmin = Inf, tmax = Inf),
    supp_errors = 'Y',
    convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Eubank
  
  mod <- nls_multstart(
    µ ~ eubank(temp, a, tpk, b),
    data = df.i,
    iter = c(4, 4, 4),
    start_lower = c(a = 0.1, tpk = 20, b = 0.1),
    start_upper = c(a = 10, tpk = 35, b = 5),
    lower = c(a = 0, tpk = 0, b = 0.01),
    upper = c(a = Inf, tpk = Inf, b = Inf),
    supp_errors = 'Y',
    convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Taylor-Sexton
  
  mod <- nls_multstart(
    µ ~ tay_sex(temp, bpk, tmin, tpk),
    data = df.i,
    iter = c(4, 4, 4),
    start_lower = c(bpk = 0.1, tmin = 0, tpk = 20),
    start_upper = c(bpk = 10, tmin = 10, tpk = 40),
    lower = c(bpk = 0, tmin = -Inf, tpk = -Inf),
    upper = c(bpk = Inf, tmin = Inf, tpk = Inf),
    supp_errors = 'Y',
    convergence_count = FALSE
  )
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Briere II
  
  start.vals <- get_start_vals(df.i$temp, df.i$µ, model_name = 'briere2_1999')
  
  mod <- nls_multstart(µ~briere2_1999(temp = temp, Tmin, Tmax, a, b),
                       data = df.i,
                       iter = c(4,4,4,4),
                       start_lower = start.vals - 10,
                       start_upper = start.vals + 10,
                       lower = get_lower_lims(df.i$temp, df.i$µ, model_name = 'briere2_1999'),
                       upper = get_upper_lims(df.i$temp, df.i$µ, model_name = 'briere2_1999'),
                       supp_errors = 'Y',
                       convergence_count = F)
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  # Thomas 1
  
  start.vals <- get_start_vals(df.i$temp, df.i$µ, model_name = 'thomas_2012')
  
  mod <- nls_multstart(µ~thomas_2012(temp = temp, a, b, c, topt),
                       data = df.i,
                       iter = c(4,4,4,4),
                       start_lower = start.vals - 10,
                       start_upper = start.vals + 10,
                       lower = get_lower_lims(df.i$temp, df.i$µ, model_name = 'thomas_2012'),
                       upper = get_upper_lims(df.i$temp, df.i$µ, model_name = 'thomas_2012'),
                       supp_errors = 'Y',
                       convergence_count = F)
  
  aicc_pop <- c(aicc_pop, AICc(mod))
  
  aicc_df[[paste0("AICc_pop", i)]] <- aicc_pop #Add this population's AICc values as a new column
  
}

mean_aicc <- rowMeans(aicc_df[, -1], na.rm = TRUE)

aicc_df <- cbind(aicc_df[, 1, drop = FALSE], 
                 Mean_AICc = mean_aicc, 
                 aicc_df[, -1])

write.csv(aicc_df, "processed-data/03_nls_rTPC_comparison_AICcs.csv") # Save model comparison data

# Bayesian modelling ------------------------------------------------------

# Now we are going to fit Lactin II TPCs to each replicate using R2jags. 
# We selected this model because it (1) was consistently among the best-performing models based on AICc scores,
# (2) allows growth rates to be negative at low and high temperatures, and (3) is simple enough to permit fitting with R2jags. 

df.µ <- df.µ %>%
  mutate(rep.id = str_c(population, str_sub(well.ID, 1, 3), sep = ".")) # Need to create a unique replicate ID

length(unique(df.µ$rep.id)) # 148, should be 148! 888/148 = 6 temperatures per replicate. 37 x 4 = 148 (pops, reps)

summary.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  
  population = character(), # population ID
  rep.ID = character(),     # rep ID
  
  T.min = numeric(),        # Minimum T (calculus)
  T.min.min = numeric(),    # Minimum T (lower HDPI)
  T.min.max = numeric(),    # Minimum T (upper HDPI)
  T.min.na = numeric(),     # %% NA returns
  
  T.max = numeric(),        # Maximum T (calculus)
  T.max.min = numeric(),    # Maximum T (lower HDPI)
  T.max.max = numeric(),    # Maximum T (upper HDPI)
  T.max.na = numeric(),     # %% NA returns
  
  T.opt = numeric(),        # Optimal T (calculus)
  T.opt.min = numeric(),    # Optimal T (lower HDPI)
  T.opt.max = numeric(),    # Optimal T (upper HDPI)
  
  r.max = numeric(),        # Maximum growth rate (calculus)
  r.max.min = numeric(),    # Maximum growth rate (lower HDPI)
  r.max.max = numeric(),    # Maximum growth rate (upper HDPI)
  
  T.br.min = numeric(),     # T breadth (calculus)
  T.br.min.min = numeric(), # T breadth (lower HDPI)
  T.br.min.max = numeric(), # T breadth (upper HDPI)
  T.br.min.na = numeric(),  # T breadth na's
  
  T.br.max = numeric(),     # T breadth (calculus)
  T.br.max.min = numeric(), # T breadth (lower HDPI)
  T.br.max.max = numeric(), # T breadth (upper HDPI)
  T.br.max.na = numeric(),  # T breadth na's
  
  a = numeric(),            # parameter: a
  b = numeric(),            # parameter: b
  tmax = numeric(),         # parameter: tmax
  d.t = numeric(),          # parameter: deltaT
  
  a.mod = numeric(),        # nls.LM parameter: a
  b.mod = numeric(),        # nls.LM parameter: b
  tmax.mod = numeric(),     # nls.LM parameter: tmax
  d.t.mod = numeric(),      # nls.LM parameter: deltaT
  
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  
  population = character(), # population ID
  rep.ID = character(),     # rep ID
  
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~6000)
  stringsAsFactors = FALSE            
)

# Set generous MCMC settings still for our models. 
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model. 

parameters.lactin2 <- c("a", "b", "tmax", "d.t", "sigma", "r.pred") # parameters to estimate

Temp.xs <- seq(0, 45, 0.05) # Temperature gradient we're interested in - upped the granularity here
N.Temp.xs <-length(Temp.xs)

inits.lactin <- function() { # The other static initial values. 
  list(
    a = runif(1, 0.05, 0.15),  # More constrained initial values
    tmax = runif(1, 37, 43),
    d.t = runif(1, 1, 5),
    b = runif(1, -2.5, -1),
    sigma = runif(1, 0.1, 2)
  )
}

rep.ids <- unique(df.µ$rep.id) # Save the unique rep ids so we can run the foor loop in chunks if needed

n <- 0 # for tracking progress

for (i in rep.ids[129:length(rep.ids)]) { # For all replicates, here starting at 129 to account for multiple coding sessions to generate all of the objects
  
  n <- n + 1
  
  df.i <- df.µ %>% 
    filter(rep.id == i)
  
  df.i <- droplevels(df.i)
  
  trait <- df.i$µ    # format the data for jags
  N.obs <- length(trait)
  temp <- df.i$temp
  
  jag.data <- list(trait = trait, N.obs = N.obs, temp = temp, Temp.xs = Temp.xs, N.Temp.xs = N.Temp.xs) # assemble the jag data
  
  lac.jag <- jags(  # run the model
    data = jag.data, 
    inits = inits.lactin,
    parameters.to.save = parameters.lactin2, 
    model.file = "lactin.txt",
    n.thin = nt.fit, 
    n.chains = nc.fit, 
    n.burnin = nb.fit, 
    n.iter = ni.fit, 
    DIC = TRUE, 
    working.directory = getwd()
  )
  
  save(lac.jag, file = paste0("R2jags-models/rep_", i, "_lactin2.RData")) # save the lactin2 model
  # This folder is listed in gitignore, because the objects are too big to load
  
  post <- as.data.frame(lac.jag$BUGSoutput$sims.matrix) # The posterior distributions
  
  post <- post %>% 
    select(a, b, d.t, tmax)
  
  post.a <- median(post$a, na.rm = T)    # Extract parameters
  post.b <- median(post$b, na.rm = T)
  post.tmax <- median(post$tmax, na.rm = T)
  post.d.t <- median(post$d.t, na.rm = T)
  
  # Calculate summary metrics for the curve across all posteriors
  
  # Topt
  calc_Topt <- function(a, b, tmax, d.t) {
    tryCatch(
      uniroot(
        function(temp) lactin2_deriv(temp, a, b, tmax, d.t),
        interval = c(-10, 45)
      )$root,
      error = function(e) NA
    )
  }
  
  T.opt <- mapply(
    calc_Topt,
    a = post$a,
    b = post$b,
    tmax = post$tmax,
    d.t = post$d.t
  )
  
  #rmax
  r.max <- mapply(
    lactin2,
    temp = T.opt,
    a = post$a,
    b = post$b,
    tmax = post$tmax,
    d.t = post$d.t
  )
  
  #Tmin
  Tmin.safe <- function(Topt, a, b, tmax, d.t, # Tmin equation
                        lower = -50) {
    
    f_low  <- lactin2(temp = lower, a = a, b = b, tmax = tmax, d.t = d.t)
    f_high <- lactin2(temp = Topt,  a = a, b = b, tmax = tmax, d.t = d.t)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2(temp, a, b, tmax, d.t),
      interval = c(lower, Topt)
    )$root
  }
  
  T.min <- mapply(
    Tmin.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t
  )
  
  #Tmax
  Tmax.safe <- function(Topt, a, b, tmax, d.t, # Tmax equation
                        upper = 50) {
    
    f_low  <- lactin2(temp = Topt, a = a, b = b, tmax = tmax, d.t = d.t)
    f_high <- lactin2(temp = upper,  a = a, b = b, tmax = tmax, d.t = d.t)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2(temp, a, b, tmax, d.t),
      interval = c(Topt, upper)
    )$root
  }
  
  T.max <- mapply(
    Tmax.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t
  )
  
  # Tbr
  r.half <- r.max/2
  
  Tmin.half.safe <- function(Topt, a, b, tmax, d.t, r_half, # Tmin equation
                             lower = -10) {
    
    f_low  <- lactin2_halfmax(temp = lower, a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    f_high <- lactin2_halfmax(temp = Topt,  a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2_halfmax(temp, a, b, tmax, d.t, r_half),
      interval = c(lower, Topt)
    )$root
  }
  
  T.min.half <- mapply(
    Tmin.half.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t,
    r_half = r.half
  )
  
  Tmax.half.safe <- function(Topt, a, b, tmax, d.t, r_half, # Tmax equation
                             upper = 50) {
    
    f_low  <- lactin2_halfmax(temp = Topt, a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    f_high <- lactin2_halfmax(temp = upper,  a = a, b = b, tmax = tmax, d.t = d.t, r_half)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2_halfmax(temp, a, b, tmax, d.t, r_half),
      interval = c(Topt, upper)
    )$root
  }
  
  T.max.half <- mapply(
    Tmax.half.safe,
    Topt = T.opt,
    a    = post$a,
    b    = post$b,
    tmax = post$tmax,
    d.t  = post$d.t,
    r_half = r.half
  )

  summary.df <- rbind(summary.df, data.frame(                                   # Add summary data
    
    population = df.i$population[1],                                            # population ID
    rep.ID = df.i$rep.id[1],                                                    # rep ID
    
    T.min = median(T.min, na.rm = T),                                           # Minimum T (calculus)
    T.min.min = hdi(T.min[!is.na(T.min)], credMass = 0.95)$CI_low,              # Minimum T (lower HDPI)
    T.min.max = hdi(T.min[!is.na(T.min)], credMass = 0.95)$CI_high,             # Minimum T (upper HDPI)
    T.min.na = mean(is.na(T.min)),                                              # % NA returns
    
    T.max = median(T.max, na.rm = T),                                           # Maximum T (calculus)
    T.max.min = hdi(T.max[!is.na(T.max)], credMass = 0.95)$CI_low,              # Maximum T (lower HDPI)
    T.max.max = hdi(T.max[!is.na(T.max)], credMass = 0.95)$CI_high,             # Maximum T (upper HDPI)
    T.max.na = mean(is.na(T.max)),                                              # % NA returns
    
    T.opt = median(T.opt, na.rm = T),                                           # Optimal T (calculus)
    T.opt.min = hdi(T.opt[!is.na(T.opt)], credMass = 0.95)$CI_low,              # Optimal T (lower HDPI)
    T.opt.max = hdi(T.opt[!is.na(T.opt)], credMass = 0.95)$CI_high,             # Optimal T (upper HDPI)
    
    r.max = median(r.max, na.rm = T),                                           # Maximum growth rate  (calculus)
    r.max.min = hdi(r.max[!is.na(r.max)], credMass = 0.95)$CI_low,              # Maximum growth rate  (lower HDPI)
    r.max.max = hdi(r.max[!is.na(r.max)], credMass = 0.95)$CI_high,             # Maximum growth rate  (upper HDPI)
    
    T.br.min = median(T.min.half, na.rm = T),                                   # T breadth (calculus)
    T.br.min.min = hdi(T.min.half[!is.na(T.min.half)], credMass = 0.95)$CI_low, # T breadth (lower HDPI)
    T.br.min.max = hdi(T.min.half[!is.na(T.min.half)], credMass = 0.95)$CI_high,# T breadth (upper HDPI)
    T.br.min.na = mean(is.na(T.min.half)),                                      # % NA returns
    
    T.br.max = median(T.max.half, na.rm = T),                                   # T breadth (calculus)
    T.br.max.min = hdi(T.max.half[!is.na(T.max.half)], credMass = 0.95)$CI_low, # T breadth (lower HDPI)
    T.br.max.max = hdi(T.max.half[!is.na(T.max.half)], credMass = 0.95)$CI_high,# T breadth (upper HDPI)
    T.br.max.na = mean(is.na(T.max.half)),                                      # % NA returns
    
    a = post.a,                                                                 # parameter: a
    b = post.b,                                                                 # parameter: b
    tmax = post.tmax,                                                           # parameter: tmax
    d.t = post.d.t,                                                             # parameter: deltaT
    
    a.mod = lac.jag$BUGSoutput$summary[1,1],                                    # Jags parameter: a
    b.mod = lac.jag$BUGSoutput$summary[2,1],                                    # Jags parameter: b
    tmax.mod = lac.jag$BUGSoutput$summary[907,1],                               # Jags parameter: tmax
    d.t.mod =  lac.jag$BUGSoutput$summary[3,1]                                  # Jags parameter: deltaT
  ))
  
  for (j in c(1:3,907)){
    fit.df <- rbind(fit.df, data.frame(                         # Model performance data
     
      population = df.i$population[1],                          # population ID
      rep.ID = df.i$rep.id[1],                                  # rep ID 
      
      Parameter = rownames(lac.jag$BUGSoutput$summary)[j],      # Model parameter (e.g. cf.a, cf.tmax, etc.)
      mean = lac.jag$BUGSoutput$summary[j,1],                   # Posterior mean
      Rhat = lac.jag$BUGSoutput$summary[j,8],                   # Rhat values
      n.eff = lac.jag$BUGSoutput$summary[j,9],                  # Sample size estimates (should be ~6000)
      stringsAsFactors = FALSE            
    ))
  }
  
  print(paste("Done", n, "of ", length(unique(df.µ$rep.id))))
  
}
 
write.csv(summary.df, "processed-data/04_TPC_summary.csv",
          row.names = FALSE) # 144 measurements

write.csv(fit.df, "processed-data/05_TPC_fits.csv",
          row.names = FALSE) # 144 measurements
