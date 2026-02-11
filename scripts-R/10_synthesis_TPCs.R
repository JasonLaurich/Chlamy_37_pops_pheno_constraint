# Jason R Laurich

# February 11th, 2026

# In this script I will fit TPCs to the synthesis data set, working with estimates from a wide range of phytoplankton and studies
  # In some cases, I will work directly with published estimates of µ, while in others will have to first estimate µ from 
  # time-series data (as in scripts 01-05.R). 

# Throughout, we will trim data sets down to species that meet certain criteria. (1) They must have 5 or more growth estimates at
  # unique temperatures, and (2) they must have at least one data point above a putative Topt (the temperature at which growth
  # peaks in the raw data). 

# Unlike with Chlamydomonas data from the experimental evolution work, we will fit these models in nls to accommodate the 
  # diversity present in these synthesis data sets, which makes fitting models in R2jags difficult. 

# Inputs: 28_Thomas_2012_raw_data.csv, 31_Bestion_2018_raw_data.csv, 34_Lewington-Pearce_2019_raw_data.csv, 38_Edwards_2016_raw_data.csv,
  # 41_Levasseur2025_l_rawdata.csv, 42_Levasseur2025_n_rawdata.csv, 43_Levasseur2025_p_rawdata.csv
# Outputs: in processed-data : 29_Thomas_2012_TPCs.csv, 30_Thomas_2012_TPCs_fits.csv, 32_Bestion_2018_TPCs.csv, 
  # 33_Bestion_2018_TPCs_fits.csv, 35_Lewington-Pearce_2019_µ_estimates_temp.csv, 36_Lewington_2019_TPCs.csv, 
  # 37_Lewington_2019_TPCs_fits.csv, 39_Edwards_2016_TPCs.csv, 40_Edwards_2016_TPCs_fits.csv,
  # 44_Levasseur_2025_µ_estimates_temp.csv, 45_Levasseur_2025_TPCs.csv, 46_Levasseur_2025_TPCs_fits.csv

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(rTPC)
library(nls.multstart)
library(Deriv)
library(car)
library(boot)
library(minpack.lm)

lactin2 <- function(temp, a, tmax, d.t, b) { 
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
} # OK we're going to modify the function to calculate T_breadth, based on a modified lactin.

# Thomas 2012 ------------------------------------------------------------

###### Load the data ######

df.t.raw <- read.csv('processed-data/28_Thomas_2012_raw_data.csv') # Thomas raw data
head(df.t.raw)

length(unique(df.t.raw$id.number))

min(df.t.raw[df.t.raw$Growth.rate>0,]$Temperature) # Tmin => -1.8
max(df.t.raw[df.t.raw$Growth.rate>0,]$Temperature) # Tmax <= 37

df.t.raw %>% # Number of observations by species
  group_by(id.number) %>%
  summarise(n = n(), .groups = "drop") %>% 
  print(n = 200)

df.t.raw %>% # 24 species are lacking sufficient data. 
  group_by(id.number) %>%
  summarise(n = n(), .groups = "drop") %>% 
  filter(n <5) %>% 
  print(n = 200)

df.t <- df.t.raw %>% # Remove species with fewer than 5 data points
  group_by(id.number) %>%
  mutate(n = n()) %>% 
  filter(n >= 5) 

length(unique(df.t$id.number)) # Down to 170 out of 194 species. 

df.t <- df.t %>% # Remove species where there are no growth data above putative Topt. 
  group_by(id.number) %>%
  filter(!(max(Temperature) %in% Temperature[Growth.rate == max(Growth.rate)])) 

length(unique(df.t$id.number)) # Down to 169 out of 194 species.

df.t <-df.t %>% 
  rename(temp = Temperature,
         mu = Growth.rate)

###### Model fitting ######

thomas.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Study = character(),      # Study
  
  RSE = numeric(),          # Raw residual square error
  RSE.mu = numeric(),       # Residual square error standardized against max µ value
  RSE.range = numeric(),    # Residual square error standardized against range of µ values
  RSE.sd = numeric(),       # Residual square error standardized against sd of µ values
  
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
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Study = character(),      # Study
  Parameter = character(),  # Model parameter (e.g. a, tmax, etc.)
  est = numeric(),          # Estimate
  se = numeric(),           # Standard error
  p = numeric(),            # p-value
  stringsAsFactors = FALSE        
)

n <-0 # progression tracker

for (i in unique(df.t$id.number[df.t$id.number >= 1])) { # for each species ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.t %>% 
    filter(id.number == i)
  df.i <- droplevels(df.i)
  
  # df.i <- add_left_zero_if_needed(df.i)
  
  lac_nls <- nls_multstart(mu ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = c(0.01, -2.5, 5, 1),
                           start_upper = c(0.19, -0.5, max(df.i$temp), 11),
                           lower = c(0, -3, min(df.i$temp), 0.1),
                           upper = c(0.5, -0.001, max(df.i$temp) + 5, 40),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  sum <- summary(lac_nls)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  lac.a <- df.nls[1,1] # Extract parameters
  lac.b <- df.nls[2,1]
  lac.tmax <- df.nls[3,1]
  lac.d.t <- df.nls[4,1]
  
  lac.LM <- nlsLM(mu ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                  data = df.i,
                  start = c(a = lac.a, b = lac.b, tmax = lac.tmax, delta_t = lac.d.t),
                  lower = c(a = max(0.001,lac.a - 0.05), b = lac.b - 0.5, tmax = max(0, lac.tmax -3), delta_t = max(0.1, lac.d.t - 2)),
                  upper = c(a = lac.a + 0.05, b = min(lac.b + 0.5, -0.001), tmax = lac.tmax + 3, delta_t = lac.d.t + 2),
                  control = nls.lm.control(maxiter=500)
  )
  
  sum <- summary(lac.LM)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  boot <- Boot(lac.LM, method = 'residual')
  
  post <- as.data.frame(boot$t) # Get the bootstrapped values
  post$d.t <- post$delta_t
  
  boot.a <- median(post$a, na.rm = T)    # Extract parameters
  boot.b <- median(post$b, na.rm = T)
  boot.tmax <- median(post$tmax, na.rm = T)
  boot.d.t <- median(post$d.t, na.rm = T)
  
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
                        lower = -100) {
    
    f_low  <- lactin2(temp = lower, a = a, b = b, tmax = tmax, d.t = d.t)
    f_high <- lactin2(temp = Topt,  a = a, b = b, tmax = tmax, d.t = d.t)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2(temp, a, tmax, d.t, b),
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
      function(temp) lactin2(temp, a, tmax, d.t, b),
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
  
  thomas.summ.df <- rbind(thomas.summ.df, data.frame(                           # Add summary data
    Sp.id = df.i$id.number[1],                                                  # Species #
    Sp.name = df.i$Species.name[1],                                             # Species name
    Study = df.i$Study[1],                                                      # Study
    
    RSE = sigma(lac.LM),                                                        # Raw residual square error
    RSE.mu = sigma(lac.LM)/max(df.i$mu),                                        # Residual square error standardized against sd of µ values
    RSE.range = sigma(lac.LM)/(max(df.i$mu) - min(df.i$mu)),                    # Residual square error standardized against range of µ values
    RSE.sd = sigma(lac.LM)/sd(df.i$mu),                                         # Residual square error standardized against sd of µ values
    
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
    
    a = boot.a,                                                                 # parameter: a
    b = boot.b,                                                                 # parameter: b
    tmax = boot.tmax,                                                           # parameter: tmax
    d.t = boot.d.t,                                                             # parameter: deltaT
    
    a.mod = df.nls[1,1],                                                        # nls.LM parameter: a
    b.mod = df.nls[2,1],                                                        # nls.LM parameter: b
    tmax.mod = df.nls[3,1],                                                     # nls.LM parameter: tmax
    d.t.mod = df.nls[4,1]                                                       # nls.LM parameter: deltaT
  ))
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(                                         # Model performance data
      Sp.id = df.i$id.number[1],                                                # Species #
      Sp.name = df.i$Species.name[1],                                           # Species name
      Study = df.i$Study[1],                                                    # Study
      Parameter = df.nls$parameter[j],                                          # Model parameter (e.g. a, b, tmax etc.)
      est = df.nls$Estimate[j],                                                 # Estimate
      se = df.nls$`Std. Error`[j],                                              # Error
      p = df.nls$`Pr(>|t|)`[j],                                                 # p-values
      stringsAsFactors = FALSE            
    ))
  }
  
  print(paste("Done", n, "of ", length(unique(df.t$id.number))))
  
}

write.csv(thomas.summ.df, "processed-data/28_Thomas_2012_TPCs.csv") # Save Thomas 2012 summary table
write.csv(fit.df, "processed-data/29_Thomas_2012_TPCs_fits.csv") # Save model fit summary table

# Bestion 2018 ------------------------------------------------------------

###### Load the data ######

df.b.raw <- read.csv("processed-data/31_Bestion_2018_raw_data.csv") # Raw growth data
head(df.b.raw)

df.b.raw <- df.b.raw[df.b.raw$Phosphate_c == 30,] # Going with 30 because there is more data for more species here. 

df.b.raw %>% # Number of observations by species
  group_by(SpeciesNb) %>%
  summarise(n = n(), .groups = "drop") %>% 
  print(n = 200) # They all have more than 

df.b <- df.b.raw %>% # Remove species with fewer than 5 data points
  group_by(SpeciesNb) %>%
  mutate(n = n()) %>% 
  filter(n >= 5) 

df.b <-df.b %>% 
  rename(temp = Temperature_c,
         id.number = SpeciesNb,
         Species.name = SpeciesName)


df.b <- df.b %>% # Remove species where there are no growth data above putative Topt. 
  group_by(id.number) %>%
  filter(!(max(temp) %in% temp[mu == max(mu)])) 

length(unique(df.b$id.number)) # Down to 5 out of 6 species.

###### Model fitting ######

bestion.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  
  RSE = numeric(),          # Raw residual square error
  RSE.mu = numeric(),       # Residual square error standardized against max µ value
  RSE.range = numeric(),    # Residual square error standardized against range of µ values
  RSE.sd = numeric(),       # Residual square error standardized against sd of µ values
  
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
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Parameter = character(),  # Model parameter (e.g. a, tmax, etc.)
  est = numeric(),          # Estimate
  se = numeric(),           # Standard error
  p = numeric(),            # p-value
  stringsAsFactors = FALSE        
)

n <-0 # progression tracker

for (i in unique(df.b$id.number[df.b$id.number >= 1])) { # for each species ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.b %>% 
    filter(id.number == i)
  df.i <- droplevels(df.i)
  
  # df.i <- add_left_zero_if_needed(df.i)
  
  lac_nls <- nls_multstart(mu ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = c(0.01, -2.5, 5, 1),
                           start_upper = c(0.19, -0.5, max(df.i$temp), 11),
                           lower = c(0, -3, min(df.i$temp), 0.1),
                           upper = c(0.5, -0.001, max(df.i$temp) + 5, 40),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  sum <- summary(lac_nls)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  lac.a <- df.nls[1,1] # Extract parameters
  lac.b <- df.nls[2,1]
  lac.tmax <- df.nls[3,1]
  lac.d.t <- df.nls[4,1]
  
  lac.LM <- nlsLM(mu ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                  data = df.i,
                  start = c(a = lac.a, b = lac.b, tmax = lac.tmax, delta_t = lac.d.t),
                  lower = c(a = max(0.001,lac.a - 0.05), b = lac.b - 0.5, tmax = max(0, lac.tmax -3), delta_t = max(0.1, lac.d.t - 2)),
                  upper = c(a = lac.a + 0.05, b = min(lac.b + 0.5, -0.001), tmax = lac.tmax + 3, delta_t = lac.d.t + 2),
                  control = nls.lm.control(maxiter=500)
  )
  
  sum <- summary(lac.LM)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  boot <- Boot(lac.LM, method = 'residual')
  
  post <- as.data.frame(boot$t) # Get the bootstrapped values
  post$d.t <- post$delta_t
  
  boot.a <- median(post$a, na.rm = T)    # Extract parameters
  boot.b <- median(post$b, na.rm = T)
  boot.tmax <- median(post$tmax, na.rm = T)
  boot.d.t <- median(post$d.t, na.rm = T)
  
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
                        lower = -100) {
    
    f_low  <- lactin2(temp = lower, a = a, b = b, tmax = tmax, d.t = d.t)
    f_high <- lactin2(temp = Topt,  a = a, b = b, tmax = tmax, d.t = d.t)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2(temp, a, tmax, d.t, b),
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
      function(temp) lactin2(temp, a, tmax, d.t, b),
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
  
  bestion.summ.df <- rbind(bestion.summ.df, data.frame(                           # Add summary data
    Sp.id = df.i$id.number[1],                                                  # Species #
    Sp.name = df.i$Species.name[1],                                             # Species name
    
    RSE = sigma(lac.LM),                                                        # Raw residual square error
    RSE.mu = sigma(lac.LM)/max(df.i$mu),                                        # Residual square error standardized against sd of µ values
    RSE.range = sigma(lac.LM)/(max(df.i$mu) - min(df.i$mu)),                    # Residual square error standardized against range of µ values
    RSE.sd = sigma(lac.LM)/sd(df.i$mu),                                         # Residual square error standardized against sd of µ values
    
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
    
    a = boot.a,                                                                 # parameter: a
    b = boot.b,                                                                 # parameter: b
    tmax = boot.tmax,                                                           # parameter: tmax
    d.t = boot.d.t,                                                             # parameter: deltaT
    
    a.mod = df.nls[1,1],                                                        # nls.LM parameter: a
    b.mod = df.nls[2,1],                                                        # nls.LM parameter: b
    tmax.mod = df.nls[3,1],                                                     # nls.LM parameter: tmax
    d.t.mod = df.nls[4,1]                                                       # nls.LM parameter: deltaT
  ))
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(                                         # Model performance data
      Sp.id = df.i$id.number[1],                                                # Species #
      Sp.name = df.i$Species.name[1],                                           # Species name
      Parameter = df.nls$parameter[j],                                          # Model parameter (e.g. a, b, tmax etc.)
      est = df.nls$Estimate[j],                                                 # Estimate
      se = df.nls$`Std. Error`[j],                                              # Error
      p = df.nls$`Pr(>|t|)`[j],                                                 # p-values
      stringsAsFactors = FALSE            
    ))
  }
  
  print(paste("Done", n, "of ", length(unique(df.b$id.number))))
  
}

write.csv(bestion.summ.df, "processed-data/32_Bestion_2018_TPCs.csv") # Save Bestion_2018 summary table
write.csv(fit.df, "processed-data/33_Bestion_2018_TPCs_fits.csv") # Save model fit summary table

# Lewington-Pearce 2019 ------------------------------------------------------------

###### Load the data ######

df <- read.csv("processed-data/34_Lewington-Pearce_2019_raw_data.csv") # Raw data file
head(df)
str(df)

df$Sp.fac <- as.factor(df$genus_species)

df$id <- paste0(df$block, ".", df$replicate)  # A unique identifier for each replicate.

df$log.fluorescence <- log(df$fluorescence + 0.001)

head(df)

###### Calculate µ across temperatures ######

df.t <- df %>% 
  filter(light_level == 105.5, nitrate_level == 1000)

temp <- as.vector(unique(df.t$temperature))# for looping through nitrate levels
ord.temp<- sort(temp)

df.r.exp <- data.frame(                # Summary dataframe for r_exp
  Sp.id = character(),                 # Species
  temperature = numeric(),             # Temperature
  id = character(),                    # Unique ID
  r.exp = numeric()                    # Thresholded r.exp
)

for (i in unique(df.t$Sp.fac)){ # for every species
  
  for (t in ord.temp){ # at every temperature
    
    df.i <- df.t[df.t$Sp.fac == i, ]
    df.it <- df.i[df.i$temperature == t, ]
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in unique(df.it$id)){ # Doing this separately for each replicate
      
      df.it.wl <- subset(df.it, as.numeric(df.it$id) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$day), ]
      df.it.wl <- df.it.wl %>% 
        mutate(N0 = fluorescence[1])
      
      t.series <- unique(df.it.wl$day) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$day <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(log.fluorescence~day, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
      s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
      
      df.it.wl.th <- df.it.wl[df.it.wl$day <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
      # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
      
      r_exp <- nls_multstart(fluorescence ~ N0 * exp(r*day),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp <- rbind(df.r.exp, data.frame(
          Sp.id = df.it.wl.th$Sp.fac[1],          
          temperature = df.it.wl$temperature[1],        
          id = df.it.wl$id[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          Sp.id = df.it.wl.th$Sp.fac[1],             # Species
          temperature = df.it.wl$temperature[1],     # Temperature
          id = df.it.wl$id[1],                       # ID
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
      }
      
    }
    
  }
  
}

write.csv(df.r.exp, "processed-data/35_Lewington-Pearce_2019_µ_estimates_temp.csv") # let's save the file.

# Can start from here

df.l.raw <- read.csv("processed-data/35_Lewington-Pearce_2019_µ_estimates_temp.csv") # growth data
head(df.l.raw)

df.l.raw %>% # Number of observations by species
  group_by(Sp.id) %>%
  summarise(n = n(), .groups = "drop") %>% 
  print(n = 200) # They all have more than 

df.l <- df.l.raw %>% # Remove species with fewer than 5 data points
  group_by(Sp.id) %>%
  mutate(n = n()) %>% 
  filter(n >= 5) 

df.l <-df.l %>% 
  ungroup() %>% 
  rename(temp = temperature,
         mu = r.exp,
         Species.name = Sp.id) %>% 
  mutate(id.number = as.numeric(as.factor(Species.name)))


df.l <- df.l %>% # Remove species where there are no growth data above putative Topt. 
  group_by(id.number) %>%
  filter(!(max(temp) %in% temp[mu == max(mu)])) 

length(unique(df.l$id.number)) # Still at 6! 

###### Model fitting ######

lewington.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  
  RSE = numeric(),          # Raw residual square error
  RSE.mu = numeric(),       # Residual square error standardized against max µ value
  RSE.range = numeric(),    # Residual square error standardized against range of µ values
  RSE.sd = numeric(),       # Residual square error standardized against sd of µ values
  
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
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Parameter = character(),  # Model parameter (e.g. a, tmax, etc.)
  est = numeric(),          # Estimate
  se = numeric(),           # Standard error
  p = numeric(),            # p-value
  stringsAsFactors = FALSE        
)

n <-0 # progression tracker

for (i in unique(df.l$id.number[df.l$id.number >= 1])) { # for each species ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.l %>% 
    filter(id.number == i) %>% 
    filter(!is.na(mu))
  
  df.i <- droplevels(df.i)
  
  # df.i <- add_left_zero_if_needed(df.i)
  
  lac_nls <- nls_multstart(mu ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = c(0.01, -2.5, 5, 1),
                           start_upper = c(0.19, -0.5, max(df.i$temp), 11),
                           lower = c(0, -3, min(df.i$temp), 0.1),
                           upper = c(0.5, -0.001, max(df.i$temp) + 5, 40),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  sum <- summary(lac_nls)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  lac.a <- df.nls[1,1] # Extract parameters
  lac.b <- df.nls[2,1]
  lac.tmax <- df.nls[3,1]
  lac.d.t <- df.nls[4,1]
  
  lac.LM <- nlsLM(mu ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                  data = df.i,
                  start = c(a = lac.a, b = lac.b, tmax = lac.tmax, delta_t = lac.d.t),
                  lower = c(a = max(0.001,lac.a - 0.05), b = lac.b - 0.5, tmax = max(0, lac.tmax -3), delta_t = max(0.1, lac.d.t - 2)),
                  upper = c(a = lac.a + 0.05, b = min(lac.b + 0.5, -0.001), tmax = lac.tmax + 3, delta_t = lac.d.t + 2),
                  control = nls.lm.control(maxiter=500)
  )
  
  sum <- summary(lac.LM)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  boot <- Boot(lac.LM, method = 'residual')
  
  post <- as.data.frame(boot$t) # Get the bootstrapped values
  post$d.t <- post$delta_t
  
  boot.a <- median(post$a, na.rm = T)    # Extract parameters
  boot.b <- median(post$b, na.rm = T)
  boot.tmax <- median(post$tmax, na.rm = T)
  boot.d.t <- median(post$d.t, na.rm = T)
  
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
                        lower = -100) {
    
    f_low  <- lactin2(temp = lower, a = a, b = b, tmax = tmax, d.t = d.t)
    f_high <- lactin2(temp = Topt,  a = a, b = b, tmax = tmax, d.t = d.t)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2(temp, a, tmax, d.t, b),
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
      function(temp) lactin2(temp, a, tmax, d.t, b),
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
  
  lewington.summ.df <- rbind(lewington.summ.df, data.frame(                     # Add summary data
    Sp.id = df.i$id.number[1],                                                  # Species #
    Sp.name = df.i$Species.name[1],                                             # Species name
    
    RSE = sigma(lac.LM),                                                        # Raw residual square error
    RSE.mu = sigma(lac.LM)/max(df.i$mu),                                        # Residual square error standardized against sd of µ values
    RSE.range = sigma(lac.LM)/(max(df.i$mu) - min(df.i$mu)),                    # Residual square error standardized against range of µ values
    RSE.sd = sigma(lac.LM)/sd(df.i$mu),                                         # Residual square error standardized against sd of µ values
    
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
    
    a = boot.a,                                                                 # parameter: a
    b = boot.b,                                                                 # parameter: b
    tmax = boot.tmax,                                                           # parameter: tmax
    d.t = boot.d.t,                                                             # parameter: deltaT
    
    a.mod = df.nls[1,1],                                                        # nls.LM parameter: a
    b.mod = df.nls[2,1],                                                        # nls.LM parameter: b
    tmax.mod = df.nls[3,1],                                                     # nls.LM parameter: tmax
    d.t.mod = df.nls[4,1]                                                       # nls.LM parameter: deltaT
  ))
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(                                         # Model performance data
      Sp.id = df.i$id.number[1],                                                # Species #
      Sp.name = df.i$Species.name[1],                                           # Species name
      Parameter = df.nls$parameter[j],                                          # Model parameter (e.g. a, b, tmax etc.)
      est = df.nls$Estimate[j],                                                 # Estimate
      se = df.nls$`Std. Error`[j],                                              # Error
      p = df.nls$`Pr(>|t|)`[j],                                                 # p-values
      stringsAsFactors = FALSE            
    ))
  }
  
  print(paste("Done", n, "of ", length(unique(df.l$id.number))))
  
}

write.csv(lewington.summ.df, "processed-data/36_Lewington_2019_TPCs.csv") # Lewington_2019 summary table
write.csv(fit.df, "processed-data/37_Lewington_2019_TPCs_fits.csv") # Save model fit summary table

# Edwards 2016 ------------------------------------------------------------

###### Load the data ######

df.e.raw <- read.csv("processed-data/38_Edwards_2016_raw_data.csv") # Raw data file
head(df.e.raw)

df.e.raw %>% # Number of observations by species
  group_by(species) %>%
  summarise(n = n(), .groups = "drop") %>% 
  print(n = 200) # They all have more than 

df.e.raw <- df.e.raw %>% # The Edwards data has some species that show up in multiple references, so we need to create a unique sp.idx ref# combo
  mutate(unique.id = paste(species, reference, sep = "_")) # So that each entry is treated separately!

df.e.raw$id.number <- as.integer(factor(df.e.raw$unique.id))

df.e.raw <- df.e.raw %>% 
  filter(irradiance<= 250) # After looking at the light models, insanely high lights are driving down growth rates. We'll cap this out at 250 to improve model fits.

df.e <- df.e.raw %>% # Remove species with fewer than 5 data points
  group_by(id.number) %>%
  filter(irradiance == max(irradiance)) %>% 
  mutate(n = n()) %>% 
  filter(n >= 5) 

df.e <-df.e %>% 
  ungroup() %>% 
  rename(temp = temperature,
         mu = growth.rate,
         Species.name = species)

df.e <- df.e %>% # Remove species where there are no growth data above putative Topt. 
  group_by(id.number) %>%
  filter(!(max(temp) %in% temp[mu == max(mu)])) 

length(unique(df.e$id.number)) # Down to 29 out of 63 species....

###### Model fitting ######

edwards.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Study = character(),      # Study
  
  RSE = numeric(),          # Raw residual square error
  RSE.mu = numeric(),       # Residual square error standardized against max µ value
  RSE.range = numeric(),    # Residual square error standardized against range of µ values
  RSE.sd = numeric(),       # Residual square error standardized against sd of µ values
  
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
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Parameter = character(),  # Model parameter (e.g. a, tmax, etc.)
  est = numeric(),          # Estimate
  se = numeric(),           # Standard error
  p = numeric(),            # p-value
  stringsAsFactors = FALSE        
)

n <-0 # progression tracker

for (i in unique(df.e$id.number[df.e$id.number >= 1])) { # for each species ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.e %>% 
    filter(id.number == i) %>% 
    filter(irradiance == max(irradiance)) %>% 
    filter(!is.na(mu))
  
  df.i <- droplevels(df.i)
  
  # df.i <- add_left_zero_if_needed(df.i)
  
  lac_nls <- nls_multstart(mu ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = c(0.01, -2.5, 5, 1),
                           start_upper = c(0.19, -0.5, max(df.i$temp), 11),
                           lower = c(0, -3, min(df.i$temp), 0.1),
                           upper = c(0.5, -0.001, max(df.i$temp) + 5, 40),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  sum <- summary(lac_nls)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  lac.a <- df.nls[1,1] # Extract parameters
  lac.b <- df.nls[2,1]
  lac.tmax <- df.nls[3,1]
  lac.d.t <- df.nls[4,1]
  
  lac.LM <- nlsLM(mu ~ lactin2_1995(temp = temp, a, b, tmax, delta_t),
                  data = df.i,
                  start = c(a = lac.a, b = lac.b, tmax = lac.tmax, delta_t = lac.d.t),
                  lower = c(a = max(0.001,lac.a - 0.05), b = lac.b - 0.5, tmax = max(0, lac.tmax -3), delta_t = max(0.1, lac.d.t - 2)),
                  upper = c(a = lac.a + 0.05, b = min(lac.b + 0.5, -0.001), tmax = lac.tmax + 3, delta_t = lac.d.t + 2),
                  control = nls.lm.control(maxiter=500)
  )
  
  sum <- summary(lac.LM)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  boot <- Boot(lac.LM, method = 'residual')
  
  post <- as.data.frame(boot$t) # Get the bootstrapped values
  post$d.t <- post$delta_t
  
  boot.a <- median(post$a, na.rm = T)    # Extract parameters
  boot.b <- median(post$b, na.rm = T)
  boot.tmax <- median(post$tmax, na.rm = T)
  boot.d.t <- median(post$d.t, na.rm = T)
  
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
                        lower = -100) {
    
    f_low  <- lactin2(temp = lower, a = a, b = b, tmax = tmax, d.t = d.t)
    f_high <- lactin2(temp = Topt,  a = a, b = b, tmax = tmax, d.t = d.t)
    
    # Check for valid root conditions
    if (!is.finite(f_low) || !is.finite(f_high)) return(NA_real_)
    if (f_low * f_high > 0) return(NA_real_)   # no sign change
    
    # Root exists → compute it
    uniroot(
      function(temp) lactin2(temp, a, tmax, d.t, b),
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
      function(temp) lactin2(temp, a, tmax, d.t, b),
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
  
  edwards.summ.df <- rbind(edwards.summ.df, data.frame(                         # Add summary data
    Sp.id = df.i$id.number[1],                                                  # Species #
    Sp.name = df.i$Species.name[1],                                             # Species name
    Study = df.i$reference[1],                                                  # Reference #
    
    RSE = sigma(lac.LM),                                                        # Raw residual square error
    RSE.mu = sigma(lac.LM)/max(df.i$mu),                                        # Residual square error standardized against sd of µ values
    RSE.range = sigma(lac.LM)/(max(df.i$mu) - min(df.i$mu)),                    # Residual square error standardized against range of µ values
    RSE.sd = sigma(lac.LM)/sd(df.i$mu),                                         # Residual square error standardized against sd of µ values
    
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
    
    a = boot.a,                                                                 # parameter: a
    b = boot.b,                                                                 # parameter: b
    tmax = boot.tmax,                                                           # parameter: tmax
    d.t = boot.d.t,                                                             # parameter: deltaT
    
    a.mod = df.nls[1,1],                                                        # nls.LM parameter: a
    b.mod = df.nls[2,1],                                                        # nls.LM parameter: b
    tmax.mod = df.nls[3,1],                                                     # nls.LM parameter: tmax
    d.t.mod = df.nls[4,1]                                                       # nls.LM parameter: deltaT
  ))
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(                                         # Model performance data
      Sp.id = df.i$id.number[1],                                                # Species #
      Sp.name = df.i$Species.name[1],                                           # Species name
      Parameter = df.nls$parameter[j],                                          # Model parameter (e.g. a, b, tmax etc.)
      est = df.nls$Estimate[j],                                                 # Estimate
      se = df.nls$`Std. Error`[j],                                              # Error
      p = df.nls$`Pr(>|t|)`[j],                                                 # p-values
      stringsAsFactors = FALSE            
    ))
  }
  
  print(paste("Done", n, "of ", length(unique(df.e$id.number))))
  
}

write.csv(edwards.summ.df, "processed-data/39_Edwards_2016_TPCs.csv") # Edwards_2016 summary table
write.csv(fit.df, "processed-data/40_Edwards_2016_TPCs_fits.csv") # Save model fit summary table

# Levasseur 2025 ------------------------------------------------------------

###### Load the data ######

df.l <- read.csv("processed-data/41_Levasseur2025_l_rawdata.csv") # Raw data file
head(df.l)
str(df.l)

df.l$sp.num <- as.integer(factor(df.l$species_updated))
df.l$log.RFU <- log(df.l$RFU + 0.001)

df.lt <- df.l %>% 
  filter(value_resource == max(value_resource, na.rm = TRUE))

df.n <- read.csv("processed-data/42_Levasseur2025_n_rawdata.csv") # Raw data file
head(df.n)
str(df.n)

df.n$sp.num <- as.integer(factor(df.n$species_updated))
df.n$log.RFU <- log(df.n$RFU + 0.001)

df.nt <- df.n %>% 
  filter(value_resource == max(value_resource, na.rm = TRUE))

df.p <- read.csv("processed-data/43_Levasseur2025_p_rawdata.csv") # Raw data file
head(df.p)
str(df.p)

df.p$sp.num <- as.integer(factor(df.p$species_updated))
df.p$log.RFU <- log(df.p$RFU + 0.001)

df.pt <- df.p %>% 
  filter(value_resource == max(value_resource, na.rm = TRUE))

df.t <- rbind(df.lt, df.nt, df.pt)

df.t <- df.t %>% # We'll concatenate well id and resource level just in case the same well appears in different experiments
  mutate(unique.id = paste(well_id, level_of_resource, sep = "_")) # So that each entry is treated seperately!

###### Estimate µ ######

temp <- as.vector(unique(df.t$temperature))# for looping through temps
ord.temp<- sort(temp)

df.r.exp <- data.frame(                # Summary dataframe for r_exp
  Sp.id = character(),                 # Species
  temperature = numeric(),             # Temperature
  id = character(),                    # Unique ID
  r.exp = numeric()                    # Thresholded r.exp
)

for (i in 1:20){ # for every species
  
  for (t in ord.temp){ # at every temperature
    
    df.i <- df.t[df.t$sp.num == i, ]
    df.it <- df.i[df.i$temperature == t, ]
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in unique(df.it$unique.id)){ # Doing this separately for each replicate
      
      df.it.wl <- subset(df.it, df.it$unique.id == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$Time_days), ]
      df.it.wl <- df.it.wl %>% 
        mutate(N0 = RFU[1])
      
      t.series <- unique(df.it.wl$Time_days) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$Time_days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(log.RFU~Time_days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
      s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
      
      df.it.wl.th <- df.it.wl[df.it.wl$Time_days <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
      # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*Time_days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp <- rbind(df.r.exp, data.frame(
          Sp.id = df.it.wl.th$species_updated[1],           
          temperature = df.it.wl$temperature[1],        
          id = df.it.wl$sp.num[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          Sp.id = df.it.wl.th$species_updated[1],    # Species
          temperature = df.it.wl$temperature[1],     # Temperature
          id = df.it.wl$sp.num[1],                   # ID
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
      }
      
    }
    
  }
  
}

write.csv(df.r.exp, "processed-data/44_Levasseur_2025_µ_estimates_temp.csv") # let's save the file.

###### Fit TPCs ######

# df.r.exp <- read.csv("processed-data/44_Levasseur_2025_µ_estimates_temp.csv") # if needed

levasseur.t.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species ID
  Sp.name = character(),    # Species name
  T.min = numeric(),        # Minimum T (calculus)
  T.max = numeric(),        # Maximum T (calculus)
  T.opt = numeric(),        # Optimal T (calculus)
  r.max = numeric(),        # Maximum growth rate (calculus)
  T.br = numeric(),         # T breadth (calculus)
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  est = numeric(),          # Estimate
  se = numeric(),           # Standard error
  p = numeric(),            # p-value
  stringsAsFactors = FALSE                  
)

for (i in 1:20){ # For each species
  
  df.i <- df.r.exp %>% 
    filter(id == i, !is.na(r.exp)) %>% 
    arrange(temperature)
  df.i <- droplevels(df.i)
  
  lac_nls <- nls_multstart(r.exp ~ lactin2_1995(temp = temperature, a, b, tmax, delta_t),
                           data = df.i,
                           iter = c(4, 4, 4, 4), 
                           start_lower = c(0.01, -2.5, 5, 1),
                           start_upper = c(0.19, -0.5, max(df.i$temperature), 11),
                           lower = c(0, -3, min(df.i$temperature), 0.1),
                           upper = c(0.5, 0, max(df.i$temperature) + 5, 40),
                           supp_errors = 'Y',
                           convergence_count = FALSE
  )
  
  sum <- summary(lac_nls)
  
  df.nls <- as.data.frame(sum$coefficients)
  df.nls$parameter <- rownames(df.nls)
  rownames(df.nls) <- NULL
  
  df.nls
  
  cf.a <- df.nls[1,1] # Extract parameters
  cf.b <- df.nls[2,1]
  cf.tmax <- df.nls[3,1]
  cf.delta_t <- df.nls[4,1]
  
  # Find the T_opt: where the derivative crosses zero
  T_opt <- uniroot(
    function(temp) lactin2_deriv(temp, cf.a, cf.b, cf.tmax, cf.delta_t),
    interval = c(-10, 45)
  )$root
  
  r_max <- lactin2(temp=T_opt, cf.a=cf.a, cf.b=cf.b, cf.tmax=cf.tmax, cf.delta_t=cf.delta_t)
  
  Tmin <- uniroot(lactin2, interval = c(-1000000000, T_opt), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t)$root
  
  Tmax <- uniroot(lactin2, interval = c(T_opt,45), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t)$root
  
  r_half <- r_max/2 # calculate half of rmax and get the roots.
  
  Tlow <- uniroot(lactin2_halfmax, interval = c(Tmin, T_opt), cf.a = cf.a, 
                  cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t, r_half = r_half)$root
  
  Thigh <- uniroot(lactin2_halfmax, interval = c(T_opt, Tmax), cf.a = cf.a, 
                   cf.b = cf.b, cf.tmax = cf.tmax, cf.delta_t = cf.delta_t, r_half = r_half)$root
  
  levasseur.t.summ.df <- rbind(levasseur.t.summ.df, data.frame( # Add summary data
    Sp.id = df.i$id[1],                               # Species #
    Sp.name = df.i$Sp.id[1],                          # Species name
    T.min = Tmin,                                     # Minimum T (calculus)
    T.max = Tmax,                                     # Maximum T (calculus)
    T.opt = T_opt,                                    # Optimal T (calculus)
    r.max = r_max,                                    # Maximum growth rate (calculus)
    T.br = Thigh - Tlow                               # T breadth (calculus)
  ))
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(                         # Model performance data
      Sp.id = df.i$id[1],                                       # Species #
      Sp.name = df.i$Sp.id[1],                                  # Species name
      Parameter = df.nls$parameter[j],                          # Model parameter (e.g. a, b, tmax etc.)
      est = df.nls$Estimate[j],                                 # Estimate
      se = df.nls$`Std. Error`[j],                              # Error
      p = df.nls$`Pr(>|t|)`[j],                                 # p-values
      stringsAsFactors = FALSE            
    ))
  }
  
}

write.csv(levasseur.t.summ.df, "processed-data/45_Levasseur_2025_TPCs.csv") # Save Levasseur 2025 TPC summary table
write.csv(fit.df, "processed-data/46_Levasseur_2025_TPCs_fits.csv") # Save model fit summary table
