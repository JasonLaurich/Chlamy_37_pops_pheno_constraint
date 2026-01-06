# Jason R Laurich

# December 28th, 2025

# We are going to reestimate Thomas 2012 TPCs using nls.multstart

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(rTPC)
library(nls.multstart)
library(Deriv)

lactin2 <- function(temp, cf.a, cf.tmax, cf.delta_t, cf.b) { 
  exp(cf.a * temp) - exp(cf.a * cf.tmax - (cf.tmax - temp) / cf.delta_t) + cf.b
} # Define the Lactin II function

lactin2_deriv <- function(temp, cf.a, cf.b, cf.tmax, cf.delta_t) {
  rho <- cf.a
  T_max <- cf.tmax
  delta_T <- cf.delta_t
  
  term1 <- rho * exp(rho * temp)
  term2 <- (1 / delta_T) * exp(rho * T_max - (T_max - temp) / delta_T)
  
  return(term1 - term2)
} # Derivative of the Lactin II function

lactin2_halfmax <- function(temp, cf.a, cf.b, cf.tmax, cf.delta_t, r_half) {
  exp(cf.a * temp) - exp(cf.a * cf.tmax - (cf.tmax - temp) / cf.delta_t) + cf.b - r_half
} # OK we're going to modify the function to calculate T_breadth, based on a modified lactin.

# Load & examine the data ------------------------------------------------------------

df.t.raw <- read.csv('data-processed/10_Thomas_2012_raw_data.csv') # Thomas raw data
head(df.t.raw)

length(unique(df.t.raw$id.number))

mat <- split(df.t.raw, df.t.raw$id.number)  # Matrix

min(df.t.raw[df.t.raw$Growth.rate>0,]$Temperature) # Tmin => -1.8
max(df.t.raw[df.t.raw$Growth.rate>0,]$Temperature) # Tmax <= 37

# Looping -----------------------------------------------------------------

thomas.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Study = character(),      # Study
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
  Study = character(),      # Study
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  est = numeric(),          # Estimate
  se = numeric(),           # Standard error
  p = numeric(),            # p-value
  stringsAsFactors = FALSE            
)

for (i in c(1:62, 64:length(mat))){ # for each species. Can't do species 63 - data does not support estimating Tmin
  
  df.i <- subset(mat[[i]])
  df.i <- droplevels(df.i)
  
  lac_nls <- nls_multstart(Growth.rate ~ lactin2_1995(temp = Temperature, a, b, tmax, delta_t),
                             data = df.i,
                             iter = c(4, 4, 4, 4), 
                           start_lower = c(0.01, -2.5, 5, 1),
                           start_upper = c(0.19, -0.5, max(df.i$Temperature), 11),
                           lower = c(0, -3, min(df.i$Temperature), 0.1),
                           upper = c(0.5, 0, max(df.i$Temperature) + 5, 40),
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
  
  thomas.summ.df <- rbind(thomas.summ.df, data.frame( # Add summary data
    Sp.id = df.i$id.number[1],                        # Species #
    Sp.name = df.i$Species.name[1],                   # Species name
    Study = df.i$Study[1],                            # Study
    T.min = Tmin,                                     # Minimum T (calculus)
    T.max = Tmax,                                     # Maximum T (calculus)
    T.opt = T_opt,                                    # Optimal T (calculus)
    r.max = r_max,                                    # Maximum growth rate (calculus)
    T.br = Thigh - Tlow                               # T breadth (calculus)
  ))
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(                         # Model performance data
      Sp.id = df.i$id.number[1],                                # Species #
      Sp.name = df.i$Species.name[1],                           # Species name
      Study = df.i$Study[1],                                    # Study
      Parameter = df.nls$parameter[j],                          # Model parameter (e.g. a, b, tmax etc.)
      est = df.nls$Estimate[j],                                 # Estimate
      se = df.nls$`Std. Error`[j],                              # Error
      p = df.nls$`Pr(>|t|)`[j],                                 # p-values
      stringsAsFactors = FALSE            
    ))
  }
  
  message(sprintf("Done %d of 194", i))
  
}

write.csv(thomas.summ.df, "data-processed/305a_Thomas_2012_TPCs_newest.csv") # Save Thomas 2012 summary table
write.csv(fit.df, "data-processed/305b_Thomas_2012_TPCs_fits_newest.csv") # Save model fit summary table
