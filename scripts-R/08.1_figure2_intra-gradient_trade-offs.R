# Jason R Laurich

# May 26th, 2026

# We are going to play around with a new (simpler) method for inferreing whether not experimental evolution resulted in shifts towards
# pareto-optimal regions of space

# AND test drive some new figure looks. 

# Inputs: 27_summary_table.csv
# Outputs: in figures-main : 02_fig_2_intra-gradient_tradeoffs.jpeg
# in figures-supp: 02_fig_s2_intra-gradient_tradeoffs_qr.jpeg, 03_fig_s3_intra-gradient_tradeoffs_qr_0.33.jpeg, 10_fig_s10_intra-gradient_evol.jpeg
# in figures-misc: 01_intra-gradient_tradeoffs_fig2.jpeg

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(quantreg)
library(sp) # For the point.in.polygon function
library(scam)
library(pracma) # For calculating polygon area
library(cowplot)

par_frt <- function(df, xvar, yvar) { # Simple Pareto front function / convex hull algorithm (one sided)
  
  df <- df[order(-df[[xvar]], df[[yvar]]), ]  
  pareto_points <- df[1, ]  # Start with the first point
  
  for (i in 2:nrow(df)) {
    if (df[i, yvar] >= tail(pareto_points[[yvar]], 1)) {  # Ensure increasing y values
      pareto_points <- rbind(pareto_points, df[i,])
    }
  }
  
  return(pareto_points)
}

find_nearest_index <- function(x, ref_vec) { # Function to find the closest point to the actual r.max value in the pred data frame. 
  which.min(abs(ref_vec - x))
} # For significance testing based on location of points relative to a curve

# Load & examine the data -------------------------------------------------