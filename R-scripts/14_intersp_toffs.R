# Jason R Laurich
# March 5, 2025

# Bringing in interspecific data

############# Packages ########################

library(dplyr)
library(ggplot2)
library(rPref)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(vegan)  # For PCA and RDA
library(ggrepel)
library(quantreg)
library(brms)

############# Upload and organize data #######################

df<-read.csv("data-processed/16a_thomas_2012_supp_data_summ.csv") # Summary file
head(df)
str(df)

par_frt <- function(df, xvar, yvar) { # Simple Pareto frontier function
  
  df <- df[order(-df[[xvar]], df[[yvar]]), ]  
  pareto_points <- df[1, ]  # Start with the first point
  
  for (i in 2:nrow(df)) {
    if (df[i, yvar] > tail(pareto_points[[yvar]], 1)) {  # Ensure increasing y values
      pareto_points <- rbind(pareto_points, df[i,])
    }
  }
  
  return(pareto_points)
}

############# Quick examination of the relationship #####################

par.res <- par_frt(df, xvar = "Optimum", yvar = "Niche.width")

T_thomas <- ggplot(df, aes(x = Optimum, y = Niche.width)) +  # Remove shape from aes() for regression
  geom_point( size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Optimal temperature", 
       y = "Thermal niche breadth", 
       title = "Thermal performance") +
  geom_line(data = par.res, aes(x = Optimum, y = Niche.width), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

T_thomas # Raw pareto front.
