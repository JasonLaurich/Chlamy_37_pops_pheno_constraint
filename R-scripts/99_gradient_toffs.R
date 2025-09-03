# Jason R Laurich
# March 5, 2025

# Modelling inter-gradient competitive/ thermal performance metrics here.

# Use same methods developed and tested in 12_gen-spec_glean-opp_toffs.R

# Start with examination of PCA, RDA data from 12
# Pick out some relationships to dive into

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

df<-read.csv("data-processed/14_summary_metric_table.csv") # Summary file

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

############# Relationships to look at? ###########################

# T_br ~ .Comps

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "L", 8, 16)) # Ls are now equivalent to 8, for later mapping

df$evol.bin <- ifelse(df$evol == "none", 'ancestral', 
                      ifelse(df$evol == "L", 'light', 'other')) # for testing regressions.

par.res.TL <- par_frt(df[df$I.comp<10,], xvar = "I.comp", yvar = "T.br")

TL_par <- ggplot(df[df$I.comp<10,], aes(x = I.comp, y = T.br)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Competitive ability (1/R*)", 
       y = "Thermal breadth", 
       title = "Temp ~ Light") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +
  geom_line(data = par.res.TL, aes(x = I.comp, y = T.br), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

TL_par # Raw pareto front.




















df$shape <- ifelse(df$evol == "none", 22, 
                  ifelse(df$evol %in% c("S", "BS"), 8, 16)) # Ss and BSs are now equivalent to 8, for later mapping

df$evol.bin <- ifelse(df$evol == "none", 'ancestral', 
                      ifelse(df$evol %in% c("S", "BS"), 'salt', 'other')) # for testing regressions.