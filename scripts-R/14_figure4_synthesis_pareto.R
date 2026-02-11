# Jason R Laurich

# February 11th, 2026

# Here I am going to import the various datasets that contain synthesis data for phytoplankton growth across
  # light, nitrogen, phosphorous, and temperature gradients.

# I will then apply analyses based on Shoval 2012 (Science) to calculate the triangularity of our data and compare this 
  # estimate to null randomizations. I will further integrate the Pareto front analyses I have been employing elsewhere
  # in this data analysis pipeline (scripts 08/09)

# Inputs: 47_Edwards_2016_light_monods.csv, 50_Lewington_2019_light_monods.csv, 53_Levasseur_2025_light_monods.csv, 66_Narwani_2015_summary.csv,
  # 68_Edwards_2015_summary.csv, 56_Lewington_2019_nit_monods.csv, 59_Levasseur_2025_nit_monods.csv, 61_Bestion_2018_phos_monods.csv,
  # 64_Levasseur_2025_phos_monods.csv, 32_Bestion_2018_TPCs.csv, 39_Edwards_2016_TPCs.csv, 36_Lewington_2019_TPCs.csv, 45_Levasseur_2025_TPCs.csv,
  # 29_Thomas_2012_TPCs.csv
# Outputs: in processed-data : 67_light_metadata.csv, 69_nit_metadata.csv, 70_phos_metadata.csv,

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(pracma) # For calculating polygon area
library(cowplot)

library(quantreg)
library(sp) # For the point.in.polygon function
library(scam)

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

# Load and compile the data sets -------------------------------------------------------------------

###### Light synthesis data ######

# Edwards 2016

df.l.e <- read.csv("processed-data/47_Edwards_2016_light_monods.csv")
head(df.l.e)

df.l.e$dataset <- "Edwards et al., 2016"

# Lewington-Pearce 2019

df.l.l <- read.csv("processed-data/50_Lewington_2019_light_monods.csv")
head(df.l.l)

df.l.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.l.lv <- read.csv("processed-data/53_Levasseur_2025_light_monods.csv")
head(df.l.lv)

df.l.lv$dataset <- "Levasseur et al., 2025"

# Narwani 2015

df.l.n <- read.csv("processed-data/66_Narwani_2015_summary.csv")
head(df.l.n)

df.l.n$dataset <- "Narwani et al., 2015"

df.l.n <- df.l.n %>% 
  filter(!is.na(Istar))

df.l <- bind_rows(
  
  df.l.e %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.l.l %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.l.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.l.n %>% 
    transmute(Species = Species.name,
              dataset = dataset,
              r.max = umax.light,
              R = Istar,
              comp = 1/Istar)
  
)

head(df.l) # 131 entries!

df.l <- df.l %>% 
  filter(R > 0) # 130 viable entries

plot(df.l$comp ~ df.l$r.max) # Some really wonky (high) estimates of 1/R* that suggest error (filter at > 3.5)

df.l <- df.l %>% 
  filter(comp < 3.5,
         r.max < 2.5,
         !is.na(comp))

write.csv(df.l, "processed-data/67_light_metadata.csv") # Summary

###### Nitrogen synthesis data ######

# Edwards 2015

df.n.e <- read.csv("processed-data/68_Edwards_2015_summary.csv")
head(df.n.e)

df.n.e <- df.n.e %>% 
  filter(!is.na(k_nit_m),
         !is.na(mu_nit)) %>% 
  mutate(R = 0.1*k_nit_m/(mu_nit - 0.1))

df.n.e$dataset <- "Edwards et al., 2015"

# Lewington-Pearce 2019

df.n.l <- read.csv("processed-data/56_Lewington_2019_nit_monods.csv")
head(df.n.l)

df.n.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.n.lv <- read.csv("processed-data/59_Levasseur_2025_nit_monods.csv")
head(df.n.lv)

df.n.lv$dataset <- "Levasseur et al., 2025"

df.n <- bind_rows(
  
  df.n.e %>% 
    transmute(Species = species,
              dataset = dataset,
              r.max = mu_nit,
              R = R,
              comp = 1/R),
  
  df.n.l %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.n.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post)
  
)

head(df.n) # 58 entries!

df.n <- df.n %>% 
  filter(R > 0) # 57 viable entries

plot(df.n$comp ~ df.n$r.max) # Some really wonky (high) estimates of 1/R* that suggest error or nitrogen fixation?

df.n <- df.n %>% 
  filter(comp < 60,
         !is.na(r.max),
         !is.na(comp))

write.csv(df.n, "processed-data/69_nit_metadata.csv") # Summary

###### Phosphorous synthesis data ######

# Bestion 2018

df.p.b <- read.csv("processed-data/61_Bestion_2018_phos_monods.csv")
head(df.p.b)

df.p.b$dataset <- "Bestion et al., 2018"

# Edwards 2015

df.p.e <- read.csv("processed-data/68_Edwards_2015_summary.csv")
head(df.p.e)

df.p.e <- df.p.e %>% 
  filter(!is.na(k_p_m),
         !is.na(mu_p)) %>% 
  mutate(R = 0.1*k_p_m/(mu_p - 0.1))

df.p.e$dataset <- "Edwards et al., 2015"

# Levasseur 2025

df.p.lv <- read.csv("processed-data/64_Levasseur_2025_phos_monods.csv")
head(df.p.lv)

df.p.lv$dataset <- "Levasseur et al., 2025"

# Narwani 2015

df.p.n <- read.csv("processed-data/66_Narwani_2015_summary.csv")
head(df.p.n)

df.p.n$dataset <- "Narwani et al., 2015"

df.p <- bind_rows(
  
  df.p.b %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.p.e %>% 
    transmute(Species = species,
              dataset = dataset,
              r.max = mu_p,
              R = R,
              comp = 1/R),
  
  df.p.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.p.n %>% 
    transmute(Species = Species.name,
              dataset = dataset,
              r.max = umax.nitrate,
              R = Nstar,
              comp = 1/Nstar)
  
)

head(df.p) # 156 entries!

df.p <- df.p %>% 
  filter(R > 0) # 124 viable entries (NAs in the narwani data set)

plot(df.p$comp ~ df.p$r.max) # Some really wonky (high) estimates of 1/R* that suggest error or nitrogen fixation?

df.p <- df.p %>% 
  filter(comp < 300,
         r.max < 2.5,
         !is.na(r.max),
         !is.na(comp))

write.csv(df.p, "processed-data/70_phos_metadata.csv") # Summary

###### Temperature synthesis data ######

# Bestion 2018

df.t.b <- read.csv("processed-data/32_Bestion_2018_TPCs.csv")
head(df.t.b)

df.t.b$dataset <- "Bestion et al., 2018"

# Edwards 2016

df.t.e <- read.csv("processed-data/39_Edwards_2016_TPCs.csv")
head(df.t.e)

df.t.e$dataset <- "Edwards et al., 2016"

# Lewington-Pearce 2019

df.t.l <- read.csv("processed-data/36_Lewington_2019_TPCs.csv")
head(df.t.l)

df.t.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.t.lv <- read.csv("processed-data/45_Levasseur_2025_TPCs.csv")
head(df.t.lv)

df.t.lv$dataset <- "Levasseur et al., 2025"

# Thomas 2012

df.t.t <- read.csv("processed-data/29_Thomas_2012_TPCs.csv")
head(df.t.t)

df.t.t$dataset <- "Thomas et al., 2012"

df.t <- bind_rows(
  
  df.t.b %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.e %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.l %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.t %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
)

head(df.t) # 229 entries!

plot(df.t$T.br ~ df.t$r.max)

df.t <- df.t %>% 
  filter(r.max < 3.5,
         T.br < 36,
         !is.na(r.max),
         !is.na(T.br))

write.csv(df.t, "processed-data/71_temp_metadata.csv") # Summary
