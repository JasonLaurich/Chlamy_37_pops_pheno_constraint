# Jason R Laurich

# February 8th, 2026

# This file will assemble summary data for all C. reinhardtii populations from our experimental evolution project.
# We incorporate niche-determining traits from scripts 01-05.R as well as pigmentation, stoichiometry, and biovolume data.

# Inputs: in processed-data : 04_TPC_summary.csv, 08_light_monod_summary.csv, 12_nit_monod_summary.csv, 16_phos_monod_summary.csv,
  # 20_salt_tol_summary.csv, 
# Outputs: in processed-data : 19_µ_estimates_salt.csv, 20_salt_tol_summary.csv, 21_salt_tol_fits.csv

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(vegan)  # For PCA

# Load & examine the data -------------------------------------------------

# Temperature

df.t <- read_csv("processed-data/04_TPC_summary.csv") # This file has the TPC data for each replicate.
head(df.t)

df.t <- df.t %>% 
  mutate(T.br = T.br.max - T.br.min) %>% 
  select(population, rep.ID, T.min, T.max, T.opt, r.max, T.br, a, b, tmax, d.t)

# Light

df.l <- read_csv("processed-data/08_light_monod_summary.csv") # This file has the light Monod curve data for each replicate.
head(df.l)

df.l <- df.l %>% 
  select(rep.ID, K.s.post, r.max.post, R.post) %>% 
  mutate(comp = 1/R.post)

# Nitrogen

df.n <- read_csv("processed-data/12_nit_monod_summary.csv") # This file has the nitrogen Monod curve data for each replicate.
head(df.n)

df.n <- df.n %>% 
  select(rep.ID, K.s.post, r.max.post, R.post) %>% 
  mutate(comp = 1/R.post)

# Phosphorous

df.p <- read_csv("processed-data/16_phos_monod_summary.csv") # This file has the phosphorous Monod curve data for each replicate.
head(df.p)

df.p <- df.p %>% 
  select(rep.ID, K.s.post, r.max.post, R.post) %>% 
  mutate(comp = 1/R.post)

# Salt 

df.s <- read_csv("processed-data/20_salt_tol_summary.csv") # This file has the salt tolerance curve data for each population. 
head(df.s)

df.s <- df.s %>% 
  select(population, r.max.post, c.mod.post) %>% 
  mutate(comp = 1/R.post)
