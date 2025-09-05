# Jason R Laurich

# September 5th, 2025

# Creating Figure 1: representations of model fits, parameter extraction and data comparions for methodological clarity

# Packages & Functions ----------------------------------------------------

library(tidyverse)

# Load & examine data -----------------------------------------------------

df <- read.csv("data-processed/20_summary_table.csv") # Summary file
str(df)
head(df)
