# Jason R Laurich

# October 2nd, 2025

# I am going to play around with quantile regression

# Load packages & specify functions -----------------------------------------------------------

library(tidyverse)
library(cowplot)
library(quantreg)

# Load & examine the data -------------------------------------------------

df <- read.csv("data-processed/20_summary_table.csv") # Summary file
str(df)
head(df)
