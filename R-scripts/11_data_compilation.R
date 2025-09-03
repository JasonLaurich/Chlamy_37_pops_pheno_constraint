# Jason R Laurich

# September 3rd, 2025

# Loading all of the disparate files that contain relevant ecological data and combining into a single massive data frame.
# Am also going to (1) re-estimate thermal breadth at relevant levels (0.56 for our study, 0.1 for comparison with other datasets) and 
# (2) calculate the 95% HDPI intervals for all estimates using the posterior distributions from our JAGS models!

# Packages & Functions ----------------------------------------------------

library(dplyr)

# Load & examine data -----------------------------------------------------

df.hist <- read.csv("data-processed/16_rfus_time_summary.csv") # This file has the history (ancestry, evolutionary treatment) for each well ID.
head(df.hist)

df.hist$anc <- as.factor(df.hist$ancestor_id) # Will port these into combined dataset, sorting by population. 
df.hist$evol <- as.factor(df.hist$treatment)
df.hist$Pop.fac <- as.factor(df.hist$population)

df.tpc <- read.csv("data-processed/05_TPCs.csv") # TPC summary data. We'll use Lactin II model fits. 
head(df.tpc)

df <- df.tpc[,c(2,5:9)] # We'll work with the analytically determined parameters (Deriv package)
head(df)

df.hist.agg <- df.hist %>% # Aggregate the history data for merging.
  group_by(Pop.fac) %>%
  summarise(
    anc = first(anc),
    evol = first(evol),
    .groups = "drop" 
  )

df <- merge(df.hist.agg, df, by = "Pop.fac", all.x = TRUE)
head(df) 
str(df)

df$Pop.fac <- as.factor(df$Pop.fac) # Looks great! Now we layer in other summary metrics.
names(df) <- c("Pop.fac", "anc", "evol", "T.min", "T.max", "T.br", "Topt", "r.max_T") # rename a few variables. 

# Bring in light data

df.I <- read.csv("data-processed/06b_Monod_light_estimates.csv") # I* data 
head(df.I)
str(df.I)

df.I.par <- df.I[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.I.par) <- c("Pop.fac", "K.s_I", "r.max_I", "I.jag", "I.mth") # rename
df.I.par$I.comp <- 1/df.I.par$I.mth

df <- merge(df, df.I.par, by = "Pop.fac", all.x = TRUE)

# Bring in nitrogen data

df.N <- read.csv("data-processed/07b_Monod_nitrogen_estimates.csv") # N* data 
head(df.N)
str(df.N)

df.N.par <- df.N[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.N.par) <- c("Pop.fac", "K.s_N", "r.max_N", "N.jag", "N.mth") # rename
df.N.par$N.comp <- 1/df.N.par$N.mth

df <- merge(df, df.N.par, by = "Pop.fac", all.x = TRUE)

# Bring in phosphorous data

df.P <- read.csv("data-processed/08b_Monod_phosphorous_estimates.csv") # P* data 
head(df.P)
str(df.P)

df.P.par <- df.P[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.P.par) <- c("Pop.fac", "K.s_P", "r.max_P", "P.jag", "P.mth") # rename
df.P.par$P.comp <- 1/df.P.par$P.mth

df <- merge(df, df.P.par, by = "Pop.fac", all.x = TRUE)

# Bring in salt data

df.S <- read.csv("data-Processed/09b_salt_tolerance_estimates.csv") # S* data 
head(df.S)
str(df.S)

df.S$Pop.fac <- gsub("Anc ", "anc", df.S$Pop.fac) # ancestors are written differently here

df.S.par <- df.S[, c("Pop.fac", "r.max", "c.mod", "c.pred")] # need these
names(df.S.par) <- c("Pop.fac", "r.max_S", "S.c.mod", "S.c.pred") # rename

df <- merge(df, df.S.par, by = "Pop.fac", all.x = TRUE)

# Calculate & store new metrics ---------------------------------------------------




write.csv(df, "data-processed/14_summary_metric_table.csv") # Save summary table

