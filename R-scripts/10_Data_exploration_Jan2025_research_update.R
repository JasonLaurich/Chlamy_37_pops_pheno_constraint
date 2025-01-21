# Jason R Laurich
# January 21, 2025

# OK so here I am going to load up all of my previously summarized datasets.
# I'll add in a few key checks (e.g. matching them up to original data to make sure my population labels are consistent across the board).
# I'll also want to incorporate ancestry and evolutionary environment as factors.
# And I'll want to calculate 1/R* (competitive ability) for N, P, and I.
# Change in R*, TPC values? From ancestral population?

# Then we'll do some PCAs (maybe accounting for ancestry? Evolutionary environment?) and generate some simple plots showcasing differences in TPC, R* shape.

############# Packages ########################

library(dplyr)
library(ggplot2)
library(rPref)

############# Upload and organize data #######################

df.hist <- read.csv("data-processed/chlamee-acute-exponential.csv") # This file has the history (ancestry, evolutionary treatment) for each well ID.
head(df.hist)
str(df.hist)
df.hist$anc <- as.factor(df.hist$ancestor_id) # Will port these into combined dataset, remember: we are sorting by population. 
df.hist$evol <- as.factor(df.hist$treatment)  # so we really only need population level - evol - anc mapping, not well-specific.
df.hist$Pop.fac <- as.factor(df.hist$population)

df.tpc <- read.csv("data-processed/09_TPC_shape_values_BayesTPC.csv") # TPC summary data. We'll use Lactin II model fits. 
head(df.tpc)
str(df.tpc)

df.tpc <- df.tpc[df.tpc$Pop.fac != 'cc1629' & df.tpc$Model == 'Lactin 2', ] # Trim off Thomas models and that one weird population.

df <- df.tpc[,c(2,5:9)]

df.hist.agg <- df.hist %>% # Aggregate the history data for merging.
  group_by(Pop.fac) %>%
  summarise(
    anc = first(anc),
    evol = first(evol),
    .groups = "drop" 
  )

df <- merge(df, df.hist.agg, by = "Pop.fac", all.x = TRUE)
head(df) 
str(df)
df$Pop.fac <- as.factor(df$Pop.fac) # Looks great! Now we layer in other summary metrics.
names(df) <- c("Pop.fac", "T.min", "T.max", "T.br", "Topt", "r.max_T", "anc", "evol") # rename a few variables. 

df.I <- read.csv("data-processed/10b_light_Monod_estimates.csv") # I* data 
head(df.I)
str(df.I)

df.I.par <- df.I[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.I.par) <- c("Pop.fac", "K.s_I", "r.max_I", "I.jag", "I.mth") # rename
df.I.par$I.comp <- 1/df.I.par$I.mth

df <- merge(df, df.I.par, by = "Pop.fac", all.x = TRUE)

df.N <- read.csv("data-processed/11b_nitrogen_Monod_estimates.csv") # N* data 
head(df.N)
str(df.N)

df.N.par <- df.N[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.N.par) <- c("Pop.fac", "K.s_N", "r.max_N", "N.jag", "N.mth") # rename
df.N.par$N.comp <- 1/df.N.par$N.mth

df <- merge(df, df.N.par, by = "Pop.fac", all.x = TRUE)

df.P <- read.csv("data-processed/12b_phosphorous_Monod_estimates.csv") # P* data 
head(df.P)
str(df.P)

df.P.par <- df.P[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.P.par) <- c("Pop.fac", "K.s_P", "r.max_P", "P.jag", "P.mth") # rename
df.P.par$P.comp <- 1/df.P.par$P.mth

df <- merge(df, df.P.par, by = "Pop.fac", all.x = TRUE)

df.S <- read.csv("data-Processed/13b_salt_tolerance_estimates.csv") # S* data 
head(df.S)
str(df.S)

df.S.par <- df.S[, c("Pop.fac", "r.max_S", "c.mod", "c.pred")] # need these
names(df.S.par) <- c("Pop.fac", "r.max_S", "S.c.mod", "S.c.pred") # rename

df <- merge(df, df.S.par, by = "Pop.fac", all.x = TRUE)

write.csv(df, "data-processed/14_summary_metric_table.csv") # Save summary table
