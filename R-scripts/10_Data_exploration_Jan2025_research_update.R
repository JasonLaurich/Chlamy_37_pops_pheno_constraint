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

########################## Visualization #####################################

df<-read.csv("data-processed/14_summary_metric_table.csv")

par_frt <- function(df, xvar, yvar) { # general pareto frontier function
  
  df <- df[order(-df[[xvar]], df[[yvar]]), ]  
  pareto_points <- df[1, ]  # Start with the first point
  
  for (i in 2:nrow(df)) {
    if (df[i, yvar] > tail(pareto_points[[yvar]], 1)) {  # Ensure increasing y values
      pareto_points <- rbind(pareto_points, df[i,])
    }
  }
  
  return(pareto_points)
}

par.res <- par_frt(df[df$I.comp < 10,], xvar = "I.comp", yvar = "T.br")

df$evol.plt <- factor(df$evol, 
                      levels = c("none", "L", "N", "P", "S", "B", "BS", "C"),
                      labels = c("Ancestral", 
                                 "Light limitation", 
                                 "Nitrogen limitation", 
                                 "Phosphorous limitation", 
                                 "Salt stress", 
                                 "Biotic depletion", 
                                 "Biotic depletion x Salt", 
                                 "Control"))


T.I <- ggplot(df[df$I.comp < 10,], aes(x = I.comp, y = T.br)) +  # removing the outlier data point
  geom_point(size = 2) +  # Original data points
  geom_line(data = par.res, aes(x = I.comp, y = T.br), color = "red", size = 1) +  # Pareto frontier line
  labs(x = "Competitive ability (1/I*)", y = "Thermal breadth (T.br)") +
  theme_classic()

T.I

T.I.evol <- T.I + 
  geom_point(aes(colour = evol.plt)) +  # Add the color aesthetic mapping
  scale_colour_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Ancestral" = "black",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Phosphorous limitation" = "firebrick",
               "Salt stress" = "blue",
               "Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen")
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.85, 0.3),  # Move legend to bottom right
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )


T.I.evol


