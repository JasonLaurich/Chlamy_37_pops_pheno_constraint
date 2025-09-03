# Jason R Laurich
# April 14, 2025

# We are going to generate "final" figures for the project here

####### Important notes.

#1. need to explore varying the buffer more - is 0.3 too high? What about the point exclusion thresholds?
      # Base this on scaled Euclidean distance? By percentage of total variation between min and max?
#2. for the Li et al null simulations - should I add in error?
      # Yes I should - need to extract values for my 6,000 (?) models and take the 95% HDPI intervals. 
#3. I should be setting a max on k for the scam models I think. 6?
#4. For Li et al randomization - null PFs inflating the number of points in PF.2? Skewing results? Work with the convex hulls directly?
      # Yes definetly! Also the null datasets (reasonably) frequently feature too few Pareto-optimal points for scam fits. 

# Load packages, specify functions  -----------------------------------------------------------

library(cowplot)
library(tidyverse)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(vegan)  # For PCA and RDA
library(ggrepel)
library(quantreg)
library(scam)
library(lme4)
library(emmeans)
library(lmerTest)
library(pracma)

par_frt <- function(df, xvar, yvar) { # Simple Pareto front function / crude convex hull algorithm (one sided)
  
  df <- df[order(-df[[xvar]], df[[yvar]]), ]  
  pareto_points <- df[1, ]  # Start with the first point
  
  for (i in 2:nrow(df)) {
    if (df[i, yvar] > tail(pareto_points[[yvar]], 1)) {  # Ensure increasing y values
      pareto_points <- rbind(pareto_points, df[i,])
    }
  }
  
  return(pareto_points)
}

find_nearest_index <- function(x, ref_vec) { # Write a function to find the closest point to the actual r.max value in the pred data frame. 
  which.min(abs(ref_vec - x))
} # For significance testing based on location of points relative to a curve

point_line_distance <- function(x0, y0, x1, y1, x2, y2) {
  numerator <- abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1)
  denominator <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  numerator / denominator
} # Function to compute perpendicular distance from point to line segment

lactin2 <- function(temp, cf.a, cf.b, cf.delta_t, cf.tmax) {
  exp(cf.a * temp) - exp(cf.a * cf.tmax - ((cf.tmax - temp) / cf.delta_t)) + cf.b
}

find_roots <- function(cf.a, cf.b, cf.delta_t, cf.tmax) {
  fn <- function(temp) lactin2(temp, cf.a, cf.b, cf.delta_t, cf.tmax) - target
  
  root1 <- tryCatch(uniroot(fn, lower = 5, upper = cf.tmax - 10)$root, error = function(e) NA)
  root2 <- tryCatch(uniroot(fn, lower = cf.tmax - 10, upper = 45)$root, error = function(e) NA)
  
  c(root1, root2)
}

buffer <- 0.3 # Buffer (in SD) for point inclusion in the pareto fronts.

# Load and examine data ---------------------------------------------------

df <- read.csv("data-processed/14_summary_metric_table.csv") # Summary file

head(df)
str(df)

levels(as.factor(df$Pop.fac)) # Character-coded, for ancestral populations as well. 

df.stoich <- read.csv("data-processed/24_stoich_data.csv") # Stoichiometry data

head(df.stoich)
str(df.stoich)

levels(as.factor(df.stoich$Name)) # These are 1 - 40 (ie I need to match them to their corresponding population names)

df.pig <- read.csv("data-processed/23a_pigment_data.csv") # Pigment data

head(df.pig)
str(df.pig)

levels(as.factor(df.pig$HPLC.Nummer)) # Hmm so these are already numbered 1-37, ie. the missing populations do not feature here. 

length(unique(df.stoich$Name))
length(unique(df.pig$HPLC.Nummer))
length(unique(df$Pop.fac))

# OK so the following #s should be missing: 12, 13, 34, anc1 (35 + 5 ancestors + 1 general anc = 41 total)
# Missing 4 makes the math work!

levels(as.factor(df.stoich$Name))
length(unique(df.stoich$Name)) # Weirdly I have stoichiometry data for 12 and 13, but not 27 and 29? I also have it for 34?
# So we have three extra reads, but are missing 2 — 38 makes sense
# I don't know how to make sense of this I will talk to Joey about this.

levels(as.factor(df.pig$HPLC.Nummer)) # Ok so we can maybe just assume these are ordered?

df.stoich.NP <- read.csv("data-processed/24a_stoich_NP.csv") # Stoichiometry data for N and P

levels(as.factor(df.stoich.NP$Name)) # Ok so we can maybe just assume these are ordered?
length(unique(df.stoich.NP$Name)) # OK that matches better

df.id <- read.csv("data-processed/22_stoich_pigment_ID_mapping.csv") # File containing the identity assignment

df.id

df.stoich.NP %>% 
  group_by(Name) %>% 
  summarize(mean.N.µg.l = mean(N.µg.l), mean.P.µg.l = mean(P.µg.l)) %>% 
  print(n=37)

df.stoich.NP.sum <- df.stoich.NP %>%    # we have two observations per population, we need to take the means
  group_by(Name) %>% 
  summarize(mean.N.µg.l = mean(N.µg.l), mean.P.µg.l = mean(P.µg.l)) %>% 
  print(n=37)

df.stoich.NP.sum <- df.stoich.NP.sum %>% # Join this with the df.id matching schema
  left_join(df.id, by = c("Name" = "sample"))

df.pig.sum <- df.pig %>%
  mutate(population = df.id$population)

df.final <- df %>% # Make a beautiful final data frame
  left_join(
    df.pig.sum %>% 
      select(chl.a, chl.b, luthein, population), by = c("Pop.fac" = "population")) %>%
  left_join(
    df.stoich.NP.sum %>% 
      select(mean.N.µg.l, mean.P.µg.l, population), by = c("Pop.fac" = "population"))

df.final$evol.plt <- factor(df$evol, 
                      levels = c("none", "L", "N", "P", "S", "B", "BS", "C"),
                      labels = c("Ancestral", 
                                 "Light limitation", 
                                 "Nitrogen limitation", 
                                 "Phosphorous limitation", 
                                 "Salt stress", 
                                 "Biotic depletion", 
                                 "Biotic depletion x Salt", 
                                 "Control"))

df.final$anc.plt <- factor(df$anc, 
                            levels = c("anc2", "anc3", "anc4", "anc5", "cc1690"),
                            labels = c("Population 2", 
                                       "Population 3", 
                                       "Population 4", 
                                       "Population 5", 
                                       "Mixed population"))

df.tpc <- read_csv('data-processed//09.5_TPC_Bayes_model_fit_stats.csv') # Load the data with the TPC shape parameters

df.tpc <- df.tpc %>% # Pivot to wide format
  filter(Model == 'Lactin 2', Pop.fac != 'cc1629') %>% 
  select(Pop.fac, Parameter, mean) %>%
  pivot_wider(names_from = Parameter, values_from = mean)

target <- 0.56

df.tpc.roots <- df.tpc %>%
  mutate(roots = pmap(list(cf.a, cf.b, cf.delta_t, cf.tmax), find_roots)) %>%
  transmute(
    Pop.fac,
    T.min.0.56 = map_dbl(roots, 1),
    T.max.0.56 = map_dbl(roots, 2)
  )

df.final <- df.final %>%
  left_join(
    df.tpc.roots %>% 
      select(Pop.fac, T.min.0.56, T.max.0.56), 
    by = "Pop.fac"
  ) %>%
  mutate(T.br.0.56 = T.max.0.56 - T.min.0.56) %>% 
  mutate(T.br.0 = T.max - T.min)

df.final <- df.final %>%
  mutate(I.mth = I.mth*2.5) %>% 
  mutate(I.comp = 1/I.mth) %>% # Convert light %s to µmol photons (consistent with other interspecific datasets)
  mutate(K.s_I = K.s_I*2.5) %>% 
  mutate(I.comp.10 = 0.1*K.s_I/(r.max_I - 0.1)) %>% 
  mutate(N.comp.10 = 0.1*K.s_N/(r.max_N - 0.1)) %>% 
  mutate(P.comp.10 = 0.1*K.s_P/(r.max_P - 0.1)) # Recalculate 1/R*s at m = 0.1 (for comparison with interspecific datasets)

# Statistical analysis ----------------------------------------------------

df.anc <- df.final %>%   # Isolate ancestral values
  filter(Pop.fac %in% unique(df.final$anc)) %>%
  select(Pop.fac, T.br.0.56, r.max_T, r.max_I, I.comp, r.max_N, N.comp, r.max_P, P.comp, r.max_S, S.c.mod) %>%
  rename(
    anc = Pop.fac,
    T.br.0.56_anc = T.br.0.56,
    r.max_T_anc = r.max_T,
    r.max_I_anc = r.max_I,
    I.comp_anc = I.comp,
    r.max_N_anc = r.max_N,
    N.comp_anc = N.comp,
    r.max_P_anc = r.max_P,
    P.comp_anc = P.comp,
    r.max_S_anc = r.max_S,
    S.c.mod_anc = S.c.mod
  )

df.delta <- df.final %>%   # Calculate parameter changes
  left_join(df.anc, by = "anc") %>%
  mutate(
    delta_T.br.0.56 = T.br.0.56 - T.br.0.56_anc,
    delta_r.max_T = r.max_T - r.max_T_anc,
    delta_r.max_I = r.max_I - r.max_I_anc,
    delta_I.comp = I.comp - I.comp_anc,
    delta_r.max_N = r.max_N - r.max_N_anc,
    delta_N.comp = N.comp - N.comp_anc,
    delta_r.max_P = r.max_P - r.max_P_anc,
    delta_P.comp = P.comp - P.comp_anc,
    delta_r.max_S = r.max_S - r.max_S_anc,
    delta_S.c.mod = S.c.mod - S.c.mod_anc
  ) %>% 
  filter(evol != 'none')

df.delta <- df.delta %>% # Convert to factors
  mutate(evol = factor(evol), anc = factor(anc))

traits <- c("delta_T.br.0.56", "delta_r.max_T", "delta_r.max_I", "delta_I.comp",
            "delta_r.max_N", "delta_N.comp", "delta_r.max_P", "delta_P.comp",
            "delta_r.max_S", "delta_S.c.mod")  # Traits to analyze

stat.sum.df <- data.frame(  # We'll create a dataframe to store model results
  trait = character(),      # trait
  evol.p = numeric(),       # p-value (evolutionary environment)
  anc.prop.var = numeric(), # proportion of variance explained by ancestry
  em.B = numeric(),         # emmean (Biotic depletion)
  em.BS = numeric(),        # emmean (Biotic depletion x Salt)
  em.C = numeric(),         # emmean (Control)
  em.L = numeric(),         # emmean (Light)
  em.N = numeric(),         # emmean (Nitrogen)
  em.P = numeric(),         # emmean (Phosphorous)
  em.S = numeric(),         # emmean (Salt)
  stringsAsFactors = FALSE  # Avoid factor conversion
)

for (i in traits){
  
  model <- lmer(as.formula(paste(i, "~ evol + (1|anc)")), data = df.delta)
  
  anova.mod <- anova(model)
  
  var.cor <- as.data.frame(VarCorr(model))
  
  emm <- as.data.frame(emmeans(model, "evol")) # emmeans
  
  stat.sum.df <- rbind(stat.sum.df, data.frame(                            # Add model outputs
    trait = i,                                                             # trait
    evol.p = anova.mod$`Pr(>F)`,                                           # p-value (evolutionary environment)
    anc.prop.var = var.cor$vcov[1]/(var.cor$vcov[1] + var.cor$vcov[2]),    # proportion of variance explained by ancestry
    em.B = emm$emmean[1],                                                  # emmean (Biotic depletion)
    em.BS = emm$emmean[2],                                                 # emmean (Biotic depletion x Salt)
    em.C = emm$emmean[3],                                                  # emmean (Control)
    em.L = emm$emmean[4],                                                  # emmean (Light)
    em.N = emm$emmean[5],                                                  # emmean (Nitrogen)
    em.P = emm$emmean[6],                                                  # emmean (Phosphorous)
    em.S = emm$emmean[7]                                                   # emmean (Salt)
  )
  )
}

write.csv(stat.sum.df, "data-processed/25_stats_effects_evol_anc_traits.csv") # Save stats results table

# Figure 1: Representation of model fits ----------------------------------

df.final %>% 
  filter(T.br.0.56 %in% c(min(T.br.0.56), max(T.br.0.56), median(T.br.0.56))) %>% 
  select(Pop.fac, T.br.0.56, r.max_T, S.c.mod, r.max_S, P.comp, r.max_P) %>% 
  print() 
  
# Ok so we are going to work with 3 populations - 4, 9, and 11. These represent a good spread
# So number conversions: 4 - 27, 9 - 32, 11 - 3 (pop.num, also how they are recorded in the jags objects!)

# Temperature

for (i in c(3, 27, 32)){      # Temperature R2jags (model fits)
  load(paste0("R2jags-objects/pop_", i, "_lactin2.RData"))
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:6),]
  df.jags.plot$temp <- seq(0, 45, 0.05)
  assign(paste0("df.T.jags", i), df.jags.plot)
}

df.r.t <- read.csv("data-processed/05_final_r_estimates.csv") # Growth data across temperatures
head(df.r.t) 

df.r.t <- df.r.t %>% 
  filter(population.number %in% c(3, 27, 32)) %>% 
  print()

p.t1 <- ggplot(df.r.t[df.r.t$population.number==3,], aes(x = temperature, y = r.exp, colour = population)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-2, 5) +
  xlim(0, 45) +
  
  geom_line(data = df.T.jags3, aes(x = temp, y= mean), colour = "forestgreen", linewidth = 1) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "A - Temperature") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) +
  
  geom_segment(aes(x = -Inf, xend = 34.16441, y = 3.968991, yend = 3.968991), linetype = "dashed", colour = "black", size = 1.2) +
  geom_segment(aes(x = 12.5, xend = 12.5 + df.final[df.final$Pop.fac == "3",]$T.br.0.56, 
                   y = 0.56, yend = 0.56),
              linetype = "dashed", colour = "black", size = 1.2) +
                    
  
  annotate("text", x = 30, y = 4.6, label = "µ max", size = 4.5, fontface = "bold") +
  annotate("text", x = 30, y = 0.56 + 0.6 , label = "Thermal breadth", size = 4.5, fontface = "bold")
  
p.t1

p.t2 <- ggplot(df.r.t[df.r.t$population.number==3,], aes(x = temperature, y = r.exp, colour = population)) +
  ylim(-2, 5) +
  xlim(0, 45) +
  
  geom_line(data = df.T.jags3, aes(x = temp, y= mean), colour = "forestgreen", size = 1) +
  geom_line(data = df.T.jags27, aes(x = temp, y= mean), colour = "magenta2", size = 1) +
  geom_line(data = df.T.jags32, aes(x = temp, y= mean), colour = "darkorange1", size = 1) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "B - Temperature") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 
  
p.t2

p.t3 <- ggplot(df.final[df.final$Pop.fac %in% c("11", "4", "9"),], aes(x = r.max_T, y = T.br, colour = Pop.fac)) +
  geom_point(size=4.5) +
  scale_colour_manual(values = c("forestgreen", "magenta", "darkorange")) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Thermal breadth (°C)",
       title = "C - Temperature") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 
  
p.t3

# Phosphorous

for (i in c(3, 27, 32)){      # Phosphorous
  load(paste0("R2jags-objects/pop_", i, "_phosphorous_monod.RData"))
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:3, 2005),]
  df.jags.plot$phos <- seq(0, 50, 0.025)
  assign(paste0("df.P.jags", i), df.jags.plot)
}

df.r.p <- read.csv("data-processed/12a_phosphorous_r_estimates.csv") # Growth data across P levels
head(df.r.p) # population is the factor, population.number is the corresponding # (e.g. 7, 18, 37)

df.r.p <- df.r.p %>% 
  filter(population.number %in% c(3, 27, 32)) %>% 
  print()

p.p1 <- ggplot(df.r.p[df.r.p$population.number==3,], aes(x = phos.lvl, y = r.exp, colour = population)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = df.P.jags3, aes(x = phos, y= mean), colour = "forestgreen", size = 1) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Phosphorous concentration (µM)", 
       y = "Exponential growth rate (µ)",
       title = "D - Phosphorous") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) +
  
  geom_segment(aes(x = -Inf, xend = 50, y = 1.452906, yend = 1.452906), linetype = "dashed", colour = "black", size = 1.2) +
  geom_segment(aes(x = 1.3431559, xend = 1.3431559, 
                   y = -Inf, yend = 0.56),
               linetype = "dashed", colour = "black", size = 1.2) +
  
  
  annotate("text", x = 25, y = 1.6, label = "µ max", size = 4.5, fontface = "bold") +
  annotate("text", x = 4, y = 0.3, label = "P*", size = 4.5, fontface = "bold")

p.p1

p.p2 <- ggplot(df.r.p[df.r.p$population.number==3,], aes(x = phos.lvl, y = r.exp, colour = population)) +
  scale_colour_manual(values = c("darkorange1", "magenta2", "forestgreen")) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = df.P.jags3, aes(x = phos, y= mean), colour = "forestgreen", size = 1) +
  geom_line(data = df.P.jags27, aes(x = phos, y= mean), colour = "magenta2", size = 1) +
  geom_line(data = df.P.jags32, aes(x = phos, y= mean), colour = "darkorange1", size = 1) +
  
  labs(x = "Phosphorous concentration (µM)", 
       y = "Exponential growth rate (µ)",
       title = "E - Phosphorous") +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.p2

p.p3 <- ggplot(df.final[df.final$Pop.fac %in% c("11", "4", "9"),], aes(x = r.max_P, y = P.comp, colour = Pop.fac)) +
  geom_point(size=4.5) +
  scale_colour_manual(values = c("forestgreen", "magenta", "darkorange")) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/P*)",
       title = "F - Phosphorous") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.p3

# Salt

for (i in c(3, 27, 32)){      # Salt
  load(paste0("R2jags-objects/pop_", i, "_salt_tolerance.RData"))
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:4, 2006),]
  df.jags.plot$salt <- seq(0, 10, 0.005)
  assign(paste0("df.S.jags", i), df.jags.plot)
}

df.r.s <- read.csv("data-processed/13a_salt_r_estimates.csv") # Growth data across salt levels
head(df.r.s) # population is the factor, population.number is the corresponding # (e.g. 7, 18, 37)

df.r.s <- df.r.s %>% 
  filter(population.number %in% c(3, 27, 32)) %>% 
  print()

p.s1 <- ggplot(df.r.s[df.r.s$population.number==3,], aes(x = salt.lvl, y = r.exp, colour = population)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = df.S.jags3, aes(x = salt, y= mean), colour = "forestgreen", size = 1) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Salt concentration (g/L)", 
       y = "Exponential growth rate (µ)",
       title = "G - Salt") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) +
  
  geom_segment(aes(x = -Inf, xend = 10, y = 1.102104, yend = 1.102104), linetype = "dashed", colour = "black", size = 1.2) +
  geom_segment(aes(x = 4.135936, xend = 4.135936, 
                   y = -Inf, yend = 1.102104/2),
               linetype = "dashed", colour = "black", size = 1.2) +
  
  
  annotate("text", x = 5, y = 1.3, label = "µ max", size = 4.5, fontface = "bold") +
  annotate("text", x = 2.1, y = 0.45, label = "Salt tolerance (c)", size = 4.5, fontface = "bold")

p.s1

p.s2 <- ggplot(df.r.p[df.r.p$population.number==3,], aes(x = salt.lvl, y = r.exp, colour = population)) +
  scale_colour_manual(values = c("forestgreen","darkorange1", "magenta2")) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = df.S.jags3, aes(x = salt, y= mean), colour = "forestgreen", size = 1) +
  geom_line(data = df.S.jags27, aes(x = salt, y= mean), colour = "magenta2", size = 1) +
  geom_line(data = df.S.jags32, aes(x = salt, y= mean), colour ="darkorange1", size = 1) +
  
  labs(x = "Salt concentration (g/L)", 
       y = "Exponential growth rate (µ)",
       title = "H - Salt") +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.s2

p.s3 <- ggplot(df.final[df.final$Pop.fac %in% c("11", "4", "9"),], aes(x = r.max_S, y = S.c.mod, colour = Pop.fac)) +
  geom_point(size=4.5) +
  scale_colour_manual(values = c("forestgreen", "magenta", "darkorange")) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Salt tolerance (c)",
       title = "I - Salt") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.s3

fig.1 <- plot_grid(p.t1, p.t2, p.t3, p.p1, p.p2, p.p3, p.s1, p.s2, p.s3, nrow = 3, align='hv', rel_widths = c(1,1,1))

fig.1

ggsave("figures/16a_fig_1_sample_models.T.0.56.jpeg", fig.1, width = 15, height = 15) # PDF was rendering weird

# Figure 2: PCAs and RDAs -----------------------------------------------------------

# We'll need the full dataset with metabolic and pigment data.

df.pca <- df.final %>% select(T.br.0.56, r.max_T, I.comp, r.max_I, N.comp, r.max_N, P.comp, r.max_P, r.max_S, 
                              S.c.mod, chl.a, chl.b, luthein, mean.N.µg.l, mean.P.µg.l, evol.plt) # Prepare the data: selecting only the relevant columns
df.pca

df.pca <- df.pca[-31,] # For now removing wonky/missing points (I.comp is way too high in row 31)

evol.fil <- df.pca$evol.plt

df.pca <- df.pca %>% select(-evol.plt)

pca.result <- prcomp(df.pca, center = TRUE, scale. = TRUE) # Perform PCA
pca.result

sdev <- pca.result$sdev  

sdev %>%
  { .^2 / sum(.^2) } %>%
  print() 
  
df.pca.res <- data.frame(pca.result$x, evol.fil)  # Add grouping factor
colnames(df.pca.res) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", 
                      "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "Evolution")

PCA <- ggplot(df.pca.res, aes(x = PC1, y = PC2, color = Evolution)) +  # PCA biplot visualization
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("PC1 (", round(pca.result$sdev[1]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pca.result$sdev[2]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       title = "PCA of Thermal Performance and Competitive Abilities") +
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  )

PCA

loadings <- as.data.frame(pca.result$rotation) # Extract PCA loadings (rotation matrix)

loadings$PC1 <- loadings$PC1 * max(abs(df.pca.res$PC1)) # Scale loadings to fit within the PCA plot (can adjust scaling factor)
loadings$PC2 <- loadings$PC2 * max(abs(df.pca.res$PC2))

loadings$variable <- rownames(loadings) # Add variable names for annotation

loadings$metric <- factor(loadings$variable, 
                            levels = c("T.br.0.56", "r.max_T", "I.comp", "r.max_I", "N.comp", "r.max_N", 
                                       "P.comp", "r.max_P", "r.max_S", "S.c.mod", "chl.a", "chl.b", "luthein",
                                       "mean.N.µg.l", "mean.P.µg.l"),
                            labels = c("Thermal~breadth", 
                                       "mu~max~(T)", 
                                       "1/I^\"*\"", 
                                       "mu~max~(I)", 
                                       "1/N^\"*\"", 
                                       "mu~max~(N)",
                                       "1/P^\"*\"", 
                                       "mu~max~(P)", 
                                       "mu~max~(Salt)", 
                                       "Salt~tolerance",
                                       "Chlorophyll~italic(a)",
                                       "Chlorophyll~italic(b)",
                                       "Luthein", 
                                       "N~content",
                                       "P~content"))

adjust.x <- c(0.8, 0.4 , 0.02, 0.7, 0.22, 0.75, 0.25, -0.05, -0.05, 0.98, 0.45, 0.97, 0.6, 0.72, 0.72) # adjustments for each label
adjust.y <- c(0.27, -0.02, -0.05, 0.1, -0.03, 0.20, -0.03, 0.2, 0.15, -0.05, 0.3, 0.05, 0, 0, 0)

loadings$var.x <- loadings$PC1 + adjust.x # Need to manually space out labels here. 
loadings$var.y <- loadings$PC2 + adjust.y

# Create PCA plot with arrows
pca_plot_arrows <- ggplot(df.pca.res, aes(x = PC1, y = PC2, color = Evolution)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("PC 1 (", round(pca.result$sdev[1]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       y = paste("PC 2 (", round(pca.result$sdev[2]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = "")) +
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +  # Use your custom colors
  # Add arrows for variable contributions
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size= 1.1) +
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
  scale_y_continuous(breaks = seq(-5, 7, by = 1)) +
  # Add variable names to the plot
  geom_text(data = loadings, aes(x = var.x, y = var.y, label = metric),
            vjust = 1, hjust = 1,  color = "black", size = 5, parse = T) +
  theme(legend.position= c(0.25, 0.25),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 12, face = "plain"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")) 

pca_plot_arrows

ggsave("figures/17_fig_2a_PCA.jpeg", pca_plot_arrows, width = 12, height = 9)

###### Evolutionary history RDA #######

response_vars <- df.pca
explanatory_vars <- model.matrix(~ evol.fil)[, -1]  # Remove intercept

rda_result_evol <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result_evol)
sum(summary(rda_result_evol)$cont$importance[2, 1:rda_result_evol$CCA$rank]) # with P and N, evolutionary environment explains 20.23616% of the variation

rda_var_explained <- summary(rda_result_evol)$cont$importance["Proportion Explained", 1:2] * 100

rda_sites_evol <- as.data.frame(scores(rda_result_evol, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_evol <- as.data.frame(scores(rda_result_evol, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_evol <- as.data.frame(scores(rda_result_evol, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_evol$Evolution <- evol.fil # Add the evolutionary treatment labels to the site scores

rda_constraints_evol$label <- rownames(rda_constraints_evol) # Assign readable labels

rda_species_evol$metric <- factor(rownames(rda_species_evol), 
                          levels = c("T.br.0.56", "r.max_T", "I.comp", "r.max_I", "N.comp", "r.max_N", 
                                     "P.comp", "r.max_P", "r.max_S", "S.c.mod", "chl.a", "chl.b", "luthein",
                                     "mean.N.µg.l", "mean.P.µg.l"),
                          labels = c("Thermal~breadth", 
                                     "mu~max~(T)", 
                                     "1/I^\"*\"", 
                                     "mu~max~(I)", 
                                     "1/N^\"*\"", 
                                     "mu~max~(N)",
                                     "1/P^\"*\"", 
                                     "mu~max~(P)", 
                                     "mu~max~(Salt)", 
                                     "Salt~tolerance",
                                     "Chlorophyll~italic(a)",
                                     "Chlorophyll~italic(b)",
                                     "Luthein", 
                                     "N~content",
                                     "P~content"))

adjust.x.rda1 <- c(0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) # adjustments for each label
adjust.y.rda1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

rda_species_evol$var.x <- rda_species_evol$RDA1 + adjust.x.rda1 # Need to manually space out labels here. 
rda_species_evol$var.y <- rda_species_evol$RDA2 + adjust.y.rda1

# Create evolutionary environment RDA plot with arrows
rda_evol_plot_arrows <- ggplot(rda_sites_evol, aes(x = RDA1, y = RDA2, color = Evolution)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste0("RDA 1 (", round(rda_var_explained[1], 2), "%)"),
       y = paste0("RDA 2 (", round(rda_var_explained[2], 2), "%)"),
       title = "A") +
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +  # Use your custom colors
  # Add arrows for variable contributions
  geom_segment(data = rda_species_evol, aes(x = 0, y = 0, xend = RDA1*2, yend = RDA2*2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size= 1.1) +
  scale_x_continuous(breaks = seq(-20, 25, by = 5)) +
  scale_y_continuous(breaks = seq(-20, 40, by = 5)) +
  # Add variable names to the plot
  geom_text(data = rda_species_evol, aes(x = var.x, y = var.y, label = metric),
            vjust = 1, hjust = 1, color = "black", size = 5, parse = T) +
  theme(legend.position= c(0.15, 0.75),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 12, face = "plain"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.03))

rda_evol_plot_arrows # N and P content is skewing everything. 

df.pca1 <- df.pca %>% 
  select(-mean.N.µg.l, -mean.P.µg.l)

response_vars <- df.pca1
explanatory_vars <- model.matrix(~ evol.fil)[, -1]  # Remove intercept

rda_result_evol <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result_evol) 
sum(summary(rda_result_evol)$cont$importance[2, 1:rda_result_evol$CCA$rank]) # without P and N, evolutionary environment explains 20.23616% of the variation

rda_sites_evol <- as.data.frame(scores(rda_result_evol, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_evol <- as.data.frame(scores(rda_result_evol, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_evol <- as.data.frame(scores(rda_result_evol, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_evol$Evolution <- evol.fil # Add the evolutionary treatment labels to the site scores

rda_constraints_evol$label <- rownames(rda_constraints_evol) # Assign readable labels for centroids

rda_species_evol$metric <- factor(rownames(rda_species_evol), 
                                  levels = c("T.br.0.56", "r.max_T", "I.comp", "r.max_I", "N.comp", "r.max_N", 
                                             "P.comp", "r.max_P", "r.max_S", "S.c.mod", "chl.a", "chl.b", "luthein"),
                                  labels = c("Thermal~breadth", 
                                             "mu~max~(T)", 
                                             "1/I^\"*\"", 
                                             "mu~max~(I)", 
                                             "1/N^\"*\"", 
                                             "mu~max~(N)",
                                             "1/P^\"*\"", 
                                             "mu~max~(P)", 
                                             "mu~max~(Salt)", 
                                             "Salt~tolerance",
                                             "Chlorophyll~italic(a)",
                                             "Chlorophyll~italic(b)",
                                             "Luthein"))

adjust.x.rda1 <- c(0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) # adjustments for each label
adjust.y.rda1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

rda_species_evol$var.x <- rda_species_evol$RDA1 + adjust.x.rda1 # Need to manually space out labels here. 
rda_species_evol$var.y <- rda_species_evol$RDA2 + adjust.y.rda1

# Create evolutionary environment RDA plot with arrows
rda_evol_plot_arrows_no_NP <- ggplot(rda_sites_evol, aes(x = RDA1, y = RDA2, color = Evolution)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("RDA 1 (", round(summary(rda_result_evol)$cont$importance[2, 1] * 100, 2), "%)", sep = ""),
       y = paste("RDA 2 (", round(summary(rda_result_evol)$cont$importance[2, 2] * 100, 2), "%)", sep = ""),
       title = "B") +
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +  # Use your custom colors
  # Add arrows for variable contributions
  geom_segment(data = rda_species_evol, aes(x = 0, y = 0, xend = RDA1*2, yend = RDA2*2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size= 1.1) +
  scale_x_continuous(breaks = seq(-5, 6, by = 1)) +
  scale_y_continuous(breaks = seq(-3, 5, by = 1)) +
  # Add variable names to the plot
  geom_text(data = rda_species_evol, aes(x = var.x, y = var.y, label = metric),
            vjust = 1, hjust = 1, color = "black", size = 5, parse = T) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.03))

rda_evol_plot_arrows_no_NP

###### Ancestry RDA ######

df.pca2 <- df.final %>% select(T.br.0.56, r.max_T, I.comp, r.max_I, N.comp, r.max_N, P.comp, r.max_P, r.max_S, 
                              S.c.mod, chl.a, chl.b, luthein, mean.N.µg.l, mean.P.µg.l, anc.plt) # Prepare the data: selecting only the relevant columns
df.pca2

df.pca2 <- df.pca2[-31,] # For now removing wonky/missing points (I.comp is way too high in row 31)

anc.fil <- df.pca2$anc.plt

df.pca2 <- df.pca2 %>% select(-anc.plt)

response_vars <- df.pca2

explanatory_vars <- model.matrix(~ anc.fil)[, -1]  # Remove intercept

rda_result_anc <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result_anc)
sum(summary(rda_result_anc)$cont$importance[2, 1:rda_result_anc$CCA$rank]) # with P and N, ancestry explains 12.96396% of the variation

rda_sites_anc <- as.data.frame(scores(rda_result_anc, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_anc <- as.data.frame(scores(rda_result_anc, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_anc <- as.data.frame(scores(rda_result_anc, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_anc$Ancestry <- anc.fil # Add the evolutionary treatment labels to the site scores

rda_constraints_anc$label <- rownames(rda_constraints_anc) # Assign readable labels for centroids

rda_species_anc$metric <- factor(rownames(rda_species_anc), 
                                  levels = c("T.br.0.56", "r.max_T", "I.comp", "r.max_I", "N.comp", "r.max_N", 
                                             "P.comp", "r.max_P", "r.max_S", "S.c.mod", "chl.a", "chl.b", "luthein",
                                             "mean.N.µg.l", "mean.P.µg.l"),
                                  labels = c("Thermal~breadth", 
                                             "mu~max~(T)", 
                                             "1/I^\"*\"", 
                                             "mu~max~(I)", 
                                             "1/N^\"*\"", 
                                             "mu~max~(N)",
                                             "1/P^\"*\"", 
                                             "mu~max~(P)", 
                                             "mu~max~(Salt)", 
                                             "Salt~tolerance",
                                             "Chlorophyll~italic(a)",
                                             "Chlorophyll~italic(b)",
                                             "Luthein", 
                                             "N~content",
                                             "P~content"))

adjust.x.rda2 <- c(0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) # adjustments for each label
adjust.y.rda2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

rda_species_anc$var.x <- rda_species_anc$RDA1 + adjust.x.rda2 # Need to manually space out labels here. 
rda_species_anc$var.y <- rda_species_anc$RDA2 + adjust.y.rda2

# Create evolutionary environment RDA plot with arrows
rda_anc_plot_arrows <- ggplot(rda_sites_anc, aes(x = RDA1, y = RDA2, color = Ancestry)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("RDA 1 (", round(summary(rda_result_anc)$cont$importance[2, 1] * 100, 2), "%)", sep = ""),
       y = paste("RDA 2 (", round(summary(rda_result_anc)$cont$importance[2, 2] * 100, 2), "%)", sep = ""),
       title = "C") +
  scale_color_manual(
    name = "Ancestry",  # Update the legend title
    values = c("Population 2" = "darkorange",
               "Population 3" = "deepskyblue1",
               "Population 4" = "forestgreen",
               "Population 5" = "gold",
               "Mixed population" = "magenta3")
  ) +  # Use custom colors
  # Add arrows for variable contributions
  geom_segment(data = rda_species_anc, aes(x = 0, y = 0, xend = RDA1*2, yend = RDA2*2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size= 1.1) +
  scale_x_continuous(breaks = seq(-35, 25, by = 5)) +
  scale_y_continuous(breaks = seq(-35, 20, by = 5)) +
  # Add variable names to the plot
  geom_text(data = rda_species_anc, aes(x = var.x, y = var.y, label = metric),
            vjust = 1, hjust = 1, color = "black", size = 5, parse = T) +
  theme(legend.position= c(0.15, 0.75),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 12, face = "plain"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.03))

rda_anc_plot_arrows # Again, N and P are skewing everything

df.pca3 <- df.pca2 %>% 
  select(-mean.N.µg.l, -mean.P.µg.l)

df.pca3

response_vars <- df.pca3

explanatory_vars <- model.matrix(~ anc.fil)[, -1]  # Remove intercept

rda_result_anc <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result_anc)
sum(summary(rda_result_anc)$cont$importance[2, 1:rda_result_anc$CCA$rank]) # without P and N, ancestry explains 3.230246% of the variation

rda_sites_anc <- as.data.frame(scores(rda_result_anc, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_anc <- as.data.frame(scores(rda_result_anc, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_anc <- as.data.frame(scores(rda_result_anc, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_anc$Ancestry <- anc.fil # Add the evolutionary treatment labels to the site scores

rda_constraints_anc$label <- rownames(rda_constraints_anc) # Assign readable labels for centroids

rda_species_anc$metric <- factor(rownames(rda_species_anc), 
                                 levels = c("T.br.0.56", "r.max_T", "I.comp", "r.max_I", "N.comp", "r.max_N", 
                                            "P.comp", "r.max_P", "r.max_S", "S.c.mod", "chl.a", "chl.b", "luthein"),
                                 labels = c("Thermal~breadth", 
                                            "mu~max~(T)", 
                                            "1/I^\"*\"", 
                                            "mu~max~(I)", 
                                            "1/N^\"*\"", 
                                            "mu~max~(N)",
                                            "1/P^\"*\"", 
                                            "mu~max~(P)", 
                                            "mu~max~(Salt)", 
                                            "Salt~tolerance",
                                            "Chlorophyll~italic(a)",
                                            "Chlorophyll~italic(b)",
                                            "Luthein"))

adjust.x.rda2 <- c(0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) # adjustments for each label
adjust.y.rda2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

rda_species_anc$var.x <- rda_species_anc$RDA1 + adjust.x.rda2 # Need to manually space out labels here. 
rda_species_anc$var.y <- rda_species_anc$RDA2 + adjust.y.rda2

# Create evolutionary environment RDA plot with arrows
rda_anc_plot_arrows_no_NP <- ggplot(rda_sites_anc, aes(x = RDA1, y = RDA2, color = Ancestry)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("RDA 1 (", round(summary(rda_result_anc)$cont$importance[2, 1] * 100, 2), "%)", sep = ""),
       y = paste("RDA 2 (", round(summary(rda_result_anc)$cont$importance[2, 2] * 100, 2), "%)", sep = ""),
       title = "D") +
  scale_color_manual(
    name = "Ancestry",  # Update the legend title
    values = c("Population 2" = "darkorange",
               "Population 3" = "deepskyblue1",
               "Population 4" = "forestgreen",
               "Population 5" = "gold",
               "Mixed population" = "magenta3")
  ) +  # Use custom colors
  # Add arrows for variable contributions
  geom_segment(data = rda_species_anc, aes(x = 0, y = 0, xend = RDA1*2, yend = RDA2*2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size= 1.1) +
  scale_x_continuous(breaks = seq(-14, 12, by = 2)) +
  scale_y_continuous(breaks = seq(-8, 14, by = 2)) +
  # Add variable names to the plot
  geom_text(data = rda_species_anc, aes(x = var.x, y = var.y, label = metric),
            vjust = 1, hjust = 1, color = "black", size = 5, parse = T) +
  theme(legend.position = "none",
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.03))

rda_anc_plot_arrows_no_NP 

df.evol.dummy.leg <- data.frame( # dummy df for making an extractable legend
  RDA1 = rep(0, 8),
  RDA2 = rep(0, 8),
  Evolution = factor(c("Biotic depletion", "Biotic depletion x Salt", "Control",
                       "Light limitation", "Nitrogen limitation", "Ancestral",
                       "Phosphorous limitation", "Salt stress"))
)

legend_evol <- ggplot(df.evol.dummy.leg, aes(x = RDA1, y = RDA2, color = Evolution)) +
  geom_point() +
  scale_color_manual(
    name = "Evolution environment",
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 12, face = "plain"))

evol_legend <- get_legend(legend_evol)

df.anc.dummy.leg <- data.frame( # dummy df for making an extractable legend
  RDA1 = rep(0, 5),
  RDA2 = rep(0, 5),
  Ancestry = factor(c("Population 2", "Population 3", "Population 4",
                       "Population 5", "Mixed population"))
)

legend_anc <- ggplot(df.anc.dummy.leg, aes(x = RDA1, y = RDA2, color = Ancestry)) +
  geom_point() +
  scale_color_manual(
    name = "Ancestry",
    values = c("Population 2" = "darkorange",
               "Population 3" = "deepskyblue1",
               "Population 4" = "forestgreen",
               "Population 5" = "gold",
               "Mixed population" = "magenta3")
  ) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 12, face = "plain"))

anc_legend <- get_legend(legend_anc)

rda_evol_plot_arrows <- rda_evol_plot_arrows + theme(legend.position = "none")
rda_anc_plot_arrows <- rda_anc_plot_arrows + theme(legend.position = "none")

rdas <- plot_grid(rda_evol_plot_arrows, rda_evol_plot_arrows_no_NP, evol_legend, rda_anc_plot_arrows, rda_anc_plot_arrows_no_NP, anc_legend,
                  align= 'hv',
                  nrow = 2)

ggsave("figures/17_fig_2b_supp_rdas.jpeg", rdas, width = 16, height = 8) # PDF was rendering weird

# Figure 3: Intra-gradient trade-offs -------------------------------------

# We'll need our full dataset for the start

###### Temperature ######

df.final$evol.bin <- ifelse(df.final$evol == "none", "ancestral", "evolved") # For binning into evolutionary treatments

df.filt <- df.final %>% 
  mutate(
    z.y = T.br.0.56,
    z.x = r.max_T
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.4)), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

T.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +

  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       color = "Evolutionary History",
       title = "E — Temperature") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("evolved" = 1,  # open circle
               "ancestral" = 5)  # diamond
  ) +
  
  ylim(24.5,29.5) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.scam  # Display the plot

# Quadrant-based statistical testing

x.thresh <- mean(df.filt$z.x)
y.thresh <- mean(df.filt$z.y)

obs.cnt <- df.filt %>%
  filter(z.x > x.thresh, z.y > y.thresh) %>%
  nrow()

null.counts <- replicate(1000, {
  shuffled.df <- df.filt %>%
    mutate(z.x = sample(z.x, rep= F),
           z.y = sample(z.y, rep = F))
  
  sum(shuffled.df$z.x > x.thresh & shuffled.df$z.y > y.thresh)
})

mean(null.counts >= obs.cnt) # p-value 0.445

# Empty space testing (Adapted from Li et al., 2019), using only the convex hull PF points — there are not consistently enough points to fit scams to null data.

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)
  
a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

# We are going to write this out as a for-loop and save results in a df

null.df <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt %>%
    mutate(z.x = sample(z.x, replace = FALSE),
           z.y = sample(z.y, replace = FALSE))
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y)[, 1],
      z.x2 = scale(z.x)[, 1]
    ) # Scale for the point addition algorithm. 
  
  x.ref <- min(shuffled.df$z.x2, na.rm = TRUE) # Min x
  y.ref <- min(shuffled.df$z.y2, na.rm = TRUE) # Min y
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 
  
  shuffled.df <- shuffled.df %>% 
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y)[, 1],
      z.x2 = scale(z.x)[, 1]
    ) # re-scale for the point addition algorithm.
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x", yvar = "z.y") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y)
  
  poly.n <- par.res.n[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x = x.max.n, z.y = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x, poly.n$z.y) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.655?

###### Light ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "L", 'light', 'other')) # for testing regressions.

df.filt <- df.final %>% 
  mutate(
    z.y = I.comp,
    z.x = r.max_I
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.4)), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

I.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +

  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "A — Light") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "light" = 16)  # filled circle
  ) +
  
  ylim(0.025, 0.225) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

I.scam  # Display the plot

# Statistical testing

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin == "light") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff == "light") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.422

# Quadrant-based statistical testing

x.thresh <- mean(df.filt$z.x)
y.thresh <- mean(df.filt$z.y)

obs.cnt <- df.filt %>%
  filter(z.x > x.thresh, z.y > y.thresh) %>%
  nrow()

null.counts <- replicate(1000, {
  shuffled.df <- df.filt %>%
    mutate(z.x = sample(z.x, rep= F),
           z.y = sample(z.y, rep = F))
  
  sum(shuffled.df$z.x > x.thresh & shuffled.df$z.y > y.thresh)
})

mean(null.counts >= obs.cnt) # p-value 0.569

# Li et al 2019 empty space testing

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

# We are going to write this out as a for-loop and save results in a df

null.df <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt %>%
    mutate(z.x = sample(z.x, replace = FALSE),
           z.y = sample(z.y, replace = FALSE))
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y)[, 1],
      z.x2 = scale(z.x)[, 1]
    ) # Scale for the point addition algorithm. 
  
  x.ref <- min(shuffled.df$z.x2, na.rm = TRUE) # Min x
  y.ref <- min(shuffled.df$z.y2, na.rm = TRUE) # Min y
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 
  
  shuffled.df <- shuffled.df %>% 
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y)[, 1],
      z.x2 = scale(z.x)[, 1]
    ) # re-scale for the point addition algorithm.
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x", yvar = "z.y") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y)
  
  poly.n <- par.res.n[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x = x.max.n, z.y = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x, poly.n$z.y) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.488

###### Nitrogen ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "N", 'nitrogen', 'other')) # for testing regressions.

df.filt <- df.final %>% 
  mutate(
    z.y = N.comp,
    z.x = r.max_N
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.4)), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

N.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "B — Nitrogen") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "nitrogen" = 16)  # filled circle
  ) +
  
  ylim(0,1.25) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

N.scam  # Display the plot

# Statistical testing

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin == "nitrogen") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff == "nitrogen") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.549

# Quadrant-based statistical testing

x.thresh <- mean(df.filt$z.x)
y.thresh <- mean(df.filt$z.y)

obs.cnt <- df.filt %>%
  filter(z.x > x.thresh, z.y > y.thresh) %>%
  nrow()

null.counts <- replicate(1000, {
  shuffled.df <- df.filt %>%
    mutate(z.x = sample(z.x, rep= F),
           z.y = sample(z.y, rep = F))
  
  sum(shuffled.df$z.x > x.thresh & shuffled.df$z.y > y.thresh)
})

mean(null.counts >= obs.cnt) # p-value 0.906

# Li et al 2019 empty space testing

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

# We are going to write this out as a for-loop and save results in a df

null.df <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt %>%
    mutate(z.x = sample(z.x, replace = FALSE),
           z.y = sample(z.y, replace = FALSE))
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y)[, 1],
      z.x2 = scale(z.x)[, 1]
    ) # Scale for the point addition algorithm. 
  
  x.ref <- min(shuffled.df$z.x2, na.rm = TRUE) # Min x
  y.ref <- min(shuffled.df$z.y2, na.rm = TRUE) # Min y
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 
  
  shuffled.df <- shuffled.df %>% 
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y)[, 1],
      z.x2 = scale(z.x)[, 1]
    ) # re-scale for the point addition algorithm.
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x", yvar = "z.y") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y)
  
  poly.n <- par.res.n[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x = x.max.n, z.y = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x, poly.n$z.y) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.33

###### Phosphorous ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "P", 'phosphorous', 'other')) # for testing regressions.

df.filt <- df.final %>% 
  mutate(
    z.y = P.comp,
    z.x = r.max_P
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.4)), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

P.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       color = "Evolutionary History",
       title = "C — Phosphorous") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "phosphorous" = 16)  # filled circle
  ) +
  
  ylim(0, 2.7) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

P.scam  # Display the plot

# Statistical testing

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin == "phosphorous") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff == "phosphorous") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.15

# Quadrant-based statistical testing

x.thresh <- mean(df.filt$z.x)
y.thresh <- mean(df.filt$z.y)

obs.cnt <- df.filt %>%
  filter(z.x > x.thresh, z.y > y.thresh) %>%
  nrow()

null.counts <- replicate(1000, {
  shuffled.df <- df.filt %>%
    mutate(z.x = sample(z.x, rep= F),
           z.y = sample(z.y, rep = F))
  
  sum(shuffled.df$z.x > x.thresh & shuffled.df$z.y > y.thresh)
})

mean(null.counts >= obs.cnt) # p-value 0.979

# Li et al 2019 empty space testing

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

# We are going to write this out as a for-loop and save results in a df

null.df <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt %>%
    mutate(z.x = sample(z.x, replace = FALSE),
           z.y = sample(z.y, replace = FALSE))
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y)[, 1],
      z.x2 = scale(z.x)[, 1]
    ) # Scale for the point addition algorithm. 
  
  x.ref <- min(shuffled.df$z.x2, na.rm = TRUE) # Min x
  y.ref <- min(shuffled.df$z.y2, na.rm = TRUE) # Min y
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 
  
  shuffled.df <- shuffled.df %>% 
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y)[, 1],
      z.x2 = scale(z.x)[, 1]
    ) # re-scale for the point addition algorithm.
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x", yvar = "z.y") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y)
  
  poly.n <- par.res.n[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x = x.max.n, z.y = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x, poly.n$z.y) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.015

# Let's create the descriptive figure for the supplement - how I detect Pareto points

df.filt <- df.filt %>%
  mutate(par.stat = if_else(Pop.fac %in% par.res.1$Pop.fac, "Y", "N"))

p.par <- ggplot(df.filt, aes(x = z.x, y = z.y, colour = par.stat)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = par.res.1, aes(x = z.x, y = z.y), colour = "black", size = 1) +  # Pareto frontier line
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       color = "Evolutionary History",
       title = "A — Raw Pareto frontier") +  # labels
  
  scale_colour_manual(
    name = "Pareto frontier",  # Update the legend title
    values = c("Y" = "magenta2",
               "N" = "black")
  ) +
  
  ylim(0, 2.7) +
  xlim(1.22, 1.64) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )


p.par

df.filt <- df.filt %>%
  mutate(par.stat = case_when(
    Pop.fac %in% par.res.1$Pop.fac & !(Pop.fac %in% par.res.2$Pop.fac) ~ "removed",
    Pop.fac %in% par.res.1$Pop.fac ~ "initial",
    Pop.fac %in% par.res.2$Pop.fac ~ "added",
    TRUE ~ "N"
  ))

p.par2 <- ggplot(df.filt, aes(x = z.x, y = z.y, colour = par.stat)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = par.res.1, aes(x = z.x, y = z.y), colour = "black", size = 1) +  # Pareto frontier line
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       color = "Evolutionary History",
       title = "B — Modified Pareto front dataset") +  # labels
  
  scale_colour_manual(
    name = "Pareto frontier",  # Update the legend title
    values = c("initial" = "magenta2",
               "removed" = "darkorange",
               "added" = "forestgreen",
               "N" = "black")
  ) +
  
  ylim(0, 2.7) +
  xlim(1.22, 1.64) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

p.par2

df.filt <- df.filt %>%
  filter(!par.stat %in% c("N", "removed"))

p.par3 <- ggplot(df.filt, aes(x = z.x, y = z.y, colour = par.stat)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       color = "Evolutionary History",
       title = "C — SCAM fit to Pareto front points") +  # labels
  
  scale_colour_manual(
    name = "Pareto frontier",  # Update the legend title
    values = c("initial" = "magenta2",
               "added" = "forestgreen")
  ) +
  
  ylim(0, 2.7) +
  xlim(1.22, 1.64) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

p.par3

plots.par <- plot_grid(p.par, p.par2, p.par3, nrow = 1)

ggsave("figures/21_supp_fig_1_pareto_front_construction.jpeg", plots.par, width = 15, height = 5)

###### Salt ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                   ifelse(df.final$evol %in% c("S", "BS"), 'salt', 'other')) # Ss and BSs are now equivalent 

df.filt <- df.final %>% 
  mutate(
    z.y = S.c.mod,
    z.x = r.max_S
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = par.res.4), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

S.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +

  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Salt tolerance (c)", 
       color = "Evolutionary History",
       title = "D — Salt") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "salt" = 16)  # filled circle
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

S.scam  # Display the plot

# Statistical testing

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin == "salt") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff == "salt") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0!

# Quadrant-based statistical testing

x.thresh <- mean(df.filt$z.x)
y.thresh <- mean(df.filt$z.y)

obs.cnt <- df.filt %>%
  filter(z.x > x.thresh, z.y > y.thresh) %>%
  nrow()

null.counts <- replicate(1000, {
  shuffled.df <- df.filt %>%
    mutate(z.x = sample(z.x, rep= F),
           z.y = sample(z.y, rep = F))
  
  sum(shuffled.df$z.x > x.thresh & shuffled.df$z.y > y.thresh)
})

mean(null.counts >= obs.cnt) # p-value 0.812

# Li et al 2019 empty space testing

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

# We are going to write this out as a for-loop and save results in a df

null.df <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt %>%
    mutate(z.x = sample(z.x, replace = FALSE),
           z.y = sample(z.y, replace = FALSE))
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y)[, 1],
      z.x2 = scale(z.x)[, 1]
    ) # Scale for the point addition algorithm. 
  
  x.ref <- min(shuffled.df$z.x2, na.rm = TRUE) # Min x
  y.ref <- min(shuffled.df$z.y2, na.rm = TRUE) # Min y
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 
  
  shuffled.df <- shuffled.df %>% 
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y)[, 1],
      z.x2 = scale(z.x)[, 1]
    ) # re-scale for the point addition algorithm.
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x", yvar = "z.y") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y)
  
  poly.n <- par.res.n[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x = x.max.n, z.y = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x, poly.n$z.y) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.088

# Figure 1 : all plots together

plots <- list(I.scam, N.scam, P.scam, S.scam, T.scam)

legend_df <- data.frame(
  x = c(1, 2, 1, 2, 1, 2, 1, 2),
  y = c(1, 1, 2, 2, 1, 1, 2, 2),
  Group = factor(c("Ancestral", "Other", "Matching", "Matching", "Ancestral", "Other", "Matching", "Matching")),
  Group2 = factor(c("Biotic depletion", "Biotic depletion x Salt", "Control", "Light limitation", "Nitrogen limitation", "Ancestral", "Phosphorous limitation", "Salt stress")),
  LineType = factor(c("Outer", "Inner (75%)", "Outer", "Inner (75%)", "Outer", "Inner (75%)", "Outer", "Inner (75%)"))
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y)) +
  geom_point(aes(shape = Group, colour = Group2), size = 3, stroke = 1.5) +
  geom_line(aes(linetype = LineType), size = 1) +
  scale_shape_manual(name = NULL,
                    values = c("Ancestral" = 5, 
                    "Other" = 1, 
                    "Matching" = 16)) +
  
  scale_color_manual(
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_linetype_manual(values = c("Outer" = "dashed", "Inner (75%)" = "solid"),
                        labels = c("Outer", "Inner (75%)"),
                        name = "Pareto front") +
  labs(linetype = "Quantile regression", color = "Evolutionary context") +
  theme_void() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "lines")
  )


legend_only <- get_legend(legend_plot)

all_plots <- c(plots, list(legend_only))

grad_toffs <- plot_grid(plotlist = all_plots,
          ncol = 2,
          align = "hv")

ggsave("figures/18_fig_3_intra-gradient_tradeoffs.jpeg", grad_toffs, width = 8, height = 12)

# Then we will bring in inter-specific datasets and plot the position of their metrics on our plots

# Figure 4: Intra-gradient, interspecific trade-offs-------------------------------------

# Bringing in interspecific data sets

df.thomas <- read.csv("data-processed/17_Thomas2012_TPCs.csv")

df.thomas <- df.thomas %>% 
  mutate(source = 'Thomas2012')

df.bestion.t <- read.csv("data-processed/18a_Bestion2018_TPCs.csv")

df.bestion.t <- df.bestion.t %>%
  mutate(source = 'Bestion2018') 

df.bestion.p <- read.csv("data-processed/18c_Bestion2018_P_Monodss.csv")

df.bestion.p <- df.bestion.p %>% 
  mutate(source = 'Bestion2018') %>%
  select(-Temp)

df.lewington.t <- read.csv("data-processed/19b_Lewington2019_TPCs.csv")

df.lewington.t <- df.lewington.t %>% 
  mutate(source = 'Lewington2019') 

df.lewington.n <- read.csv("data-processed/19e_Lewington2019_Nit_Monods.csv")

df.lewington.n <- df.lewington.n %>% 
  mutate(source = 'Lewington2019') 

df.lewington.l <- read.csv("data-processed/19h_Lewington2019_Light_Monods.csv")

df.lewington.l <- df.lewington.l %>% 
  mutate(source = 'Lewington2019') 

df.narwani <- read.csv("data-processed/20_Narwani2015_summmary.csv") # Already includes the relevant statistics (no raw data) so I did not re-estimate anything. 

df.narwani <- df.narwani %>% 
  mutate(source = 'Narwani2015') 

#df.kontopoulos <- read.csv("data-processed/21a_Kontopoulos2020_TPCs.csv") These numbers are really low, we'll avoid them for now.
#head(df.kontopoulos) # temperature, low r.maxes

df.edwards.2015 <- read.csv("data-processed/26a_Edwards_2015_growth_data.csv")

df.edwards.2015 <- df.edwards.2015 %>% 
  mutate(source = 'Edwards2015') 

df.edwards.2016.t <- read.csv("data-processed/27b_Edwards2016_TPCs.csv")

df.edwards.2016.t <- df.edwards.2016.t %>% 
  mutate(source = 'Edwards2016') 

df.edwards.2016.l <- read.csv("data-processed/27d_Edwards2016_light.csv")

df.edwards.2016.l <- df.edwards.2016.l %>% 
  mutate(source = 'Edwards2016') 

###### Temperature ######

# We have thomas, bestion, lewington and edwards 2016

df.int <- rbind(df.thomas, df.bestion.t, df.lewington.t, df.edwards.2016.t)

df.filt <- df.int %>%
  drop_na(T.max.raw, T.min.raw, r.max.raw) %>%  # remove rows with any NA in these columns
  mutate(
    z.y = T.max.raw - T.min.raw,
    z.x = r.max.raw
  )

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

# OK now generate the same Pareto front for our data!

df.filt2 <- df.final %>% 
  mutate(
    z.y = T.br.0,
    z.x = r.max_T
  ) # Specify the x and y variables

df.filt2 <- df.filt2 %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt2$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt2$z.y2, na.rm = TRUE) # Min y

df.filt2 <- df.filt2 %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt2 <- df.filt2 %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt2 <- df.filt2 %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt2$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt2$z.y2, na.rm = TRUE) # Min y

df.filt2 <- df.filt2 %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.3 <- par_frt(df.filt2[df.filt2$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt2 %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.4)), data = par.res.4) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt2$z.x), max(df.filt2$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.2 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

T.int.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, colour=source)) + # Thomas 2012
  geom_point(size = 2) +
  
  geom_point(data= df.filt2, aes(x = z.x, y = z.y), colour = "forestgreen", size = 2) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "forestgreen", size = 1.1, inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       title = "D — Temperature") +  # labels
  
  scale_color_manual(
    name = "Source of data",  # Update the legend title
    values = c("Edwards2016" = "magenta3",
               "Lewington2019" = "deepskyblue1",
               "Thomas2012" = "gold",
               "Bestion2018" = "darkorange")
  ) +
  
  ylim(5,40) +
  
  theme_classic() +
  theme(
    legend.position = c('none'),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

T.int.scam  # Display the plot

min(df.filt2$z.x)
mean(df.filt$z.x < min(df.filt2$z.x)) * 100

###### Light ######

# Narwani, Lewington, and Edwards 2016

df.narwani.l <- df.narwani %>%
  mutate(
    X = NA,
    Sp.id = Species.name.phylo,
    DIC = NA,               
    K.s = NA,
    r.max = umax.light,
    R.jag = Istar,
    R.mth = Istar,
    source = 'Narwani2015'
  ) %>%
  select(X, Sp.id, DIC, K.s, r.max, R.jag, R.mth, source)

df.int <- rbind(df.narwani.l, df.lewington.l, df.edwards.2016.l)

df.filt <- df.int %>%
  drop_na(r.max, R.mth) %>%  # remove rows with any NA in these columns
  mutate(
    z.y = 1/R.mth,
    z.x = r.max
  )

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

# OK now generate the same Pareto front for our data!

df.filt2 <- df.final %>% 
  mutate(
    z.y = I.comp.10,
    z.x = r.max_I
  ) # Specify the x and y variables

df.filt2 <- df.filt2 %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt2$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt2$z.y2, na.rm = TRUE) # Min y

df.filt2 <- df.filt2 %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt2 <- df.filt2 %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt2 <- df.filt2 %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt2$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt2$z.y2, na.rm = TRUE) # Min y

df.filt2 <- df.filt2 %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.3 <- par_frt(df.filt2[df.filt2$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt2 %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.4)), data = par.res.4) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt2$z.x), max(df.filt2$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.2 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

L.int.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, colour=source)) + 
  geom_point(size = 2) +
  
  geom_point(data= df.filt2, aes(x = z.x, y = z.y), colour = "forestgreen", size = 2) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "forestgreen", size = 1.1, inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  scale_color_manual(
    name = "Source of data",  # Update the legend title
    values = c("Edwards2016" = "magenta3",
               "Lewington2019" = "deepskyblue1",
               "Narwani2015" = "darkred")
  ) +
  
  theme_classic() +
  theme(
    legend.position = c('none'),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

L.int.scam  # Display the plot

min(df.filt2$z.x)
mean(df.filt$z.x < min(df.filt2$z.x)) * 100

###### Nitrogen ######

# So we have Edwards 2015, Narwani and Lewington

df.edwards.2015.n <- df.edwards.2015 %>%
  mutate(
    X = NA,
    Sp.id = species,
    DIC = NA,               
    K.s = k_nit_m,
    r.max = coalesce(mu_nit, mu_inf_nit),  # Prefer mu_nit, fall back on mu_inf_nit
    R.jag = 0.1*K.s/(r.max-0.1),
    R.mth = 0.1*K.s/(r.max-0.1),
    source = 'Edwards2015'
  ) %>%
  select(X, Sp.id, DIC, K.s, r.max, R.jag, R.mth, source)

df.narwani.n <- df.narwani %>%
  mutate(
    X = NA,
    Sp.id = Species.name.phylo,
    DIC = NA,               
    K.s = NA,
    r.max = umax.nitrate,
    R.jag = Nstar,
    R.mth = Nstar,
    source = 'Narwani2015'
  ) %>%
  select(X, Sp.id, DIC, K.s, r.max, R.jag, R.mth, source)

df.int <- rbind(df.lewington.n, df.edwards.2015.n)

df.filt <- df.int %>%
  drop_na(r.max, R.mth) %>%  # remove rows with any NA in these columns
  mutate(
    z.y = 1/R.mth,
    z.x = r.max
  )

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts, df.filt[c(23,25),]) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

# OK now generate the same Pareto front for our data!

df.filt2 <- df.final %>% 
  mutate(
    z.y = N.comp.10,
    z.x = r.max_N
  ) # Specify the x and y variables

df.filt2 <- df.filt2 %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt2$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt2$z.y2, na.rm = TRUE) # Min y

df.filt2 <- df.filt2 %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt2 <- df.filt2 %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt2 <- df.filt2 %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt2$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt2$z.y2, na.rm = TRUE) # Min y

df.filt2 <- df.filt2 %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.3 <- par_frt(df.filt2[df.filt2$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt2 %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.4)), data = par.res.4) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt2$z.x), max(df.filt2$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.2 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

N.int.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, colour=source)) + 
  geom_point(size = 2) +
  
  geom_point(data= df.filt2, aes(x = z.x, y = z.y), colour = "forestgreen", size = 2) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "forestgreen", size = 1.1, inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/N*)", 
       title = "B — Nitrogen") +  # labels
  
  scale_color_manual(
    name = "Source of data",  # Update the legend title
    values = c("Edwards2015" = "blue",
               "Lewington2019" = "deepskyblue1")
  ) +
  
  ylim(0,30) +
  
  theme_classic() +
  theme(
    legend.position = c('none'),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

N.int.scam  # Display the plot

min(df.filt2$z.x)
mean(df.filt$z.x < min(df.filt2$z.x)) * 100

###### Phosphorous ######

# Bestion, Edwards 2015 and Narwani

df.edwards.2015.p <- df.edwards.2015 %>%
  mutate(
    X = NA,
    Sp.id = species,
    DIC = NA,               
    K.s = k_p_m,
    r.max = coalesce(mu_p, mu_inf_p),  # Prefer mu_nit, fall back on mu_inf_nit
    R.jag = 0.1*K.s/(r.max-0.1),
    R.mth = 0.1*K.s/(r.max-0.1),
    source = 'Edwards2015'
  ) %>%
  select(X, Sp.id, DIC, K.s, r.max, R.jag, R.mth, source)

df.narwani.p <- df.narwani %>%
  mutate(
    X = NA,
    Sp.id = Species.name.phylo,
    DIC = NA,               
    K.s = NA,
    r.max = umax.phosphate,
    R.jag = Pstar,
    R.mth = Pstar,
    source = 'Narwani2015'
  ) %>%
  select(X, Sp.id, DIC, K.s, r.max, R.jag, R.mth, source)

df.int <- rbind(df.lewington.n, df.edwards.2015.n, df.bestion.p, df.narwani.p)

df.filt <- df.int %>%
  drop_na(r.max, R.mth) %>%  # remove rows with any NA in these columns
  mutate(
    z.y = 1/R.mth,
    z.x = r.max
  )

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt[df.filt$dist.sc < 2.1,]$z.x), max(df.filt[df.filt$dist.sc < 2.1,]$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

# OK now generate the same Pareto front for our data!

df.filt2 <- df.final %>% 
  mutate(
    z.y = P.comp.10,
    z.x = r.max_P
  ) # Specify the x and y variables

df.filt2 <- df.filt2 %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt2$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt2$z.y2, na.rm = TRUE) # Min y

df.filt2 <- df.filt2 %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt2 <- df.filt2 %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt2 <- df.filt2 %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt2$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt2$z.y2, na.rm = TRUE) # Min y

df.filt2 <- df.filt2 %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.3 <- par_frt(df.filt2[df.filt2$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt2 %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit <- lm(z.y ~ z.x, data = par.res.4) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt2$z.x), max(df.filt2$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.2 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

P.int.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, colour=source)) + 
  geom_point(size = 2) +
  
  geom_point(data= df.filt2, aes(x = z.x, y = z.y), colour = "forestgreen", size = 2) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "forestgreen", size = 1.1, inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       title = "C — Phosphorous") +  # labels
  
  scale_color_manual(
    name = "Source of data",  # Update the legend title
    values = c("Edwards2015" = "blue",
               "Lewington2019" = "deepskyblue1",
               "Bestion2018" = "darkorange",
               "Narwani2015" = "darkred")
  ) +
  
  ylim(0,30) +
  
  theme_classic() +
  theme(
    legend.position = c('none'),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

P.int.scam  # Display the plot

min(df.filt2$z.x)
mean(df.filt$z.x < min(df.filt2$z.x)) * 100

###### Make the figures ######

plots2 <- list(L.int.scam, N.int.scam, P.int.scam, T.int.scam)

legend_sp <- data.frame(
  x = 1:7,
  y = 1:7,
  Dataset = factor(c(
    "Bestion et al., 2018",
    "Edwards et al., 2015",
    "Edwards et al., 2016",
    "Lewington-Pearce et al., 2019",
    "Narwani et al., 2015",
    "Thomas et al., 2012",
    "Laurich et al., 2025"
  ))
)

legend_sp_plot <- ggplot(legend_sp, aes(x = x, y = y, color = Dataset)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Bestion et al., 2018" = "darkorange",
               "Edwards et al., 2015" = "blue",
               "Edwards et al., 2016" = "magenta3",
               "Lewington-Pearce et al., 2019" = "deepskyblue",
               "Narwani et al., 2015" = "darkred",
               "Thomas et al., 2012" = "gold",
               "Laurich et al., 2025" = "forestgreen"),
    name = "Data set"
  ) +
  theme_void() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )


sp_legend <- get_legend(legend_sp_plot)

all_sp_plots <- c(plots2, list(sp_legend))

grad_sp_toffs <- plot_grid(plotlist = all_sp_plots,
                        ncol = 2,
                        align = "hv")

ggsave("figures/19_fig_4a_sp_tradeoffs.jpeg", grad_sp_toffs, width = 8, height = 12)

# Figure 5: Inter-gradient trade-offs -------------------------------------

# I think for now we will limit this to competitive abilities, salt tolerance and thermal breadth.
# We'll create a half-full grid. 
# Order: light, nitrogen, phosphorous, salt, temperature

###### Light comparisons ######

###### Light v Nitrogen ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "L", 'light', ifelse(df.final$evol == 'N', 'nit', 'other'))) # for testing regressions.

df.filt <- df.final %>% 
  mutate(
    z.y = I.comp,
    z.x = N.comp
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.4)), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

LN.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Competitive ability (1/N*)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "A — Light ~ Nitrogen") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "nit" = 16,
               "light" = 16)  # filled circle
  ) +
  
  ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LN.scam.PF  # Display the plot

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin %in% c("light", "nit")) %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff %in% c("light", "nit")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.253

###### Light v Phosphorous ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "L", 'light', ifelse(df.final$evol == 'P', 'phos', 'other'))) # for testing regressions.

df.filt <- df.final %>% 
  mutate(
    z.y = I.comp,
    z.x = P.comp
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- lm(z.y ~ z.x, data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

LP.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "B — Light ~ Phosphorous") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "phos" = 16,
               "light" = 16)  # filled circle
  ) +
  
  ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LP.scam.PF  # Display the plot

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin %in% c("light", "phos")) %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff %in% c("light", "phos")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.439

###### Light v salt ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "L", 'light', ifelse(df.final$evol %in% c("S", "BS"), 'salt', 'other'))) # for testing regressions.

df.filt <- df.final %>% 
  mutate(
    z.y = I.comp,
    z.x = S.c.mod
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.4)), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

LS.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Salt tolerance (c)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "C — Light ~ Salt") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "salt" = 16,
               "light" = 16)  # filled circle
  ) +
  
  ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LS.scam.PF  # Display the plot

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin %in% c("light", "salt")) %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff %in% c("light", "salt")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.175

###### Light v temperature ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "L", 'light', 'other')) # for testing regressions.

df.filt <- df.final %>% 
  mutate(
    z.y = I.comp,
    z.x = T.br.0.56
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.4)), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

LT.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Thermal breadth (°C)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "D — Light ~ Temperature") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "light" = 16)  # filled circle
  ) +
  
  ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LT.scam.PF  # Display the plot

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin %in% c("light")) %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff %in% c("light")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.353

###### Nitrogen comparisons ######

###### N v P ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "N", 'nit', ifelse(df.final$evol == "P", 'phos', 'other'))) # for testing regressions.

df.filt <- df.final %>% 
  mutate(
    z.y = N.comp,
    z.x = P.comp
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- lm(z.y ~ z.x, data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

NP.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "E — Nitrogen ~ Phosphorous") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "phos" = 16,
               "nit" = 16)  # filled circle
  ) +
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NP.scam.PF  # Display the plot

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin %in% c("nit", "phos")) %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff %in% c("nit", "phos")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.024

###### Nit v Salt ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "N", 'nit', ifelse(df.final$evol %in% c("S", "BS"), 'salt', 'other')))

df.filt <- df.final %>% 
  mutate(
    z.y = N.comp,
    z.x = S.c.mod
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- lm(z.y ~ z.x, data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

NS.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Salt tolerance (c)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "F — Nitrogen ~ Salt") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "salt" = 16,
               "nit" = 16)  # filled circle
  ) +
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NS.scam.PF  # Display the plot

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin %in% c("salt", "nit")) %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff %in% c("salt", "nit")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.22

###### N v Temp ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "N", 'nit', 'other')) # for testing regressions.

df.filt <- df.final %>% 
  mutate(
    z.y = N.comp,
    z.x = T.br.0.56
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- lm(z.y ~ z.x, data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

NT.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Thermal breadth  (°C)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "G — Nitrogen ~ Temperature") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "nit" = 16)  # filled circle
  ) +
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NT.scam.PF  # Display the plot

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin %in% c("nit")) %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff %in% c("nit")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.68

###### Phosphorous comparisons ######

###### P v salt ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "P", 'phos', ifelse(df.final$evol %in% c("S", "BS"), 'salt', 'other')))

df.filt <- df.final %>% 
  mutate(
    z.y = P.comp,
    z.x = S.c.mod
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- lm(z.y ~ z.x, data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

PS.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Salt tolerance (c)",    
       y = "Competitive ability (1/P*)", 
       color = "Evolutionary History",
       title = "H — Phosphorous ~ Salt") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "salt" = 16,
               "phos" = 16)  # filled circle
  ) +
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

PS.scam.PF  # Display the plot

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin %in% c("salt", "phos")) %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff %in% c("salt", "phos")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.005

###### P v T ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "P", 'phos', 'other')) # for testing regressions.

df.filt <- df.final %>% 
  mutate(
    z.y = P.comp,
    z.x = T.br.0.56
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- lm(z.y ~ z.x, data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

PT.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Thermal breadth  (°C)",    
       y = "Competitive ability (1/P*)", 
       color = "Evolutionary History",
       title = "I — Phosphorous ~ Temperature") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "phos" = 16)  # filled circle
  ) +
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

PT.scam.PF  # Display the plot

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin %in% c("phos")) %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff %in% c("phos")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.351

###### Finally, salt v temp! ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol %in% c("S", "BS"), 'salt', 'other')) # for testing regressions.

df.filt <- df.final %>% 
  mutate(
    z.y = S.c.mod,
    z.x = T.br.0.56
  ) # Specify the x and y variables

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # Scale for the point addition algorithm. 

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

df.filt <- df.filt %>% 
  filter(dist.sc < 2.4) # Error trimming

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale for the point addition algorithm.

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.1) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.1$z.x2[i]
  y1 <- par.res.1$z.y2[i]
  x2 <- par.res.1$z.x2[i + 1]
  y2 <- par.res.1$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.1$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.2 <- bind_rows(par.res.1, buff.pts) %>% distinct() # Add these to my existing PF

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = nrow(par.res.2)), data = par.res.2) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

buff.pts <- data.frame() # This will hold our extra data to potentially add

for (i in 1:(nrow(par.res.3) - 1)) { # Loop over each segment (slope)
  x1 <- par.res.3$z.x2[i]
  y1 <- par.res.3$z.y2[i]
  x2 <- par.res.3$z.x2[i + 1]
  y2 <- par.res.3$z.y2[i + 1]
  
  df.cand <- df.filt %>% # Check all other points
    filter(!X %in% par.res.3$X) %>%  # Exclude existing Pareto points
    rowwise() %>%
    mutate(
      dist = point_line_distance(z.x2, z.y2, x1, y1, x2, y2),
      in_x_range = between(z.x2, min(x1, x2) - buffer, max(x1, x2) + buffer),
      in_y_range = between(z.y2, min(y1, y2) - buffer, max(y1, y2) + buffer)
    ) %>%
    filter(dist <= buffer, in_x_range, in_y_range)
  
  # Append
  buff.pts <- bind_rows(buff.pts, df.cand)
}

par.res.4 <- bind_rows(par.res.3, buff.pts) %>% distinct() # Add these to my existing PF

fit2 <- lm(z.y ~ z.x, data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

ST.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Thermal breadth  (°C)",    
       y = "Salt tolerance (c)", 
       color = "Evolutionary History",
       title = "J — Salt ~ Temperature") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "salt" = 16)  # filled circle
  ) +
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

ST.scam.PF  # Display the plot

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin %in% c("salt")) %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n_iter <- 1000
null_counts <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts[i] <- df.shuff %>%
    filter(evol.shuff %in% c("salt")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.067

# The real figure

legend_df.2 <- data.frame(
  x = c(1, 2, 1, 2, 1, 2, 1, 2),
  y = c(1, 1, 2, 2, 1, 1, 2, 2),
  Group = factor(c("Ancestral", "Other", "Matching", "Matching", "Ancestral", "Other", "Matching", "Matching")),
  Group2 = factor(c("Biotic depletion", "Biotic depletion x Salt", "Control", "Light limitation", "Nitrogen limitation", "Ancestral", "Phosphorous limitation", "Salt stress")),
  LineType = factor(c("Pareto Front", "Pareto Front", "Pareto Front", "Pareto Front", "Pareto Front", "Pareto Front", "Pareto Front", "Pareto Front"))
)

legend_plot.2 <- ggplot(legend_df.2, aes(x = x, y = y)) +
  geom_point(aes(shape = Group, colour = Group2), size = 3, stroke = 1.5) +
  geom_line(aes(linetype = LineType), size = 1) +
  
  scale_shape_manual(name = NULL,
                     values = c("Ancestral" = 5, 
                                "Other" = 1, 
                                "Matching" = 16)) +
  
  scale_color_manual(
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_linetype_manual(name = NULL, values = c("Pareto Front" = "solid")) +
  
  labs(color = "Evolutionary context") +
  theme_void() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "lines")
  )

legend_only.2 <- get_legend(legend_plot.2)

full_toffs <- plot_grid(
  LN.scam.PF, LP.scam.PF, LS.scam.PF, LT.scam.PF,
  NULL, NP.scam.PF, NS.scam.PF, NT.scam.PF,
  NULL, NULL, PS.scam.PF, PT.scam.PF,
  legend_only.2, NULL, NULL, ST.scam.PF,
  ncol = 4,
  align = "hv",
  axis = "tblr"
)

ggsave("figures/20_fig_5_full_intergradient_toffs.jpeg", full_toffs, width = 15, height = 15)

###### Done!!! ######
