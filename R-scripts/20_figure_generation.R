# Jason R Laurich
# April 14, 2025

# We are going to generate "final" figures for the project here

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

par_frt_tolerant <- function(df, xvar, yvar, tolerance = 0.05) {
  df <- df[order(-df[[xvar]], df[[yvar]]), ]
  pareto_points <- df[1, ]
  last_y <- df[1, yvar]
  
  for (i in 2:nrow(df)) {
    current_y <- df[i, yvar]
    threshold <- last_y * (1 - tolerance)  # local threshold
    if (current_y >= threshold) {
      pareto_points <- rbind(pareto_points, df[i, ])
      if (current_y > last_y) last_y <- current_y  # update "best" y value
    }
  }
  
  return(pareto_points)
}

find_nearest_index <- function(x, ref_vec) { # Write a function to find the closest point to the actual r.max value in the pred data frame. 
  which.min(abs(ref_vec - x))
} # For significance testing based on location of points relative to a curve

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

# Figure 1: Representation of model fits ----------------------------------

df.final %>% 
  filter(T.br %in% c(min(T.br), max(T.br), median(T.br))) %>% 
  select(Pop.fac, T.br, Topt, r.max_T, S.c.mod, r.max_S, P.comp, r.max_P) %>% 
  print() 
  
# Ok so we are going to work with 3 populations - 17, 27, and cc1690. These represent a good spread
# As R2jags objects, 17 is number 7, 27 is number 18, cc1690 is 38 (for T) or 37 (for salt and P)

# Temperature

for (i in c(7,18,38)){      # Temperature R2jags (model fits)
  load(paste0("R2jags-objects/pop_", i, "_lactin2.RData"))
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:6),]
  df.jags.plot$temp <- seq(0, 45, 0.05)
  assign(paste0("df.T.jags", i), df.jags.plot)
}

df.r.t <- read.csv("data-processed/05_final_r_estimates.csv") # Growth data across temperatures
head(df.r.t) # population is the factor, population.number is the corresponding # (e.g. 7, 18, 38)

df.r.t <- df.r.t %>% 
  filter(population.number %in% c(7, 18, 38)) %>% 
  print()

p.t <- ggplot(df.r.t, aes(x = temperature, y = r.exp, colour = population)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  scale_colour_manual(values = c("darkorange1", "magenta2", "forestgreen")) +
  ylim(-2, 7) +

  geom_line(data = df.T.jags7, aes(x = temp, y= mean), colour = "darkorange1", size = 1) +
  geom_line(data = df.T.jags18, aes(x = temp, y= mean), colour = "magenta2", size = 1) +
  geom_line(data = df.T.jags38, aes(x = temp, y= mean), colour = "forestgreen", size = 1) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "A") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.t

# Phosphorous

for (i in c(7,18,37)){      # Phosphorous
  load(paste0("R2jags-objects/pop_", i, "_phosphorous_monod.RData"))
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:3, 2005),]
  df.jags.plot$phos <- seq(0, 50, 0.025)
  assign(paste0("df.P.jags", i), df.jags.plot)
}

df.r.p <- read.csv("data-processed/12a_phosphorous_r_estimates.csv") # Growth data across P levels
head(df.r.p) # population is the factor, population.number is the corresponding # (e.g. 7, 18, 37)

df.r.p <- df.r.p %>% 
  filter(population.number %in% c(7, 18, 37)) %>% 
  print()

p.p <- ggplot(df.r.p, aes(x = phos.lvl, y = r.exp, colour = population)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  scale_colour_manual(values = c("darkorange1", "magenta2", "forestgreen")) +
  ylim(-0.5, 2.5) +
  
  geom_line(data = df.P.jags7, aes(x = phos, y= mean), colour = "darkorange1", size = 1) +
  geom_line(data = df.P.jags18, aes(x = phos, y= mean), colour = "magenta2", size = 1) +
  geom_line(data = df.P.jags37, aes(x = phos, y= mean), colour = "forestgreen", size = 1) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Phosphorous concentration (µM)", 
       y = "Exponential growth rate (µ)",
       title = "B") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.p

# Salt

for (i in c(7,18,37)){      # Salt
  load(paste0("R2jags-objects/pop_", i, "_salt_tolerance.RData"))
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:4, 2006),]
  df.jags.plot$salt <- seq(0, 10, 0.005)
  assign(paste0("df.S.jags", i), df.jags.plot)
}

df.r.s <- read.csv("data-processed/13a_salt_r_estimates.csv") # Growth data across salt levels
head(df.r.s) # population is the factor, population.number is the corresponding # (e.g. 7, 18, 37)

df.r.s <- df.r.s %>% 
  filter(population.number %in% c(7, 18, 37)) %>% 
  print()

p.s <- ggplot(df.r.s, aes(x = salt.lvl, y = r.exp, colour = population)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  scale_colour_manual(values = c("darkorange1", "magenta2", "forestgreen")) +
  ylim(-0.5, 2.5) +
  
  geom_line(data = df.S.jags7, aes(x = salt, y= mean), colour = "darkorange1", size = 1) +
  geom_line(data = df.S.jags18, aes(x = salt, y= mean), colour = "magenta2", size = 1) +
  geom_line(data = df.S.jags37, aes(x = salt, y= mean), colour = "forestgreen", size = 1) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Salt concentration (g/L)", 
       y = "Exponential growth rate (µ)",
       title = "C") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.s

fig.1 <- plot_grid(p.t, p.p, p.s, nrow = 1, align=hv, rel_widths = c(1,1,1))

fig.1

ggsave("figures/16_fig_1_sample_models.jpeg", fig.1, width = 15, height = 5) # PDF was rendering weird

# Figure 2: PCAs and RDAs -----------------------------------------------------------

# We'll need the full dataset with metabolic and pigment data.

df.pca <- df.final %>% select(T.br, r.max_T, I.comp, r.max_I, N.comp, r.max_N, P.comp, r.max_P, r.max_S, 
                              S.c.mod, chl.a, chl.b, luthein, mean.N.µg.l, mean.P.µg.l, evol.plt) # Prepare the data: selecting only the relevant columns
df.pca

df.pca <- df.pca[-31,] # For now removing wonky/missing points (I.comp is way too high in row 31)

evol.fil <- df.pca$evol.plt

df.pca <- df.pca %>% select(-evol.plt)

pca.result <- prcomp(df.pca, center = TRUE, scale. = TRUE) # Perform PCA
pca.result

sdev <- pca.result$sdev  # or just use the vector if you already assigned it

sdev^2 %>%
  { . / sum(.) } %>% # The amount of variation explained by each PC
  sum(.)
  
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
                            levels = c("T.br", "r.max_T", "I.comp", "r.max_I", "N.comp", "r.max_N", 
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
                                       "tolerance",
                                       "Chlorophyll~italic(a)",
                                       "Chlorophyll~italic(b)",
                                       "Luthein", 
                                       "N~content",
                                       "P~content"))

adjust.x <- c(0.55, -0.05 , 0.36, -0.05, -0.02, 0, -0.02, 0.7, 0.8, 0.3, -0.05, 0.5, -0.05, -0.05, -0.05) # adjustments for each label
adjust.y <- c(-0.1, 0.15, 0.15, 0.1, 0.15, -0.05, 0.15, 0.02, -0.01, 0.24, 0.05, -0.15, 0.05, 0.09, 0.1)

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
  theme(legend.position= c(0.15, 0.75),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 12, face = "plain"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")) +
  geom_text(aes(x = -0.885, y = 0.95),  # adjust these coords
                     label = "Salt",
                     size = 5,
                     color = "black")

pca_plot_arrows

ggsave("figures/17_fig_2a_PCA.jpeg", pca_plot_arrows, width = 12, height = 9)

###### Evolutionary history RDA #######

response_vars <- df.pca
explanatory_vars <- model.matrix(~ evol.fil)[, -1]  # Remove intercept

rda_result_evol <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result_evol)
sum(summary(rda_result_evol)$cont$importance[2, 1:rda_result_evol$CCA$rank]) # with P and N, evolutionary environment explains 20.2712% of the variation

rda_var_explained <- summary(rda_result_evol)$cont$importance["Proportion Explained", 1:2] * 100

rda_sites_evol <- as.data.frame(scores(rda_result_evol, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_evol <- as.data.frame(scores(rda_result_evol, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_evol <- as.data.frame(scores(rda_result_evol, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_evol$Evolution <- evol.fil # Add the evolutionary treatment labels to the site scores

rda_constraints_evol$label <- rownames(rda_constraints_evol) # Assign readable labels for centroids

rda_species_evol$metric <- factor(rownames(rda_species_evol), 
                          levels = c("T.br", "r.max_T", "I.comp", "r.max_I", "N.comp", "r.max_N", 
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
sum(summary(rda_result_evol)$cont$importance[2, 1:rda_result_evol$CCA$rank]) # without P and N, evolutionary environment explains 20.58758% of the variation

rda_sites_evol <- as.data.frame(scores(rda_result_evol, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_evol <- as.data.frame(scores(rda_result_evol, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_evol <- as.data.frame(scores(rda_result_evol, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_evol$Evolution <- evol.fil # Add the evolutionary treatment labels to the site scores

rda_constraints_evol$label <- rownames(rda_constraints_evol) # Assign readable labels for centroids

rda_species_evol$metric <- factor(rownames(rda_species_evol), 
                                  levels = c("T.br", "r.max_T", "I.comp", "r.max_I", "N.comp", "r.max_N", 
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

df.pca2 <- df.final %>% select(T.br, r.max_T, I.comp, r.max_I, N.comp, r.max_N, P.comp, r.max_P, r.max_S, 
                              S.c.mod, chl.a, chl.b, luthein, mean.N.µg.l, mean.P.µg.l, anc.plt) # Prepare the data: selecting only the relevant columns
df.pca2

df.pca2 <- df.pca2[-31,] # For now removing wonky/missing points (I.comp is way too high in row 31)

anc.fil <- df.pca2$anc.plt

df.pca2 <- df.pca2 %>% select(-anc.plt)

response_vars <- df.pca2

explanatory_vars <- model.matrix(~ anc.fil)[, -1]  # Remove intercept

rda_result_anc <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result_anc)
sum(summary(rda_result_anc)$cont$importance[2, 1:rda_result_anc$CCA$rank]) # with P and N, ancestry explains 12.96264% of the variation

rda_sites_anc <- as.data.frame(scores(rda_result_anc, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_anc <- as.data.frame(scores(rda_result_anc, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_anc <- as.data.frame(scores(rda_result_anc, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_anc$Ancestry <- anc.fil # Add the evolutionary treatment labels to the site scores

rda_constraints_anc$label <- rownames(rda_constraints_anc) # Assign readable labels for centroids

rda_species_anc$metric <- factor(rownames(rda_species_anc), 
                                  levels = c("T.br", "r.max_T", "I.comp", "r.max_I", "N.comp", "r.max_N", 
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
sum(summary(rda_result_anc)$cont$importance[2, 1:rda_result_anc$CCA$rank]) # without P and N, ancestry explains 2.886395% of the variation

rda_sites_anc <- as.data.frame(scores(rda_result_anc, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_anc <- as.data.frame(scores(rda_result_anc, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_anc <- as.data.frame(scores(rda_result_anc, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_anc$Ancestry <- anc.fil # Add the evolutionary treatment labels to the site scores

rda_constraints_anc$label <- rownames(rda_constraints_anc) # Assign readable labels for centroids

rda_species_anc$metric <- factor(rownames(rda_species_anc), 
                                 levels = c("T.br", "r.max_T", "I.comp", "r.max_I", "N.comp", "r.max_N", 
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

df.filt <- df.final # Filter out extreme outliers

par.res.T <- par_frt(df.filt, xvar = "r.max_T", yvar = "T.br") # Get the raw Pareto Front

fit <- scam(T.br ~ s(r.max_T, bs = "mpd", k = 6), data = par.res.T) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$r.max_T), max(df.filt$r.max_T), length.out = 100) # Generate an x sequence for plotting

pred.curve.t <- data.frame( # Get the corresponding y values
  r.max_T = x.vals,
  T.br = predict(fit, newdata = data.frame(r.max_T = x.vals))
)

buffer <- 0.05  # 5% vertical tolerance - we'll include points within this range

pred.curve.t$threshold <- pred.curve.t$T.br - pred.curve.t$T.br * buffer # calculate threshold

df.filt$closest_idx <- sapply(df.filt$r.max_T, function(x) { # Match data in our df.filt to the predicted curve
  which.min(abs(pred.curve.t$r.max_T - x))
})

df.filt$predicted_T.br <- pred.curve.t$T.br[df.filt$closest_idx] # Get the predicted T.br
df.filt$threshold <- pred.curve.t$threshold[df.filt$closest_idx] # Get the threshold point

df.filt$within_band <- df.filt$T.br >= df.filt$threshold # Label points based on proximity to curve
df.filt2 <- df.filt[df.filt$within_band, ] # smaller subset of points close to the Pareto front

fit2 <- scam(T.br ~ s(r.max_T, bs = "mpd", k = 6), data = df.filt2)

pred.curve2 <- data.frame(
  r.max_T = x.vals,
  T.br = predict(fit2, newdata = data.frame(r.max_T = x.vals))
)

T_par <- ggplot(df.final, aes(x = r.max_T, y = T.br)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  geom_line(data=pred.curve, aes(x=r.max_T, y= T.br), colour = "black", size=1) +
  geom_line(data=pred.curve2, aes(x=r.max_T, y= T.br), colour="forestgreen", size = 1) +
  labs(x = "Maximum exponential growth rate", 
       y = "Thermal breadth", 
       title = "Thermal performance") +
  scale_shape_manual(values = c(16, 3)) +  # Keep custom shapes
  geom_line(data = par.res.T, aes(x = r.max_T, y = T.br), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

T_par # Raw pareto fronts.

# We have our 2 fits, now I want to remove the top 25% of points

df.filt <- df.filt %>% # Scale the data
  mutate(
    z.x = scale(r.max_T)[, 1],
    z.y = scale(T.br)[, 1],
  )

x.ref <- min(df.filt$z.x, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y, na.rm = TRUE) # Min y

df.filt3 <- df.filt %>% # Filter out based on Euclidean distance from min
  mutate(
    distance = sqrt((z.x - x.ref)^2 + (z.y - y.ref)^2)
  ) %>%
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.T2 <- par_frt(df.filt3, xvar = "r.max_T", yvar = "T.br") # Pareto frontier on this data

fit3 <- scam(T.br ~ s(r.max_T, bs = "mpd", k = 5), data = par.res.T2) # Model fit

pred.curve3 <- data.frame( # predicted data frame
  r.max_T = x.vals,
  T.br = predict(fit3, newdata = data.frame(r.max_T = x.vals))
)

pred.curve3$threshold <- pred.curve3$T.br - pred.curve3$T.br * buffer # predicted threshold

df.filt3$closest_idx <- sapply(df.filt3$r.max_T, function(x) { # Get closest ID
  which.min(abs(pred.curve3$r.max_T - x))
})

df.filt3$predicted_T.br <- pred.curve3$T.br[df.filt3$closest_idx] # Get corresponding predicted T.br
df.filt3$threshold <- pred.curve3$threshold[df.filt3$closest_idx] # Associated threshold

df.filt3$within_band <- df.filt3$T.br >= df.filt3$threshold # Label based on proximity to PF
df.filt4 <- df.filt3[df.filt3$within_band, ] # Filter to include only close points

fit4 <- scam(T.br ~ s(r.max_T, bs = "mpd", k = 6), data = df.filt4) # Model PF on that

pred.curve4 <- data.frame( # assemble a dataframe
  r.max_T = x.vals,
  T.br = predict(fit4, newdata = data.frame(r.max_T = x.vals))
)

T.scam <- ggplot(df.filt, aes(x = r.max_T, y = T.br, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.t, aes(x = r.max_T, y = T.br), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  # geom_line(data = pred.curve2, aes(x = r.max_T, y = T.br), color = "black", size = 1.1, inherit.aes = FALSE) +
  geom_line(data = pred.curve3, aes(x = r.max_T, y = T.br), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  # geom_line(data = pred.curve4, aes(x = r.max_T, y = T.br), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
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
  
  theme_classic() +
  ylim(12,20) +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.scam  # Display the plot

# Quadrant-based statistical testing

x.thresh.t <- mean(df.final$r.max_T)
y.thresh.t <- mean(df.final$T.br)

obs.cnt.t <- df.final %>%
  filter(r.max_T > x.thresh.t, T.br > y.thresh.t) %>%
  nrow()

null.counts.t <- replicate(1000, {
  shuffled.df <- df.final %>%
    mutate(r.max_T = sample(r.max_T, rep= F),
           T.br = sample(T.br, rep = F))
  
  sum(shuffled.df$r.max_T > x.thresh.t & shuffled.df$T.br > y.thresh.t)
})

p.val.t.quad <- mean(null.counts.t >= obs.cnt.t)
p.val.t.quad # 0.896

###### Light ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "L", 'light', 'other')) # for testing regressions.

df.filt <- df.final[df.final$I.comp <10, ] # Filter out extreme outliers

par.res.I <- par_frt(df.filt, xvar = "r.max_I", yvar = "I.comp") # Get the raw Pareto Front

fit <- scam(I.comp ~ s(r.max_I, bs = "mpd", k = 6), data = par.res.I) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$r.max_I), max(df.filt$r.max_I), length.out = 100) # Generate an x sequence for plotting

pred.curve.i <- data.frame( # Get the corresponding y values
  r.max_I = x.vals,
  I.comp = predict(fit, newdata = data.frame(r.max_I = x.vals))
)

pred.curve.i$threshold <- pred.curve.i$I.comp - pred.curve.i$I.comp * buffer # calculate threshold

df.filt$closest_idx <- sapply(df.filt$r.max_I, function(x) { # Match data in our df.filt to the predicted curve
  which.min(abs(pred.curve.i$r.max_I - x))
})

df.filt$predicted_I.comp <- pred.curve.i$I.comp[df.filt$closest_idx] # Get the predicted I.comp
df.filt$threshold <- pred.curve.i$threshold[df.filt$closest_idx] # Get the threshold point

df.filt$within_band <- df.filt$I.comp >= df.filt$threshold # Label points based on proximity to curve
df.filt2 <- df.filt[df.filt$within_band, ] # smaller subset of points close to the Pareto front

fit2 <- scam(I.comp ~ s(r.max_I, bs = "mpd", k = 4), data = df.filt2)

pred.curve2 <- data.frame(
  r.max_I = x.vals,
  I.comp = predict(fit2, newdata = data.frame(r.max_I = x.vals))
)

df.filt <- df.filt %>% # Scale the data
  mutate(
    z.x = scale(r.max_I)[, 1],
    z.y = scale(I.comp)[, 1],
  )

x.ref <- min(df.filt$z.x, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y, na.rm = TRUE) # Min y

df.filt3 <- df.filt %>% # Filter out based on Euclidean distance from min
  mutate(
    distance = sqrt((z.x - x.ref)^2 + (z.y - y.ref)^2)
  ) %>%
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.I2 <- par_frt(df.filt3, xvar = "r.max_I", yvar = "I.comp") # Pareto frontier on this data

fit3 <- scam(I.comp ~ s(r.max_I, bs = "mpd", k = 4), data = par.res.I2) # Model fit

pred.curve3 <- data.frame( # predicted data frame
  r.max_I = x.vals,
  I.comp = predict(fit3, newdata = data.frame(r.max_I = x.vals))
)

pred.curve3$threshold <- pred.curve3$I.comp - pred.curve3$I.comp * buffer # predicted threshold

df.filt3$closest_idx <- sapply(df.filt3$r.max_I, function(x) { # Get closest ID
  which.min(abs(pred.curve3$r.max_I - x))
})

df.filt3$predicted_I.comp <- pred.curve3$I.comp[df.filt3$closest_idx] # Get corresponding predicted I.comp
df.filt3$threshold <- pred.curve3$threshold[df.filt3$closest_idx] # Associated threshold

df.filt3$within_band <- df.filt3$I.comp >= df.filt3$threshold # Label based on proximity to PF
df.filt4 <- df.filt3[df.filt3$within_band, ] # Filter to include only close points

fit4 <- scam(I.comp ~ s(r.max_I, bs = "mpd", k = 2), data = df.filt4) # Model PF on that

pred.curve4 <- data.frame( # assemble a dataframe
  r.max_I = x.vals,
  I.comp = predict(fit4, newdata = data.frame(r.max_I = x.vals))
)

I.scam <- ggplot(df.filt, aes(x = r.max_I, y = I.comp, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.i, aes(x = r.max_I, y = I.comp), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  # geom_line(data = pred.curve2, aes(x = r.max_I, y = I.comp), color = "black", size = 1.1) +
  geom_line(data = pred.curve3, aes(x = r.max_I, y = I.comp), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  # geom_line(data = pred.curve4, aes(x = r.max_I, y = I.comp), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
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
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

I.scam  # Display the plot

l.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin == "light") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(r.max_I, pred.curve3$r.max_I),
    I.comp.pred = pred.curve3$I.comp[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(I.comp > I.comp.pred) %>%
  nrow()

set.seed(123)  # for reproducibility
n_iter <- 1000
null_counts_l <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts_l[i] <- df.shuff %>%
    filter(evol.shuff == "light") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve3$r.max_I - r.max_I)),
      I.comp.pred = pred.curve3$I.comp[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(I.comp > I.comp.pred) %>%
    nrow()
}

p_val_l <- mean(null_counts_l >= l.75) 
p_val_l # 0.378

# Quadrant-based statistical testing

x.thresh.l <- mean(df.filt$r.max_I)
y.thresh.l <- mean(df.filt$I.comp)

obs.cnt.l <- df.filt %>%
  filter(r.max_I > x.thresh.l, I.comp > y.thresh.l) %>%
  nrow()

null.counts.l <- replicate(1000, {
  shuffled.df <- df.filt %>%
    mutate(r.max_I = sample(r.max_I, rep= F),
           I.comp = sample(I.comp, rep = F))
  
  sum(shuffled.df$r.max_I > x.thresh.l & shuffled.df$I.comp > y.thresh.l)
})

p.val.l.quad <- mean(null.counts.l >= obs.cnt.l)
p.val.l.quad # 0.58

###### Nitrogen ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "N", 'nitrogen', 'other')) # for testing regressions.

df.filt <- df.final[df.final$N.comp < 1, ] # Filter out extreme outliers

par.res.N <- par_frt(df.filt, xvar = "r.max_N", yvar = "N.comp") # Get the raw Pareto Front

fit <- scam(N.comp ~ s(r.max_N, bs = "mpd", k = 5), data = par.res.N) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$r.max_N), max(df.filt$r.max_N), length.out = 100) # Generate an x sequence for plotting

pred.curve.n <- data.frame( # Get the corresponding y values
  r.max_N = x.vals,
  N.comp = predict(fit, newdata = data.frame(r.max_N = x.vals))
)

pred.curve.n$threshold <- pred.curve.n$N.comp - pred.curve.n$N.comp * buffer # calculate threshold

df.filt$closest_idx <- sapply(df.filt$r.max_N, function(x) { # Match data in our df.filt to the predicted curve
  which.min(abs(pred.curve.n$r.max_N - x))
})

df.filt$predicted_N.comp <- pred.curve.n$N.comp[df.filt$closest_idx] # Get the predicted N.comp
df.filt$threshold <- pred.curve.n$threshold[df.filt$closest_idx] # Get the threshold point

df.filt$within_band <- df.filt$N.comp >= df.filt$threshold # Label points based on proximity to curve
df.filt2 <- df.filt[df.filt$within_band, ] # smaller subset of points close to the Pareto front

fit2 <- scam(N.comp ~ s(r.max_N, bs = "mpd", k = 3), data = df.filt2)

pred.curve2 <- data.frame(
  r.max_N = x.vals,
  N.comp = predict(fit2, newdata = data.frame(r.max_N = x.vals))
)

df.filt <- df.filt %>% # Scale the data
  mutate(
    z.x = scale(r.max_N)[, 1],
    z.y = scale(N.comp)[, 1],
  )

x.ref <- min(df.filt$z.x, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y, na.rm = TRUE) # Min y

df.filt3 <- df.filt %>% # Filter out based on Euclidean distance from min
  mutate(
    distance = sqrt((z.x - x.ref)^2 + (z.y - y.ref)^2)
  ) %>%
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.N2 <- par_frt(df.filt3, xvar = "r.max_N", yvar = "N.comp") # Pareto frontier on this data

fit3 <- scam(N.comp ~ s(r.max_N, bs = "mpd", k = 5), data = par.res.N2) # Model fit

pred.curve3 <- data.frame( # predicted data frame
  r.max_N = x.vals,
  N.comp = predict(fit3, newdata = data.frame(r.max_N = x.vals))
)

pred.curve3$threshold <- pred.curve3$N.comp - pred.curve3$N.comp * buffer # predicted threshold

df.filt3$closest_idx <- sapply(df.filt3$r.max_N, function(x) { # Get closest ID
  which.min(abs(pred.curve3$r.max_N - x))
})

df.filt3$predicted_N.comp <- pred.curve3$N.comp[df.filt3$closest_idx] # Get corresponding predicted N.comp
df.filt3$threshold <- pred.curve3$threshold[df.filt3$closest_idx] # Associated threshold

df.filt3$within_band <- df.filt3$N.comp >= df.filt3$threshold # Label based on proximity to PF
df.filt4 <- df.filt3[df.filt3$within_band, ] # Filter to include only close points

fit4 <- scam(N.comp ~ s(r.max_N, bs = "mpd", k = 6), data = df.filt4) # Model PF on that

pred.curve4 <- data.frame( # assemble a dataframe
  r.max_N = x.vals,
  N.comp = predict(fit4, newdata = data.frame(r.max_N = x.vals))
)

N.scam <- ggplot(df.filt, aes(x = r.max_N, y = N.comp, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.n, aes(x = r.max_N, y = N.comp), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  # geom_line(data = pred.curve2, aes(x = r.max_N, y = N.comp), color = "black", size = 1.1, inherit.aes = FALSE) +
  geom_line(data = pred.curve3, aes(x = r.max_N, y = N.comp), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  # geom_line(data = pred.curve4, aes(x = r.max_N, y = N.comp), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
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
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

N.scam  # Display the plot

n.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin == "nitrogen") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(r.max_N, pred.curve3$r.max_N),
    N.comp.pred = pred.curve3$N.comp[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(N.comp > N.comp.pred) %>%
  nrow()

set.seed(123)  # for reproducibility
n_iter <- 1000
null_counts_n <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts_n[i] <- df.shuff %>%
    filter(evol.shuff == "nitrogen") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve3$r.max_N - r.max_N)),
      N.comp.pred = pred.curve3$N.comp[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(N.comp > N.comp.pred) %>%
    nrow()
}

p_val_n <- mean(null_counts_n >= n.75) 
p_val_n # 0.354

# Quadrant-based statistical testing

x.thresh.n <- mean(df.filt$r.max_N)
y.thresh.n <- mean(df.filt$N.comp)

obs.cnt.n <- df.filt %>%
  filter(r.max_N > x.thresh.n, N.comp > y.thresh.n) %>%
  nrow()

null.counts.n <- replicate(1000, {
  shuffled.df <- df.filt %>%
    mutate(r.max_N = sample(r.max_N, rep= F),
           N.comp = sample(N.comp, rep = F))
  
  sum(shuffled.df$r.max_N > x.thresh.n & shuffled.df$N.comp > y.thresh.n)
})

p.val.n.quad <- mean(null.counts.n >= obs.cnt.n)
p.val.n.quad # 0.956

###### Phosphorous ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "P", 'phosphorous', 'other')) # for testing regressions.

df.filt <- df.final # Filter out extreme outliers

par.res.P <- par_frt(df.filt, xvar = "r.max_P", yvar = "P.comp") # Get the raw Pareto Front

fit <- scam(P.comp ~ s(r.max_P, bs = "mpd", k = 6), data = par.res.P) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$r.max_P), max(df.filt$r.max_P), length.out = 100) # Generate an x sequence for plotting

pred.curve.p <- data.frame( # Get the corresponding y values
  r.max_P = x.vals,
  P.comp = predict(fit, newdata = data.frame(r.max_P = x.vals))
)

pred.curve.p$threshold <- pred.curve.p$P.comp - pred.curve.p$P.comp * buffer # calculate threshold

df.filt$closest_idx <- sapply(df.filt$r.max_P, function(x) { # Match data in our df.filt to the predicted curve
  which.min(abs(pred.curve.p$r.max_P - x))
})

df.filt$predicted_P.comp <- pred.curve.p$P.comp[df.filt$closest_idx] # Get the predicted P.comp
df.filt$threshold <- pred.curve.p$threshold[df.filt$closest_idx] # Get the threshold point

df.filt$within_band <- df.filt$P.comp >= df.filt$threshold # Label points based on proximity to curve
df.filt2 <- df.filt[df.filt$within_band, ] # smaller subset of points close to the Pareto front

fit2 <- scam(P.comp ~ s(r.max_P, bs = "mpd", k = 6), data = df.filt2)

pred.curve2 <- data.frame(
  r.max_P = x.vals,
  P.comp = predict(fit2, newdata = data.frame(r.max_P = x.vals))
)

df.filt <- df.filt %>% # Scale the data
  mutate(
    z.x = scale(r.max_P)[, 1],
    z.y = scale(P.comp)[, 1],
  )

x.ref <- min(df.filt$z.x, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y, na.rm = TRUE) # Min y

df.filt3 <- df.filt %>% # Filter out based on Euclidean distance from min
  mutate(
    distance = sqrt((z.x - x.ref)^2 + (z.y - y.ref)^2)
  ) %>%
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.P2 <- par_frt(df.filt3, xvar = "r.max_P", yvar = "P.comp") # Pareto frontier on this data

fit3 <- scam(P.comp ~ s(r.max_P, bs = "mpd", k = 5), data = par.res.P2) # Model fit

pred.curve3 <- data.frame( # predicted data frame
  r.max_P = x.vals,
  P.comp = predict(fit3, newdata = data.frame(r.max_P = x.vals))
)

pred.curve3$threshold <- pred.curve3$P.comp - pred.curve3$P.comp * buffer # predicted threshold

df.filt3$closest_idx <- sapply(df.filt3$r.max_P, function(x) { # Get closest ID
  which.min(abs(pred.curve3$r.max_P - x))
})

df.filt3$predicted_P.comp <- pred.curve3$P.comp[df.filt3$closest_idx] # Get corresponding predicted P.comp
df.filt3$threshold <- pred.curve3$threshold[df.filt3$closest_idx] # Associated threshold

df.filt3$within_band <- df.filt3$P.comp >= df.filt3$threshold # Label based on proximity to PF
df.filt4 <- df.filt3[df.filt3$within_band, ] # Filter to include only close points

fit4 <- scam(P.comp ~ s(r.max_P, bs = "mpd", k = 6), data = df.filt4) # Model PF on that

pred.curve4 <- data.frame( # assemble a dataframe
  r.max_P = x.vals,
  N.comp = predict(fit4, newdata = data.frame(r.max_P = x.vals))
)

P.scam <- ggplot(df.filt, aes(x = r.max_P, y = P.comp, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.p, aes(x = r.max_P, y = P.comp), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  # geom_line(data = pred.curve2, aes(x = r.max_P, y = P.comp), color = "black", size = 1.1, inherit.aes = FALSE) +
  geom_line(data = pred.curve3, aes(x = r.max_P, y = P.comp), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  # geom_line(data = pred.curve4, aes(x = r.max_P, y = P.comp), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
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
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

P.scam  # Display the plot

p.75 <- df.filt %>% # Now calculate the number of P points above that. 
  filter(evol.bin == "phosphorous") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(r.max_P, pred.curve3$r.max_P),
    P.comp.pred = pred.curve3$P.comp[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(P.comp > P.comp.pred) %>%
  nrow()

set.seed(123)  # for reproducibility
n_iter <- 1000
null_counts_p <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts_p[i] <- df.shuff %>%
    filter(evol.shuff == "phosphorous") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve3$r.max_P - r.max_P)),
      P.comp.pred = pred.curve3$P.comp[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(P.comp > P.comp.pred) %>%
    nrow()
}

p_val_p <- mean(null_counts_p >= p.75) 
p_val_p # 0.433

# Quadrant-based statistical testing

x.thresh.p <- mean(df.filt$r.max_P)
y.thresh.p <- mean(df.filt$P.comp)

obs.cnt.p <- df.filt %>%
  filter(r.max_P > x.thresh.p, P.comp > y.thresh.p) %>%
  nrow()

null.counts.p <- replicate(1000, {
  shuffled.df <- df.filt %>%
    mutate(r.max_P = sample(r.max_P, rep= F),
           P.comp = sample(P.comp, rep = F))
  
  sum(shuffled.df$r.max_P > x.thresh.p & shuffled.df$P.comp > y.thresh.p)
})

p.val.p.quad <- mean(null.counts.p >= obs.cnt.p)
p.val.p.quad # 0.977

###### Salt ######

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                   ifelse(df.final$evol %in% c("S", "BS"), 'salt', 'other')) # Ss and BSs are now equivalent 

df.filt <- df.final # Filter out extreme outliers

par.res.S <- par_frt(df.filt, xvar = "r.max_S", yvar = "S.c.mod") # Get the raw Pareto Front

fit <- scam(S.c.mod ~ s(r.max_S, bs = "mpd", k = 5), data = par.res.S) # Fit a scam to the raw Pareto front

x.vals <- seq(min(df.filt$r.max_S), max(df.filt$r.max_S), length.out = 100) # Generate an x sequence for plotting

pred.curve.s <- data.frame( # Get the corresponding y values
  r.max_S = x.vals,
  S.c.mod = predict(fit, newdata = data.frame(r.max_S = x.vals))
)

pred.curve.s$threshold <- pred.curve.s$S.c.mod - pred.curve.s$S.c.mod * buffer # calculate threshold

df.filt$closest_idx <- sapply(df.filt$r.max_S, function(x) { # Match data in our df.filt to the predicted curve
  which.min(abs(pred.curve.s$r.max_S - x))
})

df.filt$predicted_S.c.mod <- pred.curve.s$S.c.mod[df.filt$closest_idx] # Get the predicted S.c.mod
df.filt$threshold <- pred.curve.s$threshold[df.filt$closest_idx] # Get the threshold point

df.filt$within_band <- df.filt$S.c.mod >= df.filt$threshold # Label points based on proximity to curve
df.filt2 <- df.filt[df.filt$within_band, ] # smaller subset of points close to the Pareto front

fit2 <- scam(S.c.mod ~ s(r.max_S, bs = "mpd", k = 6), data = df.filt2)

pred.curve2 <- data.frame(
  r.max_S = x.vals,
  S.c.mod = predict(fit2, newdata = data.frame(r.max_S = x.vals))
)

df.filt <- df.filt %>% # Scale the data
  mutate(
    z.x = scale(r.max_S)[, 1],
    z.y = scale(S.c.mod)[, 1],
  )

x.ref <- min(df.filt$z.x, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y, na.rm = TRUE) # Min y

df.filt3 <- df.filt %>% # Filter out based on Euclidean distance from min
  mutate(
    distance = sqrt((z.x - x.ref)^2 + (z.y - y.ref)^2)
  ) %>%
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.S2 <- par_frt(df.filt3, xvar = "r.max_S", yvar = "S.c.mod") # Pareto frontier on this data

fit3 <- scam(S.c.mod ~ s(r.max_S, bs = "mpd", k = 6), data = par.res.S2) # Model fit

pred.curve3 <- data.frame( # predicted data frame
  r.max_S = x.vals,
  S.c.mod = predict(fit3, newdata = data.frame(r.max_S = x.vals))
)

pred.curve3$threshold <- pred.curve3$S.c.mod - pred.curve3$S.c.mod * buffer # predicted threshold

df.filt3$closest_idx <- sapply(df.filt3$r.max_S, function(x) { # Get closest ID
  which.min(abs(pred.curve3$r.max_S - x))
})

df.filt3$predicted_S.c.mod <- pred.curve3$S.c.mod[df.filt3$closest_idx] # Get corresponding predicted S.c.mod
df.filt3$threshold <- pred.curve3$threshold[df.filt3$closest_idx] # Associated threshold

df.filt3$within_band <- df.filt3$S.c.mod >= df.filt3$threshold # Label based on proximity to PF
df.filt4 <- df.filt3[df.filt3$within_band, ] # Filter to include only close points

fit4 <- scam(S.c.mod ~ s(r.max_S, bs = "mpd", k = 6), data = df.filt4) # Model PF on that

pred.curve4 <- data.frame( # assemble a dataframe
  r.max_S = x.vals,
  S.c.mod = predict(fit4, newdata = data.frame(r.max_S = x.vals))
)

S.scam <- ggplot(df.filt, aes(x = r.max_S, y = S.c.mod, color = evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.s, aes(x = r.max_S, y = S.c.mod), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  # geom_line(data = pred.curve2, aes(x = r.max_S, y = S.c.mod), color = "black", size = 1.1, inherit.aes = FALSE) +
  geom_line(data = pred.curve3, aes(x = r.max_S, y = S.c.mod), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  # geom_line(data = pred.curve4, aes(x = r.max_S, y = S.c.mod), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
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

s.75 <- df.filt %>% # Now calculate the number of light points above that. 
  filter(evol.bin == "salt") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(r.max_S, pred.curve3$r.max_S),
    S.c.mod.pred = pred.curve3$S.c.mod[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(S.c.mod > S.c.mod.pred) %>%
  nrow()

set.seed(123)  # for reproducibility
n_iter <- 1000
null_counts_s <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.filt %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts_s[i] <- df.shuff %>%
    filter(evol.shuff == "salt") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve3$r.max_S - r.max_S)),
      S.c.mod.pred = pred.curve3$S.c.mod[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(S.c.mod > S.c.mod.pred) %>%
    nrow()
}

p_val_s <- mean(null_counts_s >= s.75) 
p_val_s # 0

# Quadrant-based statistical testing

x.thresh.s <- mean(df.filt$r.max_S)
y.thresh.s <- mean(df.filt$S.c.mod)

obs.cnt.s <- df.filt %>%
  filter(r.max_S > x.thresh.s, S.c.mod > y.thresh.s) %>%
  nrow()

null.counts.s <- replicate(1000, {
  shuffled.df <- df.filt %>%
    mutate(r.max_S = sample(r.max_S, rep= F),
           S.c.mod = sample(S.c.mod, rep = F))
  
  sum(shuffled.df$r.max_S > x.thresh.s & shuffled.df$S.c.mod > y.thresh.s)
})

p.val.s.quad <- mean(null.counts.s >= obs.cnt.s)
p.val.s.quad # 0.835

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
head(df.thomas) # temperature only

df.bestion.t <- read.csv("data-processed/18a_Bestion2018_TPCs.csv")
head(df.bestion.t) # temperature

df.bestion.p <- read.csv("data-processed/18c_Bestion2018_P_Monodss.csv")
head(df.bestion.p) # phosphorous

df.lewington.t <- read.csv("data-processed/19b_Lewington2019_TPCs.csv")
head(df.lewington.t) # temperature

df.lewington.n <- read.csv("data-processed/19e_Lewington2019_Nit_Monods.csv")
head(df.lewington.n) # nitrogen, some negative R.mth scores here

df.lewington.l <- read.csv("data-processed/19h_Lewington2019_Light_Monods.csv")
head(df.lewington.l) # light, some negative R.mth scores here

df.narwani <- read.csv("data-processed/20_Narwani2015_summmary.csv")
head(df.narwani) # light, nitrogen, phosphate

df.kontopoulos <- read.csv("data-processed/21a_Kontopoulos2020_TPCs.csv")
head(df.kontopoulos) # temperature, low r.maxes

###### Temperature ######

# We have thomas, bestion, lewington and kontopoulos

sp.t <- ggplot(df.thomas, aes(x = r.max.raw, y = T.br.raw, colour = "Thomas et al., 2012")) + # Thomas 2012
  geom_point(size = 2) +
  geom_point(data= df.bestion.t, aes(x = r.max.raw, y = T.br.raw, colour = "Bestion et al., 2018"), size = 2) + # Bestion 2018
  geom_point(data= df.lewington.t, aes(x = r.max.raw, y = T.br.raw, colour = "Lewington-Pearce et al., 2019"), size = 2) + # Lewington 2019
  geom_point(data= df.kontopoulos, aes(x = r.max.raw, y = T.br.raw, colour = "Kontopoulos et al., 2020"), size = 2) + # Kontopoulos 2020
  
  geom_line(data = pred.t, aes(x = r.max_T, y = T.br.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       title = "D") +  # labels
  
  scale_color_manual(values = c("chartreuse2", "maroon1", "gold1", "cornflowerblue"), 
                     name = "Data set",
                     labels = c("Bestion et al., 2018", "Kontopoulos et al., 2020", "Lewington-Pearce et al., 2019", "Thomas et al., 2012")) +  
  theme_classic() +
  theme(
    legend.position = c(0.7, 0.8),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

sp.t  # Display the plot

df.final.T.scaled <- df.final %>%
  mutate(r.max_T.scaled = scale(r.max_T),
         T.br.scaled = scale(T.br))

quant.T.scaled <- rq(T.br.scaled ~ poly(r.max_T.scaled, 2),
                     data = df.final.T.scaled, tau = 0.95)

mean_rmax <- mean(df.final.T.scaled$r.max_T, na.rm = TRUE)
sd_rmax   <- sd(df.final.T.scaled$r.max_T, na.rm = TRUE)

pred.t2 <- pred.t %>%
  mutate(r.max_T.scaled = (r.max_T - mean_rmax) / sd_rmax)

pred.t2$T.br.scaled <- predict(quant.T.scaled, newdata = pred.t2)

sp.t.scaled <- ggplot(df.thomas, aes(x = scale(r.max.raw), y = scale(T.br.raw), colour = "Thomas et al., 2012")) + # Thomas 2012
  geom_point(size = 2) +
  geom_point(data= df.bestion.t, aes(x = scale(r.max.raw), y = scale(T.br.raw), colour = "Bestion et al., 2018"), size = 2) + # Bestion 2018
  geom_point(data= df.lewington.t, aes(x = scale(r.max.raw), y = scale(T.br.raw), colour = "Lewington-Pearce et al., 2019"), size = 2) + # Lewington 2019
  geom_point(data= df.kontopoulos, aes(x = scale(r.max.raw), y = scale(T.br.raw), colour = "Kontopoulos et al., 2020"), size = 2) + # Kontopoulos 2020
  
  geom_line(data = pred.t2, aes(x = r.max_T.scaled, y = T.br.scaled), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       title = "D") +  # labels
  
  scale_color_manual(values = c("chartreuse2", "maroon1", "gold1", "cornflowerblue"), 
                     name = "Data set",
                     labels = c("Bestion et al., 2018", "Kontopoulos et al., 2020", "Lewington-Pearce et al., 2019", "Thomas et al., 2012")) +  
  theme_classic() +
  theme(
    legend.position = c(0.2, 0.8),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

sp.t.scaled

###### Light ######

# Only the Narwani data and Lewington

sp.l <- ggplot(df.narwani, aes(x = umax.nitrate, y = Nstar, colour = "Narwani et al., 2015")) + # Narwani 2015
  geom_point(size = 2) +
  
  geom_point(data= df.lewington.l, aes(x = r.max, y = R.mth, colour = "Lewington-Pearce et al., 2019"), size = 2) + # Lewington 2019
  
  geom_line(data = pred.l, aes(x = r.max_I, y = I.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/R*)", , 
       title = "A") +  # labels
  
  scale_color_manual(values = c("gold1","darkgreen"), 
                     name = "Data set",
                     labels = c("Lewington-Pearce et al., 2019", "Narwani et al., 2015")) +  
  theme_classic() +
  theme(
    legend.position = c(0.9, 0.9),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

sp.l  # Display the plot

df.final.I.scaled <- df.final %>%
  filter(I.comp < 10) %>% 
  mutate(r.max_I.scaled = scale(r.max_I),
         I.comp.scaled = scale(I.comp))

quant.I.scaled <- rq(I.comp.scaled ~ poly(r.max_I.scaled, 2),
                     data = df.final.I.scaled, tau = 0.95)

mean_rmax <- mean(df.final.I.scaled$r.max_I, na.rm = TRUE)
sd_rmax   <- sd(df.final.I.scaled$r.max_I, na.rm = TRUE)

pred.l2 <- pred.l %>%
  mutate(r.max_I.scaled = (r.max_I - mean_rmax) / sd_rmax)

pred.l2$I.comp.scaled <- predict(quant.I.scaled, newdata = pred.l2)

pred.l2 <- pred.l2 %>%
  mutate(I.comp.scaled = ifelse(row_number() <= which.max(pred.l2$I.comp.scaled) & I.comp.scaled < max(I.comp.scaled), NA, I.comp.scaled))

sp.l.scaled <- ggplot(df.narwani, aes(x = scale(umax.nitrate), y = scale(Nstar), colour = "Narwani et al., 2015")) + # Narwani 2015
  geom_point(size = 2) +
  
  geom_point(data= df.lewington.l, aes(x = scale(r.max), y = scale(R.mth), colour = "Lewington-Pearce et al., 2019"), size = 2) + # Lewington 2019
  
  geom_line(data = pred.l2, aes(x = r.max_I.scaled, y = I.comp.scaled), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/R*)", , 
       title = "A") +  # labels
  
  scale_color_manual(values = c("gold1","darkgreen"), 
                     name = "Data set",
                     labels = c("Lewington-Pearce et al., 2019", "Narwani et al., 2015")) +  
  theme_classic() +
  theme(
    legend.position = c(0.25, 0.8),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

sp.l.scaled  # Display the plot

###### Nitrogen ######

# So we have Narwani and Lewington

sp.n <- ggplot(df.narwani, aes(x = umax.light, y = Istar, colour = "Narwani et al., 2015")) + # Narwani 2015
  geom_point(size = 2) +
  
  geom_point(data= df.lewington.n, aes(x = r.max, y = R.mth, colour = "Lewington-Pearce et al., 2019"), size = 2) + # Lewington 2019
  
  geom_line(data = pred.n, aes(x = r.max_N, y = N.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/R*)", , 
       title = "B") +  # labels
  
  scale_color_manual(values = c("gold1", "darkgreen"), 
                     name = "Data set",
                     labels = c("Lewington-Pearce et al., 2019", "Narwani et al., 2015")) +  
  
  theme_classic() +
  theme(
    legend.position = c(0.7, 0.5),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

sp.n  # Display the plot

df.final.N.scaled <- df.final %>%
  mutate(r.max_N.scaled = scale(r.max_N),
         N.comp.scaled = scale(N.comp))

quant.N.scaled <- rq(N.comp.scaled ~ poly(r.max_N.scaled, 2),
                     data = df.final.N.scaled, tau = 0.95)

mean_rmax <- mean(df.final.N.scaled$r.max_N, na.rm = TRUE)
sd_rmax   <- sd(df.final.N.scaled$r.max_N, na.rm = TRUE)

pred.n2 <- pred.n %>%
  mutate(r.max_N.scaled = (r.max_N - mean_rmax) / sd_rmax)

pred.n2$N.comp.scaled <- predict(quant.N.scaled, newdata = pred.n2)

sp.n.scaled <- ggplot(df.narwani, aes(x = scale(umax.light), y = scale(Istar), colour = "Narwani et al., 2015")) + # Narwani 2015
  geom_point(size = 2) +
  
  geom_point(data= df.lewington.n, aes(x = scale(r.max), y = scale(R.mth), colour = "Lewington-Pearce et al., 2019"), size = 2) + # Lewington 2019
  
  geom_line(data = pred.n2, aes(x = r.max_N.scaled, y = N.comp.scaled), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/R*)", , 
       title = "B") +  # labels
  
  scale_color_manual(values = c("gold1", "darkgreen"), 
                     name = "Data set",
                     labels = c("Lewington-Pearce et al., 2019", "Narwani et al., 2015")) +  
  
  theme_classic() +
  theme(
    legend.position = c(0.7, 0.5),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

sp.n.scaled  # Display the plot

###### Phosphorous ######

# Bestion and Narwani

sp.p <- ggplot(df.narwani, aes(x = umax.phosphate, y = Pstar, colour = "Narwani et al., 2015")) + # Narwani 2015
  geom_point(size = 2) +
  
  geom_point(data= df.bestion.p, aes(x = r.max, y = R.mth, colour = "Bestion et al., 2018"), size = 2) + # Bestion 2018
  
  geom_line(data = pred.p, aes(x = r.max_P, y = P.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/R*)", , 
       title = "C") +  # labels
  
  scale_color_manual(values = c("chartreuse2", "darkgreen"), 
                     name = "Data set",
                     labels = c("Lewington-Pearce et al., 2019", "Narwani et al., 2015")) +  
  
  theme_classic() +
  theme(
    legend.position = c(0.3, 0.8),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

sp.p  # Display the plot

df.final.P.scaled <- df.final %>%
  mutate(r.max_P.scaled = scale(r.max_P),
         P.comp.scaled = scale(P.comp))

quant.P.scaled <- rq(P.comp.scaled ~ poly(r.max_P.scaled, 2),
                     data = df.final.P.scaled, tau = 0.95)

mean_rmax <- mean(df.final.P.scaled$r.max_P, na.rm = TRUE)
sd_rmax   <- sd(df.final.P.scaled$r.max_P, na.rm = TRUE)

pred.p2 <- pred.p %>%
  mutate(r.max_P.scaled = (r.max_P - mean_rmax) / sd_rmax)

pred.p2$P.comp.scaled <- predict(quant.P.scaled, newdata = pred.p2)

sp.p.scaled <- ggplot(df.narwani, aes(x = scale(umax.phosphate), y = scale(Pstar), colour = "Narwani et al., 2015")) + # Narwani 2015
  geom_point(size = 2) +
  
  geom_point(data= df.bestion.p, aes(x = scale(r.max), y = scale(R.mth), colour = "Bestion et al., 2018"), size = 2) + # Bestion 2018
  
  geom_line(data = pred.p2, aes(x = r.max_P.scaled, y = P.comp.scaled), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/R*)", , 
       title = "C") +  # labels
  
  scale_color_manual(values = c("chartreuse2", "darkgreen"), 
                     name = "Data set",
                     labels = c("Lewington-Pearce et al., 2019", "Narwani et al., 2015")) +  
  
  theme_classic() +
  theme(
    legend.position = c(0.3, 0.8),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

sp.p.scaled  # Display the plot

###### Make the figures ######

plots2 <- list(sp.l, sp.n, sp.p, sp.t)

plots2_nolegend <- lapply(plots2, function(p) p + theme(legend.position = "none"))

legend_sp <- data.frame(
  x = 1:5,
  y = 1:5,
  Dataset = factor(c(
    "Bestion et al., 2018",
    "Kontopoulos et al., 2020",
    "Lewington-Pearce et al., 2019",
    "Narwani et al., 2015",
    "Thomas et al., 2012"
  ))
)

legend_sp_plot <- ggplot(legend_sp, aes(x = x, y = y, color = Dataset)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Bestion et al., 2018" = "chartreuse2",
               "Kontopoulos et al., 2020" = "maroon1",
               "Lewington-Pearce et al., 2019" = "gold1",
               "Narwani et al., 2015" = "darkgreen",
               "Thomas et al., 2012" = "cornflowerblue"),
    name = "Data set"
  ) +
  theme_void() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right"
  )


sp_legend <- get_legend(legend_sp_plot)

all_sp_plots <- c(plots2_nolegend, list(sp_legend))

grad_sp_toffs <- plot_grid(plotlist = all_sp_plots,
                        ncol = 2,
                        align = "hv")

ggsave("figures/19_fig_4a_sp_tradeoffs.jpeg", grad_sp_toffs, width = 8, height = 12)

plots2.scaled <- list(sp.l.scaled, sp.n.scaled, sp.p.scaled, sp.t.scaled)

plots2_nolegend.scaled <- lapply(plots2.scaled, function(p) p + theme(legend.position = "none"))

all_sp_plots.scaled <- c(plots2_nolegend.scaled, list(sp_legend))

grad_sp_toffs.scaled <- plot_grid(plotlist = all_sp_plots.scaled,
                           ncol = 2,
                           align = "hv")

ggsave("figures/19_fig_4b_sp_tradeoffs_scaled.jpeg", grad_sp_toffs.scaled, width = 8, height = 12)

# Figure 5: Inter-gradient trade-offs -------------------------------------

# I think for now we will limit this to competitive abilities, salt tolerance and thermal breadth.
# We'll create a half-full grid. 
# Order: light, nitrogen, phosphorous, salt, temperature

###### Light comparisons ######

# Light v Nitrogen

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "L", 'light', ifelse(df.final$evol == 'N', 'nit', 'other'))) # for testing regressions.

pred.ln <- data.frame(N.comp = seq(min(df$N.comp), max(df$N.comp), length.out = 100)) # Dataframe to collect quantile info in

quant.IN.95 <- rq(I.comp ~ poly(N.comp, 2), data = df.final[df.final$I.comp<10,], tau = 0.95) 
pred.ln$I.comp.95 <- predict(quant.IN.95, newdata = pred.ln)

quant.IN.75 <- rq(I.comp ~ poly(N.comp, 2), data = df.final[df.final$I.comp<10,], tau = 0.75) 
pred.ln$I.comp.75 <- predict(quant.IN.75, newdata = pred.ln)

df.final$evol.bin <- factor(df.final$evol.bin, levels = c("ancestral", "other", "light", "nit"))

LN.qrs <- ggplot(df.final[df.final$I.comp<10,], aes(x = N.comp, y = I.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.ln, aes(x = N.comp, y = I.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.ln, aes(x = N.comp, y = I.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Competitive ability (1/N*)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "A") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3", "dodgerblue"), 
                     labels = c("Ancestral", "Other", "Light", "Nitrogen")) +  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff # theme stuff
  )

LN.qrs  # Display the plot

# Light v Phosphorous

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "L", 'light', ifelse(df.final$evol == 'P', 'phos', 'other'))) # for testing regressions.

pred.lp <- data.frame(P.comp = seq(min(df$P.comp), max(df$P.comp), length.out = 100)) # Dataframe to collect quantile info in

quant.IP.95 <- rq(I.comp ~ poly(P.comp, 2), data = df.final[df.final$I.comp<10,], tau = 0.95) 
pred.lp$I.comp.95 <- predict(quant.IP.95, newdata = pred.lp)

quant.IP.75 <- rq(I.comp ~ poly(P.comp, 2), data = df.final[df.final$I.comp<10,], tau = 0.75) 
pred.lp$I.comp.75 <- predict(quant.IP.75, newdata = pred.lp)

df.final$evol.bin <- factor(df.final$evol.bin, levels = c("ancestral", "other", "light", "phos"))

LP.qrs <- ggplot(df.final[df.final$I.comp<10,], aes(x = P.comp, y = I.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.lp, aes(x = P.comp, y = I.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.lp, aes(x = P.comp, y = I.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "B") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3", "dodgerblue"), 
                     labels = c("Ancestral", "Other", "Light", "Phosphorous")) +  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff # theme stuff
  )

LP.qrs  # Display the plot

par.res.LP <- par_frt(df.final[df.final$I.comp<10,], xvar = "P.comp", yvar = "I.comp")

par.res.LP2 <- par_frt_tolerant(df.final[df.final$I.comp<10,], xvar = "P.comp", yvar = "I.comp")

pred.lp.pf <- data.frame(P.comp = seq((min(par.res.LP$P.comp)-sd(par.res.LP$P.comp)/2), (max(par.res.LP$P.comp)+sd(par.res.LP$P.comp)/2), length.out = 100)) # New df for the par front data

quant.IP.pf <- lm(I.comp ~ poly(P.comp, 2), data = par.res.LP) 
pred.lp.pf$I.comp.pf <- predict(quant.IP.pf, newdata = pred.lp.pf)

smooth.model <- smooth.spline(par.res.LP$P.comp, par.res.LP$I.comp, spar = 0.5)

x.vals <- seq(min(par.res.LP$P.comp), max(par.res.LP$P.comp), length.out = 100)
pred.curve <- data.frame(
  P.comp = x.vals,
  I.comp = predict(smooth.model, x = x.vals)$y
)

fit <- scam(I.comp ~ s(P.comp, bs = "mpd", k = 6), data = par.res.LP2)

x.vals <- seq(min(par.res.LP$P.comp), max(par.res.LP$P.comp), length.out = 100)
pred.curve <- data.frame(
  P.comp = x.vals,
  I.comp = predict(fit, newdata = data.frame(P.comp = x.vals))
)

pred.lp.pf$I.comp.pf2 <- predict(smooth.model, newdata = pred.lp.pf)

LP.qrs2 <- ggplot(df.final[df.final$I.comp<10,], aes(x = P.comp, y = I.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve, aes(x = P.comp, y = I.comp), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.lp, aes(x = P.comp, y = I.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "B") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3", "dodgerblue"), 
                     labels = c("Ancestral", "Other", "Light", "Phosphorous")) +  
  theme_classic() +
  theme(
    legend.position = "none",  # Move legend inside the plot
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff # theme stuff
  )

LP.qrs2

examp.par.front <- plot_grid(LP.qrs, LP.qrs2)

examp.par.front

ggsave("figures/20b_pareto_front_options.jpeg", examp.par.front, width = 12, height = 8)

# Light v salt

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "L", 'light', ifelse(df.final$evol %in% c("S", "BS"), 'salt', 'other'))) # for testing regressions.

pred.ls <- data.frame(S.c.mod = seq(min(df$S.c.mod), max(df$S.c.mod), length.out = 100)) # Dataframe to collect quantile info in

quant.IS.95 <- rq(I.comp ~ poly(S.c.mod, 2), data = df.final[df.final$I.comp<10,], tau = 0.95) 
pred.ls$I.comp.95 <- predict(quant.IS.95, newdata = pred.ls)

quant.IS.75 <- rq(I.comp ~ poly(S.c.mod, 2), data = df.final[df.final$I.comp<10,], tau = 0.75) 
pred.ls$I.comp.75 <- predict(quant.IS.75, newdata = pred.ls)

df.final$evol.bin <- factor(df.final$evol.bin, levels = c("ancestral", "other", "light", "phos"))

LS.qrs <- ggplot(df.final[df.final$I.comp<10,], aes(x = S.c.mod, y = I.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.ls, aes(x = S.c.mod, y = I.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.ls, aes(x = S.c.mod, y = I.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Salt tolerance (c)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "C") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3", "dodgerblue"), 
                     labels = c("Ancestral", "Other", "Light", "Salt or Biotic x Salt")) +  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff # theme stuff
  )

LS.qrs  # Display the plot

# Light v temperature

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "L", 'light', 'other')) # for testing regressions.

pred.lt <- data.frame(T.br = seq(min(df$T.br), max(df$T.br), length.out = 100)) # Dataframe to collect quantile info in

quant.IT.95 <- rq(I.comp ~ poly(T.br, 2), data = df.final[df.final$I.comp<10,], tau = 0.95) 
pred.lt$I.comp.95 <- predict(quant.IT.95, newdata = pred.lt)

quant.IT.75 <- rq(I.comp ~ poly(T.br, 2), data = df.final[df.final$I.comp<10,], tau = 0.75) 
pred.lt$I.comp.75 <- predict(quant.IT.75, newdata = pred.lt)

df.final$evol.bin <- factor(df.final$evol.bin, levels = c("ancestral", "other", "light"))

LT.qrs <- ggplot(df.final[df.final$I.comp<10,], aes(x = T.br, y = I.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.lt, aes(x = T.br, y = I.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.lt, aes(x = T.br, y = I.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Thermal breadth (°C)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "D") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3", "dodgerblue"), 
                     labels = c("Ancestral", "Other", "Light", "Temperature")) +  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff # theme stuff
  )

LT.qrs  # Display the plot

###### Nitrogen comparisons ######

# N v P

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "N", 'nit', ifelse(df.final$evol == 'P', 'phos', 'other'))) # for testing regressions.

pred.np <- data.frame(P.comp = seq(min(df$P.comp), max(df$P.comp), length.out = 100)) # Dataframe to collect quantile info in

quant.NP.95 <- rq(N.comp ~ poly(P.comp, 2), data = df.final, tau = 0.95) 
pred.np$N.comp.95 <- predict(quant.NP.95, newdata = pred.np)

quant.NP.75 <- rq(N.comp ~ poly(P.comp, 2), data = df.final, tau = 0.75) 
pred.np$N.comp.75 <- predict(quant.NP.75, newdata = pred.np)

df.final$evol.bin <- factor(df.final$evol.bin, levels = c("ancestral", "other", "nit", "phos"))

NP.qrs <- ggplot(df.final, aes(x = P.comp, y = N.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.np, aes(x = P.comp, y = N.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.np, aes(x = P.comp, y = N.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "E") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3", "dodgerblue"), 
                     labels = c("Ancestral", "Other", "Nitrogen", "Phosphorous")) +  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff # theme stuff
  )

NP.qrs  # Display the plot

# Nit v Salt

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "N", 'nit', ifelse(df.final$evol %in% c("S", "BS"), 'salt', 'other')))

pred.ns <- data.frame(S.c.mod = seq(min(df$S.c.mod), max(df$S.c.mod), length.out = 100)) # Dataframe to collect quantile info in

quant.NS.95 <- rq(N.comp ~ poly(S.c.mod, 2), data = df.final, tau = 0.95) 
pred.ns$N.comp.95 <- predict(quant.NP.95, newdata = pred.np)

quant.NS.75 <- rq(N.comp ~ poly(S.c.mod, 2), data = df.final, tau = 0.75) 
pred.ns$N.comp.75 <- predict(quant.NP.75, newdata = pred.np)

df.final$evol.bin <- factor(df.final$evol.bin, levels = c("ancestral", "other", "nit", "salt"))

NS.qrs <- ggplot(df.final, aes(x = S.c.mod, y = N.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.ns, aes(x = S.c.mod, y = N.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.ns, aes(x = S.c.mod, y = N.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Salt tolerance (c)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "F") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3", "dodgerblue"), 
                     labels = c("Ancestral", "Other", "Nitrogen", "Salt or Biotic x Salt")) +  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff # theme stuff
  )

NS.qrs  # Display the plot

# Nitrogen v temp

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "N", 'nit', 'other')) # for testing regressions.

pred.nt <- data.frame(T.br = seq(min(df$T.br), max(df$T.br), length.out = 100)) # Dataframe to collect quantile info in

quant.NT.95 <- rq(N.comp ~ poly(T.br, 2), data = df.final[df.final$N.comp<10,], tau = 0.95) 
pred.nt$N.comp.95 <- predict(quant.NT.95, newdata = pred.nt)

quant.NT.75 <- rq(N.comp ~ poly(T.br, 2), data = df.final[df.final$N.comp<10,], tau = 0.75) 
pred.nt$N.comp.75 <- predict(quant.NT.75, newdata = pred.nt)

df.final$evol.bin <- factor(df.final$evol.bin, levels = c("ancestral", "other", "light"))

NT.qrs <- ggplot(df.final[df.final$N.comp<10,], aes(x = T.br, y = N.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.nt, aes(x = T.br, y = N.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.nt, aes(x = T.br, y = N.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Thermal breadth (°C)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "G") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Nitrogen")) +  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff # theme stuff
  )

NT.qrs  # Display the plot

###### Phosphorous comparisons ######

# P v salt

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "P", 'phos', ifelse(df.final$evol %in% c("S", "BS"), 'salt', 'other')))

pred.ps <- data.frame(S.c.mod = seq(min(df$S.c.mod), max(df$S.c.mod), length.out = 100)) # Dataframe to collect quantile info in

quant.PS.95 <- rq(P.comp ~ poly(S.c.mod, 2), data = df.final, tau = 0.95) 
pred.ps$P.comp.95 <- predict(quant.PS.95, newdata = pred.ns)

quant.PS.75 <- rq(P.comp ~ poly(S.c.mod, 2), data = df.final, tau = 0.75) 
pred.ps$P.comp.75 <- predict(quant.PS.75, newdata = pred.ns)

df.final$evol.bin <- factor(df.final$evol.bin, levels = c("ancestral", "other", "phos", "salt"))

PS.qrs <- ggplot(df.final, aes(x = S.c.mod, y = P.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.ps, aes(x = S.c.mod, y = P.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.ps, aes(x = S.c.mod, y = P.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Salt tolerance (c)",    
       y = "Competitive ability (1/P*)", 
       color = "Evolutionary History",
       title = "H") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3", "dodgerblue"), 
                     labels = c("Ancestral", "Other", "Phosphorous", "Salt or Biotic x Salt")) +  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff # theme stuff
  )

PS.qrs  # Display the plot

plot_grid(LN.qrs, LP.qrs, LS.qrs, LT.qrs, NP.qrs, NS.qrs, NT.qrs, PS.qrs)

# P v T

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "P", 'phos', 'other')) # for testing regressions.

pred.pt <- data.frame(T.br = seq(min(df$T.br), max(df$T.br), length.out = 100)) # Dataframe to collect quantile info in

quant.PT.95 <- rq(N.comp ~ poly(T.br, 2), data = df.final[df.final$N.comp<10,], tau = 0.95) 
pred.pt$N.comp.95 <- predict(quant.PT.95, newdata = pred.pt)

quant.NT.75 <- rq(N.comp ~ poly(T.br, 2), data = df.final[df.final$N.comp<10,], tau = 0.75) 
pred.pt$N.comp.75 <- predict(quant.NT.75, newdata = pred.pt)

df.final$evol.bin <- factor(df.final$evol.bin, levels = c("ancestral", "other", "light"))

NT.qrs <- ggplot(df.final[df.final$N.comp<10,], aes(x = T.br, y = N.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.pt, aes(x = T.br, y = N.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.pt, aes(x = T.br, y = N.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Thermal breadth (°C)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "G") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Phosphorous")) +  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff # theme stuff
  )

NT.qrs  # Display the plot

