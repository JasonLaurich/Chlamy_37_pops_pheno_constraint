# Jason R Laurich
# April 14, 2025

# We are going to generate "final" figures for the project here

# Load packages -----------------------------------------------------------

library(cowplot)
library(tidyverse)
library(R2jags) # Fits Bayesian models
library(mcmcplots) # Diagnostic plots for fits
library(gridExtra)
library(vegan)  # For PCA and RDA
library(ggrepel)
library(quantreg)

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

###### Temperature ######

df.final$shape <- ifelse(df.final$evol == "none", 22, 16) # I want to add a shape column to the dataframe that I will update
# The idea is to label un-evolved populations with a plus, and then later (not for T) relevant experimental evolution nutrient conditions with a star

df.final$evol.bin <- ifelse(df.final$evol == "none", "ancestral", "evolved") # For regressions

par.res.T <- par_frt(df, xvar = "r.max_T", yvar = "T.br")

T_par <- ggplot(df.final, aes(x = r.max_T, y = T.br)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Thermal breadth", 
       title = "Thermal performance") +
  scale_shape_manual(values = c(16, 3)) +  # Keep custom shapes
  geom_line(data = par.res.T, aes(x = r.max_T, y = T.br), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

T_par # Raw pareto front.

pred.t <- data.frame(r.max_T = seq(min(df$r.max_T), max(df$r.max_T), length.out = 100)) # Dataframe to collect quantile info in

quant.T.95 <- rq(T.br ~ poly(r.max_T, 2), data = df, tau = 0.95) 
pred.t$T.br.95 <- predict(quant.T.95, newdata = pred.t)

quant.T.75 <- rq(T.br ~ poly(r.max_T, 2), data = df, tau = 0.75) 
pred.t$T.br.75 <- predict(quant.T.75, newdata = pred.t)

T.qrs <- ggplot(df.final, aes(x = r.max_T, y = T.br, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.t, aes(x = r.max_T, y = T.br.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.t, aes(x = r.max_T, y = T.br.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       color = "Evolutionary History",
       title = "E") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1"), 
                     labels = c("Ancestral", "Evolved")) +  
  theme_classic() +
  theme(
    legend.position = c(0.8, 0.85),  # Move legend inside the plot
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 12, face = "plain"),  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.qrs  # Display the plot

###### Light ######

df.final$shape <- ifelse(df.final$evol == "none", 22, 
                   ifelse(df.final$evol == "L", 8, 16)) # Ls are now equivalent to 8, for later mapping

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                      ifelse(df.final$evol == "L", 'light', 'other')) # for testing regressions.

par.res.L <- par_frt(df.final[df.final$I.comp<10,], xvar = "r.max_I", yvar = "I.comp")

L_par <- ggplot(df.final[df.final$I.comp<10,], aes(x = r.max_I, y = I.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Light limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res.L, aes(x = r.max_I, y = I.comp), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

L_par # Raw pareto front.

pred.l <- data.frame(r.max_I = seq(min(df$r.max_I), max(df$r.max_I), length.out = 100)) # Dataframe to collect quantile info in

quant.I.95 <- rq(I.comp ~ poly(r.max_I, 2), data = df.final[df.final$I.comp<10,], tau = 0.95) 
pred.l$I.comp.95 <- predict(quant.I.95, newdata = pred.l)

pred.l <- pred.l %>%
  mutate(I.comp.95 = ifelse(row_number() <= which.max(pred.l$I.comp.95) & I.comp.95 < max(I.comp.95), NA, I.comp.95))

quant.I.75 <- rq(I.comp ~ poly(r.max_I, 2), data = df.final[df.final$I.comp<10,], tau = 0.75) 
pred.l$I.comp.75 <- predict(quant.I.75, newdata = pred.l)

df.final$evol.bin <- factor(df.final$evol.bin, levels = c("ancestral", "other", "light"))

L.qrs <- ggplot(df.final[df.final$I.comp<10,], aes(x = r.max_I, y = I.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.l, aes(x = r.max_I, y = I.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.l, aes(x = r.max_I, y = I.comp.75), color = "black", size = 1.2, linetype = "dashed") +  

  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/R*)", 
       color = "Evolutionary History",
       title = "A") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Light")) +  
  theme_classic() +
  theme(
    legend.position = c(0.25, 0.85),  # Move legend inside the plot
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 12, face = "plain"),  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff # theme stuff
  )

L.qrs  # Display the plot

find_nearest_index <- function(x, ref_vec) { # Write a function to find the closest point to the actual r.max value in the pred data frame. 
  which.min(abs(ref_vec - x))
}

l.75 <- df.final %>% # Now calculate the number of light points above that. 
  filter(evol.bin == "light", I.comp < 10) %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(r.max_I, pred.l$r.max_I),
    I.comp.75.pred = pred.l$I.comp.75[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(I.comp > I.comp.75.pred) %>%
  nrow()

set.seed(123)  # for reproducibility
n_iter <- 1000
null_counts_l <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.final %>%
    filter(I.comp < 10) %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts_l[i] <- df.shuff %>%
    filter(evol.shuff == "light") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.l$r.max_I - r.max_I)),
      I.comp.75.pred = pred.l$I.comp.75[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(I.comp > I.comp.75.pred) %>%
    nrow()
}

p_val_l <- mean(null_counts_l >= l.75) # 0.726

###### Nitrogen ######

df.final$shape <- ifelse(df.final$evol == "none", 22, 
                   ifelse(df.final$evol == "N", 8, 16)) # Ns are now equivalent to 8, for later mapping

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                      ifelse(df.final$evol == "N", 'nit', 'other')) # for testing regressions.

par.res.N <- par_frt(df.final, xvar = "r.max_N", yvar = "N.comp")

N_par <- ggplot(df.final, aes(x = r.max_N, y = N.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Nitrogen limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res.N, aes(x = r.max_N, y = N.comp), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

N_par # Raw pareto front.

pred.n <- data.frame(r.max_N = seq(min(df.final$r.max_N), max(df.final$r.max_N), length.out = 100)) # Dataframe to collect quantile info in

quant.N.95 <- rq(N.comp ~ poly(r.max_N, 2), data = df.final, tau = 0.95) 
pred.n$N.comp.95 <- predict(quant.N.95, newdata = pred.n)

quant.N.75 <- rq(N.comp ~ poly(r.max_N, 2), data = df.final, tau = 0.75) 
pred.n$N.comp.75 <- predict(quant.N.75, newdata = pred.n)

df.final$evol.bin <- factor(df.final$evol.bin, levels = c("ancestral", "other", "nit"))

N.qrs <- ggplot(df.final, aes(x = r.max_N, y = N.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.n, aes(x = r.max_N, y = N.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.n, aes(x = r.max_N, y = N.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/R*)", 
       color = "Evolutionary History",
       title = "B") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Nitrogen")) +  
  theme_classic() +
  theme(
    legend.position = c(0.35, 0.85),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

N.qrs  # Display the plot

n.75 <- df.final %>% # Now calculate the number of light points above that. 
  filter(evol.bin == "nit") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(r.max_N, pred.n$r.max_N),
    N.comp.75.pred = pred.n$N.comp.75[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(N.comp > N.comp.75.pred) %>%
  nrow()

null_counts_n <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.final %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts_n[i] <- df.shuff %>%
    filter(evol.shuff == "nit") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.n$r.max_N - r.max_N)),
      N.comp.75.pred = pred.n$N.comp.75[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(N.comp > N.comp.75.pred) %>%
    nrow()
}

p_val_n <- mean(null_counts_n >= n.75) # 0.354

###### Phosphorous ######

df.final$shape <- ifelse(df.final$evol == "none", 22, 
                         ifelse(df.final$evol == "P", 8, 16)) # Ps are now equivalent to 8, for later mapping

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                            ifelse(df.final$evol == "P", 'phos', 'other')) # for testing regressions.

par.res.P <- par_frt(df.final, xvar = "r.max_P", yvar = "P.comp")

P_par <- ggplot(df.final, aes(x = r.max_P, y = P.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res.P, aes(x = r.max_P, y = P.comp), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_par # Raw pareto front.

pred.p <- data.frame(r.max_P = seq(min(df.final$r.max_P), max(df.final$r.max_P), length.out = 100)) # Dataframe to collect quantile info in

quant.P.95 <- rq(P.comp ~ poly(r.max_P, 2), data = df.final, tau = 0.95) 
pred.p$P.comp.95 <- predict(quant.P.95, newdata = pred.p)

quant.P.75 <- rq(P.comp ~ poly(r.max_P, 2), data = df.final, tau = 0.75) 
pred.p$P.comp.75 <- predict(quant.P.75, newdata = pred.p)

P.qrs <- ggplot(df.final, aes(x = r.max_P, y = P.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.p, aes(x = r.max_P, y = P.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.p, aes(x = r.max_P, y = P.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/R*)", 
       color = "Evolutionary History",
       title = "C") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Phosphorous")) +  
  theme_classic() +
  theme(
    legend.position = c(0.35, 0.85),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

P.qrs  # Display the plot

p.75 <- df.final %>% # Now calculate the number of light points above that. 
  filter(evol.bin == "phos") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(r.max_P, pred.p$r.max_P),
    P.comp.75.pred = pred.p$P.comp.75[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(P.comp > P.comp.75.pred) %>%
  nrow()

null_counts_p <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.final %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts_p[i] <- df.shuff %>%
    filter(evol.shuff == "phos") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.p$r.max_P - r.max_P)),
      P.comp.75.pred = pred.p$P.comp.75[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(P.comp > P.comp.75.pred) %>%
    nrow()
}

p_val_p <- mean(null_counts_p >= p.75) # 0.116

###### Salt ######

df.final$shape <- ifelse(df.final$evol == "none", 22, 
                   ifelse(df.final$evol %in% c("S", "BS"), 8, 16)) # Ss and BSs are now equivalent to 8, for later mapping

df.final$evol.bin <- ifelse(df.final$evol == "none", 'ancestral', 
                      ifelse(df.final$evol %in% c("S", "BS"), 'salt', 'other')) # for testing regressions.

par.res.S <- par_frt(df.final, xvar = "r.max_S", yvar = "S.c.mod")

S_par <- ggplot(df.final, aes(x = r.max_S, y = S.c.mod)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Salt stress") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res.S, aes(x = r.max_S, y =S.c.mod), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

S_par # Raw pareto front.

pred.s <- data.frame(r.max_S = seq(min(df.final$r.max_S), max(df.final$r.max_S), length.out = 100)) # Dataframe to collect quantile info in

quant.S.95 <- rq(S.c.mod ~ poly(r.max_S, 2), data = df.final, tau = 0.95) 
pred.s$S.comp.95 <- predict(quant.S.95, newdata = pred.s)

quant.S.75 <- rq(S.c.mod ~ poly(r.max_S, 2), data = df.final, tau = 0.75) 
pred.s$S.comp.75 <- predict(quant.S.75, newdata = pred.s)

S.qrs <- ggplot(df.final, aes(x = r.max_S, y = S.c.mod, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 3) +  # Scatter plot of raw data
  
  geom_line(data = pred.s, aes(x = r.max_S, y = S.comp.95), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.s, aes(x = r.max_S, y = S.comp.75), color = "black", size = 1.2, linetype = "dashed") +  
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Salt tolerance (c)", 
       color = "Evolutionary History",
       title = "D") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Salt or Biotic x Salt")) +  
  theme_classic() +
  theme(
    legend.position = c(0.35, 0.85),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 12, face = "plain"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

S.qrs  # Display the plot

s.75 <- df.final %>% # Now calculate the number of light points above that. 
  filter(evol.bin == 'salt') %>% 
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(r.max_S, pred.s$r.max_S),
    S.comp.75.pred = pred.s$S.comp.75[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(S.c.mod > S.comp.75.pred) %>%
  nrow()

null_counts_s <- numeric(n_iter)

for (i in 1:n_iter) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  df.shuff <- df.final %>%
    mutate(evol.shuff = sample(evol.bin, rep = F))
  
  null_counts_s[i] <- df.shuff %>%
    filter(evol.shuff == "salt") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.s$r.max_S - r.max_S)),
      S.comp.75.pred = pred.s$S.comp.75[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(S.c.mod > S.comp.75.pred) %>%
    nrow()
}

p_val_s <- mean(null_counts_s >= s.75) # 0.001!

# Figure 1 : all plots together

plots <- list(L.qrs, N.qrs, P.qrs, S.qrs, T.qrs)

plots_nolegend <- lapply(plots, function(p) p + theme(legend.position = "none"))

legend_df <- data.frame(
  x = c(1, 2, 1, 2),
  y = c(1, 1, 2, 2),
  Group = factor(c("Ancestral", "Other", "Matching", "Matching")),
  LineType = factor(c("QR 95%", "QR 75%", "QR 95%", "QR 75%"))
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y)) +
  geom_point(aes(color = Group), size = 3) +
  geom_line(aes(linetype = LineType), size = 1) +
  scale_color_manual(values = c("Ancestral" = "black", 
                                "Other" = "gold", 
                                "Matching" = "orchid")) +
  scale_linetype_manual(values = c("QR 95%" = "solid", "QR 75%" = "dashed"),
                        labels = c("95th", "75th"),
                        name = "Quantile regression") +
  labs(linetype = "Quantile regression", color = "Evolutionary context") +
  theme_void() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "lines")
  )


legend_only <- get_legend(legend_plot)

all_plots <- c(plots_nolegend, list(legend_only))

grad_toffs <- plot_grid(plotlist = all_plots,
          ncol = 2,
          align = "hv")

ggsave("figures/18_fig_3_intra-gradient_tradeoffs.jpeg", grad_toffs, width = 8, height = 12) # Still looks off?.

# Then we will bring in inter-specific datasets and plot the position of their metrics on our plots

# Figure 4: Intra-gradient, interspecific trade-offs-------------------------------------

# Bringing in interspecific data sets

###### Temperature ######

# Figure 5: Inter-gradient trade-offs -------------------------------------



