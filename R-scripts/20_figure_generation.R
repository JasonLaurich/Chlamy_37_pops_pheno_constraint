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

# Load and examine data ---------------------------------------------------

df<-read.csv("data-processed/14_summary_metric_table.csv") # Summary file

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

# Figure 1: PCAs and RDAs -----------------------------------------------------------

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
                                       "Salt~tolerance",
                                       "Chlorophyll~italic(a)",
                                       "Chlorophyll~italic(b)",
                                       "Luthein", 
                                       "N~content",
                                       "P~content"))

adjust.x <- c(0.55, -0.05 , 0.36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) # adjustments for each label
adjust.y <- c(-0.1, 0.15, 0.11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

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
            vjust = 1, hjust = 1, color = "black", size = 5, parse = T) +
  theme(legend.position= c(0.15, 0.75),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 12, face = "plain"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"))

pca_plot_arrows

ggsave("figures/16_fig_1a_PCA.pdf", pca_plot_arrows, width = 12, height = 9)

# OK let's move onto the RDAs

# Figure 2: Intra-gradient trade-offs -------------------------------------

# We'll need our full dataset for the start

# Then we will bring in inter-specific datasets and plot the position of their metrics on our plots

# Figure 3: Inter-gradient trade-offs -------------------------------------

# We'll need our full dataset again

# Then we will bring in inter-specific datasets? Maybe?

