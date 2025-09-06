# Jason R Laurich

# September 5th, 2025

# Creating Figure 2: PCA and RDA analyses (RDA to be included in the supplement)

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(cowplot)
library(vegan)  # For PCA and RDA

# Load & examine the data -------------------------------------------------

df <- read.csv("data-processed/20_summary_table.csv") # Summary file
str(df)
head(df)

df$Evol.plt <- factor(df$Evol, 
                      levels = c("none", "L", "N", "P", "S", "B", "BS", "C"),
                      labels = c("Ancestral",
                                 "Light limitation",
                                 "Nitrogen limitation", 
                                 "Phosphorous limitation", 
                                 "Salt stress", 
                                 "Biotic depletion", 
                                 "Biotic depletion x Salt", 
                                 "Control"))

df$Anc.plt <- factor(df$Anc, 
                     levels = c("anc2", "anc3", "anc4", "anc5", "cc1690"),
                     labels = c("Population 2", 
                                "Population 3", 
                                "Population 4", 
                                "Population 5", 
                                "Mixed population"))

# Run the PCA and RDAs ----------------------------------------------------

df.pca <- df %>% select(T.br.0.56, T.µ.max, I.comp.0.56, I.µ.max, N.comp.0.56, N.µ.max, P.comp.0.56, P.µ.max, S.µ.max, 
                              S.c, chl.a, chl.b, luthein, mean.N.µg.l, mean.P.µg.l, Evol.plt) # Prepare the data: selecting only the relevant columns
df.pca

df.pca <- df.pca[-31,] # For now removing wonky/missing points (I.comp is way too high in row 31)

evol.pca <- df.pca$Evol.plt

df.pca <- df.pca %>% select(-Evol.plt)

pca.result <- prcomp(df.pca, center = TRUE, scale. = TRUE) # Perform PCA
pca.result

sdev <- pca.result$sdev  

sdev %>%
  { .^2 / sum(.^2) } %>%
  print() 

df.pca.res <- data.frame(pca.result$x, evol.pca)  # Add grouping factor
colnames(df.pca.res) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", 
                          "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "Evolution")

PCA <- ggplot(df.pca.res, aes(x = PC1, y = PC2, color = evol.pca)) +  # PCA biplot visualization
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
                          levels = c("T.br.0.56", "T.µ.max","I.comp.0.56", "I.µ.max", "N.comp.0.56", "N.µ.max",
                                     "P.comp.0.56", "P.µ.max", "S.µ.max", "S.c", "chl.a", "chl.b", "luthein",
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

adjust.x <- c(1.2, 0.4, 0.05, 0.67, 0.25, 0.75, 0.25, -0.05, -0.05, 0.8, 0.95, 0.45, 0.55, 0.67, 0.67) # adjustments for each label
adjust.y <- c(0.05, 0, -0.02, 0.11, -0.01, 0.15, -0.01, 0.1, 0.1, -0.02, 0.05, 0.2, -0.05, -0.05, -0.05)

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
  ) + 

  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size= 1.1) + # Add arrows for variable contributions
  scale_x_continuous(breaks = seq(-5, 5, by = 1)) +
  scale_y_continuous(breaks = seq(-5, 7, by = 1)) +
  
  geom_text(data = loadings, aes(x = var.x, y = var.y, label = metric),
            vjust = 1, hjust = 1,  color = "black", size = 5, parse = T) + # Add variable names to the plot
  
  theme(legend.position= c(0.25, 0.25),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 12, face = "plain"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")) 

pca_plot_arrows

ggsave("figures/09_fig2_PCA.jpeg", pca_plot_arrows, width = 12, height = 9)

###### Evolutionary history RDA #######

response_vars <- df.pca
explanatory_vars <- model.matrix(~ evol.pca)[, -1]  # Remove intercept

rda_result_evol <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result_evol)
sum(summary(rda_result_evol)$cont$importance[2, 1:rda_result_evol$CCA$rank]) # with P and N, evolutionary environment explains 20.27051% of the variation

rda_var_explained <- summary(rda_result_evol)$cont$importance["Proportion Explained", 1:2] * 100

rda_sites_evol <- as.data.frame(scores(rda_result_evol, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_evol <- as.data.frame(scores(rda_result_evol, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_evol <- as.data.frame(scores(rda_result_evol, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_evol$Evolution <- evol.pca # Add the evolutionary treatment labels to the site scores

rda_constraints_evol$label <- rownames(rda_constraints_evol) # Assign readable labels

rda_species_evol$metric <- factor(rownames(rda_species_evol), 
                                  levels = c("T.br.0.56", "T.µ.max","I.comp.0.56", "I.µ.max", "N.comp.0.56", "N.µ.max",
                                             "P.comp.0.56", "P.µ.max", "S.µ.max", "S.c", "chl.a", "chl.b", "luthein",
                                             "mean.N.µg.l", "mean.P.µg.l"),
                                  labels = c(" ", 
                                             " ", 
                                             " ", 
                                             " ", 
                                             " ", 
                                             " ",
                                             " ", 
                                             " ", 
                                             " ", 
                                             " ",
                                             " ",
                                             " ",
                                             " ", 
                                             "N~content",
                                             "P~content")) # Only N and P content can really be made out on this plot

adjust.x.rda1 <- c(0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -5.2, -2.9) # adjustments for each label
adjust.y.rda1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.75, -0.2)

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
  ) +
  
  geom_segment(data = rda_species_evol, aes(x = 0, y = 0, xend = RDA1*2, yend = RDA2*2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size= 1.1) +
  scale_x_continuous(breaks = seq(-20, 25, by = 5)) +
  scale_y_continuous(breaks = seq(-20, 40, by = 5)) + # Add arrows for variable contributions
  
  
  geom_text(data = rda_species_evol, aes(x = var.x, y = var.y, label = metric),
            vjust = 1, hjust = 1, color = "black", size = 5, parse = T) + # Add variable names to the plot
  
  theme(legend.position= c(0.15, 0.45),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 12, face = "plain"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.03)) 

rda_evol_plot_arrows # N and P content is skewing everything. 

df.pca1 <- df.pca %>% 
  select(-mean.N.µg.l, -mean.P.µg.l)

response_vars <- df.pca1
explanatory_vars <- model.matrix(~ evol.pca)[, -1]  # Remove intercept

rda_result_evol <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result_evol) 
sum(summary(rda_result_evol)$cont$importance[2, 1:rda_result_evol$CCA$rank]) # without P and N, evolutionary environment explains 20.23622% of the variation

rda_sites_evol <- as.data.frame(scores(rda_result_evol, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_evol <- as.data.frame(scores(rda_result_evol, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_evol <- as.data.frame(scores(rda_result_evol, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_evol$Evolution <- evol.pca # Add the evolutionary treatment labels to the site scores

rda_constraints_evol$label <- rownames(rda_constraints_evol) # Assign readable labels for centroids

rda_species_evol$metric <- factor(rownames(rda_species_evol), 
                                  levels = c("T.br.0.56", "T.µ.max","I.comp.0.56", "I.µ.max", "N.comp.0.56", "N.µ.max",
                                             "P.comp.0.56", "P.µ.max", "S.µ.max", "S.c", "chl.a", "chl.b", "luthein"),
                                  labels = c(" ", 
                                             "mu~max~(T)", 
                                             " ", 
                                             " ", 
                                             " ", 
                                             " ",
                                             "1/P^\"*\"", 
                                             " ", 
                                             " ", 
                                             "Salt~tolerance",
                                             "Chlorophyll~italic(a)",
                                             "Chlorophyll~italic(b)",
                                             " "))

adjust.x.rda1 <- c(1.15, -0.2, 0, 0, 0, 0, 0.1, 0, 0, 1.2, 2.75, 3.75, 0.7) # adjustments for each label
adjust.y.rda1 <- c(0, 0, 0, 0, 0, 0, -0.4, 0, 0, -0.85, -0.15, 0.5, 0.25)

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

df.pca2 <- df %>% select(T.br.0.56, T.µ.max, I.comp.0.56, I.µ.max, N.comp.0.56, N.µ.max, P.comp.0.56, P.µ.max, S.µ.max, 
                                S.c, chl.a, chl.b, luthein, mean.N.µg.l, mean.P.µg.l, Anc.plt) # Prepare the data: selecting only the relevant columns
df.pca2

df.pca2 <- df.pca2[-31,] # For now removing wonky/missing points (I.comp is way too high in row 31)

anc.pca <- df.pca2$Anc.plt

df.pca2 <- df.pca2 %>% select(-Anc.plt)

response_vars <- df.pca2

explanatory_vars <- model.matrix(~ anc.pca)[, -1]  # Remove intercept

rda_result_anc <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result_anc)
sum(summary(rda_result_anc)$cont$importance[2, 1:rda_result_anc$CCA$rank]) # with P and N, ancestry explains 12.96415% of the variation

rda_sites_anc <- as.data.frame(scores(rda_result_anc, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_anc <- as.data.frame(scores(rda_result_anc, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_anc <- as.data.frame(scores(rda_result_anc, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_anc$Ancestry <- anc.pca # Add the evolutionary treatment labels to the site scores

rda_constraints_anc$label <- rownames(rda_constraints_anc) # Assign readable labels for centroids

rda_species_anc$metric <- factor(rownames(rda_species_anc), 
                                 levels = c("T.br.0.56", "T.µ.max","I.comp.0.56", "I.µ.max", "N.comp.0.56", "N.µ.max",
                                            "P.comp.0.56", "P.µ.max", "S.µ.max", "S.c", "chl.a", "chl.b", "luthein",
                                            "mean.N.µg.l", "mean.P.µg.l"),
                                 labels = c(" ", 
                                            " ", 
                                            " ", 
                                            " ", 
                                            " ", 
                                            " ",
                                            " ", 
                                            " ", 
                                            " ", 
                                            " ",
                                            " ",
                                            " ",
                                            " ", 
                                            "N~content",
                                            "P~content"))

adjust.x.rda2 <- c(0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20.5, 6.5) # adjustments for each label
adjust.y.rda2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, -1.5)

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

explanatory_vars <- model.matrix(~ anc.pca)[, -1]  # Remove intercept

rda_result_anc <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result_anc)
sum(summary(rda_result_anc)$cont$importance[2, 1:rda_result_anc$CCA$rank]) # without P and N, ancestry explains 3.322815% of the variation

rda_sites_anc <- as.data.frame(scores(rda_result_anc, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species_anc <- as.data.frame(scores(rda_result_anc, display = "species")) # Extract RDA species (trait arrows)

rda_constraints_anc <- as.data.frame(scores(rda_result_anc, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites_anc$Ancestry <- anc.pca # Add the evolutionary treatment labels to the site scores

rda_constraints_anc$label <- rownames(rda_constraints_anc) # Assign readable labels for centroids

rda_species_anc$metric <- factor(rownames(rda_species_anc), 
                                 levels = c("T.br.0.56", "T.µ.max","I.comp.0.56", "I.µ.max", "N.comp.0.56", "N.µ.max",
                                             "P.comp.0.56", "P.µ.max", "S.µ.max", "S.c", "chl.a", "chl.b", "luthein"),
                                 labels = c("Thermal~breadth", 
                                            " ", 
                                            " ", 
                                            " ", 
                                            " ", 
                                            " ",
                                            " ", 
                                            " ", 
                                            " ", 
                                            "Salt~tolerance",
                                            "Chlorophyll~italic(a)",
                                            "Chlorophyll~italic(b)",
                                            " "))

adjust.x.rda2 <- c(7.85, 0 , 0, 0, 0, 0, 0.5, 0, 0, -0.2, 6.85, 6.65, 0) # adjustments for each label
adjust.y.rda2 <- c(1.25, 0, 0, 0, 0, 0, -0.1, 0, 0.5, 0, -0.15, 0.45, 0)

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

# Compile the supplemental figure -----------------------------------------

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

ggsave("figures/10_supp_fig2_rdas.jpeg", rdas, width = 16, height = 8) # PDF was rendering weird
