# Jason R Laurich

# February 12th, 2026

# Here I am going to going to run PCAs and RDAs on the experimental evolution Chlamydomonas reinhardtii data, and output some supplemental figures. 

# Inputs: 27_summary_table.csv
# Outputs: in figures-supplemental : 01_fig_s1_model_fitting.jpeg

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(vegan)  # For PCA and RDA

# Load & examine the data -------------------------------------------------

df <- read.csv("processed-data/27_summary_table.csv") # Summary file
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
                     labels = c("Ancestor 2", 
                                "Ancestor 3", 
                                "Ancestor 4", 
                                "Ancestor 5", 
                                "Mixed ancestry"))

# Run the PCA and RDAs ----------------------------------------------------

df.pca <- df %>% select(T.br, T.µ.max, I.comp, I.µ.max, N.comp, N.µ.max, P.comp, P.µ.max, S.µ.max, 
                        S.c, chl.a, chl.b, luthein, mean.N.µg.l, mean.P.µg.l, bio.vol, Evol.plt) # Prepare the data: selecting only the relevant columns
df.pca

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
                          "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "Evolution")

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

loadings$PC1 <- loadings$PC1 * 1.5
loadings$PC2 <- loadings$PC2 * 1.5

loadings$variable <- rownames(loadings) # Add variable names for annotation

loadings$metric <- factor(loadings$variable, 
                          levels = c("T.br", "T.µ.max","I.comp", "I.µ.max", "N.comp", "N.µ.max",
                                     "P.comp", "P.µ.max", "S.µ.max", "S.c", "chl.a", "chl.b", "luthein",
                                     "mean.N.µg.l", "mean.P.µg.l", "bio.vol"),
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
                                     "Lutein", 
                                     "N~content",
                                     "P~content",
                                     "Biovolume"))

adjust.x <- c(-0.05, 0.45, 0.35, 0.45, 0.46, 0.1, 0.45, 0.80, 0.50, 0.0, 0.60, 0.70, 0.60, 0.1, 0.4, -0.05) # adjustments for each label
adjust.y <- c(0.06, 0.00, 0.18, 0.25, 0.18, -0.1, 0.15, 0.25, 0.22, -0.02, 0.25, 0.02, -0.08, -0.02, -0.02, 0.06)

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
  
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text  = element_text(size = 12, face = "plain"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold")) 

pca_plot_arrows

ggsave("figures-supplemental/04_figs4_PCA.jpeg", pca_plot_arrows, width = 12, height = 9)
