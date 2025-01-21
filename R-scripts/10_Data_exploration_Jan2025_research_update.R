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
library(gridExtra)
library(cowplot)
library(ggpubr)
library(vegan)  # For PCA and RDA
library(ggrepel)

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

T.I.evol.noleg <- T.I + 
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
    legend.position = 'none',  # Move legend to bottom right
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )


T.I.evol.noleg

par.res <- par_frt(df, xvar = "N.comp", yvar = "T.br")

T.N <- ggplot(df, aes(x = N.comp, y = T.br)) +  # removing the outlier data point
  geom_point(size = 2) +  # Original data points
  geom_line(data = par.res, aes(x = N.comp, y = T.br), color = "red", size = 1) +  # Pareto frontier line
  labs(x = "Competitive ability (1/N*)", y = "Thermal breadth (T.br)") +
  theme_classic()

T.N

T.N.evol <- T.N + 
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
    legend.position = 'none',  # no legend
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )


T.N.evol

par.res <- par_frt(df, xvar = "P.comp", yvar = "T.br")

T.P <- ggplot(df, aes(x = P.comp, y = T.br)) +  # removing the outlier data point
  geom_point(size = 2) +  # Original data points
  geom_line(data = par.res, aes(x = P.comp, y = T.br), color = "red", size = 1) +  # Pareto frontier line
  labs(x = "Competitive ability (1/P*)", y = "Thermal breadth (T.br)") +
  theme_classic()

T.P

T.P.evol <- T.P + 
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
    legend.position = 'none',  # no legend
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )


T.P.evol

par.res <- par_frt(df, xvar = "S.c.mod", yvar = "T.br")

T.S <- ggplot(df, aes(x = S.c.mod, y = T.br)) +  # removing the outlier data point
  geom_point(size = 2) +  # Original data points
  geom_line(data = par.res, aes(x = S.c.mod, y = T.br), color = "red", size = 1) +  # Pareto frontier line
  labs(x = "Salt tolerance (c)", y = "Thermal breadth (T.br)") +
  theme_classic()

T.S

T.S.evol <- T.S + 
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
    legend.position = 'none',  # no legend
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )


T.S.evol

par.res <- par_frt(df[df$I.comp < 10,], xvar = "N.comp", yvar = "I.comp")

I.N <- ggplot(df[df$I.comp < 10,], aes(x = N.comp, y = I.comp)) +  # removing the outlier data point
  geom_point(size = 2) +  # Original data points
  geom_line(data = par.res, aes(x = N.comp, y = I.comp), color = "red", size = 1) +  # Pareto frontier line
  labs(x = "Competitive ability (1/N*)", y = "Competitive ability (1/I*)") +
  theme_classic()

I.N

I.N.evol <- I.N + 
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
    legend.position = 'none',  # no legend
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )


I.N.evol

par.res <- par_frt(df[df$I.comp < 10,], xvar = "P.comp", yvar = "I.comp")

I.P <- ggplot(df[df$I.comp < 10,], aes(x = P.comp, y = I.comp)) +  # removing the outlier data point
  geom_point(size = 2) +  # Original data points
  geom_line(data = par.res, aes(x = P.comp, y = I.comp), color = "red", size = 1) +  # Pareto frontier line
  labs(x = "Competitive ability (1/P*)", y = "Competitive ability (1/I*)") +
  theme_classic()

I.P

I.P.evol <- I.P + 
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
    legend.position = 'none',  # no legend
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )


I.P.evol

par.res <- par_frt(df[df$I.comp < 10,], xvar = "S.c.mod", yvar = "I.comp")

I.S <- ggplot(df[df$I.comp < 10,], aes(x = S.c.mod, y = I.comp)) +  # removing the outlier data point
  geom_point(size = 2) +  # Original data points
  geom_line(data = par.res, aes(x = S.c.mod, y = I.comp), color = "red", size = 1) +  # Pareto frontier line
  labs(x = "Salt tolerance (c)", y = "Competitive ability (1/I*)") +
  theme_classic()

I.S

I.S.evol <- I.S + 
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
    legend.position = 'none',  # no legend
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )


I.S.evol

par.res <- par_frt(df, xvar = "P.comp", yvar = "N.comp")

N.P <- ggplot(df, aes(x = P.comp, y = N.comp)) +  # removing the outlier data point
  geom_point(size = 2) +  # Original data points
  geom_line(data = par.res, aes(x = P.comp, y = N.comp), color = "red", size = 1) +  # Pareto frontier line
  labs(x = "Competitive ability (1/P*)", y = "Competitive ability (1/N*)") +
  theme_classic()

N.P

N.P.evol <- N.P + 
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
    legend.position = 'none',  # no legend
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )


N.P.evol

par.res <- par_frt(df, xvar = "S.c.mod", yvar = "N.comp")

N.S <- ggplot(df, aes(x = S.c.mod, y = N.comp)) +  # removing the outlier data point
  geom_point(size = 2) +  # Original data points
  geom_line(data = par.res, aes(x = S.c.mod, y = N.comp), color = "red", size = 1) +  # Pareto frontier line
  labs(x = "Salt tolerance (c)", y = "Competitive ability (1/N*)") +
  theme_classic()

N.S

N.S.evol <- N.S + 
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
    legend.position = 'none',  # no legend
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )


N.S.evol

par.res <- par_frt(df[df$P.comp < 4,], xvar = "S.c.mod", yvar = "P.comp")

P.S <- ggplot(df[df$P.comp < 4,], aes(x = S.c.mod, y = P.comp)) +  # removing the outlier data point
  geom_point(size = 2) +  # Original data points
  geom_line(data = par.res, aes(x = S.c.mod, y = P.comp), color = "red", size = 1) +  # Pareto frontier line
  labs(x = "Salt tolerance (c)", y = "Competitive ability (1/P*)") +
  theme_classic()

P.S

P.S.evol <- P.S + 
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
    legend.position = 'none',  # no legend
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )


P.S.evol

T.I.evol.leg <- T.I + 
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
    text = element_text(face = "bold"),  # Make all text bold
    axis.title = element_text(face = "bold"),  # Bold axis titles
    axis.text = element_text(face = "bold"),  # Bold axis tick labels
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.text = element_text(face = "bold")   # Bold legend labels
  )

T.I.evol.leg

leg <- get_legend(T.I.evol.leg)

as_ggplot(leg)

p_list <- list(T.I.evol.noleg, T.N.evol, T.P.evol, T.S.evol, # OK let's plot them together
               NULL, I.N.evol, I.P.evol, I.S.evol,
               NULL, NULL, N.P.evol, N.S.evol,
               NULL, NULL, NULL, P.S.evol)

empty_plot <- ggplot() + theme_void() # Convert NULL values to empty plots to maintain the layout

p_list <- lapply(p_list, function(x) if(is.null(x)) empty_plot else x) # Replace NULL values with empty_plot to avoid errors in grid arrangement

pareto.plot <- grid.arrange(grobs = p_list, ncol = 4, nrow = 4) # Arrange plots in a 4x4 grid

ggsave("figures/12_pareto_grid_plot.pdf", pareto.plot, width = 20, height = 16)

# Let's do a PCA

df.pca <- df %>% select(T.br, I.comp, N.comp, P.comp, S.c.mod, evol.plt) # Prepare the data: selecting only the relevant columns
df.pca

df.pca <- df.pca[-c(31, 33:36),] # For now removing wonky/missing points

evol.fil <- df.pca$evol.plt

df.pca <- df.pca %>% select(-evol.plt)

pca.result <- prcomp(df.pca, center = TRUE, scale. = TRUE) # Perform PCA
pca.result

pca.df <- data.frame(pca.result$x, evol.fil)  # Add grouping factor
colnames(pca.df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "Evolution")

screeplot(pca.result, type = "lines") # Scree plot to check explained variance

PCA <- ggplot(pca.df, aes(x = PC1, y = PC2, color = Evolution)) +  # PCA biplot visualization
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

# Extract PCA loadings (rotation matrix)
loadings <- as.data.frame(pca.result$rotation)

# Scale loadings to fit within the PCA plot (adjust the scaling factor if needed)
loadings$PC1 <- loadings$PC1 * max(abs(pca.df$PC1))
loadings$PC2 <- loadings$PC2 * max(abs(pca.df$PC2))

# Add variable names for annotation
loadings$variable <- rownames(loadings)

# Create PCA plot with arrows
pca_plot_arrows <- ggplot(pca.df, aes(x = PC1, y = PC2, color = Evolution)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("PC1 (", round(pca.result$sdev[1]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pca.result$sdev[2]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       title = "PCA of Thermal Performance and Competitive Abilities") +
  scale_color_manual(values = evol_colors) +  # Use your custom colors
  # Add arrows for variable contributions
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  # Add variable names to the plot
  geom_text(data = loadings, aes(x = PC1, y = PC2, label = variable),
            vjust = 1, hjust = 1, color = "black", size = 5)

# Show the plot
pca_plot_arrows

ggsave("figures/13_PCA.pdf", pca_plot_arrows, width = 10, height = 8)

# Redundancy analysis

response_vars <- df.pca
explanatory_vars <- model.matrix(~ evol.fil)[, -1]  # Remove intercept

rda_result <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result)

# Extract RDA site scores (sample coordinates)
rda_sites <- as.data.frame(scores(rda_result, display = "sites"))

# Extract RDA species (trait arrows)
rda_species <- as.data.frame(scores(rda_result, display = "species"))

# Extract explanatory variable centroids (e.g., treatment centroids)
rda_constraints <- as.data.frame(scores(rda_result, display = "bp"))

# Add the evolutionary treatment labels to the site scores
rda_sites$Evolution <- evol.fil

# Prepare RDA site data (samples)
rda_sites <- rda_sites %>% 
  rename(RDA1 = RDA1, RDA2 = RDA2)

# Prepare RDA trait data (arrows)
rda_species <- rda_species %>%
  rename(RDA1 = RDA1, RDA2 = RDA2)

# Prepare centroids of treatment conditions
rda_constraints <- rda_constraints %>%
  rename(RDA1 = RDA1, RDA2 = RDA2)

# Assign readable labels for centroids
rda_constraints$label <- rownames(rda_constraints)

# Define custom colors for evolution environments
evol_colors <- c(
  "Biotic depletion" = "darkorange",
  "Biotic depletion x Salt" = "deepskyblue1",
  "Control" = "forestgreen",
  "Light limitation" = "gold",
  "Nitrogen limitation" = "magenta3",
  "Ancestral" = "black",
  "Phosphorous limitation" = "firebrick",  
  "Salt stress" = "blue"
)

# Create the RDA biplot
rda_plot <- ggplot() +
  geom_point(data = rda_sites, aes(x = RDA1, y = RDA2, color = Evolution), size = 3, alpha = 0.8) +
  
  # Retain trait arrows
  geom_segment(data = rda_species, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  
  # Retain trait labels only
  geom_text_repel(data = rda_species, aes(x = RDA1, y = RDA2, label = rownames(rda_species)), size = 5) +
  
  # Customize color legend
  scale_color_manual(name = "Evolution environment", values = evol_colors) +
  
  # Axes labels
  labs(x = paste("RDA1 (", round(summary(rda_result)$cont$importance[2, 1] * 100, 2), "%)", sep = ""),
       y = paste("RDA2 (", round(summary(rda_result)$cont$importance[2, 2] * 100, 2), "%)", sep = ""),
       title = "Redundancy Analysis (RDA) of Traits by Evolution Environment") +
  
  theme_classic() +
  theme(text = element_text(size = 14, face = "bold"),
        legend.position = "right",
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

rda_plot

ggsave("figures/14_RDA.pdf", rda_plot, width = 10, height = 8)

# Let's pull up my R2jags objects to look at some plots.

pops <- sample(1:37, 5, replace = FALSE)

file.names <- paste0("R2jags-objects/pop_", pops, "_lactin2.RData")

for (i in pops){
  load(paste0("R2jags-objects/pop_", i, "_lactin2.RData"))
  assign(paste0("tpc", i), lac_jag)
}

tpc2$BUGSoutput$summary[1:6,] # Get estimates

df.jags <- data.frame(tpc2$BUGSoutput$summary)
df.jags.plot <- df.jags[-c(1:6),]
df.jags.plot$temp <- seq(0, 45, 0.05)

lact.jag.plot <- ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gold", alpha = 0.5) +# Add shaded uncertainty region (LCL to UCL)
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  scale_x_continuous(limits = c(0, 45)) + 
  scale_y_continuous(limits = c(-2, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "TPC, Lactin II"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

df.jags2 <- data.frame(tpc9$BUGSoutput$summary[-c(1:6),])
df.jags2$temp <- seq(0, 45, 0.05)

df.jags3 <- data.frame(tpc18$BUGSoutput$summary[-c(1:6),])
df.jags3$temp <- seq(0, 45, 0.05)

df.jags4 <- data.frame(tpc26$BUGSoutput$summary[-c(1:6),])
df.jags4$temp <- seq(0, 45, 0.05)

df.jags5 <- data.frame(tpc30$BUGSoutput$summary[-c(1:6),])
df.jags5$temp <- seq(0, 45, 0.05)

lact.jag.plot2 <- ggplot(data = df.jags.plot, aes(x = temp)) +
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  
  geom_line(data = df.jags2, aes(x = temp, y = mean), color = "darkorange", size = 1) + 
  
  geom_line(data = df.jags3, aes(x = temp, y = mean), color = "forestgreen", size = 1) + 
  
  geom_line(data = df.jags4, aes(x = temp, y = mean), color = "goldenrod1", size = 1) + 
  
  geom_line(data = df.jags5, aes(x = temp, y = mean), color = "dodgerblue3", size = 1) + 
  
  scale_x_continuous(limits = c(0, 45)) + 
  scale_y_continuous(limits = c(-2, 5.5)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Temperature (", degree, "C)")),
    y = "Growth rate",
    title = "TPCs, Lactin II"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

lact.jag.plot2

for (i in pops){
  load(paste0("R2jags-objects/pop_", i, "_salt_tolerance.RData"))
  assign(paste0("salt", i), monod_jag)
}

df.jags1 <- data.frame(salt2$BUGSoutput$summary)[-c(1:4,2006),]   # generate the sequence of r.pred values
df.jags1$salt <- seq(0, 10, 0.005)

df.jags2 <- data.frame(salt9$BUGSoutput$summary)[-c(1:4,2006),]   # generate the sequence of r.pred values
df.jags2$salt <- seq(0, 10, 0.005)

df.jags3 <- data.frame(salt18$BUGSoutput$summary)[-c(1:4,2006),]   # generate the sequence of r.pred values
df.jags3$salt <- seq(0, 10, 0.005)

df.jags4 <- data.frame(salt26$BUGSoutput$summary)[-c(1:4,2006),]   # generate the sequence of r.pred values
df.jags4$salt <- seq(0, 10, 0.005)

df.jags5 <- data.frame(salt30$BUGSoutput$summary)[-c(1:4,2006),]   # generate the sequence of r.pred values
df.jags5$salt <- seq(0, 10, 0.005)

salt.jag.plot <- ggplot(data = df.jags1, aes(x = salt)) +
  
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  
  geom_line(data = df.jags2, aes(x = salt, y = mean), color = "darkorange", size = 1) + 
  
  geom_line(data = df.jags3, aes(x = salt, y = mean), color = "forestgreen", size = 1) + 
  
  geom_line(data = df.jags4, aes(x = salt, y = mean), color = "goldenrod1", size = 1) + 
  
  geom_line(data = df.jags5, aes(x = salt, y = mean), color = "dodgerblue3", size = 1) + 
  
  scale_x_continuous(limits = c(0, 10)) + 
  scale_y_continuous(limits = c(-0.25, 2.25)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Salt (concentration)")),
    y = "Growth rate",
    title = "Logistic growth curve, salt stress"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

salt.jag.plot

for (i in pops){
  load(paste0("R2jags-objects/pop_", i, "_phosphorous_monod.RData"))
  assign(paste0("phos", i), monod_jag)
}

df.jags1 <- data.frame(phos2$BUGSoutput$summary)[-c(1:3,2005),]   # generate the sequence of r.pred values
df.jags1$phos <- seq(0, 50, 0.025)

df.jags2 <- data.frame(phos9$BUGSoutput$summary)[-c(1:3,2005),]   # generate the sequence of r.pred values
df.jags2$phos <- seq(0, 50, 0.025)

df.jags3 <- data.frame(phos18$BUGSoutput$summary)[-c(1:3,2005),]   # generate the sequence of r.pred values
df.jags3$phos <- seq(0, 50, 0.025)

df.jags4 <- data.frame(phos26$BUGSoutput$summary)[-c(1:3,2005),]   # generate the sequence of r.pred values
df.jags4$phos <- seq(0, 50, 0.025)

df.jags5 <- data.frame(phos30$BUGSoutput$summary)[-c(1:3,2005),]   # generate the sequence of r.pred values
df.jags5$phos <- seq(0, 50, 0.025)

phos.jag.plot <- ggplot(data = df.jags1, aes(x = phos)) +
  
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  
  geom_line(data = df.jags2, aes(x = phos, y = mean), color = "darkorange", size = 1) + 
  
  geom_line(data = df.jags3, aes(x = phos, y = mean), color = "forestgreen", size = 1) + 
  
  geom_line(data = df.jags4, aes(x = phos, y = mean), color = "goldenrod1", size = 1) + 
  
  geom_line(data = df.jags5, aes(x = phos, y = mean), color = "dodgerblue3", size = 1) + 
  
  scale_x_continuous(limits = c(0, 10)) + 
  scale_y_continuous(limits = c(-0.25, 2.25)) + # Customize the axes and labels +
  labs(
    x = expression(paste("Phosphorous (concentration)")),
    y = "Growth rate",
    title = "Monod function, phosphorous limitation"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

phos.jag.plot
