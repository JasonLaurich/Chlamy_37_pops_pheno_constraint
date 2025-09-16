# Jason R Laurich

# September 9th, 2025

# Creating an introductory figure to lay out hypotheses and predictions/ serve as a roadmap for the paper

# 3 panels?

# Idea - working off of the Metcalfe 2025 JEB paper: one panel could be at the level of organism? The three ways in which trade-offs
# can emerge as a result of energy ceilings? (3) plastic ceiling - location of PF could change within a single org.
# How does this relate to our experiments? Physiological changes (epigenetic) or adaptive evolution of this frontier?

# Idea - working of the Y-model resource-dependent trade off idea (e.g. Johnson and Nassrullah, 2025) - trade-offs could be detected 
# with quantile regressions? What about at the bleeding edge (PFs)?!?
# Also more variation in resource acquisition (e.g. 1/R*s) leads to weaker/ positive correlations? Could test this with my data
# Both across the whole dataset and among Pareto-optimal points?

# Question - what is the formal definition of Pareto fronts in biology? Would they apply at the level of individuals, or just within sp. and across them?

# Garland et al 2022: (1) allocation constraints, (2) functional conflicts, (3) shared biochemical pathways, (4) antagonistic pleiotropy,
# (5) ecological circumstances/ selection regime, (6) sexual selection. 
# 1-3 all for sure involved in our study. 

# Idea - do I go back and fit models to each replicate within each population? Would allow me to assess variation in resource acquisition within experimental
# populations and its effects on weaker negative correlations? OR positive? (e.g. Mark Johnston's question). Also the Haave-Audet meta-analysis'.... point 
# about looking at trade-offs within individuals (here really within genetically similar populations) and among them. 

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(cowplot)

# Create panel A (inset) — variation w/in individuals  --------------------

df.A2 <- data.frame(
  x = c(3, 7),
  y = c(23, 27)
)

poly.band <- data.frame(
  x = c(2, 10, 5, 0),
  y = c(20, 25, 30, 22)
)

p.A2 <- ggplot(df.A2, aes(x = x, y = y)) +
  geom_point(size= 3) +
  
  scale_y_continuous(limits = c(20, 30), breaks = c(20, 25, 30)) +
  scale_x_continuous(limits = c(0, 10), breaks = c(0, 5, 10)) +
  
  geom_segment(aes(x = 2, xend = 10, y = 20, yend = 25), colour = "black", size = 1.1) +
  geom_segment(aes(x = 0, xend = 5, y = 22, yend = 30), colour = "black", size = 1.1) +
  
  geom_segment(
    aes(x = 3.2, y = 23.2, xend = 6.8, yend = 26.8),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 1.5,
    arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 3.2,  y = 22.8,
           xend = 3.8,   yend = 21.4,
           curvature = -0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 2.8,  y = 23.2,
           xend = 1.4,   yend = 23.8,
           curvature = 0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 7.2,  y = 26.8,
           xend = 8.3,   yend = 24.2,
           curvature = -0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 6.8,  y = 27.2,
           xend = 4.4,   yend = 28.3,
           curvature = 0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  geom_polygon(data = poly.band, aes(x, y),
               fill = "grey80", alpha = 0.35, colour = NA,
               inherit.aes = FALSE) +
  
  theme_classic() +
  
  labs(x = "Salt tolerance (c)",    
       y = "Thermal breadth (°C)", 
       title = "Inset — individual plasiticity") +  # labels
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

p.A2
