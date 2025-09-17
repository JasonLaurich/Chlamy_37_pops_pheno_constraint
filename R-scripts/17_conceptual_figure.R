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

# Create panel A — variation across individuals/ genotypes ----------------

df.A <- data.frame(
  x = c(0.7, 1.2, 2.3, 2.4, 3.5, 2.8, 3.8, 4.1, 2.9, 4.5, 5.8, 5.1, 6.2, 6.2, 6.4, 6.3, 7.0, 8.9, 8.55, 8.1, 5.7, 9.55, 4.1, 7.5, 6.9, 9.45, 8.0),
  y = c(0.8, 3.1, 1.5, 4.3, 1.3, 3.0, 7.0, 2.5, 5.6, 4.4, 8.2, 5.1, 3.5, 6.7, 2.4, 4.7, 8.6, 5.1, 7.3, 8.8, 9.55, 6.1, 9.1, 5.55, 9.45, 6.9, 3.7)
)

poly.band <- data.frame(
  x = c(0,2, 10, 4, 0),
  y = c(0, 0, 4, 10, 2)
)

# Create a quadratic Bezier function to capture the curved pareto frontier

M  <- c(8.5, 8.5)          # point the curve should pass through
P0 <- c(10, 4)             # right endpoint
P2 <- c(4, 10)             # top endpoint

t0 <- 0.5
P1 <- (M - (1 - t0)^2 * P0 - t0^2 * P2) / (2 * (1 - t0) * t0)  # (10,10)

t <- seq(0, 1, length.out = 80)          # curve from P0 -> P2 (right -> top)
curve_pts <- data.frame(
  x = (1 - t)^2 * P0[1] + 2 * (1 - t) * t * P1[1] + t^2 * P2[1],
  y = (1 - t)^2 * P0[2] + 2 * (1 - t) * t * P1[2] + t^2 * P2[2]
)

poly.crv <- rbind(
  poly.band[1:3,],                         # (0,0)->(2,0)->(10,4)
  curve_pts[-c(1, length(t)), ],           # curved edge (avoid duplicated ends)
  poly.band[4:5,]                          # (4,10)->(0,2)
)

p.A <- ggplot(df.A, aes(x = x, y = y)) +
  geom_point(size= 3) +
  
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 5, 10)) +
  scale_x_continuous(limits = c(0, 10), breaks = c(0, 5, 10)) +
  
  theme_classic() +
  
  geom_segment(
    aes(x = 0.1, y = 0.1, xend = 8.4, yend = 8.4),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 1.5) +
  
  annotate(
    "rect",
    xmin = 7, xmax = 9,
    ymin = 2.7, ymax = 4.7,
    fill = NA, colour = "black", linewidth = 1.5
  ) +
  
  geom_polygon(data = poly.crv, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  labs(x = "Trait 1 (e.g. maximum growth rate)",    
       y = "Trait 2 (e.g. salt tolerance)", 
       title = "A — variation among genotypes") +  # labels
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  coord_cartesian(xlim = c(0,10), ylim = c(0,10), expand = FALSE)

p.A


# Create panel A (inset) — variation w/in individuals  --------------------

df.A2 <- data.frame(
  y = c(3.7, 4.2),
  x = c(8, 8.5)
)

p.A2 <- ggplot(df.A2, aes(x = x, y = y)) +
  geom_point(size= 3) +
  
  geom_polygon(data = poly.crv, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  geom_segment(aes(x = 7.3, xend = 9, y = 2.7, yend = 3.7), colour = "black", size = 1.1) +
  geom_segment(aes(x = 7, xend = 8, y = 3.0, yend = 4.7), colour = "black", size = 1.1) +
  
  theme_classic() +

  
  labs(x = "",    
       y = "", 
       title = "Inset — individual plasiticity") +  # labels
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  geom_segment(
    aes(x = 8.05, y = 3.75, xend = 8.45, yend = 4.15),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 1.5,
    arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 8.05,  y = 3.65,
           xend = 8.20,   yend = 3.3,
           curvature = -0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 7.95,  y = 3.75,
           xend = 7.6,   yend = 3.90,
           curvature = 0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 8.55,  y = 4.15,
           xend = 8.80,   yend = 3.67,
           curvature = -0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 8.45,  y = 4.25,
           xend = 7.97,   yend = 4.50,
           curvature = 0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  coord_cartesian(xlim = c(7,9), ylim = c(2.7,4.7), expand = FALSE)

p.A2

