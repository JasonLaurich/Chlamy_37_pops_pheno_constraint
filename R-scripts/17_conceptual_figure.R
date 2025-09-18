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
  x = c(5.7, 6.2, 7.3, 7.4, 8.5, 7.8, 8.8, 9.1, 7.9, 9.5, 10.8, 10.1, 10.2, 11.2, 11.4, 11.3, 12.0, 13.9, 13.55, 13.1, 10.7, 14.55, 9.1, 12.5, 11.9, 14.45, 13.0),
  y = c(5.8, 8.1, 6.5, 9.3, 6.3, 8.0, 12.0, 7.5, 10.6, 9.4, 8.2, 10.1, 12.9, 11.7, 7.4, 9.7, 13.6, 10.1, 12.3, 13.8, 14.55, 11.1, 14.1, 10.55, 14.45, 11.9, 8.7)
)

poly.band <- data.frame(
  x = c(5, 7, 15, 9, 5),
  y = c(5, 5, 9, 15, 7)
)

# Create a quadratic Bezier function to capture the curved pareto frontier

M  <- c(13.5, 13.5)          # point the curve should pass through
P0 <- c(15, 9)               # right endpoint
P2 <- c(9, 15)               # top endpoint

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
  
  scale_y_continuous(limits = c(5, 15), breaks = c(5, 10, 15)) +
  scale_x_continuous(limits = c(5, 15), breaks = c(5, 10, 15)) +
  
  theme_classic() +
  
  geom_segment(
    aes(x = 5.1, y = 5.1, xend = 13.4, yend = 13.4),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 1.5) +
  
  annotate(
    "rect",
    xmin = 12, xmax = 14,
    ymin = 7.7, ymax = 9.7,
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
  
  coord_cartesian(xlim = c(5,15), ylim = c(5,15), expand = FALSE)

p.A


# Create panel A (inset) — variation w/in individuals  --------------------

df.A2 <- data.frame(
  y = c(8.7, 9.2),
  x = c(13, 13.5)
)

p.A2 <- ggplot(df.A2, aes(x = x, y = y)) +
  geom_point(size= 3) +
  
  geom_polygon(data = poly.crv, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  geom_segment(aes(x = 12.3, xend = 14, y = 7.7, yend = 8.7), colour = "black", size = 1.1) +
  geom_segment(aes(x = 12, xend = 13, y = 8.0, yend = 9.7), colour = "black", size = 1.1) +
  
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
    aes(x = 13.05, y = 8.75, xend = 13.45, yend = 9.15),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 1.5,
    arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 13.05,  y = 8.65,
           xend = 13.20,   yend = 8.3,
           curvature = -0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 12.95,  y = 8.75,
           xend = 12.6,   yend = 8.90,
           curvature = 0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 13.55,  y = 9.15,
           xend = 13.80,   yend = 8.67,
           curvature = -0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  annotate("curve",
           x = 13.45,  y = 9.25,
           xend = 12.97,   yend = 9.50,
           curvature = 0.15, ncp = 50,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  coord_cartesian(xlim = c(12,14), ylim = c(7.7,9.7), expand = FALSE)

p.A2

p.top <- plot_grid(
  p.A, p.A2,
  ncol = 2,
  rel_widths = c(1, 0.4),   # tweak this ratio to taste
  rel_heights = c(1, 0.4),
  align = ""              # align panels nicely
)

p.top

p.A2.zoom <- p.A2 +
  coord_equal(xlim = c(12, 14), ylim = c(7.7, 9.7), expand = FALSE) +
  theme(text = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.margin = margin(4, 4, 4, 4))


right_col <- plot_grid(
  NULL,
  p.A2.zoom,
  NULL,
  ncol = 1,
  rel_heights = c(0.3, 1, 0.3)   # <- controls inset height in the column
)

right_col

p.top <- plot_grid(
  p.A, right_col,
  ncol = 2,
  rel_widths = c(1, 1),            # <- inset column ~40% of the width; tweak
  align = "h"
)

p.top

# Create panel B - interspecific variation ---------------------------------------

df.B <- data.frame(
  x = c(1.9, 1.3, 2.8, 4.2, 2.2, 3.1, 7.5, 5.0, 6.0, 7.0, 7.0, 8.0, 9.0, 9.0, 10.0, 11.0, 11.0, 12.0, 13.0, 13.0, 14.0, 15.0, 15.0, 16.0, 17.0, 17.0, 18.0, 19.0, 19.0),
  y = c(1.2, 2.1, 2.1, 1.3, 3.8, 8.1, 4.2, 5.0, 6.0, 7.0, 7.0, 8.0, 9.0, 9.0, 10.0, 11.0, 11.0, 12.0, 13.0, 13.0, 14.0, 15.0, 15.0, 16.0, 17.0, 17.0, 18.0, 19.0, 19.0)
)

df.B <- df.B %>% 
  mutate(x.ci = runif(n(),0.2,2), y.ci = runif(n(),0.2,2)) # Randomly assign spread values representing variation in ellipse shape and orientation

poly.band.B <- data.frame(
  x = c(0, 5, 20, 8, 0),
  y = c(0, 0, 8, 20, 5)
)

# Create a quadratic Bezier function to capture the curved pareto frontier

M.B  <- c(17, 17)          # point the curve should pass through
P0.B <- c(20, 8)              # right endpoint
P2.B <- c(8, 20)              # top endpoint

t0 <- 0.5
P1.B <- (M.B - (1 - t0)^2 * P0.B - t0^2 * P2.B) / (2 * (1 - t0) * t0)  # (10,10)

t <- seq(0, 1, length.out = 80)          # curve from P0 -> P2 (right -> top)
curve_pts.B <- data.frame(
  x = (1 - t)^2 * P0.B[1] + 2 * (1 - t) * t * P1.B[1] + t^2 * P2.B[1],
  y = (1 - t)^2 * P0.B[2] + 2 * (1 - t) * t * P1.B[2] + t^2 * P2.B[2]
)

poly.crv.B <- rbind(
  poly.band.B[1:3,],                         # (0,0)->(2,0)->(10,4)
  curve_pts.B[-c(1, length(t)), ],         # curved edge (avoid duplicated ends)
  poly.band.B[4:5,]                          # (4,10)->(0,2)
)

df.B <- data.frame(
  x = c(1.9, 1.3, 3.3, 4.2, 2.1, 3.1, 7.5, 4.8, 5.7, 8.2, 11.2, 5.6, 8.1, 9.2, 13.1, 15.1, 16.7, 7.7, 12.1, 10.3, 9.7, 18.8, 14.5, 19.5, 11.1, 14.1, 17.6, 16.9, 8.3, 17.5),
  y = c(1.2, 2.1, 2.6, 1.3, 4.3, 8.1, 4.2, 5.1, 9.1, 7.1, 5.5, 12.0, 10.3, 8.8, 6.7, 11.9, 9.3, 14.3, 10.7, 12.8, 16.6, 12.3, 14.9, 9.3, 19.3, 17.6, 15.4, 16.8, 19.7, 13.4)
)

df.B <- df.B %>% 
  mutate(x.ci = runif(n(),0.2,2), y.ci = runif(n(),0.2,2)) # Randomly assign spread values representing variation in ellipse shape and orientation


p.B <- ggplot(df.B, aes(x = x, y = y)) +
  geom_point(size= 3) +
  
  scale_y_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  
  geom_spoke(aes(angle = pi/4,        radius = y.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = pi/4 + pi,   radius = y.ci), linewidth = 0.6) +
  
  geom_spoke(aes(angle = -pi/4,       radius = x.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = -pi/4 + pi,  radius = x.ci), linewidth = 0.6) +
  
  geom_polygon(data = poly.crv.B, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  labs(x = "Trait 1 (e.g. maximum growth rate)",    
       y = "Trait 2 (e.g. salt tolerance)", 
       title = "B — variation among species") +  # labels
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  theme_classic() +
  
  coord_cartesian(xlim = c(0,21), ylim = c(0,21), expand = FALSE)

p.B

p.B.null <- plot_grid(
  p.B, NULL,
  ncol = 2,
  rel_heights = c(1,1)   # <- controls inset height in the column
)

p.final <- plot_grid(
  p.top, p.B.null,
  nrow = 2,
  rel_widths = c(1, 1),            # <- inset column ~40% of the width; tweak
  align = "h"
)

p.final

ggsave("figures/15_fig_conceptual.jpeg", p.final, width = 15, height = 15)
