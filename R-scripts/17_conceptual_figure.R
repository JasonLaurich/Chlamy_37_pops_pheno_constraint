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
library(patchwork)

# Create panel A — variation w/in individuals  --------------------

df.A <- data.frame(
  y = c(8.7, 9.2),
  x = c(13, 13.5)
)

p.A <- ggplot(df.A, aes(x = x, y = y)) +
  geom_point(size= 3) +
  
  geom_polygon(data = poly.crv, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  geom_segment(aes(x = 12.3, xend = 14, y = 7.7, yend = 8.7), colour = "black", size = 1.1) +
  geom_segment(aes(x = 12, xend = 13, y = 8.0, yend = 9.7), colour = "black", size = 1.1) +
  
  theme_classic() +
  
  
  labs(x = "",    
       y = "", 
       title = "A — plastic variation") +  # labels
  
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

p.A

# Create panel B — variation across individuals/ genotypes ----------------

df.B <- data.frame(
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

p.B <- ggplot(df.B, aes(x = x, y = y)) +
  
  scale_y_continuous(limits = c(5, 16), breaks = c(5, 10, 15)) +
  scale_x_continuous(limits = c(5, 16), breaks = c(5, 10, 15)) +
  
  theme_classic() +
  
  geom_segment(
    aes(x = 5.1, y = 5.1, xend = 13.4, yend = 13.4),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 1.5,
    arrow = arrow(length = unit(2.1, "mm"), ends = "both", type = "closed")) +
  
  geom_line(data = curve_pts, aes(x = x, y = y), color = "orchid2", size = 1.5, inherit.aes = FALSE) +  # Pareto front
  
  annotate(
    "rect",
    xmin = 12, xmax = 14,
    ymin = 7.7, ymax = 9.7,
    fill = NA, colour = "black", linewidth = 1.5
  ) +
  
  geom_polygon(data = poly.crv, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  geom_point(size= 3) +
  
  labs(x = "Trait 1 (e.g. nutrient acquisition)",    
       y = "Trait 2 (e.g. stress tolerance)", 
       title = "B — genotypic variation") +  # labels
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  annotate("text", x = 13.8, y = 13.8, label = "Pareto front", size = 5.1, fontface = "bold", colour = "orchid2", angle = -45) +
  
  coord_cartesian(xlim = c(5,16), ylim = c(5,16), expand = FALSE) 

p.B

# Create panel C — variation across individuals/ genotypes with a gap ----------------

df.C <- data.frame(
  x = c(5.7, 6.2, 7.3, 7.4, 8.5, 7.8, 9.1, 7.90, 9.5, 10.8, 12.0, 14.1, 13.55, 13.1, 10.70, 14.55, 9.10, 11.90, 14.45, 10.1, 14.7, 13.1),
  y = c(5.8, 8.1, 6.5, 9.3, 6.3, 8.0, 7.5, 10.6, 9.4, 8.20, 13.6, 10.2, 12.30, 13.8, 14.55, 11.10, 14.1, 14.45, 11.90, 14.7, 9.10, 13.0)
)

p.C <- ggplot(df.C, aes(x = x, y = y)) +
  
  scale_y_continuous(limits = c(5, 16), breaks = c(5, 10, 15)) +
  scale_x_continuous(limits = c(5, 16), breaks = c(5, 10, 15)) +
  
  theme_classic() +
  
  geom_segment(
    aes(x = 5.1, y = 5.1, xend = 13.4, yend = 13.4),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 1.5,
    arrow = arrow(length = unit(2.1, "mm"), ends = "both", type = "closed")) +
  
  geom_line(data = curve_pts, aes(x = x, y = y), color = "orchid2", size = 1.5, inherit.aes = FALSE) +  # Pareto front
  
  geom_polygon(data = poly.crv, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  geom_point(size= 3) +
  
  labs(x = "Trait 1 (e.g. nutrient acquisition)",    
       y = "Trait 2 (e.g. stress tolerance)", 
       title = "C — genotypic variation") +  # labels
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  annotate("text", x = 13.8, y = 13.8, label = "Pareto front", size = 5.1, fontface = "bold", colour = "orchid2", angle = -45) +
  
  coord_cartesian(xlim = c(5,16), ylim = c(5,16), expand = FALSE) 

p.C

# Create panel C.2 — variation across individuals/ genotypes with a gap, PF extended into new space ----------------

df.C2 <- data.frame(
  x = c(5.7, 6.2, 7.3, 7.4, 8.5, 7.8, 8.8, 9.1, 7.9, 9.5, 10.8, 10.1, 10.2, 11.2, 11.4, 11.3, 12.0, 13.9, 13.55, 13.1, 10.7, 14.55, 9.1, 12.5, 11.9, 14.45, 13.0, 15.25, 10.3, 12.4, 16.1, 14.1, 16.8, 16.5),
  y = c(5.8, 8.1, 6.5, 9.3, 6.3, 8.0, 12.0, 7.5, 10.6, 9.4, 8.2, 10.1, 12.9, 11.7, 7.4, 9.7, 13.6, 10.1, 12.3, 13.8, 14.55, 11.1, 14.1, 10.55, 14.45, 11.9, 8.70, 15.10, 17.0, 16.7, 14.2, 15.9, 11.2, 12.3)
)

# Create a second larger quadratic Bezier function to capture the curved pareto frontier

M.2  <- c(15.3, 15.3)      # point the curve should pass through
P0.2 <- c(17.1, 9.6)       # right endpoint
P2.2 <- c(9.6, 17.1)       # top endpoint

P1.2 <- (M.2 - (1 - t0)^2 * P0.2 - t0^2 * P2.2) / (2 * (1 - t0) * t0)  # (10,10)

t <- seq(0, 1, length.out = 80)          # curve from P0 -> P2 (right -> top)
curve_pts.2 <- data.frame(
  x = (1 - t)^2 * P0.2[1] + 2 * (1 - t) * t * P1.2[1] + t^2 * P2.2[1],
  y = (1 - t)^2 * P0.2[2] + 2 * (1 - t) * t * P1.2[2] + t^2 * P2.2[2]
)

p.C2 <- ggplot(df.C2, aes(x = x, y = y)) +
  
  scale_y_continuous(limits = c(5, 19), breaks = c(5, 10, 15)) +
  scale_x_continuous(limits = c(5, 19), breaks = c(5, 10, 15)) +
  
  theme_classic() +
  
  geom_segment(
    aes(x = 5.1, y = 5.1, xend = 15.2, yend = 15.2),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 1.5,
    arrow = arrow(length = unit(2.1, "mm"), ends = "both", type = "closed")) +
  
  geom_line(data = curve_pts, aes(x = x, y = y), color = "orchid2", size = 1.5, inherit.aes = FALSE) +  # Pareto front
  
  geom_line(data = curve_pts.2, aes(x = x, y = y), color = "orchid2", size = 1.5, inherit.aes = FALSE) +  # 2nd Pareto front
  
  geom_polygon(data = poly.crv, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  geom_point(size= 3) +
  
  labs(x = "Trait 1 (e.g. nutrient acquisition)",    
       y = "Trait 2 (e.g. stress tolerance)", 
       title = "C — genotypic variation") +  # labels
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  annotate("text", x = 15.8, y = 15.8, label = "Pareto fronts", size = 5.1, fontface = "bold", colour = "orchid2", angle = -45) +
  
  coord_cartesian(xlim = c(5,19), ylim = c(5,19), expand = FALSE) 

p.C2

# Create panel D - interspecific variation (broad, positive)---------------------------------------

poly.band.D <- data.frame(
  x = c(2, 20, 9.5),
  y = c(2, 9.5, 20)
)

df.D <- data.frame(
  x = c(3.30, 4.7, 7.5, 4.8, 12.2, 11.2, 6.6, 7.3, 10.2, 13.1, 15.1, 16.7, 7.7, 12.8, 10.3, 9.7, 14.2, 19.5, 10.1),
  y = c(2.75, 8.1, 4.9, 5.1, 12.4, 6.5, 12.0, 9.5, 8.80, 6.70, 11.9, 9.3, 14.3, 9.70, 12.8, 16.6, 14.6, 9.6, 19.3)
)

df.D <- df.D %>% 
  mutate(x.ci = runif(n(),0.2,2.5), y.ci = runif(n(),0.2,2)) # Randomly assign spread values representing variation in ellipse shape and orientation


p.D <- ggplot(df.D, aes(x = x, y = y)) +
  geom_point(size= 3) +
  
  scale_y_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  
  geom_spoke(aes(angle = pi/4,        radius = y.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = pi/4 + pi,   radius = y.ci), linewidth = 0.6) +
  
  geom_spoke(aes(angle = -pi/4,       radius = x.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = -pi/4 + pi,  radius = x.ci), linewidth = 0.6) +
  
  geom_polygon(data = poly.band.D, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  labs(x = "Trait 1 (e.g. nutrient acquisition)",    
       y = "Trait 2 (e.g. stress tolerance)", 
       title = "D — variation among species") +  # labels
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  coord_cartesian(xlim = c(0,21), ylim = c(0,21), expand = FALSE)

p.D

# Create panel E - interspecific variation (negative, narrower) ---------------------

poly.band.E <- data.frame(             # Create a polygon! We'll set the lower limits at y and x 7.5. Slope is 7.5/18.
  x = c(17, 20, 9.5, 7.833333),        # solve for x = 16? ▲y = 4*(7.5/18) = 1.666667. y = 7.833333. Reverse for top left
  y = c(7.833333, 9.5, 20, 17)
)

df.E <- df.D %>% 
  filter(x >= 14.5 | y >= 14.5)

df.E2 <- data.frame(
  x = c(10.5, 12.5, 17.1, 12.5, 15.9, 13.6),
  y = c(15.1, 16.3, 10.6, 14.9, 13.1, 12.4)
)

df.E2 <- df.E2 %>% 
  mutate(x.ci = runif(n(),0.2,2.5), y.ci = runif(n(),0.2,2)) # Randomly assign spread values representing variation in ellipse shape and orientation

df.E <- rbind(df.E, df.E2)
  
p.E <- ggplot(df.E, aes(x = x, y = y)) +
  geom_point(size= 3) +
  
  scale_y_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  
  geom_spoke(aes(angle = pi/4,        radius = y.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = pi/4 + pi,   radius = y.ci), linewidth = 0.6) +
  
  geom_spoke(aes(angle = -pi/4,       radius = x.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = -pi/4 + pi,  radius = x.ci), linewidth = 0.6) +
  
  geom_polygon(data = poly.band.E, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  labs(x = "Trait 1 (e.g. nutrient acquisition)",    
       y = "Trait 2 (e.g. stress tolerance)", 
       title = "E — variation among species") +  # labels
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  coord_cartesian(xlim = c(0,21), ylim = c(0,21), expand = FALSE)

p.E

top   <- plot_spacer() + p.A + plot_spacer() + plot_layout(widths = c(0.25, 0.5, 0.25))
p.final <- top / (p.B | p.C2) / (p.D | p.E)
p.final

ggsave("figures/15_fig_conceptual.jpeg", p.final, width = 15, height = 21)
ggsave("figures/15_fig_conceptual.pdf", p.final, width = 15, height = 21)
