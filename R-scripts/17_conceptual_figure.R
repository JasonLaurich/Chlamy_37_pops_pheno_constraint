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
  
  geom_segment(
    aes(x = 13.05, y = 8.65, xend = 13.30, yend = 8.4),
    inherit.aes = FALSE,
    colour = "blue3", linewidth = 1.5,
    arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  geom_segment(
           x = 12.95,  y = 8.75,
           xend = 12.7,   yend = 9.0,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  geom_segment(
           x = 13.55,  y = 9.15,
           xend = 13.9,   yend = 8.8,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  geom_segment(
           x = 13.45,  y = 9.25,
           xend = 13.1,   yend = 9.6,
           colour = "blue3", linewidth = 1.5,
           arrow = arrow(length = unit(3, "mm"), type = "closed")
  ) +
  
  coord_cartesian(xlim = c(12.5,14), ylim = c(8.2,9.7), expand = FALSE)

p.A

# Create panel B — variation across individuals/ genotypes ----------------

df.B <- data.frame(
  x = c(5.7, 6.2, 7.3, 7.4, 8.5, 7.8, 8.8, 9.1, 7.9, 9.5, 10.8, 9.9, 10.2, 11.6, 11.4, 11.3, 13.9, 9.1, 12.5, 13.0, 10.8),
  y = c(5.8, 8.1, 6.5, 9.3, 6.3, 8.0, 12.0, 7.5, 10.6, 9.4, 8.2, 10.8, 12.9, 12.1, 7.4, 9.7, 10.1, 14.1, 10.55, 8.7, 13.1)
)

poly.band <- data.frame(
  x = c(5, 7, 15, 9, 5),
  y = c(5, 5, 9, 15, 7)
)

p.B <- ggplot(df.B, aes(x = x, y = y)) +
  
  geom_polygon(data = poly.band, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  scale_y_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  scale_x_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  
  theme_classic() +
  
  geom_segment(
    aes(x = 5.1, y = 5.1, xend = 11.9, yend = 11.9),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 1.5,
    arrow = arrow(length = unit(2.1, "mm"), ends = "both", type = "closed")) +
  
  geom_segment(
    x = 15,  y = 9,
    xend = 9,   yend = 15,
    colour = "goldenrod2", linewidth = 1.5) +
  
  # annotate(
  #  "rect",
  #  xmin = 12.5, xmax = 14,
  #  ymin = 8.2, ymax = 9.7,
  #  fill = NA, colour = "black", linewidth = 1.5
  #) +
  
  geom_point(size= 3) +
  
  labs(x = "Trait 1 (e.g. maximum growth rate)",    
       y = "Trait 2 (e.g. stress tolerance)", 
       title = "A — genotypic variation") +  # labels
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
#  annotate("text", x = 12.5, y = 12.5, label = "Pareto front", size = 5.1, fontface = "bold", colour = "goldenrod2", angle = -45) +
  
  coord_cartesian(xlim = c(5,18), ylim = c(5,18), expand = FALSE) 

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
       title = "B — genotypic variation") +  # labels
  
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
  x = c(5.7, 6.2, 7.3, 7.4, 8.5, 7.8, 8.8, 9.1, 7.9, 9.5, 10.8, 9.9, 10.2, 11.6, 11.4, 11.3, 13.9, 9.1, 12.5, 13.0, 10.8, 13.55, 11.25, 9.9, 15.2, 12.85, 16.6, 15.8, 12.5, 10.5, 14.0),
  y = c(5.8, 8.1, 6.5, 9.3, 6.3, 8.0, 12.0, 7.5, 10.6, 9.4, 8.2, 10.8, 12.9, 12.1, 7.4, 9.7, 10.1, 14.1, 10.55, 8.7, 13.1, 13.3, 15.85, 16.7, 11.9, 14.3, 10.1, 10.7, 13.4, 15.1, 12.1)
)

p.C2 <- ggplot(df.C2, aes(x = x, y = y)) +
  
  geom_polygon(data = poly.band, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  scale_y_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  scale_x_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  
  theme_classic() +
  
  geom_segment(
    x = 15,  y = 9,
    xend = 9,   yend = 15,
    colour = "goldenrod2", linewidth = 1.5) +
  
  geom_segment(
    x = 17,  y = 10,
    xend = 10,   yend = 17,
    colour = "goldenrod2", linewidth = 1.5) +
  
  geom_segment(
    aes(x = 5.1, y = 5.1, xend = 13.4, yend = 13.4),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 1.5,
    arrow = arrow(length = unit(2.1, "mm"), ends = "both", type = "closed")) +
  
  geom_point(size= 3) +
  
  labs(x = "Trait 1 (e.g. maximum growth rate)",    
       y = "Trait 2 (e.g. stress tolerance)",  
       title = "B — genotypic variation") +  # labels
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  # annotate("text", x = 14, y = 14, label = "Pareto fronts", size = 5.1, fontface = "bold", colour = "goldenrod2", angle = -45) +
  
  coord_cartesian(xlim = c(5,18), ylim = c(5,18), expand = FALSE) 

p.C2

# Create panel E - interspecific variation (broad, positive)---------------------------------------

poly.band.E <- data.frame(
  x = c(2, 20, 9.5),
  y = c(2, 9.5, 20)
)

df.E <- data.frame(
  x = c(3.30, 4.7, 7.5, 4.8, 12.2, 11.2, 6.6, 7.3, 10.2, 13.1, 15.1, 16.7, 7.7, 12.8, 10.3, 9.7, 14.2, 19.5, 10.1),
  y = c(2.75, 8.1, 4.9, 5.1, 12.4, 6.5, 12.0, 9.5, 8.80, 6.70, 11.9, 9.3, 14.3, 9.70, 12.8, 16.6, 14.6, 9.6, 19.3)
)

df.E <- df.E %>% 
  mutate(x.ci = runif(n(),0.2,2.5), y.ci = runif(n(),0.2,2)) # Randomly assign spread values representing variation in ellipse shape and orientation


p.E <- ggplot(df.E, aes(x = x, y = y)) +
  
  geom_polygon(data = poly.band.E, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  geom_segment(
    x = 20,  y = 9.5,
    xend = 9.5,   yend = 20,
    colour = "goldenrod2", linewidth = 1.5) +
  
  geom_point(size= 3) +
  
  scale_y_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  
  geom_spoke(aes(angle = pi/4,        radius = y.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = pi/4 + pi,   radius = y.ci), linewidth = 0.6) +
  
  geom_spoke(aes(angle = -pi/4,       radius = x.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = -pi/4 + pi,  radius = x.ci), linewidth = 0.6) +
  
  labs(x = "Trait 1 (e.g. maximum growth rate)",    
       y = "Trait 2 (e.g. stress tolerance)",  
       title = "D — variation among species") +  # labels
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  # annotate("text", x = 15.25, y = 15.25, label = "Pareto front", size = 5.1, fontface = "bold", colour = "goldenrod2", angle = -45) +
  
  coord_cartesian(xlim = c(0,21), ylim = c(0,21), expand = FALSE)

p.E

# Create panel D - interspecific variation (negative, narrower) ---------------------

poly.band.D <- data.frame(             # Create a polygon! We'll set the lower limits at y and x 7.5. Slope is 7.5/18.
  x = c(17, 20, 9.5, 7.833333),        # solve for x = 16? ▲y = 4*(7.5/18) = 1.666667. y = 7.833333. Reverse for top left
  y = c(7.833333, 9.5, 20, 17)
)

df.D <- df.E %>% 
  filter(x >= 14.5 | y >= 14.5)

df.D2 <- data.frame(
  x = c(10.5, 12.5, 17.1, 12.5, 15.9, 13.6),
  y = c(15.1, 16.3, 10.6, 14.9, 13.1, 12.4)
)

df.D2 <- df.D2 %>% 
  mutate(x.ci = runif(n(),0.2,2.5), y.ci = runif(n(),0.2,2)) # Randomly assign spread values representing variation in ellipse shape and orientation

df.D <- rbind(df.D, df.D2)
  
p.D <- ggplot(df.D, aes(x = x, y = y)) +
  geom_point(size= 3) +
  
  geom_polygon(data = poly.band.D, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  geom_segment(
    x = 20,  y = 9.5,
    xend = 9.5,   yend = 20,
    colour = "goldenrod2", linewidth = 1.5) +
  
  scale_y_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  
  geom_spoke(aes(angle = pi/4,        radius = y.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = pi/4 + pi,   radius = y.ci), linewidth = 0.6) +
  
  geom_spoke(aes(angle = -pi/4,       radius = x.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = -pi/4 + pi,  radius = x.ci), linewidth = 0.6) +
  
  labs(x = "Trait 1 (e.g. maximum growth rate)",    
       y = "Trait 2 (e.g. stress tolerance)", 
       title = "C — variation among species") +  # labels
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  coord_cartesian(xlim = c(0,21), ylim = c(0,21), expand = FALSE) 
  
  # annotate("text", x = 15.25, y = 15.25, label = "Pareto front", size = 5.1, fontface = "bold", colour = "goldenrod2", angle = -45)

p.D

top   <- plot_spacer() + p.A + plot_spacer() + plot_layout(widths = c(0.25, 0.5, 0.25))
p.final <- top / (p.B | p.C2) / (p.D | p.E)
p.final

p.final2 <- (p.B | p.C2) / (p.D | p.E)
p.final2

ggsave("figures/15_fig_conceptual.jpeg", p.final, width = 15, height = 21)
ggsave("figures/15_fig_conceptual.pdf", p.final, width = 15, height = 21)

ggsave("figures/15_fig_conceptual_v2.jpeg", p.final2, width = 15, height = 15)
