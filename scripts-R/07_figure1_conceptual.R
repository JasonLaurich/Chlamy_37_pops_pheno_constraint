# Jason R Laurich

# February 8th, 2026

# This script will generate Figure 1, a conceptual figure mapping out Pareto fronts and the concepts of multivariate 
  # phenotypic constraint among genotypes and species.

# Inputs: none
# Outputs: in figures-main : 01_conceptual_figure.jpeg

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(cowplot)
library(patchwork)

# Create panel A — genetic variation --------------------

df.A <- data.frame(
  x = c(5.7, 6.2, 7.3, 7.4, 8.5, 7.8, 8.8, 9.1, 7.9, 9.5, 10.8, 9.9, 10.2, 11.6, 11.4, 11.3, 13.9, 9.1, 12.5, 13.0, 10.8),
  y = c(5.8, 8.1, 6.5, 9.3, 6.3, 8.0, 12.0, 7.5, 10.6, 9.4, 8.2, 10.8, 12.9, 12.1, 7.4, 9.7, 10.1, 14.1, 10.55, 8.7, 13.1)
)

poly.band <- data.frame(
  x = c(5, 7, 15, 9, 5),
  y = c(5, 5, 9, 15, 7)
)

p.A <- ggplot(df.A, aes(x = x, y = y)) +
  
  geom_polygon(data = poly.band, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  scale_y_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  scale_x_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  
  theme_classic() +
  
  geom_segment(
    x = 15,  y = 9,
    xend = 9,   yend = 15,
    colour = "goldenrod2", linewidth = 0.75) +
  
  geom_segment(
    aes(x = 5.1, y = 5.1, xend = 11.9, yend = 11.9),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 0.75,
    arrow = arrow(length = unit(0.9, "mm"), ends = "both", type = "closed")) +
  
  geom_point(size= 1) +
  
  labs(x = "Trait 1 (e.g. growth rate)",    
       y = "Trait 2 (e.g. stress tolerance)", 
       title = "A — genotypes") +  # labels
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  coord_cartesian(xlim = c(5,18), ylim = c(5,18), expand = FALSE) 

p.A

# Create panel B — variation across individuals/ genotypes with a gap ----------------

df.B <- data.frame(
  x = c(5.7, 6.2, 7.3, 7.4, 8.5, 7.8, 8.8, 9.1, 7.9, 9.5, 10.8, 9.9, 10.2, 11.6, 11.4, 11.3, 13.9, 9.1, 12.5, 13.0, 10.8, 13.55, 11.25, 9.9, 15.2, 12.85, 16.6, 15.8, 12.5, 10.5, 14.0),
  y = c(5.8, 8.1, 6.5, 9.3, 6.3, 8.0, 12.0, 7.5, 10.6, 9.4, 8.2, 10.8, 12.9, 12.1, 7.4, 9.7, 10.1, 14.1, 10.55, 8.7, 13.1, 13.3, 15.85, 16.7, 11.9, 14.3, 10.1, 10.7, 13.4, 15.1, 12.1)
)

p.B <- ggplot(df.B, aes(x = x, y = y)) +
  
  geom_polygon(data = poly.band, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  scale_y_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  scale_x_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  
  theme_classic() +
  
  geom_segment(
    x = 15,  y = 9,
    xend = 9,   yend = 15,
    colour = "goldenrod2", linewidth = 0.75) +
  
  geom_segment(
    x = 17,  y = 10,
    xend = 10,   yend = 17,
    colour = "goldenrod2", linewidth = 0.75) +
  
  geom_segment(
    aes(x = 5.1, y = 5.1, xend = 13.4, yend = 13.4),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 0.75,
    arrow = arrow(length = unit(0.9, "mm"), ends = "both", type = "closed")) +
  
  geom_point(size= 1) +
  
  labs(x = "Trait 1 (e.g. growth rate)",    
       y = "Trait 2 (e.g. stress tolerance)",  
       title = "B — genotypes") +  # labels
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  # annotate("text", x = 14, y = 14, label = "Pareto fronts", size = 5.1, fontface = "bold", colour = "goldenrod2", angle = -45) +
  
  coord_cartesian(xlim = c(5,18), ylim = c(5,18), expand = FALSE) 

p.B

# Create panel C - interspecific variation (negative, narrower) ---------------------

poly.band.C <- data.frame(             # Create a polygon! We'll set the lower limits at y and x 7.5. Slope is 7.5/18.
  x = c(17, 20, 9.5, 7.833333),        # solve for x = 16? ▲y = 4*(7.5/18) = 1.666667. y = 7.833333. Reverse for top left
  y = c(7.833333, 9.5, 20, 17)
)

df.D <- data.frame(
  x = c(3.30, 4.7, 7.5, 4.8, 12.2, 11.2, 6.6, 7.3, 10.2, 13.1, 15.1, 16.7, 7.7, 12.8, 10.3, 9.7, 14.2, 19.5, 10.1),
  y = c(2.75, 8.1, 4.9, 5.1, 12.4, 6.5, 12.0, 9.5, 8.80, 6.70, 11.9, 9.3, 14.3, 9.70, 12.8, 16.6, 14.6, 9.6, 19.3)
)

df.D <- df.D %>% 
  mutate(x.ci = runif(n(),0.2,2.5), y.ci = runif(n(),0.2,2)) # Randomly assign spread values representing variation in ellipse shape and orientation


df.C <- df.D %>% 
  filter(x >= 14.5 | y >= 14.5)

df.C2 <- data.frame(
  x = c(10.5, 12.5, 17.1, 12.5, 15.9, 13.6),
  y = c(15.1, 16.3, 10.6, 14.9, 13.1, 12.4)
)

df.C2 <- df.C2 %>% 
  mutate(x.ci = runif(n(),0.2,2.5), y.ci = runif(n(),0.2,2)) # Randomly assign spread values representing variation in ellipse shape and orientation

df.C <- rbind(df.C, df.C2)

p.C <- ggplot(df.C, aes(x = x, y = y)) +
  
  geom_polygon(data = poly.band.C, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  geom_point(size= 1) +
  
  geom_segment(
    x = 20,  y = 9.5,
    xend = 9.5,   yend = 20,
    colour = "goldenrod2", linewidth = 0.75) +
  
  scale_y_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  
  geom_spoke(aes(angle = pi/4,        radius = y.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = pi/4 + pi,   radius = y.ci), linewidth = 0.6) +
  
  geom_spoke(aes(angle = -pi/4,       radius = x.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = -pi/4 + pi,  radius = x.ci), linewidth = 0.6) +
  
  labs(x = "Trait 1 (e.g. growth rate)",    
       y = "Trait 2 (e.g. stress tolerance)", 
       title = "C — species") +  # labels
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  coord_cartesian(xlim = c(0,21), ylim = c(0,21), expand = FALSE) 

p.C

# Create panel E - interspecific variation (broad, positive)---------------------------------------

poly.band.D <- data.frame(
  x = c(2, 20, 9.5),
  y = c(2, 9.5, 20)
)

p.D <- ggplot(df.D, aes(x = x, y = y)) +
  
  geom_polygon(data = poly.band.D, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  geom_segment(
    x = 20,  y = 9.5,
    xend = 9.5,   yend = 20,
    colour = "goldenrod2", linewidth = 0.75) +
  
  geom_point(size= 1) +
  
  scale_y_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  
  geom_spoke(aes(angle = pi/4,        radius = y.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = pi/4 + pi,   radius = y.ci), linewidth = 0.6) +
  
  geom_spoke(aes(angle = -pi/4,       radius = x.ci), linewidth = 0.6) +
  geom_spoke(aes(angle = -pi/4 + pi,  radius = x.ci), linewidth = 0.6) +
  
  labs(x = "Trait 1 (e.g. growth rate)",    
       y = "Trait 2 (e.g. stress tolerance)",  
       title = "D — species") +  # labels
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  # annotate("text", x = 15.25, y = 15.25, label = "Pareto front", size = 5.1, fontface = "bold", colour = "goldenrod2", angle = -45) +
  
  coord_cartesian(xlim = c(0,21), ylim = c(0,21), expand = FALSE)

p.D

# Assemble the figure -----------------------------------------------------

p <- (p.A | p.B) / (p.C | p.D)
p

ggsave("figures-main/01_conceptual_figure.jpeg", p, width = 5.5, height = 5.5) # aiming for ~ 2/3 of a page in width
