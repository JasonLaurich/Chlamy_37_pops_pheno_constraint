# Jason R Laurich

# June 7th, 2026

# lightly edited for re-submission on August 16th, 2026.

# This script will generate a new Figure 1, a conceptual figure mapping out Pareto fronts and the concepts of multivariate 
  # phenotypic constraint among genotypes and species.

# Inputs: none
# Outputs: in figures-main : 01_conceptual_figure.jpeg

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(cowplot)
library(patchwork)
library(scam)
library(smatr)

nearest_pf <- function(x0, y0){ # calculate closest point on the pf
  
  df.pf.dense %>%
    mutate(
      dist = sqrt((x - x0)^2 + (y - y0)^2)
    ) %>%
    slice_min(dist, n = 1, with_ties = FALSE) %>%
    select(x.pf = x,
           y.pf = y)
  
}

# Create panel A — negative correlation w/ Pareto front --------------------

df.pf <- data.frame(   # Create a fake Pareto front for visualization
  x = c(9.1, 9.75, 10.4, 10.5, 11.2, 11.9, 13.4, 14.4),
  y = c(15.1, 13.8, 13.1, 12.4, 11.2, 10.3, 10.1, 9.0)
)

fit <- scam(y ~ s(x, bs = "mpd", k = 6), data = df.pf)

x.vals <- seq(min(df.pf$x), max(df.pf$x), length.out = 100)

df.pred <- data.frame(x = x.vals) %>%
  mutate(y = predict(fit, newdata = .))

fit <- scam(y ~ s(x, bs = "mpd", k = 6), data = df.pf) # Fit a scam to the PF

x.vals <- seq(min(df.pf$x), max(df.pf$x), length.out = 100) # Generate an x sequence for plotting

df.pred <- data.frame(x = x.vals) %>%
  mutate(y = predict(fit, newdata = .))

p.pf <- ggplot() +
  
  geom_point(
    data = df.pf,
    aes(x = x, y = y),
    size = 2
  ) +
  
  geom_line(
    data = df.pred,
    aes(x = x, y = y),
    colour = "goldenrod2",
    linewidth = 1
  ) +
  
  scale_x_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  scale_y_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  
  labs(
    x = "Trait 1",
    y = "Trait 2"
  ) +
  
  theme_classic() +
  coord_cartesian(xlim = c(5, 18), ylim = c(5, 18), expand = FALSE)

p.pf

df.A <- data.frame(
  x = c(
    7.4, 7.7, 8.0, 8.3, 8.6,
    8.9, 9.2, 9.5, 9.8, 10.1,
    10.4, 10.7, 11.0, 11.3, 11.6,
    11.9, 12.2, 12.5, 12.8, 13.1,
    13.4, 13.7, 14.0, 14.3
  ),
  y = c(
    14, 13.2, 14.5, 15.0, 14.4,
    13.2, 13.9, 12.9, 11.5, 11.9,
    10.8, 8.9, 9.1, 10.3, 8.6,
    9.2, 9.9, 7.6, 8.9, 8.5,
    9.1, 9.5, 8.3, 8.6
  )
)

df.A.all <- bind_rows(df.pf, df.A)

sma.A <- sma(y ~ x, data = df.A.all, method = "SMA")

xr <- range(df.A.all$x, na.rm = TRUE)
seg.A <- data.frame(
  x    = xr[1],
  xend = xr[2],
  y    = coef(sma.A)[1] + coef(sma.A)[2] * xr[1],
  yend = coef(sma.A)[1] + coef(sma.A)[2] * xr[2]
)

p.A <- ggplot() +
  
  geom_point(
    data = df.pf,
    aes(x = x, y = y),
    size = 4
  ) +
  
  geom_point(
    data = df.A,
    aes(x = x, y = y),
    size = 2
  ) +
  
  geom_line(
    data = df.pred,
    aes(x = x, y = y),
    colour = "goldenrod2",
    linewidth = 1
  ) +
  
  geom_segment(
    data = seg.A,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "red3", linewidth = 0.75
  ) +
  
  scale_x_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  scale_y_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  
  labs(
    x = "Trait 1 (e.g. growth rate)",    
    y = "Trait 2 (e.g. stress tolerance)", 
    title = "A) Classical trade-off"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text = element_text(size = 10, face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) +
  
  coord_cartesian(xlim = c(5, 15.5), ylim = c(5, 15.5), expand = FALSE)

p.A

# Create panel B — positive correlation w/ Pareto front --------------------

df.B <- data.frame(
  x = c(
    6.2, 6.8, 7.4, 8.0, 8.5,
    8.9, 9.2, 9.5, 9.8, 10.1,
    10.4, 10.7, 11.0, 11.3, 11.6,
    11.9, 12.2, 12.8,
    13.4, 13.7, 14.3, 6.5, 7.1,
    7.3, 7.1, 9.1, 9, 8.1,
    7.5, 10.1, 9.5, 7.8, 8.0,
    7.4, 7.8, 8.3, 8.7, 9,
    7.5, 6.3, 7.1, 9.9, 8
  ),
  y = c(
    5.8, 6.3, 6.1, 7.0, 7.6,
    8.4, 9.1, 8.0, 10.5, 9.3,
    11.7, 8.6, 10.1, 9.0, 8.8, 
    10.1, 9.5, 9.2,
    9.8, 8.9, 8.7, 7.5, 9.1,
    6.9, 7.1, 12.5, 14, 11.1,
    8, 13.2, 12.4, 8, 10.5,
    9.1, 9.5, 8.1, 11.1, 10.1,
    13.1, 6.5, 7.5, 12.1, 9
  )
)

df.B.all <- bind_rows(df.pf, df.B)

sma.B <- sma(y ~ x, data = df.B.all, method = "SMA")

xr <- range(df.B.all$x, na.rm = TRUE)
seg.B <- data.frame(
  x    = xr[1],
  xend = xr[2],
  y    = coef(sma.B)[1] + coef(sma.B)[2] * xr[1],
  yend = coef(sma.B)[1] + coef(sma.B)[2] * xr[2]
)

p.B <- ggplot() +
  
  geom_point(
    data = df.pf,
    aes(x = x, y = y),
    size = 4
  ) +
  
  geom_point(
    data = df.B,
    aes(x = x, y = y),
    size = 2
  ) +
  
  geom_line(
    data = df.pred,
    aes(x = x, y = y),
    colour = "goldenrod2",
    linewidth = 1
  ) +
  
  geom_segment(
    data = seg.B,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "red3", linewidth = 0.75
  ) +
  
  scale_x_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  scale_y_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  
  labs(
    x = "Trait 1 (e.g. growth rate)",    
    y = "Trait 2 (e.g. stress tolerance)", 
    title = "B) Masked trade-off"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text = element_text(size = 10, face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) +
  
  coord_cartesian(xlim = c(5, 15.5), ylim = c(5, 15.5), expand = FALSE)

p.B

# Create panel C — statistical testing of Pareto fronts and evolutionary optimization ----------------

df.pf.match <- data.frame(   # Create a fake Pareto front for visualization
  x = c(9.75, 10.4, 13.4),
  y = c(13.8, 13.1, 10.1)
)

df.pf.other <- data.frame(   # Create a fake Pareto front for visualization
  x = c(9.1, 10.5, 11.2, 11.9, 14.4),
  y = c(15.1, 12.4, 11.2, 10.3, 9.0)
)

df.C.match <- data.frame(  # break df.B into matching, ancestral and not. Here matching
  x = c(
    10.1,
    11.0, 11.6,
    12.2,
    14.3,
    9.1, 7.5, 9
  ),
  y = c(
    9.3,
    10.1, 8.8, 
    9.5,
    8.7,
    12.5, 12.5, 10.1
  )
)

df.C.oth <- data.frame(  # break df.B into matching, ancestral and not. Here other
  x = c(
    8.0, 6.2, 7.4,
    9.2, 9.5, 9.8,
    10.4, 10.7, 11.3,
    11.9, 12.8,
    13.4, 13.7, 
    7.1, 9, 8.1,
    7.5, 10.1, 9.5, 7.8, 8.0,
    7.4, 7.8, 8.3, 8.7,
    6.3, 9.9
  ),
  y = c(
    7.0, 5.8, 6.1,
    9.1, 8.0, 10.5,
    11.7, 8.6, 9.0,
    10.1, 9.2,
    9.8, 8.9, 
    7.1, 14, 11.1,
    8, 13.2, 12.4, 8, 10.5,
    9.1, 9.5, 8.1, 11.1,
    6.5, 12.1
  )
)

df.C.anc <- data.frame(  # break df.B into matching, ancestral and not. Here ancestral
  x = c(
    6.8, 8.5,
    8.9,
    6.5, 7.1,
    7.3, 
    7.1, 8
  ),
  y = c(
    6.3, 7.6,
    8.4,
    7.5, 9.1,
    6.9, 
    7.5, 9
  )
)

# Calculate the closest point on the Pareto front for a couple of points

pf.approx <- with( # Dense version of actual segmented PF
  df.pf %>% arrange(x),
  approx(x = x, y = y, xout = seq(min(x), max(x), length.out = 1000))
)

df.pf.dense <- data.frame(
  x = pf.approx$x,
  y = pf.approx$y
)

df.dist <- data.frame(
  x = c(7.1, 8.5, 9.1, 11.6),
  y = c(7.5, 7.6, 12.5, 8.8),
  type = c("ancestral", "ancestral",
           "matching", "matching")
)

df.lines <- df.dist %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf(x, y))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

# Create a polygon for plotting the phenotpyic space rendered inaccesible by the Pareto front

x.max <- max(df.pf$x)
y.max <- max(df.pf$y)

df.inaccessible <- bind_rows(
  data.frame(x = min(df.pf.dense$x), y = y.max),
  data.frame(x = x.max, y = y.max),
  data.frame(x = x.max, y = df.pf.dense$y[df.pf.dense$x == max(df.pf.dense$x)]),
  df.pf.dense %>%
    arrange(desc(x))
)

p.C <- ggplot() +
  
  geom_polygon(
    data = df.inaccessible,
    aes(x = x, y = y),
    fill = "grey60",
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  geom_segment(
    data = subset(df.lines, type == "ancestral"),
    aes(x = x, y = y,
        xend = x.pf,
        yend = y.pf),
    colour = "black",
    alpha = 0.6,
    linewidth = 0.8,
    linetype = "31"
  ) +
  
  geom_segment(
    data = subset(df.lines, type == "matching"),
    aes(x = x, y = y,
        xend = x.pf,
        yend = y.pf),
    colour = "blue3",
    alpha = 0.6,
    linewidth = 0.8,
    linetype = "31"
  ) +
  
  geom_point(
    data = df.pf.other,
    aes(x = x, y = y),
    size = 4,
    shape = 16,
    alpha = 0.6
  ) +
  
  geom_point(
    data = df.pf.match,
    aes(x = x, y = y),
    shape = 21,
    size = 4,
    stroke = 1.5,
    fill = scales::alpha("blue3", 0.6)
  ) +
  
  geom_point(
    data = df.C.match,
    aes(x = x, y = y),
    shape = 21,
    size = 2,
    stroke = 1.5,
    fill = scales::alpha("blue3", 0.6)
  ) +
  
  geom_point(
    data = df.C.oth,
    aes(x = x, y = y),
    size = 2,
    shape = 16,
    alpha = 0.6
  ) +
  
  geom_point(
    data = df.C.anc,
    aes(x = x, y = y),
    shape = 5, 
    size = 2,
    colour = "black",
    stroke = 1,
    alpha = 0.6
  ) +
  
  geom_path(
    data = df.pf,
    aes(x = x, y = y),
    colour = "goldenrod2",
    linewidth = 0.8,
    linetype = "31"
  ) +
  
  scale_x_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  scale_y_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  
  labs(
    x = "Trait 1 (e.g. growth rate)",    
    y = "Trait 2 (e.g. stress tolerance)", 
    title = "C) Evolution towards Pareto fronts"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text = element_text(size = 10, face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) +
  
  coord_cartesian(xlim = c(5, 15.5), ylim = c(5, 15.5), expand = FALSE)

p.C

# Create panel D - interspecific variation (broad, positive)---------------------------------------

df.D <- data.frame(
  x = c(6.8, 7.2, 7.5, 7.8, 8.0, 8.4, 8.8, 9.1, 9.5,
        9.8, 10.2, 10.6, 11.0, 11.5, 12.0, 12.4, 12.8, 13.2,
        8.2, 9.1, 5.9, 6.6, 11),
  y = c(13.8, 11.2, 10.8, 13.0, 9.5, 11.5, 13.2, 10.2, 9.8,
        8.9, 11.2, 9.4, 10.8, 8.5, 9.1, 7.8, 8.6, 7.4,
        7.5, 6.6, 5.7, 8.6, 7.5)
)

df.pf.D <- data.frame(   # PF anchor points — should sit at/near the outer edge of df.D
  x = c(7.2, 8.8, 9.8, 11.0, 13.2, 14.5, 13.5),
  y = c(14.0, 13.5, 13.1, 11.5, 8.7, 6.1, 8.5)
)

ggplot() +
  geom_point(data = df.D, aes(x, y), size = 1.5, alpha = 0.6, colour = "black") +
  geom_point(data = df.pf.D, aes(x, y), colour = "goldenrod2", size = 3) +
  scale_x_continuous(limits = c(5, 15)) +
  scale_y_continuous(limits = c(5, 15)) +
  theme_classic()

x.max <- max(df.pf.D$x)
y.max <- max(df.pf.D$y)

pf.approx.D <- with( # Dense version of actual segmented PF
  df.pf.D %>% arrange(x),
  approx(x = x, y = y, xout = seq(min(x), max(x), length.out = 1000))
)

df.pf.dense <- data.frame(
  x = pf.approx.D$x,
  y = pf.approx.D$y
)

df.inaccessible <- bind_rows(
  data.frame(x = min(df.pf.dense$x), y = y.max),
  data.frame(x = x.max, y = y.max),
  data.frame(x = x.max, y = df.pf.dense$y[df.pf.dense$x == max(df.pf.dense$x)]),
  df.pf.dense %>%
    arrange(desc(x))
)

sma.D <- sma(y ~ x, data = df.D, method = "SMA")

xr <- range(df.D$x, na.rm = TRUE)
seg.D <- data.frame(
  x    = xr[1],
  xend = xr[2],
  y    = coef(sma.D)[1] + coef(sma.D)[2] * xr[1],
  yend = coef(sma.D)[1] + coef(sma.D)[2] * xr[2]
)

p.D <- ggplot(df.D, aes(x = x, y = y)) +
  
  geom_polygon(data = df.inaccessible, aes(x, y),
               fill = "grey60", colour = NA,
               inherit.aes = FALSE) +
  
  geom_path(data = df.pf.dense, aes(x = x, y = y),
            colour = "goldenrod2", linewidth = 0.75,
            inherit.aes = FALSE) +
  
  geom_point(size = 1.5, alpha = 0.6, colour = "black") +
  
  geom_point(data = df.pf.D, 
             size = 1.5, alpha =0.6, colour = "black") +
  
  scale_y_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  scale_x_continuous(limits = c(0, 21), breaks = c(0, 5, 10, 15, 20)) +
  
  labs(x = "Trait 1 (e.g. growth rate)",    
       y = "Trait 2 (e.g. stress tolerance)",  
       title = "D) Species-level Pareto fronts") +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text  = element_text(size = 10, face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) +
  
  coord_cartesian(xlim = c(3.5, 15), ylim = c(3.5, 15), expand = FALSE)

p.D

# Assemble the figure -----------------------------------------------------

p <- plot_grid(p.A, p.B, p.C, p.D, nrow = 2, align ='hv')
p

ggsave("figures-main/01_fig_1_conceptual_figure.jpeg", p, width = 7, height = 7) # aiming for ~ 2/3 of a page in width
