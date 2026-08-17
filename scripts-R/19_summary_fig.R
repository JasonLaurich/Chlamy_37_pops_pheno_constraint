# Jason R Laurich

# March 17, 2026

# This script will create a final summary figure that communicates the paper's key results in one single figure.

# We'll import a table with P values (PF and Evol) then assess 50th quantile corrs on z-standardized data, and plot relationships between the two. 

#Using Pearson's r now

# Updated with new results from the overhauled evolutionary optimization approach.

# Inputs: 27_summary_table.csv, 79_summary_fig_data.csv
# Outputs: in figures-main : 05_fig_5_summary.pdf

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(cowplot)
library(quantreg)
library(ggforce) # for half-circle plotting

# Load & examine the data -------------------------------------------------

df.raw <- read.csv("processed-data/27_summary_table.csv") # Summary data - variation in niche-determining traits
head(df.raw)

df.sum <- read.csv("processed-data/79_summary_fig_data.csv") # P-value summaries for comparison (Pareto fronts and evolutionary optimization)
head(df.sum)

# Estimate Pearson correlations -------------------------------------------

corrs <- numeric(nrow(df.sum))

for (i in 1:nrow(df.sum)){
  
  x <- df.sum$trait.x[i]
  y <- df.sum$trait.y[i]
  
  corrs[i] <- cor(df.raw[[x]], df.raw[[y]],
                  method = "pearson")
}

df.sum$corrs <- corrs

# Plotting ----------------------------------------------------------------

df.sum$colour.x <- c("goldenrod2", "plum3", "brown4", "navyblue", "black", "plum3", "brown4", "navyblue", "black", "brown4", "navyblue", "black", "navyblue", "black", "black")
df.sum$colour.y <- c("goldenrod2", "plum3", "brown4", "navyblue", "black", "goldenrod2", "goldenrod2", "goldenrod2", "goldenrod2", "plum3", "plum3", "plum3", "brown4", "brown4", "navyblue")

# scale the p values

df.sum <- df.sum %>% 
  mutate(p.par = -log10(P.pareto),
         p.evol = -log10(P.evol))

df.plot <- df.sum %>%
  mutate(id = row_number()) %>%
  select(id, everything())

df.arc <- bind_rows(
  df.plot %>%
    mutate(start = 0, end = pi, fill = colour.x),
  
  df.plot %>%
    mutate(start = pi, end = 2*pi, fill = colour.y)
)

# Panel A: Pareto fronts and quantile regressions. 

r <- 0.06

#Jitter the p-values, where equal to 3 (too much overlap)

df.sum <- df.sum %>%
  group_by(p.par) %>%
  mutate(
    offset = (row_number() - mean(row_number())) * 0.03,
    p.par.jit = p.par + ifelse(p.par > 2.8, offset, 0)
  ) %>%
  ungroup()

df.plot2 <- df.sum %>%
  mutate(id = row_number()) %>%
  select(id, everything())

df.arc2 <- bind_rows(
  df.plot2 %>%
    mutate(start = 0, end = pi, fill = colour.x),
  
  df.plot2 %>%
    mutate(start = pi, end = 2*pi, fill = colour.y)
)

x.mid <- mean(range(df.arc$corrs, na.rm = TRUE))
y.mid <- mean(range(df.arc$p.par, na.rm = TRUE))

x.span <- diff(range(df.arc$corrs, na.rm = TRUE))
y.span <- diff(range(df.arc$p.par, na.rm = TRUE))

span <- max(x.span, y.span)
pad  <- 0.05 * span
span <- span + 2 * pad

p.A <- ggplot(df.arc2) +
  
  geom_arc_bar(
    aes(
      x0 = corrs,
      y0 = p.par.jit,
      r0 = 0,
      r = r,
      start = start,
      end = end,
      fill = fill,
      group = interaction(id, start)
    ),
    
    colour = "black",
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  
  annotate(
    "rect",
    xmin = -1, xmax = 1,
    ymin = 0, ymax = 1.2,
    fill = "white",
    alpha = 0.6
  ) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  coord_fixed(
    ratio = 1,
    xlim = c(x.mid - span/2, x.mid + span/2),
    ylim = c(y.mid - span/2, y.mid + span/2)
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  scale_fill_identity() +
  
  labs(
    x = "Pearson's correlation (r)",
    y = expression("Pareto front significance: "-log[10](italic("P")))
  )

p.A

# We need to stretch out the x axis without disrupting the circles

stretch <- 2.7

df.arc2 <- df.arc2 %>%
  mutate(corrs.plot = corrs * stretch)

df.arc2 <- df.arc2 %>%
  mutate(alpha.par = ifelse(P.pareto < 0.05, 0.6, 0.3))

r = 0.075

x.leg <- -1.65
y.leg <- 0.6
dy <- 0.18

# Some specific jittering to make the figure more legible

df.arc2 <- df.arc2 %>%
  mutate(
    p.par.jit = ifelse(
      trait.x == "I.µ.max" & trait.y == "I.comp",
      p.par.jit + 0.15,
      p.par.jit
    )
  )

df.arc2 <- df.arc2 %>%
  mutate(
    p.par.jit = ifelse(
      trait.x == "N.µ.max" & trait.y == "N.comp",
      p.par.jit - 0.15,
      p.par.jit
    )
  )

df.arc2 <- df.arc2 %>%
  mutate(
    p.par.jit = ifelse(
      trait.x == "S.c" & trait.y == "N.comp",
      p.par.jit + 0.05,
      p.par.jit
    )
  )

df.arc2 <- df.arc2 %>%
  mutate(
    p.par.jit = ifelse(
      trait.x == "T.br" & trait.y == "I.comp",
      p.par.jit - 0.10,
      p.par.jit
    )
  )

df.arc2 <- df.arc2 %>%
  mutate(
    p.par.jit = ifelse(
      trait.x == "T.br" & trait.y == "S.c",
      p.par.jit - 0.20,
      p.par.jit
    )
  )

p.A <- ggplot(df.arc2) +
  
  geom_arc_bar(
    aes(
      x0 = corrs.plot,
      y0 = p.par.jit,
      r0 = 0,
      r = r*2.1,
      start = start,
      end = end,
      fill = fill,
      alpha = alpha.par,                 # <- added
      group = interaction(id, start)
    ),
    colour = "black",
    linewidth = 0.3,
    show.legend = FALSE                  # <- dropped `alpha = 0.6,`
  ) +
  
  scale_fill_identity() +
  scale_alpha_identity() +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  coord_fixed(ratio = 1) +
  
  scale_x_continuous(
    limits = c(-0.75 * stretch, 0.5 * stretch),
    breaks = stretch * seq(-1, 1, by = 0.25),
    labels = seq(-1, 1, by = 0.25)
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) +
  
  scale_fill_identity() +
  
  labs(
    x = "Trait correlation (Pearson's r)",
    y = expression("Strength of Pareto front (" -log[10](italic("P)"))),
    title = "A"
  ) +
  
  annotate("text", x = -1.4, y = 3.3, label = "Trade-off", size = 3, fontface = "bold") +
  
  annotate("text", x = 0.75, y = 3.2, label = "Pareto \nconstraint", size = 3, fontface = "bold") +
  
  annotate("point", x = x.leg, y = y.leg,     colour = "goldenrod2", size = 3, alpha = 0.6) +
  annotate("point", x = x.leg, y = y.leg-dy,  colour = "plum3", size = 3, alpha = 0.6) +
  annotate("point", x = x.leg, y = y.leg-2*dy,colour = "brown4", size = 3, alpha = 0.6) +
  annotate("point", x = x.leg, y = y.leg-3*dy,colour = "navyblue", size = 3, alpha = 0.6) +
  annotate("point", x = x.leg, y = y.leg-4*dy,colour = "black", size = 3, alpha = 0.6) +
  
  annotate("text", x = x.leg + 0.1, y = y.leg,     
           label = "Light", hjust = 0, size = 3) +
  
  annotate("text", x = x.leg + 0.1, y = y.leg-dy,  
           label = "Nitrogen", hjust = 0, size = 3) +
  
  annotate("text", x = x.leg + 0.1, y = y.leg-2*dy,
           label = "Phosphorus", hjust = 0, size = 3) +
  
  annotate("text", x = x.leg + 0.1, y = y.leg-3*dy,
           label = "Salt", hjust = 0, size = 3) +
  
  annotate("text", x = x.leg + 0.1, y = y.leg-4*dy,
           label = "Temperature", hjust = 0, size = 3) +
  
  annotate("text", x = x.leg, y = y.leg + 0.18,
           label = "Niche axis",
           hjust = 0, fontface = "bold", size = 3)

p.A

# Panel B: Evolutionary optimization and quantile regressions. 

#Jitter 

df.sum <- df.sum %>%
  group_by(p.evol) %>%
  mutate(
    offset = (row_number() - mean(row_number())) * 0.03,
    p.evol.jit = p.evol + ifelse(p.evol > 2.8, offset, 0)
  ) %>%
  ungroup()

df.plot3 <- df.sum %>%
  mutate(id = row_number()) %>%
  select(id, everything())

df.arc3 <- bind_rows(
  df.plot3 %>%
    mutate(start = 0, end = pi, fill = colour.x),
  
  df.plot3 %>%
    mutate(start = pi, end = 2*pi, fill = colour.y)
)

df.arc3 <- df.arc3 %>%
  mutate(
    corrs.plot = corrs * stretch,                     # was qrs.plot = qrs * stretch
    alpha.evol = ifelse(P.evol < 0.05, 0.6, 0.3)
  )

# Need to bump down nitrogen v phos a little, so the label fits in
# Also need to space out the all salt and phos a little - too much overlap. Let's do temp ~ phos and salt too.

df.arc3 <- df.arc3 %>%
  mutate(
    p.evol.jit = ifelse(
      trait.x == "S.µ.max" & trait.y == "S.c",
      p.evol.jit - 0.1,
      p.evol.jit
    )
  )

df.arc3 <- df.arc3 %>%
  mutate(
    p.evol.jit = ifelse(
      trait.x == "S.c" & trait.y == "I.comp",
      p.evol.jit - 0.15,
      p.evol.jit
    )
  )

df.arc3 <- df.arc3 %>%
  mutate(
    p.evol.jit = ifelse(
      trait.x == "P.comp" & trait.y == "N.comp",
      p.evol.jit - 0.20,
      p.evol.jit
    )
  )

df.arc3 <- df.arc3 %>%
  mutate(
    p.evol.jit = ifelse(
      trait.x == "P.comp" & trait.y == "I.comp",
      p.evol.jit - 0.10,
      p.evol.jit
    )
  )

p.B <- ggplot(df.arc3) +
  
  geom_arc_bar(
    aes(
      x0 = corrs.plot,                    # was qrs.plot
      y0 = p.evol.jit,
      r0 = 0,
      r = r*2.1,
      start = start,
      end = end,
      fill = fill,
      alpha = alpha.evol,                 # <- added
      group = interaction(id, start)
    ),
    colour = "black",                     # <- was `ccolour` (typo); dropped `alpha = 0.6`
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  
  scale_fill_identity() +
  scale_alpha_identity() +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  coord_fixed(ratio = 1) +
  
  scale_x_continuous(
    limits = c(-0.75 * stretch, 0.5 * stretch),
    breaks = stretch * seq(-1, 1, by = 0.25),
    labels = seq(-1, 1, by = 0.25)
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) +
  scale_fill_identity() +
  
  labs(
    x = "Trait correlation (Pearson's r)",
    y = expression("Strength of evolutionary shift (" -log[10](italic("P)"))),
    title = "B"
  )+
  
  annotate("text", x = -1.4, y = 3.3, label = "Modularity", size = 3, fontface = "bold") +
  
  annotate("text", x = 0.75, y = 3.2, label = "Correlated \nevolution", size = 3, fontface = "bold") 

p.B

# Panel C

# create micro and macro columns (significance of these Pfs across scales)

df.sum$micro <- c(1,1,1,1,1, 0,1,0,1, 1,1,0, 1,0, 1)

df.sum2 <- df.sum %>%  # Filter out salt data (not present in macro)
  filter(trait.x != "S.c", trait.y != "S.c")

df.sum2$macro <- c(1,0,1,1, 1,1,1, 0,1, 0)

grad.order <- c("I", "N", "P", "T")

grad.labels <- c("Light", "Nit", "Phos", "Temp")

df.C <- df.sum2 %>%
  mutate(
    grad.x = sub("\\..*$", "", trait.x),
    grad.y.raw = sub("\\..*$", "", trait.y),
    x.index = match(grad.x, grad.order),
    y.index = match(grad.y.raw, grad.order)
  ) %>%
  mutate(
    grad.x = factor(grad.x, levels = grad.order),
    grad.y = factor(grad.y.raw, levels = rev(grad.order)),
    x = as.numeric(grad.x),
    y = as.numeric(grad.y),
    cell.id = paste(grad.x, grad.y.raw, sep = "_")
  ) %>%
  select(grad.x, grad.y, grad.y.raw, x, y, cell.id, micro, macro)

tri.micro <- df.C %>%
  rowwise() %>%
  do({
    tibble(
      grad.x  = .$grad.x,
      grad.y  = .$grad.y,
      x       = .$x,
      y       = .$y,
      sig     = .$micro,
      part    = "micro",
      poly.id = paste(.$cell.id, "micro", sep = "_"),
      px      = c(.$x - 0.45, .$x - 0.45, .$x + 0.45),
      py      = c(.$y + 0.45, .$y - 0.45, .$y + 0.45)
    )
  }) %>%
  ungroup()

tri.macro <- df.C %>%
  rowwise() %>%
  do({
    tibble(
      grad.x  = .$grad.x,
      grad.y  = .$grad.y,
      x       = .$x,
      y       = .$y,
      sig     = .$macro,
      part    = "macro",
      poly.id = paste(.$cell.id, "macro", sep = "_"),
      px      = c(.$x + 0.45, .$x - 0.45, .$x + 0.45),
      py      = c(.$y - 0.45, .$y - 0.45, .$y + 0.45)
    )
  }) %>%
  ungroup()

tri.C <- bind_rows(tri.micro, tri.macro) %>%
  mutate(
    fill = case_when(
      part == "micro" & sig == 1 ~ "firebrick",
      part == "macro" & sig == 1 ~ "dodgerblue3",
      TRUE ~ "white"
    )
  )

square.C <- df.C

p.C <- ggplot() +
  
  geom_tile(
    data = square.C,
    aes(x = x, y = y),
    width = 0.9,
    height = 0.9,
    fill = NA,
    color = "black",
    linewidth = 0.6
  ) +
  
  geom_polygon(
    data = tri.C,
    aes(x = px, y = py, group = poly.id, fill = fill),
    color = "black",
    linewidth = 0.4
  ) +
  
  scale_fill_identity() +
  
  scale_x_continuous(
    breaks = seq_along(grad.order),
    labels = grad.labels,
    expand = c(0, 0)
    ) +
  
  scale_y_continuous(
    breaks = seq_along(rev(grad.order)),
    labels = rev(grad.labels),
    expand = c(0, 0)
  ) +
  
  coord_fixed() +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  labs(
    x = "Niche axis 1",
    y = "Niche axis 2",
    title = "C"
  ) +
  
  annotate(
    "polygon",
    x = c(0.75, 0.75, 1.05),
    y = c(1.45, 1.15, 1.45),
    fill = "firebrick",
    color = "black",
    linewidth = 0.3
  ) +
  
  # blue legend triangle
  annotate(
    "polygon",
    x = c(1.05, 0.75, 1.05),
    y = c(0.95, 0.95, 1.25),
    fill = "dodgerblue3",
    color = "black",
    linewidth = 0.3
  ) +
  
  # red legend text
  annotate(
    "text",
    x = 1.15,
    y = 1.4,
    label = "Microevolutionary scale",
    hjust = 0,
    size = 3,
    fontface = "bold"
  ) +
  
  # blue legend text
  annotate(
    "text",
    x = 1.15,
    y = 1.05,
    label = "Macroevolutionary scale",
    hjust = 0,
    size = 3,
    fontface = "bold"
  )

p.C

# Compile the figure ------------------------------------------------------

summ.fig <- plot_grid(p.A, p.B, p.C, ncol = 1)

ggsave("figures-main/05_fig_5_summary.pdf", summ.fig, width = 4, height = 12)
