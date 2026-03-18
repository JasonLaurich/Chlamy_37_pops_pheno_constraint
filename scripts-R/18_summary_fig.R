# Jason R Laurich

# March 17, 2026

# This script will create a final summary figure that communicates the paper's key results in one single figure.

# We'll import a table with P values (PF and Evol) then assess 50th quantile QRs on z-standardized data, and plot relationships between the two. 

# Inputs: 27_summary_table.csv, 73_summary_fig_data.csv
# Outputs: in figures-main : 05_summary_fig

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(cowplot)
library(quantreg)
library(ggforce) # for half-circle plotting

# Load & examine the data -------------------------------------------------

df.raw <- read.csv("processed-data/27_summary_table.csv") # Summary data - variation in niche-determining traits
head(df.raw)

df.sum <- read.csv("processed-data/73_summary_fig_data.csv") # P-value summaries for comparison (Pareto fronts and evolutionary optimization)
head(df.sum)

# Standardize and fit 50th qr ---------------------------------------------

qrs <- numeric(nrow(df.sum)) # initialize a quantile regression vector

for (i in 1:nrow(df.sum)){
  
  x <- df.sum$trait.x[i]
  y <- df.sum$trait.y[i]
  
  x.var <- df.raw %>% 
    select(x)
  
  z.x <- scale(x.var)
  
  y.var <- df.raw %>% 
    select(y)
  
  z.y <- scale(y.var)
  
  df <- cbind(z.x, z.y)
  
  df <- as.data.frame(df)
  
  q50  <- rq(df[,2] ~ df[,1], tau = 0.50, data = df) 
  
  qrs[i] <- q50$coefficients[2]

}

qrs <- as.data.frame(qrs)
df.sum <- cbind(df.sum, qrs)

# Plotting ----------------------------------------------------------------

df.sum$colour.x <- c("gold", "magenta3", "firebrick", "blue", "black", "magenta3", "firebrick", "blue", "black", "firebrick", "blue", "black", "blue", "black", "black")
df.sum$colour.y <- c("gold", "magenta3", "firebrick", "blue", "black", "gold", "gold", "gold", "gold", "magenta3", "magenta3", "magenta3", "firebrick", "firebrick", "blue")

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

p.A <- ggplot(df.arc2) +
  
  geom_arc_bar(
    aes(
      x0 = qrs,
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
    x = "Standardized 50th quantile regression slope",
    y = expression("Pareto front significance: "-log[10](italic("P")))
  )

p.A

# We need to stretch out the x axis without disrupting the circles

stretch <- 2.4

df.arc2 <- df.arc2 %>%
  mutate(qrs.plot = qrs * stretch)

r = 0.06

p.A <- ggplot(df.arc2) +
  
  geom_arc_bar(
    aes(
      x0 = qrs.plot,
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
    xmin = -Inf, xmax = Inf,
    ymin = 0, ymax = 1.2,
    fill = "white",
    alpha = 0.6
  ) +
  
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
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)
  ) +
  
  scale_fill_identity() +
  
  labs(
    x = "Standardized 50th quantile regression slope",
    y = expression("Pareto front significance: " -log[10](italic("P"))),
    title = "A"
  )

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
  mutate(qrs.plot = qrs * stretch)

p.B <- ggplot(df.arc3) +
  
  geom_arc_bar(
    aes(
      x0 = qrs.plot,
      y0 = p.evol.jit,
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
    xmin = -Inf, xmax = Inf,
    ymin = 0, ymax = 1.2,
    fill = "white",
    alpha = 0.6
  ) +
  
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
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)
  ) +
  scale_fill_identity() +
  
  labs(
    x = "Standardized 50th quantile regression slope",
    y = expression("Evolutionary optimization significance: " -log[10](italic("P"))),
    title = "B"
  )

p.B

# Panel C

# create micro and macro columns (significance of these Pfs across scales

df.sum$micro <- c(1,1,1,1,1, 0,1,0,1, 1,1,0, 1,0, 1)

df.sum2 <- df.sum %>%  # Filter out salt data (not present in macro)
  filter(trait.x != "S.c", trait.y != "S.c")

df.sum2$macro <- c(1,0,1,1, 1,1,1, 0,1, 0)

df.C <- df.sum2 %>%
  mutate(
    grad.x = sub("\\..*$", "", trait.x),
    grad.y = sub("\\..*$", "", trait.y)
  ) %>%
  select(grad.x, grad.y, micro, macro)

grad.order <- c("I", "N", "P", "T")

df.C <- df.C %>%
  mutate(
    grad.x = factor(grad.x, levels = grad.order),
    grad.y = factor(grad.y, levels = rev(grad.order))
  )

df.C.long <- df.C %>%
  pivot_longer(
    cols = c(micro, macro),
    names_to = "scale",
    values_to = "sig"
  ) %>%
  mutate(
    x.num = as.numeric(grad.x),
    y.num = as.numeric(grad.y),
    y.plot = ifelse(scale == "micro", y.num + 0.25, y.num - 0.25)
  )

p.C <- ggplot(df.C.long) +
  
  geom_tile(
    aes(
      x = x.num,
      y = y.plot,
      fill = factor(sig)
    ),
    width = 0.9,
    height = 0.45,
    color = "black",
    linewidth = 0.4
  ) +
  
  scale_fill_manual(
    values = c("0" = "white", "1" = "black")
  ) +
  
  scale_x_continuous(
    breaks = seq_along(grad.order),
    labels = grad.order,
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = seq_along(rev(grad.order)),
    labels = rev(grad.order),
    expand = c(0, 0)
  ) +
  
  coord_fixed() +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.ticks = element_blank()
  )

df.C <- df.sum2 %>%
  mutate(
    grad.x = sub("\\..*$", "", trait.x),
    grad.y = sub("\\..*$", "", trait.y)
  ) %>%
  select(grad.x, grad.y, micro, macro) %>%
  mutate(
    grad.x = factor(grad.x, levels = grad.order),
    grad.y = factor(grad.y, levels = rev(grad.order)),
    x = as.numeric(grad.x),
    y = as.numeric(grad.y)
  )

tri.micro <- df.C %>%
  rowwise() %>%
  do({
    tibble(
      grad.x = .$grad.x,
      grad.y = .$grad.y,
      x = .$x,
      y = .$y,
      sig = .$micro,
      part = "micro",
      px = c(.$x - 0.45, .$x - 0.45, .$x + 0.45),
      py = c(.$y + 0.45, .$y - 0.45, .$y + 0.45)
    )
  }) %>%
  ungroup()

tri.macro <- df.C %>%
  rowwise() %>%
  do({
    tibble(
      grad.x = .$grad.x,
      grad.y = .$grad.y,
      x = .$x,
      y = .$y,
      sig = .$macro,
      part = "macro",
      px = c(.$x + 0.45, .$x - 0.45, .$x + 0.45),
      py = c(.$y - 0.45, .$y - 0.45, .$y + 0.45)
    )
  }) %>%
  ungroup()

tri.C <- bind_rows(tri.micro, tri.macro) %>%
  mutate(
    fill = case_when(
      part == "micro" & sig == 1 ~ "firebrick",
      part == "macro" & sig == 1 ~ "dodgerblue3",
      TRUE ~ "white"
    ),
    id = paste(grad.x, grad.y, part, row_number(), sep = "_")
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
    aes(px, py, group = id, fill = fill),
    color = "black",
    linewidth = 0.4
  ) +
  
  scale_fill_identity() +
  
  scale_x_continuous(
    breaks = seq_along(grad.order),
    labels = grad.order,
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = seq_along(rev(grad.order)),
    labels = rev(grad.order),
    expand = c(0, 0)
  ) +
  
  coord_fixed() +
  
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.ticks = element_blank(),
    legend.position = "none"
  )

p.C

df.C <- df.sum2 %>%
  mutate(
    grad.x = sub("\\..*$", "", trait.x),
    grad.y = sub("\\..*$", "", trait.y)
  ) %>%
  select(grad.x, grad.y, micro, macro) %>%
  mutate(
    grad.x = factor(grad.x, levels = grad.order),
    grad.y = factor(grad.y, levels = rev(grad.order)),
    x = as.numeric(grad.x),
    y = as.numeric(grad.y),
    cell.id = paste(grad.x, grad.y, sep = "_")
  )

tri.micro <- df.C %>%
  rowwise() %>%
  do({
    tibble(
      grad.x = .$grad.x,
      grad.y = .$grad.y,
      x = .$x,
      y = .$y,
      sig = .$micro,
      part = "micro",
      poly.id = paste(.$cell.id, "micro", sep = "_"),
      px = c(.$x - 0.45, .$x - 0.45, .$x + 0.45),
      py = c(.$y + 0.45, .$y - 0.45, .$y + 0.45)
    )
  }) %>%
  ungroup()

tri.macro <- df.C %>%
  rowwise() %>%
  do({
    tibble(
      grad.x = .$grad.x,
      grad.y = .$grad.y,
      x = .$x,
      y = .$y,
      sig = .$macro,
      part = "macro",
      poly.id = paste(.$cell.id, "macro", sep = "_"),
      px = c(.$x + 0.45, .$x - 0.45, .$x + 0.45),
      py = c(.$y - 0.45, .$y - 0.45, .$y + 0.45)
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
    labels = grad.order,
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    breaks = seq_along(rev(grad.order)),
    labels = rev(grad.order),
    expand = c(0, 0)
  ) +
  
  coord_fixed() +
  
  theme_classic() +
  
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  ) +
  
  labs(
    x = "Niche axis 1",
    y = "Niche axis 2",
    title = "C"
  )

p.C

summ.fig <- plot_grid(
  p.A, p.B, p.C,
  ncol = 3,
  align = "hv",
  axis = "tblr"
)
