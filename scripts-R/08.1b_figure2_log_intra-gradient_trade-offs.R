# Jason R Laurich

# June 2026

# This script replicates 08.1_figure2_intra-gradient_trade-offs.R on LOG-TRANSFORMED data.
# Log transformation is applied selectively to niche/competitive traits that show pronounced
# right-skewed distributions (1/I*, 1/N*, 1/P*, salt tolerance c) — log is justified by
# distributional skew for these ratio-scale measurements. Thermal breadth (T.br) and all
# µ.max variables are left on the raw scale as they are approximately symmetric.
# This addresses reviewer concerns about unit choice while remaining data-driven.

# Inputs: 27_summary_table.csv
# Outputs: figures-misc/02b_fig_2_log_intra-gradient_tradeoffs.jpeg
# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(quantreg)
library(sp) # For the point.in.polygon function
library(scam)
library(pracma) # For calculating polygon area
library(cowplot)
library(lme4)
library(emmeans)
library(gridExtra)
library(grid)
library(smatr)

par_frt <- function(df, xvar, yvar) { # Simple Pareto front function / convex hull algorithm (one sided)
  
  df <- df[order(-df[[xvar]], df[[yvar]]), ]  
  pareto_points <- df[1, ]  # Start with the first point
  
  for (i in 2:nrow(df)) {
    if (df[i, yvar] >= tail(pareto_points[[yvar]], 1)) {  # Ensure increasing y values
      pareto_points <- rbind(pareto_points, df[i,])
    }
  }
  
  return(pareto_points)
}

closest_on_segment <- function(px, py, ax, ay, bx, by) { # find the pf front segment closest to each point
  
  apx <- px - ax
  apy <- py - ay
  abx <- bx - ax
  aby <- by - ay
  
  ab2 <- abx^2 + aby^2
  
  t <- (apx * abx + apy * aby) / ab2
  t <- pmax(0, pmin(1, t))  # constrain to segment
  
  x.closest <- ax + t * abx
  y.closest <- ay + t * aby
  
  dist <- sqrt((px - x.closest)^2 + (py - y.closest)^2)
  
  data.frame(
    x.pf = x.closest,
    y.pf = y.closest,
    dist.to.pf = dist
  )
}

nearest_pf_segment <- function(px, py) { # calculate distance to closest point on nearest pf segment. 
  
  df.pf.segments %>%
    rowwise() %>%
    mutate(
      closest = list(
        closest_on_segment(px, py, ax, ay, bx, by)
      )
    ) %>%
    tidyr::unnest(closest) %>%
    ungroup() %>%
    slice_min(dist.to.pf, n = 1, with_ties = FALSE) %>%
    select(x.pf, y.pf, dist.to.pf)
}

# Load & examine the data -------------------------------------------------

df <- read.csv("processed-data/27_summary_table.csv") # Summary file
head(df)

df$Evol.plt <- factor(df$Evol,   # Create nicer labels for plotting
                      levels = c("none", "L", "N", "P", "S", "B", "BS", "C"),
                      labels = c("Ancestral",
                                 "Light limitation",
                                 "Nitrogen limitation", 
                                 "Phosphorus limitation", 
                                 "Salt stress", 
                                 "Biotic depletion", 
                                 "Biotic depletion x Salt", 
                                 "Control"))

df$Anc.plt <- factor(df$Anc,     # Create nicer labels for plotting
                     levels = c("anc2", "anc3", "anc4", "anc5", "cc1690"),
                     labels = c("Ancestor 2", 
                                "Ancestor 3", 
                                "Ancestor 4", 
                                "Ancestor 5", 
                                "Mixed ancestry"))

# Distribution check: µ.max and niche traits (raw vs log) --------------------

vars_gr <- c("T.µ.max", "I.µ.max", "N.µ.max", "P.µ.max", "S.µ.max") # Reshape to long for faceting
vars_niche <- c("T.br", "I.comp", "N.comp", "P.comp", "S.c")

labels_gr <- c(
  "T.µ.max" = "Temperature µmax (day⁻¹)",
  "I.µ.max" = "Light µmax (day⁻¹)",
  "N.µ.max" = "Nitrogen µmax (day⁻¹)",
  "P.µ.max" = "Phosphorus µmax (day⁻¹)",
  "S.µ.max" = "Salt µmax (day⁻¹)"
)
labels_niche <- c(
  "T.br"   = "Thermal breadth T_br (°C)",
  "I.comp" = "Light comp. 1/I*",
  "N.comp" = "Nitrogen comp. 1/N*",
  "P.comp" = "Phosphorus comp. 1/P*",
  "S.c"    = "Salt tolerance c (g/L)"
)

make_dist_grid <- function(vars, labels, df, title) {
  
  long <- df %>%
    select(all_of(vars)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    filter(!is.na(value), value > 0) %>%
    mutate(
      log_value = log(value),
      label     = labels[variable],
      label     = factor(label, levels = labels)
    )
  
  p_raw <- ggplot(long, aes(x = value)) +
    geom_histogram(bins = 15, fill = "steelblue", colour = "white", alpha = 0.85) +
    facet_wrap(~ label, scales = "free", nrow = 1) +
    labs(title = paste(title, "— raw"), x = NULL, y = "Count") +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 8),
          plot.title = element_text(size = 10, face = "bold"))
  
  p_log <- ggplot(long, aes(x = log_value)) +
    geom_histogram(bins = 15, fill = "darkorange", colour = "white", alpha = 0.85) +
    facet_wrap(~ label, scales = "free", nrow = 1) +
    labs(title = paste(title, "— log-transformed"), x = NULL, y = "Count") +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 8),
          plot.title = element_text(size = 10, face = "bold"))
  
  grid.arrange(p_raw, p_log, nrow = 2)
}

make_dist_grid(vars_gr,    labels_gr,    df, "Growth rates (µmax)")
make_dist_grid(vars_niche, labels_niche, df, "Niche/competitive traits")

# Temperature: all logged -------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral", "evolved") # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = log(T.br),
    z.x = log(T.µ.max)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # Looks fine!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # Temp: r -0.6677927, p < 2.2e-16

sma.T <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.T)      # slope -05880105
coef(sma.T)   

T.qp2 <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  
  geom_point(
    data = subset(df.filt, evol.bin =="evolved" & pareto.opt == "N"),
    shape = 16, 
    size = 2,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="evolved" & pareto.opt == "Y"),
    shape = 16, 
    size = 4,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="ancestral"),
    shape = 5, 
    size = 2,
    colour = "black",
    stroke = 1,
    alpha = 0.6
  ) +
  
  geom_abline(intercept = coef(sma.T)[1], slope = coef(sma.T)[2], lwd = 0.6, linetype = "dashed") +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 0.6, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = expression(atop("log(Maximum growth rate)",
                           log(italic(μ)[italic(max)]))),
       y = expression(atop("log(Thermal breadth)",
                           log(italic("T")[italic(br)]))),
       color = "Evolutionary History",
       title = "E) Temperature") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "chocolate3",
               "Biotic depletion x Salt" = "skyblue",
               "Control" = "olivedrab4",
               "Light limitation" = "goldenrod2",
               "Nitrogen limitation" = "plum3",
               "Ancestral" = "black",
               "Phosphorus limitation" = "brown4",  
               "Salt stress" = "navyblue")
  ) +
  
  ylim(2.6, 3.1) +
  xlim(0.8, 1.7) + 
  
  annotate(
    "text",
    x = 0.8 + 0.6 * (1.7 - 0.8),
    y = 2.6 + 0.95 * (3.1 - 2.6),
    label = "PF\nr = -0.668*",                          # removing QR
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

T.qp2 # Display the plot

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(ax = z.x2, ay = z.y2, bx = lead(z.x2), by = lead(z.y2)) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(nearest = list(nearest_pf_segment(z.x2, z.y2))) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # none significantly closer

# Light: all logged -------------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "L", 'light', 'other')) # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = log(I.comp),
    z.x = log(I.µ.max)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # Light: r -0.6640078, p < 2.2e-16

sma.L <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.L)      # slope -4.655543
coef(sma.L)      

I.qp2 <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  
  geom_point(
    data = subset(df.filt, evol.bin =="other" & pareto.opt == "N"),
    shape = 16, 
    size = 2,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="other" & pareto.opt == "Y"),
    shape = 16, 
    size = 4,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="ancestral" & pareto.opt == "N"),
    shape = 5, 
    size = 2,
    colour = "black",
    stroke = 1,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="ancestral" & pareto.opt == "Y"),
    shape = 5, 
    size = 4,
    colour = "black",
    stroke = 1,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin == "light" & pareto.opt == "N"),
    shape = 21,
    size = 2,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("goldenrod2", 0.6)
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin == "light" & pareto.opt == "Y"),
    shape = 21,
    size = 4,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("goldenrod2", 0.6)
  ) +
  
  geom_abline(intercept = coef(sma.L)[1], slope = coef(sma.L)[2], lwd = 0.6, linetype = "dashed", inherit.aes = FALSE) +
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 0.6, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = expression(atop("log(Maximum growth rate)",
                           log(italic(μ)[italic(max)]))),
       y = expression(atop("log(Low-L tolerance)",
                           log(1/italic(L)^"*"))),
       color = "Evolutionary History",
       title = "A) Light") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "chocolate3",
               "Biotic depletion x Salt" = "skyblue",
               "Control" = "olivedrab4",
               "Light limitation" = "goldenrod2",
               "Nitrogen limitation" = "plum3",
               "Ancestral" = "black",
               "Phosphorus limitation" = "brown4",  
               "Salt stress" = "navyblue")
  ) +
  
  ylim(-3.2, -0.4) +
  xlim(0.3, 0.85) +
  
  annotate(
    "text",
    x = 0.3 + 0.6 * (0.85 - 0.3),
    y = -3.2 + 0.95 * (-0.4 - -3.2),
    label = "PF\nr = -0.664*", # removed QR
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

I.qp2 # Display the plot

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # light v anc: -0.1118 (0.8405)

# Nitrogen: all logged ----------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "N", 'nit', 'other')) # for testing significance of matching evolutionary conditions.

df.filt <- df %>% 
  mutate(
    z.y = log(N.comp),
    z.x = log(N.µ.max)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # Looks good

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Quantile regression ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # nit: r -0.4356368, p 3.135e-05

sma.N <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.N)      # slope -7.337128,
coef(sma.N) 

N.qp2 <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  
  geom_point(
    data = subset(df.filt, evol.bin =="other" & pareto.opt == "N"),
    shape = 16, 
    size = 2,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="other" & pareto.opt == "Y"),
    shape = 16, 
    size = 4,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="ancestral" & pareto.opt == "N"),
    shape = 5, 
    size = 2,
    colour = "black",
    stroke = 1,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="ancestral" & pareto.opt == "Y"),
    shape = 5, 
    size = 4,
    colour = "black",
    stroke = 1,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin == "nit" & pareto.opt == "N"),
    shape = 21, 
    size = 2,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("plum3", 0.6)
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin == "nit" & pareto.opt == "Y"),
    shape = 21, 
    size = 4,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("plum3", 0.6)
  ) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 0.6, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_abline(intercept = coef(sma.N)[1], slope = coef(sma.N)[2], lwd = 0.6, linetype = "dashed") +
  
  labs(x = expression(atop("log(Maximum growth rate)",
                           log(italic(μ)[italic(max)]))),
       y = expression(atop("log(Low-N tolerance)",
                           log(1/italic(N)^"*"))),
       color = "Evolutionary History",
       title = "B) Nitrogen") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "chocolate3",
               "Biotic depletion x Salt" = "skyblue",
               "Control" = "olivedrab4",
               "Light limitation" = "goldenrod2",
               "Nitrogen limitation" = "plum3",
               "Ancestral" = "black",
               "Phosphorus limitation" = "brown4",  
               "Salt stress" = "navyblue")
  ) +
  
  ylim(-2.75, 0.35) +
  xlim(0.25, 0.65) +
  
  annotate(
    "text",
    x = 0.25 + 0.6 * (0.65 - 0.25),
    y = -2.75 + 0.95 * (0.35 - -2.75),
    label = "PF\nr = -0.436*", # removed QR
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

N.qp2  # Display the plot

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # nit v anc: -0.547 (0.0807), phos: -0.658 (0.0236)

# Phosphorus: all logged -------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "P", 'phos', 'other')) # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = log(P.comp),
    z.x = log(P.µ.max)
  ) # Specify the x and y variables and their 95% CIs

plot(z.y ~ z.x, data= df.filt) # Looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Quantile regression ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # phos: r -0.2113107, p 0.009935

sma.P <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.P)      # slope -8.151480
coef(sma.P)   

P.qp2 <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  
  geom_point(
    data = subset(df.filt, evol.bin =="other" & pareto.opt == "N"),
    shape = 16, 
    size = 2,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="other" & pareto.opt == "Y"),
    shape = 16, 
    size = 4,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="ancestral" & pareto.opt == "N"),
    shape = 5, 
    size = 2,
    colour = "black",
    stroke = 1,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="ancestral" & pareto.opt == "Y"),
    shape = 5, 
    size = 4,
    colour = "black",
    stroke = 1,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin == "phos" & pareto.opt == "N"),
    shape = 21, 
    size = 2,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("brown4", 0.6)
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin == "phos" & pareto.opt == "Y"),
    shape = 21, 
    size = 4,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("brown4", 0.6)
  ) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 0.6, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_abline(intercept = coef(sma.P)[1], slope = coef(sma.P)[2], lwd = 0.6, linetype = "dashed") +
  
  labs(x = expression(atop("log(Maximum growth rate)",
                           log(italic(μ)[italic(max)]))),
       y = expression(atop("log(Low-P tolerance)",
                           log(1/italic(P)^"*"))),
       color = "Evolutionary History",
       title = "C) Phosphorus") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "chocolate3",
               "Biotic depletion x Salt" = "skyblue",
               "Control" = "olivedrab4",
               "Light limitation" = "goldenrod2",
               "Nitrogen limitation" = "plum3",
               "Ancestral" = "black",
               "Phosphorus limitation" = "brown4",  
               "Salt stress" = "navyblue")
  ) +
  
  ylim(-1.5, 1.85) +
  xlim(0.2, 0.6) +
  
  annotate(
    "text",
    x = 0.2 + 0.6 * (0.6 - 0.2),
    y = 1.5 + 0.95 * (1.85 - 1.5),
    label = "PF\nEvo\nr = -0.211*", # removed QR
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

P.qp2  # Display the plot

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # phos v anc: -0.962 (0.0002), BS: -0.860 (0.0022)

# Salt: logged all --------------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol %in% c("S", "BS"), 'salt', 'other')) # Ss and BSs are treated as equivalent

df.filt <- df %>% 
  mutate(
    z.y = log(S.c),
    z.x = log(S.µ.max)
  ) # Specify the x and y variables and their 95% CIs

plot(z.y ~ z.x, data= df.filt) # Looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # salt: r --0.1941612, p 0.01805

sma.S <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.S)      # slope -2.446888
coef(sma.S)    

S.qp2 <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  
  geom_point(
    data = subset(df.filt, evol.bin =="other" & pareto.opt == "N"),
    shape = 16, 
    size = 2,
    alpha = 0.21
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="other" & pareto.opt == "Y"),
    shape = 16, 
    size = 4,
    alpha = 0.21
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="ancestral" & pareto.opt == "N"),
    shape = 5, 
    size = 2,
    colour = "black",
    stroke = 1,
    alpha = 0.21
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="ancestral" & pareto.opt == "Y"),
    shape = 5, 
    size = 4,
    colour = "black",
    stroke = 1,
    alpha = 0.21
  ) +
  
  geom_point(
    data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & pareto.opt == "N"),
    shape = 21, 
    size = 2,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("skyblue", 0.21)
  ) +
  
  geom_point(
    data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & pareto.opt == "Y"),
    shape = 21, 
    size = 4,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("skyblue", 0.21)
  ) +
  
  geom_point(
    data = subset(df.filt, Evol.plt == "Salt stress" & pareto.opt == "N"),
    shape = 21, 
    size = 2,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("navyblue", 0.21)
  ) +
  
  geom_point(
    data = subset(df.filt, Evol.plt == "Salt stress" & pareto.opt == "Y"),
    shape = 21, 
    size = 4,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("navyblue", 0.21)
  ) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 0.6, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_abline(intercept = coef(sma.S)[1], slope = coef(sma.S)[2], lwd = 0.6, linetype = "dashed") +
  
  labs(x = expression(atop("log(Maximum growth rate)",
                           log(italic(μ)[italic(max)]))),
       y = expression(atop("log(Salt tolerance)",
                           log(italic(S)))),
       color = "Evolutionary History",
       title = "D) Salt") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "chocolate3",
               "Biotic depletion x Salt" = "skyblue",
               "Control" = "olivedrab4",
               "Light limitation" = "goldenrod2",
               "Nitrogen limitation" = "plum3",
               "Ancestral" = "black",
               "Phosphorus limitation" = "brown4",  
               "Salt stress" = "navyblue")
  ) +
  
  ylim(0.5, 2.3) +
  xlim(-0.15, 0.55) +
  
  annotate(
    "text",
    x = -0.15 + 0.6 * (0.55 - -0.15),
    y = 0.5 + 0.95 * (2.3 - 0.5),
    label = "PF\nEvo\nr = -0.194*",
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.03)# theme stuff
  )

S.qp2  # Display the plot

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0.011

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # S v anc: -0.9990 (<0.0001), BS: -0.8946 (<0.0001)

# Light ~ Nitrogen -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = log(I.comp),
    z.x = log(N.comp)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0.545

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # L v anc: 0.0485 (0.9974), N: -0.4594 (0.2666), P sig (-1.5699, <0.0001)

# Light ~ Phosphorus -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = log(I.comp),
    z.x = log(P.comp)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0.013

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # L v anc: 0.73675 (1), P: -1.0246 (0.0001)

# Light ~ Salt -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = log(I.comp),
    z.x = log(S.c)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0.337

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # L v anc: 0.191 (1), S: -1.511 (<0.0001), BS: -1.394 (<0.0001), P too

# Light ~ Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = log(I.comp),
    z.x = log(T.br)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # L v anc: 0.3473 (1)

# Nitrogen ~ Phosphorus -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = log(N.comp),
    z.x = log(P.comp)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0.07

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # nit v anc: -0.2925 (0.6157), P: -1.5092 (<0.0001)

# Nitrogen ~ Salt -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = log(N.comp),
    z.x = log(S.c)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0.025

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # nit v anc: 0.0582 (0.9987), S: -1.1134 (<0.0001), BS: -1.1924 (<0.0001), P too

# Nitrogen ~ Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = log(N.comp),
    z.x = log(T.br)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF
fit <- lm(z.y~z.x, data = par.res.1)

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0.246

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # nit v anc: -0.451 (0.342)

# Phosphorus ~ Salt -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = log(P.comp),
    z.x = log(S.c)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0.086

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # p v anc: -0.991 (<0.0001), S -1.633 (<0.0001) BS: -1.887 (<0.0001)

# Phosphorus ~ Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = log(P.comp),
    z.x = log(T.br)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0.021

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # p v anc: 0.4531 (1)

# Salt ~ Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = log(S.c),
    z.x = log(T.br)
  ) # Specify the x and y variables

plot(z.y ~ z.x, data= df.filt) # looks good!

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

pf.ids <- par.res.1 %>% # extract Id labels for the Pareto-optimal points
  distinct(rep.ID) %>%
  mutate(pareto.opt = "Y")  

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt)) # Add a 'Y' to pareto-optimal points in the corresponding new column

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n())

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){

  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),
    )

  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>%
    arrange(z.x.sim)

  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)

  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)

  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)

  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF = nrow(par.res.n)
  ))

}

mean(null.df3$a.emp.n >= a.emp) # 0

###### Evolutionary optimization test ######

df.pf.segments <- par.res.1 %>%
  arrange(z.x2) %>%
  transmute(
    ax = z.x2,
    ay = z.y2,
    bx = lead(z.x2),
    by = lead(z.y2)
  ) %>%
  filter(!is.na(bx), !is.na(by))

df.filt <- df.filt %>%
  rowwise() %>%
  mutate(
    nearest = list(nearest_pf_segment(z.x2, z.y2))
  ) %>%
  tidyr::unnest(nearest) %>%
  ungroup()

mod <- lmer(dist.to.pf ~ Evol.plt + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<") # s v anc: -0.778 (0.0097), BS: -0.977 (0.0017)

# Compile & save the figures -----------------------------------------------

###### Figure 2, logged ######

plots.qp <- list(I.qp2, N.qp2, P.qp2, S.qp2, T.qp2)

legend_df3 <- data.frame(
  x = c(1, 2, 1, 2, 1, 2, 1, 2),
  y = c(1, 1, 2, 2, 1, 1, 2, 2),
  Group = factor(c("Ancestral", "Other", "Matching", "Matching", "Ancestral", "Other", "Matching", "Matching")),
  Group2 = factor(c("Biotic depletion", "Biotic depletion x Salt", "Control", "Light limitation", "Nitrogen limitation", "Ancestral", "Phosphorus limitation", "Salt stress")),
  LineType = factor(c("Pareto front", "SMA regression", "Pareto front", "SMA regression", "Pareto front", "SMA regression", "Pareto front", "SMA regression"))
)

legend_plot3 <- ggplot() +
  
  geom_point(
    data = legend_df3,
    aes(x = x, y = y, colour = Group2),
    shape = 16,
    size = 2,
    alpha = 0.6
  ) +
  
  # Shape legend: Ancestral + Matching + Sub-optimal + Pareto-optimal
  geom_point(
    data = data.frame(
      x     = c(1, 1, 1, 1),
      y     = c(1, 1, 1, 1),
      Group = factor(c("Ancestral", "Matching", "Sub-optimal", "Pareto-optimal"),
                     levels = c("Ancestral", "Matching", "Sub-optimal", "Pareto-optimal"))
    ),
    aes(x = x, y = y, shape = Group),
    colour = "black",
    size = 2,
    stroke = 1,
    alpha = 0.6
  ) +
  
  geom_line(
    data = legend_df3,
    aes(x = x, y = y, linetype = LineType),
    colour = "black",
    size = 0.6,
    show.legend = TRUE
  ) +
  
  scale_shape_manual(
    name = NULL,
    values = c("Ancestral" = 5, "Matching" = 21, "Sub-optimal" = 21, "Pareto-optimal" = 21)
  ) +
  
  scale_color_manual(
    name = "Evolution environment",
    values = c(
      "Biotic depletion"        = "chocolate3",
      "Biotic depletion x Salt" = "skyblue",
      "Control"                 = "olivedrab4",
      "Light limitation"        = "goldenrod2",
      "Nitrogen limitation"     = "plum3",
      "Phosphorus limitation"   = "brown4",
      "Salt stress"             = "navyblue"
    )
  ) +
  
  scale_linetype_manual(
    name = "Line",
    values = c(
      "Pareto front" = "solid",
      "SMA regression" = "dashed"
    ),
    labels = c(
      "Pareto front" = "Pareto front",
      "SMA regression" = "SMA regression")
  ) +
  
  guides(
    colour = guide_legend(
      order = 1,
      override.aes = list(
        shape    = 16,
        size     = 2,
        alpha    = 0.6,
        linetype = 0
      )
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        # One value per legend row: Ancestral, Matching, Sub-optimal, Pareto-optimal
        shape    = c(5,       21,      21,        21),
        colour   = c("black", "black", NA,        NA),
        fill     = c(NA,      "grey75", "grey75",  "grey75"),
        size     = c(2,       2,       2.5,         4),
        stroke   = c(1,       1,       0,         0),
        alpha    = c(0.6,     1,       1,         1),
        linetype = 0
      )
    ),
    linetype = guide_legend(
      order = 3,
      override.aes = list(
        colour    = "black",
        linewidth = 0.6,
        shape     = NA
      )
    )
  ) +
  
  theme_void() +
  theme(
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "lines"),
    legend.key.width = unit(1.2, "lines")
  )

legend_only3 <- get_legend(legend_plot3)

all_plots_qp <- c(plots.qp, list(legend_only3))

qp_toffs <- plot_grid(plotlist = all_plots_qp,
                      ncol = 2,
                      align = "hv")

ggsave("figures-misc/01.1_fig_2b_log_td_intra-gradient_tradeoffs.jpeg", qp_toffs, width = 5.5, height = 8.25)
