# Jason R Laurich

# May 26th, 2026

# edited for resubmission on August 16th
  # quantile regressions replaced with Pearson correlations and SMA plotting. 

# This script will assess the existence of Pareto fronts within abiotic gradients, through the implementation of Li et al.
# (2019)'s analysis. It will also evaluate whether experimental evolution shifted populations towards Pareto-optimal trait values

# It will generate Figure 2. 

# Inputs: 27_summary_table.csv
# Outputs: in figures-main : 02_fig_2_intra-gradient_tradeoffs.pdf
# in figures-supp: 03_fig_s3_evol_opt_test.pdf

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(sp) # For the point.in.polygon function
library(scam)
library(pracma) # For calculating polygon area
library(cowplot)
library(lme4)
library(emmeans)
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

# Temperature -------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral", "evolved") # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = T.br,
    z.x = T.µ.max
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

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # Temp: r -0.6641438, p < 2.2e-16

sma.T <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.T)      # slope -2.677501
coef(sma.T)      

T.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  
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
  
  geom_abline(intercept = coef(sma.T)[1], slope = coef(sma.T)[2],
              lwd = 0.6, linetype = "dashed") +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 0.6, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Thermal breadth",
                           italic("T")[italic(br)] * " (°C)")),
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
  
  ylim(14, 22) +
  xlim(2.25, 5.25) + 
  
  annotate(
    "text",
    x = 2.25 + 0.6 * (5.25 - 2.25),
    y = 14 + 0.95 * (22 - 14),
    label = "PF\nr = -0.664*",                          # removing QR
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
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)# theme stuff
  )

T.qp # Display the plot

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

contrast(em, method = "trt.vs.ctrl", side = "<") 

# Light -------------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "L", 'light', 'other')) # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = I.comp,
    z.x = I.µ.max
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

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # Light: r -0.5029304, p 1.142e-11

sma.L <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.L)      # slope -0.3761067
coef(sma.L)      

I.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  
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
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Low-L tolerance",
                           1/italic(L)^"*" ~ (mu*mol^-1 ~ m^2 ~ s))),
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
  
  ylim(0, 0.6) +
  xlim(1.4, 2.35) +
  
  annotate(
    "text",
    x = 1.4 + 0.6 * (2.35 - 1.4),
    y = 0.0 + 0.95 * (0.60 - 0.0),
    label = "PF\nr = -0.503*", # removed QR
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
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)# theme stuff
  )

I.qp # Display the plot

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

contrast(em, method = "trt.vs.ctrl", side = "<") # light v anc: -0.1216 (0.676)

# Nitrogen ----------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "N", 'nit', 'other')) # for testing significance of matching evolutionary conditions.

df.filt <- df %>% 
  mutate(
    z.y = N.comp,
    z.x = N.µ.max
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

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # nit: r -0.3169419, p 8.682e-05

sma.N <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.N)      # slope -1.754382,
coef(sma.N)         

N.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  
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
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Low-N tolerance",
                           1/italic(N)^"*" ~ (mu*mol^-1))),
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
  
  ylim(0, 1.35) +
  xlim(1.275, 1.95) +
  
  annotate(
    "text",
    x = 1.275 + 0.6 * (1.95 - 1.275),
    y = 0.0 + 0.95 * (1.35 - 0.0),
    label = "PF\nr = -0.317*", # removed QR
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
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)# theme stuff
  )

N.qp  # Display the plot

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

contrast(em, method = "trt.vs.ctrl", side = "<") # nit v anc: -0.466 (0.1965)

# Phosphorus -------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "P", 'phos', 'other')) # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = P.comp,
    z.x = P.µ.max
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

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # phos: r -0.2892904, p 0.0003623

sma.P <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.P)      # slope -9.269761
coef(sma.P)        

P.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  
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
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Low-P tolerance",
                           1/italic(P)^"*" ~ (mu*mol^-1))),
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
  
  ylim(0, 6) +
  xlim(1.2, 1.8) +
  
  annotate(
    "text",
    x = 1.2 + 0.6 * (1.8 - 1.2),
    y = 0.0 + 0.95 * (6 - 0.0),
    label = "PF\nEvo\nr = -0.289*", # removed QR
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
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)# theme stuff
  )

P.qp  # Display the plot

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

contrast(em, method = "trt.vs.ctrl", side = "<") # phos v anc: -0.927 (0.0001), SB -0.938 (0.0004)

# Supplemental figure 2, panel A: full data with Pareto front

supp2a <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  
  
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
  
  labs(x = expression("Maximum growth rate, " * italic("\u03bc")[italic(max)] * " (day"^-1 * ")"),
       y = expression("Low-P tolerance, 1/" * italic(P) * "* (" * mu * "mol"^-1 * ")"),
       color = "Evolutionary History",
       title = "A) Pareto front with full data") +  # labels
  
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
  
  ylim(0, 6) +
  xlim(1.2, 1.8) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)# theme stuff
  )

supp2a  # Display the plot

# Supplement figure 2, panel B

df.B <- df.filt %>%
  filter(evol.bin %in% c("ancestral", "phos"))

pf.pts <- par.res.1 %>% 
  arrange(z.x2) %>%
  select(z.x2, z.y2) # PF line in z-score space (piecewise linear)

supp2b <- ggplot() +
  
  geom_line(data = pf.pts,       # Piecewise linear Pareto front
            aes(x = z.x2, y = z.y2),
            color = "goldenrod2", 
            linewidth = 1,
            linetype = "31") +
  
  geom_segment(data = df.B, # Connecting lines from each point to its nearest location on the PF
               aes(x = z.x2, y = z.y2,
                   xend = x.pf, yend = y.pf,
                   color = evol.bin),
               linewidth = 0.8,
               linetype = "31",
               alpha = 0.6) +
  
  geom_point(
    data = subset(df.filt, evol.bin == "phos" & pareto.opt == "N"),
    aes(x = z.x2, y = z.y2),
    shape = 21, 
    size = 2,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("brown4", 0.6)
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin == "phos" & pareto.opt == "Y"),
    aes(x = z.x2, y = z.y2),
    shape = 21, 
    size = 4,
    stroke = 1,
    colour = "black",
    fill = scales::alpha("brown4", 0.6)
  ) +

  geom_point(
    data = subset(df.filt, evol.bin =="ancestral" & pareto.opt == "N"),
    aes(x = z.x2, y = z.y2),
    shape = 5, 
    size = 2,
    colour = "black",
    stroke = 1,
    alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.filt, evol.bin =="ancestral" & pareto.opt == "Y"),
    aes(x = z.x2, y = z.y2),
    shape = 5, 
    size = 4,
    colour = "black",
    stroke = 1,
    alpha = 0.6
  ) +
  
  scale_color_manual(values = c("ancestral" = "black", "phos" = "brown4"),
                     guide = "none") +
  
  labs(x = expression("Maximum growth rate, " * italic("\u03bc")[italic(max)] * " (z-score)"),
       y = expression("Low-P tolerance, 1/" * italic(P) * "* (z-score)"),
       title = "B) distances to Pareto front") +
  
  theme_classic() +
  theme(
    axis.title = element_text(size = 10),
    axis.text  = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  )

supp2b

x_min <- -1.5; x_max <- 1.5
y_min <- -0.25; y_max <- 2.75

df.B.inset <- df.B %>%
  filter(z.x2 >= x_min, z.x2 <= x_max,
         z.y2 >= y_min, z.y2 <= y_max)

inset <- ggplot() +
  
  geom_line(data = pf.pts,
            aes(x = z.x2, y = z.y2),
            color = "goldenrod2",
            linewidth = 1,
            linetype = "31") +
  
  geom_segment(data = df.B.inset,
               aes(x = z.x2, y = z.y2,
                   xend = x.pf, yend = y.pf,
                   color = evol.bin),
               linewidth = 0.6,
               linetype = "31",
               alpha = 0.8) +
  
  geom_point(
    data = subset(df.B.inset, evol.bin == "phos" & pareto.opt == "N"),
    aes(x = z.x2, y = z.y2),
    shape = 21, size = 1.5, stroke = 1,
    colour = "black", fill = scales::alpha("brown4", 0.6)
  ) +
  
  geom_point(
    data = subset(df.B.inset, evol.bin == "phos" & pareto.opt == "Y"),
    aes(x = z.x2, y = z.y2),
    shape = 21, size = 3, stroke = 1,
    colour = "black", fill = scales::alpha("brown4", 0.6)
  ) +
  
  geom_point(
    data = subset(df.B.inset, evol.bin == "ancestral" & pareto.opt == "N"),
    aes(x = z.x2, y = z.y2),
    shape = 5, size = 1.5, stroke = 0.8,
    colour = "black", alpha = 0.6
  ) +
  
  geom_point(
    data = subset(df.B.inset, evol.bin == "ancestral" & pareto.opt == "Y"),
    aes(x = z.x2, y = z.y2),
    shape = 5, size = 3, stroke = 0.8,
    colour = "black", alpha = 0.6
  ) +
  
  scale_color_manual(values = c("ancestral" = "black", "phos" = "brown4"),
                     guide = "none") +
  
  coord_fixed(ratio = 1, xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
  
  theme_classic() +
  theme(
    axis.title       = element_blank(),
    axis.text        = element_text(size = 6),
    plot.background  = element_rect(color = "black", linewidth = 0.5, fill = "white"),
    panel.background = element_rect(fill = "white")
  )

inset

panel.B.inset <- ggdraw(supp2b) +
  draw_plot(inset, x = 0.54, y = 0.54, width = 0.39, height = 0.39)

panel.B.inset

# Supplement figure 2, panel C

evol_order <- c("Ancestral", "Phosphorus limitation",
                "Light limitation", "Nitrogen limitation", "Salt stress",
                "Biotic depletion x Salt", "Biotic depletion",
                "Control"
                )

df.filt <- df.filt %>%
  mutate(
    Evol.plt  = factor(Evol.plt, levels = evol_order),
    highlight = ifelse(evol.bin %in% c("ancestral", "phos"), "yes", "no")
  )

# Set bracket height just above max dist.to.pf
y.max  <- max(df.filt$dist.to.pf, na.rm = TRUE)
y.br   <- y.max * 1.05   # bracket height
y.text <- y.max * 1.10   # p-value text height

supp2c <- ggplot(df.filt, aes(x = Evol.plt, y = dist.to.pf,
                              fill = Evol.plt, alpha = highlight)) +
  
  geom_boxplot(outlier.shape = NA, colour = "grey30") +
  
  geom_jitter(aes(colour = Evol.plt, alpha = highlight),
              width = 0.2, size = 1.5, shape = 16) +
  
  # p-value bracket between Ancestral (x=1) and Phosphorus limitation (x=2)
  annotate("segment", x = 1, xend = 2, y = y.br, yend = y.br,
           colour = "black", linewidth = 0.5) +
  annotate("segment", x = 1, xend = 1, y = y.br, yend = y.br * 0.97,
           colour = "black", linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = y.br, yend = y.br * 0.97,
           colour = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = y.text,
           label = expression(italic(P) * " < 0.001"), size = 3.5) +
  
  scale_fill_manual(values = c(
    "Ancestral"               = "grey40",
    "Phosphorus limitation"   = "brown4",
    "Biotic depletion"        = "chocolate3",
    "Biotic depletion x Salt" = "skyblue",
    "Control"                 = "olivedrab4",
    "Light limitation"        = "goldenrod2",
    "Nitrogen limitation"     = "plum3",
    "Salt stress"             = "navyblue"
  ), guide = "none") +
  
  scale_colour_manual(values = c(
    "Ancestral"               = "grey40",
    "Phosphorus limitation"   = "brown4",
    "Biotic depletion"        = "chocolate3",
    "Biotic depletion x Salt" = "skyblue",
    "Control"                 = "olivedrab4",
    "Light limitation"        = "goldenrod2",
    "Nitrogen limitation"     = "plum3",
    "Salt stress"             = "navyblue"
  ), guide = "none") +
  
  scale_alpha_manual(values = c("yes" = 0.6, "no" = 0.30), guide = "none") +
  
  labs(x = NULL,
       y = "Distance to Pareto front (z-score)",
       title = "C) evolutionary optimization test") +
  
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    axis.title.y = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    plot.title   = element_text(size = 10, face = "bold", hjust = 0.03)
  )

supp2c

supp2 <- plot_grid(supp2a, panel.B.inset, supp2c, 
                   nrow = 1,
                   rel_widths = c(1, 1, 1))

ggsave("figures-supplemental/03_fig_s3_evol_opt_test.pdf", supp2, width = 12, height = 5)

# Salt --------------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol %in% c("S", "BS"), 'salt', 'other')) # Ss and BSs are treated as equivalent

df.filt <- df %>% 
  mutate(
    z.y = S.c,
    z.x = S.µ.max
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

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # salt: r -0.1751515, p 0.03324

sma.S <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.S)      # slope -9.198362
coef(sma.S)         

S.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  
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
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Salt tolerance",
                           italic(S) ~ " (g/L)")),
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
  
  ylim(1, 9.5) +
  xlim(0.85, 1.75) +
  
  annotate(
    "text",
    x = 0.85 + 0.6 * (1.75 - 0.85),
    y = 1 + 0.95 * (9.5 - 1),
    label = "PF\nEvo\nr = -0.175*", # removed QR
    hjust = 0,
    vjust = 0.63,
    size = 3,
    fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 10, face = "plain"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)# theme stuff
  )

S.qp  # Display the plot

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

contrast(em, method = "trt.vs.ctrl", side = "<") # salt v anc: -1.0371 (<0.0001), SB -0.9409 (<0.0001), ctl: -1.219 (<0.0001)

# Compile & save Figure 2 -----------------------------------------------

plots.qp <- list(I.qp, N.qp, P.qp, S.qp, T.qp)

legend_df <- data.frame(
  x = c(1, 2, 1, 2, 1, 2, 1, 2),
  y = c(1, 1, 2, 2, 1, 1, 2, 2),
  Group = factor(c("Ancestral", "Other", "Matching", "Matching", "Ancestral", "Other", "Matching", "Matching")),
  Group2 = factor(c("Biotic depletion", "Biotic depletion x Salt", "Control", "Light limitation", "Nitrogen limitation", "Ancestral", "Phosphorus limitation", "Salt stress")),
  LineType = factor(c("Pareto front", "SMA regression", "Pareto front", "SMA regression", "Pareto front", "SMA regression", "Pareto front", "SMA regression"))
)

legend_plot <- ggplot() +
  
  geom_point(
    data = legend_df,
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
    data = legend_df,
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

legend_only <- get_legend(legend_plot)

all_plots_qp <- c(plots.qp, list(legend_only))

qp_toffs <- plot_grid(plotlist = all_plots_qp,
                      ncol = 2,
                      align = "hv")

ggsave("figures-main/02_fig_2_intra-gradient_tradeoffs.pdf", qp_toffs, width = 5.5, height = 8.25)
