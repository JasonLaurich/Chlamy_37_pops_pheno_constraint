# Jason R Laurich

# June 2026

# edited for re-submission on August 16th, 2026. 

# This script assesses Pareto fronts between niche-determining ecological traits (I*, N*, P*, S, Tbr)
# and secondary biological traits (N content, P content, pigmentation, biovolume). It implements
# Li et al. (2019) polygonal empty space analysis and evaluates whether experimental evolution shifted
# populations towards Pareto-optimal trait values using a distance-to-PF mixed model approach
# (lmer + emmeans).

# It generates Figure S7.

# Inputs: 27_summary_table.csv
# Outputs: figures-supplemental/07_fig_s7_extra_inter-gradient_toffs.pdf

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(sp)
library(scam)
library(pracma)
library(cowplot)
library(lme4)
library(emmeans)
library(smatr)

par_frt <- function(df, xvar, yvar) {
  df <- df[order(-df[[xvar]], df[[yvar]]), ]
  pareto_points <- df[1, ]
  for (i in 2:nrow(df)) {
    if (df[i, yvar] >= tail(pareto_points[[yvar]], 1)) {
      pareto_points <- rbind(pareto_points, df[i, ])
    }
  }
  return(pareto_points)
}

closest_on_segment <- function(px, py, ax, ay, bx, by) {
  apx <- px - ax; apy <- py - ay
  abx <- bx - ax; aby <- by - ay
  ab2 <- abx^2 + aby^2
  t <- (apx * abx + apy * aby) / ab2
  t <- pmax(0, pmin(1, t))
  x.closest <- ax + t * abx
  y.closest <- ay + t * aby
  dist <- sqrt((px - x.closest)^2 + (py - y.closest)^2)
  data.frame(x.pf = x.closest, y.pf = y.closest, dist.to.pf = dist)
}

nearest_pf_segment <- function(px, py) {
  df.pf.segments %>%
    rowwise() %>%
    mutate(closest = list(closest_on_segment(px, py, ax, ay, bx, by))) %>%
    tidyr::unnest(closest) %>%
    ungroup() %>%
    slice_min(dist.to.pf, n = 1, with_ties = FALSE) %>%
    select(x.pf, y.pf, dist.to.pf)
}


# Shared colour palette ---------------------------------------------------

evol_colours <- c(
  "Biotic depletion"        = "chocolate3",
  "Biotic depletion x Salt" = "skyblue",
  "Control"                 = "olivedrab4",
  "Light limitation"        = "goldenrod2",
  "Nitrogen limitation"     = "plum3",
  "Ancestral"               = "black",
  "Phosphorus limitation"   = "brown4",
  "Salt stress"             = "navyblue"
)


# Load & examine the data -------------------------------------------------

df <- read.csv("processed-data/27_summary_table.csv")
head(df)

df$Evol.plt <- factor(df$Evol,
                      levels = c("none", "L", "N", "P", "S", "B", "BS", "C"),
                      labels = c("Ancestral", "Light limitation", "Nitrogen limitation",
                                 "Phosphorus limitation", "Salt stress",
                                 "Biotic depletion", "Biotic depletion x Salt", "Control"))

df$Anc.plt <- factor(df$Anc,
                     levels = c("anc2", "anc3", "anc4", "anc5", "cc1690"),
                     labels = c("Ancestor 2", "Ancestor 3", "Ancestor 4",
                                "Ancestor 5", "Mixed ancestry"))


# Axis label expressions --------------------------------------------------

xlab.I   <- expression(atop("Low-L tolerance",
                             1/italic(L)^"*" ~ (mu*mol^-1 ~ m^2 ~ s)))
xlab.N   <- expression(atop("Low-N tolerance",
                             1/italic(N)^"*" ~ (mu*mol^-1)))
xlab.P   <- expression(atop("Low-P tolerance",
                             1/italic(P)^"*" ~ (mu*mol^-1)))
xlab.S   <- expression(atop("Salt tolerance",
                             italic(S) ~ "(g/L)"))
xlab.T   <- expression(atop("Thermal breadth",
                             italic(T)[italic(br)] ~ "(°C)"))

xlab.Nc  <- expression(atop("Nitrogen content",
                             paste("(", mu, "g/L)")))
xlab.Pc  <- expression(atop("Phosphorus content",
                             paste("(", mu, "g/L)")))
xlab.Pig <- expression(atop("Pigmentation", "PC1"))
xlab.Bv  <- expression(atop("Biovolume",
                             paste("(", mu, m^3, ")")))


# N content comparisons ---------------------------------------------------

# A — N content v light ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "L", "light", "other"))

df.filt <- df %>% mutate(z.y = I.comp, z.x = mean.N.µg.l)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # light v nit cont: r 0.01364102, p 0.8693

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)      # slope 0.0005833191,
coef(sma)         

NcL.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "light" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "light" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  labs(x = xlab.Nc, y = xlab.I, title = "A) Light—N cont.") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) + xlim(145, 635) +
  
  annotate(
    "text",
    x = 145 + 0.6 * (635 - 145),
    y = 0.0 + 0.95 * (0.60 - 0.0),
    label = "r = 0.0136", # removed QR
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

NcL.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<") # light v anc: -0.00737 (0.9907), salt: -0.62115 (0.0250)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.660

# E — N content v nitrogen ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "N", "nit", "other"))

df.filt <- df %>% mutate(z.y = N.comp, z.x = mean.N.µg.l)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # nit v nit cont: r 0.1136184, p 0.1691

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)    
coef(sma)         

NcN.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  labs(x = xlab.Nc, y = xlab.N, title = "E) Nitrogen—N cont.") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 1.35) + xlim(145, 635) +
  
  annotate(
    "text",
    x = 145 + 0.6 * (635 - 145),
    y = 0.0 + 0.95 * (1.35 - 0.0),
    label = "r = 0.114", # removed QR
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

NcN.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # nit v anc: -0.360 (0.5529)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.235

# I — N content v phosphorus ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "P", "phos", "other"))

df.filt <- df %>% mutate(z.y = P.comp, z.x = mean.N.µg.l)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # phos v nit cont: r 0.1212686, p 0.142

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)    
coef(sma)         

NcP.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  labs(x = xlab.Nc, y = xlab.P, title = "I) Phosphorus—N cont.") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 6) + xlim(145, 635) +
  
  annotate(
    "text",
    x = 145 + 0.6 * (635 - 145),
    y = 0.0 + 0.95 * (6 - 0.0),
    label = "r = 0.121", # removed QR
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

NcP.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # phos v anc: -0.3636 (0.517)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.655

# M — N content v salt ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol %in% c("S", "BS"), "salt", "other"))  # Salt stress and Biotic depletion x Salt treated as equivalent

df.filt <- df %>% mutate(z.y = S.c, z.x = mean.N.µg.l)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # salt v nit cont: r 0.1163256, p 0.1591

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)     
coef(sma)         

NcS.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, Evol == "S" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "S" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "BS" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "BS" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  annotate("text",
           x     = 145 + 0.25 * (635 - (145)),
           y     = 1 + 0.95 * (9.5 - (1)),
           label = "Evo\nr = 0.116",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +

  labs(x = xlab.Nc, y = xlab.S, title = "M) Salt—N cont.") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(1, 9.5) + xlim(145, 635) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

NcS.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # Salt v anc: -1.5120 (<0.0001), SB: -0.3736 (0.4767)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.132

# Q — N content v temperature ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral", "evolved")

df.filt <- df %>% mutate(z.y = T.br, z.x = mean.N.µg.l)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # temp v nit cont: r 0.003053605, p 0.9706

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)    
coef(sma)         

NcT.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "evolved" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "evolved" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  labs(x = xlab.Nc, y = xlab.T, title = "Q) Temperature—N cont.") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(14, 22) + xlim(145, 635) +
  
  annotate("text",
           x     = 145 + 0.6 * (635 - (145)),
           y     = 14 + 0.95 * (22 - 14),
           label = "r = 0.00305",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

NcT.qp

# Evolutionary optimization test (evo opt not applicable: no temperature evolution treatment)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.068

# P content comparisons ---------------------------------------------------

# B — P content v light ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "L", "light", "other"))

df.filt <- df %>% mutate(z.y = I.comp, z.x = mean.P.µg.l)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # light v phos cont: r 0.05666032, p 0.494

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)      
coef(sma)         

PcL.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "light" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "light" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  labs(x = xlab.Pc, y = xlab.I, title = "B) Light—P cont.") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) + xlim(30, 175) +
  
  annotate("text",
           x     = 30 + 0.6 * (175 - 30),
           y     = 0 + 0.95 * (0.6 - 0),
           label = "r = 0.0567",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

PcL.qp

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

mod <- lmer(dist.to.pf ~ Evol.plt2+ (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt)

contrast(em, method = "trt.vs.ctrl", side = "<")  # light v anc: + -0.2060 (0.8120), s: -1.512 (<0.0001)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.731

# F — P content v nitrogen ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "N", "nit", "other"))

df.filt <- df %>% mutate(z.y = N.comp, z.x = mean.P.µg.l)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # nit v nit cont: r 0.1904284, p 0.02043

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)    
coef(sma)         

PcN.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_abline(intercept = coef(sma)[1], slope = coef(sma)[2], lwd = 0.6, linetype = "dashed") +
  
  labs(x = xlab.Pc, y = xlab.N, title = "F) Nitrogen—P cont.") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 1.35) + xlim(30, 175) +
  
  annotate("text",
           x     = 30 + 0.6 * (175 - 30),
           y     = 0 + 0.95 * (1.35 - 0),
           label = "r = 0.190*",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +

  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

PcN.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # nit v anc: - 0.2276 (0.7917), P: -1.1112 (0.0002)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.102

# J — P content v phosphorus ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "P", "phos", "other"))

df.filt <- df %>% mutate(z.y = P.comp, z.x = mean.P.µg.l)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # phos v phos cont: r 0.3196471, p 7.492e-05

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)    
coef(sma)         

PcP.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_abline(intercept = coef(sma)[1], slope = coef(sma)[2], lwd = 0.6, linetype = "dashed") +

  labs(x = xlab.Pc, y = xlab.P, title = "J) Phosphorus—P cont.") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 6) + xlim(30, 175) +
  
  annotate("text",
           x     = 30 + 0.6 * (175 - 30),
           y     = 0 + 0.95 * (6 - 0),
           label = "r = 0.320*",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

PcP.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # phos v anc: - 0.3806 (0.5483)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.896

# N — P content v salt ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol %in% c("S", "BS"), "salt", "other"))  # Salt stress and Biotic depletion x Salt treated as equivalent

df.filt <- df %>% mutate(z.y = S.c, z.x = mean.P.µg.l)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # nit v nit cont: r 0.08881984, p 0.283

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)  
coef(sma)         

PcS.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, Evol == "S" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "S" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "BS" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "BS" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +

  labs(x = xlab.Pc, y = xlab.S, title = "N) Salt—P cont.") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(1, 9.5) + xlim(30, 175) +
  
  annotate("text",
           x     = 30 + 0.6 * (175 - (30)),
           y     = 1 + 0.95 * (9.5 - (1)),
           label = "PF\nEvo\nr = 0.0888",
           hjust = 0,
           vjust = 1,
           size = 3,
           fontface = "bold"
  ) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

PcS.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # Salt v anc: - 1.2753 (< 0.0001), BS: -0.5206 (0.1250)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.002

# R — P content v temperature ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral", "evolved")

df.filt <- df %>% mutate(z.y = T.br, z.x = mean.P.µg.l)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # temp v nit cont: r -0.02940671, p 0.7227

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)      
coef(sma)         

PcT.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "evolved" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "evolved" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  
  annotate("text",
           x     = 30 + 0.6 * (175 - (30)),
           y     = 14 + 0.95 * (22 - (14)),
           label = "PF\nr = -0.0294",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +

  labs(x = xlab.Pc, y = xlab.T, title = "R) Temperature—P cont.") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(14, 22) + xlim(30, 175) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

PcT.qp

# Evolutionary optimization test (evo opt not applicable: no temperature evolution treatment)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.016

# Pigmentation comparisons ------------------------------------------------

# C — Pigmentation v light ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "L", "light", "other"))

df.filt <- df %>% mutate(z.y = I.comp, z.x = pig.PC)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # light v pig: r -0.2039153, p 0.01292

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)      
coef(sma)         

PigL.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "light" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "light" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_abline(intercept = coef(sma)[1], slope = coef(sma)[2], lwd = 0.6, linetype = "dashed") +
  
  annotate("text",
           x     = -3.6 + 0.6 * (3.3 - (-3.6)),
           y     = 0 + 0.95 * (0.6 - (0)),
           label = "PF\nr = -0.204*",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +

  labs(x = xlab.Pig, y = xlab.I, title = "C) Light—Pig. PC1") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) + xlim(-3.6, 3.3) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

PigL.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # light v anc:  + 0.04137 (light further from PF!)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.007

# G — Pigmentation v nitrogen ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "N", "nit", "other"))

df.filt <- df %>% mutate(z.y = N.comp, z.x = pig.PC)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # nit v pig: r -0.04987084, p 0.5472

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)    
coef(sma)         

PigN.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dashed", color = "grey75") +
  
  annotate("text",
           x     = -3.6 + 0.6 * (3.3 - (-3.6)),
           y     = 0 + 0.95 * (1.35 - (0)),
           label = "PF\nr = -0.0499",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +

  labs(x = xlab.Pig, y = xlab.N, title = "G) Nitrogen—Pig. PC1") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 1.35) + xlim(-3.6, 3.3) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

PigN.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # nit v anc: + 0.44907 (nit further from PF!)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.005

# K — Pigmentation v phosphorus ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "P", "phos", "other"))

df.filt <- df %>% mutate(z.y = P.comp, z.x = pig.PC)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # phos v pig: r -0.08083937, p 0.3287

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)      # slope 0.0005833191,
coef(sma)         

PigP.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  labs(x = xlab.Pig, y = xlab.P, title = "K) Phosphorus—Pig. PC1") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 6) + xlim(-3.6, 3.3) +
  
  annotate("text",
           x     = -3.6 + 0.6 * (3.3 - (-3.6)),
           y     = 0 + 0.95 * (6 - 0),
           label = "r = -0.0808",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

PigP.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # + 1.1253 (phos further from PF!)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.807

# O — Pigmentation v salt ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol %in% c("S", "BS"), "salt", "other"))  # Salt stress and Biotic depletion x Salt treated as equivalent

df.filt <- df %>% mutate(z.y = S.c, z.x = pig.PC)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # salt v pig: r 0.02136372, p 0.7966

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)    
coef(sma)         

PigS.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, Evol == "S" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "S" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "BS" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "BS" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  
  annotate("text",
           x     = -3.6 + 0.6 * (3.3 - (-3.6)),
           y     = 1 + 0.95 * (9.5 - (1)),
           label = "PF\nr = 0.0214",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +

  labs(x = xlab.Pig, y = xlab.S, title = "O) Salt—Pig. PC1") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(1, 9.5) + xlim(-3.6, 3.3) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

PigS.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # Salt v anc: - 0.3133 ( 0.5834), BS: 0.2329 (closer to pf)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.004

# S — Pigmentation v temperature ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral", "evolved")

df.filt <- df %>% mutate(z.y = T.br, z.x = pig.PC)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # nit v nit cont: r 0.2564488, p 0.001654

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)  
coef(sma)         

PigT.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "evolved" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "evolved" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_abline(intercept = coef(sma)[1], slope = coef(sma)[2], lwd = 0.6, linetype = "dashed") +

  labs(x = xlab.Pig, y = xlab.T, title = "S) Temperature—Pig. PC1") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(14, 22) + xlim(-3.6, 3.3) +
  
  annotate("text",
           x     = -3.6 + 0.6 * (3.3 - (-3.6)),
           y     = 14 + 0.95 * (22 - 14),
           label = "r = 0.256*",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

PigT.qp

# Evolutionary optimization test (evo opt not applicable: no temperature evolution treatment)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.319

# Biovolume comparisons ---------------------------------------------------

# D — Biovolume v light ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "L", "light", "other"))

df.filt <- df %>% mutate(z.y = I.comp, z.x = bio.vol)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # light v bv: r -0.4955415, p 1.52e-10

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)      
coef(sma)         

BvL.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "light" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "light" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_abline(intercept = coef(sma)[1], slope = coef(sma)[2], lwd = 0.6, linetype = "dashed") +
  
  annotate("text",
           x     = 230 + 0.6 * (545 - (230)),
           y     = 0 + 0.95 * (0.6 - (0)),
           label = "PF\nr = -0.496*",
           hjust = 0,
           vjust = 1,
           size = 3,
           fontface = "bold") +

  labs(x = xlab.Bv, y = xlab.I, title = "D) Light—Biovolume") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) + xlim(230, 545) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

BvL.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # light v anc:  - 0.2097 (0.7316)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.032

# H — Biovolume v nitrogen ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "N", "nit", "other"))

df.filt <- df %>% mutate(z.y = N.comp, z.x = bio.vol)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # nit v bv: r -0.4021031, p 4.067e-07

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)      
coef(sma)         

BvN.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_abline(intercept = coef(sma)[1], slope = coef(sma)[2], lwd = 0.6, linetype = "dashed") +
  
  annotate("text",
           x     = 230 + 0.6 * (545 - (230)),
           y     = 0 + 0.95 * (1.35 - (0)),
           label = "PF\nr = -0.402*",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +

  labs(x = xlab.Bv, y = xlab.N, title = "H) Nitrogen—Biovolume") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 1.35) + xlim(230, 545) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

BvN.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # nit v anc:  + 0.13441 (nit further from PF)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.004

# L — Biovolume v phosphorus ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "P", "phos", "other"))

df.filt <- df %>% mutate(z.y = P.comp, z.x = bio.vol)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # phos v bv: r -0.2104781, p 0.01024

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)  
coef(sma)         

BvP.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_abline(intercept = coef(sma)[1], slope = coef(sma)[2], lwd = 0.6, linetype = "dashed") +
  
  labs(x = xlab.Bv, y = xlab.P, title = "L) Phosphorus—Biovolume") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 6) + xlim(230, 545) +
  
  annotate("text",
           x     = 230 + 0.6 * (545 - 230),
           y     = 0 + 0.95 * (6 - 0),
           label = "r = -0.210*",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

BvP.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # phos v anc: + 1.33505 (phos further from PF)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.068

# P — Biovolume v salt ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol %in% c("S", "BS"), "salt", "other"))  # Salt stress and Biotic depletion x Salt treated as equivalent

df.filt <- df %>% mutate(z.y = S.c, z.x = bio.vol)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # salt v bv: r 0.02663687, p 0.7479

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)  
coef(sma)         

BvS.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, Evol == "S" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "S" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "BS" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol == "BS" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  
  annotate("text",
           x     = 230 + 0.2 * (545 - (230)),
           y     = 1 + 0.95 * (9.5 - (1)),
           label = "PF\nEvo\nr = 0.0266",
           hjust = 0,
           vjust = 1,
           size = 3,
           fontface = "bold"
  ) +

  labs(x = xlab.Bv, y = xlab.S, title = "P) Salt—Biovolume") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(1, 9.5) + xlim(230, 545) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

BvS.qp

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

contrast(em, method = "trt.vs.ctrl", side = "<")  # Salt v anc: - 1.436 (< 0.0001), BS: -0.914 (<0.0001)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.021

# T — Biovolume v temperature ------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral", "evolved")

df.filt <- df %>% mutate(z.y = T.br, z.x = bio.vol)

plot(z.y ~ z.x, data = df.filt)

df.filt <- df.filt %>%
  mutate(z.y2 = scale(z.y)[, 1], z.x2 = scale(z.x)[, 1])

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
         dist.sc  = scale(distance)[, 1]) %>%
  arrange(distance)

df.filt %>% {
  bind_rows(arrange(., desc(z.y)) %>% slice_head(n = 3),
            arrange(., z.y)       %>% slice_head(n = 3),
            arrange(., desc(z.x)) %>% slice_head(n = 3),
            arrange(., z.x)       %>% slice_head(n = 3))
} %>% print()

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")

df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Correlations ######

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # temp v bv: r -0.04593963, p 0.5793

sma <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma)    
coef(sma)         

BvT.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "evolved" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "evolved" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  
  annotate("text",
           x     = 230 + 0.6 * (545 - (230)),
           y     = 14 + 0.95 * (22 - (14)),
           label = "PF\nr = -0.0459",
           size  = 3, hjust = 0, vjust = 1, fontface = "bold") +

  labs(x = xlab.Bv, y = xlab.T, title = "T) Temperature—Biovolume") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(14, 22) + xlim(230, 545) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.03))

BvT.qp

# Evolutionary optimization test (evo opt not applicable: no temperature evolution treatment)

###### Top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

x.max <- max(df.filt$z.x); y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y)

null.df3 <- data.frame(a.emp.n = numeric(), n.PF = numeric(), stringsAsFactors = FALSE)

for (i in 1:1000) {
  
  shuffled.df <- df.filt3 %>%
    mutate(z.x.sim = sample(df.filt3$z.x, replace = FALSE),
           z.y.sim = sample(df.filt3$z.y, replace = FALSE))
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim); y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  null.df3 <- rbind(null.df3, data.frame(a.emp.n = polyarea(poly.n$z.x.sim, poly.n$z.y.sim),
                                         n.PF = nrow(par.res.n)))

}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.002

# Compile & save figure S7 --------------------------------------------------

extra_toffs <- plot_grid(
  NcL.qp,  PcL.qp,  PigL.qp,  BvL.qp,
  NcN.qp,  PcN.qp,  PigN.qp,  BvN.qp,
  NcP.qp,  PcP.qp,  PigP.qp,  BvP.qp,
  NcS.qp,  PcS.qp,  PigS.qp,  BvS.qp,
  NcT.qp,  PcT.qp,  PigT.qp,  BvT.qp,
  ncol = 4, align = "hv", axis = "tblr"
)

ggsave("figures-supplemental/07_fig_s7_extra_inter-gradient_toffs.pdf",
       extra_toffs, width = 11, height = 13.75)