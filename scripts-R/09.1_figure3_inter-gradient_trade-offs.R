# Jason R Laurich

# June 2026

# This script assesses Pareto fronts across abiotic gradients (inter-gradient niche trait pairs),
# implements Li et al. (2019) polygonal empty space analysis, and evaluates whether experimental
# evolution shifted populations towards Pareto-optimal trait values using a distance-to-PF
# mixed model approach (lmer + emmeans).

# It generates Figure 3 and associated supplementals.

# Inputs: 27_summary_table.csv
# Outputs: figures-main/03_fig_3_inter-gradient_toffs.jpeg
#          figures-supplemental/05.1_fig_s5_inter-gradient_qrs.jpeg
#          figures-supplemental/06.1_fig_s6_inter-gradient_qrs_33.jpeg


# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(quantreg)
library(sp)
library(scam)
library(pracma)
library(cowplot)
library(lme4)
library(emmeans)

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

xlab.I  <- expression(atop("Low-L tolerance",
                            1/italic(L)^"*" ~ (mu*mol^-1 ~ m^2 ~ s)))
xlab.N  <- expression(atop("Low-N tolerance",
                            1/italic(N)^"*" ~ (mu*mol^-1)))
xlab.P  <- expression(atop("Low-P tolerance",
                            1/italic(P)^"*" ~ (mu*mol^-1)))
xlab.S  <- expression(atop("Salt tolerance",
                            italic(S) ~ "(g/L)"))
xlab.T  <- expression(atop("Thermal breadth",
                            italic(T)[italic(br)] ~ "(°C)"))

# A — Light v Nitrogen ----------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "L", "light",
                             ifelse(df$Evol == "N", "nit", "other")))

df.filt <- df %>% mutate(z.y = I.comp, z.x = N.comp)

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

###### Quantile regression ######

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000) # 0.07297, p 0.00696
summary(q75, se = "boot", R = 1000) # 0.09192, p 0.19607
summary(q90, se = "boot", R = 1000) # 0.30813, p 0.01439

LN.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "light"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.N, y = xlab.I, title = "A — Light & Nit") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LN.qr

LN.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
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
  
  geom_point(data = subset(df.filt, evol.bin == "nit" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dashed") +
  
  labs(x = xlab.N, y = xlab.I, title = "A — Light & Nit") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) + xlim(0, 1.35) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LN.qp

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

pairs(em, adjust = "none") # light v anc: - 0.00961 (0.9732), nit v anc: 0.72110 (0.0127)

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)

summary(q50, se = "boot", R = 1000) # -0.05828, p 0.06328
summary(q75, se = "boot", R = 1000) # -0.02992, p 0.64586
summary(q90, se = "boot", R = 1000) # 0.24567, p 0.11803

LN.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt3, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt3, evol.bin == "light"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "nit"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("plum3", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.N, y = xlab.I, title = "A — Light & Nit") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LN.qr2

# B — Light v Phosphorus --------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "L", "light",
                             ifelse(df$Evol == "P", "phos", "other")))

df.filt <- df %>% mutate(z.y = I.comp, z.x = P.comp)

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

###### Quantile regression ######

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000) # 0.01259, P 0.01783
summary(q75, se = "boot", R = 1000) # 0.03210, P 0.00033 
summary(q90, se = "boot", R = 1000) # 0.02415, P 0.29173

LP.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "light"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("brown4", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.P, y = xlab.I, title = "B — Light & Phos") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LP.qr

LP.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
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
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dashed") +
  
  labs(x = xlab.P, y = xlab.I, title = "B — Light & Phos") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) +
  xlim(0, 6) +
  
  annotate(
    "text",
    x = 0 + 0.85 * (6 - 0),
    y = 0 + 0.95 * (0.6 - 0),
    label = "PF",
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +

  theme_classic() +
  
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LP.qp

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

pairs(em, adjust = "none") # light v anc: - 0.5147 (0.0039) (closer to PF!), phos : 0.7493 (<0.0001)

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)

summary(q50, se = "boot", R = 1000) # -0.01226, P 0.02652
summary(q75, se = "boot", R = 1000) # -0.01812, P 0.06453
summary(q90, se = "boot", R = 1000) # -0.04438, P 0.02503

LP.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt3, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt3, evol.bin == "light"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "phos"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("brown4", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.P, y = xlab.I, title = "B — Light & Phos") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LP.qr2

# C — Light v Salt --------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "L", "light",
                             ifelse(df$Evol %in% c("S", "BS"), "salt", "other")))

df.filt <- df %>% mutate(z.y = I.comp, z.x = S.c)

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

###### Quantile regression ######

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000) # 0.00417, P 0.15350
summary(q75, se = "boot", R = 1000) # 0.00260, P 0.38833
summary(q90, se = "boot", R = 1000) # -0.00028, P 0.94739

LS.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "light"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & evol.bin == "salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & evol.bin == "salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("skyblue", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.S, y = xlab.I, title = "C — Light & Salt") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LS.qr

LS.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
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
  
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "grey75", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dashed", colour = "grey75") +
  
  labs(x = xlab.S, y = xlab.I, title = "C — Light & Salt") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) +
  xlim(1, 9.5) +
  
  annotate(
    "text",
    x = 1 + 0.85 * (9.5 - 1),
    y = 0 + 0.95 * (0.6 - 0),
    label = "Evo",
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +

  theme_classic() +
  
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LS.qp

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

df.filt <- df.filt %>%
  mutate(Evol.plt2 = as.character(Evol.plt),
         Evol.plt2 = ifelse(Evol.plt2 %in% c("Salt stress", "Biotic depletion x Salt"), 
                            "match", 
                            Evol.plt2),
         Evol.plt2 = factor(Evol.plt2))

mod <- lmer(dist.to.pf ~ Evol.plt2 + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt2)

pairs(em, adjust = "none") # light v anc: - 0.10586 (0.5850), further from PF. salt: 1.33669 (<0.0001)

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)

summary(q50, se = "boot", R = 1000) # -0.01456, P 0.00091
summary(q75, se = "boot", R = 1000) # -0.01298, P 0.15688 
summary(q90, se = "boot", R = 1000) # -0.03065, P 0.03027

LS.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt3, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt3, evol.bin == "light"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_point(data = subset(df.filt3, Evol.plt == "Salt stress"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt3, Evol.plt == "Biotic depletion x Salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("skyblue", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.S, y = xlab.I, title = "C — Light & Salt") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LS.qr2

# D — Light v Temperature -------------------------------------------------
# No temperature-evolved populations, so evolutionary optimization test for light only.

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "L", "light", "other"))

df.filt <- df %>% mutate(z.y = I.comp, z.x = T.br)

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

###### Quantile regression ######

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000) # -0.00008, P 0.97406
summary(q75, se = "boot", R = 1000) # -0.00547, P 0.41853
summary(q90, se = "boot", R = 1000) # -0.01590, P 0.18730

LT.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "light"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.T, y = xlab.I, title = "D — Light & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LT.qr

LT.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
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
  
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dashed", colour = "grey75") +
  
  labs(x = xlab.T, y = xlab.I, title = "D — Light & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) +
  
  xlim(13.5, 22) +
  
  annotate(
    "text",
    x = 13.5 + 0.85 * (22 - 14),
    y = 0 + 0.95 * (0.6 - 0),
    label = "PF",
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +

  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LT.qp

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

pairs(em, adjust = "none") # light v anc: 0.2409 (0.2958) - light further than ancestral from the PF

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)

summary(q50, se = "boot", R = 1000) # -0.06374, P 0.00004
summary(q75, se = "boot", R = 1000) # -0.06586, P 0.02418 
summary(q90, se = "boot", R = 1000) # -0.04010, P 0.25878

LT.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt3, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt3, evol.bin == "light"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("goldenrod2", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.T, y = xlab.I, title = "D — Light & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 0.6) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

LT.qr2

# E — Nitrogen v Phosphorus -----------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "N", "nit",
                             ifelse(df$Evol == "P", "phos", "other")))

df.filt <- df %>% mutate(z.y = N.comp, z.x = P.comp)

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

###### Quantile regression ######

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000) # 0.07814, P 0.00058
summary(q75, se = "boot", R = 1000) # 0.09402, P 0.00059
summary(q90, se = "boot", R = 1000) # 0.22667, P 0.00160

NP.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("plum3", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("brown4", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.P, y = xlab.N, title = "E — Nit & Phos") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

NP.qr

NP.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
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
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("brown4", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dashed") +
  
  labs(x = xlab.P, y = xlab.N, title = "E — Nit & Phos") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 1.35) +
  xlim(0, 6) +
  
  annotate(
    "text",
    x = 0 + 0.85 * (6 - 0),
    y = 0 + 0.95 * (1.35 - 0),
    label = "PF
Evo",
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +

  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

NP.qp

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

pairs(em, adjust = "none") # nit v anc: - 0.5394 (0.0473). phos: 1.6805 (<0.0001)

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)

summary(q50, se = "boot", R = 1000) # -0.03569, P 0.31026
summary(q75, se = "boot", R = 1000) # -0.06268, P 0.10884
summary(q90, se = "boot", R = 1000) # -0.17067, P 0.06514

NP.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt3, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt3, evol.bin == "nit"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("plum3", 0.6)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "phos"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("brown4", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.P, y = xlab.N, title = "E — Nit & Phos") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

NP.qr2

# F — Nitrogen v Salt -----------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "N", "nit",
                             ifelse(df$Evol %in% c("S", "BS"), "salt", "other")))

df.filt <- df %>% mutate(z.y = N.comp, z.x = S.c)

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

###### Quantile regression ######

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000) # -0.01062, P 0.01168
summary(q75, se = "boot", R = 1000) # -0.02768, P 0.00149
summary(q90, se = "boot", R = 1000) # -0.06527, P 0

NS.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("plum3", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & evol.bin == "salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & evol.bin == "salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("skyblue", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.S, y = xlab.N, title = "F — Nit & Salt") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

NS.qr

NS.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
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
  
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
  
  labs(x = xlab.S, y = xlab.N, title = "F — Nit & Salt") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 1.35) +
  xlim(1, 9.5) +
  
  annotate(
    "text",
    x = 1 + 0.85 * (9.5 - 1),
    y = 0 + 0.95 * (1.35 - 0),
    label = "PF
Evo",
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +

  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

NS.qp

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

df.filt <- df.filt %>%
  mutate(Evol.plt2 = as.character(Evol.plt),
         Evol.plt2 = ifelse(Evol.plt2 %in% c("Salt stress", "Biotic depletion x Salt"), 
                            "match", 
                            Evol.plt2),
         Evol.plt2 = factor(Evol.plt2))

mod <- lmer(dist.to.pf ~ Evol.plt2 + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt2)

pairs(em, adjust = "none") # nit v anc: - 0.1704 (0.2566), salt: 1.0306 (<0.0001)

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)

summary(q50, se = "boot", R = 1000) # -0.07887, P 0
summary(q75, se = "boot", R = 1000) # -0.09032, P 0.00013
summary(q90, se = "boot", R = 1000) # -0.13925, P 0

NS.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt3, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt3, evol.bin == "nit"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("plum3", 0.6)) +
  
  geom_point(data = subset(df.filt3, Evol.plt == "Salt stress"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt3, Evol.plt == "Biotic depletion x Salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("skyblue", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.S, y = xlab.N, title = "F — Nit & Salt") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

NS.qr2

# G — Nitrogen v Temperature ----------------------------------------------
# No temperature-evolved populations; evolutionary optimization test for nit only.

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "N", "nit", "other"))

df.filt <- df %>% mutate(z.y = N.comp, z.x = T.br)

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
fit <- lm(z.y~z.x, data = par.res.1)

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)
pred.curve.1 <- data.frame(z.x = x.vals,
                            z.y = predict(fit, newdata = data.frame(z.x = x.vals)))

pf.ids <- par.res.1 %>% distinct(rep.ID) %>% mutate(pareto.opt = "Y")
df.filt <- df.filt %>%
  left_join(pf.ids, by = "rep.ID") %>%
  mutate(pareto.opt = if_else(is.na(pareto.opt), "N", pareto.opt))

###### Quantile regression ######

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000) # -0.00950, P 0.54289
summary(q75, se = "boot", R = 1000) # 0.00464, P 0.85966
summary(q90, se = "boot", R = 1000) # -0.00888, P 0.91023

NT.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "nit"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("plum3", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.T, y = xlab.N, title = "G — Nit & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

NT.qr

NT.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
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
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dashed", colour = "grey75") +
  
  labs(x = xlab.T, y = xlab.N, title = "G — Nit & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 1.35) +
  xlim(13.5, 22) +

  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

NT.qp

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

pairs(em, adjust = "none") # nit vs ancestral: -0.44538 (0.1587)

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)

summary(q50, se = "boot", R = 1000) # -0.20319, P 0.00005
summary(q75, se = "boot", R = 1000) # -0.16831, P 0.01284
summary(q90, se = "boot", R = 1000) # -0.18959, P 0.04263

NT.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt3, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt3, evol.bin == "nit"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("plum3", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.T, y = xlab.N, title = "G — Nit & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

NT.qr2

# H — Phosphorus v Salt ---------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "P", "phos",
                             ifelse(df$Evol %in% c("S", "BS"), "salt", "other")))

df.filt <- df %>% mutate(z.y = P.comp, z.x = S.c)

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

###### Quantile regression ######

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000) # 0.00290, P 0.92615
summary(q75, se = "boot", R = 1000) # -0.05791, P 0.33835
summary(q90, se = "boot", R = 1000) # 0.00233, P 0.98304

PS.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("brown4", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & evol.bin == "salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & evol.bin == "salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("skyblue", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.S, y = xlab.P, title = "H — Phos & Salt") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

PS.qr

PS.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
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
  
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dashed", colour = "grey75") +
  
  labs(x = xlab.S, y = xlab.P, title = "H — Phos & Salt") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 6) +
  xlim(1, 9.5) +
  
  annotate(
    "text",
    x = 1 + 0.85 * (9.5 - 1),
    y = 0 + 0.95 * (6 - 0),
    label = "PF
Evo",
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +

  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

PS.qp

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

df.filt <- df.filt %>%
  mutate(Evol.plt2 = as.character(Evol.plt),
         Evol.plt2 = ifelse(Evol.plt2 %in% c("Salt stress", "Biotic depletion x Salt"), 
                            "match", 
                            Evol.plt2),
         Evol.plt2 = factor(Evol.plt2))

mod <- lmer(dist.to.pf ~ Evol.plt2 + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt2)

pairs(em, adjust = "none") # phos v anc: - 0.834 (<0.0001); salt -1.714 (<0.0001)

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)

summary(q50, se = "boot", R = 1000) # -0.27609, P 0.00001
summary(q75, se = "boot", R = 1000) # -0.32403, P 0.00008
summary(q90, se = "boot", R = 1000) # -0.50261, P 0.00095

PS.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt3, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt3, evol.bin == "phos"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("brown4", 0.6)) +
  
  geom_point(data = subset(df.filt3, Evol.plt == "Salt stress"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt3, Evol.plt == "Biotic depletion x Salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("skyblue", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.S, y = xlab.P, title = "H — Phos & Salt") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

PS.qr2

# I — Phosphorus v Temperature --------------------------------------------
# No temperature-evolved populations; evolutionary optimization test for phos only.

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol == "P", "phos", "other"))

df.filt <- df %>% mutate(z.y = P.comp, z.x = T.br)

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

###### Quantile regression ######

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000) # -0.08255, P 0.02211
summary(q75, se = "boot", R = 1000) # -0.11350, P 0.09711
summary(q90, se = "boot", R = 1000) # -0.00321, P 0.99155

PT.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "phos"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("brown4", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.T, y = xlab.P, title = "I — Phos & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

PT.qr

PT.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
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
  
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dashed") +
  
  labs(x = xlab.T, y = xlab.P, title = "I — Phos & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(0, 6) +
  xlim(13.5, 22) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

PT.qp

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

pairs(em, adjust = "none") # phos v ancestral: + 0.4531 (0.0957) - phos further from PF

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)

summary(q50, se = "boot", R = 1000) # -0.76422, P 0
summary(q75, se = "boot", R = 1000) # -0.55033, P 0.02090
summary(q90, se = "boot", R = 1000) # -0.60035, P 0.08213

PT.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt3, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt3, evol.bin == "phos"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("brown4", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.T, y = xlab.P, title = "I — Phos & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

PT.qr2

# J — Salt v Temperature --------------------------------------------------
# No temperature-evolved populations; evolutionary optimization test for salt (combined) only.

df$evol.bin <- ifelse(df$Evol == "none", "ancestral",
                      ifelse(df$Evol %in% c("S", "BS"), "salt", "other"))

df.filt <- df %>% mutate(z.y = S.c, z.x = T.br)

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

###### Quantile regression ######

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000) # 0.07003, P 0.45758
summary(q75, se = "boot", R = 1000) # 0.73870, P 0.02528
summary(q90, se = "boot", R = 1000) # -0.05454, P 0.90342

ST.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & evol.bin == "salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & evol.bin == "salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("skyblue", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.T, y = xlab.S, title = "J — Salt & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(1, 9.5) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

ST.qr

ST.qp <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "N"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "other" & pareto.opt == "Y"),
             shape = 16, size = 4, alpha = 0.5) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "N"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, evol.bin == "ancestral" & pareto.opt == "Y"),
             shape = 5, size = 4, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Salt stress" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & pareto.opt == "N"),
             shape = 21, size = 2, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_point(data = subset(df.filt, Evol.plt == "Biotic depletion x Salt" & pareto.opt == "Y"),
             shape = 21, size = 4, stroke = 1, colour = "black",
             fill = scales::alpha("skyblue", 0.6)) +
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
            color = "black", linewidth = 0.6, inherit.aes = FALSE) +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dashed", colour = "grey75") +
  
  labs(x = xlab.T, y = xlab.S, title = "J — Salt & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  ylim(1, 9.5) +
  xlim(13.5, 22) +
  
  annotate(
    "text",
    x = 13.5 + 0.85 * (22 - 14),
    y = 1 + 0.95 * (9.5 - 1),
    label = "PF
Evo",
    hjust = 0,
    vjust = 1,
    size = 3,
    fontface = "bold"
  ) +

  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

ST.qp

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

df.filt <- df.filt %>%
  mutate(Evol.plt2 = as.character(Evol.plt),
         Evol.plt2 = ifelse(Evol.plt2 %in% c("Salt stress", "Biotic depletion x Salt"), 
                            "match", 
                            Evol.plt2),
         Evol.plt2 = factor(Evol.plt2))

mod <- lmer(dist.to.pf ~ Evol.plt2 + (1|Anc), data = df.filt)

em  <- emmeans(mod, ~ Evol.plt2)

pairs(em, adjust = "none") # salt v anc: - 0.8219 (0.0004)

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%
  select(-distance)

q50 <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3)
q75 <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3)
q90 <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)

summary(q50, se = "boot", R = 1000) # -1.56632, P 0.00395
summary(q75, se = "boot", R = 1000) # -1.40537, P 0.00026
summary(q90, se = "boot", R = 1000) # -1.54548, P 0

ST.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = Evol.plt)) +
  
  geom_point(data = subset(df.filt3, evol.bin == "other"),
             shape = 16, size = 2, alpha = 0.5) +
  
  geom_point(data = subset(df.filt3, evol.bin == "ancestral"),
             shape = 5, size = 2, colour = "black", stroke = 1, alpha = 0.6) +
  
  geom_point(data = subset(df.filt3, Evol.plt == "Salt stress" & evol.bin == "salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("navyblue", 0.6)) +
  
  geom_point(data = subset(df.filt3, Evol.plt == "Biotic depletion x Salt" & evol.bin == "salt"),
             shape = 21, size = 2, stroke = 1, colour = "black", fill = scales::alpha("skyblue", 0.6)) +
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 0.6) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 0.6, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 0.6, linetype = "dotted") +
  
  labs(x = xlab.T, y = xlab.S, title = "J — Salt & Temp") +
  
  scale_color_manual(name = "Evolution environment", values = evol_colours) +
  
  theme_classic() +
  
  theme(legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text  = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.03))

ST.qr2

# Compile figure 3 and other supplementals --------------------------------------------------------

legend_df2 <- data.frame(
  x = c(1, 2, 1, 2, 1, 2, 1, 2),
  y = c(1, 1, 2, 2, 1, 1, 2, 2),
  Group = factor(c("Ancestral (N = 5)", "Other", "Matching", "Matching", "Ancestral (N = 5)", "Other", "Matching", "Matching")),
  Group2 = factor(c("Biotic depletion (N = 5)", "Biotic depletion x Salt (N = 4)", "Control (N = 3)", "Light limitation (N = 5)", "Nitrogen limitation (N = 5)", "Ancestral (N = 5)", "Phosphorus limitation (N = 5)", "Salt stress (N = 5)")),
  LineType = factor(c("50th", "75th", "90th", "50th", "75th", "90th", "50th", "75th"))
)

legend_plot2 <- ggplot() +
  
  geom_point(
    data = legend_df2,
    aes(x = x, y = y, colour = Group2),
    shape = 16,
    size = 2,
    alpha = 0.6
  ) +
  
  geom_point(
    data = data.frame(
      x     = c(1, 1),
      y     = c(1, 1),
      Group = factor(c("Ancestral (N = 5)", "Matching"),
                     levels = c("Ancestral (N = 5)", "Matching"))
    ),
    aes(x = x, y = y, shape = Group),
    colour = "black",
    size = 2,
    stroke = 1,
    alpha = 0.6
  ) +
  
  geom_line(
    data = legend_df2,
    aes(x = x, y = y, linetype = LineType),
    colour = "black",
    size = 0.6,
    show.legend = TRUE
  ) +
  
  scale_shape_manual(
    name = NULL,
    values = c("Ancestral (N = 5)" = 5, "Matching" = 21)
  ) +
  
  scale_color_manual(
    name = "Evolution environment",
    values = c(
      "Biotic depletion (N = 5)"        = "chocolate3",
      "Biotic depletion x Salt (N = 4)" = "skyblue",
      "Control (N = 3)"                 = "olivedrab4",
      "Light limitation (N = 5)"        = "goldenrod2",
      "Nitrogen limitation (N = 5)"     = "plum3",
      "Phosphorus limitation (N = 5)"   = "brown4",
      "Salt stress (N = 5)"             = "navyblue"
    )
  ) +
  
  scale_linetype_manual(
    name = "Quantile regression",
    values = c(
      "50th" = "dotted",
      "75th" = "dashed",
      "90th" = "solid"
    )
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
        shape    = c(5,       21),
        colour   = c("black", "black"),
        fill     = c(NA,      "grey75"),
        size     = c(2,       2),
        stroke   = c(1,       1),
        alpha    = c(0.6,     1),
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
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 10),
    legend.key.size  = unit(0.6, "lines"),
    legend.key.width = unit(1.2, "lines")
  )

legend_only2 <- get_legend(legend_plot2)

full_qrs <- plot_grid(
  LN.qr, LP.qr, LS.qr, LT.qr,
  NULL, NP.qr, NS.qr, NT.qr,
  legend_only2, NULL, PS.qr, PT.qr,
  NULL, NULL, NULL, ST.qr,
  ncol = 4,
  align = "hv",
  axis = "tblr"
)

ggsave("figures-supplemental/06.1_fig_s6_inter-gradient_qrs.jpeg", full_qrs, width = 10, height = 10)

full_qrs2 <- plot_grid(
  LN.qr2, LP.qr2, LS.qr2, LT.qr2,
  NULL, NP.qr2, NS.qr2, NT.qr2,
  legend_only2, NULL, PS.qr2, PT.qr2,
  NULL, NULL, NULL, ST.qr2,
  ncol = 4,
  align = "hv",
  axis = "tblr"
)

ggsave("figures-supplemental/07.1_fig_s7_inter-gradient_qrs_33.jpeg", full_qrs2, width = 10, height = 10)

###### Figure 3 ######

legend_df3 <- data.frame(
  x = c(1, 2, 1, 2, 1, 2, 1, 2),
  y = c(1, 1, 2, 2, 1, 1, 2, 2),
  Group = factor(c("Ancestral (N = 5)", "Other", "Matching", "Matching", "Ancestral (N = 5)", "Other", "Matching", "Matching")),
  Group2 = factor(c("Biotic depletion (N = 5)", "Biotic depletion x Salt (N = 4)", "Control (N = 3)", "Light limitation (N = 5)", "Nitrogen limitation (N = 5)", "Ancestral (N = 5)", "Phosphorus limitation (N = 5)", "Salt stress (N = 5)")),
  LineType = factor(c("Pareto front", "50th quantile regression", "Pareto front", "50th quantile regression", "Pareto front", "50th quantile regression", "Pareto front", "50th quantile regression"))
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
      Group = factor(c("Ancestral (N = 5)", "Matching", "Sub-optimal", "Pareto-optimal"),
                     levels = c("Ancestral (N = 5)", "Matching", "Sub-optimal", "Pareto-optimal"))
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
    values = c("Ancestral (N = 5)" = 5, "Matching" = 21, "Sub-optimal" = 21, "Pareto-optimal" = 21)
  ) +
  
  scale_color_manual(
    name = "Evolution environment",
    values = c(
      "Biotic depletion (N = 5)"        = "chocolate3",
      "Biotic depletion x Salt (N = 4)" = "skyblue",
      "Control (N = 3)"                 = "olivedrab4",
      "Light limitation (N = 5)"        = "goldenrod2",
      "Nitrogen limitation (N = 5)"     = "plum3",
      "Phosphorus limitation (N = 5)"   = "brown4",
      "Salt stress (N = 5)"             = "navyblue"
    )
  ) +
  
  scale_linetype_manual(
    name = "Line",
    values = c(
      "Pareto front" = "solid",
      "50th quantile regression" = "dashed"
    ),
    labels = c(
      "Pareto front" = "Pareto front",
      "50th quantile regression" =
        expression(paste("Quantile regression (", tau, " = 0.5)"))
    )
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

full_qps <- plot_grid(
  LN.qp, LP.qp, LS.qp, LT.qp,
  NULL, NP.qp, NS.qp, NT.qp,
  legend_only3, NULL, PS.qp, PT.qp,
  NULL, NULL, NULL, ST.qp,
  ncol = 4,
  align = "hv",
  axis = "tblr"
)

ggsave("figures-main/03.1_fig_3_inter-gradient_toffs.jpeg", full_qps, width = 10, height = 10)
