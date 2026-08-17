# Jason R Laurich
# June 2026

# lightly edited for re-submission on August 16th, 2026. 

# We're going to show that the Li et al. (2019) polygonal empty-space test is
# invariant to linear scaling of trait axes (z-scoring), and explore how a
# non-linear transformation (log) can change results.
#

# Structure:
#   Section 1 — Side-by-side plots: raw vs z-scored (same PF, different axes)
#   Section 2 — Empty space + 100 permutations: raw vs z-scored (same p-value)
#   Section 3 — Non-linear transformation: log(P.comp) — does it change things?

# Inputs: 27_summary_table.csv

# Packages & Functions --------------------------------------------------------

library(tidyverse)
library(pracma)   # polyarea
library(cowplot)
library(scam)

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

empty_space_test <- function(x, y, n_perm = 100) {
  # Returns observed area and vector of null areas
  df <- data.frame(x = x, y = y)

  pf <- par_frt(df, xvar = "x", yvar = "y") %>% arrange(x)

  x.max <- max(df$x); y.max <- max(df$y)
  poly <- pf[, c("x", "y")] %>% add_row(x = x.max, y = y.max)
  a.obs <- polyarea(poly$x, poly$y)

  null_areas <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    sh <- df %>% mutate(x = sample(x, replace = FALSE),
                        y = sample(y, replace = FALSE))
    pf_n <- par_frt(sh, xvar = "x", yvar = "y") %>% arrange(x)
    x.max.n <- max(sh$x); y.max.n <- max(sh$y)
    poly_n <- pf_n[, c("x", "y")] %>% add_row(x = x.max.n, y = y.max.n)
    null_areas[i] <- polyarea(poly_n$x, poly_n$y)
  }

  list(a.obs = a.obs, null_areas = null_areas,
       p = mean(null_areas >= a.obs))
}


# Load data -------------------------------------------------------------------

df <- read.csv("processed-data/27_summary_table.csv")

df$Evol.plt <- factor(df$Evol,
                      levels = c("none", "L", "N", "P", "S", "B", "BS", "C"),
                      labels = c("Ancestral", "Light limitation",
                                 "Nitrogen limitation", "Phosphorus limitation",
                                 "Salt stress", "Biotic depletion",
                                 "Biotic depletion x Salt", "Control"))

# Phosphorus panel
df.p <- df %>%
  mutate(raw.x = P.µ.max,
         raw.y = P.comp) %>%
  filter(!is.na(raw.x), !is.na(raw.y)) %>%
  mutate(z.x   = scale(raw.x)[, 1],
         z.y   = scale(raw.y)[, 1],
         log.y = log(raw.y),          # non-linear transform of y (P.comp = 1/P*)
         evol.bin = ifelse(Evol == "none", "ancestral",
                           ifelse(Evol == "P", "phos", "other")))

pf_raw  <- par_frt(df.p, xvar = "raw.x", yvar = "raw.y") %>% arrange(raw.x)
pf_z    <- par_frt(df.p, xvar = "z.x",   yvar = "z.y")   %>% arrange(z.x)
pf_log  <- par_frt(df.p, xvar = "raw.x", yvar = "log.y") %>% arrange(raw.x)

# SCAM curves for display
fit_raw <- scam(raw.y ~ s(raw.x, bs = "mpd", k = min(6, nrow(pf_raw))), data = pf_raw)
fit_z   <- scam(z.y   ~ s(z.x,   bs = "mpd", k = min(6, nrow(pf_z))),   data = pf_z)
fit_log <- scam(log.y ~ s(raw.x, bs = "mpd", k = min(6, nrow(pf_log))), data = pf_log)

curve_raw <- data.frame(raw.x = seq(min(df.p$raw.x), max(df.p$raw.x), length.out = 100)) %>%
  mutate(raw.y = predict(fit_raw, newdata = data.frame(raw.x = raw.x)))
curve_z   <- data.frame(z.x = seq(min(df.p$z.x), max(df.p$z.x), length.out = 100)) %>%
  mutate(z.y = predict(fit_z, newdata = data.frame(z.x = z.x)))
curve_log <- data.frame(raw.x = seq(min(df.p$raw.x), max(df.p$raw.x), length.out = 100)) %>%
  mutate(log.y = predict(fit_log, newdata = data.frame(raw.x = raw.x)))

# PF-optimal flags
df.p <- df.p %>%
  left_join(pf_raw %>% distinct(rep.ID) %>% mutate(pf.raw = "Y"), by = "rep.ID") %>%
  left_join(pf_z   %>% distinct(rep.ID) %>% mutate(pf.z   = "Y"), by = "rep.ID") %>%
  left_join(pf_log %>% distinct(rep.ID) %>% mutate(pf.log = "Y"), by = "rep.ID") %>%
  mutate(across(c(pf.raw, pf.z, pf.log), ~ if_else(is.na(.), "N", .)))

# Section 1 — Side-by-side scatter plots ------------------------------------

pt_size_raw <- ifelse(df.p$pf.raw == "Y", 4, 2)
pt_size_z   <- ifelse(df.p$pf.z   == "Y", 4, 2)
pt_size_log <- ifelse(df.p$pf.log == "Y", 4, 2)

col_map <- c("ancestral" = "black", "phos" = "brown4", "other" = "grey60")

p1.raw <- ggplot(df.p, aes(x = raw.x, y = raw.y)) +
  geom_point(aes(colour = evol.bin,
                 size   = pf.raw,
                 shape  = evol.bin)) +
  geom_line(data = curve_raw, aes(x = raw.x, y = raw.y),
            colour = "black", linewidth = 0.7, inherit.aes = FALSE) +
  scale_colour_manual(values = col_map, guide = "none") +
  scale_shape_manual(values = c("ancestral" = 5, "phos" = 21, "other" = 16),
                     guide = "none") +
  scale_size_manual(values = c("N" = 2, "Y" = 4), guide = "none") +
  labs(title = "A — Raw units",
       x = expression(italic(mu)[max] ~ (day^-1)),
       y = expression(1/italic(P)^"*" ~ (mu*mol^-1))) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.03),
        axis.title = element_text(size = 10))

p1.z <- ggplot(df.p, aes(x = z.x, y = z.y)) +
  geom_point(aes(colour = evol.bin,
                 size   = pf.z,
                 shape  = evol.bin)) +
  geom_line(data = curve_z, aes(x = z.x, y = z.y),
            colour = "black", linewidth = 0.7, inherit.aes = FALSE) +
  scale_colour_manual(values = col_map, guide = "none") +
  scale_shape_manual(values = c("ancestral" = 5, "phos" = 21, "other" = 16),
                     guide = "none") +
  scale_size_manual(values = c("N" = 2, "Y" = 4), guide = "none") +
  labs(title = "B — Z-scored",
       x = expression(italic(mu)[max] ~ (z-score)),
       y = expression(1/italic(P)^"*" ~ (z-score))) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.03),
        axis.title = element_text(size = 10))

p1.log <- ggplot(df.p, aes(x = raw.x, y = log.y)) +
  geom_point(aes(colour = evol.bin,
                 size   = pf.log,
                 shape  = evol.bin)) +
  geom_line(data = curve_log, aes(x = raw.x, y = log.y),
            colour = "black", linewidth = 0.7, inherit.aes = FALSE) +
  scale_colour_manual(values = col_map, guide = "none") +
  scale_shape_manual(values = c("ancestral" = 5, "phos" = 21, "other" = 16),
                     guide = "none") +
  scale_size_manual(values = c("N" = 2, "Y" = 4), guide = "none") +
  labs(title = "C — Log(1/P*), raw µmax",
       x = expression(italic(mu)[max] ~ (day^-1)),
       y = expression(log(1/italic(P)^"*"))) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.03),
        axis.title = element_text(size = 10))

fig1 <- plot_grid(p1.raw, p1.z, p1.log, ncol = 3)
fig1

# Section 2 — Empty space test: raw vs z-scored ----------------------------

set.seed(42)
res_raw <- empty_space_test(df.p$raw.x, df.p$raw.y, n_perm = 100)
res_z   <- empty_space_test(df.p$z.x,   df.p$z.y,   n_perm = 100)

# Section 2 plot — null distributions with observed area marked ------------

null_df <- bind_rows(
  data.frame(area = res_raw$null_areas, type = "Raw"),
  data.frame(area = res_z$null_areas,   type = "Z-scored")
)

obs_df <- data.frame(
  type = c("Raw", "Z-scored"),
  area = c(res_raw$a.obs, res_z$a.obs),
  p    = c(res_raw$p,     res_z$p)
)

# Scale raw areas to z-scored range for visual comparison
scale_factor <- res_z$a.obs / res_raw$a.obs

null_df_scaled <- bind_rows(
  data.frame(area = res_raw$null_areas * scale_factor,
             area_orig = res_raw$null_areas,
             type = "Raw (scaled for comparison)"),
  data.frame(area = res_z$null_areas,
             area_orig = res_z$null_areas,
             type = "Z-scored")
)

obs_df_scaled <- data.frame(
  type = c("Raw (scaled for comparison)", "Z-scored"),
  area = c(res_raw$a.obs * scale_factor, res_z$a.obs),
  p    = c(res_raw$p,                    res_z$p)
)

p2 <- ggplot(null_df_scaled, aes(x = area, fill = type)) +
  geom_histogram(bins = 20, alpha = 0.5, position = "identity", colour = "white") +
  geom_vline(data = obs_df_scaled,
             aes(xintercept = area, colour = type),
             linewidth = 1, linetype = "dashed") +
  geom_text(data = obs_df_scaled,
            aes(x = area, y = Inf,
                label = paste0("p = ", sprintf("%.3f", p)),
                colour = type),
            vjust = 1.5, hjust = -0.1, size = 3.5, show.legend = FALSE) +
  scale_fill_manual(values = c("Raw (scaled for comparison)" = "steelblue",
                                "Z-scored" = "tomato")) +
  scale_colour_manual(values = c("Raw (scaled for comparison)" = "steelblue4",
                                  "Z-scored" = "tomato4")) +
  labs(title = "Null distributions overlap when raw areas are scaled",
       subtitle = "Raw null areas multiplied by (z-obs / raw-obs) for visual alignment",
       x = "Empty space area (aligned units)",
       y = "Count",
       fill = "Data type", colour = "Data type") +
  theme_classic() +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        legend.position = "bottom")

p2

# Section 3 — Non-linear transformation: does log change the test? ---------

set.seed(42)
res_log <- empty_space_test(df.p$raw.x, df.p$log.y, n_perm = 100)

sw_raw <- shapiro.test(df.p$raw.y)
sw_log <- shapiro.test(df.p$log.y)

# Distribution comparison plot
p3a <- ggplot(df.p, aes(x = raw.y)) +
  geom_histogram(bins = 15, fill = "brown4", alpha = 0.7, colour = "white") +
  labs(title = "Raw P.comp",
       subtitle = sprintf("Shapiro p = %.3f", sw_raw$p.value),
       x = expression(1/italic(P)^"*"), y = "Count") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey40"))

p3b <- ggplot(df.p, aes(x = log.y)) +
  geom_histogram(bins = 15, fill = "brown4", alpha = 0.7, colour = "white") +
  labs(title = "Log(P.comp)",
       subtitle = sprintf("Shapiro p = %.3f", sw_log$p.value),
       x = expression(log(1/italic(P)^"*")), y = "Count") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey40"))

# Summary: p-values across all three approaches
summary_df <- data.frame(
  Transformation = c("Raw", "Z-scored\n(linear)", "Log(y)\n(non-linear)"),
  p_value        = c(res_raw$p, res_z$p, res_log$p),
  obs_area       = c(res_raw$a.obs, res_z$a.obs, res_log$a.obs),
  ratio          = c(res_raw$a.obs / mean(res_raw$null_areas),
                     res_z$a.obs   / mean(res_z$null_areas),
                     res_log$a.obs / mean(res_log$null_areas))
)

print(summary_df)

p3c <- ggplot(summary_df, aes(x = Transformation, y = p_value,
                               fill = Transformation)) +
  geom_col(width = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0.6, y = 0.052, label = "α = 0.05",
           hjust = 0, size = 3, colour = "grey40") +
  scale_fill_manual(values = c("Raw" = "steelblue",
                                "Z-scored\n(linear)" = "tomato",
                                "Log(y)\n(non-linear)" = "goldenrod3"),
                    guide = "none") +
  labs(title = "p-value by transformation",
       subtitle = "Linear transforms → identical p; non-linear may differ",
       x = NULL, y = "p-value (Li et al. test)") +
  theme_classic() +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey40"))

p3d <- ggplot(summary_df, aes(x = Transformation, y = ratio,
                               fill = Transformation)) +
  geom_col(width = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("Raw" = "steelblue",
                                "Z-scored\n(linear)" = "tomato",
                                "Log(y)\n(non-linear)" = "goldenrod3"),
                    guide = "none") +
  labs(title = "Obs/mean(null) ratio by transformation",
       subtitle = "Linear transforms → identical ratio",
       x = NULL, y = "Observed area / mean(null)") +
  theme_classic() +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey40"))

fig3 <- plot_grid(p3a, p3b, p3c, p3d, ncol = 2)
fig3
