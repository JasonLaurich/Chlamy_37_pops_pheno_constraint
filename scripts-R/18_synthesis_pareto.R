# Jason R Laurich

# Edited for re-submission on August 16th, 2026

# Updated from 14_figure4_synthesis_pareto.R
# Key changes from script 14:
#  (1) Removed 75% inner Pareto front (X.scam) sections
#  (2) Li et al. PF significance test only on thresholded 33% data
#  (3) Triangularity and polygonality removed
#  (4) Qrs replaced with Pearson corrs and SMA

# Inputs: 37_Edwards_2016_light_monods.csv, 40_Lewington_2019_light_monods.csv, 43_Levasseur_2025_light_monods.csv, 
#   68_Narwani_2015_summary.csv, 69_Edwards_2015_summary.csv, 
#   46_Lewington_2019_nit_monods.csv, 49_Levasseur_2025_nit_monods.csv, 
#   51_Bestion_2018_phos_monods.csv, 54_Levasseur_2025_phos_monods.csv, 
#   56_Edwards_2016_TPCs.csv, 59_Levasseur_2025_TPCs.csv, 62_Lewington_2019_TPCs.csv, 
#   64_Bestion_2018_TPCs.csv, 66_Thomas_2012_TPCs.csv, 
#   71_light_metadata_sp.csv, 73_nit_metadata_sp.csv, 75_phos_metadata_sp.csv, 77_temp_metadata_sp.csv,
#   27_summary_table.csv

# Outputs: in figures-main : 04_fig4_synthesis_comps.pdf
#   in figures-supp : 08_fig_s8_synthesis_with_Laurich_2026_data.pdf
#   in data-processed : 70_light_metadata.csv, 72_nit_metadata.csv, 74_phos_metadata.csv,
#       # 76_temp_metadata.csv, 78_synthesis_metadata_sp_summary

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(pracma) # For calculating polygon area
library(cowplot)

library(sp) # For the point.in.polygon function
library(scam)
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

find_nearest_index <- function(x, ref_vec) { # Function to find the closest point to the actual r.max value in the pred data frame. 
  which.min(abs(ref_vec - x))
} # For significance testing based on location of points relative to a curve

# Load and compile the data sets -------------------------------------------------------------------

###### Light synthesis data ######

# Edwards 2016

df.l.e <- read.csv("processed-data/37_Edwards_2016_light_monods.csv")
head(df.l.e)

df.l.e$dataset <- "Edwards et al., 2016"

# Lewington-Pearce 2019

df.l.l <- read.csv("processed-data/40_Lewington_2019_light_monods.csv")
head(df.l.l)

df.l.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.l.lv <- read.csv("processed-data/43_Levasseur_2025_light_monods.csv")
head(df.l.lv)

df.l.lv$dataset <- "Levasseur et al., 2025"

# Narwani 2015

df.l.n <- read.csv("processed-data/68_Narwani_2015_summary.csv")
head(df.l.n)

df.l.n$dataset <- "Narwani et al., 2015"

df.l.n <- df.l.n %>% 
  filter(!is.na(Istar))

df.l <- bind_rows(
  
  df.l.e %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.l.l %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.l.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.l.n %>% 
    transmute(Species = Species.name,
              dataset = dataset,
              r.max = umax.light,
              R = Istar,
              comp = 1/Istar)
  
)

head(df.l) # 131 entries!

df.l <- df.l %>% 
  filter(R > 0) # 130 viable entries

plot(df.l$comp ~ df.l$r.max) # Some really wonky (high) estimates of 1/R* that suggest error (filter at > 3.5)

df.l <- df.l %>% 
  filter(comp < 3.5,
         r.max < 2.5,
         !is.na(comp))

write.csv(df.l, "processed-data/70_light_metadata.csv") # Summary

###### Nitrogen synthesis data ######

# Edwards 2015

df.n.e <- read.csv("processed-data/69_Edwards_2015_summary.csv")
head(df.n.e)

df.n.e <- df.n.e %>% 
  filter(!is.na(k_nit_m),
         !is.na(mu_nit)) %>% 
  mutate(R = 0.1*k_nit_m/(mu_nit - 0.1))

df.n.e$dataset <- "Edwards et al., 2015"

# Lewington-Pearce 2019

df.n.l <- read.csv("processed-data/46_Lewington_2019_nit_monods.csv")
head(df.n.l)

df.n.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.n.lv <- read.csv("processed-data/49_Levasseur_2025_nit_monods.csv")
head(df.n.lv)

df.n.lv$dataset <- "Levasseur et al., 2025"

df.n <- bind_rows(
  
  df.n.e %>% 
    transmute(Species = species,
              dataset = dataset,
              r.max = mu_nit,
              R = R,
              comp = 1/R),
  
  df.n.l %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.n.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post)
  
)

head(df.n) # 58 entries!

df.n <- df.n %>% 
  filter(R > 0) # 57 viable entries

plot(df.n$comp ~ df.n$r.max)

df.n <- df.n %>% 
  filter(comp < 60,
         !is.na(r.max),
         !is.na(comp))

write.csv(df.n, "processed-data/72_nit_metadata.csv") # Summary

###### Phosphorus synthesis data ######

# Bestion 2018

df.p.b <- read.csv("processed-data/51_Bestion_2018_phos_monods.csv")
head(df.p.b)

df.p.b$dataset <- "Bestion et al., 2018"

# Edwards 2015

df.p.e <- read.csv("processed-data/69_Edwards_2015_summary.csv")
head(df.p.e)

df.p.e <- df.p.e %>% 
  filter(!is.na(k_p_m),
         !is.na(mu_p)) %>% 
  mutate(R = 0.1*k_p_m/(mu_p - 0.1))

df.p.e$dataset <- "Edwards et al., 2015"

# Levasseur 2025

df.p.lv <- read.csv("processed-data/54_Levasseur_2025_phos_monods.csv")
head(df.p.lv)

df.p.lv$dataset <- "Levasseur et al., 2025"

# Narwani 2015

df.p.n <- read.csv("processed-data/68_Narwani_2015_summary.csv")
head(df.p.n)

df.p.n$dataset <- "Narwani et al., 2015"

df.p <- bind_rows(
  
  df.p.b %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.p.e %>% 
    transmute(Species = species,
              dataset = dataset,
              r.max = mu_p,
              R = R,
              comp = 1/R),
  
  df.p.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.p.n %>% 
    transmute(Species = Species.name,
              dataset = dataset,
              r.max = umax.nitrate,
              R = Nstar,
              comp = 1/Nstar)
  
)

head(df.p) # 156 entries!

df.p <- df.p %>% 
  filter(R > 0)

plot(df.p$comp ~ df.p$r.max)

df.p <- df.p %>% 
  filter(comp < 300,
         r.max < 2.5,
         !is.na(r.max),
         !is.na(comp))

write.csv(df.p, "processed-data/74_phos_metadata.csv") # Summary

###### Temperature synthesis data ######

# Bestion 2018

df.t.b <- read.csv("processed-data/64_Bestion_2018_TPCs.csv")
head(df.t.b)

df.t.b$dataset <- "Bestion et al., 2018"

# Edwards 2016

df.t.e <- read.csv("processed-data/56_Edwards_2016_TPCs.csv")
head(df.t.e)

df.t.e$dataset <- "Edwards et al., 2016"

# Lewington-Pearce 2019

df.t.l <- read.csv("processed-data/62_Lewington_2019_TPCs.csv")
head(df.t.l)

df.t.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.t.lv <- read.csv("processed-data/59_Levasseur_2025_TPCs.csv")
head(df.t.lv)

df.t.lv$dataset <- "Levasseur et al., 2025"

# Thomas 2012

df.t.t <- read.csv("processed-data/66_Thomas_2012_TPCs.csv")
head(df.t.t)

df.t.t$dataset <- "Thomas et al., 2012"

df.t <- bind_rows(
  
  df.t.b %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.e %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.l %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.t %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
)

head(df.t) # 229 entries!

plot(df.t$T.br ~ df.t$r.max)

df.t <- df.t %>% 
  filter(r.max < 3.5,
         T.br < 36,
         !is.na(r.max),
         !is.na(T.br))

write.csv(df.t, "processed-data/76_temp_metadata.csv") # Summary

###### Get sp-level summary data (means) ######

# For comparing across gradients, where we cannot be confident that observations on the x and y axis were obtained from the same study.

df.l <- read.csv("processed-data/71_light_metadata_sp.csv")
head(df.l)

length(unique(df.l$Sp.name)) # 99 species
sort(unique(df.l$Sp.name))

# Nitrogen

df.n <- read.csv("processed-data/73_nit_metadata_sp.csv")
head(df.n)

length(unique(df.n$Sp.name)) #30 species
sort(unique(df.n$Sp.name))

# Phosphorus 

df.p <- read.csv("processed-data/75_phos_metadata_sp.csv")
head(df.p)

length(unique(df.p$Sp.name)) # 57 species
sort(unique(df.p$Sp.name))

# Temperature

df.t <- read.csv("processed-data/77_temp_metadata_sp.csv")
head(df.t)

length(unique(df.t$Sp.name)) # 117 species
sort(unique(df.t$Sp.name))

# Combine the data frames

df.t.1 <- df.t %>%
  group_by(Sp.name) %>%
  filter(Sp.name != "NA") %>% 
  summarise(
    r.max.T = mean(r.max, na.rm = TRUE),
    T.br    = mean(T.br,  na.rm = TRUE),
    .groups = "drop"
  )

df.n.1 <- df.n %>%
  group_by(Sp.name) %>%
  filter(Sp.name != "NA") %>% 
  summarise(
    r.max.N = mean(r.max, na.rm = TRUE),
    comp.N  = mean(comp,  na.rm = TRUE),
    .groups = "drop"
  )

df.p.1 <- df.p %>%
  group_by(Sp.name) %>%
  filter(Sp.name != "NA") %>% 
  summarise(
    r.max.P = mean(r.max, na.rm = TRUE),
    comp.P  = mean(comp,  na.rm = TRUE),
    .groups = "drop"
  )

df.l.1 <- df.l %>%
  group_by(Sp.name) %>%
  filter(Sp.name != "NA") %>% 
  summarise(
    r.max.L = mean(r.max, na.rm = TRUE),
    comp.L  = mean(comp,  na.rm = TRUE),
    .groups = "drop"
  )

df <- df.t.1 %>%
  full_join(df.n.1, by = "Sp.name") %>%
  full_join(df.p.1, by = "Sp.name") %>%
  full_join(df.l.1, by = "Sp.name")

write.csv(df, "processed-data/78_synthesis_metadata_sp_summary.csv") # Summary

###### Add in our data ######

df.cr <- read.csv("processed-data/27_summary_table.csv") # Summary file
head(df.cr)

# Light v Light -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.L,
    z.x = r.max.L
  ) %>%
  filter(!is.na(z.y), !is.na(z.x))

plot(z.y ~ z.x, data = df.filt)

p.l <- ggplot(df.filt, aes(z.x, z.y)) +
  
  geom_point(size = 1.5, alpha = 0.6, colour = "black") +
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Low-L tolerance",
                           1/italic(L)^"*" ~ (mu*mol^-1 ~ m^2 ~ s))),
       title = "A) Light") +  # labels
  
  ylim(-0.10, 3.25) +
  xlim(-0.30, 3) +
  
  annotate(
    "text",
    x = -0.3 + 0.65 * (3 - -0.3),
    y = -0.1 + 0.95 * (3.25 - -0.1),
    label = "PF\nr = 0.107",
    hjust = 0, vjust = 1, size = 3, fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10, face = "plain"),
    axis.text  = element_text(size = 9,  face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  )

p.l

###### Pareto fronts ######

df.filt <- df.filt %>%
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc  = scale(distance)
  ) %>%
  arrange(distance)

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)

pred.curve.1 <- data.frame(
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # Light: r 0.1066047, p 0.2936

sma.L <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.L)      # slope 1.0532809
coef(sma.L)

p.l3 <- p.l + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
                color = "black", size = 0.6, inherit.aes = FALSE) # PF 33% p = 0 (significant)
p.l3

p.l4 <- p.l3 + geom_point(data = df.cr, aes(x = I.µ.max, y = I.comp),
                size = 1.5, shape = 8, inherit.aes = FALSE)
p.l4

###### Polygonal empty space analysis (Li et al 2019) — top 33% of data ######

x.max <- max(df.filt$z.x)
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>%
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Reference area from full Pareto front

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%  # keep top 33%
  select(-distance)

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF    = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE)
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim")
  par.res.n <- par.res.n %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)
  
  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF    = nrow(par.res.n)
  ))
}

mean(null.df3$a.emp.n >= a.emp) # 0.002

# Light v Nitrogen -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.L,
    z.x = comp.N
  ) %>%
  filter(!is.na(z.y), !is.na(z.x))

plot(z.y ~ z.x, data = df.filt)

x90 <- quantile(df.filt$z.x, 0.90, na.rm = TRUE)
y90 <- quantile(df.filt$z.y, 0.90, na.rm = TRUE)

df.filt.x <- df.filt %>% 
  filter(!(z.x >= x90 & z.y >= y90))

p.ln <- ggplot(df.filt, aes(z.x, z.y)) +
  
  geom_point(size = 1.5, alpha = 0.6, colour = "black") +
  
  labs(x = expression(atop("Low-N tolerance",
                           1/italic(N)^"*" ~ (mu*mol^-1))),
       y = expression(atop("Low-L tolerance",
                           1/italic(L)^"*" ~ (mu*mol^-1 ~ m^2 ~ s))),
       title = "B) Light—Nitrogen") +  # labels
  
  ylim(-0.10, 2) +
  xlim(-0.30, 52) +
  
  annotate(
    "text",
    x = -0.3 + 0.65 * (52 - -0.3),
    y = -0.1 + 0.95 * (2 - -0.1),
    label = "PF\nr = 0.437",
    hjust = 0, vjust = 1, size = 3, fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10, face = "plain"),
    axis.text  = element_text(size = 9,  face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  )

p.ln

###### Pareto fronts ######

df.filt.x <- df.filt.x %>%
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )

x.ref <- min(df.filt.x$z.x2, na.rm = TRUE)
y.ref <- min(df.filt.x$z.y2, na.rm = TRUE)

df.filt.x <- df.filt.x %>%
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc  = scale(distance)
  ) %>%
  arrange(distance)

par.res.1 <- par_frt(df.filt.x, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)

x.vals <- seq(min(df.filt.x$z.x), max(df.filt.x$z.x), length.out = 100)

pred.curve.1 <- data.frame(
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # Light n vit: r 0.4369932, p 0.05403

sma.LN <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.LN)      # slope 0.03934010
coef(sma.LN)

p.ln3 <- p.ln + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
                color = "black", size = 0.6, inherit.aes = FALSE) # PF 33% p = 0.002 (significant)
p.ln3

p.ln4 <- p.ln3 + geom_point(data = df.cr, aes(x = N.comp, y = I.comp),
                size = 1.5, shape = 8, inherit.aes = FALSE)
p.ln4

###### Polygonal empty space analysis (Li et al 2019) — top 33% of data ######

x.max <- max(df.filt.x$z.x)
y.max <- max(df.filt.x$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>%
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Reference area from full Pareto front

df.filt3 <- df.filt.x %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%  # keep top 33%
  select(-distance)

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF    = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE)
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim")
  par.res.n <- par.res.n %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)
  
  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF    = nrow(par.res.n)
  ))
}

mean(null.df3$a.emp.n >= a.emp) # 0.002

# Light v Phosphorus -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.L,
    z.x = comp.P
  ) %>%
  filter(!is.na(z.y), !is.na(z.x))

plot(z.y ~ z.x, data = df.filt)

x90 <- quantile(df.filt$z.x, 0.90, na.rm = TRUE)
y90 <- quantile(df.filt$z.y, 0.90, na.rm = TRUE)

df.filt.x <- df.filt %>% 
  filter(!(z.x >= x90 & z.y >= y90))

p.lp <- ggplot(df.filt, aes(z.x, z.y)) +
  
  geom_point(size = 1.5, alpha = 0.6, colour = "black") +
  
  labs(x = expression(atop("Low-P tolerance",
                           1/italic(P)^"*" ~ (mu*mol^-1))),
       y = expression(atop("Low-L tolerance",
                           1/italic(L)^"*" ~ (mu*mol^-1 ~ m^2 ~ s))),
       title = "C) Light—Phosphorus") +  # labels
  
  ylim(0, 3) +
  xlim(0, 295) +
  
  annotate(
    "text",
    x = 0 + 0.65 * (295 - 0),
    y = 0 + 0.95 * (3 - 0),
    label = "PF\nr = 0.178",
    hjust = 0, vjust = 1, size = 3, fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10, face = "plain"),
    axis.text  = element_text(size = 9,  face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  )

p.lp

###### Pareto fronts ######

df.filt.x <- df.filt.x %>%
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )

x.ref <- min(df.filt.x$z.x2, na.rm = TRUE)
y.ref <- min(df.filt.x$z.y2, na.rm = TRUE)

df.filt.x <- df.filt.x %>%
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc  = scale(distance)
  ) %>%
  arrange(distance)

par.res.1 <- par_frt(df.filt.x, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)

x.vals <- seq(min(df.filt.x$z.x), max(df.filt.x$z.x), length.out = 100)

pred.curve.1 <- data.frame(
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # Light v phos: r 0.1776382, p 0.286

sma.LP <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.LP)      # slope 0.008792447
coef(sma.LP)

p.lp3 <- p.lp + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
                color = "black", size = 0.6, inherit.aes = FALSE) # PF 33% p = 0.010 (significant)
p.lp3

p.lp4 <- p.lp3 + geom_point(data = df.cr, aes(x = P.comp, y = I.comp),
                size = 1.5, shape = 8, inherit.aes = FALSE)
p.lp4

###### Polygonal empty space analysis (Li et al 2019) — top 33% of data ######

x.max <- max(df.filt.x$z.x)
y.max <- max(df.filt.x$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>%
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Reference area from full Pareto front

df.filt3 <- df.filt.x %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%  # keep top 33%
  select(-distance)

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF    = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE)
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim")
  par.res.n <- par.res.n %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)
  
  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF    = nrow(par.res.n)
  ))
}

mean(null.df3$a.emp.n >= a.emp) # 0.010

# Light v Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.L,
    z.x = T.br
  ) %>%
  filter(!is.na(z.y), !is.na(z.x))

plot(z.y ~ z.x, data = df.filt)

x90 <- quantile(df.filt$z.x, 0.90, na.rm = TRUE)
y90 <- quantile(df.filt$z.y, 0.90, na.rm = TRUE)

df.filt.x <- df.filt %>% 
  filter(!(z.x >= x90 & z.y >= y90))

p.lt <- ggplot(df.filt, aes(z.x, z.y)) +
  
  geom_point(size = 1.5, alpha = 0.6, colour = "black") +
  
  labs(x = expression(atop("Thermal breadth",
                           italic(T)[italic(br)] ~ "(\u00b0C)")),
       y = expression(atop("Low-L tolerance",
                           1/italic(L)^"*" ~ (mu*mol^-1 ~ m^2 ~ s))),
       title = "D) Light—Temperature") +  # labels
  
  ylim(0, 2.1) +
  xlim(5, 32) +
  
  annotate(
    "text",
    x = 5 + 0.65 * (32 - 5),
    y = 0 + 0.95 * (2.1 - 0),
    label = "PF\nr = 0.150",
    hjust = 0, vjust = 1, size = 3, fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10, face = "plain"),
    axis.text  = element_text(size = 9,  face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  )

p.lt

###### Pareto fronts ######

df.filt.x <- df.filt.x %>%
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )

x.ref <- min(df.filt.x$z.x2, na.rm = TRUE)
y.ref <- min(df.filt.x$z.y2, na.rm = TRUE)

df.filt.x <- df.filt.x %>%
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc  = scale(distance)
  ) %>%
  arrange(distance)

par.res.1 <- par_frt(df.filt.x, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)

x.vals <- seq(min(df.filt.x$z.x), max(df.filt.x$z.x), length.out = 100)

pred.curve.1 <- data.frame(
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # Light v temp: r 0.1496693, p 0.3208

sma.LT <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.LT)      # slope 0.07358704
coef(sma.LT)

p.lt3 <- p.lt + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
                color = "black", size = 0.6, inherit.aes = FALSE) # PF 33% p = 0.005 (significant)
p.lt3

p.lt4 <- p.lt3 + geom_point(data = df.cr, aes(x = T.br, y = I.comp),
                size = 1.5, shape = 8, inherit.aes = FALSE)
p.lt4

###### Polygonal empty space analysis (Li et al 2019) — top 33% of data ######

x.max <- max(df.filt.x$z.x)
y.max <- max(df.filt.x$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>%
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Reference area from full Pareto front

df.filt3 <- df.filt.x %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%  # keep top 33%
  select(-distance)

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF    = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE)
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim")
  par.res.n <- par.res.n %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)
  
  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF    = nrow(par.res.n)
  ))
}

mean(null.df3$a.emp.n >= a.emp) # 0.005

# Nitrogen v Nitrogen -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.N,
    z.x = r.max.N
  ) %>%
  filter(!is.na(z.y), !is.na(z.x))

plot(z.y ~ z.x, data = df.filt)

p.n <- ggplot(df.filt, aes(z.x, z.y)) +
  
  geom_point(size = 1.5, alpha = 0.6, colour = "black") +
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Low-N tolerance",
                           1/italic(N)^"*" ~ (mu*mol^-1))),
       title = "E) Nitrogen") +  # labels
  
  
  ylim(-0.1, 55) +
  xlim(0, 2.5) +
  
  annotate(
    "text",
    x = 0 + 0.65 * (2.5 - 0),
    y = -0.1 + 0.95 * (55 - -0.1),
    label = "r = 0.472*",
    hjust = 0, vjust = 1, size = 3, fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10, face = "plain"),
    axis.text  = element_text(size = 9,  face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  )

p.n

###### Pareto fronts ######

df.filt <- df.filt %>%
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc  = scale(distance)
  ) %>%
  arrange(distance)

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)
fit <- lm(z.y ~ z.x, data = par.res.1) # Override with lm (too few PF points for scam)

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)

pred.curve.1 <- data.frame(
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # Nitrogen: r 0.4720871, p 0.00844

sma.N <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.N)      # slope 24.58714
coef(sma.N)

p.n3 <- p.n + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
                color = "grey75", size = 0.6, inherit.aes = FALSE) + # PF 33% p = 0.236 (not significant) 

              geom_abline(intercept = coef(sma.N)[1], slope = coef(sma.N)[2], lwd = 0.6, linetype = "dashed") 


p.n3

p.n4 <- p.n3 + geom_point(data = df.cr, aes(x = N.µ.max, y = N.comp),
                size = 1.5, shape = 8, inherit.aes = FALSE)
p.n4

###### Polygonal empty space analysis (Li et al 2019) — top 33% of data ######

x.max <- max(df.filt$z.x)
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>%
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Reference area from full Pareto front

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%  # keep top 33%
  select(-distance)

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF    = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE)
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim")
  par.res.n <- par.res.n %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)
  
  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF    = nrow(par.res.n)
  ))
}

mean(null.df3$a.emp.n >= a.emp) # 0.054

# Nitrogen v Phosphorus -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.N,
    z.x = comp.P
  ) %>%
  filter(!is.na(z.y), !is.na(z.x))

plot(z.y ~ z.x, data = df.filt)

p.np <- ggplot(df.filt, aes(z.x, z.y)) +
  
  geom_point(size = 1.5, alpha = 0.6, colour = "black") +
  
  labs(x = expression(atop("Low-P tolerance",
                           1/italic(P)^"*" ~ (mu*mol^-1))),
       y = expression(atop("Low-N tolerance",
                           1/italic(N)^"*" ~ (mu*mol^-1))),
       title = "F) Nitrogen—Phosphorus") +  # labels
  
  ylim(-0.1, 60) +
  xlim(-5, 185) +
  
  annotate(
    "text",
    x = -5 + 0.65 * (185 - -5),
    y = -0.1 + 0.95 * (60 - -0.1),
    label = "r = -0.0113",
    hjust = 0, vjust = 1, size = 3, fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10, face = "plain"),
    axis.text  = element_text(size = 9,  face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  )

p.np

###### Pareto fronts ######

df.filt <- df.filt %>%
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc  = scale(distance)
  ) %>%
  arrange(distance)

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

# fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1) # Too few points
fit <- lm(z.y ~ z.x, data = par.res.1) # Override with lm (too few PF points for scam)

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)

pred.curve.1 <- data.frame(
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # N v P: r -0.01127607, p 0.9603

sma.NP <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.NP)      # slope -0.2779173
coef(sma.NP)

p.np3 <- p.np + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
                color = "grey75", size = 0.6, inherit.aes = FALSE) # PF 33% p = 0.288 (not significant)
p.np3

p.np4 <- p.np3 + geom_point(data = df.cr, aes(x = P.comp, y = N.comp),
                size = 1.5, shape = 8, inherit.aes = FALSE)
p.np4

###### Polygonal empty space analysis (Li et al 2019) — top 33% of data ######

x.max <- max(df.filt$z.x)
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>%
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Reference area from full Pareto front

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%  # keep top 33%
  select(-distance)

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF    = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE)
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim")
  par.res.n <- par.res.n %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)
  
  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF    = nrow(par.res.n)
  ))
}

mean(null.df3$a.emp.n >= a.emp) # 0.288

# Nitrogen v Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.N,
    z.x = T.br
  ) %>%
  filter(!is.na(z.y), !is.na(z.x))

plot(z.y ~ z.x, data = df.filt)

p.nt <- ggplot(df.filt, aes(z.x, z.y)) +
  
  geom_point(size = 1.5, alpha = 0.6, colour = "black") +
  
  labs(x = expression(atop("Thermal breadth",
                           italic(T)[italic(br)] ~ "(\u00b0C)")),
       y = expression(atop("Low-N tolerance",
                           1/italic(N)^"*" ~ (mu*mol^-1))),
       title = "G) Nitrogen—Temperature") +  # labels
  
  ylim(-1, 45) +
  xlim(12, 28.5) +
  
  annotate(
    "text",
    x = 12 + 0.65 * (28.5 - 12),
    y = -1 + 0.95 * (45 - -1),
    label = "PF\nr = -0.205",
    hjust = 0, vjust = 1, size = 3, fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10, face = "plain"),
    axis.text  = element_text(size = 9,  face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  )

p.nt

###### Pareto fronts ######

df.filt <- df.filt %>%
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc  = scale(distance)
  ) %>%
  arrange(distance)

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)

pred.curve.1 <- data.frame(
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # N and T: r -0.2048489, p 0.3731

sma.NT <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.NT)      # slope -2.921529
coef(sma.NT)

p.nt3 <- p.nt + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
                color = "black", size = 0.6, inherit.aes = FALSE) # PF 33% p = 0 (significant)
p.nt3

p.nt4 <- p.nt3 + geom_point(data = df.cr, aes(x = T.br, y = N.comp),
                size = 1.5, shape = 8, inherit.aes = FALSE)
p.nt4

###### Polygonal empty space analysis (Li et al 2019) — top 33% of data ######

x.max <- max(df.filt$z.x)
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>%
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Reference area from full Pareto front

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%  # keep top 33%
  select(-distance)

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF    = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE)
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim")
  par.res.n <- par.res.n %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)
  
  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF    = nrow(par.res.n)
  ))
}

mean(null.df3$a.emp.n >= a.emp) # 0

# Phosphorus v Phosphorus -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.P,
    z.x = r.max.P
  ) %>%
  filter(!is.na(z.y), !is.na(z.x))

plot(z.y ~ z.x, data = df.filt)

p.p <- ggplot(df.filt, aes(z.x, z.y)) +
  
  geom_point(size = 1.5, alpha = 0.6, colour = "black") +
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Low-P tolerance",
                           1/italic(P)^"*" ~ (mu*mol^-1))),
       title = "H) Phosphorus") +  # labels
  
  ylim(-0.1, 290) +
  xlim(0, 1.8) +
  
  annotate(
    "text",
    x = 0 + 0.65 * (1.8 - 0),
    y = -0.1 + 0.95 * (290 - -0.1),
    label = "PF\nr = -0.219",
    hjust = 0, vjust = 1, size = 3, fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10, face = "plain"),
    axis.text  = element_text(size = 9,  face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  )

p.p

###### Pareto fronts ######

df.filt <- df.filt %>%
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc  = scale(distance)
  ) %>%
  arrange(distance)

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)

pred.curve.1 <- data.frame(
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # Phos: r -0.2192673, p 0.1013

sma.P <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.P)      # slope -204.524
coef(sma.P)

p.p3 <- p.p + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
                color = "black", size = 0.6, inherit.aes = FALSE) # PF 33% p = 0 (significant)
p.p3

p.p4 <- p.p3 + geom_point(data = df.cr, aes(x = P.µ.max, y = P.comp),
                size = 1.5, shape = 8, inherit.aes = FALSE)
p.p4

###### Polygonal empty space analysis (Li et al 2019) — top 33% of data ######

x.max <- max(df.filt$z.x)
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>%
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Reference area from full Pareto front

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%  # keep top 33%
  select(-distance)

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF    = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE)
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim")
  par.res.n <- par.res.n %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)
  
  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF    = nrow(par.res.n)
  ))
}

mean(null.df3$a.emp.n >= a.emp) # 0

# Phosphorus v Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.P,
    z.x = T.br
  ) %>%
  filter(!is.na(z.y), !is.na(z.x))

plot(z.y ~ z.x, data = df.filt)

x90 <- quantile(df.filt$z.x, 0.90, na.rm = TRUE)
y90 <- quantile(df.filt$z.y, 0.90, na.rm = TRUE)

df.filt.x <- df.filt %>% 
  filter(!(z.x >= x90 & z.y >= y90))

p.pt <- ggplot(df.filt, aes(z.x, z.y)) +
  
  geom_point(size = 1.5, alpha = 0.6, colour = "black") +
  
  labs(x = expression(atop("Thermal breadth",
                           italic(T)[italic(br)] ~ "(\u00b0C)")),
       y = expression(atop("Low-P tolerance",
                           1/italic(P)^"*" ~ (mu*mol^-1))),
       title = "I) Phosphorus—Temperature") +  # labels
  
  ylim(0, 190) +
  xlim(11, 28) +
  
  annotate(
    "text",
    x = 11 + 0.65 * (28 - 11),
    y = 0 + 0.95 * (190 - 0),
    label = "r = 0.354",
    hjust = 0, vjust = 1, size = 3, fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10, face = "plain"),
    axis.text  = element_text(size = 9,  face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  )

p.pt

###### Pareto fronts ######

df.filt.x <- df.filt.x %>%
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )

x.ref <- min(df.filt.x$z.x2, na.rm = TRUE)
y.ref <- min(df.filt.x$z.y2, na.rm = TRUE)

df.filt.x <- df.filt.x %>%
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc  = scale(distance)
  ) %>%
  arrange(distance)

par.res.1 <- par_frt(df.filt.x, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)

x.vals <- seq(min(df.filt.x$z.x), max(df.filt.x$z.x), length.out = 100)

pred.curve.1 <- data.frame(
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # P & T: r 0.3538382, p 0.1062

sma.PT <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.PT)      # slope 12.959682
coef(sma.PT)

p.pt3 <- p.pt + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
                color = "grey75", size = 0.6, inherit.aes = FALSE) # PF 33% p = 0.127 (not significant)
p.pt3

p.pt4 <- p.pt3 + geom_point(data = df.cr, aes(x = T.br, y = P.comp),
                size = 1.5, shape = 8, inherit.aes = FALSE)
p.pt4

###### Polygonal empty space analysis (Li et al 2019) — top 33% of data ######

x.max <- max(df.filt.x$z.x)
y.max <- max(df.filt.x$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>%
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Reference area from full Pareto front

df.filt3 <- df.filt.x %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%  # keep top 33%
  select(-distance)

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF    = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE)
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim")
  par.res.n <- par.res.n %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)
  
  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF    = nrow(par.res.n)
  ))
}

mean(null.df3$a.emp.n >= a.emp) # 0.127

# Temperature v Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = T.br,
    z.x = r.max.T
  ) %>%
  filter(!is.na(z.y), !is.na(z.x))

plot(z.y ~ z.x, data = df.filt)

p.t <- ggplot(df.filt, aes(z.x, z.y)) +
  
  geom_point(size = 1.5, alpha = 0.6, colour = "black") +
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Thermal breadth",
                           italic(T)[italic(br)] ~ "(\u00b0C)")),
       title = "J) Temperature") +  # labels
  
  ylim(3, 33) +
  xlim(0, 4) +
  
  annotate(
    "text",
    x = 0 + 0.65 * (4 - 0),
    y = 4 + 0.95 * (33 - 4),
    label = "PF\nr = 0.223*",
    hjust = 0, vjust = 1, size = 3, fontface = "bold"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 10, face = "plain"),
    axis.text  = element_text(size = 9,  face = "plain"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  )

p.t

###### Pareto fronts ######

df.filt <- df.filt %>%
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )

x.ref <- min(df.filt$z.x2, na.rm = TRUE)
y.ref <- min(df.filt$z.y2, na.rm = TRUE)

df.filt <- df.filt %>%
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc  = scale(distance)
  ) %>%
  arrange(distance)

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y")

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6, nrow(par.res.1))), data = par.res.1)

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100)

pred.curve.1 <- data.frame(
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

cor.test(df.filt$z.y, df.filt$z.x, method = "pearson")   # Light: r 0.222619, p 0.01584

sma.T <- sma(z.y ~ z.x, data = df.filt, method = "SMA") # RMA / SMA line for plotting (both axes error-prone, scale-invariant)
summary(sma.T)      # slope 7.287684
coef(sma.T)

p.t3 <- p.t + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y),
                color = "black", size = 0.6, inherit.aes = FALSE) + 
  
  geom_abline(intercept = coef(sma.T)[1], slope = coef(sma.T)[2], lwd = 0.6, linetype = "dashed") 

p.t3

p.t4 <- p.t3 + geom_point(data = df.cr, aes(x = T.µ.max, y = T.br),
                size = 1.5, shape = 8, inherit.aes = FALSE)
p.t4

###### Polygonal empty space analysis (Li et al 2019)) Top 33% of data ######

x.max <- max(df.filt$z.x)
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>%
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Reference area from full Pareto front

df.filt3 <- df.filt %>%
  arrange(distance) %>%
  slice((floor(0.6667 * n()) + 1):n()) %>%  # keep top 33%
  select(-distance)

null.df3 <- data.frame(
  a.emp.n = numeric(),
  n.PF    = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),
      z.y.sim = sample(df.filt3$z.y, replace = FALSE)
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim")
  par.res.n <- par.res.n %>% arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim)
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>%
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)
  
  null.df3 <- rbind(null.df3, data.frame(
    a.emp.n = a.emp.n,
    n.PF    = nrow(par.res.n)
  ))
}

mean(null.df3$a.emp.n >= a.emp) # 0.022

# Assemble Figures 4 and S8 ------------------------------------------------------

meta_cross_PF <- plot_grid(p.l3, p.ln3, p.lp3, p.lt3,
                           NULL, p.n3, p.np3, p.nt3,
                           NULL, NULL, p.p3, p.pt3,
                           NULL, NULL, NULL, p.t3,
                           ncol = 4, align = "hv", axis = "tblr")

ggsave("figures-main/04_fig_4_synthesis_comps.pdf", meta_cross_PF, width = 11.5, height = 11.5)

# Supplemental figure 8 — with Laurich 2026 data

meta_cross_ourdata <- plot_grid(p.l4, p.ln4, p.lp4, p.lt4,
                                NULL, p.n4, p.np4, p.nt4,
                                NULL, NULL, p.p4, p.pt4,
                                NULL, NULL, NULL, p.t4,
                                ncol = 4, align = "hv", axis = "tblr")

ggsave("figures-supplemental/08_fig_s8_synthesis_with_Laurich_2026_data.pdf", meta_cross_ourdata, width = 11.5, height = 11.5)
