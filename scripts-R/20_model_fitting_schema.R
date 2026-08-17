# Jason R Laurich

# February 12th, 2026

# lightly edited for re-submission on August 16th, 2026.

# Here I am going to going to create figure S1, which illustrates the models we are fitting and how we extract parameters from them.

# Adjusted labels and axis titles.

# Inputs: 27_summary_table.csv, 10_phosphorus_µ_estimates.csv, 14_salt_µ_estimates.csv, 18_temperature_µ_estimates.csv

# Outputs: in figures-supplemental : 01_fig_s1_model_fitting.pdf

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(cowplot)

pred_lact <- function(temp, a, b, delta_t, tmax) {
  exp(a * temp) - exp(a * tmax - ((tmax - temp) / delta_t)) + b
}

pred_mon <- function(res, r.max, k.s) {
  r.max * res / (k.s + res)
}

pred_salt <- function(salt, a, b, c) {
  a / (1 + exp(b * (salt - c)))
}

# Load & examine data -----------------------------------------------------

df <- read.csv('processed-data/27_summary_table.csv') # Summary file
head(df)

# Temperature

df %>%
  slice_min(T.br, n = 1) %>%                                        # lowest T.br
  bind_rows(slice_max(df, T.br, n = 1)) %>%                         # highest T.br
  bind_rows(df %>% slice(which.min(abs(T.br - median(T.br))))) %>%  # closest to median
  select(population, rep.ID, T.br, T.µ.max) %>% 
  print()                                                           # rep.IDs 4.B03, 27.C07, 2.E04

df.t.mu <- read.csv('processed-data/18_temperature_µ_estimates.csv') # Growth data, temperature
head(df.t.mu)

df.t.mu <- df.t.mu %>%
  mutate(rep.id = str_c(population, str_sub(well.ID, 1, 3), sep = ".")) # Need to create a unique replicate ID

# Phosphorus

df %>%
  slice_min(P.comp, n = 1) %>%                                          # lowest T.br
  bind_rows(slice_max(df, P.comp, n = 1)) %>%                           # highest T.br
  bind_rows(df %>% slice(which.min(abs(P.comp - median(P.comp))))) %>%  # closest to median
  select(population, rep.ID, P.comp, P.µ.max) %>% 
  print()                                                               # rep.IDs 22.D03, 35.D07, 32.D07

df.p.mu <- read.csv('processed-data/10_phosphorus_µ_estimates.csv') # Growth data, phosphorus
head(df.p.mu)

df.p.mu <- df.p.mu %>%
  mutate(rep.id = str_c(population, str_sub(well.ID, 1, 3), sep = ".")) # Need to create a unique replicate ID

# Salt

df %>%
  slice_min(S.c, n = 1) %>%                                          # lowest T.br
  bind_rows(slice_max(df, S.c, n = 1)) %>%                           # highest T.br
  bind_rows(df %>% slice(which.min(abs(S.c - median(S.c))))) %>%     # closest to median
  select(population, S.c, S.µ.max) %>%                                  # Working with population-level estimates here
  print()                                                            # populations 1, 15, 10 

df.s.mu <- read.csv('processed-data/14_salt_µ_estimates.csv') # Growth data, salt
head(df.s.mu)

# Plotting ----------------------------------------------------------------

###### Temperature ######

# Panel A

t.tpc <- df %>% filter(rep.ID == "2.E04")

t.opt <- t.tpc$T.opt[1]
t.mu.max <- t.tpc$T.µ.max[1]
t.br.min <- t.tpc$T.br.min[1]
t.br.max <- t.tpc$T.br.max[1]

curve.med <- tibble::tibble(
  res  = seq(0, 45, length.out = 200),
  rate = pred_lact(res, a = t.tpc$T.a, b = t.tpc$T.b, delta_t = t.tpc$T.d.t, tmax = t.tpc$T.tmax)
)

p.t1 <- ggplot(df.t.mu[df.t.mu$rep.id == "2.E04",], aes(x = temp, y = µ)) +
  
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-2, 5) +
  xlim(0, 45) +
  
  geom_line(data = curve.med, aes(x = res, y= rate), colour = "forestgreen", linewidth = 0.6) +
  
  geom_jitter(size = 2,width = 0.5, height = 0) + # small horizontal shift 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Temperature (°C)", 
       y = expression("Exponential growth rate (" * italic("\u03bc") * ")"),
       title = "A) Temperature") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 10, face = "plain"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) +
  
  geom_segment(aes(x = -Inf, xend = t.opt, y = t.mu.max - 0.1, yend = t.mu.max - 0.1), linetype = "dashed", colour = "black", size = 0.6) +
  geom_segment(aes(x = t.br.min, xend = t.br.max, 
                   y = t.mu.max/2, yend = t.mu.max/2),
               linetype = "dashed", colour = "black", size = 0.6) +
  
  
  annotate("text", x = t.opt - 5, y = t.mu.max + 0.15, label = expression(italic("\u03bc")[italic(max)]), size = 3, fontface = "plain") +
  annotate("text", x = t.opt - 3.0, y = t.mu.max/2 - 0.25 , label = "Thermal breadth", size = 3, fontface = "plain")

p.t1

# Panel B

t.tpc <- df %>% filter(rep.ID == "4.B03")

curve.min <- tibble::tibble(
  res  = seq(0, 45, length.out = 200),
  rate = pred_lact(res, a = t.tpc$T.a, b = t.tpc$T.b, delta_t = t.tpc$T.d.t, tmax = t.tpc$T.tmax)
)

t.tpc <- df %>% filter(rep.ID == "27.C07")

curve.max <- tibble::tibble(
  res  = seq(0, 45, length.out = 200),
  rate = pred_lact(res, a = t.tpc$T.a, b = t.tpc$T.b, delta_t = t.tpc$T.d.t, tmax = t.tpc$T.tmax)
)

p.t2 <- ggplot(df.t.mu[df.t.mu$unique.id == "2.E04",], aes(x = temp, y = µ, colour = population)) +
  ylim(-2, 5) +
  xlim(0, 45) +
  
  geom_line(data = curve.min, aes(x = res, y= rate), colour = "magenta2", size = 0.6) +
  geom_line(data = curve.med, aes(x = res, y= rate), colour = "forestgreen", size = 0.6) +
  geom_line(data = curve.max, aes(x = res, y= rate), colour = "darkorange1", size = 0.6) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Temperature (°C)", 
       y = expression("Exponential growth rate (" * italic("\u03bc") * ")"),
       title = "B) Temperature") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 10, face = "plain"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) 

p.t2

# Panel C

p.t3 <- ggplot(df[df$rep.ID %in% c("2.E04", "4.B03", "27.C07"),], aes(x = T.µ.max, y = T.br, colour = population)) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = c("forestgreen", "darkorange", "magenta2")) +
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Thermal breadth",
                           italic("T")[italic(br)] * " (°C)")),
       title = "C) Temperature") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 10, face = "plain"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) 

p.t3

###### Phosphorus ######

# Panel D

p.monod <- df %>% filter(rep.ID == "32.D07")

p.mu.max <- p.monod$P.µ.max[1]
p.comp <- p.monod$P.comp[1]

curve.med.p <- tibble::tibble(
  res  = seq(0, 50, length.out = 200),
  rate = pred_mon(res, r.max = p.monod$P.µ.max, k.s = p.monod$P.K.s)
)

p.p1 <- ggplot(df.p.mu[df.p.mu$rep.id == "32.D07",], aes(x = phos, y = µ)) +
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = curve.med.p, aes(x = res, y= rate), colour = "forestgreen", size = 0.6) +
  
  geom_jitter(size = 2,width = 0.5, height = 0) + # small horizontal shift 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Phosphorus concentration (µM)", 
       y = expression("Exponential growth rate (" * italic("\u03bc") * ")"),
       title = "D) Phosphorus") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 10, face = "plain"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) +
  
  geom_segment(aes(x = -Inf, xend = 50, y = p.mu.max, yend = p.monod$p.mu.max), linetype = "dashed", colour = "black", size = 0.6) +
  geom_segment(aes(x = 1/p.comp, xend = 1/p.comp, 
                   y = 0, yend = 0.56),
               linetype = "dashed", colour = "black", size = 0.6) +
  
  
  annotate("text", x = 25, y = p.monod$P.µ.max + 0.1, label = expression(italic("\u03bc")[italic(max)]), size = 3, fontface = "plain") +
  annotate("text", x = 1/p.monod$P.comp + 2, y = 0.28, label = expression(italic("P*")), size = 3, fontface = "plain")

p.p1

#Panel E

p.monod <- df %>% filter(rep.ID == "22.D03")

curve.min.p <- tibble::tibble(
  res  = seq(0, 50, length.out = 200),
  rate = pred_mon(res, r.max = p.monod$P.µ.max, k.s = p.monod$P.K.s)
)

p.monod <- df %>% filter(rep.ID == "35.D07")

curve.max.p <- tibble::tibble(
  res  = seq(0, 50, length.out = 200),
  rate = pred_mon(res, r.max = p.monod$P.µ.max, k.s = p.monod$P.K.s)
)

p.p2 <- ggplot(df.p.mu[df.p.mu$rep.id == "35.D07",], aes(x = phos, y = µ, colour = population)) +
  
  ylim(-0.2, 1.7) +
  
  geom_line(data = curve.min.p, aes(x = res, y= rate), colour = "magenta2", size = 0.6) +
  geom_line(data = curve.med.p, aes(x = res, y= rate), colour = "forestgreen", size = 0.6) +
  geom_line(data = curve.max.p, aes(x = res, y= rate), colour = "darkorange1", size = 0.6) +
  
  labs(x = "Phosphorus concentration (µM)", 
       y = expression("Exponential growth rate (" * italic("\u03bc") * ")"),
       title = "E) Phosphorus") +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 10, face = "plain"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) 

p.p2

# Panel F

p.p3 <- ggplot(df[df$rep.ID %in% c("22.D03", "35.D07", "32.D07"),], aes(x = P.µ.max, y = P.comp, colour = population)) +
  geom_point(size=2.5) +
  scale_colour_manual(values = c("magenta2", "forestgreen", "darkorange")) +
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Low-P tolerance",
                           1/italic(P)^"*" ~ (mu*mol^-1))),
       title = "F) Phosphorus") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 10, face = "plain"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) 

p.p3

###### Salt ######

# Panel G

s.tol <- df %>% filter(population == 10)

s.mu.max <- s.tol$S.µ.max[1]
s.c <- s.tol$S.c[1]

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = s.tol$S.µ.max, b = s.tol$S.b, c = s.tol$S.c)
)

p.s1 <- ggplot(df.s.mu[df.s.mu$population==10,], aes(x = salt, y = µ)) +
  
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = curve.s, aes(x = salt, y= rate), colour = "forestgreen", size = 0.6) +
  
  geom_jitter(size = 2,width = 0.5, height = 0) + # small horizontal shift 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Salt concentration (g/L)", 
       y = expression("Exponential growth rate (" * italic("\u03bc") * ")"),
       title = "G) Salt") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 10, face = "plain"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) +
  
  geom_segment(aes(x = -Inf, xend = 10, y = s.mu.max, yend = s.mu.max), linetype = "dashed", colour = "black", size = 0.6) +
  geom_segment(aes(x = s.c, xend = s.c, 
                   y = 0, yend = s.mu.max/2),
               linetype = "dashed", colour = "black", size = 0.6) +
  
  
  annotate("text", x = 5, y = s.mu.max + 0.1, label = expression(italic("\u03bc")[italic(max)]), size = 3, fontface = "plain") +
  annotate("text", x = s.c - 1.6, y = s.mu.max/4 + 0.06, 
           label = "Salt", size = 3) +
  annotate("text", x = s.c - 1.6, y = s.mu.max/4 - 0.06, 
           label = expression("tolerance ("*italic(S)*")"), size = 3)

p.s1

# Panel H

s.tol <- df %>% filter(population == 1)

curve.min.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = s.tol$S.µ.max, b = s.tol$S.b, c = s.tol$S.c)
)

s.tol <- df %>% filter(population == 15)

curve.max.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = s.tol$S.µ.max, b = s.tol$S.b, c = s.tol$S.c)
)

p.s2 <- ggplot(df.s.mu[df.s.mu$population==10,], aes(x = salt, y = µ, colour = population)) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = curve.min.s, aes(x = salt, y= rate), colour = "magenta2", size = 0.6) +
  geom_line(data = curve.s, aes(x = salt, y= rate), colour = "forestgreen", size = 0.6) +
  geom_line(data = curve.max.s, aes(x = salt, y= rate), colour = "darkorange1", size = 0.6) +
  
  labs(x = "Salt concentration (g/L)", 
       y = expression("Exponential growth rate (" * italic("\u03bc") * ")"),
       title = "H) Salt") +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 10, face = "plain"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) 

p.s2

# Panel I

p.s3 <- ggplot(df[df$population %in% c("1", "10", "15"),], aes(x = S.µ.max, y = S.c, colour = population)) +
  geom_point(size=2.5) +
  scale_colour_manual(values = c("magenta2", "forestgreen", "darkorange")) +
  
  labs(x = expression(atop("Maximum growth rate",
                           italic("\u03bc")[italic(max)] ~ (day^-1))),
       y = expression(atop("Salt tolerance",
                           italic(S) ~ " (g/L)")),
       title = "I) Salt") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 10, face = "plain"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.03)
  ) 

p.s3

fig.S1 <- plot_grid(p.t1, p.t2, p.t3, p.p1, p.p2, p.p3, p.s1, p.s2, p.s3, nrow = 3, align='h', rel_widths = c(1,1,1))

fig.S1

ggsave("figures-supplemental/01_fig_s1_model_fitting.pdf", fig.S1, width = 10, height = 10)
