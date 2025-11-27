# Jason R Laurich

# November 27th, 2025

# Creating Figure S1: representations of model fits, parameter extraction and data comparions for methodological clarity

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

df <- read.csv("data-processed/304_summary_table_final.csv") # Summary file
str(df)
head(df)

###### Temperature ######

df %>%
  slice_min(T.br, n = 1) %>%                                        # lowest T.br
  bind_rows(slice_max(df, T.br, n = 1)) %>%                         # highest T.br
  bind_rows(df %>% slice(which.min(abs(T.br - median(T.br))))) %>%  # closest to median
  select(Pop.fac, unique.id, T.br, T.µ.max) %>% 
  print()                                                           # Unique.ids : 27.B03, 18.B07, 7.EO7 

df.t <- read.csv("data-processed/303b_TPC_fits_final.csv") # Get the model fits — actual parameter estimates needed for plotting

df.t <- df.t %>%
  select(Pop.fac, unique.id, Parameter, mean) %>%
  pivot_wider(names_from = Parameter, values_from = mean)

head(df.t)

df.t.mu <- read.csv('data-processed/204_µ_estimates_temp_new.csv') # Growth data, temperature
head(df.t.mu)

df.t.mu <- df.t.mu %>%
  mutate(unique.id = str_c(population.number, str_sub(well.ID, 1, 3), sep = "."))

###### Phosphorous ######

df %>%
  slice_min(P.comp, n = 1) %>%                                          # lowest T.br
  bind_rows(slice_max(df, P.comp, n = 1)) %>%                           # highest T.br
  bind_rows(df %>% slice(which.min(abs(P.comp - median(P.comp))))) %>%  # closest to median
  select(Pop.fac, unique.id, P.comp, P.µ.max) %>% 
  print()                                                               # Unique.ids : 13.F06, 30.D10, 26.D08 

df.p <- read.csv("data-processed/302b_Monod_phosphorous_fits_new.csv") # Get the model fits — actual parameter estimates needed for plotting

df.p <- df.p %>%
  select(Pop.fac, unique.id, Parameter, mean) %>%
  pivot_wider(names_from = Parameter, values_from = mean)

head(df.p)

df.p.mu <- read.csv('data-processed/200_µ_estimates_phosphorous.2-7.csv') # Growth data, phosphorous
head(df.p.mu)

df.p.mu <- df.p.mu %>%
  mutate(unique.id = str_c(population.number, str_sub(well.ID, 1, 3), sep = "."))

###### Salt ######

df %>%
  slice_min(S.c, n = 1) %>%                                          # lowest T.br
  bind_rows(slice_max(df, S.c, n = 1)) %>%                           # highest T.br
  bind_rows(df %>% slice(which.min(abs(S.c - median(S.c))))) %>%     # closest to median
  select(Pop.fac, S.c, S.µ.max) %>%                                  # Working with population-level estimates here
  print()                                                            # Pop.facs : 1, 2, 10 

df.s <- read.csv("data-processed/217_salt_pops_fits_final.csv") # Get the model fits — actual parameter estimates needed for plotting

df.s <- df.s %>%
  select(Pop.fac, Parameter, mean) %>%
  pivot_wider(names_from = Parameter, values_from = mean)

head(df.s)

df.s.mu <- read.csv('data-processed/203_µ_estimates_salt_new.csv') # Growth data, salt
head(df.s.mu)

# Plotting ----------------------------------------------------------------

###### Temperature ######

# Panel A

t.tpc <- df.t %>% filter(unique.id == "7.E07")

curve.med <- tibble::tibble(
  res  = seq(0, 45, length.out = 200),
  rate = pred_lact(res, a = t.tpc$cf.a, b = t.tpc$cf.b, delta_t = t.tpc$cf.delta_t, tmax = t.tpc$cf.tmax)
)

p.t1 <- ggplot(df.t.mu[df.t.mu$unique.id == "7.E07",], aes(x = temperature, y = r.exp)) +
  
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-2, 5) +
  xlim(0, 45) +
  
  geom_line(data = curve.med, aes(x = res, y= rate), colour = "forestgreen", linewidth = 1) +
  
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "A - Temperature") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) +
  
  geom_segment(aes(x = -Inf, xend = 33.467337, y = 4.273670, yend = 4.273670), linetype = "dashed", colour = "black", size = 1.2) +
  geom_segment(aes(x = 21.75, xend = 21.75 + df[df$unique.id == "7.E07",]$T.br, 
                   y = 4.273670/2, yend = 4.273670/2),
               linetype = "dashed", colour = "black", size = 1.2) +
  
  
  annotate("text", x = 28, y = 4.273670 + 0.45, label = "µ max", size = 4.5, fontface = "bold") +
  annotate("text", x = 22.5 + df[df$unique.id == "7.E07",]$T.br/2, y = 4.273670/2 + 0.45 , label = "Thermal breadth", size = 4.5, fontface = "bold")

p.t1

# Panel B

t.tpc <- df.t %>% filter(unique.id == "27.B03")

curve.min <- tibble::tibble(
  res  = seq(0, 45, length.out = 200),
  rate = pred_lact(res, a = t.tpc$cf.a, b = t.tpc$cf.b, delta_t = t.tpc$cf.delta_t, tmax = t.tpc$cf.tmax)
)

t.tpc <- df.t %>% filter(unique.id == "18.B07")

curve.max <- tibble::tibble(
  res  = seq(0, 45, length.out = 200),
  rate = pred_lact(res, a = t.tpc$cf.a, b = t.tpc$cf.b, delta_t = t.tpc$cf.delta_t, tmax = t.tpc$cf.tmax)
)

p.t2 <- ggplot(df.t.mu[df.t.mu$unique.id == "7.E07",], aes(x = temperature, y = r.exp, colour = population)) +
  ylim(-2, 5) +
  xlim(0, 45) +
  
  geom_line(data = curve.min, aes(x = res, y= rate), colour = "magenta2", size = 1) +
  geom_line(data = curve.med, aes(x = res, y= rate), colour = "forestgreen", linewidth = 1) +
  geom_line(data = curve.max, aes(x = res, y= rate), colour = "darkorange1", size = 1) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "B - Temperature") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.t2

# Panel C

p.t3 <- ggplot(df[df$unique.id %in% c("27.B03", "7.E07", "18.B07"),], aes(x = T.µ.max, y = T.br, colour = Pop.fac)) +
  geom_point(size=4.5) +
  scale_colour_manual(values = c("forestgreen", "darkorange", "magenta2")) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Thermal breadth (°C)",
       title = "C - Temperature") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.t3

###### Phosphorous ######

# Panel D

p.monod <- df.p %>% filter(unique.id == "30.D10")

curve.med.p <- tibble::tibble(
  res  = seq(0, 50, length.out = 200),
  rate = pred_mon(res, r.max = p.monod$r_max, k.s = p.monod$K_s)
)

p.p1 <- ggplot(df.p.mu[df.p.mu$unique.id == "30.D10",], aes(x = phos.lvl, y = r.exp)) +
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = curve.med.p, aes(x = res, y= rate), colour = "forestgreen", linewidth = 1) +
  
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Phosphorous concentration (µM)", 
       y = "Exponential growth rate (µ)",
       title = "D - Phosphorous") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) +
  
  geom_segment(aes(x = -Inf, xend = 50, y = 1.439323, yend = 1.439323), linetype = "dashed", colour = "black", size = 1.2) +
  geom_segment(aes(x = 0.8061286, xend = 0.8061286, 
                   y = -Inf, yend = 0.56),
               linetype = "dashed", colour = "black", size = 1.2) +
  
  
  annotate("text", x = 27, y = 1.439323 + 0.15, label = "µ max", size = 4.5, fontface = "bold") +
  annotate("text", x = 3.7, y = 0.28, label = "P*", size = 4.5, fontface = "bold")

p.p1

#Panel E

p.monod <- df.p %>% filter(unique.id == "13.F06")

curve.min.p <- tibble::tibble(
  res  = seq(0, 50, length.out = 200),
  rate = pred_mon(res, r.max = p.monod$r_max, k.s = p.monod$K_s)
)

p.monod <- df.p %>% filter(unique.id == "26.D07")

curve.max.p <- tibble::tibble(
  res  = seq(0, 50, length.out = 200),
  rate = pred_mon(res, r.max = p.monod$r_max, k.s = p.monod$K_s)
)

p.p2 <- ggplot(df.p.mu[df.p.mu$unique.id == "30.D10",], aes(x = phos.lvl, y = r.exp, colour = population)) +
  
  ylim(-0.2, 1.7) +
  
  geom_line(data = curve.min.p, aes(x = res, y= rate), colour = "magenta2", size = 1) +
  geom_line(data = curve.med.p, aes(x = res, y= rate), colour = "forestgreen", linewidth = 1) +
  geom_line(data = curve.max.p, aes(x = res, y= rate), colour = "darkorange1", size = 1) +
  
  labs(x = "Phosphorous concentration (µM)", 
       y = "Exponential growth rate (µ)",
       title = "E - Phosphorous") +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.p2

# Panel F

p.p3 <- ggplot(df[df$unique.id %in% c("13.F06", "30.D10", "26.D07"),], aes(x = P.µ.max, y = P.comp, colour = Pop.fac)) +
  geom_point(size=4.5) +
  scale_colour_manual(values = c("magenta2", "darkorange", "forestgreen")) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/P*)",
       title = "F - Phosphorous") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.p3

###### Salt ######

# Panel G

p.monod <- df.p %>% filter(unique.id == "30.D10")

curve.med.p <- tibble::tibble(
  res  = seq(0, 50, length.out = 200),
  rate = pred_mon(res, r.max = p.monod$r_max, k.s = p.monod$K_s)
)

s.tol <- df.s %>% filter(Pop.fac == 10)

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = s.tol$a, b = s.tol$b, c = s.tol$c)
)

p.s1 <- ggplot(df.s.mu[df.s.mu$population==10,], aes(x = salt, y = r.exp)) +
  
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = curve.s, aes(x = salt, y= rate), colour = "forestgreen", size = 1) +
  
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  labs(x = "Salt concentration (g/L)", 
       y = "Exponential growth rate (µ)",
       title = "G - Salt") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) +
  
  geom_segment(aes(x = -Inf, xend = 10, y = 1.245147, yend = 1.245147), linetype = "dashed", colour = "black", size = 1.2) +
  geom_segment(aes(x = 3.223195, xend = 3.223195, 
                   y = -Inf, yend = 1.245147/2),
               linetype = "dashed", colour = "black", size = 1.2) +
  
  
  annotate("text", x = 5, y = 1.245147 + 0.15, label = "µ max", size = 4.5, fontface = "bold") +
  annotate("text", x = 1.6, y = 1.245147/4, label = "Salt \n tolerance (c)", size = 4.5, fontface = "bold")

p.s1

# Panel H

s.tol <- df.s %>% filter(Pop.fac == 1)

curve.min.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = s.tol$a, b = s.tol$b, c = s.tol$c)
)

s.tol <- df.s %>% filter(Pop.fac == 2)

curve.max.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = s.tol$a, b = s.tol$b, c = s.tol$c)
)

p.s2 <- ggplot(df.s.mu[df.s.mu$population==10,], aes(x = salt, y = r.exp, colour = population)) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = curve.min.s, aes(x = salt, y= rate), colour = "magenta2", size = 1) +
  geom_line(data = curve.s, aes(x = salt, y= rate), colour = "forestgreen", size = 1) +
  geom_line(data = curve.max.s, aes(x = salt, y= rate), colour = "darkorange1", size = 1) +
  
  labs(x = "Salt concentration (g/L)", 
       y = "Exponential growth rate (µ)",
       title = "H - Salt") +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.s2

# Panel I

p.s3 <- ggplot(df[df$Pop.fac %in% c("1", "2", "10"),], aes(x = S.µ.max, y = S.c, colour = Pop.fac)) +
  geom_point(size=4.5) +
  scale_colour_manual(values = c("magenta2", "forestgreen", "darkorange")) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Salt tolerance (c)",
       title = "I - Salt") +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  ) 

p.s3

fig.S1 <- plot_grid(p.t1, p.t2, p.t3, p.p1, p.p2, p.p3, p.s1, p.s2, p.s3, nrow = 3, align='hv', rel_widths = c(1,1,1))

fig.S1

ggsave("figures/100_figs1_model_fits.jpeg", fig.S1, width = 15, height = 15)
