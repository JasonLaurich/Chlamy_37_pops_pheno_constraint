# Jason R Laurich

# September 5th, 2025

# Creating Figure 1: representations of model fits, parameter extraction and data comparions for methodological clarity

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(cowplot)

# Load & examine data -----------------------------------------------------

df <- read.csv("data-processed/20_summary_table.csv") # Summary file
str(df)
head(df)

df %>%    
  filter(T.br.0.56 %in% c(min(T.br.0.56), max(T.br.0.56), median(T.br.0.56))) %>% 
  select(Pop.fac, T.br.0.56, T.µ.max, S.c, S.µ.max, P.comp.0.56, P.µ.max) %>% 
  print() # Identify populations that vary a lot in T.br

# Ok so we are going to work with 3 populations - Pop.fac 4, anc2, and anc5. These represent a good spread

df.tpc <- read.csv("data-processed/05_TPCs.csv") # TPC file - just to match up numbers (R2jags objects) to the Pop.facs

df.tpc %>%    
  filter(Pop.fac %in% c("4", "anc2", "anc5")) %>% 
  select(Pop.fac, Pop.num) %>% 
  print() # Identify populations that vary a lot in T.br

# So number conversions: 4 - 27, anc2 - 33, anc5 - 36 (pop.num, also how they are recorded in the jags objects!)

# Load the R2jags objects : Temperature

for (i in c(27, 33, 36)){      # Temperature R2jags (model fits)
  load(paste0("R2jags-objects/pop_", i, "_lactin.RData"))
  df.jags <- data.frame(lac_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:6),]
  df.jags.plot$temp <- seq(0, 45, 0.05)
  assign(paste0("df.T.jags", i), df.jags.plot)
}

df.µ.t <- read.csv("data-processed/01_µ_estimates_temp.csv") # Growth data across temperatures
head(df.µ.t) 

df.µ.t <- df.µ.t %>% 
  filter(population.number %in% c(27, 33, 36)) %>% 
  print()

# Load the R2jags objects : Phosphorous

for (i in c(27, 33, 36)){      # Phosphorous Monods
  load(paste0("R2jags-objects/pop_", i, "_phosphorous_monod.RData"))
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:3, 2005),]
  df.jags.plot$phos <- seq(0, 50, 0.025)
  assign(paste0("df.P.jags", i), df.jags.plot)
}

df.µ.p <- read.csv("data-processed/08a_µ_estimates_phosphorous.csv") # Growth data across P levels
head(df.µ.p) # population is the factor, population.number is the corresponding # (e.g. 27, 33, 36)

df.µ.p <- df.µ.p %>% 
  filter(population.number %in% c(27, 33, 36)) %>% 
  print()

# Load the R2jags objects : Salt

for (i in c(27, 33, 36)){      # Salt
  load(paste0("R2jags-objects/pop_", i, "_salt_tolerance.RData"))
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:4, 2006),]
  df.jags.plot$salt <- seq(0, 10, 0.005)
  assign(paste0("df.S.jags", i), df.jags.plot)
}

df.µ.s <- read.csv("data-processed/09a_µ_estimates_salt.csv") # Growth data across salt levels
head(df.µ.s) # population is the factor, population.number is the corresponding # (e.g. 27, 33, 36)

df.µ.s <- df.µ.s %>% 
  filter(population.number %in% c(27, 33, 36)) %>% 
  print()

# Make the figure ---------------------------------------------------------

# Temperature panels

p.t1 <- ggplot(df.µ.t[df.µ.t$population.number==33,], aes(x = temperature, y = r.exp, colour = population)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-2, 5) +
  xlim(0, 45) +
  
  geom_line(data = df.T.jags33, aes(x = temp, y= mean), colour = "forestgreen", linewidth = 1) +
  
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
  
  geom_segment(aes(x = -Inf, xend = 32.7, y = 3.604304, yend = 3.604304), linetype = "dashed", colour = "black", size = 1.2) +
  geom_segment(aes(x = 12, xend = 12 + df[df$Pop.fac == "33",]$T.br.0.56, 
                   y = 0.56, yend = 0.56),
               linetype = "dashed", colour = "black", size = 1.2) +
  
  
  annotate("text", x = 28, y = 3.604304 + 0.45, label = "µ max", size = 4.5, fontface = "bold") +
  annotate("text", x = 28, y = 0.56 + 0.45 , label = "Thermal breadth", size = 4.5, fontface = "bold")

p.t1

p.t2 <- ggplot(df.µ.t[df.µ.t$population.number==27,], aes(x = temperature, y = r.exp, colour = population)) +
  ylim(-2, 5) +
  xlim(0, 45) +
  
  geom_line(data = df.T.jags27, aes(x = temp, y= mean), colour = "magenta2", size = 1) +
  geom_line(data = df.T.jags33, aes(x = temp, y= mean), colour = "forestgreen", size = 1) +
  geom_line(data = df.T.jags36, aes(x = temp, y= mean), colour = "darkorange1", size = 1) +
  
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

p.t3 <- ggplot(df[df$Pop.fac %in% c("4", "anc2", "anc5"),], aes(x = T.µ.max, y = T.br.0.56, colour = Pop.fac)) +
  geom_point(size=4.5) +
  scale_colour_manual(values = c("magenta2", "forestgreen", "darkorange")) +
  
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

# Phosphorous panels

p.p1 <- ggplot(df.µ.p[df.µ.p$population.number==33,], aes(x = phos.lvl, y = r.exp, colour = population)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-0.2, 1.6) +
  
  geom_line(data = df.P.jags27, aes(x = phos, y= mean), colour = "forestgreen", size = 1) +
  
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
  
  geom_segment(aes(x = -Inf, xend = 50, y = 1.324788, yend = 1.324788), linetype = "dashed", colour = "black", size = 1.2) +
  geom_segment(aes(x = 0.7421387, xend = 0.7421387, 
                   y = -Inf, yend = 0.56),
               linetype = "dashed", colour = "black", size = 1.2) +
  
  
  annotate("text", x = 27, y = 1.324788 + 0.15, label = "µ max", size = 4.5, fontface = "bold") +
  annotate("text", x = 3.7, y = 0.28, label = "P*", size = 4.5, fontface = "bold")

p.p1

p.p2 <- ggplot(df.µ.p[df.µ.p$population.number==33,], aes(x = phos.lvl, y = r.exp, colour = population)) +
  scale_colour_manual(values = c("darkorange1", "magenta2", "forestgreen")) +
  ylim(-0.2, 1.6) +
  
  geom_line(data = df.P.jags27, aes(x = phos, y= mean), colour = "magenta2", size = 1) +
  geom_line(data = df.P.jags33, aes(x = phos, y= mean), colour = "forestgreen", size = 1) +
  geom_line(data = df.P.jags36, aes(x = phos, y= mean), colour = "darkorange1", size = 1) +
  
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

p.p3 <- ggplot(df[df$Pop.fac %in% c("4", "anc2", "anc5"),], aes(x = P.µ.max, y = P.comp.0.56, colour = Pop.fac)) +
  geom_point(size=4.5) +
  scale_colour_manual(values = c("magenta2", "forestgreen", "darkorange")) +
  
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

# Salt panels

p.s1 <- ggplot(df.µ.s[df.µ.s$population.number==33,], aes(x = salt.lvl, y = r.exp, colour = population)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) + # small horizontal shift 
  scale_colour_manual(values = c("forestgreen")) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = df.S.jags33, aes(x = salt, y= mean), colour = "forestgreen", size = 1) +
  
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
  
  geom_segment(aes(x = -Inf, xend = 10, y = 1.374038, yend = 1.374038), linetype = "dashed", colour = "black", size = 1.2) +
  geom_segment(aes(x = 4.18795, xend = 4.18795, 
                   y = -Inf, yend = 1.374038/2),
               linetype = "dashed", colour = "black", size = 1.2) +
  
  
  annotate("text", x = 5.3, y = 1.374038 + 0.15, label = "µ max", size = 4.5, fontface = "bold") +
  annotate("text", x = 2, y = 1.374038/4, label = "Salt \n tolerance (c)", size = 4.5, fontface = "bold")

p.s1

p.s2 <- ggplot(df.µ.p[df.µ.p$population.number==3,], aes(x = salt.lvl, y = r.exp, colour = population)) +
  scale_colour_manual(values = c("forestgreen","darkorange1", "magenta2")) +
  ylim(-0.2, 1.7) +
  
  geom_line(data = df.S.jags27, aes(x = salt, y= mean), colour = "magenta2", size = 1) +
  geom_line(data = df.S.jags33, aes(x = salt, y= mean), colour = "forestgreen", size = 1) +
  geom_line(data = df.S.jags36, aes(x = salt, y= mean), colour ="darkorange1", size = 1) +
  
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

p.s3 <- ggplot(df[df$Pop.fac %in% c("4", "anc2", "anc5"),], aes(x = S.µ.max, y = S.c, colour = Pop.fac)) +
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

fig.1 <- plot_grid(p.t1, p.t2, p.t3, p.p1, p.p2, p.p3, p.s1, p.s2, p.s3, nrow = 3, align='hv', rel_widths = c(1,1,1))

fig.1

ggsave("figures/08_fig1_model_fits.jpeg", fig.1, width = 15, height = 15)
