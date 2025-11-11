# Jason R Laurich

# November 10th, 2025

# Moving towards the end of this project — submission soon hopefully! I want to make sure our data makes sense in the light
# of previous papers that used the same data (metanalyses and Joey's 2020 paper)!

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(cowplot)
library(nls.multstart)

# Load & examine the data -------------------------------------------------

df.hpd <- read.csv("data-processed/20_summary_table.csv") # Summary file, with the confidence intervals
str(df.hpd)
head(df.hpd)

# (a) Evolutionary changes in R* and salt tolerance -----------------------

# We need to calculate changes in values, and then plot them. Also determine whether or not estimate diverge significantly from ancestral counterparts.

df.sum <- df.hpd %>%
  mutate(
    I.star = 1 / I.comp.0.56,                            # calculate I* for each population
    I.anc = 1 / I.comp.0.56[ match(Anc, Pop.fac) ],      # Extract I* of the ancestor
    d.I.star = I.star - I.anc,                           # Calculate difference for the mean.
    d.I.star.L = (1 / I.comp.0.56.L) - I.anc,            # Calculate CIs for the difference
    d.I.star.U = (1 / I.comp.0.56.U) - I.anc,
    
    N.star = 1 / N.comp.0.56,                            # same for N*
    N.anc = 1 / N.comp.0.56[ match(Anc, Pop.fac) ],      
    d.N.star = N.star - N.anc,                          
    d.N.star.L = (1 / N.comp.0.56.L) - N.anc,            
    d.N.star.U = (1 / N.comp.0.56.U) - N.anc,
    
    P.star = 1 / P.comp.0.56,                            # and P*
    P.anc = 1 / P.comp.0.56[ match(Anc, Pop.fac) ],      
    d.P.star = P.star - P.anc,                          
    d.P.star.L = (1 / P.comp.0.56.L) - P.anc,            
    d.P.star.U = (1 / P.comp.0.56.U) - P.anc,
    
    S.c.anc = S.c[ match(Anc, Pop.fac) ],                # and salt tolerance (c) 
    d.S.c = S.c - S.c.anc,                          
    d.S.c.L = S.c.L - S.c.anc,            
    d.S.c.U = S.c.U - S.c.anc
  ) %>%
  select(Pop.fac, Anc, Evol, I.star, d.I.star, d.I.star.L, d.I.star.U, I.µ.max, 
         N.star, d.N.star, d.N.star.L, d.N.star.U, N.µ.max,
         P.star, d.P.star, d.P.star.L, d.P.star.U, P.µ.max,
         S.c, d.S.c, d.S.c.L, d.S.c.U, S.µ.max)

# Recreate Figure 2 - one panel at a time

df.sum <- df.sum %>%
  mutate(Evol = if_else(Evol == "none", "A", Evol)) %>% 
  mutate(Evol = factor(Evol, levels = c("A","C","L","N","P","B","S","BS")))

# Panel A - change in P* ~ evol condition. 

p.P <- ggplot(df.sum, aes(x = Evol, y = d.P.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * P^"*" ~ "(µM P)")  # or I^"*", N^"*", etc.
  ) +
  
  ggtitle("P, old µ data") +
  
  ylim(-1,2) +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.P

df.p <- read.csv("data-processed/210_Monod_phosphorous_pops_estimates.csv") # Summary file, with the confidence intervals
str(df.p)
head(df.p)

df.p <- df.p %>% 
  left_join(df.hpd %>% 
              select(Anc, Evol, Pop.fac), by = "Pop.fac")

df.p <- df.p %>%
  mutate(Evol = if_else(Evol == "none", "A", Evol)) %>% 
  mutate(Evol = factor(Evol, levels = c("A","C","L","N","P","B","S","BS")))

df.p <- df.p %>%
  mutate(
    P.star = R.mth,                                      # P* for each population
    P.anc = P.star[ match(Anc, Pop.fac) ],               # Extract P* of the ancestor
    d.P.star = P.star - P.anc                         # Calculate difference for the mean.
  ) %>%
  select(Pop.fac, Anc, Evol, P.star, d.P.star, r.max)

# Panel A - change in P* ~ evol condition. 

p.P.new <- ggplot(df.p, aes(x = Evol, y = d.P.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * P^"*" ~ "(µM P)")  # or I^"*", N^"*", etc.
  ) +
  
  ylim(-1,2) +
  
  ggtitle("P, new µ data") +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.P.new

fig.P <- plot_grid(p.P, p.P.new, nrow = 1, align='hv', rel_widths = c(1,1))
fig.P

# OK, nitrogen

p.N <- ggplot(df.sum, aes(x = Evol, y = d.N.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * N^"*" ~ "(µM N)")  # or I^"*", N^"*", etc.
  ) +
  
  ggtitle("N, old µ data") +
  
  ylim(-4,4)+
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.N

df.n <- read.csv("data-processed/212_Monod_nitrogen_pops_estimates.csv") # Summary file, with the confidence intervals
str(df.n)
head(df.n)

df.n <- df.n %>% 
  left_join(df.hpd %>% 
              select(Anc, Evol, Pop.fac), by = "Pop.fac")

df.n <- df.n %>%
  mutate(Evol = if_else(Evol == "none", "A", Evol)) %>% 
  mutate(Evol = factor(Evol, levels = c("A","C","L","N","P","B","S","BS")))

df.n <- df.n %>%
  mutate(
    N.star = R.mth,                                      # N* for each population
    N.anc = N.star[ match(Anc, Pop.fac) ],               # Extract N* of the ancestor
    d.N.star = N.star - N.anc                         # Calculate difference for the mean.
  ) %>%
  select(Pop.fac, Anc, Evol, N.star, d.N.star, r.max)

# Panel B - change in N* ~ evol condition. 

p.N.new <- ggplot(df.n, aes(x = Evol, y = d.N.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * P^"*" ~ "(µM P)")  # or I^"*", N^"*", etc.
  ) +
  
  ylim(-4,4) +
  
  ggtitle("N, new µ data") +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.N.new

fig.N <- plot_grid(p.N, p.N.new, nrow = 1, align='hv', rel_widths = c(1,1))
fig.N

# Now light

# Panel C - change in I* ~ evol condition. 

p.I <- ggplot(df.sum, aes(x = Evol, y = d.I.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * I^"*" ~ "(µmol m-2s-1)")  # or I^"*", N^"*", etc.
  ) +
  
  ggtitle("L, old µ data") +
  
  ylim(-14,6) +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.I

df.l <- read.csv("data-processed/214_Monod_light_pops_estimates.csv") # Summary file
str(df.l)
head(df.l)

df.l <- df.l %>% 
  left_join(df.hpd %>% 
              select(Anc, Evol, Pop.fac), by = "Pop.fac")

df.l <- df.l %>%
  mutate(Evol = if_else(Evol == "none", "A", Evol)) %>% 
  mutate(Evol = factor(Evol, levels = c("A","C","L","N","P","B","S","BS")))

df.l <- df.l %>%
  mutate(
    L.star = R.mth,                                      # L* for each population
    L.anc = L.star[ match(Anc, Pop.fac) ],               # Extract L* of the ancestor
    d.L.star = L.star - L.anc                            # Calculate difference for the mean.
  ) %>%
  select(Pop.fac, Anc, Evol, L.star, d.L.star, r.max)

# Panel C - change in I* ~ evol condition. 

p.L.new <- ggplot(df.l, aes(x = Evol, y = d.L.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * P^"*" ~ "(µM P)")  # or I^"*", N^"*", etc.
  ) +
  
  ggtitle("L, new µ data") +
  
  ylim(-5,3) +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.L.new

fig.L <- plot_grid(p.I, p.L.new, nrow = 1, align='hv', rel_widths = c(1,1))
fig.L

# OK salt

p.S <- ggplot(df.sum, aes(x = Evol, y = d.S.c, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in salt tolerance" ~ "(g/L)")  # or I^"*", N^"*", etc.
  ) +
  
  ylim(-2,4) +
  
  ggtitle("S, old µ data")+
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.S

df.s <- read.csv("data-processed/216_salt_pops_estimates_new.csv") # Summary file
str(df.s)
head(df.s)

df.s$Pop.fac[33] <- "anc2"
df.s$Pop.fac[34] <- "anc3"
df.s$Pop.fac[35] <- "anc4"
df.s$Pop.fac[36] <- "anc5"

df.s <- df.s %>% 
  left_join(df.hpd %>% 
              select(Anc, Evol, Pop.fac), by = "Pop.fac")

df.s <- df.s %>%
  mutate(Evol = if_else(Evol == "none", "A", Evol)) %>% 
  mutate(Evol = factor(Evol, levels = c("A","C","L","N","P","B","S","BS")))

df.s <- df.s %>%
  mutate(
    c.pop = c.mod,                                      # c for each population
    c.anc = c.pop[ match(Anc, Pop.fac) ],               # Extract c of the ancestor
    d.c = c.pop - c.anc                            # Calculate difference for the mean.
  ) %>%
  select(Pop.fac, Anc, Evol, c.pop, d.c, r.max)

# Panel C - change in I* ~ evol condition. 

p.S.new <- ggplot(df.s, aes(x = Evol, y = d.c, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in salt tolerance" ~ "(g/L)")  # or I^"*", N^"*", etc.
  ) +
  
  ylim(-2,4) +
  
  ggtitle("S, new µ data")+
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.S.new

fig.S <- plot_grid(p.S, p.S.new, nrow = 1, align='hv', rel_widths = c(1,1))
fig.S

# Figure 4 ----------------------------------------------------------------

p.P2 <- ggplot(df.sum, aes(x = P.star, y = P.µ.max)) +
  
  geom_point() +
  
  ggtitle("P, old µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.P2 

p.P2.new <- ggplot(df.p, aes(x = P.star, y = r.max)) +
  
  geom_point() +
  
  ggtitle("P, new µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.P2.new  

fig.p2 <- plot_grid(p.P2, p.P2.new, nrow = 1, align='hv', rel_widths = c(1,1))
fig.p2

# N

p.N2 <- ggplot(df.sum, aes(x = N.star, y = N.µ.max)) +
  
  geom_point() +
  
  ggtitle("N, old µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.N2  

p.N2.new <- ggplot(df.n, aes(x = N.star, y = r.max)) +
  
  geom_point() +
  
  ggtitle("N, new µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.N2.new  

fig.n2 <- plot_grid(p.N2, p.N2.new, nrow = 1, align='hv', rel_widths = c(1,1))
fig.n2

# light

p.I2 <- ggplot(df.sum, aes(x = I.star, y = I.µ.max)) +
  
  geom_point() +
  
  ggtitle("L, old µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.I2 

p.I2.new <- ggplot(df.l, aes(x = L.star, y = r.max)) +
  
  geom_point() +
  
  ggtitle("L, new µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.I2.new 

fig.l2 <- plot_grid(p.I2, p.I2.new, nrow = 1, align='hv', rel_widths = c(1,1))
fig.l2

# Figure 5 ----------------------------------------------------------------

p.PI <- ggplot(df.sum, aes(x = d.I.star, y = d.P.star)) +
  
  geom_point() +
  
  ggtitle ("P~I, old µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.PI 

df.new <- df.p %>% 
  left_join(df.l, by = "Pop.fac") %>% 
  left_join(df.n, by = "Pop.fac")

p.PI.new <- ggplot(df.new, aes(x = d.L.star, y = d.P.star)) +
  
  geom_point() +
  
  ggtitle ("P~I, new µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.PI.new

fig.pl2 <- plot_grid(p.PI, p.PI.new, nrow = 1, align='hv', rel_widths = c(1,1))
fig.pl2

p.NI <- ggplot(df.sum, aes(x = d.I.star, y = d.N.star)) +
  
  geom_point() +
  
  ggtitle("N~I, old µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.NI 

p.NI.new <- ggplot(df.new, aes(x = d.L.star, y = d.N.star)) +
  
  geom_point() +
  
  ggtitle("N~I, new µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.NI.new

fig.nl2 <- plot_grid(p.NI, p.NI.new, nrow = 1, align='hv', rel_widths = c(1,1))
fig.nl2

p.PN <- ggplot(df.sum, aes(x = d.N.star, y = d.P.star)) +
  
  geom_point() +
  
  ggtitle("P~N, old µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.PN

p.PN.new <- ggplot(df.new, aes(x = d.N.star, y = d.P.star)) +
  
  geom_point() +
  
  ggtitle("P~N, new µ data") +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.PN.new

fig.pn2 <- plot_grid(p.PN, p.PN.new, nrow = 1, align='hv', rel_widths = c(1,1))
fig.pn2

# let's compare them directly? To Joey's data -----------------------------

### We'll start with growth data

# P

df.joey.p <- read.csv("data-processed/100_P_growth_Joey.csv")  # Joey's estimates of µ
head(df.joey.p)

df.joey.p$phos.lvl <- df.joey.p$phosphate_concentration 
df.joey.p$well.ID <- df.joey.p$well_plate

df.jason.p <- read.csv("data-processed/08a_µ_estimates_phosphorous.csv")  # My estimates of µ
head(df.jason.p)

df.growth <- df.jason.p %>% 
  left_join(df.joey.p %>% 
              select(population, well.ID, phos.lvl, estimate),
            by = c("population", "well.ID", "phos.lvl"))


df.sum$population <- df.sum$Pop.fac

df.growth2 <- df.growth %>% 
  left_join(df.sum %>% 
              select(Evol, population), by = "population")

reg <- ggplot(df.growth2, aes(x= r.exp, y =estimate, colour = phos.lvl)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

reg

df.joey.n <- read.csv("data-processed/101_N_growth_Joey.csv")  # Joey's estimates of µ
head(df.joey.n)

df.joey.n$nitrate.lvl <- df.joey.n$nitrate_concentration 
df.joey.n$well.ID <- df.joey.n$well_plate

df.jason.n <- read.csv("data-processed/07a_µ_estimates_nitrogen_old.csv")  # My estimates of µ
head(df.jason.n)

df.growth.n <- df.jason.n %>% 
  left_join(df.joey.n %>% 
              select(population, well.ID, nitrate.lvl, estimate),
            by = c("population", "well.ID", "nitrate.lvl"))

reg2 <- ggplot(df.growth.n, aes(x= r.exp, y =estimate, colour = nitrate.lvl)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

reg2

nold <- lm(estimate ~ r.exp, data= df.growth.n)
summary(nold)

# OK let's look at R*s and µmaxes now! ------------------------------------

df.tot <- df.sum %>% 
  left_join(df.new, by = "Pop.fac")

head(df.tot)

df.joey.p <- read.csv("data-processed/100_P_growth_Joey.csv")  # Joey's estimates of µ
head(df.joey.p)

df.joey.p$phos.lvl <- df.joey.p$phosphate_concentration 
df.joey.p$well.ID <- df.joey.p$well_plate

df.jason.p <- read.csv("data-processed/08a_µ_estimates_phosphorous.csv")  # My estimates of µ
head(df.jason.p)

df.growth <- df.jason.p %>% 
  left_join(df.joey.p %>% 
              select(population, well.ID, phos.lvl, estimate),
            by = c("population", "well.ID", "phos.lvl"))