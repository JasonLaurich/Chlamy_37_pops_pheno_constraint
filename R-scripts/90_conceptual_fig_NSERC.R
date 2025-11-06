# Jason R Laurich

# October 11th, 2025

# Creating a conceptual figure for my NSERC application.

# Variation in growth with and without microbes (TPC, salt)
# 4 panels? Then the representation of two G matrices?

# Figure 2: Start with 20 TPCs? and then extract at a given temperature or two the relationship between Tbr and r?
# Fitness regression/ strength of selection on standardize traits data (with relative fitness)

# Load packages & specify functions -----------------------------------------------------------

library(tidyverse)
library(cowplot)

pred_lact <- function(temp, a, b, delta_t, tmax) {     
  exp(a * temp) - exp(a * tmax - ((tmax - temp) / delta_t)) + b
} # For fitting a lactin curve after modifying variables

pred_salt <- function(salt, a, b, c) {
  a / (1 + exp(b * (salt - c)))
} # And for salt

# Load & examine the data -------------------------------------------------

df.t <- read_csv('data-processed/05a_TPC_fits.csv')
head(df.t)

df.t <- df.t %>%
  filter(Model == 'Lactin 2') %>% 
  select(Pop.fac, Parameter, mean) %>%
  pivot_wider(names_from = Parameter, values_from = mean)


df.s <- read_csv('data-processed/09c_salt_tolerance_fits.csv')
head(df.s)

df.s <- df.s %>%
  select(Pop.fac, Parameter, mean) %>%
  pivot_wider(names_from = Parameter, values_from = mean) %>% 
  mutate(Pop.fac = case_when(
    Pop.fac == "Anc 2" ~ "anc2",
    Pop.fac == "Anc 3" ~ "anc3",
    Pop.fac == "Anc 4" ~ "anc4",
    Pop.fac == "Anc 5" ~ "anc5",
    TRUE ~ Pop.fac
  ))

head(df.s)

# Let's make panel 1 and 2 : 3 genotypes across temp and salt grad --------

###### Temperature gradients — first no microbes ######

pop <- df.t %>% filter(Pop.fac == "1") # Choose your population
pop1_df <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.t %>% filter(Pop.fac == "2") # Choose your population
pop2_df <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.t %>% filter(Pop.fac == "3") # Choose your population
pop3_df <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.t %>% filter(Pop.fac == "4") # Choose your population
pop4_df <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.t %>% filter(Pop.fac == "5") # Choose your population
pop5_df <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.t %>% filter(Pop.fac == "6") # Choose your population
pop6_df <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

p.t <- ggplot(pop1_df, aes(x = temp, y =rate)) +
  theme_classic() +
  
  geom_line(data = pop1_df, aes(x = temp, y= rate), colour = "black", linewidth = 1) +
  geom_line(data = pop2_df, aes(x = temp, y= rate), colour = "gray", linewidth = 1) +
  geom_line(data = pop3_df, aes(x = temp, y= rate), colour = "yellow", linewidth = 1) +
  geom_line(data = pop4_df, aes(x = temp, y= rate), colour = "red", linewidth = 1) +
  geom_line(data = pop5_df, aes(x = temp, y= rate), colour = "blue", linewidth = 1) +
  geom_line(data = pop6_df, aes(x = temp, y= rate), colour = "purple", linewidth = 1) +
  
  labs(
    x = "Temperature (°C)",
    y = "Growth rate",
    title = "TPC"
  ) +
  
  ylim(-1,5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2)

p.t
  
# We'll go with 3, 4, and 6

p.t <- ggplot(pop1_df, aes(x = temp, y =rate)) +
  theme_classic() +
  
  geom_line(data = pop3_df, aes(x = temp, y= rate), colour = "gold", linewidth = 1) +
  geom_line(data = pop4_df, aes(x = temp, y= rate), colour = "firebrick2", linewidth = 1) +
  geom_line(data = pop6_df, aes(x = temp, y= rate), colour = "blue3", linewidth = 1) +
  
  labs(
    x = "Temperature (°C)",
    y = "Predicted exponential growth rate",
    title = "TPC"
  ) +
  
  ylim(-1,5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2)

p.t

pop <- df.t %>% filter(Pop.fac == "6") # Choose your population
pop <- pop %>% 
  mutate(cf.tmax = 42, delta_t = cf.b *1.1, )

pop6_df.mod <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

ggplot(pop1_df, aes(x = temp, y =rate)) +
  theme_classic() +
  
  geom_line(data = pop6_df.mod, aes(x = temp, y= rate), colour = "black", linewidth = 1) +
  geom_line(data = pop6_df, aes(x = temp, y= rate), colour = "blue3", linewidth = 1) +
  
  labs(
    x = "Temperature (°C)",
    y = "Predicted exponential growth rate",
    title = "TPC"
  ) +
  
  ylim(-1,5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2)

pop <- df.t %>% filter(Pop.fac == "4") # Choose your population
pop <- pop %>% 
  mutate(cf.tmax = 38)

pop4_df.mod <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

ggplot(pop1_df, aes(x = temp, y =rate)) +
  theme_classic() +
  
  geom_line(data = pop4_df.mod, aes(x = temp, y= rate), colour = "black", linewidth = 1) +
  geom_line(data = pop4_df, aes(x = temp, y= rate), colour = "firebrick2", linewidth = 1) +
  
  labs(
    x = "Temperature (°C)",
    y = "Predicted exponential growth rate",
    title = "TPC"
  ) +
  
  ylim(-1,5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2)

pop <- df.t %>% filter(Pop.fac == "3") # Choose your population
pop <- pop %>% 
  mutate(cf.a = cf.a*1.10)

pop3_df.mod <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

p.t.1 <- ggplot(pop1_df, aes(x = temp, y =rate)) +
  theme_classic() +
  
  geom_line(data = pop4_df.mod, aes(x = temp, y= rate), colour = "firebrick2", linewidth = 1.5) +
  geom_line(data = pop6_df.mod, aes(x = temp, y= rate), colour = "blue3", linewidth = 1.5) +
  geom_line(data = pop3_df.mod, aes(x = temp, y= rate), colour = "gold", linewidth = 1.5) +
  
  labs(
    x = "Temperature (°C)",
    y = "Growth rate",
    title = "A"
  ) +
  
  ylim(-1,6) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2)

p.t.1

pop <- df.t %>% filter(Pop.fac == "3") # Choose your population
pop <- pop %>% 
  mutate(cf.a = cf.a*1.12, cf.tmax = 42.5)

pop3_df.mod2 <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.t %>% filter(Pop.fac == "4") # Choose your population
pop <- pop %>% 
  mutate(cf.tmax = 38.5, cf.a = cf.a*1.21)

pop4_df.mod2 <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.t %>% filter(Pop.fac == "6") # Choose your population
pop <- pop %>% 
  mutate(cf.tmax = 42.5, cf.b = cf.b*0.6)

pop6_df.mod2 <-
  seq(5, 45, length.out = 200) %>%
  tibble(temp = .) %>%
  mutate(
    rate = pred_lact(
      temp,
      a       = as.numeric(pop$cf.a),
      b       = as.numeric(pop$cf.b),
      delta_t = as.numeric(pop$cf.delta_t),
      tmax    = as.numeric(pop$cf.tmax)
    ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

p.t.2 <- ggplot(pop1_df, aes(x = temp, y =rate)) +
  theme_classic() +
  
  geom_line(data = pop4_df.mod2, aes(x = temp, y= rate), colour = "firebrick2", linewidth = 1.5, linetype = 'dashed') +
  geom_line(data = pop4_df.mod, aes(x = temp, y= rate), colour = "firebrick2", linewidth = 0.5) +
  geom_line(data = pop6_df.mod2, aes(x = temp, y= rate), colour = "blue3", linewidth = 1.5, linetype = 'dashed') +
  geom_line(data = pop6_df.mod, aes(x = temp, y= rate), colour = "blue3", linewidth = 0.5) +
  geom_line(data = pop3_df.mod2, aes(x = temp, y= rate), colour = "gold", linewidth = 1.5, linetype='dashed') +
  geom_line(data = pop3_df.mod, aes(x = temp, y= rate), colour = "gold", linewidth = 0.5) +
  
  labs(
    x = "Temperature (°C)",
    y = "Growth rate",
    title = "B"
  ) +
  
  ylim(-1,6) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2)

p.t.2

plot_grid(p.t.1, p.t.2, nrow = 2, align='hv')

###### Salt ######

pop <- df.s %>% filter(Pop.fac == "1") # Choose your population
pop <- pop %>% 
  mutate(a=a*1.3)

pop1_df.mod <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
    Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.s %>% filter(Pop.fac == "2") # Choose your population
pop <- pop %>% 
  mutate(b=b*1.2)

pop2_df.mod <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
  Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.s %>% filter(Pop.fac == "3") # Choose your population
pop <- pop %>% 
  mutate(c = c + 2.1, b = b*0.84, a = a + 0.3)

pop3_df.mod <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
  Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

p.s.1 <- ggplot(pop1_df, aes(x = salt, y =rate)) +
  theme_classic() +
  
  geom_line(data = pop1_df.mod, aes(x = salt, y= rate), colour = "firebrick2", linewidth = 1.5) +
  geom_line(data = pop2_df.mod, aes(x = salt, y= rate), colour = "blue3", linewidth = 1.5) +
  geom_line(data = pop3_df.mod, aes(x = salt, y= rate), colour = "gold", linewidth = 1.5) +
  
  labs(
    x = "Salinity (g/L)",
    y = "Growth rate",
    title = ""
  ) +
  
  ylim(-0.5,3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2)

p.s.1

pop <- df.s %>% filter(Pop.fac == "1") # Choose your population
pop <- pop %>% 
  mutate(a=a*1.8, c = c + 0.5, b = b *0.9)

pop1_df.mod2 <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
  Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.s %>% filter(Pop.fac == "2") # Choose your population
pop <- pop %>% 
  mutate(c= c +0.3, a = a*1.1)

pop2_df.mod2 <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
  Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.s %>% filter(Pop.fac == "3") # Choose your population
pop <- pop %>% 
  mutate(c = c + 1, b = b*0.6, a = a + 1.2)

pop3_df.mod2 <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
  Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

p.s.2 <- ggplot(pop1_df, aes(x = salt, y =rate)) +
  theme_classic() +
  
  geom_line(data = pop1_df.mod, aes(x = salt, y= rate), colour = "firebrick2", linewidth = 0.5) +
  geom_line(data = pop1_df.mod2, aes(x = salt, y= rate), colour = "firebrick2", linewidth = 1.5, linetype='dashed') +
  geom_line(data = pop2_df.mod, aes(x = salt, y= rate), colour = "blue3", linewidth = 0.5) +
  geom_line(data = pop2_df.mod2, aes(x = salt, y= rate), colour = "blue3", linewidth = 1.5, linetype='dashed') +
  geom_line(data = pop3_df.mod, aes(x = salt, y= rate), colour = "gold", linewidth = 0.5) +
  geom_line(data = pop3_df.mod2, aes(x = salt, y= rate), colour = "gold", linewidth = 1.5, linetype='dashed') +
  
  labs(
    x = "Salinity (g/L)",
    y = "Growth rate",
    title = ""
  ) +
  
  ylim(-0.5,3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.2)

p.s.2

# G matrices --------------------------------------------------------------

###### Let's create the genetic variance plots ######
dsplitnorm0 <- function(x, sd_left = 1, sd_right = 1) {
  # location shift so E[X]=0: mu = -sqrt(2/pi)*(sd_right - sd_left)
  mu <- -sqrt(2/pi) * (sd_right - sd_left)
  c  <- sqrt(2/pi) / (sd_left + sd_right)        # normalizing constant
  ifelse(
    x < mu,
    c * exp(-0.5 * ((x - mu)^2) / sd_left^2),
    c * exp(-0.5 * ((x - mu)^2) / sd_right^2)
  )
}

# ----- make objects -----
A <- list(name = "A", sd_left = 0.8, sd_right = 1.6)
B <- list(name = "B", sd_left = 2, sd_right = 2.2)
C <- list(name = "C", sd_left = 1.0, sd_right = 1.0) 
D <- list(name = "D", sd_left = 0.8, sd_right = 1.1)

# ----- attach data frames -----
z <- seq(-4, 4, length.out = 1200)

A$df <- tibble(z, dens = dsplitnorm0(z, A$sd_left, A$sd_right))
B$df <- tibble(z, dens = dsplitnorm0(z, B$sd_left, B$sd_right))
C$df <- tibble(z, dens = dsplitnorm0(z, C$sd_left, C$sd_right))
D$df <- tibble(z, dens = dsplitnorm0(z, D$sd_left, D$sd_right))

p.v.A <- ggplot(A$df, aes(z, dens)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "C",
       x = "Thermal breadth", y = "Density") +
  ylim(0,0.45) +
  theme_classic()

p.v.A

p.v.C <- ggplot(C$df, aes(z, dens)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "D",
       x = "Thermal breadth", y = "Density") +
  ylim(0,0.45) +
  theme_classic()

p.v.C

p.v.B <- ggplot(B$df, aes(z, dens)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "",
       x = "Salt tolerance", y = "Density") +
  ylim(0,0.45) +
  theme_classic()

p.v.B

p.v.D <- ggplot(D$df, aes(z, dens)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(title = "",
       x = "Salt tolerance", y = "Density") +
  ylim(0,0.45) +
  theme_classic()

p.v.D

cov2 <- function(sd1, sd2, rho){
  matrix(c(sd1^2, rho*sd1*sd2,
           rho*sd1*sd2, sd2^2), 2, 2)
}

## First off-diagonal
sd1  <- 1.2      # SD of trait 1
sd2  <- 0.8      # SD of trait 2
rho  <- 0.9     # correlation between traits (-1..1)
level <- 0.95    # ellipse confidence (e.g., 0.50, 0.95)
mu   <- c(0, 0)  # center (use c(0,0) for mean-standardized traits)

Sigma <- matrix(c(sd1^2, rho*sd1*sd2,
                  rho*sd1*sd2, sd2^2), 2, 2)

ev  <- eigen(Sigma)
rad <- sqrt(qchisq(level, df = 2))
theta <- seq(0, 2*pi, length.out = 400)
circle <- rad * rbind(cos(theta), sin(theta))
A <- ev$vectors %*% diag(sqrt(ev$values))
E <- t(A %*% circle)
ell_df <- data.frame(x = E[,1] + mu[1], y = E[,2] + mu[2])

p.AB.cv <- ggplot(ell_df, aes(x, y)) +
  geom_path(linewidth = 1.3) +
  coord_equal() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Salt tolerance",
       y = "Thermal breadth",
       title = "") +
  xlim(-3.5,3.5) +
  theme_classic()

p.AB.cv


## Second off-diagonal
sd1  <- 1      # SD of trait 1
sd2  <- 1.3      # SD of trait 2
rho  <- -0.4     # correlation between traits (-1..1)
level <- 0.95    # ellipse confidence (e.g., 0.50, 0.95)
mu   <- c(0, 0)  # center (use c(0,0) for mean-standardized traits)

Sigma <- matrix(c(sd1^2, rho*sd1*sd2,
                  rho*sd1*sd2, sd2^2), 2, 2)

ev  <- eigen(Sigma)
rad <- sqrt(qchisq(level, df = 2))
theta <- seq(0, 2*pi, length.out = 400)
circle <- rad * rbind(cos(theta), sin(theta))
A <- ev$vectors %*% diag(sqrt(ev$values))
E <- t(A %*% circle)
ell_df <- data.frame(x = E[,1] + mu[1], y = E[,2] + mu[2])

p.CD.cv <- ggplot(ell_df, aes(x, y)) +
  geom_path(linewidth = 1.3) +
  coord_equal() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Salt tolerance",
       y = "Thermal breadth",
       title = "") +
  xlim(-3.5,3.5) +
  
  theme_classic()

p.CD.cv

###### ok now I will make panels with only covariances?

#the first two work still? NO let's do it fresh

## First off-diagonal
sd1  <- 1.2      # SD of trait 1
sd2  <- 0.8      # SD of trait 2
rho  <- 0.9     # correlation between traits (-1..1)
level <- 0.95    # ellipse confidence (e.g., 0.50, 0.95)
mu   <- c(0, 0)  # center (use c(0,0) for mean-standardized traits)

Sigma <- matrix(c(sd1^2, rho*sd1*sd2,
                  rho*sd1*sd2, sd2^2), 2, 2)

ev  <- eigen(Sigma)
rad <- sqrt(qchisq(level, df = 2))
theta <- seq(0, 2*pi, length.out = 400)
circle <- rad * rbind(cos(theta), sin(theta))
A <- ev$vectors %*% diag(sqrt(ev$values))
E <- t(A %*% circle)
ell_df <- data.frame(x = E[,1] + mu[1], y = E[,2] + mu[2])

c1 <- ggplot(ell_df, aes(x, y)) +
  geom_path(linewidth = 1.3) +
  coord_equal() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "",
       y = "",
       title = "C") +
  xlim(-3.5,3.5) +
  ylim(-3.5,3.5) +
  theme_classic()

c1

### 2

## Second off-diagonal
sd1  <- 1      # SD of trait 1
sd2  <- 1.3      # SD of trait 2
rho  <- -0.4     # correlation between traits (-1..1)
level <- 0.95    # ellipse confidence (e.g., 0.50, 0.95)
mu   <- c(0, 0)  # center (use c(0,0) for mean-standardized traits)

Sigma <- matrix(c(sd1^2, rho*sd1*sd2,
                  rho*sd1*sd2, sd2^2), 2, 2)

ev  <- eigen(Sigma)
rad <- sqrt(qchisq(level, df = 2))
theta <- seq(0, 2*pi, length.out = 400)
circle <- rad * rbind(cos(theta), sin(theta))
A <- ev$vectors %*% diag(sqrt(ev$values))
E <- t(A %*% circle)
ell_df <- data.frame(x = E[,1] + mu[1], y = E[,2] + mu[2])

c2 <- ggplot(ell_df, aes(x, y)) +
  geom_path(linewidth = 1.3) +
  coord_equal() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "",
       y = "",
       title = "") +
  xlim(-3.5,3.5) +
  ylim(-3.5,3.5) +
  theme_classic()

c2

### 3
sd1  <- 1      # SD of trait 1
sd2  <- 1.1      # SD of trait 2
rho  <- -0.1     # correlation between traits (-1..1)
level <- 0.95    # ellipse confidence (e.g., 0.50, 0.95)
mu   <- c(0, 0)  # center (use c(0,0) for mean-standardized traits)

Sigma <- matrix(c(sd1^2, rho*sd1*sd2,
                  rho*sd1*sd2, sd2^2), 2, 2)

ev  <- eigen(Sigma)
rad <- sqrt(qchisq(level, df = 2))
theta <- seq(0, 2*pi, length.out = 400)
circle <- rad * rbind(cos(theta), sin(theta))
A <- ev$vectors %*% diag(sqrt(ev$values))
E <- t(A %*% circle)
ell_df <- data.frame(x = E[,1] + mu[1], y = E[,2] + mu[2])

c3 <- ggplot(ell_df, aes(x, y)) +
  geom_path(linewidth = 1.3) +
  coord_equal() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "",
       y = "",
       title = "") +
  xlim(-3.5,3.5) +
  ylim(-3.5,3.5) +
  theme_classic()

c3

#### 4

sd1  <- 1.4      # SD of trait 1
sd2  <- 1.2      # SD of trait 2
rho  <- 0.45     # correlation between traits (-1..1)
level <- 0.95    # ellipse confidence (e.g., 0.50, 0.95)
mu   <- c(0, 0)  # center (use c(0,0) for mean-standardized traits)

Sigma <- matrix(c(sd1^2, rho*sd1*sd2,
                  rho*sd1*sd2, sd2^2), 2, 2)

ev  <- eigen(Sigma)
rad <- sqrt(qchisq(level, df = 2))
theta <- seq(0, 2*pi, length.out = 400)
circle <- rad * rbind(cos(theta), sin(theta))
A <- ev$vectors %*% diag(sqrt(ev$values))
E <- t(A %*% circle)
ell_df <- data.frame(x = E[,1] + mu[1], y = E[,2] + mu[2])

c4 <- ggplot(ell_df, aes(x, y)) +
  geom_path(linewidth = 1.3) +
  coord_equal() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "",
       y = "",
       title = "") +
  xlim(-3.5,3.5) +
  ylim(-3.5,3.5) +
  theme_classic()

c4

#### 5

sd1  <- 1.1      # SD of trait 1
sd2  <- 1.3      # SD of trait 2
rho  <- 0.7     # correlation between traits (-1..1)
level <- 0.95    # ellipse confidence (e.g., 0.50, 0.95)
mu   <- c(0, 0)  # center (use c(0,0) for mean-standardized traits)

Sigma <- matrix(c(sd1^2, rho*sd1*sd2,
                  rho*sd1*sd2, sd2^2), 2, 2)

ev  <- eigen(Sigma)
rad <- sqrt(qchisq(level, df = 2))
theta <- seq(0, 2*pi, length.out = 400)
circle <- rad * rbind(cos(theta), sin(theta))
A <- ev$vectors %*% diag(sqrt(ev$values))
E <- t(A %*% circle)
ell_df <- data.frame(x = E[,1] + mu[1], y = E[,2] + mu[2])

c5 <- ggplot(ell_df, aes(x, y)) +
  geom_path(linewidth = 1.3) +
  coord_equal() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "",
       y = "",
       title = "") +
  xlim(-3.5,3.5) +
  ylim(-3.5,3.5) +
  theme_classic()

c5

#### 6

sd1  <- 1.2      # SD of trait 1
sd2  <- 1.35      # SD of trait 2
rho  <- -0.2     # correlation between traits (-1..1)
level <- 0.95    # ellipse confidence (e.g., 0.50, 0.95)
mu   <- c(0, 0)  # center (use c(0,0) for mean-standardized traits)

Sigma <- matrix(c(sd1^2, rho*sd1*sd2,
                  rho*sd1*sd2, sd2^2), 2, 2)

ev  <- eigen(Sigma)
rad <- sqrt(qchisq(level, df = 2))
theta <- seq(0, 2*pi, length.out = 400)
circle <- rad * rbind(cos(theta), sin(theta))
A <- ev$vectors %*% diag(sqrt(ev$values))
E <- t(A %*% circle)
ell_df <- data.frame(x = E[,1] + mu[1], y = E[,2] + mu[2])

c6 <- ggplot(ell_df, aes(x, y)) +
  geom_path(linewidth = 1.3) +
  coord_equal() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "",
       y = "",
       title = "") +
  xlim(-3.5,3.5) +
  ylim(-3.5,3.5) +
  theme_classic()

c6

gc1 <- set_panel_size(c1,   width = unit(4, "cm"), height = unit(4, "cm"))
gc2 <- set_panel_size(c2,   width = unit(4, "cm"), height = unit(4, "cm"))
gc3 <- set_panel_size(c3,   width = unit(4, "cm"), height = unit(4, "cm"))
gc4 <- set_panel_size(c4,   width = unit(4, "cm"), height = unit(4, "cm"))
gc5 <- set_panel_size(c5,   width = unit(4, "cm"), height = unit(4, "cm"))
gc6 <- set_panel_size(c6,   width = unit(4, "cm"), height = unit(4, "cm"))


G3 <- plot_grid(
  gc1, gc2, gc3,
  NULL, gc4, gc5,
  NULL, NULL, gc6,
  ncol = 3,
  align = "hv",
  axis = "tblr"
)

G4 <- plot_grid(c1, c2, c3, NULL, c4, c5, NULL, NULL, c6, ncol = 3, align = "hv")

ggsave("figures/66_G4.jpeg", G4, width = 6.5, height = 6.5)

#### And finally, the Pareto front panel?

df.B <- data.frame(
  x = c(5.7, 6.2, 7.3, 7.4, 8.5, 7.8, 8.8, 9.1, 7.9, 9.5, 10.8, 9.9, 10.2, 11.6, 11.4, 11.3, 13.9, 9.1, 12.5, 13.0, 10.8),
  y = c(5.8, 8.1, 6.5, 9.3, 6.3, 8.0, 12.0, 7.5, 10.6, 9.4, 8.2, 10.8, 12.9, 12.1, 7.4, 9.7, 10.1, 14.1, 10.55, 8.7, 13.1)
)

poly.band <- data.frame(
  x = c(5, 7, 15, 9, 5),
  y = c(5, 5, 9, 15, 7)
)

p.B <- ggplot(df.B, aes(x = x, y = y)) +
  
  geom_polygon(data = poly.band, aes(x, y),
               fill = "grey60", alpha = 0.3, colour = NA,
               inherit.aes = FALSE) +
  
  scale_y_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  scale_x_continuous(limits = c(5, 18), breaks = c(5, 10, 15)) +
  
  theme_classic() +
  
  geom_segment(
    aes(x = 5.1, y = 5.1, xend = 11.9, yend = 11.9),
    inherit.aes = FALSE,
    colour = "red3", linewidth = 1.5,
    arrow = arrow(length = unit(2.1, "mm"), ends = "both", type = "closed")) +
  
  geom_segment(
    x = 15,  y = 9,
    xend = 9,   yend = 15,
    colour = "goldenrod2", linewidth = 1.5) +
  
  geom_point(size= 2.5) +
  
  labs(x = "Trait 1 (e.g. thermal breadth)",    
       y = "Trait 2 (e.g. salt tolerance)", 
       title = "D") +  # labels
  
  theme(legend.position = "none") +
  
  annotate("text", x = 12.5, y = 12.5, label = "Pareto front", size = 5.1, fontface = "bold", colour = "goldenrod2", angle = -45) +
  
  coord_cartesian(xlim = c(5,16), ylim = c(5,16), expand = FALSE) 

p.B

ggsave("figures/67_PF.pdf", p.B, width = 6.5, height = 6.5)



















gA <- set_panel_size(p.v.A,   width = unit(4, "cm"), height = unit(4, "cm"))
gB <- set_panel_size(p.CD.cv, width = unit(4, "cm"), height = unit(4, "cm"))
gC <- set_panel_size(p.v.B,   width = unit(4, "cm"), height = unit(4, "cm"))

G1 <- plot_grid(
  gA, gB,
  NULL, gC,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

G1

gD <- set_panel_size(p.v.C,   width = unit(4, "cm"), height = unit(4, "cm"))
gE <- set_panel_size(p.AB.cv, width = unit(4, "cm"), height = unit(4, "cm"))
gF <- set_panel_size(p.v.D,   width = unit(4, "cm"), height = unit(4, "cm"))

G2 <- plot_grid(
  gD, gE,
  NULL, gF,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

G2

gT1 <- set_panel_size(p.t.1,   width = unit(6, "cm"), height = unit(4, "cm"))
gS1 <- set_panel_size(p.s.1, width = unit(6, "cm"), height = unit(4, "cm"))


C1 <- plot_grid(gT1, gS1, ncol = 1, align = "hv")
C1

gT2 <- set_panel_size(p.t.2,   width = unit(6, "cm"), height = unit(4, "cm"))
gS2 <- set_panel_size(p.s.2, width = unit(6, "cm"), height = unit(4, "cm"))

C2 <- plot_grid(gT2, gS2, ncol = 1, align = "hv")
C2

plot_grid(C1, G1,
          C2, G2,
          ncol = 2,
          align = "hv")


ggsave("figures/61_G1.jpeg", G1, width = 4.5, height = 4.5)
ggsave("figures/62_G2.jpeg", G2, width = 4.5, height = 4.5)

ggsave("figures/63_C1.jpeg", C1, width = 6.5, height = 4.5)
ggsave("figures/64_C2.jpeg", C2, width = 6.5, height = 4.5)










G1 <- plot_grid(
  p.v.A, p.AB.cv,
  NULL, p.v.B,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

G1

ggsave("figures/61_G1.jpeg", G1, width = 9, height = 9)


G2 <- plot_grid(
  p.v.C, p.CD.cv,
  NULL, p.v.D,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

G2

C1 <- plot_grid(p.t.1, p.s.1, ncol = 1, align = "hv")
C1

C2 <- plot_grid(p.t.2, p.s.2, ncol = 1, align = "hv")
C2

plot_grid(C1, G1,
          C2, G2,
          ncol = 2,
          align = "hv")


# Salt gradient plot ------------------------------------------------------

pop <- df.s %>% filter(Pop.fac == "1") # Choose your population
pop <- pop %>% 
  mutate(a=a*1.3)

pop1_df.mod <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
  Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.s %>% filter(Pop.fac == "2") # Choose your population
pop <- pop %>% 
  mutate(b=b*1.2)

pop2_df.mod <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
  Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.s %>% filter(Pop.fac == "3") # Choose your population
pop <- pop %>% 
  mutate(c = c + 2.1, b = b*0.84, a = a + 0.3)

pop3_df.mod <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
  Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.s %>% filter(Pop.fac == "4") # Choose your population
pop <- pop %>% 
  mutate(a=a*1)

pop4_df.mod <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
  Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.s %>% filter(Pop.fac == "5") # Choose your population
pop <- pop %>% 
  mutate(b=b*2)

pop5_df.mod <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
  Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )

pop <- df.s %>% filter(Pop.fac == "6") # Choose your population
pop <- pop %>% 
  mutate(b=b*1.1, c = c*1.25)

pop6_df.mod <-
  seq(0, 10, length.out = 200) %>%
  tibble(salt = .) %>%
  mutate(rate = pred_salt(
    salt,
    a = pop$a,
    b = pop$b,
    c = pop$c
  ),
  Pop = as.character(if ("Pop.fac" %in% names(pop)) pop$Pop.fac else "pop")
  )


p.s.1 <- ggplot(pop1_df, aes(x = salt, y =rate)) +
  
  theme_classic() +
  
  geom_line(data = pop1_df.mod, aes(x = salt, y= rate), colour = "black", linewidth = 1.2) +
  geom_line(data = pop2_df.mod, aes(x = salt, y= rate), colour = "black", linewidth = 1.2) +
  geom_line(data = pop3_df.mod, aes(x = salt, y= rate), colour = "black", linewidth = 1.2) +
  geom_line(data = pop4_df.mod, aes(x = salt, y= rate), colour = "black", linewidth = 1.2) +
  geom_line(data = pop5_df.mod, aes(x = salt, y= rate), colour = "black", linewidth = 1.2) +
  geom_line(data = pop6_df.mod, aes(x = salt, y= rate), colour = "black", linewidth = 1.2) +
  
  labs(
    x = "Salinity (g/L)",
    y = "Growth rate",
    title = "A"
  ) +
  
  ylim(-0.5,2.75) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.9) +
  geom_vline(xintercept = 2, linetype = "dashed", size = 0.9)

p.s.1

res <- data.frame(
  Pop = character(),
  rate_at2 = numeric(),
  rate_at0 = numeric(),
  salt_at_halfmax = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:6) {
  nm <- paste0("pop", i, "_df.mod")
  df <- get(nm)
  
  # labels
  pop_label <- if ("Pop" %in% names(df)) unique(df$Pop)[1] else nm
  
  # rate at salt = 2  (linear interpolation on your grid)
  rate2 <- approx(x = df$salt, y = df$rate, xout = 2, rule = 2)$y
  
  # rate at salt = 0  (take this as the max)
  rate0 <- approx(x = df$salt, y = df$rate, xout = 0, rule = 2)$y
  half_target <- 0.5 * rate0
  
  # Find salt where rate = half_target by "inverting" the curve.
  # Ensure monotonic order for interpolation (high rate -> low rate as salt increases)
  ord <- order(df$rate, decreasing = TRUE)
  salt_sorted <- df$salt[ord]
  rate_sorted <- df$rate[ord]
  
  salt_half <- approx(
    x = rate_sorted, y = salt_sorted,
    xout = half_target, rule = 2
  )$y
  
  # add row
  res <- rbind(
    res,
    data.frame(
      Pop = pop_label,
      rate_at2 = rate2,
      rate_at0 = rate0,
      salt_at_halfmax = salt_half,
      stringsAsFactors = FALSE
    )
  )
}

res

res <- res %>%
  mutate(z.c = scale(salt_at_halfmax, center = TRUE, scale = TRUE)[,1]) %>% 
  mutate(w = rate_at2 / mean(rate_at2, na.rm = TRUE))

res2 <- res %>% 
  mutate(salt_at_halfmax = c(0.51, 0.82, 0.95, 1.12, 1.24, 1.39), rate_at2 = c(0.72, 1.25, 1.4, 1.81, 2.09, 2.4))

res2 <- res2 %>%
  mutate(z.c = scale(salt_at_halfmax, center = TRUE, scale = TRUE)[,1]) %>% 
  mutate(w = rate_at2 / mean(rate_at2, na.rm = TRUE))

res <- rbind(res,res2)

res$mic <- c(rep('N', 6), rep('Y', 6)) 

p.s.2 <- ggplot(res, aes(x = z.c, y = w, colour = mic)) +
  
  theme_classic() +
  
  geom_smooth(aes(colour = mic, group = mic), method = "lm", se = FALSE, fullrange = TRUE) +
  
  geom_point(size = 1.5) +
  
  scale_colour_manual(
    name   = "",                  
    breaks = c("N","Y"),                      # legend order
    labels = c("Plants only","With microbiome"),       # legend text
    values = c("N" = "black", "Y" = "red3")  # line + point colors
  ) +

  labs(
    x = "Mean-standardized salt tolerance (c)",
    y = "Relative fitness",
    title = "B"
  ) +
  
  theme(legend.position = c(0.7,0.3)) +
  
  ylim(0,1.75) 

p.s.2

p2 <- plot_grid(p.s.1, p.s.2, align = 'hv')

p2

ggsave("figures/65_fig2.jpeg", p2, width = 9, height = 5)