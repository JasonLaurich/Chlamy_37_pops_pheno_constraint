# Jason R Laurich
# November 25th, 2025

# OK we're chasing down the final few models that are still not performing.  

# Packages & functions ----------------------------------------------------

library(tidyr)
library(dplyr)
library(R2jags)
library(mcmcplots)
library(cowplot)

pred_salt <- function(salt, a, b, c) {
  a / (1 + exp(b * (salt - c)))
}

# Load & examine the data ------------------------------------------------------------

df.r <- read.csv("data-processed/203_µ_estimates_salt_new.csv")

head(df.r)
str(df.r)

df.r$pop.fac <- as.factor(df.r$population)
df.r$pop.num <- as.numeric(df.r$population.number)
df.r$well.ID <- as.factor(df.r$well.ID)
df.r$salt <- as.numeric(df.r$salt)

salt <- as.vector(as.numeric(as.character(unique(df.r$salt)))) # for looping through salt levels
ord.salt<- sort(salt)

mat <- split(df.r, df.r$pop.num)  # Matrixify the data!

# Fit tolerance curves ----------------------------------------------------

inits.salt <- function() { # Original
  list(
    a = runif(1, 0.1, 3.5),  # Initial guess for a
    b = runif(1, 0.1, 5),  # Initial guess for b
    c = runif(1, 0.1, max(df.i$salt)),  # Initial guess for c
    sigma = runif(1, 0.1, 2)  # Initial guess for error
  )
}

inits.salt.smalla <- function() { # Smaller a (prior 0.5-> 2)
  list(
    a = runif(1, 0.7, 1.8),  # Smaller window for a
    b = runif(1, 0.1, 5),  # Initial guess for b
    c = runif(1, 0.1, max(df.i$salt)),  # Initial guess for c
    sigma = runif(1, 0.1, 2)  # Initial guess for error
  )
}

parameters.salt <- c("a", "b", "c", "sigma", "r_pred_new") # Repeats

S.pred <- seq(0, 10, 0.005) # Repeats
N.S.pred <-length(S.pred)

# Repeated here
ni.fit <- 330000   # iterations / chain
nb.fit <- 30000    # burn in periods for each chain
nt.fit <- 300      # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 6        # number of chains, total of 6,000 estimates for each model.

# Population 9 ------------------------------------------------------------
# number 32

i<- 32

###### Original (load it) #######

load(paste0("R2jags-objects/pop_", i, "_salt_newer.RData"))
assign(paste0("monod_jag", 9), monod_jag)

head(monod_jag.smalla$BUGSoutput$summary)

###### Older - I think c was still at 10 here? (load it) #######

load(paste0("R2jags-objects/pop_", i, "_salt_new.RData"))
assign(paste0("monod_jag", 9, "oldc"), monod_jag)

head(monod_jag9oldc$BUGSoutput$summary)

###### Smaller a ######

df.i <- subset(mat[[i]])

df.i <- droplevels(df.i)
  
trait <- df.i$r.exp     # format the data for jags
N.obs <- length(trait)
salt <- df.i$salt
  
jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)
  
monod_jag.smalla <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.smalla.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag.smalla$BUGSoutput$summary)

###### Intermediate a ######

monod_jag.meda <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.meda.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag.meda$BUGSoutput$summary)

###### Smaller a, medium c ######

monod_jag.smalla.medc <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.smalla.medc.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag.smalla.medc$BUGSoutput$summary)

###### Smaller a, lower c ######

monod_jag.smalla.lowc <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.smalla.lowc.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag.smalla.lowc$BUGSoutput$summary)

###### Smaller a, small c ######

monod_jag.smalla.smallc <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.smalla.smallc.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag.smalla.smallc$BUGSoutput$summary)

###### Plotting ######

s.mu <- df.r %>% filter(population.number == i) # Raw data

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag9$BUGSoutput$summary[1,1], b = monod_jag9$BUGSoutput$summary[2,1], c = monod_jag9$BUGSoutput$summary[3,1])
)

p.og <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "A: original priors"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.og

### Older one - c still cut off at 10 I think?

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag9oldc$BUGSoutput$summary[1,1], b = monod_jag9oldc$BUGSoutput$summary[2,1], c = monod_jag9oldc$BUGSoutput$summary[3,1])
)

p.ogold <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "B: original priors, c < 10"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.ogold

### Smaller a

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag.smalla$BUGSoutput$summary[1,1], b = monod_jag.smalla$BUGSoutput$summary[2,1], c = monod_jag.smalla$BUGSoutput$summary[3,1])
)

p.smalla <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "C: smaller a (0.5-2)"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.smalla

### Intermediate a

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag.meda$BUGSoutput$summary[1,1], b = monod_jag.meda$BUGSoutput$summary[2,1], c = monod_jag.meda$BUGSoutput$summary[3,1])
)

p.meda <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "D: intermediate a (0.5-2.5)"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.meda

### Smaller a, medium c

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag.smalla.medc$BUGSoutput$summary[1,1], b = monod_jag.smalla.medc$BUGSoutput$summary[2,1], c = monod_jag.smalla.medc$BUGSoutput$summary[3,1])
)

p.smalla.medc <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "E: smaller a (0.5-2), medium c (0-15)"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.smalla.medc

### Smaller a, low c

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag.smalla.lowc$BUGSoutput$summary[1,1], b = monod_jag.smalla.lowc$BUGSoutput$summary[2,1], c = monod_jag.smalla.lowc$BUGSoutput$summary[3,1])
)

p.smalla.lowc <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "F: smaller a (0.5-2), lower c (0-12)"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.smalla.lowc

### Smaller a, small c

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag.smalla.smallc$BUGSoutput$summary[1,1], b = monod_jag.smalla.smallc$BUGSoutput$summary[2,1], c = monod_jag.smalla.smallc$BUGSoutput$summary[3,1])
)

p.smalla.smallc <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "G: smaller a (0.5-2), low c (0-10)"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.smalla.smallc

plot_grid(p.og, p.ogold, p.smalla, p.meda, p.smalla.medc, p.smalla.lowc, p.smalla.smallc, nrow = 3, align='hv', rel_widths = c(1,1,1))

# Model output comparison -------------------------------------------------

head(monod_jag$BUGSoutput$summary)
head(monod_jag.smalla$BUGSoutput$summary)
head(monod_jag.meda$BUGSoutput$summary)
head(monod_jag.smalla.medc$BUGSoutput$summary)

# DIC comparison ----------------------------------------------------------

monod_jag9$BUGSoutput$DIC # 18.4055
monod_jag9oldc$BUGSoutput$DIC # 8.745938
monod_jag.smalla$BUGSoutput$DIC # 21.2328
monod_jag.meda$BUGSoutput$DIC # 19.98966
monod_jag.smalla.medc$BUGSoutput$DIC # 20.20253
monod_jag.smalla.lowc$BUGSoutput$DIC # 17.49764
monod_jag.smalla.smallc$BUGSoutput$DIC # 14.49764


# Population 2 ------------------------------------------------------------
# number 10

i<- 10

###### Original (load it) #######

load(paste0("R2jags-objects/pop_", i, "_salt_newer.RData"))
assign(paste0("monod_jag", 2), monod_jag)

head(monod_jag2$BUGSoutput$summary)

###### Older - I think c was still at 10 here? (load it) #######

load(paste0("R2jags-objects/pop_", i, "_salt_new.RData"))
assign(paste0("monod_jag", 2, "oldc"), monod_jag)

head(monod_jag2oldc$BUGSoutput$summary)

###### Smaller a ######

df.i <- subset(mat[[i]])

df.i <- droplevels(df.i)

trait <- df.i$r.exp     # format the data for jags
N.obs <- length(trait)
salt <- df.i$salt

jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)

monod_jag.2.smalla <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.smalla.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag.2.smalla$BUGSoutput$summary)

###### Intermediate a ######

monod_jag.2.meda <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.meda.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag.2.meda$BUGSoutput$summary)

###### Smaller a, medium c ######

monod_jag.2.smalla.medc <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.smalla.medc.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag.2.smalla.medc$BUGSoutput$summary)

###### Smaller a, lower c ######

monod_jag.2.smalla.lowc <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.smalla.lowc.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag.2.smalla.lowc$BUGSoutput$summary)

###### Smaller a, small c ######

monod_jag.2.smalla.smallc <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.smalla.smallc.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag.2.smalla.smallc$BUGSoutput$summary)


###### Plotting ######

s.mu <- df.r %>% filter(population.number == i) # Raw data

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag2$BUGSoutput$summary[1,1], b = monod_jag2$BUGSoutput$summary[2,1], c = monod_jag2$BUGSoutput$summary[3,1])
)

p.og2 <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "A: original priors"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.og2

### Older one - c still cut off at 10 I think?

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag2oldc$BUGSoutput$summary[1,1], b = monod_jag2oldc$BUGSoutput$summary[2,1], c = monod_jag2oldc$BUGSoutput$summary[3,1])
)

p.ogold2 <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "B: original priors, c < 10"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.ogold2


### Smaller a

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag.2.smalla$BUGSoutput$summary[1,1], b = monod_jag.2.smalla$BUGSoutput$summary[2,1], c = monod_jag.2.smalla$BUGSoutput$summary[3,1])
)

p.smalla2 <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "C: smaller a (0.5-2)"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.smalla2

### Intermediate a

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag.2.meda$BUGSoutput$summary[1,1], b = monod_jag.2.meda$BUGSoutput$summary[2,1], c = monod_jag.2.meda$BUGSoutput$summary[3,1])
)

p.meda2 <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "D: intermediate a (0.5-2.5)"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.meda2

### Smaller a, medium c

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag.2.smalla.medc$BUGSoutput$summary[1,1], b = monod_jag.2.smalla.medc$BUGSoutput$summary[2,1], c = monod_jag.2.smalla.medc$BUGSoutput$summary[3,1])
)

p.smalla.medc2 <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "E: smaller a (0.5-2), medium c (0-15)"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.smalla.medc2

### Smaller a, low c

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag.2.smalla.lowc$BUGSoutput$summary[1,1], b = monod_jag.2.smalla.lowc$BUGSoutput$summary[2,1], c = monod_jag.2.smalla.lowc$BUGSoutput$summary[3,1])
)

p.smalla.lowc2 <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "F: smaller a (0.5-2), lower c (0-12)"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.smalla.lowc2

plot_grid(p.og2, p.smalla2, p.meda2, p.smalla.medc2, p.smalla.lowc2, nrow = 3, align='hv', rel_widths = c(1,1))

### Smaller a, small c

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag.2.smalla.smallc$BUGSoutput$summary[1,1], b = monod_jag.2.smalla.smallc$BUGSoutput$summary[2,1], c = monod_jag.2.smalla.smallc$BUGSoutput$summary[3,1])
)

p.smalla.smallc2 <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "G: smaller a (0.5-2), low c (0-10)"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p.smalla.smallc2

plot_grid(p.og2, p.ogold2, p.smalla2, p.meda2, p.smalla.medc2, p.smalla.lowc2, p.smalla.smallc2, nrow = 3, align='hv', rel_widths = c(1,1,1))

# Model output comparison -------------------------------------------------

head(monod_jag2$BUGSoutput$summary)
head(monod_jag2oldc$BUGSoutput$summary)
head(monod_jag.2.smalla$BUGSoutput$summary)
head(monod_jag.2.meda$BUGSoutput$summary)
head(monod_jag.2.smalla.medc$BUGSoutput$summary)
head(monod_jag.2.smalla.lowc$BUGSoutput$summary)
head(monod_jag.2.smalla.smallc$BUGSoutput$summary)

# DIC comparison ----------------------------------------------------------

monod_jag2$BUGSoutput$DIC # 18.4055
monod_jag2oldc$BUGSoutput$DIC # 18.4055
monod_jag.2.smalla$BUGSoutput$DIC # 21.2328
monod_jag.2.meda$BUGSoutput$DIC # 19.98966
monod_jag.2.smalla.medc$BUGSoutput$DIC # 20.20253
monod_jag.2.smalla.lowc$BUGSoutput$DIC # 20.20253
monod_jag.2.smalla.smallc$BUGSoutput$DIC # 20.20253

# Try out the final model on other populations? ---------------------------

### 15

i <- 5

df.i <- subset(mat[[i]])

df.i <- droplevels(df.i)

trait <- df.i$r.exp     # format the data for jags
N.obs <- length(trait)
salt <- df.i$salt

jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)

monod_jag15.smalla.smallc <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.smalla.smallc.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag15.smalla.smallc$BUGSoutput$summary)

s.mu <- df.r %>% filter(population.number == i) # Raw data

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag15.smalla.smallc$BUGSoutput$summary[1,1], b = monod_jag15.smalla.smallc$BUGSoutput$summary[2,1], c = monod_jag15.smalla.smallc$BUGSoutput$summary[3,1])
)

p15 <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "15: small a, small c"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p15

### 14

i <- 4

df.i <- subset(mat[[i]])

df.i <- droplevels(df.i)

trait <- df.i$r.exp     # format the data for jags
N.obs <- length(trait)
salt <- df.i$salt

jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)

monod_jag14.smalla.smallc <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.smalla.smallc.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag14.smalla.smallc$BUGSoutput$summary)

s.mu <- df.r %>% filter(population.number == i) # Raw data

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag14.smalla.smallc$BUGSoutput$summary[1,1], b = monod_jag14.smalla.smallc$BUGSoutput$summary[2,1], c = monod_jag14.smalla.smallc$BUGSoutput$summary[3,1])
)

p14 <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "14: small a, small c"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p14

### 18

i <- 8

df.i <- subset(mat[[i]])

df.i <- droplevels(df.i)

trait <- df.i$r.exp     # format the data for jags
N.obs <- length(trait)
salt <- df.i$salt

jag.data <- list(trait = trait, N.obs = N.obs, S = salt, S.pred = S.pred, N.S.pred = N.S.pred)

monod_jag18.smalla.smallc <- jags( # Run the salt logistic growth curve function. 
  data = jag.data,
  inits = inits.salt.smalla,
  parameters.to.save = parameters.salt,
  model.file = "salt.tol.smalla.smallc.txt",
  n.thin = nt.fit,
  n.chains = nc.fit,
  n.burnin = nb.fit,
  n.iter = ni.fit,
  DIC = TRUE,
  working.directory = getwd()
)

head(monod_jag18.smalla.smallc$BUGSoutput$summary)

s.mu <- df.r %>% filter(population.number == i) # Raw data

curve.s <- tibble::tibble(
  salt  = seq(0, 10, length.out = 200),
  rate = pred_salt(salt, a = monod_jag18.smalla.smallc$BUGSoutput$summary[1,1], b = monod_jag18.smalla.smallc$BUGSoutput$summary[2,1], c = monod_jag18.smalla.smallc$BUGSoutput$summary[3,1])
)

p18 <- ggplot(curve.s, aes(x = salt, y = rate)) +
  geom_line(size = 1.5, colour = "orchid2") +
  geom_point(data = s.mu,
             aes(x = salt, y = r.exp),
             inherit.aes = FALSE, size = 2) +
  labs(
    x = "Salt concentration (g/L)",
    y = "Exponential growth rate",
    title = "18: small a, small c"
  ) +
  theme_classic() +
  ylim(-0.1, 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1.2)

p18

plot_grid(p14, p18, p15, nrow = 1, align='hv', rel_widths = c(1,1,1))
