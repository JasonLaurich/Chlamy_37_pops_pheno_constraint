# Jason R Laurich

# October 1st, 2025

# I am going to investigate possible autocorrelations between model parameters

# Load packages & specify functions -----------------------------------------------------------

library(tidyverse)
library(cowplot)

pred_lact <- function(temp, a, b, delta_t, tmax) {     
  exp(a * temp) - exp(a * tmax - ((tmax - temp) / delta_t)) + b
} # For fitting a lactin curve after modifying variables

pred_mon <- function(nit, r.max, k.s) {
  r.max * nit / (k.s + nit)
} # Ditto for monods

pred_salt <- function(salt, a, b, c) {
  a / (1 + exp(b * (salt - c)))
} # And for salt

# Load & examine the data -------------------------------------------------

df <- read.csv("data-processed/01_µ_estimates_temp.csv")

head(df)
str(df)

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$population.number)
df$well.ID <- as.factor(df$well.ID)
df$Temp <- df$temperature

mat <- split(df, df$pop.num)  # Matrixify the data!

df.med <- mat[[33]] # Going to pull out the median population in terms of thermal breadth at 0.56

head(df.med)

load("R2jags-objects/pop_33_lactin.RData") # The median population for thermal breadth

df.jags <- data.frame(lac_jag$BUGSoutput$summary)

head(df.jags)

df.med.n <- read.csv("data-processed/07a_µ_estimates_nitrogen.csv")

head(df.med.n)
df.med.n$nit <- df.med.n$nitrate_concentration

mat.n <- split(df.med.n, df.med.n$population.number)  # Matrixify the data!

df.med.n <- mat.n[[33]] # Going to pull out the median population in terms of thermal breadth at 0.56

head(df.med.n)

df.s <- read.csv("data-processed/09a_µ_estimates_salt.csv")

head(df.s)

mat.s <- split(df.s, df.s$population.number)  # Matrixify the data!

df.med.s <- mat.s[[33]] # Going to pull out the median population in terms of thermal breadth at 0.56

head(df.med.s)

# Temperature ----------------------------------------------------

###### let's start by modifying the underlying estimates directly  ######

df.0 <- df.jags[-c(1:6),]
df.0$temp <- seq(0, 45, 0.05)

df.1 <- df.0%>% 
  mutate(mean = mean*2)

df.2 <- df.0%>% 
  mutate(mean = mean*3)

df.3 <- df.0%>% 
  mutate(mean = mean*0.5)

p.mod.Tbr <- ggplot(df.med, aes(x = temperature, y = r.exp)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) +
  
  geom_line(data = df.0, aes(x = temp, y= mean), colour = "black", linewidth = 1.5) +
  geom_line(data = df.1, aes(x = temp, y= mean), colour = "goldenrod2", linewidth = 1.5) +
  geom_line(data = df.2, aes(x = temp, y= mean), colour = "red3", linewidth = 1.5) +
  geom_line(data = df.3, aes(x = temp, y= mean), colour = "blue3", linewidth = 1.5) +
  
  ylim(-1,12) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "Comparison of models with varying µmax") +
  
  geom_hline(yintercept = 0.56, linetype = "dashed",
             colour = "black", linewidth = 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.mod.Tbr

dfs <- list(df.0, df.1, df.2, df.3)

summ.df.tbr <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  r.max.raw = numeric(),    # Maximum growth rate (Jags raw)
  T.br.raw.056 = numeric(),     # T breadth (Jags raw)
  T.br.raw.0 = numeric(),     # T breadth (Jags raw)
  T.br.raw.half = numeric(),     # T breadth at half µmax
  stringsAsFactors = FALSE  # Avoid factor conversion
)

for (i in 1:4){
  
  df.i <- dfs[[i]]
  
  summ.df.tbr <- rbind(summ.df.tbr, data.frame(                          
    r.max.raw = max(df.i$mean),                                                                      # Maximum growth rate
    T.br.raw.056 = df.i$temp[max(which(df.i$mean > 0.56))] - df.i$temp[min(which(df.i$mean > 0.56))],     # T breadth
    T.br.raw.0 = df.i$temp[max(which(df.i$mean > 0))] - df.i$temp[min(which(df.i$mean > 0))],
    T.br.raw.half = df.i$temp[max(which(df.i$mean > (max(df.i$mean) / 2)))] - 
      df.i$temp[min(which(df.i$mean > (max(df.i$mean) / 2)))]                                            # T breadth at half the max growth rate
  ))
  
} 

lm.auto.temp.1 <- lm(T.br.raw.056 ~ r.max.raw, data = summ.df.tbr)
summary(lm.auto.temp.1) # adj. R^2 0.7429, p 0.08974

p.auto.temp.1.056 <- ggplot(summ.df.tbr, aes(x = r.max.raw, y = T.br.raw.056)) +
  geom_point(size= 3) +
  geom_smooth(method = 'lm') +
  
  theme_classic() +
  
  labs(x = "Maximum growth rate (µ)", 
       y = "Thermal breadth (C)",
       title = "A - effects of growth modification on T.br at 0.56") +
  
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.auto.temp.1.056

p.auto.temp.1.0 <- ggplot(summ.df.tbr, aes(x = r.max.raw, y = T.br.raw.0)) +
  geom_point(size= 3) +
  geom_smooth(method = 'lm') +
  
  theme_classic() +
  
  labs(x = "Maximum growth rate (µ)", 
       y = "Thermal breadth (C)",
       title = "B - effects of growth modification on T.br at 0") +
  
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.auto.temp.1.0

p.auto.temp.1.half <- ggplot(summ.df.tbr, aes(x = r.max.raw, y = T.br.raw.half)) +
  geom_point(size= 3) +
  geom_smooth(method = 'lm') +
  
  theme_classic() +
  
  labs(x = "Maximum growth rate (µ)", 
       y = "Thermal breadth (C)",
       title = "C - effects of growth modification on T.br at half µmax") +
  
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.auto.temp.1.half

###### modify cf.a ######

df.a.0 <- df.jags[c(1:6),]

df.a.1 <- df.a.0%>% 
  mutate(mean = replace(mean, 1L, mean[1] * 1.25))

df.a.1.pred <- seq(0, 45, length.out = 901) %>%
  tibble(temp = .) %>%
  mutate(rate = pred_lact(
     temp,
    a = df.a.1$mean[1],
    b = df.a.1$mean[2],
    delta_t = df.a.1$mean[3],
    tmax = df.a.1$mean[5]
    ))

df.a.2 <- df.a.0%>% 
  mutate(mean = replace(mean, 1L, mean[1] * 1.5))

df.a.2.pred <- seq(0, 45, length.out = 901) %>%
  tibble(temp = .) %>%
  mutate(rate = pred_lact(
    temp,
    a = df.a.2$mean[1],
    b = df.a.2$mean[2],
    delta_t = df.a.2$mean[3],
    tmax = df.a.2$mean[5]
  ))

df.a.3 <- df.a.0%>% 
  mutate(mean = replace(mean, 1L, mean[1] * 0.75))

df.a.3.pred <- seq(0, 45, length.out = 901) %>%
  tibble(temp = .) %>%
  mutate(rate = pred_lact(
    temp,
    a = df.a.3$mean[1],
    b = df.a.3$mean[2],
    delta_t = df.a.3$mean[3],
    tmax = df.a.3$mean[5]
  ))

p.mod.Tbr.a <- ggplot(df.med, aes(x = temperature, y = r.exp)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) +
  
  geom_line(data = df.0, aes(x = temp, y= mean), colour = "black", linewidth = 1.5) +
  geom_line(data = df.a.1.pred, aes(x = temp, y= rate), colour = "goldenrod2", linewidth = 1.5) +
  geom_line(data = df.a.2.pred, aes(x = temp, y= rate), colour = "red3", linewidth = 1.5) +
  geom_line(data = df.a.3.pred, aes(x = temp, y= rate), colour = "blue3", linewidth = 1.5) +
  
  ylim(-1,12) +
  
  labs(x = "Temperature (°C)", 
       y = "Exponential growth rate (µ)",
       title = "Comparison of models with varying µmax, modifying a") +
  
  geom_hline(yintercept = 0.56, linetype = "dashed",
             colour = "black", linewidth = 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.mod.Tbr.a

df.a.1.pred <- as.data.frame(df.a.1.pred) %>% 
  rename(mean = rate)

df.a.2.pred <- as.data.frame(df.a.2.pred) %>% 
  rename(mean = rate)

df.a.3.pred <- as.data.frame(df.a.3.pred) %>% 
  rename(mean = rate)

dfs.a <- list(df.0, df.a.1.pred, df.a.2.pred, df.a.3.pred)

summ.df.tbr.a <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  r.max.raw = numeric(),    # Maximum growth rate (Jags raw)
  T.br.raw.056 = numeric(),     # T breadth (Jags raw)
  T.br.raw.0 = numeric(),     # T breadth (Jags raw)
  T.br.raw.half = numeric(),     # T breadth at half µmax
  stringsAsFactors = FALSE  # Avoid factor conversion
)

for (i in 1:4){
  
  df.i <- dfs.a[[i]]
  
  summ.df.tbr.a <- rbind(summ.df.tbr.a, data.frame(                          
    r.max.raw = max(df.i$mean),                                                                      # Maximum growth rate
    T.br.raw.056 = df.i$temp[max(which(df.i$mean > 0.56))] - df.i$temp[min(which(df.i$mean > 0.56))],     # T breadth
    T.br.raw.0 = df.i$temp[max(which(df.i$mean > 0))] - df.i$temp[min(which(df.i$mean > 0))],
    T.br.raw.half = df.i$temp[max(which(df.i$mean > (max(df.i$mean) / 2)))] - 
      df.i$temp[min(which(df.i$mean > (max(df.i$mean) / 2)))]                                            # T breadth at half the max growth rate
  ))
  
} 

lm.auto.temp.a.1 <- lm(T.br.raw.056 ~ r.max.raw, data = summ.df.tbr.a)
summary(lm.auto.temp.a.1) # adj. R^2 0.7429, p 0.08974

p.auto.temp.1.056.a <- ggplot(summ.df.tbr.a, aes(x = r.max.raw, y = T.br.raw.056)) +
  geom_point(size= 3) +
  geom_smooth(method = 'lm') +
  
  theme_classic() +
  
  labs(x = "Maximum growth rate (µ)", 
       y = "Thermal breadth (C)",
       title = "A - effects of growth modification on T.br at 0.56") +
  
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.auto.temp.1.056.a

p.auto.temp.1.0.a <- ggplot(summ.df.tbr.a, aes(x = r.max.raw, y = T.br.raw.0)) +
  geom_point(size= 3) +
  geom_smooth(method = 'lm') +
  
  theme_classic() +
  
  labs(x = "Maximum growth rate (µ)", 
       y = "Thermal breadth (C)",
       title = "B - effects of growth modification on T.br at 0") +
  
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.auto.temp.1.0.a

p.auto.temp.1.half.a <- ggplot(summ.df.tbr.a, aes(x = r.max.raw, y = T.br.raw.half)) +
  geom_point(size= 3) +
  geom_smooth(method = 'lm') +
  
  theme_classic() +
  
  labs(x = "Maximum growth rate (µ)", 
       y = "Thermal breadth (C)",
       title = "C - effects of growth modification on T.br at half µmax") +
  
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.auto.temp.1.half.a

# Nitrogen ----------------------------------------------------------------

load("R2jags-objects/pop_33_nitrogen_monod.RData") # The median population for thermal breadth

df.n <- data.frame(monod_jag$BUGSoutput$summary)

head(df.n)

df.n.0 <- df.n[c(1:3,2005),]

df.n.0.pred <- seq(0, 1000, length.out = 2001) %>%
  tibble(nit = .) %>%
  mutate(rate = pred_mon(
    nit,
    r.max = df.n.0$mean[3],
    k.s = df.n.0$mean[1]
  ))

df.n.1 <- df.n.0%>% 
  mutate(mean = replace(mean, 3L, mean[3] * 2))

df.n.1.pred <- seq(0, 1000, length.out = 2001) %>%
  tibble(nit = .) %>%
  mutate(rate = pred_mon(
    nit,
    r.max = df.n.1$mean[3],
    k.s = df.n.1$mean[1]
  ))

df.n.2 <- df.n.0%>% 
  mutate(mean = replace(mean, 3L, mean[3] * 3))

df.n.2.pred <- seq(0, 1000, length.out = 2001) %>%
  tibble(nit = .) %>%
  mutate(rate = pred_mon(
    nit,
    r.max = df.n.2$mean[3],
    k.s = df.n.2$mean[1]
  ))


df.n.3 <- df.n.0%>% 
  mutate(mean = replace(mean, 3L, mean[3] * 0.5))

df.n.3.pred <- seq(0, 1000, length.out = 2001) %>%
  tibble(nit = .) %>%
  mutate(rate = pred_mon(
    nit,
    r.max = df.n.3$mean[3],
    k.s = df.n.3$mean[1]
  ))

p.mod.n <- ggplot(df.med.n, aes(x = nitrate.lvl, y = r.exp)) +
  geom_jitter(size = 2.5,width = 0.5, height = 0) +
  
  geom_line(data = df.n.0.pred, aes(x = nit, y= rate), colour = "black", linewidth = 1.5) +
  geom_line(data = df.n.1.pred, aes(x = nit, y= rate), colour = "goldenrod2", linewidth = 1.5) +
  geom_line(data = df.n.2.pred, aes(x = nit, y= rate), colour = "red3", linewidth = 1.5) +
  geom_line(data = df.n.3.pred, aes(x = nit, y= rate), colour = "blue3", linewidth = 1.5) +
  
  ylim(-1,5) +
  
  labs(x = "Nitrogen concentration", 
       y = "Exponential growth rate (µ)",
       title = "Comparison of models with varying µmax") +
  
  geom_hline(yintercept = 0.56, linetype = "dashed",
             colour = "black", linewidth = 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.mod.n

dfs.n <- list(df.n.0.pred, df.n.1.pred, df.n.2.pred, df.n.3.pred)

summ.df.tbr.n <- data.frame(  # We'll create a dataframe to store the data as we fit models.
  r.max.raw = numeric(),      # Maximum growth rate (Jags raw)
  R = numeric(),              # R*
  stringsAsFactors = FALSE    # Avoid factor conversion
)

for (i in 1:4){
  
  df.i <- dfs.n[[i]]
  
  summ.df.tbr.n <- rbind(summ.df.tbr.n, data.frame(                          
    r.max.raw = max(df.i$rate),                        # Maximum growth rate
    R = df.i$nit[min(which(df.i$rate > 0.56))]        # R
  ))
  
} 

lm.auto.nit <- lm(R ~ r.max.raw, data = summ.df.tbr.n)
summary(lm.auto.nit) # adj. R^2 0.3689, p 0.2389

p.auto.nit <- ggplot(summ.df.tbr.n, aes(x = r.max.raw, y = R)) +
  geom_point(size= 3) +
  geom_smooth(method = 'lm') +
  
  theme_classic() +
  
  labs(x = "Maximum growth rate (µ)", 
       y = "R* (nitrogen)",
       title = "B - effects of growth modification on N*") +
  
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.auto.nit

lm.auto.nit.2 <- lm(R ~ r.max.raw, data = summ.df.tbr.n[1:3,])
summary(lm.auto.nit.2) # adj. R^2 0.7857, p 0.2123

p.auto.nit.2 <- ggplot(summ.df.tbr.n[1:3,], aes(x = r.max.raw, y = R)) +
  geom_point(size= 3) +
  geom_smooth(method = 'lm') +
  
  theme_classic() +
  
  labs(x = "Maximum growth rate (µ)", 
       y = "R* (nitrogen)",
       title = "B - effects of growth modification on N*") +
  
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.auto.nit.2

p.auto.nit.3 <- ggplot(summ.df.tbr.n, aes(x = r.max.raw, y = 1/R)) +
  geom_point(size= 3) +
  geom_smooth(method = 'lm') +
  
  theme_classic() +
  
  labs(x = "Maximum growth rate (µ)", 
       y = "R* (nitrogen)",
       title = "B - effects of growth modification on 1/N*") +
  
  theme(
    legend.position = "none",  # delete legend
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)
  )

p.auto.nit.3

# Salt --------------------------------------------------------------------


