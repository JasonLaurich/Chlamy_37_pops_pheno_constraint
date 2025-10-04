# Jason R Laurich

# October 2nd, 2025

# I am going to play around with quantile regression

# Load packages & specify functions -----------------------------------------------------------

library(tidyverse)
library(cowplot)
library(quantreg)

# Load & examine the data -------------------------------------------------

df <- read.csv("data-processed/20_summary_table.csv") # Summary file
str(df)
head(df)


df$Evol.plt <- factor(df$Evol, 
                      levels = c("none", "L", "N", "P", "S", "B", "BS", "C"),
                      labels = c("Ancestral",
                                 "Light limitation",
                                 "Nitrogen limitation", 
                                 "Phosphorous limitation", 
                                 "Salt stress", 
                                 "Biotic depletion", 
                                 "Biotic depletion x Salt", 
                                 "Control"))

df$Anc.plt <- factor(df$Anc, 
                     levels = c("anc2", "anc3", "anc4", "anc5", "cc1690"),
                     labels = c("Population 2", 
                                "Population 3", 
                                "Population 4", 
                                "Population 5", 
                                "Mixed population"))

df.p <- read.csv("data-processed/48b_Monod_phosphorous_estimates.csv") # Models fit to replicates, not populations
head(df.p)

# Temperature -------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral", "evolved") # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = T.br.0.56,
    z.x = T.µ.max
  ) # Specify the x and y variables and their 95% CIs

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000)
summary(q75, se = "boot", R = 1000)
summary(q90, se = "boot", R = 1000)

T.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       color = "Evolutionary History",
       title = "E — Temperature") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("evolved" = 1,  # open circle
               "ancestral" = 5)  # diamond
  ) +
  
  ylim(24.5,29.5) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.qr # Display the plot

# Light -------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "L", 'light', 'other')) # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = I.comp.0.56,
    z.x = I.µ.max
  ) # Specify the x and y variables and their 95% CIs

df.filt <- df.filt %>% 
  filter(z.y < 5)

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000)
summary(q75, se = "boot", R = 1000)
summary(q90, se = "boot", R = 1000)

I.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "A — Light") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "light" = 16)  # filled circle
  ) +
  
  ylim(0.025, 0.225) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

I.qr # Display the plot

# Nitrogen ----------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "N", 'nitrogen', 'other')) # for testing significance of matching evolutionary conditions.

df.filt <- df %>% 
  mutate(
    z.y = N.comp.0.56,
    z.x = N.µ.max
  ) # Specify the x and y variables and their 95% CIs

df.filt <- df.filt %>% 
  filter(z.y < 1)

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000)
summary(q75, se = "boot", R = 1000)
summary(q90, se = "boot", R = 1000)

N.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "B — Nitrogen") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "nitrogen" = 16)  # filled circle
  ) +
  
  ylim(0, 0.7) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

N.qr # Display the plot

# Phosphorous -------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "P", 'phosphorous', 'other')) # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = P.comp.0.56,
    z.x = P.µ.max
  ) # Specify the x and y variables and their 95% CIs

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000)
summary(q75, se = "boot", R = 1000)
summary(q90, se = "boot", R = 1000)

P.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       color = "Evolutionary History",
       title = "C — Phosphorous") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "phosphorous" = 16)  # filled circle
  ) +
  
  ylim(0, 5.5) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

P.qr # Display the plot

### Let's run the data fit to replicates

df.p <- df.p %>% 
  filter(1/R.mth <= 10)

q50.ind  <- rq(1/R.mth ~ r.max, tau = 0.50, data = df.p)
q75.ind  <- rq(1/R.mth ~ r.max, tau = 0.75, data = df.p)
q90.ind  <- rq(1/R.mth ~ r.max, tau = 0.90, data = df.p)


summary(q50.ind, se = "boot", R = 1000)
summary(q75.ind, se = "boot", R = 1000)
summary(q90.ind, se = "boot", R = 1000)

P.qr.ind <- ggplot(df.p, aes(x = r.max, y = 1/R.mth)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 1.5, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90.ind)[1], slope = coef(q90.ind)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75.ind)[1], slope = coef(q75.ind)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50.ind)[1], slope = coef(q50.ind)[2], lwd = 1.1, linetype = "dotted") +
  
  ylim(0,4.5) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       color = "Evolutionary History",
       title = "C — Phosphorous") +  # labels
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

P.qr.ind # Display the plot

# Salt --------------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol %in% c("S", "BS"), 'salt', 'other')) # Ss and BSs are treated as equivalent

df.filt <- df %>% 
  mutate(
    z.y = S.c,
    z.x = S.µ.max
  ) # Specify the x and y variables and their 95% CIs

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt)
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt)
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)

summary(q50, se = "boot", R = 1000)
summary(q75, se = "boot", R = 1000)
summary(q90, se = "boot", R = 1000)

null.df <- data.frame(      # Null model results
  m50 = numeric(),          # slopes (q50)
  m75 = numeric(),          # slopes (q75)
  m90 = numeric(),          # slopes (q90)
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt %>%
    
    mutate(
      z.x.sim = sample(df.filt$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt$z.y, replace = FALSE),     # Separately reassign y
    )
  
  q50.n  <- rq(z.y.sim ~ z.x.sim, tau = 0.50, data = shuffled.df)
  q75.n  <- rq(z.y.sim ~ z.x.sim, tau = 0.75, data = shuffled.df)
  q90.n  <- rq(z.y.sim ~ z.x.sim, tau = 0.90, data = shuffled.df)
  
  sum.50.n <- summary(q50.n, se = "boot", R = 1000)
  sum.75.n <- summary(q75.n, se = "boot", R = 1000)
  sum.90.n <- summary(q90.n, se = "boot", R = 1000)
  
  null.df <- rbind(null.df, data.frame(         
    m50 = sum.50.n$coefficients[2,1],  
    m75 = sum.75.n$coefficients[2,1], 
    m90 = sum.90.n$coefficients[2,1] 
  ))
  
}

mean(null.df$m90 <= -1.66817) # p-value 0.084

S.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Salt tolerance (c)", 
       color = "Evolutionary History",
       title = "D — Salt") +  # labels
  
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_shape_manual(
    name = "Evolutionary status",
    values = c("other" = 1,  # open circle
               "ancestral" = 5, # diamond
               "salt" = 16)  # filled circle
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

S.qr # Display the plot

# Compile & save the figures -----------------------------------------------

plots <- list(I.qr, N.qr, P.qr, S.qr, T.qr)

legend_df <- data.frame(
  x = c(1, 2, 1, 2, 1, 2, 1, 2),
  y = c(1, 1, 2, 2, 1, 1, 2, 2),
  Group = factor(c("Ancestral", "Other", "Matching", "Matching", "Ancestral", "Other", "Matching", "Matching")),
  Group2 = factor(c("Biotic depletion", "Biotic depletion x Salt", "Control", "Light limitation", "Nitrogen limitation", "Ancestral", "Phosphorous limitation", "Salt stress")),
  LineType = factor(c("Outer", "Inner (75%)", "Outer", "Inner (75%)", "Outer", "Inner (75%)", "Outer", "Inner (75%)"))
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y)) +
  geom_point(aes(shape = Group, colour = Group2), size = 3, stroke = 1.5) +
  geom_line(aes(linetype = LineType), size = 1) +
  scale_shape_manual(name = NULL,
                     values = c("Ancestral" = 5, 
                                "Other" = 1, 
                                "Matching" = 16)) +
  
  scale_color_manual(
    values = c("Biotic depletion" = "darkorange",
               "Biotic depletion x Salt" = "deepskyblue1",
               "Control" = "forestgreen",
               "Light limitation" = "gold",
               "Nitrogen limitation" = "magenta3",
               "Ancestral" = "black",
               "Phosphorous limitation" = "firebrick",  
               "Salt stress" = "blue")
  ) +
  
  scale_linetype_manual(values = c("Outer" = "dashed", "Inner (75%)" = "solid"),
                        labels = c("Outer", "Inner (75%)"),
                        name = "Pareto front") +
  labs(linetype = "Quantile regression", color = "Evolutionary context") +
  theme_void() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "lines")
  )


legend_only <- get_legend(legend_plot)

all_plots <- c(plots, list(legend_only))

grad_toffs <- plot_grid(plotlist = all_plots,
                        ncol = 2,
                        align = "hv")

ggsave("figures/36_quantreg_intra-gradient_tradeoffs.jpeg", grad_toffs, width = 8, height = 12)
