# Jason R Laurich
# March 1, 2025

# Moving the analysis portion of script 11 here. It was getting unwieldy
# Modelling gleaner-opportunist, specialist-generalist trade-offs
# Adding Pareto frontier estimation to correlation plots. 

# Conceptual approach - (1) use modified convex Hull algorithm (upper right only) to identify potential Pareto frontier.
# (2) Model regression of full data set, evolved populations, ancestral using lmer/brms - are the slopes of the evolved populations different? Does the intercept differ?
# (3) Model regression using non-linear regression - same question - differences in shape (especially where evolved populations form the PF) would be strong evidence.

# So partition out ancestors, evolutionary conditions particular to the nutrient in question, and other evolutionary treatments?

# Another angle for confirming the Pareto Frontier would be bootstrapping the data, and asking whether, on average, the Pareto front is more likely
# to contain points that evolved under relevant conditions?

# Or can I estimate the Pareto Frontier using quantile regression at various levels of sensitivity for this analysis?
# And then ask whether the evolved populations relevant to a particular nutrient are above some threshold QR?
# E.g. If half of the data points are evolutionary relevant, they should be above the 50% QR?
# Confirm by randomizing and evaluating the distribution of evoluationarily relevant points under a null distribution?

############# Packages ########################

library(dplyr)
library(ggplot2)
library(rPref)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(vegan)  # For PCA and RDA
library(ggrepel)
library(quantreg)
library(brms)

############# Upload and organize data #######################

df<-read.csv("data-processed/14_summary_metric_table.csv") # Summary file

par_frt <- function(df, xvar, yvar) { # Simple Pareto frontier function
  
  df <- df[order(-df[[xvar]], df[[yvar]]), ]  
  pareto_points <- df[1, ]  # Start with the first point
  
  for (i in 2:nrow(df)) {
    if (df[i, yvar] > tail(pareto_points[[yvar]], 1)) {  # Ensure increasing y values
      pareto_points <- rbind(pareto_points, df[i,])
    }
  }
  
  return(pareto_points)
}

############# Temperature ####################################

df$shape <- ifelse(df$evol == "none", 22, 16) # I want to add a shape column to the dataframe that I will update
# The idea is to label un-evolved populations with a square, and then later (not for T) relevant experimental evolution nutrient conditions with a star

df$evol.bin <- ifelse(df$evol == "none", "ancestral", "evolved") # For regressions

par.res.T <- par_frt(df, xvar = "r.max_T", yvar = "T.br")

T_par <- ggplot(df, aes(x = r.max_T, y = T.br)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Thermal breadth", 
       title = "Thermal performance") +
  scale_shape_manual(values = c(16, 3)) +  # Keep custom shapes
  geom_line(data = par.res.T, aes(x = r.max_T, y = T.br), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

T_par # Raw pareto front.

mod.T <- lm(T.br~r.max_T*evol.bin, data=df)  # test significance of Pareto frontier
summary(mod.T)

T.regs <- ggplot(df, aes(x = r.max_T, y = T.br, color = evol.bin)) +  
  geom_point(size = 2) +  # Scatter plot of raw data
  geom_smooth(method = "lm", se = FALSE, aes(color = evol.bin)) +  # Separate regression lines
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Evolutionary Effects on Thermal Performance Trade-offs",
       color = "Evolutionary History") +  # Change legend title) +
  scale_color_manual(values = c("black", "goldenrod1"), 
                     labels = c("Ancestral", "Evolved")) +  # Custom colors per group
  theme_classic() +
  theme(
    legend.position = c(0.8, 0.85),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 10, face = "bold"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10)
  ) 

T.regs

pred.t <- data.frame(r.max_T = seq(min(df$r.max_T), max(df$r.max_T), length.out = 100)) # Dataframe to collect quantile info in

quant.T.100 <- rq(T.br ~ poly(r.max_T, 2), data = df, tau = 1.00) 
pred.t$T.br.100 <- predict(quant.T.100, newdata = pred.t)

quant.T.95 <- rq(T.br ~ poly(r.max_T, 2), data = df, tau = 0.95) 
pred.t$T.br.95 <- predict(quant.T.95, newdata = pred.t)

quant.T.90 <- rq(T.br ~ poly(r.max_T, 2), data = df, tau = 0.90) 
pred.t$T.br.90 <- predict(quant.T.90, newdata = pred.t)

quant.T.75 <- rq(T.br ~ poly(r.max_T, 2), data = df, tau = 0.75) 
pred.t$T.br.75 <- predict(quant.T.75, newdata = pred.t)

quant.T.50 <- rq(T.br ~ poly(r.max_T, 2), data = df, tau = 0.50) 
pred.t$T.br.50 <- predict(quant.T.50, newdata = pred.t)

T.qrs <- ggplot(df, aes(x = r.max_T, y = T.br, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.t, aes(x = r.max_T, y = T.br.100), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.t, aes(x = r.max_T, y = T.br.90), color = "black", size = 1.2, linetype = "dashed") +  
  geom_line(data = pred.t, aes(x = r.max_T, y = T.br.75), color = "black", size = 1.2, linetype = "dotted") +  
  geom_line(data = pred.t, aes(x = r.max_T, y = T.br.50), color = "black", size = 1.2, linetype = "dotdash") +  
  
  geom_text_repel(data = pred.t[which.min(abs(df$r.max_T - median(df$r.max_T))), ],  
                  aes(x = r.max_T, y = T.br.100, label = "100"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.t[which.min(abs(df$r.max_T - median(df$r.max_T))), ], 
                  aes(x = r.max_T, y = T.br.90, label = "90"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.t[which.min(abs(df$r.max_T - median(df$r.max_T))), ], 
                  aes(x = r.max_T, y = T.br.75, label = "75"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.t[which.min(abs(df$r.max_T - median(df$r.max_T))), ], 
                  aes(x = r.max_T, y = T.br.50, label = "50"), 
                  color = "black", size = 5, fontface = "bold") + # Adding labels where lines intersect data points
  
  labs(x = "Maximum exponential growth rate",    
       y = "Thermal breadth", 
       title = "Quantile Regressions of Temperature Performance Trade-offs",
       color = "Evolutionary History") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1"), 
                     labels = c("Ancestral", "Evolved")) +  
  theme_classic() +
  theme(
    legend.position = c(0.8, 0.85),  # Move legend inside the plot
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10, face = "bold"),  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10) # theme stuff
  )

T.qrs  # Display the plot

############# Light ##########################################

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "L", 8, 16)) # Ls are now equivalent to 8, for later mapping

df$evol.bin <- ifelse(df$evol == "none", 'ancestral', 
                      ifelse(df$evol == "L", 'light', 'other')) # for testing regressions.

par.res.L <- par_frt(df[df$I.comp<10,], xvar = "r.max_I", yvar = "I.comp")

L_par <- ggplot(df[df$I.comp<10,], aes(x = r.max_I, y = I.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Light limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res.L, aes(x = r.max_I, y = I.comp), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

L_par # Raw pareto front.

mod.L <- lm(I.comp~r.max_I*evol.bin, data=df[df$I.comp<10,])  # test significance of Pareto frontier
summary(mod.L)

L.regs <- ggplot(df[df$I.comp<10,], aes(x = r.max_I, y = I.comp, color = evol.bin)) +  
  geom_point(size = 2) +  # Scatter plot of raw data
  geom_smooth(method = "lm", se = FALSE, aes(color = evol.bin)) +  # Separate regression lines
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Evolutionary Effects on Light-Competition Trade-offs",
       color = "Evolutionary History") +  # Change legend title) +
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Light")) +  # Custom colors per group
  theme_classic() +
  theme(
    legend.position = c(0.85, 0.85),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 10, face = "bold"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10)
  ) 

L.regs

pred.l <- data.frame(r.max_I = seq(min(df$r.max_I), max(df$r.max_I), length.out = 100)) # Dataframe to collect quantile info in

quant.I.100 <- rq(I.comp ~ poly(r.max_I, 2), data = df[df$I.comp<10,], tau = 1.00) 
pred.l$I.comp.100 <- predict(quant.I.100, newdata = pred.l)

quant.I.95 <- rq(I.comp ~ poly(r.max_I, 2), data = df[df$I.comp<10,], tau = 0.95) 
pred.l$I.comp.95 <- predict(quant.I.95, newdata = pred.l)

quant.I.90 <- rq(I.comp ~ poly(r.max_I, 2), data = df[df$I.comp<10,], tau = 0.90) 
pred.l$I.comp.90 <- predict(quant.I.90, newdata = pred.l)

quant.I.75 <- rq(I.comp ~ poly(r.max_I, 2), data = df[df$I.comp<10,], tau = 0.75) 
pred.l$I.comp.75 <- predict(quant.I.75, newdata = pred.l)

quant.I.50 <- rq(I.comp ~ poly(r.max_I, 2), data = df[df$I.comp<10,], tau = 0.50) 
pred.l$I.comp.50 <- predict(quant.I.50, newdata = pred.l)

L.qrs <- ggplot(df[df$I.comp<10,], aes(x = r.max_I, y = I.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.l, aes(x = r.max_I, y = I.comp.100), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.l, aes(x = r.max_I, y = I.comp.90), color = "black", size = 1.2, linetype = "dashed") +  
  geom_line(data = pred.l, aes(x = r.max_I, y = I.comp.75), color = "black", size = 1.2, linetype = "dotted") +  
  geom_line(data = pred.l, aes(x = r.max_I, y = I.comp.50), color = "black", size = 1.2, linetype = "dotdash") +  
  
  geom_text_repel(data = pred.l[which.min(abs(df$r.max_I - median(df$r.max_I))), ],  
                  aes(x = r.max_I, y = I.comp.100, label = "100"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.l[which.min(abs(df$r.max_I - median(df$r.max_I))), ], 
                  aes(x = r.max_I, y = I.comp.90, label = "90"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.l[which.min(abs(df$r.max_I - median(df$r.max_I))), ], 
                  aes(x = r.max_I, y = I.comp.75, label = "75"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.l[which.min(abs(df$r.max_I - median(df$r.max_I))), ], 
                  aes(x = r.max_I, y = I.comp.50, label = "50"), 
                  color = "black", size = 5, fontface = "bold") + # Adding labels where lines intersect data points
  
  labs(x = "Maximum exponential growth rate",    
       y = "Competitive ability (1/R*)", 
       title = "Quantile Regressions of Light-Competition Trade-offs",
       color = "Evolutionary History") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Light")) +  
  theme_classic() +
  theme(
    legend.position = c(0.85, 0.85),  # Move legend inside the plot
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10, face = "bold"),  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10) # theme stuff
  )

L.qrs  # Display the plot

############# Nitrogen #######################################

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "N", 8, 16)) # Ns are now equivalent to 8, for later mapping

df$evol.bin <- ifelse(df$evol == "none", 'ancestral', 
                      ifelse(df$evol == "N", 'nit', 'other')) # for testing regressions.

par.res.N <- par_frt(df, xvar = "r.max_N", yvar = "N.comp")

N_par <- ggplot(df, aes(x = r.max_N, y = N.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Nitrogen limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res.N, aes(x = r.max_N, y = N.comp), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

N_par # Raw pareto front.

mod.N <- lm(N.comp~r.max_N*evol.bin, data=df)  # test significance of Pareto frontier
summary(mod.N)

N.regs <- ggplot(df, aes(x = r.max_N, y = N.comp, color = evol.bin)) +  
  geom_point(size = 2) +  # Scatter plot of raw data
  geom_smooth(method = "lm", se = FALSE, aes(color = evol.bin)) +  # Separate regression lines
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Evolutionary Effects on N-Competition Trade-offs",
       color = "Evolutionary History") +  # Change legend title) +
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Nitrogen")) +  # Custom colors per group
  theme_classic() +
  theme(
    legend.position = c(0.35, 0.85),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 10, face = "bold"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10)
  ) 

N.regs

pred.n <- data.frame(r.max_N = seq(min(df$r.max_N), max(df$r.max_N), length.out = 100)) # Dataframe to collect quantile info in

quant.N.100 <- rq(N.comp ~ poly(r.max_N, 2), data = df, tau = 1.00) 
pred.n$N.comp.100 <- predict(quant.N.100, newdata = pred.n)

quant.N.95 <- rq(N.comp ~ poly(r.max_N, 2), data = df, tau = 0.95) 
pred.n$N.comp.95 <- predict(quant.N.95, newdata = pred.n)

quant.N.90 <- rq(N.comp ~ poly(r.max_N, 2), data = df, tau = 0.90) 
pred.n$N.comp.90 <- predict(quant.N.90, newdata = pred.n)

quant.N.75 <- rq(N.comp ~ poly(r.max_N, 2), data = df, tau = 0.75) 
pred.n$N.comp.75 <- predict(quant.N.75, newdata = pred.n)

quant.N.50 <- rq(N.comp ~ poly(r.max_N, 2), data = df, tau = 0.50) 
pred.n$N.comp.50 <- predict(quant.N.50, newdata = pred.n)


N.qrs <- ggplot(df, aes(x = r.max_N, y = N.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.n, aes(x = r.max_N, y = N.comp.100), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.n, aes(x = r.max_N, y = N.comp.90), color = "black", size = 1.2, linetype = "dashed") +  
  geom_line(data = pred.n, aes(x = r.max_N, y = N.comp.75), color = "black", size = 1.2, linetype = "dotted") +  
  geom_line(data = pred.n, aes(x = r.max_N, y = N.comp.50), color = "black", size = 1.2, linetype = "dotdash") +  
  
  geom_text_repel(data = pred.n[which.min(abs(df$r.max_N - median(df$r.max_N))), ],  
                  aes(x = r.max_N, y = N.comp.100, label = "100"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.n[which.min(abs(df$r.max_N - median(df$r.max_N))), ], 
                  aes(x = r.max_N, y = N.comp.90, label = "90"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.n[which.min(abs(df$r.max_N - median(df$r.max_N))), ], 
                  aes(x = r.max_N, y = N.comp.75, label = "75"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.n[which.min(abs(df$r.max_N - median(df$r.max_N))), ], 
                  aes(x = r.max_N, y = N.comp.50, label = "50"), 
                  color = "black", size = 5, fontface = "bold") + # Adding labels where lines intersect data points
  
  labs(x = "Maximum exponential growth rate",    
       y = "Competitive ability (1/R*)", 
       title = "Quantile Regressions of N-Competition Trade-offs",
       color = "Evolutionary History") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Nitrogen")) +  
  theme_classic() +
  theme(
    legend.position = c(0.85, 0.85),  # Move legend inside the plot
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10, face = "bold"),  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10) # theme stuff
  )

N.qrs  # Display the plot

############# Phosphorous ####################################

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "P", 8, 16)) # Ps are now equivalent to 8, for later mapping

df$evol.bin <- ifelse(df$evol == "none", 'ancestral', 
                      ifelse(df$evol == "P", 'phos', 'other')) # for testing regressions.

par.res.P <- par_frt(df, xvar = "r.max_P", yvar = "P.comp")

P_par <- ggplot(df, aes(x = r.max_P, y = P.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res.P, aes(x = r.max_P, y = P.comp), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_par # Raw pareto front.

mod.P <- lm(P.comp~r.max_P*evol.bin, data=df)  # test significance of Pareto frontier
summary(mod.P)

P.regs <- ggplot(df, aes(x = r.max_P, y = P.comp, color = evol.bin)) +  
  geom_point(size = 2) +  # Scatter plot of raw data
  geom_smooth(method = "lm", se = FALSE, aes(color = evol.bin)) +  # Separate regression lines
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Evolutionary Effects on P-Competition Trade-offs",
       color = "Evolutionary History") +  # Change legend title) +
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Phosphorous")) +  # Custom colors per group
  theme_classic() +
  theme(
    legend.position = c(0.85, 0.85),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 10, face = "bold"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10)
  ) 

P.regs

pred.p <- data.frame(r.max_P = seq(min(df$r.max_P), max(df$r.max_P), length.out = 100)) # Dataframe to collect quantile info in

quant.P.100 <- rq(P.comp ~ poly(r.max_P, 2), data = df, tau = 1.00) 
pred.p$P.comp.100 <- predict(quant.P.100, newdata = pred.p)

quant.P.95 <- rq(P.comp ~ poly(r.max_P, 2), data = df, tau = 0.95) 
pred.p$P.comp.95 <- predict(quant.P.95, newdata = pred.p)

quant.P.90 <- rq(P.comp ~ poly(r.max_P, 2), data = df, tau = 0.90) 
pred.p$P.comp.90 <- predict(quant.P.90, newdata = pred.p)

quant.P.75 <- rq(P.comp ~ poly(r.max_P, 2), data = df, tau = 0.75) 
pred.p$P.comp.75 <- predict(quant.P.75, newdata = pred.p)

quant.P.50 <- rq(P.comp ~ poly(r.max_P, 2), data = df, tau = 0.50) 
pred.p$P.comp.50 <- predict(quant.P.50, newdata = pred.p)


P.qrs <- ggplot(df, aes(x = r.max_P, y = P.comp, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 2) +  # Scatter plot of raw data

  geom_line(data = pred.p, aes(x = r.max_P, y = P.comp.100), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.p, aes(x = r.max_P, y = P.comp.90), color = "black", size = 1.2, linetype = "dashed") +  
  geom_line(data = pred.p, aes(x = r.max_P, y = P.comp.75), color = "black", size = 1.2, linetype = "dotted") +  
  geom_line(data = pred.p, aes(x = r.max_P, y = P.comp.50), color = "black", size = 1.2, linetype = "dotdash") +  
  
  geom_text_repel(data = pred.p[which.min(abs(df$r.max_P - median(df$r.max_P))), ],  
                  aes(x = r.max_P, y = P.comp.100, label = "100"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.p[which.min(abs(df$r.max_P - median(df$r.max_P))), ], 
                  aes(x = r.max_P, y = P.comp.90, label = "90"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.p[which.min(abs(df$r.max_P - median(df$r.max_P))), ], 
                  aes(x = r.max_P, y = P.comp.75, label = "75"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.p[which.min(abs(df$r.max_P - median(df$r.max_P))), ], 
                  aes(x = r.max_P, y = P.comp.50, label = "50"), 
                  color = "black", size = 5, fontface = "bold") + # Adding labels where lines intersect data points
  
  labs(x = "Maximum exponential growth rate",    
       y = "Competitive ability (1/R*)", 
       title = "Quantile Regressions of P-Competition Trade-offs",
       color = "Evolutionary History") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Phosphorous")) +  
  theme_classic() +
  theme(
    legend.position = c(0.85, 0.85),  # Move legend inside the plot
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10, face = "bold"),  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10) # theme stuff
  )

P.qrs  # Display the plot

############# Salt ###########################################

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol %in% c("S", "BS"), 8, 16)) # Ss and BSs are now equivalent to 8, for later mapping

df$evol.bin <- ifelse(df$evol == "none", 'ancestral', 
                      ifelse(df$evol %in% c("S", "BS"), 'salt', 'other')) # for testing regressions.

par.res.S <- par_frt(df, xvar = "r.max_S", yvar = "S.c.mod")

S_par <- ggplot(df, aes(x = r.max_S, y = S.c.mod)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Salt tolerance (c)", 
       title = "Salt stress") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res.S, aes(x = r.max_S, y = S.c.mod), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

S_par # Raw pareto front.

mod.S <- lm(S.c.mod~r.max_S*evol.bin, data=df)  # test significance of Pareto frontier
summary(mod.S)

S.regs <- ggplot(df, aes(x = r.max_S, y = S.c.mod, color = evol.bin)) +  
  geom_point(size = 2) +  # Scatter plot of raw data
  geom_smooth(method = "lm", se = FALSE, aes(color = evol.bin)) +  # Separate regression lines
  labs(x = "Maximum exponential growth rate", 
       y = "Salt tolerance (c)", 
       title = "Evolutionary Effects on S-Tolerance Trade-offs",
       color = "Evolutionary History") +  # Change legend title) +
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Salt or Biotic Salt")) +  # Custom colors per group
  theme_classic() +
  theme(
    legend.position = c(0.75, 0.85),  # Moves legend inside the plot (x, y) in [0,1] scale
    legend.title = element_text(size = 12, face = "bold"),  # Adjust title size
    legend.text = element_text(size = 10, face = "bold"),  # Adjust text size
    axis.title = element_text(size = 12, face = "bold"),  # Bold & larger axis titles
    axis.text = element_text(size = 10)
  ) 

S.regs

pred.S <- data.frame(r.max_S = seq(min(df$r.max_S), max(df$r.max_S), length.out = 100)) # Dataframe to collect quantile info in

quant.S.100 <- rq(S.c.mod ~ poly(r.max_S, 2), data = df, tau = 1.00) 
pred.S$S.c.mod.100 <- predict(quant.S.100, newdata = pred.S)

quant.S.95 <- rq(S.c.mod ~ poly(r.max_S, 2), data = df, tau = 0.95) 
pred.S$S.c.mod.95 <- predict(quant.S.95, newdata = pred.S)

quant.S.90 <- rq(S.c.mod ~ poly(r.max_S, 2), data = df, tau = 0.90) 
pred.S$S.c.mod.90 <- predict(quant.S.90, newdata = pred.S)

quant.S.75 <- rq(S.c.mod ~ poly(r.max_S, 2), data = df, tau = 0.75) 
pred.S$S.c.mod.75 <- predict(quant.S.75, newdata = pred.S)

quant.S.50 <- rq(S.c.mod ~ poly(r.max_S, 2), data = df, tau = 0.50) 
pred.S$S.c.mod.50 <- predict(quant.S.50, newdata = pred.S)

S.qrs <- ggplot(df, aes(x = r.max_S, y = S.c.mod, color = evol.bin)) +  # Quantiles plot
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.S, aes(x = r.max_S, y = S.c.mod.100), color = "black", size = 1.2) +  # Adding all quantile regression lines as black lines
  geom_line(data = pred.S, aes(x = r.max_S, y = S.c.mod.90), color = "black", size = 1.2, linetype = "dashed") +  
  geom_line(data = pred.S, aes(x = r.max_S, y = S.c.mod.75), color = "black", size = 1.2, linetype = "dotted") +  
  geom_line(data = pred.S, aes(x = r.max_S, y = S.c.mod.50), color = "black", size = 1.2, linetype = "dotdash") +  
  
  geom_text_repel(data = pred.S[which.min(abs(df$r.max_S - median(df$r.max_S))), ],  
                  aes(x = r.max_S, y = S.c.mod.100, label = "100"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.S[which.min(abs(df$r.max_S - median(df$r.max_S))), ], 
                  aes(x = r.max_S, y = S.c.mod.90, label = "90"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.S[which.min(abs(df$r.max_S - median(df$r.max_S))), ], 
                  aes(x = r.max_S, y = S.c.mod.75, label = "75"), 
                  color = "black", size = 5, fontface = "bold") +
  geom_text_repel(data = pred.S[which.min(abs(df$r.max_S - median(df$r.max_S))), ], 
                  aes(x = r.max_S, y = S.c.mod.50, label = "50"), 
                  color = "black", size = 5, fontface = "bold") + # Adding labels where lines intersect data points
  
  labs(x = "Maximum exponential growth rate",    
       y = "Salt tolerance (c)", 
       title = "Quantile Regressions of S-Tolerance Trade-offs",
       color = "Evolutionary History") +  # labels
  
  scale_color_manual(values = c("black", "goldenrod1", "mediumorchid3"), 
                     labels = c("Ancestral", "Other", "Salt or Biotic Salt")) +  
  theme_classic() +
  theme(
    legend.position = c(0.8, 0.85),  # Move legend inside the plot
    legend.title = element_text(size = 12, face = "bold"),  
    legend.text = element_text(size = 10, face = "bold"),  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10) # theme stuff
  )

S.qrs  # Display the plot

plot_grid(L_par, N_par, P_par, S_par, T_par) # Summary plots

plot_grid(L.qrs, N.qrs, P.qrs, S.qrs, T.qrs)

L.qrs1 <- L.qrs + labs(x = "Maximum exponential growth rate",    
                       y = "Light competition (1/R*)", 
                       title = "Light") + theme(legend.position = 'none') # for summary plotting

N.qrs1 <- N.qrs + labs(x = "Maximum exponential growth rate",    
                       y = "Nitrogen competition (1/R*)", 
                       title = "Nitrogen") + theme(legend.position = 'none') # for summary plotting

P.qrs1 <- P.qrs + labs(x = "Maximum exponential growth rate",    
                       y = "Phosphorous competition (1/R*)", 
                       title = "Phosphorous") + theme(legend.position = 'none') # for summary plotting

S.qrs1 <- S.qrs + labs(x = "Maximum exponential growth rate",    
                       y = "Salt tolerance (c)", 
                       title = "Salt") + theme(legend.position = 'none') # for summary plotting

T.qrs1 <- T.qrs + labs(x = "Maximum exponential growth rate",    
                       y = "Thermal breadth", 
                       title = "Temperature") + theme(legend.position = 'none') # for summary plotting
      
plot_grid(L.qrs1, N.qrs1, P.qrs1, S.qrs1, T.qrs1)

####################### Multivariate analysis ################################

# PCA

df$evol.plt <- factor(df$evol, 
                      levels = c("none", "L", "N", "P", "S", "B", "BS", "C"),
                      labels = c("Ancestral", 
                                 "Light limitation", 
                                 "Nitrogen limitation", 
                                 "Phosphorous limitation", 
                                 "Salt stress", 
                                 "Biotic depletion", 
                                 "Biotic depletion x Salt", 
                                 "Control")) # For plotting

df.pca <- df %>% select(T.br, I.comp, N.comp, P.comp, S.c.mod, evol.plt) # Prepare the data: selecting only the relevant columns
df.pca

df.pca <- df.pca[df.pca$I.comp<10,] # For now removing wonky/missing points

evol.fil <- df.pca$evol.plt

df.pca <- df.pca %>% select(-evol.plt)

pca.result <- prcomp(df.pca, center = TRUE, scale. = TRUE) # Perform PCA
pca.result

pca.df <- data.frame(pca.result$x, evol.fil)  # Add grouping factor
colnames(pca.df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "Evolution")
pca.df$Evolution <- factor(pca.df$Evolution)


screeplot(pca.result, type = "lines") # Scree plot to check explained variance

PCA <- ggplot(pca.df, aes(x = PC1, y = PC2, color = Evolution)) +  # PCA biplot visualization
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("PC1 (", round(pca.result$sdev[1]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pca.result$sdev[2]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       title = "PCA of Thermal Performance and Competitive Abilities") +
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
  )

PCA

loadings <- as.data.frame(pca.result$rotation) # Extract PCA loadings (rotation matrix)

loadings$PC1 <- loadings$PC1 * max(abs(pca.df$PC1)) # Scale loadings to fit within the PCA plot
loadings$PC2 <- loadings$PC2 * max(abs(pca.df$PC2))

loadings$variable <- rownames(loadings) # names

# Create PCA plot with arrows
pca_plot_arrows <- ggplot(pca.df, aes(x = PC1, y = PC2, color = Evolution)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("PC1 (", round(pca.result$sdev[1]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pca.result$sdev[2]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       title = "PCA of Thermal Performance and Competitive Abilities") +
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
  # Add arrows for variable contributions
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  # Add variable names to the plot
  geom_text(data = loadings, aes(x = PC1, y = PC2, label = variable),
            vjust = 1, hjust = 1, color = "black", size = 5)

pca_plot_arrows

df.pca2 <- df %>% select(T.br, r.max_T, I.comp, r.max_I, N.comp, r.max_N, P.comp, r.max_P, S.c.mod, r.max_S, evol.plt) # Prepare the data: selecting only the relevant columns
df.pca2

df.pca2 <- df.pca2[df.pca2$I.comp<10,] # For now removing wonky/missing points

df.pca2 <- df.pca2 %>% select(-evol.plt)

pca.result2 <- prcomp(df.pca2, center = TRUE, scale. = TRUE) # Perform PCA
pca.result2

pca.df2 <- data.frame(pca.result2$x, evol.fil)  # Add grouping factor
colnames(pca.df2) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Evolution")
pca.df2$Evolution <- factor(pca.df2$Evolution)

PCA2 <- ggplot(pca.df2, aes(x = PC1, y = PC2, color = Evolution)) +  # PCA biplot visualization
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("PC1 (", round(pca.result2$sdev[1]^2 / sum(pca.result2$sdev^2) * 100, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pca.result2$sdev[2]^2 / sum(pca.result2$sdev^2) * 100, 2), "%)", sep = ""),
       title = "PCA of Thermal Performance and Nutrient Utilization") +
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
  )

PCA2

loadings2 <- as.data.frame(pca.result2$rotation) # loadings

loadings2$PC1 <- loadings2$PC1 * max(abs(pca.df2$PC1)) # Scale loadings to fit within the PCA plot
loadings2$PC2 <- loadings2$PC2 * max(abs(pca.df2$PC2))

loadings2$variable <- rownames(loadings2)

pca_plot_arrows2 <- ggplot(pca.df2, aes(x = PC1, y = PC2, color = Evolution)) + # Create PCA plot with arrows
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("PC1 (", round(pca.result2$sdev[1]^2 / sum(pca.result2$sdev^2) * 100, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pca.result2$sdev[2]^2 / sum(pca.result2$sdev^2) * 100, 2), "%)", sep = ""),
       title = "PCA of Thermal Performance and Nutrient Utilization") +
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
  # Add arrows for variable contributions
  geom_segment(data = loadings2, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  # Add variable names to the plot
  geom_text(data = loadings2, aes(x = PC1, y = PC2, label = variable),
            vjust = 1, hjust = 1, color = "black", size = 5) 

pca_plot_arrows2

# Redundancy analysis

response_vars <- df.pca
explanatory_vars <- model.matrix(~ evol.fil)[, -1]  # Remove intercept

rda_result <- rda(response_vars ~ ., data = as.data.frame(explanatory_vars)) # run the RDA
summary(rda_result)

rda_sites <- as.data.frame(scores(rda_result, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species <- as.data.frame(scores(rda_result, display = "species")) # Extract RDA species (trait arrows)

rda_constraints <- as.data.frame(scores(rda_result, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites$Evolution <- evol.fil # Add the evolutionary treatment labels to the site scores

# Prepare RDA site data (samples)
rda_sites <- rda_sites %>% 
  rename(RDA1 = RDA1, RDA2 = RDA2)

# Prepare RDA trait data (arrows)
rda_species <- rda_species %>%
  rename(RDA1 = RDA1, RDA2 = RDA2)

# Prepare centroids of treatment conditions
rda_constraints <- rda_constraints %>%
  rename(RDA1 = RDA1, RDA2 = RDA2)

# Assign readable labels for centroids
rda_constraints$label <- rownames(rda_constraints)

# Create the RDA biplot
rda_plot <- ggplot() +
  geom_point(data = rda_sites, aes(x = RDA1, y = RDA2, color = Evolution), size = 3, alpha = 0.8) +
  
  # Retain trait arrows
  geom_segment(data = rda_species, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  
  # Retain trait labels only
  geom_text_repel(data = rda_species, aes(x = RDA1, y = RDA2, label = rownames(rda_species)), 
                  size = 5, color = "black") +
  
  # Color legend to match PCA
  scale_color_manual(
    name = "Evolution environment",  # Update the legend title
    values = c(
      "Biotic depletion" = "darkorange",
      "Biotic depletion x Salt" = "deepskyblue1",
      "Control" = "forestgreen",
      "Light limitation" = "gold",
      "Nitrogen limitation" = "magenta3",
      "Ancestral" = "black",
      "Phosphorous limitation" = "firebrick",  
      "Salt stress" = "blue"
    )
  ) +
  
  # Axes labels
  labs(
    x = paste("RDA1 (", round(summary(rda_result)$cont$importance[2, 1] * 100, 2), "%)", sep = ""),
    y = paste("RDA2 (", round(summary(rda_result)$cont$importance[2, 2] * 100, 2), "%)", sep = ""),
    title = "Redundancy Analysis (RDA) of Traits by Evolution Environment"
  ) +
  
  theme_classic() +
  theme(
    text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Display the plot
rda_plot

#OK now working with ancestry as the explanatory variable?

df.pca1.1 <- df %>% select(T.br, I.comp, N.comp, P.comp, S.c.mod, anc) # Prepare the data: selecting only the relevant columns
df.pca1.1

df.pca1.1 <- df.pca1.1[df.pca1.1$I.comp<10,] # For now removing wonky/missing points

anc.fil <- df.pca1.1$anc

df.pca1.1 <- df.pca1.1 %>% select(-anc)

response_vars1.1 <- df.pca1.1
explanatory_vars1.1 <- model.matrix(~ anc.fil)[, -1]  # Remove intercept

rda_result1.1 <- rda(response_vars1.1 ~ ., data = as.data.frame(explanatory_vars1.1)) # run the RDA
summary(rda_result1.1)

rda_sites1.1 <- as.data.frame(scores(rda_result1.1, display = "sites")) # Extract RDA site scores (sample coordinates)

rda_species1.1 <- as.data.frame(scores(rda_result1.1, display = "species")) # Extract RDA species (trait arrows)

rda_constraints1.1 <- as.data.frame(scores(rda_result1.1, display = "bp")) # Extract explanatory variable centroids (e.g., treatment centroids)

rda_sites1.1$anc <- anc.fil # Add the evolutionary treatment labels to the site scores

# Prepare RDA site data (samples)
rda_sites1.1 <- rda_sites1.1 %>% 
  rename(RDA1 = RDA1, RDA2 = RDA2)

# Prepare RDA trait data (arrows)
rda_species1.1 <- rda_species1.1 %>%
  rename(RDA1 = RDA1, RDA2 = RDA2)

# Prepare centroids of treatment conditions
rda_constraints1.1 <- rda_constraints1.1 %>%
  rename(RDA1 = RDA1, RDA2 = RDA2)

# Assign readable labels for centroids
rda_constraints1.1$label <- rownames(rda_constraints1.1)

rda_plot1.1 <- ggplot() +
  geom_point(data = rda_sites1.1, aes(x = RDA1, y = RDA2, color = anc), size = 3, alpha = 0.8) +
  
  # Retain trait arrows
  geom_segment(data = rda_species1.1, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.25, "cm")), color = "black") +
  
  # Retain trait labels only
  geom_text_repel(data = rda_species1.1, aes(x = RDA1, y = RDA2, label = rownames(rda_species1.1)), 
                  size = 5, color = "black") +
  
  # Color legend to match PCA
  scale_color_manual(
    name = "Ancestry",  # Update the legend title
    values = c(
      "anc2" = "darkorange",
      "anc3" = "gold",
      "anc4" = "deepskyblue1",
      "anc5" = "forestgreen",
      "cc1690" = "magenta3"
    )
  ) +
  
  # Axes labels
  labs(
    x = paste("RDA1 (", round(summary(rda_result1.1)$cont$importance[2, 1] * 100, 2), "%)", sep = ""),
    y = paste("RDA2 (", round(summary(rda_result1.1)$cont$importance[2, 2] * 100, 2), "%)", sep = ""),
    title = "Redundancy Analysis (RDA) of Traits by Ancestry"
  ) +
  
  theme_classic() +
  theme(
    text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Display the plot
rda_plot1.1

plot_grid(pca_plot_arrows, rda_plot)
pca_plot_arrows
