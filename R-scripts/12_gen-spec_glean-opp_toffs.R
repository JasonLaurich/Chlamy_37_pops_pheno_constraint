# Jason R Laurich
# March 1, 2025

# Moving the analysis portion of script 11 here. It was getting unwieldy
# Modelling gleaner-opportunist, specialist-generalist trade-offs
# Adding Pareto frontier estimation to correlation plots. 

# Conceptual approach - (1) use Convex Hull algorithm to identify potential Pareto frontier.
# (2) Model regression of full data set, evolved populations, ancestral using brms - are the slopes of the evolved populations different? Does the intercept differ?
# (3) Model regression using non-linear regression - same question - differences in shape (especially where evolved populations form the PF) would be strong evidence.

# So partition out ancestors, evolutionary conditions particular to the nutrient in question, and other evolutionary treatments?

# Another angle for confirming the Pareto Frontier would be bootstrapping the data, and asking whether, on average, the Pareto front is more likely
# to contain points that evolved under relevant conditions?

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

# Let's define a Pareto front function
par_frt <- function(df, xvar, yvar) { # general pareto frontier function
  
  df <- df[order(-df[[xvar]], df[[yvar]]), ]  
  pareto_points <- df[1, ]  # Start with the first point
  
  for (i in 2:nrow(df)) {
    if (df[i, yvar] > tail(pareto_points[[yvar]], 1)) {  # Ensure increasing y values
      pareto_points <- rbind(pareto_points, df[i,])
    }
  }
  
  return(pareto_points)
}

par.res <- par_frt(df, xvar = "r.max_T", yvar = "T.br")

T_par <- ggplot(df, aes(x = r.max_T, y = T.br)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Thermal breadth", 
       title = "Thermal performance") +
  scale_shape_manual(values = c(16, 3)) +  # Keep custom shapes
  geom_line(data = par.res, aes(x = r.max_T, y = T.br), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

T_par # Raw pareto front.

mod_qr <- rq(T.br ~ poly(r.max_T, 2), tau = 0.95, data = df) # Frequentist quantile regression for plotting Pareto frontiers. 

new.data <- data.frame(r.max_T = seq(min(df$r.max_T), max(df$r.max_T), length.out = 100))
new.data$pred_br <- predict(mod_qr, new.data)

T_par_qr <- ggplot(df, aes(x = r.max_T, y = T.br)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  geom_line(data = new.data, aes(x = r.max_T, y = pred_br), color = "blue", size = 1) + # quantile regression in blue
  labs(x = "Maximum exponential growth rate", 
       y = "Thermal breadth", 
       title = "Thermal performance") +
  scale_shape_manual(values = c(16, 3)) +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

T_par_qr # Quantile regression pareto front.

brm_sp <- brm(
  bf(T.br ~ s(r.max_T), quantile = 0.95),
  data = par.res,  # Use only the Pareto front subset
  family = asym_laplace(),
  prior = prior(sigma ~ normal(0, 5)),  # Optional priors
  iter = 4000, chains = 4, cores = 4
)

brm_qd <- brm(
  bf(T.br ~ poly(r.max_T, 2), quantile = 0.95),  # Quadratic term instead of spline
  data = par.res,  
  family = asym_laplace(),
  iter = 4000, chains = 4, cores = 4
) # Quadratic Bayesian function

new.data.bqr <- data.frame(r.max_T = seq(min(df$r.max_T), max(df$r.max_T), length.out = 100))

preds <- fitted(brm_sp, newdata = new.data.bqr, probs = c(0.05, 0.95))  

new.data.bqr$pred_br <- preds[, "Estimate"] # Add predictions to new_data
new.data.bqr$lower <- preds[, "Q5"]   # 5% CI
new.data.bqr$upper <- preds[, "Q95"]  # 95% CI

T_par_bqr <- ggplot(df, aes(x = r.max_T, y = T.br)) +
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  geom_ribbon(data = new.data.bqr, aes(x = r.max_T, ymin = lower, ymax = upper), 
              fill = "grey75", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = new.data.bqr, aes(x = r.max_T, y = pred_br), color = "blue", size = 1) + # Bayesian quantile regression
  labs(x = "Maximum exponential growth rate", 
       y = "Thermal breadth", 
       title = "Thermal performance with Bayesian quantile regression (95% HDPI)") +
  scale_shape_manual(values = c(16, 3)) +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

T_par_bqr  # Display the plot

############# Light ##########################################

############# Nitrogen #######################################

############# Phosphorous ####################################

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "P", 8, 16)) # Ps are now equivalent to 8, for later mapping

par.res <- par_frt(df, xvar = "r.max_P", yvar = "P.comp")

P_par <- ggplot(df, aes(x = r.max_P, y = P.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res, aes(x = r.max_P, y = P.comp), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_par # Raw pareto front.

brm_qd_P <- brm(
  bf(P.comp ~ poly(r.max_P, 2), quantile = 0.95),  # Quadratic term instead of spline
  data = par.res,  
  family = asym_laplace(),
  iter = 4000, chains = 4, cores = 4
)

brm_sp_P <- brm(
  bf(P.comp ~ s(r.max_P, bs = "gp"), quantile = 0.95),  # Use a Gaussian Process spline
  data = par.res,  
  family = asym_laplace(),
  prior = c(prior(normal(0, 5), class = "b")), # Regularization to prevent overfitting
  iter = 4000, chains = 4, cores = 4
)

brm_sp_P2 <- brm(
  bf(P.comp ~ s(r.max_P, k = 5, bs = "tp"), quantile = 0.95),  # Thin-plate regression spline
  data = par.res,  
  family = asym_laplace(),
  iter = 4000, chains = 4, cores = 4
)

new_data <- data.frame(r.max_P = seq(min(par.res$r.max_P), max(par.res$r.max_P), length.out = 100))
preds <- as.data.frame(fitted(brm_sp_P2, newdata = new_data, probs = c(0.05, 0.95))) 

# Add predictions to new_data
new_data$pred_P.comp <- preds[, "Estimate"]
new_data$lower <- preds[, "Q5"]   # 5% CI
new_data$upper <- preds[, "Q95"]  # 95% CI

preds$r.max_P <- new_data$r.max_P

P_par_bsp <- ggplot(df, aes(x = r.max_P, y = P.comp)) +  
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  geom_ribbon(data = new_data, aes(x = r.max_P, ymin = lower, ymax = upper), 
              fill = "grey75", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = preds, aes(x = r.max_P, y = Estimate), color = "blue", size = 1) +  # Bayesian quantile regression in blue
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation, Bayesian quantile regression (95% HDPI)") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_par_bsp # Display the plot

P_par_par <- ggplot(par.res, aes(x = r.max_P, y = P.comp)) +  
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  geom_ribbon(data = new_data, aes(x = r.max_P, ymin = lower, ymax = upper), 
              fill = "grey75", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = preds, aes(x = r.max_P, y = Estimate), color = "blue", size = 1) +  # Bayesian quantile regression in blue
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation, Bayesian quantile regression (95% HDPI)") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_par_par # Just the pareto points

brm_sp_P2.1 <- brm(
  bf(P.comp ~ s(r.max_P, k = 5, bs = "tp"), quantile = 0.95),  # Thin-plate regression spline
  data = df,  
  family = asym_laplace(),
  iter = 4000, chains = 4, cores = 4
)

new_data <- data.frame(r.max_P = seq(min(par.res$r.max_P), max(par.res$r.max_P), length.out = 100))
preds <- as.data.frame(fitted(brm_sp_P2.1, newdata = new_data, probs = c(0.05, 0.95))) 

# Add predictions to new_data
new_data$pred_P.comp <- preds[, "Estimate"]
new_data$lower <- preds[, "Q5"]   # 5% CI
new_data$upper <- preds[, "Q95"]  # 95% CI

preds$r.max_P <- new_data$r.max_P

P_par_bq <- ggplot(df, aes(x = r.max_P, y = P.comp)) +  
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  geom_ribbon(data = new_data, aes(x = r.max_P, ymin = lower, ymax = upper), 
              fill = "grey75", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = preds, aes(x = r.max_P, y = Estimate), color = "blue", size = 1) +  # Bayesian quantile regression in blue
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation, Bayesian quantile regression (95% HDPI)") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_par_bq # Display the plot

#OK we shouldn't need the quantile regression, just the spline. The data is already subsetted to the Pareto front.
brm_sp_P3 <- brm(
  bf(P.comp ~ s(r.max_P, k = 3, bs = "tp")), # thin-plate spline
  data = par.res,  
  family = gaussian(),  # Regular Gaussian likelihood instead of asymmetric Laplace
  iter = 4000, chains = 4, cores = 4
)

new_data <- data.frame(r.max_P = seq(min(par.res$r.max_P), max(par.res$r.max_P), length.out = 100))
preds <- as.data.frame(fitted(brm_sp_P3, newdata = new_data, probs = c(0.025, 0.975))) 

# Add predictions to new_data
new_data$pred_P.comp <- preds[, "Estimate"]
new_data$lower <- preds[, "Q2.5"]   # 5% CI
new_data$upper <- preds[, "Q97.5"]  # 95% CI

preds$r.max_P <- new_data$r.max_P

P_par_bsp2 <- ggplot(df, aes(x = r.max_P, y = P.comp)) +  
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  geom_ribbon(data = new_data, aes(x = r.max_P, ymin = lower, ymax = upper), 
              fill = "grey75", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = preds, aes(x = r.max_P, y = Estimate), color = "blue", size = 1) +  # Bayesian quantile regression in blue
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation, Bayesian quantile regression (95% HDPI)") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_par_bsp2 # Display the plot

brm_gp <- brm(
  bf(P.comp ~ gp(r.max_P)),  
  data = par.res,
  family = gaussian(),
  iter = 6000, chains = 4, cores = 4
)

new_data <- data.frame(r.max_P = seq(min(par.res$r.max_P), max(par.res$r.max_P), length.out = 100))
preds <- as.data.frame(fitted(brm_gp, newdata = new_data, probs = c(0.05, 0.95))) 

# Add predictions to new_data
new_data$pred_P.comp <- preds[, "Estimate"]
new_data$lower <- preds[, "Q5"]   # 5% CI
new_data$upper <- preds[, "Q95"]  # 95% CI

preds$r.max_P <- new_data$r.max_P

P_par_gp <- ggplot(df, aes(x = r.max_P, y = P.comp)) +  
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  geom_ribbon(data = new_data, aes(x = r.max_P, ymin = lower, ymax = upper), 
              fill = "grey75", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = preds, aes(x = r.max_P, y = Estimate), color = "blue", size = 1) +  # Bayesian quantile regression in blue
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation, Bayesian quantile regression (95% HDPI)") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_par_gp # Display the plot

brm_ns <- brm(
  bf(P.comp ~ ns(r.max_P, df = 4)),  # Natural spline with 4 degrees of freedom
  data = par.res,
  family = gaussian(),
  iter = 4000, chains = 4, cores = 4
)

new_data <- data.frame(r.max_P = seq(min(par.res$r.max_P), max(par.res$r.max_P), length.out = 100))
preds <- as.data.frame(fitted(brm_ns, newdata = new_data, probs = c(0.05, 0.95))) 

# Add predictions to new_data
new_data$pred_P.comp <- preds[, "Estimate"]
new_data$lower <- preds[, "Q5"]   # 5% CI
new_data$upper <- preds[, "Q95"]  # 95% CI

preds$r.max_P <- new_data$r.max_P

P_par_ns <- ggplot(df, aes(x = r.max_P, y = P.comp)) +  
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  geom_ribbon(data = new_data, aes(x = r.max_P, ymin = lower, ymax = upper), 
              fill = "grey75", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = preds, aes(x = r.max_P, y = Estimate), color = "blue", size = 1) +  # Bayesian quantile regression in blue
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation, Bayesian quantile regression (95% HDPI)") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_par_ns # Display the plot

#OK let's try a polynomial function again?

brm_poly <- brm(
  bf(P.comp ~ poly(r.max_P, 2)),  # Quadratic regression (degree 2)
  data = par.res,
  family = gaussian(),
  iter = 4000, chains = 4, cores = 4
)

new_data <- data.frame(r.max_P = seq(min(par.res$r.max_P), max(par.res$r.max_P), length.out = 100))
preds <- as.data.frame(fitted(brm_poly, newdata = new_data, probs = c(0.05, 0.95))) 

# Add predictions to new_data
new_data$pred_P.comp <- preds[, "Estimate"]
new_data$lower <- preds[, "Q5"]   # 5% CI
new_data$upper <- preds[, "Q95"]  # 95% CI

preds$r.max_P <- new_data$r.max_P

P_par_poly <- ggplot(df, aes(x = r.max_P, y = P.comp)) +  
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  geom_ribbon(data = new_data, aes(x = r.max_P, ymin = lower, ymax = upper), 
              fill = "grey75", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = preds, aes(x = r.max_P, y = Estimate), color = "blue", size = 1) +  # Bayesian quantile regression in blue
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation, Bayesian quantile regression (95% HDPI)") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_par_poly # Display the plot

############# Salt ###########################################



# Nitrogen

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "N", 8, 16)) # Ps are now equivalent to 8, for later mapping

par.res <- par_frt(df, xvar = "r.max_N", yvar = "N.comp")

N_par <- ggplot(df, aes(x = r.max_N, y = N.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Nitrogen limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res, aes(x = r.max_N, y = N.comp), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

N_par # Raw pareto front.

# Light

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "L", 8, 16)) # Ls are now equivalent to 8, for later mapping

par.res <- par_frt(df[df$I.comp<10,], xvar = "r.max_I", yvar = "I.comp")

L_par <- ggplot(df[df$I.comp<10,], aes(x = r.max_I, y = I.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Light limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res, aes(x = r.max_I, y = I.comp), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

L_par # Raw pareto front.

# Salt

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "S", 8, 16)) # Ls are now equivalent to 8, for later mapping

par.res <- par_frt(df, xvar = "r.max_S", yvar = "S.c.mod")

S_par <- ggplot(df, aes(x = r.max_S, y = S.c.mod)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Salt tolerance (c)", 
       title = "Salt stress") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res, aes(x = r.max_S, y = S.c.mod), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

S_par # Raw pareto front.