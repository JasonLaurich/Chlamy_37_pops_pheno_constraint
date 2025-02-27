# Jason R Laurich
# February 24, 2025

# OK so here I am going upload data for TPC shapes, Light, Nitrogen, and Phosphorous Monod curves, and Salt tolerance data.
# I'll look at some generalist-specialist and gleaner-opportunist trade-offs here.

# Additionally, I'll look at the data closely and troubleshoot (e.g. looking for outliers, figuring out why µmax is so much higher for TPCs)
# And I'll do a PCA and RDA of these data, generating some plots.

# Finally, I will want to model Pareto fronts on these data more explicitly — more on that later, still figuring this out. 

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

df.hist <- read.csv("data-processed/chlamee-acute-exponential.csv") # This file has the history (ancestry, evolutionary treatment) for each well ID.
head(df.hist)
str(df.hist)
df.hist$anc <- as.factor(df.hist$ancestor_id) # Will port these into combined dataset, remember: we are sorting by population. 
df.hist$evol <- as.factor(df.hist$treatment)  # so we really only need population level - evol - anc mapping, not well-specific.
df.hist$Pop.fac <- as.factor(df.hist$population)

df.tpc <- read.csv("data-processed/09_TPC_shape_values_BayesTPC.csv") # TPC summary data. We'll use Lactin II model fits. 
head(df.tpc)
str(df.tpc)

df.tpc <- df.tpc[df.tpc$Pop.fac != 'cc1629' & df.tpc$Model == 'Lactin 2', ] # Trim off Thomas models and that one weird population.

df <- df.tpc[,c(2,5:9)] # We'll work with the analytically determined parameters (Deriv package)

df.hist.agg <- df.hist %>% # Aggregate the history data for merging.
  group_by(Pop.fac) %>%
  summarise(
    anc = first(anc),
    evol = first(evol),
    .groups = "drop" 
  )

df <- merge(df, df.hist.agg, by = "Pop.fac", all.x = TRUE)
head(df) 
str(df)
df$Pop.fac <- as.factor(df$Pop.fac) # Looks great! Now we layer in other summary metrics.
names(df) <- c("Pop.fac", "T.min", "T.max", "T.br", "Topt", "r.max_T", "anc", "evol") # rename a few variables. 

df.I <- read.csv("data-processed/10b_light_Monod_estimates.csv") # I* data 
head(df.I)
str(df.I)

df.I.par <- df.I[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.I.par) <- c("Pop.fac", "K.s_I", "r.max_I", "I.jag", "I.mth") # rename
df.I.par$I.comp <- 1/df.I.par$I.mth

df <- merge(df, df.I.par, by = "Pop.fac", all.x = TRUE)

df.N <- read.csv("data-processed/11b_nitrogen_Monod_estimates.csv") # N* data 
head(df.N)
str(df.N)

df.N.par <- df.N[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.N.par) <- c("Pop.fac", "K.s_N", "r.max_N", "N.jag", "N.mth") # rename
df.N.par$N.comp <- 1/df.N.par$N.mth

df <- merge(df, df.N.par, by = "Pop.fac", all.x = TRUE)

df.P <- read.csv("data-processed/12b_phosphorous_Monod_estimates.csv") # P* data 
head(df.P)
str(df.P)

df.P.par <- df.P[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.P.par) <- c("Pop.fac", "K.s_P", "r.max_P", "P.jag", "P.mth") # rename
df.P.par$P.comp <- 1/df.P.par$P.mth

df <- merge(df, df.P.par, by = "Pop.fac", all.x = TRUE)

df.S <- read.csv("data-Processed/13b_salt_tolerance_estimates.csv") # S* data 
head(df.S)
str(df.S)

df.S$Pop.fac <- gsub("Anc ", "anc", df.S$Pop.fac) # ancestors are written differently here

df.S.par <- df.S[, c("Pop.fac", "r.max", "c.mod", "c.pred")] # need these
names(df.S.par) <- c("Pop.fac", "r.max_S", "S.c.mod", "S.c.pred") # rename

df <- merge(df, df.S.par, by = "Pop.fac", all.x = TRUE)

write.csv(df, "data-processed/14_summary_metric_table.csv") # Save summary table

############# Let's visualize some trade-offs #######################

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


# Phosphorous

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

new_data <- data.frame(r.max_T = seq(min(par.res$r.max_T), max(par.res$r.max_T), length.out = 100))


new_data <- data.frame(r.max_P = seq(min(par.res$r.max_P), max(par.res$r.max_P), length.out = 100))
preds <- as.data.frame(fitted(mod_bayes_qr, newdata = new_data, probs = c(0.05, 0.95))) 

preds3 <- as.data.frame(fitted(mod_bayes_qr3, newdata = new_data, probs = c(0.05, 0.95))) 

plot(preds3$Estimate ~ preds3$r.max_P)

# Add predictions to new_data
new_data$pred_br <- preds3[, "Estimate"]
new_data$lower <- preds3[, "Q5"]   # 5% CI
new_data$upper <- preds3[, "Q95"]  # 95% CI
preds3$r.max_P <- new_data$r.max_P  # Attach r.max_P back

P_par_bqr <- ggplot(df, aes(x = r.max_P, y = P.comp)) +  
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  geom_line(data = preds3, aes(x = r.max_P, y = Estimate), color = "blue", size = 1) +  # Bayesian quantile regression in blue
  geom_line(data = preds3, aes(x = r.max_P, y = Q5), color = "purple", size = 1) +
  geom_line(data = preds3, aes(x = r.max_P, y = Q95), color = "purple", size = 1) +
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation, Bayesian quantile regression (95% HDPI)") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_par_bqr # Display the plot

# Nitrogen

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "N", 8, 16)) # Ps are now equivalent to 8, for later mapping

par.res <- par_frt(df, xvar = "r.max_N", yvar = "N.comp")

N_par1 <- ggplot(df, aes(x = r.max_N, y = N.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Nitrogen limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res, aes(x = r.max_N, y = N.comp), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

N_par1 # Raw pareto front.

# Light

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "L", 8, 16)) # Ls are now equivalent to 8, for later mapping

par.res <- par_frt(df[df$I.comp<10,], xvar = "r.max_I", yvar = "I.comp")

L_par1 <- ggplot(df[df$I.comp<10,], aes(x = r.max_I, y = I.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Light limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res, aes(x = r.max_I, y = I.comp), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

L_par1 # Raw pareto front.

# Salt

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "S", 8, 16)) # Ls are now equivalent to 8, for later mapping

par.res <- par_frt(df, xvar = "r.max_S", yvar = "S.c.mod")

S_par1 <- ggplot(df[df$I.comp<10,], aes(x = r.max_S, y = S.c.mod)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Salt tolerance (c)", 
       title = "Salt stress") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  geom_line(data = par.res, aes(x = r.max_S, y = S.c.mod), color = "blue", size = 1) +  # Pareto frontier line
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

S_par1 # Raw pareto front.
