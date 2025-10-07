# Jason R Laurich

# October 7th, 2025

# Creating Figure 3, a final version
# Will be based on 4x the data (model fitting for replicates)

# Will not be using error estimates in fitting Pareto fronts

# Significance testing: 

# Hyp 1 (existence of a Pareto front):
  # Li et al method
  # Quantile regression

# Hyp 2 (evolution moves PFs) :
  # 75th QR PF comparison

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(quantreg)
library(sp) # For the point.in.polygon function
library(scam)
library(pracma) # For calculating polygon area

par_frt <- function(df, xvar, yvar) { # Simple Pareto front function / convex hull algorithm (one sided)
  
  df <- df[order(-df[[xvar]], df[[yvar]]), ]  
  pareto_points <- df[1, ]  # Start with the first point
  
  for (i in 2:nrow(df)) {
    if (df[i, yvar] > tail(pareto_points[[yvar]], 1)) {  # Ensure increasing y values
      pareto_points <- rbind(pareto_points, df[i,])
    }
  }
  
  return(pareto_points)
}

find_nearest_index <- function(x, ref_vec) { # Function to find the closest point to the actual r.max value in the pred data frame. 
  which.min(abs(ref_vec - x))
} # For significance testing based on location of points relative to a curve

# Load & examine the data -------------------------------------------------

df <- read.csv("data-processed/53_summary_table.csv") # Summary file
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

# Temperature -------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral", "evolved") # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = T.br,
    z.x = T.µ.max
  ) # Specify the x and y variables and their 95% CIs

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.2 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

T.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
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
  
  #ylim(24.5,29.5) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.scam  # Display the plot

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

# Randomize across the whole dataset

null.df <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt %>%
    
    mutate(
      z.x.sim = sample(df.filt$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0!

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # -0.38886, p 0.00005
summary(q75, se = "boot", R = 1000) # -0.49592, p 0
summary(q90, se = "boot", R = 1000) # -0.64116, p 0.00006

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
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.qr # Display the plot

###### 75th quantile PF testing ######

# Cannot do for temperature data, as no populations evolved under temperature stress

# Light -------------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "L", 'light', 'other')) # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = I.comp,
    z.x = I.µ.max
  ) # Specify the x and y variables and their 95% CIs

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.2 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

I.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
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
  
  ylim(0.025, 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

I.scam  # Display the plot

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

# Randomize across the whole dataset

null.df <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt %>%
    
    mutate(
      z.x.sim = sample(df.filt$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.727

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 0.08277, p 0
summary(q75, se = "boot", R = 1000) # 0.10484, p 0
summary(q90, se = "boot", R = 1000) # 0.16664, p 0

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
  
  ylim(0.025, 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

I.qr # Display the plot

###### 75th quantile PF testing ######

n.75 <- df.filt %>% # Now calculate the number of light points above the 75th % quantile Pareto front. 
  filter(evol.bin == "light") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

null_counts <- numeric(1000)

for (i in 1:1000) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  shuffled.df <- df.filt %>%
    
    mutate(
      z.x.sim = sample(df.filt$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt$z.y, replace = FALSE)     # Seperately reassign y
    )
  
  null_counts[i] <- shuffled.df %>%
    filter(evol.bin == "light") %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x.sim)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y.sim > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.787

# Nitrogen ----------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "N", 'nitrogen', 'other')) # for testing significance of matching evolutionary conditions.

df.filt <- df %>% 
  mutate(
    z.y = N.comp,
    z.x = N.µ.max
  ) # Specify the x and y variables and their 95% CIs

df.filt <- df.filt %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.2 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

N.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
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
  
  ylim(0,1.25) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

N.scam  # Display the plot

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt$z.x) # Extract the max values for x and y.
y.max <- max(df.filt$z.y)

par.res.1 <- par.res.1 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
  arrange(z.x)

poly <- par.res.1[, c("z.x", "z.y")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
  add_row(z.x = x.max, z.y = y.max)

a.emp <- polyarea(poly$z.x, poly$z.y) # Calculate the area enclosed by these vertices

# Randomize across the whole dataset

null.df <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt %>%
    
    mutate(
      z.x.sim = sample(df.filt$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.727

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 0.08277, p 0
summary(q75, se = "boot", R = 1000) # 0.10484, p 0
summary(q90, se = "boot", R = 1000) # 0.16664, p 0

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
  
  ylim(0.025, 0.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

I.qr # Display the plot



