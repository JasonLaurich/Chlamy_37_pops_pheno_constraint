# Jason R Laurich

# September 8th, 2025

# Creating Figure 5: Inter-gradient Pareto front panels (light, nitrogen, phosphorous, salt, and temperature)
# Am also going to include relevant statistical testing here. 

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(cowplot)
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

# Light comparisons -------------------------------------------------------

###### Light v Nitrogen ######

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                            ifelse(df$Evol == "L", 'light', ifelse(df$Evol == 'N', 'nit', 'other'))) # for testing.

df.filt <- df %>% 
  mutate(
    z.y = I.comp.0.56,
    z.y.l = I.comp.0.56.L,
    z.y.u = I.comp.0.56.U,
    z.x = N.comp.0.56,
    z.x.l = N.comp.0.56.L,
    z.x.u = N.comp.0.56.U
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
    dist.sc = distance/mean(distance)
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

df.filt <- df.filt %>% # Trim dramatic outliers as errors
  filter(dist.sc < 2.4) 

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale if necessary

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

par.res.1 <- par.res.1 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.1 <- rbind(                               # Build a polygon enclosing the space above the pareto front
  cbind(par.res.1$z.x, par.res.1$z.y),             # along the frontier
  c(max(par.res.1$z.x)*3, min(par.res.1$z.y)),     # over at the right end
  c(max(par.res.1$z.x)*3, max(par.res.1$z.y)*3),   # up to the top
  c(min(par.res.1$z.x), max(par.res.1$z.y)*3),     # over to the left end
  c(par.res.1$z.x[1], par.res.1$z.y[1])            # close it
)

par.res.2 <- df.filt %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt$z.x + (df.filt$z.x.u - df.filt$z.x)/3, df.filt$z.y + (df.filt$z.y.u - df.filt$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.1[,1], poly.par.1[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Fit a scam to the adjusted PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

par.res.3 <- par.res.3 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.2 <- rbind(                                   # Build a polygon enclosing the space above the pareto front
  cbind(par.res.3$z.x, par.res.3$z.y),                 # along the frontier
  c(max(par.res.3$z.x)*3, min(par.res.3$z.y)),         # over at the right end
  c(max(par.res.3$z.x)*3, max(par.res.3$z.y)*3),       # up to the top
  c(min(par.res.3$z.x), max(par.res.3$z.y)*3),         # over to the left end
  c(par.res.3$z.x[1], par.res.3$z.y[1])                # close it
)

par.res.4 <- df.filt2 %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt2$z.x + (df.filt2$z.x.u - df.filt2$z.x)/3, df.filt2$z.y + (df.filt2$z.y.u - df.filt2$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.2[,1], poly.par.2[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.4))), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

LN.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Competitive ability (1/N*)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "A — Light ~ Nitrogen") +  # labels
  
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
               "nit" = 16,
               "light" = 16)  # filled circle
  ) +
  
  ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LN.scam.PF  # Display the plot

###### 75th quantile PF testing ######

n.75 <- df.filt %>% # Now calculate the number of light and nitrogen points above the 75th % quantile Pareto front. 
  filter(evol.bin %in% c("light", "nit")) %>%
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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  null_counts[i] <- shuffled.df %>%
    filter(evol.bin %in% c("light", "nit")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x.sim)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y.sim > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.909

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # Scale for point exclusion determination 
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # re-scale before we fit the pareto front (to remove skew caused by drastic outliers)
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.865

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(      # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.2 %>%
    
    mutate(
      n.x = sample(nrow(par.res.2), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(par.res.2),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(par.res.2), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(par.res.2),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  par.res.n2 <- par_frt(shuffled.df2, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n2 <- par.res.n2 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n2 <- max(shuffled.df2$z.x.sim) # Extract the max values for x and y
  y.max.n2 <- max(shuffled.df2$z.y.sim)
  
  poly.n2 <- par.res.n2[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n2, z.y.sim = y.max.n2)
  
  a.emp.n2 <- polyarea(poly.n2$z.x.sim, poly.n2$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df2 <- rbind(null.df2, data.frame(  # Save the data
    a.emp.n = a.emp.n2,                    # Area above the curve 
    n.PF = nrow(par.res.n2)                # Number of data points in the PF
  ))
  
}

mean(null.df2$a.emp.n >= a.emp) # p-value 0.794

###### Light v Phosphorous ######

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "L", 'light', ifelse(df$Evol == 'P', 'phos', 'other'))) # for testing.

df.filt <- df %>% 
  mutate(
    z.y = I.comp.0.56,
    z.y.l = I.comp.0.56.L,
    z.y.u = I.comp.0.56.U,
    z.x = P.comp.0.56,
    z.x.l = P.comp.0.56.L,
    z.x.u = P.comp.0.56.U
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
    dist.sc = distance/mean(distance)
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

df.filt <- df.filt %>% # Trim dramatic outliers as errors
  filter(dist.sc < 2.4) 

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale if necessary

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

par.res.1 <- par.res.1 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.1 <- rbind(                               # Build a polygon enclosing the space above the pareto front
  cbind(par.res.1$z.x, par.res.1$z.y),             # along the frontier
  c(max(par.res.1$z.x)*3, min(par.res.1$z.y)),     # over at the right end
  c(max(par.res.1$z.x)*3, max(par.res.1$z.y)*3),   # up to the top
  c(min(par.res.1$z.x), max(par.res.1$z.y)*3),     # over to the left end
  c(par.res.1$z.x[1], par.res.1$z.y[1])            # close it
)

par.res.2 <- df.filt %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt$z.x + (df.filt$z.x.u - df.filt$z.x)/3, df.filt$z.y + (df.filt$z.y.u - df.filt$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.1[,1], poly.par.1[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Fit a scam to the adjusted PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

par.res.3 <- par.res.3 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.2 <- rbind(                                   # Build a polygon enclosing the space above the pareto front
  cbind(par.res.3$z.x, par.res.3$z.y),                 # along the frontier
  c(max(par.res.3$z.x)*3, min(par.res.3$z.y)),         # over at the right end
  c(max(par.res.3$z.x)*3, max(par.res.3$z.y)*3),       # up to the top
  c(min(par.res.3$z.x), max(par.res.3$z.y)*3),         # over to the left end
  c(par.res.3$z.x[1], par.res.3$z.y[1])                # close it
)

par.res.4 <- df.filt2 %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt2$z.x + (df.filt2$z.x.u - df.filt2$z.x)/3, df.filt2$z.y + (df.filt2$z.y.u - df.filt2$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.2[,1], poly.par.2[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.4))), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

LP.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "B — Light ~ Phosphorous") +  # labels
  
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
               "phos" = 16,
               "light" = 16)  # filled circle
  ) +
  
  ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LP.scam.PF  # Display the plot

###### 75th quantile PF testing ######

n.75 <- df.filt %>% # Now calculate the number of light and phosphorous points above the 75th % quantile Pareto front. 
  filter(evol.bin %in% c("light", "phos")) %>%
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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  null_counts[i] <- shuffled.df %>%
    filter(evol.bin %in% c("light", "phos")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x.sim)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y.sim > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.785

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # Scale for point exclusion determination 
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # re-scale before we fit the pareto front (to remove skew caused by drastic outliers)
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.961

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(      # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.2 %>%
    
    mutate(
      n.x = sample(nrow(par.res.2), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(par.res.2),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(par.res.2), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(par.res.2),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  par.res.n2 <- par_frt(shuffled.df2, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n2 <- par.res.n2 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n2 <- max(shuffled.df2$z.x.sim) # Extract the max values for x and y
  y.max.n2 <- max(shuffled.df2$z.y.sim)
  
  poly.n2 <- par.res.n2[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n2, z.y.sim = y.max.n2)
  
  a.emp.n2 <- polyarea(poly.n2$z.x.sim, poly.n2$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df2 <- rbind(null.df2, data.frame(  # Save the data
    a.emp.n = a.emp.n2,                    # Area above the curve 
    n.PF = nrow(par.res.n2)                # Number of data points in the PF
  ))
  
}

mean(null.df2$a.emp.n >= a.emp) # p-value 0.720

###### Light v Salt ######

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "L", 'light', ifelse(df$Evol %in% c("S", "BS"), 'salt', 'other'))) # for testing.

df.filt <- df %>% 
  mutate(
    z.y = I.comp.0.56,
    z.y.l = I.comp.0.56.L,
    z.y.u = I.comp.0.56.U,
    z.x = S.c,
    z.x.l = S.c.L,
    z.x.u = S.c.U
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
    dist.sc = distance/mean(distance)
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

df.filt <- df.filt %>% # Trim dramatic outliers as errors
  filter(dist.sc < 2.4) 

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale if necessary

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

par.res.1 <- par.res.1 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.1 <- rbind(                               # Build a polygon enclosing the space above the pareto front
  cbind(par.res.1$z.x, par.res.1$z.y),             # along the frontier
  c(max(par.res.1$z.x)*3, min(par.res.1$z.y)),     # over at the right end
  c(max(par.res.1$z.x)*3, max(par.res.1$z.y)*3),   # up to the top
  c(min(par.res.1$z.x), max(par.res.1$z.y)*3),     # over to the left end
  c(par.res.1$z.x[1], par.res.1$z.y[1])            # close it
)

par.res.2 <- df.filt %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt$z.x + (df.filt$z.x.u - df.filt$z.x)/3, df.filt$z.y + (df.filt$z.y.u - df.filt$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.1[,1], poly.par.1[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Fit a scam to the adjusted PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

par.res.3 <- par.res.3 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.2 <- rbind(                                   # Build a polygon enclosing the space above the pareto front
  cbind(par.res.3$z.x, par.res.3$z.y),                 # along the frontier
  c(max(par.res.3$z.x)*3, min(par.res.3$z.y)),         # over at the right end
  c(max(par.res.3$z.x)*3, max(par.res.3$z.y)*3),       # up to the top
  c(min(par.res.3$z.x), max(par.res.3$z.y)*3),         # over to the left end
  c(par.res.3$z.x[1], par.res.3$z.y[1])                # close it
)

par.res.4 <- df.filt2 %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt2$z.x + (df.filt2$z.x.u - df.filt2$z.x)/3, df.filt2$z.y + (df.filt2$z.y.u - df.filt2$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.2[,1], poly.par.2[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.4))), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

LS.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Salt tolerance (c)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "C — Light ~ Salt") +  # labels
  
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
               "salt" = 16,
               "light" = 16)  # filled circle
  ) +
  
  ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LS.scam.PF  # Display the plot

###### 75th quantile PF testing ######

n.75 <- df.filt %>% # Now calculate the number of light and phosphorous points above the 75th % quantile Pareto front. 
  filter(evol.bin %in% c("light", "salt")) %>%
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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  null_counts[i] <- shuffled.df %>%
    filter(evol.bin %in% c("light", "salt")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x.sim)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y.sim > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.785

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # Scale for point exclusion determination 
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # re-scale before we fit the pareto front (to remove skew caused by drastic outliers)
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.870

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(      # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.2 %>%
    
    mutate(
      n.x = sample(nrow(par.res.2), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(par.res.2),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(par.res.2), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(par.res.2),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  par.res.n2 <- par_frt(shuffled.df2, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n2 <- par.res.n2 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n2 <- max(shuffled.df2$z.x.sim) # Extract the max values for x and y
  y.max.n2 <- max(shuffled.df2$z.y.sim)
  
  poly.n2 <- par.res.n2[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n2, z.y.sim = y.max.n2)
  
  a.emp.n2 <- polyarea(poly.n2$z.x.sim, poly.n2$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df2 <- rbind(null.df2, data.frame(  # Save the data
    a.emp.n = a.emp.n2,                    # Area above the curve 
    n.PF = nrow(par.res.n2)                # Number of data points in the PF
  ))
  
}

mean(null.df2$a.emp.n >= a.emp) # p-value 0.720

###### Light v Temperature ######

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "L", 'light', 'other')) # for testing.

df.filt <- df %>% 
  mutate(
    z.y = I.comp.0.56,
    z.y.l = I.comp.0.56.L,
    z.y.u = I.comp.0.56.U,
    z.x = T.br.0.56,
    z.x.l = T.br.0.56.L,
    z.x.u = T.br.0.56.U
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
    dist.sc = distance/mean(distance)
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

df.filt <- df.filt %>% # Trim dramatic outliers as errors
  filter(dist.sc < 2.4) 

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale if necessary

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

par.res.1 <- par.res.1 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.1 <- rbind(                               # Build a polygon enclosing the space above the pareto front
  cbind(par.res.1$z.x, par.res.1$z.y),             # along the frontier
  c(max(par.res.1$z.x)*3, min(par.res.1$z.y)),     # over at the right end
  c(max(par.res.1$z.x)*3, max(par.res.1$z.y)*3),   # up to the top
  c(min(par.res.1$z.x), max(par.res.1$z.y)*3),     # over to the left end
  c(par.res.1$z.x[1], par.res.1$z.y[1])            # close it
)

par.res.2 <- df.filt %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt$z.x + (df.filt$z.x.u - df.filt$z.x)/3, df.filt$z.y + (df.filt$z.y.u - df.filt$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.1[,1], poly.par.1[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Fit a scam to the adjusted PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

par.res.3 <- par.res.3 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.2 <- rbind(                                   # Build a polygon enclosing the space above the pareto front
  cbind(par.res.3$z.x, par.res.3$z.y),                 # along the frontier
  c(max(par.res.3$z.x)*3, min(par.res.3$z.y)),         # over at the right end
  c(max(par.res.3$z.x)*3, max(par.res.3$z.y)*3),       # up to the top
  c(min(par.res.3$z.x), max(par.res.3$z.y)*3),         # over to the left end
  c(par.res.3$z.x[1], par.res.3$z.y[1])                # close it
)

par.res.4 <- df.filt2 %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt2$z.x + (df.filt2$z.x.u - df.filt2$z.x)/3, df.filt2$z.y + (df.filt2$z.y.u - df.filt2$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.2[,1], poly.par.2[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.4))), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

LT.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Thermal breadth (°C)",    
       y = "Competitive ability (1/I*)", 
       color = "Evolutionary History",
       title = "D — Light ~ Temperature") +  # labels
  
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
  
  ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LT.scam.PF  # Display the plot

###### 75th quantile PF testing ######

n.75 <- df.filt %>% # Now calculate the number of light and phosphorous points above the 75th % quantile Pareto front. 
  filter(evol.bin == "light") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n.75 # There is something wrong here....

null_counts <- numeric(1000)

for (i in 1:1000) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  shuffled.df <- df.filt %>%
    
    mutate(
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  null_counts[i] <- shuffled.df %>%
    filter(evol.bin %in% c("light")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x.sim)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y.sim > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.785

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # Scale for point exclusion determination 
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # re-scale before we fit the pareto front (to remove skew caused by drastic outliers)
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.870

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(      # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.2 %>%
    
    mutate(
      n.x = sample(nrow(par.res.2), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(par.res.2),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(par.res.2), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(par.res.2),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  par.res.n2 <- par_frt(shuffled.df2, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n2 <- par.res.n2 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n2 <- max(shuffled.df2$z.x.sim) # Extract the max values for x and y
  y.max.n2 <- max(shuffled.df2$z.y.sim)
  
  poly.n2 <- par.res.n2[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n2, z.y.sim = y.max.n2)
  
  a.emp.n2 <- polyarea(poly.n2$z.x.sim, poly.n2$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df2 <- rbind(null.df2, data.frame(  # Save the data
    a.emp.n = a.emp.n2,                    # Area above the curve 
    n.PF = nrow(par.res.n2)                # Number of data points in the PF
  ))
  
}

mean(null.df2$a.emp.n >= a.emp) # p-value 0.720

# Nitrogen comparisons ----------------------------------------------------

###### Nitrogen v Phosphorous ######

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "N", 'nit', ifelse(df$Evol == 'P', 'phos', 'other'))) # for testing.

df.filt <- df %>% 
  mutate(
    z.y = N.comp.0.56,
    z.y.l = N.comp.0.56.L,
    z.y.u = N.comp.0.56.U,
    z.x = P.comp.0.56,
    z.x.l = P.comp.0.56.L,
    z.x.u = P.comp.0.56.U
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
    dist.sc = distance/mean(distance)
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

df.filt <- df.filt %>% # Trim dramatic outliers as errors
  filter(dist.sc < 2.4) 

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale if necessary

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

par.res.1 <- par.res.1 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.1 <- rbind(                               # Build a polygon enclosing the space above the pareto front
  cbind(par.res.1$z.x, par.res.1$z.y),             # along the frontier
  c(max(par.res.1$z.x)*3, min(par.res.1$z.y)),     # over at the right end
  c(max(par.res.1$z.x)*3, max(par.res.1$z.y)*3),   # up to the top
  c(min(par.res.1$z.x), max(par.res.1$z.y)*3),     # over to the left end
  c(par.res.1$z.x[1], par.res.1$z.y[1])            # close it
)

par.res.2 <- df.filt %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt$z.x + (df.filt$z.x.u - df.filt$z.x)/3, df.filt$z.y + (df.filt$z.y.u - df.filt$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.1[,1], poly.par.1[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Fit a scam to the adjusted PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

par.res.3 <- par.res.3 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.2 <- rbind(                                   # Build a polygon enclosing the space above the pareto front
  cbind(par.res.3$z.x, par.res.3$z.y),                 # along the frontier
  c(max(par.res.3$z.x)*3, min(par.res.3$z.y)),         # over at the right end
  c(max(par.res.3$z.x)*3, max(par.res.3$z.y)*3),       # up to the top
  c(min(par.res.3$z.x), max(par.res.3$z.y)*3),         # over to the left end
  c(par.res.3$z.x[1], par.res.3$z.y[1])                # close it
)

par.res.4 <- df.filt2 %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt2$z.x + (df.filt2$z.x.u - df.filt2$z.x)/3, df.filt2$z.y + (df.filt2$z.y.u - df.filt2$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.2[,1], poly.par.2[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.4))), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

NP.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "E — Nitrogen ~ Phosphorous") +  # labels
  
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
               "phos" = 16,
               "nit" = 16)  # filled circle
  ) +
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NP.scam.PF  # Display the plot

###### 75th quantile PF testing ######

n.75 <- df.filt %>% # Now calculate the number of light and phosphorous points above the 75th % quantile Pareto front. 
  filter(evol.bin %in% c("nit", "phos")) %>%
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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  null_counts[i] <- shuffled.df %>%
    filter(evol.bin %in% c("light", "phos")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x.sim)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y.sim > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.000

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # Scale for point exclusion determination 
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # re-scale before we fit the pareto front (to remove skew caused by drastic outliers)
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.961

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(      # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.2 %>%
    
    mutate(
      n.x = sample(nrow(par.res.2), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(par.res.2),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(par.res.2), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(par.res.2),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  par.res.n2 <- par_frt(shuffled.df2, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n2 <- par.res.n2 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n2 <- max(shuffled.df2$z.x.sim) # Extract the max values for x and y
  y.max.n2 <- max(shuffled.df2$z.y.sim)
  
  poly.n2 <- par.res.n2[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n2, z.y.sim = y.max.n2)
  
  a.emp.n2 <- polyarea(poly.n2$z.x.sim, poly.n2$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df2 <- rbind(null.df2, data.frame(  # Save the data
    a.emp.n = a.emp.n2,                    # Area above the curve 
    n.PF = nrow(par.res.n2)                # Number of data points in the PF
  ))
  
}

mean(null.df2$a.emp.n >= a.emp) # p-value 0.720

###### Nitrogen v Salt ######

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "N", 'nit', ifelse(df$Evol %in% c("S", "BS"), 'salt', 'other'))) # for testing.

df.filt <- df %>% 
  mutate(
    z.y = N.comp.0.56,
    z.y.l = N.comp.0.56.L,
    z.y.u = N.comp.0.56.U,
    z.x = S.c,
    z.x.l = S.c.L,
    z.x.u = S.c.U
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
    dist.sc = distance/mean(distance)
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

df.filt <- df.filt %>% # Trim dramatic outliers as errors
  filter(dist.sc < 2.4) 

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale if necessary

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

par.res.1 <- par.res.1 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.1 <- rbind(                               # Build a polygon enclosing the space above the pareto front
  cbind(par.res.1$z.x, par.res.1$z.y),             # along the frontier
  c(max(par.res.1$z.x)*3, min(par.res.1$z.y)),     # over at the right end
  c(max(par.res.1$z.x)*3, max(par.res.1$z.y)*3),   # up to the top
  c(min(par.res.1$z.x), max(par.res.1$z.y)*3),     # over to the left end
  c(par.res.1$z.x[1], par.res.1$z.y[1])            # close it
)

par.res.2 <- df.filt %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt$z.x + (df.filt$z.x.u - df.filt$z.x)/3, df.filt$z.y + (df.filt$z.y.u - df.filt$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.1[,1], poly.par.1[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Fit a scam to the adjusted PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

par.res.3 <- par.res.3 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.2 <- rbind(                                   # Build a polygon enclosing the space above the pareto front
  cbind(par.res.3$z.x, par.res.3$z.y),                 # along the frontier
  c(max(par.res.3$z.x)*3, min(par.res.3$z.y)),         # over at the right end
  c(max(par.res.3$z.x)*3, max(par.res.3$z.y)*3),       # up to the top
  c(min(par.res.3$z.x), max(par.res.3$z.y)*3),         # over to the left end
  c(par.res.3$z.x[1], par.res.3$z.y[1])                # close it
)

par.res.4 <- df.filt2 %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt2$z.x + (df.filt2$z.x.u - df.filt2$z.x)/3, df.filt2$z.y + (df.filt2$z.y.u - df.filt2$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.2[,1], poly.par.2[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.4))), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

NS.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Salt tolerance (c)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "F — Nitrogen ~ Salt") +  # labels
  
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
               "salt" = 16,
               "nit" = 16)  # filled circle
  ) +
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NS.scam.PF  # Display the plot

###### 75th quantile PF testing ######

n.75 <- df.filt %>% # Now calculate the number of light and phosphorous points above the 75th % quantile Pareto front. 
  filter(evol.bin %in% c("nit", "salt")) %>%
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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  null_counts[i] <- shuffled.df %>%
    filter(evol.bin %in% c("nit", "salt")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x.sim)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y.sim > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.785

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # Scale for point exclusion determination 
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # re-scale before we fit the pareto front (to remove skew caused by drastic outliers)
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.870

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(      # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.2 %>%
    
    mutate(
      n.x = sample(nrow(par.res.2), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(par.res.2),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(par.res.2), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(par.res.2),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  par.res.n2 <- par_frt(shuffled.df2, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n2 <- par.res.n2 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n2 <- max(shuffled.df2$z.x.sim) # Extract the max values for x and y
  y.max.n2 <- max(shuffled.df2$z.y.sim)
  
  poly.n2 <- par.res.n2[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n2, z.y.sim = y.max.n2)
  
  a.emp.n2 <- polyarea(poly.n2$z.x.sim, poly.n2$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df2 <- rbind(null.df2, data.frame(  # Save the data
    a.emp.n = a.emp.n2,                    # Area above the curve 
    n.PF = nrow(par.res.n2)                # Number of data points in the PF
  ))
  
}

mean(null.df2$a.emp.n >= a.emp) # p-value 0.720

###### Nitrogen v Temperature ######

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "N", 'nit', 'other')) # for testing.

df.filt <- df %>% 
  mutate(
    z.y = N.comp.0.56,
    z.y.l = N.comp.0.56.L,
    z.y.u = N.comp.0.56.U,
    z.x = T.br.0.56,
    z.x.l = T.br.0.56.L,
    z.x.u = T.br.0.56.U
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
    dist.sc = distance/mean(distance)
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

df.filt <- df.filt %>% # Trim dramatic outliers as errors
  filter(dist.sc < 2.4) 

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale if necessary

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

par.res.1 <- par.res.1 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.1 <- rbind(                               # Build a polygon enclosing the space above the pareto front
  cbind(par.res.1$z.x, par.res.1$z.y),             # along the frontier
  c(max(par.res.1$z.x)*3, min(par.res.1$z.y)),     # over at the right end
  c(max(par.res.1$z.x)*3, max(par.res.1$z.y)*3),   # up to the top
  c(min(par.res.1$z.x), max(par.res.1$z.y)*3),     # over to the left end
  c(par.res.1$z.x[1], par.res.1$z.y[1])            # close it
)

par.res.2 <- df.filt %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt$z.x + (df.filt$z.x.u - df.filt$z.x)/3, df.filt$z.y + (df.filt$z.y.u - df.filt$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.1[,1], poly.par.1[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Fit a scam to the adjusted PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

par.res.3 <- par.res.3 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.2 <- rbind(                                   # Build a polygon enclosing the space above the pareto front
  cbind(par.res.3$z.x, par.res.3$z.y),                 # along the frontier
  c(max(par.res.3$z.x)*3, min(par.res.3$z.y)),         # over at the right end
  c(max(par.res.3$z.x)*3, max(par.res.3$z.y)*3),       # up to the top
  c(min(par.res.3$z.x), max(par.res.3$z.y)*3),         # over to the left end
  c(par.res.3$z.x[1], par.res.3$z.y[1])                # close it
)

par.res.4 <- df.filt2 %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt2$z.x + (df.filt2$z.x.u - df.filt2$z.x)/3, df.filt2$z.y + (df.filt2$z.y.u - df.filt2$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.2[,1], poly.par.2[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.4))), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

NT.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Thermal breadth  (°C)",    
       y = "Competitive ability (1/N*)", 
       color = "Evolutionary History",
       title = "G — Nitrogen ~ Temperature") +  # labels
  
  ylim(0,1.7) +
  
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
               "nit" = 16)  # filled circle
  ) +
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NT.scam.PF  # Display the plot

###### 75th quantile PF testing ######

n.75 <- df.filt %>% # Now calculate the number of light and phosphorous points above the 75th % quantile Pareto front. 
  filter(evol.bin == "nit") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n.75 # There is something wrong here....

null_counts <- numeric(1000)

for (i in 1:1000) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  shuffled.df <- df.filt %>%
    
    mutate(
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  null_counts[i] <- shuffled.df %>%
    filter(evol.bin %in% c("nit")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x.sim)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y.sim > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.785

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # Scale for point exclusion determination 
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # re-scale before we fit the pareto front (to remove skew caused by drastic outliers)
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.870

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(      # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.2 %>%
    
    mutate(
      n.x = sample(nrow(par.res.2), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(par.res.2),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(par.res.2), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(par.res.2),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  par.res.n2 <- par_frt(shuffled.df2, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n2 <- par.res.n2 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n2 <- max(shuffled.df2$z.x.sim) # Extract the max values for x and y
  y.max.n2 <- max(shuffled.df2$z.y.sim)
  
  poly.n2 <- par.res.n2[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n2, z.y.sim = y.max.n2)
  
  a.emp.n2 <- polyarea(poly.n2$z.x.sim, poly.n2$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df2 <- rbind(null.df2, data.frame(  # Save the data
    a.emp.n = a.emp.n2,                    # Area above the curve 
    n.PF = nrow(par.res.n2)                # Number of data points in the PF
  ))
  
}

mean(null.df2$a.emp.n >= a.emp) # p-value 0.720

# Phosphorous comparisons -------------------------------------------------

###### Phosphorous v Salt ######

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "P", 'phos', ifelse(df$Evol %in% c("S", "BS"), 'salt', 'other'))) # for testing.

df.filt <- df %>% 
  mutate(
    z.y = P.comp.0.56,
    z.y.l = P.comp.0.56.L,
    z.y.u = P.comp.0.56.U,
    z.x = S.c,
    z.x.l = S.c.L,
    z.x.u = S.c.U
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
    dist.sc = distance/mean(distance)
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

df.filt <- df.filt %>% # Trim dramatic outliers as errors
  filter(dist.sc < 2.4) 

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale if necessary

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

par.res.1 <- par.res.1 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.1 <- rbind(                               # Build a polygon enclosing the space above the pareto front
  cbind(par.res.1$z.x, par.res.1$z.y),             # along the frontier
  c(max(par.res.1$z.x)*3, min(par.res.1$z.y)),     # over at the right end
  c(max(par.res.1$z.x)*3, max(par.res.1$z.y)*3),   # up to the top
  c(min(par.res.1$z.x), max(par.res.1$z.y)*3),     # over to the left end
  c(par.res.1$z.x[1], par.res.1$z.y[1])            # close it
)

par.res.2 <- df.filt %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt$z.x + (df.filt$z.x.u - df.filt$z.x)/3, df.filt$z.y + (df.filt$z.y.u - df.filt$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.1[,1], poly.par.1[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Fit a scam to the adjusted PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

par.res.3 <- par.res.3 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.2 <- rbind(                                   # Build a polygon enclosing the space above the pareto front
  cbind(par.res.3$z.x, par.res.3$z.y),                 # along the frontier
  c(max(par.res.3$z.x)*3, min(par.res.3$z.y)),         # over at the right end
  c(max(par.res.3$z.x)*3, max(par.res.3$z.y)*3),       # up to the top
  c(min(par.res.3$z.x), max(par.res.3$z.y)*3),         # over to the left end
  c(par.res.3$z.x[1], par.res.3$z.y[1])                # close it
)

par.res.4 <- df.filt2 %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt2$z.x + (df.filt2$z.x.u - df.filt2$z.x)/3, df.filt2$z.y + (df.filt2$z.y.u - df.filt2$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.2[,1], poly.par.2[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.4))), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

PS.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Salt tolerance (c)",    
       y = "Competitive ability (1/P*)", 
       color = "Evolutionary History",
       title = "H — Phosphorous ~ Salt") +  # labels
  
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
               "salt" = 16,
               "phos" = 16)  # filled circle
  ) +
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

PS.scam.PF  # Display the plot

###### 75th quantile PF testing ######

n.75 <- df.filt %>% # Now calculate the number of light and phosphorous points above the 75th % quantile Pareto front. 
  filter(evol.bin %in% c("phos", "salt")) %>%
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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  null_counts[i] <- shuffled.df %>%
    filter(evol.bin %in% c("phos", "salt")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x.sim)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y.sim > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.785

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # Scale for point exclusion determination 
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # re-scale before we fit the pareto front (to remove skew caused by drastic outliers)
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.870

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(      # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.2 %>%
    
    mutate(
      n.x = sample(nrow(par.res.2), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(par.res.2),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(par.res.2), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(par.res.2),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  par.res.n2 <- par_frt(shuffled.df2, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n2 <- par.res.n2 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n2 <- max(shuffled.df2$z.x.sim) # Extract the max values for x and y
  y.max.n2 <- max(shuffled.df2$z.y.sim)
  
  poly.n2 <- par.res.n2[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n2, z.y.sim = y.max.n2)
  
  a.emp.n2 <- polyarea(poly.n2$z.x.sim, poly.n2$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df2 <- rbind(null.df2, data.frame(  # Save the data
    a.emp.n = a.emp.n2,                    # Area above the curve 
    n.PF = nrow(par.res.n2)                # Number of data points in the PF
  ))
  
}

mean(null.df2$a.emp.n >= a.emp) # p-value 0.720

###### Phosphorous v Temperature ######

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol == "P", 'phos', 'other')) # for testing.

df.filt <- df %>% 
  mutate(
    z.y = P.comp.0.56,
    z.y.l = P.comp.0.56.L,
    z.y.u = P.comp.0.56.U,
    z.x = T.br.0.56,
    z.x.l = T.br.0.56.L,
    z.x.u = T.br.0.56.U
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
    dist.sc = distance/mean(distance)
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

df.filt <- df.filt %>% # Trim dramatic outliers as errors
  filter(dist.sc < 2.4) 

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale if necessary

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

par.res.1 <- par.res.1 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.1 <- rbind(                               # Build a polygon enclosing the space above the pareto front
  cbind(par.res.1$z.x, par.res.1$z.y),             # along the frontier
  c(max(par.res.1$z.x)*3, min(par.res.1$z.y)),     # over at the right end
  c(max(par.res.1$z.x)*3, max(par.res.1$z.y)*3),   # up to the top
  c(min(par.res.1$z.x), max(par.res.1$z.y)*3),     # over to the left end
  c(par.res.1$z.x[1], par.res.1$z.y[1])            # close it
)

par.res.2 <- df.filt %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt$z.x + (df.filt$z.x.u - df.filt$z.x)/3, df.filt$z.y + (df.filt$z.y.u - df.filt$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.1[,1], poly.par.1[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Fit a scam to the adjusted PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

par.res.3 <- par.res.3 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.2 <- rbind(                                   # Build a polygon enclosing the space above the pareto front
  cbind(par.res.3$z.x, par.res.3$z.y),                 # along the frontier
  c(max(par.res.3$z.x)*3, min(par.res.3$z.y)),         # over at the right end
  c(max(par.res.3$z.x)*3, max(par.res.3$z.y)*3),       # up to the top
  c(min(par.res.3$z.x), max(par.res.3$z.y)*3),         # over to the left end
  c(par.res.3$z.x[1], par.res.3$z.y[1])                # close it
)

par.res.4 <- df.filt2 %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt2$z.x + (df.filt2$z.x.u - df.filt2$z.x)/3, df.filt2$z.y + (df.filt2$z.y.u - df.filt2$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.2[,1], poly.par.2[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.4))), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

PT.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Thermal breadth  (°C)",    
       y = "Competitive ability (1/P*)", 
       color = "Evolutionary History",
       title = "I — Phosphorous ~ Temperature") +  # labels
  
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
               "phos" = 16)  # filled circle
  ) +
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

PT.scam.PF  # Display the plot

###### 75th quantile PF testing ######

n.75 <- df.filt %>% # Now calculate the number of light and phosphorous points above the 75th % quantile Pareto front. 
  filter(evol.bin == "phos") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n.75 # There is something wrong here....

null_counts <- numeric(1000)

for (i in 1:1000) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  shuffled.df <- df.filt %>%
    
    mutate(
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  null_counts[i] <- shuffled.df %>%
    filter(evol.bin %in% c("phos")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x.sim)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y.sim > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.785

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # Scale for point exclusion determination 
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # re-scale before we fit the pareto front (to remove skew caused by drastic outliers)
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.870

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(      # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.2 %>%
    
    mutate(
      n.x = sample(nrow(par.res.2), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(par.res.2),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(par.res.2), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(par.res.2),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  par.res.n2 <- par_frt(shuffled.df2, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n2 <- par.res.n2 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n2 <- max(shuffled.df2$z.x.sim) # Extract the max values for x and y
  y.max.n2 <- max(shuffled.df2$z.y.sim)
  
  poly.n2 <- par.res.n2[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n2, z.y.sim = y.max.n2)
  
  a.emp.n2 <- polyarea(poly.n2$z.x.sim, poly.n2$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df2 <- rbind(null.df2, data.frame(  # Save the data
    a.emp.n = a.emp.n2,                    # Area above the curve 
    n.PF = nrow(par.res.n2)                # Number of data points in the PF
  ))
  
}

mean(null.df2$a.emp.n >= a.emp) # p-value 0.720

# Salt comparisons --------------------------------------------------------

###### Salt v Temperature ######

df$evol.bin <- ifelse(df$Evol == "none", 'ancestral', 
                      ifelse(df$Evol %in% c("S", "BS"), 'salt', 'other')) # for testing.

df.filt <- df %>% 
  mutate(
    z.y = S.c,
    z.y.l = S.c.L,
    z.y.u = S.c.U,
    z.x = T.br.0.56,
    z.x.l = T.br.0.56.L,
    z.x.u = T.br.0.56.U
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
    dist.sc = distance/mean(distance)
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

df.filt <- df.filt %>% # Trim dramatic outliers as errors
  filter(dist.sc < 2.4) 

df.filt <- df.filt %>% 
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  ) # re-scale if necessary

x.ref <- min(df.filt$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt$z.y2, na.rm = TRUE) # Min y

df.filt <- df.filt %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = distance/mean(distance)
  ) %>%
  arrange(distance) # Recalculate after removing errors

par.res.1 <- par_frt(df.filt[df.filt$dist.sc < 2.1,], xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

par.res.1 <- par.res.1 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.1 <- rbind(                               # Build a polygon enclosing the space above the pareto front
  cbind(par.res.1$z.x, par.res.1$z.y),             # along the frontier
  c(max(par.res.1$z.x)*3, min(par.res.1$z.y)),     # over at the right end
  c(max(par.res.1$z.x)*3, max(par.res.1$z.y)*3),   # up to the top
  c(min(par.res.1$z.x), max(par.res.1$z.y)*3),     # over to the left end
  c(par.res.1$z.x[1], par.res.1$z.y[1])            # close it
)

par.res.2 <- df.filt %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt$z.x + (df.filt$z.x.u - df.filt$z.x)/3, df.filt$z.y + (df.filt$z.y.u - df.filt$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.1[,1], poly.par.1[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.2))), data = par.res.2) # Fit a scam to the adjusted PF

x.vals <- seq(min(df.filt$z.x), max(df.filt$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt %>% # Filter out based on Euclidean distance from min
  arrange(distance) %>%
  slice(1:floor(0.75 * n())) %>%  # keep the closest 75%
  select(-distance)

par.res.3 <- par_frt(df.filt2, xvar = "z.x", yvar = "z.y") # Pareto frontier on this data

par.res.3 <- par.res.3 %>%                  # Organize the order of points for polygon construction
  arrange(z.x, desc(z.y)) 

poly.par.2 <- rbind(                                   # Build a polygon enclosing the space above the pareto front
  cbind(par.res.3$z.x, par.res.3$z.y),                 # along the frontier
  c(max(par.res.3$z.x)*3, min(par.res.3$z.y)),         # over at the right end
  c(max(par.res.3$z.x)*3, max(par.res.3$z.y)*3),       # up to the top
  c(min(par.res.3$z.x), max(par.res.3$z.y)*3),         # over to the left end
  c(par.res.3$z.x[1], par.res.3$z.y[1])                # close it
)

par.res.4 <- df.filt2 %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt2$z.x + (df.filt2$z.x.u - df.filt2$z.x)/3, df.filt2$z.y + (df.filt2$z.y.u - df.filt2$z.y)/3,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.2[,1], poly.par.2[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.4))), data = par.res.4) # Model fit

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

ST.scam.PF <- ggplot(df.filt, aes(x = z.x, y = z.y, color = Evol.plt, shape = evol.bin)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 3, stroke = 1.5) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  
  labs(x = "Thermal breadth  (°C)",    
       y = "Salt tolerance (c)", 
       color = "Evolutionary History",
       title = "J — Salt ~ Temperature") +  # labels
  
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
  
  #ylim(0.04,0.21) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

ST.scam.PF  # Display the plot

###### 75th quantile PF testing ######

n.75 <- df.filt %>% # Now calculate the number of light and phosphorous points above the 75th % quantile Pareto front. 
  filter(evol.bin == "salt") %>%
  rowwise() %>%
  mutate(
    nearest_idx = find_nearest_index(z.x, pred.curve.2$z.x),
    z.y.pred = pred.curve.2$z.y[nearest_idx]
  ) %>%
  ungroup() %>%
  filter(z.y > z.y.pred) %>%
  nrow()

n.75 # There is something wrong here....

null_counts <- numeric(1000)

for (i in 1:1000) { # Now we'll randomize and see how many light points should fall above the 75th qr by chance alone
  shuffled.df <- df.filt %>%
    
    mutate(
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  null_counts[i] <- shuffled.df %>%
    filter(evol.bin %in% c("salt")) %>%
    rowwise() %>%
    mutate(
      nearest_idx = which.min(abs(pred.curve.2$z.x - z.x.sim)),
      z.y.pred = pred.curve.2$z.y[nearest_idx]
    ) %>%
    ungroup() %>%
    filter(z.y.sim > z.y.pred) %>%
    nrow()
}

mean(null_counts>= n.75) # p = 0.785

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.x) # Extract the max values for x and y. For now threshold at the same 2.1 level?
y.max <- max(df.filt[df.filt$dist.sc < 2.1,]$z.y)

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
      n.x = sample(nrow(df.filt), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(df.filt),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(df.filt), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(df.filt),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # Scale for point exclusion determination 
  
  shuffled.df <- shuffled.df %>% # Calculate Euclidean distance from min
    mutate(
      distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
      dist.sc = distance/mean(distance)
    ) %>%
    filter(dist.sc < 2.4) # Error trimming
  
  shuffled.df <- shuffled.df %>% 
    mutate(
      z.y2 = scale(z.y.sim)[, 1],
      z.x2 = scale(z.x.sim)[, 1]
    ) # re-scale before we fit the pareto front (to remove skew caused by drastic outliers)
  
  par.res.n <- par_frt(shuffled.df[shuffled.df$dist.sc < 2.1, ], xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df[shuffled.df$dist.sc < 2.1,]$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df <- rbind(null.df, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df$a.emp.n >= a.emp) # p-value 0.870

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(      # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.2 %>%
    
    mutate(
      n.x = sample(nrow(par.res.2), replace = FALSE),     # Randomly assign row
      z.x.sim = runif(
        nrow(par.res.2),                                  
        z.x[n.x] - (z.x[n.x] - z.x.l[n.x]) / 2,         # Lower CI cut in half
        z.x[n.x] + (z.x.u[n.x] - z.x[n.x]) / 2          # Upper CI cut in half
      ),
      
      n.y = sample(nrow(par.res.2), replace = FALSE),     # Seperately reassign row for y
      z.y.sim = runif(
        nrow(par.res.2),                                  
        z.y[n.y] - (z.y[n.y] - z.y.l[n.y]) / 2,         # Lower CI cut in half
        z.y[n.y] + (z.y.u[n.y] - z.y[n.y]) / 2          # Upper CI cut in half
      )
    )
  
  par.res.n2 <- par_frt(shuffled.df2, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n2 <- par.res.n2 %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n2 <- max(shuffled.df2$z.x.sim) # Extract the max values for x and y
  y.max.n2 <- max(shuffled.df2$z.y.sim)
  
  poly.n2 <- par.res.n2[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n2, z.y.sim = y.max.n2)
  
  a.emp.n2 <- polyarea(poly.n2$z.x.sim, poly.n2$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df2 <- rbind(null.df2, data.frame(  # Save the data
    a.emp.n = a.emp.n2,                    # Area above the curve 
    n.PF = nrow(par.res.n2)                # Number of data points in the PF
  ))
  
}

mean(null.df2$a.emp.n >= a.emp) # p-value 0.720

# Compile & save the figure -----------------------------------------------

legend_df.2 <- data.frame(
  x = c(1, 2, 1, 2, 1, 2, 1, 2),
  y = c(1, 1, 2, 2, 1, 1, 2, 2),
  Group = factor(c("Ancestral", "Other", "Matching", "Matching", "Ancestral", "Other", "Matching", "Matching")),
  Group2 = factor(c("Biotic depletion", "Biotic depletion x Salt", "Control", "Light limitation", "Nitrogen limitation", "Ancestral", "Phosphorous limitation", "Salt stress")),
  LineType = factor(c("Pareto Front", "Pareto Front", "Pareto Front", "Pareto Front", "Pareto Front", "Pareto Front", "Pareto Front", "Pareto Front"))
)

legend_plot.2 <- ggplot(legend_df.2, aes(x = x, y = y)) +
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
  
  scale_linetype_manual(name = NULL, values = c("Pareto Front" = "solid")) +
  
  labs(color = "Evolutionary context") +
  theme_void() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "lines")
  )

legend_only.2 <- get_legend(legend_plot.2)

full_toffs <- plot_grid(
  LN.scam.PF, LP.scam.PF, LS.scam.PF, LT.scam.PF,
  NULL, NP.scam.PF, NS.scam.PF, NT.scam.PF,
  NULL, NULL, PS.scam.PF, PT.scam.PF,
  legend_only.2, NULL, NULL, ST.scam.PF,
  ncol = 4,
  align = "hv",
  axis = "tblr"
)

ggsave("figures/14_fig_5_inter-gradient_toffs.jpeg", full_toffs, width = 15, height = 15)
