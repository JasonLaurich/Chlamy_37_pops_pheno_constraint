# Jason R Laurich

# September 6th, 2025

# Creating Figure 3: Intra-gradient Pareto front panels (light, nitrogen, phosphorous, salt, and temperature)
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

# Temperature -------------------------------------------------------------

df$evol.bin <- ifelse(df$Evol == "none", "ancestral", "evolved") # For binning into evolutionary treatments for plotting purposes.

df.filt <- df %>% 
  mutate(
    z.y = T.br.0.56,
    z.y.l = T.br.0.56.L,
    z.y.u = T.br.0.56.U,
    z.x = T.µ.max,
    z.x.l = T.µ.max.L,
    z.x.u = T.µ.max.U
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

poly.par.1 <- rbind(                        # Build a polygon enclosing the space above the pareto front
  cbind(par.res.1$z.x, par.res.1$z.y),      # along the frontier
  c(10, min(par.res.1$z.y)),                # over at the right end
  c(10, 40),                                # up to the top
  c(min(par.res.1$z.x), 40),                # over to the left end
  c(par.res.1$z.x[1], par.res.1$z.y[1])     # close it
)

par.res.2 <- df.filt %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt$z.x.u - (df.filt$z.x.u - df.filt$z.x)/2, df.filt$z.y.u - (df.filt$z.y.u - df.filt$z.y)/2,            # Create a column that identifies points whose upper CIs fall in this polygon
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

poly.par.2 <- rbind(                        # Build a polygon enclosing the space above the pareto front
  cbind(par.res.3$z.x, par.res.3$z.y),      # along the frontier
  c(10, min(par.res.3$z.y)),                # over at the right end
  c(10, 40),                                # up to the top
  c(min(par.res.3$z.x), 40),                # over to the left end
  c(par.res.3$z.x[1], par.res.3$z.y[1])     # close it
)

par.res.4 <- df.filt2 %>%
  mutate(pf.inc = sp::point.in.polygon(df.filt2$z.x.u - (df.filt2$z.x.u - df.filt2$z.x)/2, df.filt2$z.y.u - (df.filt2$z.y.u - df.filt2$z.y)/2,            # Create a column that identifies points whose upper CIs fall in this polygon
                                       poly.par.2[,1], poly.par.2[,2])) %>% 
  filter(pf.inc == 1)                                                           # Filter based on results

fit2 <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.4))), data = par.res.4) # Model fit

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
  
  ylim(24.5,29.5) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.scam  # Display the plot

###### 75th quantile PF testing ######

# Cannot do for temperature data, as no populations evolved under temperature stress

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

mean(null.df$a.emp.n >= a.emp) # p-value 0.643

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(       # Null model results
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

mean(null.df2$a.emp.n >= a.emp) # p-value 0.121
    
  