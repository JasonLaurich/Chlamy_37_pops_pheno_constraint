# Jason R Laurich

# January 21st, 2025

# We are going to analyze Shoval and Pareto limitation across gradients for our synthesis data. Will also perform
# intra-gradient comparisons, treating sp. as the level (e.g. excluding unknown entities and collapsing replicates)
# and strains to a single mean).

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(pracma) # For calculating polygon area
library(cowplot)

library(quantreg)
library(sp) # For the point.in.polygon function
library(scam)

par_frt <- function(df, xvar, yvar) { # Simple Pareto front function / convex hull algorithm (one sided)
  
  df <- df[order(-df[[xvar]], df[[yvar]]), ]  
  pareto_points <- df[1, ]  # Start with the first point
  
  for (i in 2:nrow(df)) {
    if (df[i, yvar] >= tail(pareto_points[[yvar]], 1)) {  # Ensure increasing y values
      pareto_points <- rbind(pareto_points, df[i,])
    }
  }
  
  return(pareto_points)
}

find_nearest_index <- function(x, ref_vec) { # Function to find the closest point to the actual r.max value in the pred data frame. 
  which.min(abs(ref_vec - x))
} # For significance testing based on location of points relative to a curve

# Load and organize the data ----------------------------------------------

# Light

df.l <- read.csv("data-processed/504a_light_metadata.csv")
head(df.l)

Encoding(df.l$Sp.name) <- "bytes" # make sure R treats the strings as bytes (prevents “invalid UTF-8” headaches)

nbsp.byte <- rawToChar(as.raw(0xA0)) # build a literal 0xA0 byte and replace it with a normal space
df.l$Sp.name <- gsub(nbsp.byte, " ", df.l$Sp.name, fixed = TRUE, useBytes = TRUE)

df.l$Sp.name <- enc2utf8(df.l$Sp.name) # now safely convert to UTF-8 for downstream joins
df.l$Sp.name <- trimws(df.l$Sp.name)

unique(df.l$Sp.name)

nrow(df.l %>% 
  filter(Sp.name!="NA")) # 112 species. 

# Nitrogen

df.n <- read.csv("data-processed/504b_nit_metadata.csv")
head(df.n)

Encoding(df.n$Sp.name) <- "bytes" # make sure R treats the strings as bytes (prevents “invalid UTF-8” headaches)

df.n$Sp.name <- gsub(nbsp.byte, " ", df.n$Sp.name, fixed = TRUE, useBytes = TRUE)

df.n$Sp.name <- enc2utf8(df.n$Sp.name) # now safely convert to UTF-8 for downstream joins
df.n$Sp.name <- trimws(df.n$Sp.name)

unique(df.n$Sp.name)

nrow(df.n %>% 
       filter(Sp.name!="NA")) # 43 species. 

# Phosphorous 

df.p <- read.csv("data-processed/504c_phos_metadata.csv")
head(df.p)

Encoding(df.p$Sp.name) <- "bytes" # make sure R treats the strings as bytes (prevents “invalid UTF-8” headaches)

df.p$Sp.name <- gsub(nbsp.byte, " ", df.p$Sp.name, fixed = TRUE, useBytes = TRUE)

df.p$Sp.name <- enc2utf8(df.p$Sp.name) # now safely convert to UTF-8 for downstream joins
df.p$Sp.name <- trimws(df.p$Sp.name)

unique(df.p$Sp.name)

nrow(df.p %>% 
       filter(Sp.name!="NA")) # 91 species.

# Temperature

df.t <- read.csv("data-processed/504d_temp_metadata.csv")
head(df.t)

Encoding(df.t$Sp.name) <- "bytes" # make sure R treats the strings as bytes (prevents “invalid UTF-8” headaches)

df.t$Sp.name <- gsub(nbsp.byte, " ", df.t$Sp.name, fixed = TRUE, useBytes = TRUE)

df.t$Sp.name <- enc2utf8(df.t$Sp.name) # now safely convert to UTF-8 for downstream joins
df.t$Sp.name <- trimws(df.t$Sp.name)

unique(df.t$Sp.name)

nrow(df.t %>% 
       filter(Sp.name!="NA")) # 197 species.

# OK so we have r.max and either T.br or comp for all. We need to combine these into a single data frame. 
# Combine the data frames

df.t.1 <- df.t %>%
  group_by(Sp.name) %>%
  summarise(
    r.max.T = mean(r.max, na.rm = TRUE),
    T.br    = mean(T.br,  na.rm = TRUE),
    .groups = "drop"
  )

df.n.1 <- df.n %>%
  group_by(Sp.name) %>%
  summarise(
    r.max.N = mean(r.max, na.rm = TRUE),
    comp.N  = mean(comp,  na.rm = TRUE),
    .groups = "drop"
  )

df.p.1 <- df.p %>%
  group_by(Sp.name) %>%
  summarise(
    r.max.P = mean(r.max, na.rm = TRUE),
    comp.P  = mean(comp,  na.rm = TRUE),
    .groups = "drop"
  )

df.l.1 <- df.l %>%
  group_by(Sp.name) %>%
  summarise(
    r.max.L = mean(r.max, na.rm = TRUE),
    comp.L  = mean(comp,  na.rm = TRUE),
    .groups = "drop"
  )

df <- df.t.1 %>%
  full_join(df.n.1, by = "Sp.name") %>%
  full_join(df.p.1, by = "Sp.name") %>%
  full_join(df.l.1, by = "Sp.name")

write.csv(df, "data-processed/505_synthesis_metadata_sp_summary.csv") # Summary

# Analyses ----------------------------------------------------------------

# Light v Light -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.L,
    z.x = r.max.L
  ) %>% # Specify the x and y variables
  filter(!is.na(z.y)) %>% 
  filter(!is.na(z.x)) # 99 observations

plot(z.y ~ z.x, data= df.filt) # look at the data

###### Shoval et al triangular analysis ######

hull.id <- chull(df.filt$z.x, df.filt$z.y)

df.hull.full <- df.filt[hull.id,] %>% 
  bind_rows(df.filt[hull.id[1],]) # Close the polygon

ggplot(df.filt, aes(z.x, z.y)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(z.x, z.y, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Find the smallest triangle ######

df.hull <- df.filt[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri <- data.frame(
  anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
  
  apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
  baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
  baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
  
  area = numeric(),                                                             # Triangle area
  
  poly.area = numeric(),                                                        # Polygon size for comparison
  
  stringsAsFactors = FALSE
)


for (i in 1:n){
  
  a <- i # select the vertices to cut open the polygon on
  b <- ifelse(i<n, (i + 1), 1)
  
  idx <- a
  k <- a
  
  while (k != b) {
    k <- ifelse(k == 1, n, k - 1)
    idx <- c(idx, k)
  }
  
  df.cut <- df.hull[idx, ] # OK our loop has succesfully cut open the polygon along the edge between 1 and 2!
  
  # Step 2 - find the antipode E' of E
  
  A <- as.numeric(df.cut[1,c("z.x", "z.y")])                      # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("z.x", "z.y")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x", "z.y")]) # Points to test (calculating distance from E)
  
  AB <- B-A 
  num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
  dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
  
  k <- which.max(dist) # Identify antipode (max distance vertex)
  
  E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
  
  # Split the polygon into left and right sides
  left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
  right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
  
  # Fit data frame of these as edge sequences
  L.edges <- data.frame(
    x1 = left.verts$z.x[-nrow(left.verts)],
    y1 = left.verts$z.y[-nrow(left.verts)],
    x2 = left.verts$z.x[-1],
    y2 = left.verts$z.y[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$z.x[-nrow(right.verts)],
    y1 = right.verts$z.y[-nrow(right.verts)],
    x2 = right.verts$z.x[-1],
    y2 = right.verts$z.y[-1]
  )
  
  t.hat <- AB/ sqrt(sum(AB^2))
  n.hat <- c(-t.hat[2], t.hat[1])
  
  cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
  
  best.area <- Inf
  best <- NULL
  
  for (l in 1:nrow(L.edges)) {
    k <- L.edges[l,]
    
    p <- as.numeric(c(k[1], k[2]))
    q <- as.numeric(c(k[3], k[4]))
    d.L <- q - p
    d.L <- d.L / sqrt(sum(d.L^2))
    
    nL <- c(-d.L[2], d.L[1])
    vals <- df.hull$z.x * nL[1] + df.hull$z.y * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$z.x[idx], df.hull$z.y[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$z.x * nR[1] + df.hull$z.y * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$z.x[idx], df.hull$z.y[idx])
      
      M <- cbind(d.L, -d.R) # calculate the apex
      if (abs(det(M)) < 1e-12) next  # skip near-parallel
      
      sol <- solve(M, x0.R - x0.L)
      apex <- x0.L + sol[1] * d.L
      
      apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
      if (apex.h <= cE) next
      
      M1 <- cbind(AB, -d.L)
      sol1 <- solve(M1, x0.L - A)
      base.L <- A + sol1[1] * AB
      
      M2 <- cbind(AB, -d.R)
      sol2 <- solve(M2, x0.R - A)
      base.R <- A + sol2[1] * AB
      
      v1 <- base.L - apex
      v2 <- base.R - apex
      area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
      
      if (area < best.area) {
        best.area <- area
        best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
      }
    }
  }
  
  tri <- rbind(tri,data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$z.x, df.hull.full$z.y))               # Polygon area
  ))
  
}

best.tri <- tri %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly <- data.frame(
  x = c(
    best.tri$apex.x,
    best.tri$baseL.x,
    best.tri$baseR.x,
    best.tri$apex.x   # close polygon
  ),
  y = c(
    best.tri$apex.y,
    best.tri$baseL.y,
    best.tri$baseR.y,
    best.tri$apex.y
  )
)

p.l <- ggplot(df.filt, aes(z.x, z.y)) +
  geom_point(size = 2) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  geom_polygon(
    data = tri.poly,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  ylim(-0.10, 3.25) +
  xlim(-0.30, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.l

tri

best.tri$area/best.tri$poly.area # 1.11781

###### Randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.filt %>%          # scramble the data
    
    mutate(
      z.x = sample(df.filt$z.x, replace = FALSE),             # Randomly assign rmax
      z.y = sample(df.filt$z.y, replace = FALSE)              # Separately reassign z.y
    )
  
  hull.id.null <- chull(df.null$z.x, df.null$z.y)    # Get the polygon
  
  df.hull.full.null <- df.null[hull.id.null,] %>% 
    bind_rows(df.null[hull.id.null[1],]) # Close the polygon
  
  df.hull.null <- df.null[hull.id.null,]  # Don't close the polygon
  
  n <- nrow(df.hull.null) # need the number of vertices in the convex hull dataset
  
  df.hull.null <- df.hull.null[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 
  
  tri.null <- data.frame(
    anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
    
    apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
    baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
    baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
    
    area = numeric(),                                                             # Triangle area
    
    poly.area = numeric(),                                                        # Polygon size for z.yarison
    
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n){
    
    a <- i # select the vertices to cut open the polygon on
    b <- ifelse(i<n, (i + 1), 1)
    
    idx <- a
    k <- a
    
    while (k != b) {
      k <- ifelse(k == 1, n, k - 1)
      idx <- c(idx, k)
    }
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("z.x", "z.y")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("z.x","z.y")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x","z.y")]) # Points to test (calculating distance from E)
    
    AB <- B-A 
    num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
    dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
    
    k <- which.max(dist) # Identify antipode (max distance vertex)
    
    E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
    
    # Split the polygon into left and right sides
    left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
    right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
    
    # Fit data frame of these as edge sequences
    L.edges <- data.frame(
      x1 = left.verts$z.x[-nrow(left.verts)],
      y1 = left.verts$z.y[-nrow(left.verts)],
      x2 = left.verts$z.x[-1],
      y2 = left.verts$z.y[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$z.x[-nrow(right.verts)],
      y1 = right.verts$z.y[-nrow(right.verts)],
      x2 = right.verts$z.x[-1],
      y2 = right.verts$z.y[-1]
    )
    
    t.hat <- AB/ sqrt(sum(AB^2))
    n.hat <- c(-t.hat[2], t.hat[1])
    
    cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
    
    best.area <- Inf
    best <- NULL
    
    for (l in 1:nrow(L.edges)) {
      k <- L.edges[l,]
      
      p <- as.numeric(c(k[1], k[2]))
      q <- as.numeric(c(k[3], k[4]))
      d.L <- q - p
      d.L <- d.L / sqrt(sum(d.L^2))
      
      nL <- c(-d.L[2], d.L[1])
      vals <- df.hull.null$z.x * nL[1] + df.hull.null$z.y * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$z.x * nR[1] + df.hull.null$z.y * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
        
        M <- cbind(d.L, -d.R) # calculate the apex
        if (abs(det(M)) < 1e-12) next  # skip near-parallel
        
        sol <- solve(M, x0.R - x0.L)
        apex <- x0.L + sol[1] * d.L
        
        apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
        if (apex.h <= cE) next
        
        M1 <- cbind(AB, -d.L)
        sol1 <- solve(M1, x0.L - A)
        base.L <- A + sol1[1] * AB
        
        M2 <- cbind(AB, -d.R)
        sol2 <- solve(M2, x0.R - A)
        base.R <- A + sol2[1] * AB
        
        v1 <- base.L - apex
        v2 <- base.R - apex
        area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
        
        if (area < best.area) {
          best.area <- area
          best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
        }
      }
    }
    
    if (is.null(best) || !is.finite(best.area)) next
    
    tri.null <- rbind(tri.null,data.frame(
      anchor.edge = i,                                                            # Which segment is E (leading vertex)
      
      apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
      baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
      baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
      
      area = best.area,                                                           # Triangle area
      poly.area = abs(polyarea(df.hull.full.null$z.x, df.hull.full.null$z.y))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri$poly.area[1]) # 0.017

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.057

###### Pareto fronts ######

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

I.scam <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  ylim(-0.10, 3.25) +
  xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

I.scam  # Display the plot

p.l3 <- p.l + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE)  # Adding scam PF fits
p.l3

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

mean(null.df$a.emp.n >= a.emp) # p-value 0.273

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 0.09410, p 0.43441
summary(q75, se = "boot", R = 1000) # -0.03627, p 0.79101
summary(q90, se = "boot", R = 1000) # 0.32783, p 0.31632 

I.qr <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  ylim(-0.10, 3.25) +
  xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

I.qr # Display the plot

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%                           # smallest → largest
  slice((floor(0.6667 * n()) + 1):n()) %>%          # keep *second* half
  select(-distance)

I.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  ylim(-0.10, 3.25) +
  xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

I.scam2  # Display the plot

null.df3 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df3 <- rbind(null.df3, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 0

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -0.78720, p 0.06861
summary(q75, se = "boot", R = 1000) # -1.37943, p 0.00730
summary(q90, se = "boot", R = 1000) # -1.41332, p 0.00597

I.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  ylim(-0.10, 3.25) +
  xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

I.qr2 # Display the plot

# Light v Nitrogen -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.L,
    z.x = comp.N
  ) %>% # Specify the x and y variables
  filter(!is.na(z.y)) %>% 
  filter(!is.na(z.x)) # Only 20 observations

plot(z.y ~ z.x, data= df.filt) # look at the data 
# A really high data point in the upper right

###### Shoval et al triangular analysis ######

hull.id <- chull(df.filt$z.x, df.filt$z.y)

df.hull.full <- df.filt[hull.id,] %>% 
  bind_rows(df.filt[hull.id[1],]) # Close the polygon

ggplot(df.filt, aes(z.x, z.y)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(z.x, z.y, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Find the smallest triangle ######

df.hull <- df.filt[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri <- data.frame(
  anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
  
  apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
  baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
  baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
  
  area = numeric(),                                                             # Triangle area
  
  poly.area = numeric(),                                                        # Polygon size for comparison
  
  stringsAsFactors = FALSE
)


for (i in 1:n){
  
  a <- i # select the vertices to cut open the polygon on
  b <- ifelse(i<n, (i + 1), 1)
  
  idx <- a
  k <- a
  
  while (k != b) {
    k <- ifelse(k == 1, n, k - 1)
    idx <- c(idx, k)
  }
  
  df.cut <- df.hull[idx, ] # OK our loop has succesfully cut open the polygon along the edge between 1 and 2!
  
  # Step 2 - find the antipode E' of E
  
  A <- as.numeric(df.cut[1,c("z.x", "z.y")])                      # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("z.x", "z.y")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x", "z.y")]) # Points to test (calculating distance from E)
  
  AB <- B-A 
  num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
  dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
  
  k <- which.max(dist) # Identify antipode (max distance vertex)
  
  E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
  
  # Split the polygon into left and right sides
  left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
  right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
  
  # Fit data frame of these as edge sequences
  L.edges <- data.frame(
    x1 = left.verts$z.x[-nrow(left.verts)],
    y1 = left.verts$z.y[-nrow(left.verts)],
    x2 = left.verts$z.x[-1],
    y2 = left.verts$z.y[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$z.x[-nrow(right.verts)],
    y1 = right.verts$z.y[-nrow(right.verts)],
    x2 = right.verts$z.x[-1],
    y2 = right.verts$z.y[-1]
  )
  
  t.hat <- AB/ sqrt(sum(AB^2))
  n.hat <- c(-t.hat[2], t.hat[1])
  
  cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
  
  best.area <- Inf
  best <- NULL
  
  for (l in 1:nrow(L.edges)) {
    k <- L.edges[l,]
    
    p <- as.numeric(c(k[1], k[2]))
    q <- as.numeric(c(k[3], k[4]))
    d.L <- q - p
    d.L <- d.L / sqrt(sum(d.L^2))
    
    nL <- c(-d.L[2], d.L[1])
    vals <- df.hull$z.x * nL[1] + df.hull$z.y * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$z.x[idx], df.hull$z.y[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$z.x * nR[1] + df.hull$z.y * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$z.x[idx], df.hull$z.y[idx])
      
      M <- cbind(d.L, -d.R) # calculate the apex
      if (abs(det(M)) < 1e-12) next  # skip near-parallel
      
      sol <- solve(M, x0.R - x0.L)
      apex <- x0.L + sol[1] * d.L
      
      apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
      if (apex.h <= cE) next
      
      M1 <- cbind(AB, -d.L)
      sol1 <- solve(M1, x0.L - A)
      base.L <- A + sol1[1] * AB
      
      M2 <- cbind(AB, -d.R)
      sol2 <- solve(M2, x0.R - A)
      base.R <- A + sol2[1] * AB
      
      v1 <- base.L - apex
      v2 <- base.R - apex
      area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
      
      if (area < best.area) {
        best.area <- area
        best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
      }
    }
  }
  
  tri <- rbind(tri,data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$z.x, df.hull.full$z.y))               # Polygon area
  ))
  
}

best.tri <- tri %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly <- data.frame(
  x = c(
    best.tri$apex.x,
    best.tri$baseL.x,
    best.tri$baseR.x,
    best.tri$apex.x   # close polygon
  ),
  y = c(
    best.tri$apex.y,
    best.tri$baseL.y,
    best.tri$baseR.y,
    best.tri$apex.y
  )
)

p.ln <- ggplot(df.filt, aes(z.x, z.y)) +
  geom_point(size = 2) +
  
  labs(x = "Competitive ability (1/N*)",    
       y = "Competitive ability (1/I*)", 
       title = "B — Light ~ Nitrogen") +  # labels
  
  geom_polygon(
    data = tri.poly,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  #ylim(-0.10, 3.25) +
  #xlim(-0.30, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.ln

tri

best.tri$area/best.tri$poly.area # 2.106811

###### Randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.filt %>%          # scramble the data
    
    mutate(
      z.x = sample(df.filt$z.x, replace = FALSE),             # Randomly assign rmax
      z.y = sample(df.filt$z.y, replace = FALSE)              # Separately reassign z.y
    )
  
  hull.id.null <- chull(df.null$z.x, df.null$z.y)    # Get the polygon
  
  df.hull.full.null <- df.null[hull.id.null,] %>% 
    bind_rows(df.null[hull.id.null[1],]) # Close the polygon
  
  df.hull.null <- df.null[hull.id.null,]  # Don't close the polygon
  
  n <- nrow(df.hull.null) # need the number of vertices in the convex hull dataset
  
  df.hull.null <- df.hull.null[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 
  
  tri.null <- data.frame(
    anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
    
    apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
    baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
    baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
    
    area = numeric(),                                                             # Triangle area
    
    poly.area = numeric(),                                                        # Polygon size for z.yarison
    
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n){
    
    a <- i # select the vertices to cut open the polygon on
    b <- ifelse(i<n, (i + 1), 1)
    
    idx <- a
    k <- a
    
    while (k != b) {
      k <- ifelse(k == 1, n, k - 1)
      idx <- c(idx, k)
    }
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("z.x", "z.y")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("z.x","z.y")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x","z.y")]) # Points to test (calculating distance from E)
    
    AB <- B-A 
    num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
    dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
    
    k <- which.max(dist) # Identify antipode (max distance vertex)
    
    E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
    
    # Split the polygon into left and right sides
    left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
    right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
    
    # Fit data frame of these as edge sequences
    L.edges <- data.frame(
      x1 = left.verts$z.x[-nrow(left.verts)],
      y1 = left.verts$z.y[-nrow(left.verts)],
      x2 = left.verts$z.x[-1],
      y2 = left.verts$z.y[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$z.x[-nrow(right.verts)],
      y1 = right.verts$z.y[-nrow(right.verts)],
      x2 = right.verts$z.x[-1],
      y2 = right.verts$z.y[-1]
    )
    
    t.hat <- AB/ sqrt(sum(AB^2))
    n.hat <- c(-t.hat[2], t.hat[1])
    
    cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
    
    best.area <- Inf
    best <- NULL
    
    for (l in 1:nrow(L.edges)) {
      k <- L.edges[l,]
      
      p <- as.numeric(c(k[1], k[2]))
      q <- as.numeric(c(k[3], k[4]))
      d.L <- q - p
      d.L <- d.L / sqrt(sum(d.L^2))
      
      nL <- c(-d.L[2], d.L[1])
      vals <- df.hull.null$z.x * nL[1] + df.hull.null$z.y * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$z.x * nR[1] + df.hull.null$z.y * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
        
        M <- cbind(d.L, -d.R) # calculate the apex
        if (abs(det(M)) < 1e-12) next  # skip near-parallel
        
        sol <- solve(M, x0.R - x0.L)
        apex <- x0.L + sol[1] * d.L
        
        apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
        if (apex.h <= cE) next
        
        M1 <- cbind(AB, -d.L)
        sol1 <- solve(M1, x0.L - A)
        base.L <- A + sol1[1] * AB
        
        M2 <- cbind(AB, -d.R)
        sol2 <- solve(M2, x0.R - A)
        base.R <- A + sol2[1] * AB
        
        v1 <- base.L - apex
        v2 <- base.R - apex
        area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
        
        if (area < best.area) {
          best.area <- area
          best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
        }
      }
    }
    
    if (is.null(best) || !is.finite(best.area)) next
    
    tri.null <- rbind(tri.null,data.frame(
      anchor.edge = i,                                                            # Which segment is E (leading vertex)
      
      apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
      baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
      baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
      
      area = best.area,                                                           # Triangle area
      poly.area = abs(polyarea(df.hull.full.null$z.x, df.hull.full.null$z.y))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri$poly.area[1]) # 0.689

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.944

###### Pareto fronts ######

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
fit2 <- lm(z.y~z.x, data = par.res.2)

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

LN.scam <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  # geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Competitive ability (1/N*)",    
       y = "Competitive ability (1/I*)", 
       title = "B — Light ~ Nitrogen") +  # labels
  
  ylim(-0.10, 2.25) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LN.scam  # Display the plot

p.ln3 <- p.ln  # Adding scam PF fits
p.ln3

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

mean(null.df$a.emp.n >= a.emp) # p-value 1

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 0.02385, p 0.09510
summary(q75, se = "boot", R = 1000) # 0.02983, p 0.15108
summary(q90, se = "boot", R = 1000) # 0.00593, p 0.85549

LN.qr <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Competitive ability (1/N*)",    
       y = "Competitive ability (1/I*)", 
       title = "B — Light ~ Nitrogen") +  # labels
  
  ylim(-0.10, 2.25) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LN.qr # Display the plot

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%                           # smallest → largest
  slice((floor(0.6667 * n()) + 1):n()) %>%          # keep *second* half
  select(-distance)

LN.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  # geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Competitive ability (1/N*)",    
       y = "Competitive ability (1/I*)", 
       title = "B — Light ~ Nitrogen") +  # labels
  
  ylim(-0.10, 2.25) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LN.scam2  # Display the plot

null.df3 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df3 <- rbind(null.df3, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 1

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -0.03021, p 0.34870
summary(q75, se = "boot", R = 1000) # 0.00348, p 0.91420
summary(q90, se = "boot", R = 1000) # 0.00348, p 0.89562

LN.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Competitive ability (1/N*)",    
       y = "Competitive ability (1/I*)", 
       title = "B — Light ~ Nitrogen") +  # labels
  
  ylim(-0.10, 2.25) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LN.qr2 # Display the plot

# Light v Phosphorous -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.L,
    z.x = comp.P
  ) %>% # Specify the x and y variables
  filter(!is.na(z.y)) %>% 
  filter(!is.na(z.x)) # Only 38 observations

plot(z.y ~ z.x, data= df.filt) # look at the data 
# Another high point in the upper right that will drive relationships here. 

###### Shoval et al triangular analysis ######

hull.id <- chull(df.filt$z.x, df.filt$z.y)

df.hull.full <- df.filt[hull.id,] %>% 
  bind_rows(df.filt[hull.id[1],]) # Close the polygon

ggplot(df.filt, aes(z.x, z.y)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(z.x, z.y, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Find the smallest triangle ######

df.hull <- df.filt[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri <- data.frame(
  anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
  
  apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
  baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
  baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
  
  area = numeric(),                                                             # Triangle area
  
  poly.area = numeric(),                                                        # Polygon size for comparison
  
  stringsAsFactors = FALSE
)


for (i in 1:n){
  
  a <- i # select the vertices to cut open the polygon on
  b <- ifelse(i<n, (i + 1), 1)
  
  idx <- a
  k <- a
  
  while (k != b) {
    k <- ifelse(k == 1, n, k - 1)
    idx <- c(idx, k)
  }
  
  df.cut <- df.hull[idx, ] # OK our loop has succesfully cut open the polygon along the edge between 1 and 2!
  
  # Step 2 - find the antipode E' of E
  
  A <- as.numeric(df.cut[1,c("z.x", "z.y")])                      # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("z.x", "z.y")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x", "z.y")]) # Points to test (calculating distance from E)
  
  AB <- B-A 
  num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
  dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
  
  k <- which.max(dist) # Identify antipode (max distance vertex)
  
  E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
  
  # Split the polygon into left and right sides
  left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
  right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
  
  # Fit data frame of these as edge sequences
  L.edges <- data.frame(
    x1 = left.verts$z.x[-nrow(left.verts)],
    y1 = left.verts$z.y[-nrow(left.verts)],
    x2 = left.verts$z.x[-1],
    y2 = left.verts$z.y[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$z.x[-nrow(right.verts)],
    y1 = right.verts$z.y[-nrow(right.verts)],
    x2 = right.verts$z.x[-1],
    y2 = right.verts$z.y[-1]
  )
  
  t.hat <- AB/ sqrt(sum(AB^2))
  n.hat <- c(-t.hat[2], t.hat[1])
  
  cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
  
  best.area <- Inf
  best <- NULL
  
  for (l in 1:nrow(L.edges)) {
    k <- L.edges[l,]
    
    p <- as.numeric(c(k[1], k[2]))
    q <- as.numeric(c(k[3], k[4]))
    d.L <- q - p
    d.L <- d.L / sqrt(sum(d.L^2))
    
    nL <- c(-d.L[2], d.L[1])
    vals <- df.hull$z.x * nL[1] + df.hull$z.y * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$z.x[idx], df.hull$z.y[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$z.x * nR[1] + df.hull$z.y * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$z.x[idx], df.hull$z.y[idx])
      
      M <- cbind(d.L, -d.R) # calculate the apex
      if (abs(det(M)) < 1e-12) next  # skip near-parallel
      
      sol <- solve(M, x0.R - x0.L)
      apex <- x0.L + sol[1] * d.L
      
      apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
      if (apex.h <= cE) next
      
      M1 <- cbind(AB, -d.L)
      sol1 <- solve(M1, x0.L - A)
      base.L <- A + sol1[1] * AB
      
      M2 <- cbind(AB, -d.R)
      sol2 <- solve(M2, x0.R - A)
      base.R <- A + sol2[1] * AB
      
      v1 <- base.L - apex
      v2 <- base.R - apex
      area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
      
      if (area < best.area) {
        best.area <- area
        best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
      }
    }
  }
  
  tri <- rbind(tri,data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$z.x, df.hull.full$z.y))               # Polygon area
  ))
  
}

best.tri <- tri %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly <- data.frame(
  x = c(
    best.tri$apex.x,
    best.tri$baseL.x,
    best.tri$baseR.x,
    best.tri$apex.x   # close polygon
  ),
  y = c(
    best.tri$apex.y,
    best.tri$baseL.y,
    best.tri$baseR.y,
    best.tri$apex.y
  )
)

p.lp <- ggplot(df.filt, aes(z.x, z.y)) +
  geom_point(size = 2) +
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/I*)", 
       title = "C — Light ~ Phosphorous") +  # labels
  
  geom_polygon(
    data = tri.poly,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  ylim(-0.65, 3.25) +
  # xlim(0, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.lp

tri

best.tri$area/best.tri$poly.area # 1.789155

###### Randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.filt %>%          # scramble the data
    
    mutate(
      z.x = sample(df.filt$z.x, replace = FALSE),             # Randomly assign rmax
      z.y = sample(df.filt$z.y, replace = FALSE)              # Separately reassign z.y
    )
  
  hull.id.null <- chull(df.null$z.x, df.null$z.y)    # Get the polygon
  
  df.hull.full.null <- df.null[hull.id.null,] %>% 
    bind_rows(df.null[hull.id.null[1],]) # Close the polygon
  
  df.hull.null <- df.null[hull.id.null,]  # Don't close the polygon
  
  n <- nrow(df.hull.null) # need the number of vertices in the convex hull dataset
  
  df.hull.null <- df.hull.null[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 
  
  tri.null <- data.frame(
    anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
    
    apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
    baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
    baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
    
    area = numeric(),                                                             # Triangle area
    
    poly.area = numeric(),                                                        # Polygon size for z.yarison
    
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n){
    
    a <- i # select the vertices to cut open the polygon on
    b <- ifelse(i<n, (i + 1), 1)
    
    idx <- a
    k <- a
    
    while (k != b) {
      k <- ifelse(k == 1, n, k - 1)
      idx <- c(idx, k)
    }
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("z.x", "z.y")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("z.x","z.y")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x","z.y")]) # Points to test (calculating distance from E)
    
    AB <- B-A 
    num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
    dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
    
    k <- which.max(dist) # Identify antipode (max distance vertex)
    
    E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
    
    # Split the polygon into left and right sides
    left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
    right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
    
    # Fit data frame of these as edge sequences
    L.edges <- data.frame(
      x1 = left.verts$z.x[-nrow(left.verts)],
      y1 = left.verts$z.y[-nrow(left.verts)],
      x2 = left.verts$z.x[-1],
      y2 = left.verts$z.y[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$z.x[-nrow(right.verts)],
      y1 = right.verts$z.y[-nrow(right.verts)],
      x2 = right.verts$z.x[-1],
      y2 = right.verts$z.y[-1]
    )
    
    t.hat <- AB/ sqrt(sum(AB^2))
    n.hat <- c(-t.hat[2], t.hat[1])
    
    cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
    
    best.area <- Inf
    best <- NULL
    
    for (l in 1:nrow(L.edges)) {
      k <- L.edges[l,]
      
      p <- as.numeric(c(k[1], k[2]))
      q <- as.numeric(c(k[3], k[4]))
      d.L <- q - p
      d.L <- d.L / sqrt(sum(d.L^2))
      
      nL <- c(-d.L[2], d.L[1])
      vals <- df.hull.null$z.x * nL[1] + df.hull.null$z.y * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$z.x * nR[1] + df.hull.null$z.y * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
        
        M <- cbind(d.L, -d.R) # calculate the apex
        if (abs(det(M)) < 1e-12) next  # skip near-parallel
        
        sol <- solve(M, x0.R - x0.L)
        apex <- x0.L + sol[1] * d.L
        
        apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
        if (apex.h <= cE) next
        
        M1 <- cbind(AB, -d.L)
        sol1 <- solve(M1, x0.L - A)
        base.L <- A + sol1[1] * AB
        
        M2 <- cbind(AB, -d.R)
        sol2 <- solve(M2, x0.R - A)
        base.R <- A + sol2[1] * AB
        
        v1 <- base.L - apex
        v2 <- base.R - apex
        area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
        
        if (area < best.area) {
          best.area <- area
          best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
        }
      }
    }
    
    if (is.null(best) || !is.finite(best.area)) next
    
    tri.null <- rbind(tri.null,data.frame(
      anchor.edge = i,                                                            # Which segment is E (leading vertex)
      
      apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
      baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
      baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
      
      area = best.area,                                                           # Triangle area
      poly.area = abs(polyarea(df.hull.full.null$z.x, df.hull.full.null$z.y))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri$poly.area[1]) # 0.848

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.974

###### Pareto fronts ######

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
fit <- lm(z.y~z.x, data= par.res.1)

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

LP.scam <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/I*)", 
       title = "C — Light ~ Phosphorous") +  # labels
  
  ylim(-0.10, 3.25) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LP.scam  # Display the plot

p.lp3 <- p.lp  + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) # Adding scam PF fits
p.lp3

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

mean(null.df$a.emp.n >= a.emp) # p-value 0.941

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # -0.00008, p 0.96822
summary(q75, se = "boot", R = 1000) # -0.00121, p 0.80472
summary(q90, se = "boot", R = 1000) # 0.00409, p 0.38041

LP.qr <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/I*)", 
       title = "C — Light ~ Phosphorous") +  # labels
  
  ylim(-0.10, 3.25) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LP.qr # Display the plot

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%                           # smallest → largest
  slice((floor(0.6667 * n()) + 1):n()) %>%          # keep *second* half
  select(-distance)

LP.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/I*)", 
       title = "C — Light ~ Phosphorous") +  # labels
  
  ylim(-0.10, 3.25) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LP.scam2  # Display the plot

null.df3 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df3 <- rbind(null.df3, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.812

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -0.00503, p 0.16908
summary(q75, se = "boot", R = 1000) # -0.00400, p 0.41506
summary(q90, se = "boot", R = 1000) # 0.00301, p 0.56896

LP.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Competitive ability (1/P*)",    
       y = "Competitive ability (1/I*)", 
       title = "C — Light ~ Phosphorous") +  # labels
  
  ylim(-0.10, 3.25) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LP.qr2 # Display the plot

# Light v Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.L,
    z.x = T.br
  ) %>% # Specify the x and y variables
  filter(!is.na(z.y)) %>% 
  filter(!is.na(z.x)) # Only 46 observations

plot(z.y ~ z.x, data= df.filt) # look at the data 
# Another high point in the upper right that will drive relationships here. 

###### Shoval et al triangular analysis ######

hull.id <- chull(df.filt$z.x, df.filt$z.y)

df.hull.full <- df.filt[hull.id,] %>% 
  bind_rows(df.filt[hull.id[1],]) # Close the polygon

ggplot(df.filt, aes(z.x, z.y)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(z.x, z.y, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Find the smallest triangle ######

df.hull <- df.filt[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri <- data.frame(
  anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
  
  apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
  baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
  baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
  
  area = numeric(),                                                             # Triangle area
  
  poly.area = numeric(),                                                        # Polygon size for comparison
  
  stringsAsFactors = FALSE
)


for (i in 1:n){
  
  a <- i # select the vertices to cut open the polygon on
  b <- ifelse(i<n, (i + 1), 1)
  
  idx <- a
  k <- a
  
  while (k != b) {
    k <- ifelse(k == 1, n, k - 1)
    idx <- c(idx, k)
  }
  
  df.cut <- df.hull[idx, ] # OK our loop has succesfully cut open the polygon along the edge between 1 and 2!
  
  # Step 2 - find the antipode E' of E
  
  A <- as.numeric(df.cut[1,c("z.x", "z.y")])                      # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("z.x", "z.y")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x", "z.y")]) # Points to test (calculating distance from E)
  
  AB <- B-A 
  num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
  dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
  
  k <- which.max(dist) # Identify antipode (max distance vertex)
  
  E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
  
  # Split the polygon into left and right sides
  left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
  right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
  
  # Fit data frame of these as edge sequences
  L.edges <- data.frame(
    x1 = left.verts$z.x[-nrow(left.verts)],
    y1 = left.verts$z.y[-nrow(left.verts)],
    x2 = left.verts$z.x[-1],
    y2 = left.verts$z.y[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$z.x[-nrow(right.verts)],
    y1 = right.verts$z.y[-nrow(right.verts)],
    x2 = right.verts$z.x[-1],
    y2 = right.verts$z.y[-1]
  )
  
  t.hat <- AB/ sqrt(sum(AB^2))
  n.hat <- c(-t.hat[2], t.hat[1])
  
  cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
  
  best.area <- Inf
  best <- NULL
  
  for (l in 1:nrow(L.edges)) {
    k <- L.edges[l,]
    
    p <- as.numeric(c(k[1], k[2]))
    q <- as.numeric(c(k[3], k[4]))
    d.L <- q - p
    d.L <- d.L / sqrt(sum(d.L^2))
    
    nL <- c(-d.L[2], d.L[1])
    vals <- df.hull$z.x * nL[1] + df.hull$z.y * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$z.x[idx], df.hull$z.y[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$z.x * nR[1] + df.hull$z.y * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$z.x[idx], df.hull$z.y[idx])
      
      M <- cbind(d.L, -d.R) # calculate the apex
      if (abs(det(M)) < 1e-12) next  # skip near-parallel
      
      sol <- solve(M, x0.R - x0.L)
      apex <- x0.L + sol[1] * d.L
      
      apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
      if (apex.h <= cE) next
      
      M1 <- cbind(AB, -d.L)
      sol1 <- solve(M1, x0.L - A)
      base.L <- A + sol1[1] * AB
      
      M2 <- cbind(AB, -d.R)
      sol2 <- solve(M2, x0.R - A)
      base.R <- A + sol2[1] * AB
      
      v1 <- base.L - apex
      v2 <- base.R - apex
      area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
      
      if (area < best.area) {
        best.area <- area
        best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
      }
    }
  }
  
  tri <- rbind(tri,data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$z.x, df.hull.full$z.y))               # Polygon area
  ))
  
}

best.tri <- tri %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly <- data.frame(
  x = c(
    best.tri$apex.x,
    best.tri$baseL.x,
    best.tri$baseR.x,
    best.tri$apex.x   # close polygon
  ),
  y = c(
    best.tri$apex.y,
    best.tri$baseL.y,
    best.tri$baseR.y,
    best.tri$apex.y
  )
)

p.lt <- ggplot(df.filt, aes(z.x, z.y)) +
  geom_point(size = 2) +
  
  labs(x = "Thermal breadth (°C)",    
       y = "Competitive ability (1/I*)", 
       title = "D — Light ~ Temperature") +  # labels
  
  geom_polygon(
    data = tri.poly,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  ylim(-0.1, 3) +
  # xlim(0, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.lt

tri

best.tri$area/best.tri$poly.area # 1.271504

###### Randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.filt %>%          # scramble the data
    
    mutate(
      z.x = sample(df.filt$z.x, replace = FALSE),             # Randomly assign rmax
      z.y = sample(df.filt$z.y, replace = FALSE)              # Separately reassign z.y
    )
  
  hull.id.null <- chull(df.null$z.x, df.null$z.y)    # Get the polygon
  
  df.hull.full.null <- df.null[hull.id.null,] %>% 
    bind_rows(df.null[hull.id.null[1],]) # Close the polygon
  
  df.hull.null <- df.null[hull.id.null,]  # Don't close the polygon
  
  n <- nrow(df.hull.null) # need the number of vertices in the convex hull dataset
  
  df.hull.null <- df.hull.null[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 
  
  tri.null <- data.frame(
    anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
    
    apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
    baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
    baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
    
    area = numeric(),                                                             # Triangle area
    
    poly.area = numeric(),                                                        # Polygon size for z.yarison
    
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n){
    
    a <- i # select the vertices to cut open the polygon on
    b <- ifelse(i<n, (i + 1), 1)
    
    idx <- a
    k <- a
    
    while (k != b) {
      k <- ifelse(k == 1, n, k - 1)
      idx <- c(idx, k)
    }
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("z.x", "z.y")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("z.x","z.y")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x","z.y")]) # Points to test (calculating distance from E)
    
    AB <- B-A 
    num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
    dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
    
    k <- which.max(dist) # Identify antipode (max distance vertex)
    
    E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
    
    # Split the polygon into left and right sides
    left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
    right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
    
    # Fit data frame of these as edge sequences
    L.edges <- data.frame(
      x1 = left.verts$z.x[-nrow(left.verts)],
      y1 = left.verts$z.y[-nrow(left.verts)],
      x2 = left.verts$z.x[-1],
      y2 = left.verts$z.y[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$z.x[-nrow(right.verts)],
      y1 = right.verts$z.y[-nrow(right.verts)],
      x2 = right.verts$z.x[-1],
      y2 = right.verts$z.y[-1]
    )
    
    t.hat <- AB/ sqrt(sum(AB^2))
    n.hat <- c(-t.hat[2], t.hat[1])
    
    cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
    
    best.area <- Inf
    best <- NULL
    
    for (l in 1:nrow(L.edges)) {
      k <- L.edges[l,]
      
      p <- as.numeric(c(k[1], k[2]))
      q <- as.numeric(c(k[3], k[4]))
      d.L <- q - p
      d.L <- d.L / sqrt(sum(d.L^2))
      
      nL <- c(-d.L[2], d.L[1])
      vals <- df.hull.null$z.x * nL[1] + df.hull.null$z.y * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$z.x * nR[1] + df.hull.null$z.y * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
        
        M <- cbind(d.L, -d.R) # calculate the apex
        if (abs(det(M)) < 1e-12) next  # skip near-parallel
        
        sol <- solve(M, x0.R - x0.L)
        apex <- x0.L + sol[1] * d.L
        
        apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
        if (apex.h <= cE) next
        
        M1 <- cbind(AB, -d.L)
        sol1 <- solve(M1, x0.L - A)
        base.L <- A + sol1[1] * AB
        
        M2 <- cbind(AB, -d.R)
        sol2 <- solve(M2, x0.R - A)
        base.R <- A + sol2[1] * AB
        
        v1 <- base.L - apex
        v2 <- base.R - apex
        area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
        
        if (area < best.area) {
          best.area <- area
          best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
        }
      }
    }
    
    if (is.null(best) || !is.finite(best.area)) next
    
    tri.null <- rbind(tri.null,data.frame(
      anchor.edge = i,                                                            # Which segment is E (leading vertex)
      
      apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
      baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
      baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
      
      area = best.area,                                                           # Triangle area
      poly.area = abs(polyarea(df.hull.full.null$z.x, df.hull.full.null$z.y))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri$poly.area[1]) # 0.924

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.752

###### Pareto fronts ######

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
fit <- lm(z.y~z.x, data= par.res.1)

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

LT.scam <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Thermal breadth (°C)",    
       y = "Competitive ability (1/I*)", 
       title = "D — Light ~ Temperature") +  # labels
  
   ylim(-0.10, 2.25) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LT.scam  # Display the plot

p.lt3 <- p.lt  + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) # Adding scam PF fits
p.lt3

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

mean(null.df$a.emp.n >= a.emp) # p-value 0.836

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 0.01047, p 0.33973
summary(q75, se = "boot", R = 1000) # 0.01314, p 0.67341
summary(q90, se = "boot", R = 1000) # 0.00230, p 0.96289

LT.qr <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Thermal breadth (°C)",    
       y = "Competitive ability (1/I*)", 
       title = "D — Light ~ Temperature") +  # labels
  
  ylim(-0.10, 2.25) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LT.qr # Display the plot

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%                           # smallest → largest
  slice((floor(0.6667 * n()) + 1):n()) %>%          # keep *second* half
  select(-distance)

LT.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Thermal breadth (°C)",    
       y = "Competitive ability (1/I*)", 
       title = "D — Light ~ Temperature") +  # labels
  
  ylim(-0.10, 2.25) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LT.scam2  # Display the plot

null.df3 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df3 <- rbind(null.df3, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.492

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -0.04617, p 0.02528
summary(q75, se = "boot", R = 1000) # -0.04198, p 0.43375
summary(q90, se = "boot", R = 1000) # 0.00787, p 0.90230

LT.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Thermal breadth (°C)",    
       y = "Competitive ability (1/I*)", 
       title = "D — Light ~ Temperature") +  # labels
  
  ylim(-0.10, 2.25) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

LT.qr2 # Display the plot

# Nitrogen v Nitrogen -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.N,
    z.x = r.max.N
  ) %>% # Specify the x and y variables
  filter(!is.na(z.y)) %>% 
  filter(!is.na(z.x)) # Only 31 observations

plot(z.y ~ z.x, data= df.filt) # look at the data 

###### Shoval et al triangular analysis ######

hull.id <- chull(df.filt$z.x, df.filt$z.y)

df.hull.full <- df.filt[hull.id,] %>% 
  bind_rows(df.filt[hull.id[1],]) # Close the polygon

ggplot(df.filt, aes(z.x, z.y)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(z.x, z.y, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Find the smallest triangle ######

df.hull <- df.filt[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri <- data.frame(
  anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
  
  apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
  baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
  baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
  
  area = numeric(),                                                             # Triangle area
  
  poly.area = numeric(),                                                        # Polygon size for comparison
  
  stringsAsFactors = FALSE
)


for (i in 1:n){
  
  a <- i # select the vertices to cut open the polygon on
  b <- ifelse(i<n, (i + 1), 1)
  
  idx <- a
  k <- a
  
  while (k != b) {
    k <- ifelse(k == 1, n, k - 1)
    idx <- c(idx, k)
  }
  
  df.cut <- df.hull[idx, ] # OK our loop has succesfully cut open the polygon along the edge between 1 and 2!
  
  # Step 2 - find the antipode E' of E
  
  A <- as.numeric(df.cut[1,c("z.x", "z.y")])                      # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("z.x", "z.y")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x", "z.y")]) # Points to test (calculating distance from E)
  
  AB <- B-A 
  num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
  dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
  
  k <- which.max(dist) # Identify antipode (max distance vertex)
  
  E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
  
  # Split the polygon into left and right sides
  left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
  right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
  
  # Fit data frame of these as edge sequences
  L.edges <- data.frame(
    x1 = left.verts$z.x[-nrow(left.verts)],
    y1 = left.verts$z.y[-nrow(left.verts)],
    x2 = left.verts$z.x[-1],
    y2 = left.verts$z.y[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$z.x[-nrow(right.verts)],
    y1 = right.verts$z.y[-nrow(right.verts)],
    x2 = right.verts$z.x[-1],
    y2 = right.verts$z.y[-1]
  )
  
  t.hat <- AB/ sqrt(sum(AB^2))
  n.hat <- c(-t.hat[2], t.hat[1])
  
  cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
  
  best.area <- Inf
  best <- NULL
  
  for (l in 1:nrow(L.edges)) {
    k <- L.edges[l,]
    
    p <- as.numeric(c(k[1], k[2]))
    q <- as.numeric(c(k[3], k[4]))
    d.L <- q - p
    d.L <- d.L / sqrt(sum(d.L^2))
    
    nL <- c(-d.L[2], d.L[1])
    vals <- df.hull$z.x * nL[1] + df.hull$z.y * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$z.x[idx], df.hull$z.y[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$z.x * nR[1] + df.hull$z.y * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$z.x[idx], df.hull$z.y[idx])
      
      M <- cbind(d.L, -d.R) # calculate the apex
      if (abs(det(M)) < 1e-12) next  # skip near-parallel
      
      sol <- solve(M, x0.R - x0.L)
      apex <- x0.L + sol[1] * d.L
      
      apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
      if (apex.h <= cE) next
      
      M1 <- cbind(AB, -d.L)
      sol1 <- solve(M1, x0.L - A)
      base.L <- A + sol1[1] * AB
      
      M2 <- cbind(AB, -d.R)
      sol2 <- solve(M2, x0.R - A)
      base.R <- A + sol2[1] * AB
      
      v1 <- base.L - apex
      v2 <- base.R - apex
      area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
      
      if (area < best.area) {
        best.area <- area
        best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
      }
    }
  }
  
  tri <- rbind(tri,data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$z.x, df.hull.full$z.y))               # Polygon area
  ))
  
}

best.tri <- tri %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly <- data.frame(
  x = c(
    best.tri$apex.x,
    best.tri$baseL.x,
    best.tri$baseR.x,
    best.tri$apex.x   # close polygon
  ),
  y = c(
    best.tri$apex.y,
    best.tri$baseL.y,
    best.tri$baseR.y,
    best.tri$apex.y
  )
)

p.n <- ggplot(df.filt, aes(z.x, z.y)) +
  geom_point(size = 2) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/N*)", 
       title = "E — Nitrogen") +  # labels
  
  geom_polygon(
    data = tri.poly,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  ylim(-0.1, 55) +
  # xlim(0, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.n

tri

best.tri$area/best.tri$poly.area # 1.051291

###### Randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.filt %>%          # scramble the data
    
    mutate(
      z.x = sample(df.filt$z.x, replace = FALSE),             # Randomly assign rmax
      z.y = sample(df.filt$z.y, replace = FALSE)              # Separately reassign z.y
    )
  
  hull.id.null <- chull(df.null$z.x, df.null$z.y)    # Get the polygon
  
  df.hull.full.null <- df.null[hull.id.null,] %>% 
    bind_rows(df.null[hull.id.null[1],]) # Close the polygon
  
  df.hull.null <- df.null[hull.id.null,]  # Don't close the polygon
  
  n <- nrow(df.hull.null) # need the number of vertices in the convex hull dataset
  
  df.hull.null <- df.hull.null[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 
  
  tri.null <- data.frame(
    anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
    
    apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
    baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
    baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
    
    area = numeric(),                                                             # Triangle area
    
    poly.area = numeric(),                                                        # Polygon size for z.yarison
    
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n){
    
    a <- i # select the vertices to cut open the polygon on
    b <- ifelse(i<n, (i + 1), 1)
    
    idx <- a
    k <- a
    
    while (k != b) {
      k <- ifelse(k == 1, n, k - 1)
      idx <- c(idx, k)
    }
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("z.x", "z.y")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("z.x","z.y")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x","z.y")]) # Points to test (calculating distance from E)
    
    AB <- B-A 
    num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
    dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
    
    k <- which.max(dist) # Identify antipode (max distance vertex)
    
    E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
    
    # Split the polygon into left and right sides
    left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
    right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
    
    # Fit data frame of these as edge sequences
    L.edges <- data.frame(
      x1 = left.verts$z.x[-nrow(left.verts)],
      y1 = left.verts$z.y[-nrow(left.verts)],
      x2 = left.verts$z.x[-1],
      y2 = left.verts$z.y[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$z.x[-nrow(right.verts)],
      y1 = right.verts$z.y[-nrow(right.verts)],
      x2 = right.verts$z.x[-1],
      y2 = right.verts$z.y[-1]
    )
    
    t.hat <- AB/ sqrt(sum(AB^2))
    n.hat <- c(-t.hat[2], t.hat[1])
    
    cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
    
    best.area <- Inf
    best <- NULL
    
    for (l in 1:nrow(L.edges)) {
      k <- L.edges[l,]
      
      p <- as.numeric(c(k[1], k[2]))
      q <- as.numeric(c(k[3], k[4]))
      d.L <- q - p
      d.L <- d.L / sqrt(sum(d.L^2))
      
      nL <- c(-d.L[2], d.L[1])
      vals <- df.hull.null$z.x * nL[1] + df.hull.null$z.y * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$z.x * nR[1] + df.hull.null$z.y * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
        
        M <- cbind(d.L, -d.R) # calculate the apex
        if (abs(det(M)) < 1e-12) next  # skip near-parallel
        
        sol <- solve(M, x0.R - x0.L)
        apex <- x0.L + sol[1] * d.L
        
        apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
        if (apex.h <= cE) next
        
        M1 <- cbind(AB, -d.L)
        sol1 <- solve(M1, x0.L - A)
        base.L <- A + sol1[1] * AB
        
        M2 <- cbind(AB, -d.R)
        sol2 <- solve(M2, x0.R - A)
        base.R <- A + sol2[1] * AB
        
        v1 <- base.L - apex
        v2 <- base.R - apex
        area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
        
        if (area < best.area) {
          best.area <- area
          best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
        }
      }
    }
    
    if (is.null(best) || !is.finite(best.area)) next
    
    tri.null <- rbind(tri.null,data.frame(
      anchor.edge = i,                                                            # Which segment is E (leading vertex)
      
      apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
      baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
      baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
      
      area = best.area,                                                           # Triangle area
      poly.area = abs(polyarea(df.hull.full.null$z.x, df.hull.full.null$z.y))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri$poly.area[1]) # 0.138

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.099

###### Pareto fronts ######

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
fit2 <- lm(z.y~z.x, data = par.res.2)

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

N.scam <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/N*)", 
       title = "E — Nitrogen") +  # labels
  
  ylim(-0.10, 55) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

N.scam  # Display the plot

p.n3 <- p.n  + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) # Adding scam PF fits
p.n3

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

mean(null.df$a.emp.n >= a.emp) # p-value 0.635

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 9.76601, p 0.041863
summary(q75, se = "boot", R = 1000) # 13.90848, p 0.15975
summary(q90, se = "boot", R = 1000) # 34.01186, p 0.00298

N.qr <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/N*)", 
       title = "E — Nitrogen") +  # labels
  
  ylim(-0.10, 55) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

N.qr # Display the plot

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%                           # smallest → largest
  slice((floor(0.6667 * n()) + 1):n()) %>%          # keep *second* half
  select(-distance)

N.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/N*)", 
       title = "E — Nitrogen") +  # labels
  
  ylim(-0.10, 55) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

N.scam2  # Display the plot

null.df3 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df3 <- rbind(null.df3, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.101

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -13.41769, p 0.53543
summary(q75, se = "boot", R = 1000) # -27.60364, p 0.33548
summary(q90, se = "boot", R = 1000) # -26.32663, p 0.45281

N.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/N*)", 
       title = "E — Nitrogen") +  # labels
  
  ylim(-0.10, 55) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

N.qr2 # Display the plot

# Nitrogen v Phosphorous -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.N,
    z.x = comp.P
  ) %>% # Specify the x and y variables
  filter(!is.na(z.y)) %>% 
  filter(!is.na(z.x)) # Only 23 observations

plot(z.y ~ z.x, data= df.filt) # look at the data 

###### Shoval et al triangular analysis ######

hull.id <- chull(df.filt$z.x, df.filt$z.y)

df.hull.full <- df.filt[hull.id,] %>% 
  bind_rows(df.filt[hull.id[1],]) # Close the polygon

ggplot(df.filt, aes(z.x, z.y)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(z.x, z.y, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Find the smallest triangle ######

df.hull <- df.filt[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri <- data.frame(
  anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
  
  apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
  baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
  baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
  
  area = numeric(),                                                             # Triangle area
  
  poly.area = numeric(),                                                        # Polygon size for comparison
  
  stringsAsFactors = FALSE
)


for (i in 1:n){
  
  a <- i # select the vertices to cut open the polygon on
  b <- ifelse(i<n, (i + 1), 1)
  
  idx <- a
  k <- a
  
  while (k != b) {
    k <- ifelse(k == 1, n, k - 1)
    idx <- c(idx, k)
  }
  
  df.cut <- df.hull[idx, ] # OK our loop has succesfully cut open the polygon along the edge between 1 and 2!
  
  # Step 2 - find the antipode E' of E
  
  A <- as.numeric(df.cut[1,c("z.x", "z.y")])                      # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("z.x", "z.y")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x", "z.y")]) # Points to test (calculating distance from E)
  
  AB <- B-A 
  num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
  dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
  
  k <- which.max(dist) # Identify antipode (max distance vertex)
  
  E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
  
  # Split the polygon into left and right sides
  left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
  right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
  
  # Fit data frame of these as edge sequences
  L.edges <- data.frame(
    x1 = left.verts$z.x[-nrow(left.verts)],
    y1 = left.verts$z.y[-nrow(left.verts)],
    x2 = left.verts$z.x[-1],
    y2 = left.verts$z.y[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$z.x[-nrow(right.verts)],
    y1 = right.verts$z.y[-nrow(right.verts)],
    x2 = right.verts$z.x[-1],
    y2 = right.verts$z.y[-1]
  )
  
  t.hat <- AB/ sqrt(sum(AB^2))
  n.hat <- c(-t.hat[2], t.hat[1])
  
  cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
  
  best.area <- Inf
  best <- NULL
  
  for (l in 1:nrow(L.edges)) {
    k <- L.edges[l,]
    
    p <- as.numeric(c(k[1], k[2]))
    q <- as.numeric(c(k[3], k[4]))
    d.L <- q - p
    d.L <- d.L / sqrt(sum(d.L^2))
    
    nL <- c(-d.L[2], d.L[1])
    vals <- df.hull$z.x * nL[1] + df.hull$z.y * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$z.x[idx], df.hull$z.y[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$z.x * nR[1] + df.hull$z.y * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$z.x[idx], df.hull$z.y[idx])
      
      M <- cbind(d.L, -d.R) # calculate the apex
      if (abs(det(M)) < 1e-12) next  # skip near-parallel
      
      sol <- solve(M, x0.R - x0.L)
      apex <- x0.L + sol[1] * d.L
      
      apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
      if (apex.h <= cE) next
      
      M1 <- cbind(AB, -d.L)
      sol1 <- solve(M1, x0.L - A)
      base.L <- A + sol1[1] * AB
      
      M2 <- cbind(AB, -d.R)
      sol2 <- solve(M2, x0.R - A)
      base.R <- A + sol2[1] * AB
      
      v1 <- base.L - apex
      v2 <- base.R - apex
      area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
      
      if (area < best.area) {
        best.area <- area
        best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
      }
    }
  }
  
  tri <- rbind(tri,data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$z.x, df.hull.full$z.y))               # Polygon area
  ))
  
}

best.tri <- tri %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly <- data.frame(
  x = c(
    best.tri$apex.x,
    best.tri$baseL.x,
    best.tri$baseR.x,
    best.tri$apex.x   # close polygon
  ),
  y = c(
    best.tri$apex.y,
    best.tri$baseL.y,
    best.tri$baseR.y,
    best.tri$apex.y
  )
)

p.np <- ggplot(df.filt, aes(z.x, z.y)) +
  geom_point(size = 2) +
  
  labs(x = "Competitive ability (1/P*)", 
       y = "Competitive ability (1/N*)", 
       title = "F — Nitrogen ~ Phosphorous") +  # labels
  
  geom_polygon(
    data = tri.poly,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  ylim(-0.1, 85) +
  # xlim(0, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.np

tri

best.tri$area/best.tri$poly.area # 1.268966

###### Randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.filt %>%          # scramble the data
    
    mutate(
      z.x = sample(df.filt$z.x, replace = FALSE),             # Randomly assign rmax
      z.y = sample(df.filt$z.y, replace = FALSE)              # Separately reassign z.y
    )
  
  hull.id.null <- chull(df.null$z.x, df.null$z.y)    # Get the polygon
  
  df.hull.full.null <- df.null[hull.id.null,] %>% 
    bind_rows(df.null[hull.id.null[1],]) # Close the polygon
  
  df.hull.null <- df.null[hull.id.null,]  # Don't close the polygon
  
  n <- nrow(df.hull.null) # need the number of vertices in the convex hull dataset
  
  df.hull.null <- df.hull.null[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 
  
  tri.null <- data.frame(
    anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
    
    apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
    baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
    baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
    
    area = numeric(),                                                             # Triangle area
    
    poly.area = numeric(),                                                        # Polygon size for z.yarison
    
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n){
    
    a <- i # select the vertices to cut open the polygon on
    b <- ifelse(i<n, (i + 1), 1)
    
    idx <- a
    k <- a
    
    while (k != b) {
      k <- ifelse(k == 1, n, k - 1)
      idx <- c(idx, k)
    }
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("z.x", "z.y")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("z.x","z.y")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x","z.y")]) # Points to test (calculating distance from E)
    
    AB <- B-A 
    num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
    dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
    
    k <- which.max(dist) # Identify antipode (max distance vertex)
    
    E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
    
    # Split the polygon into left and right sides
    left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
    right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
    
    # Fit data frame of these as edge sequences
    L.edges <- data.frame(
      x1 = left.verts$z.x[-nrow(left.verts)],
      y1 = left.verts$z.y[-nrow(left.verts)],
      x2 = left.verts$z.x[-1],
      y2 = left.verts$z.y[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$z.x[-nrow(right.verts)],
      y1 = right.verts$z.y[-nrow(right.verts)],
      x2 = right.verts$z.x[-1],
      y2 = right.verts$z.y[-1]
    )
    
    t.hat <- AB/ sqrt(sum(AB^2))
    n.hat <- c(-t.hat[2], t.hat[1])
    
    cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
    
    best.area <- Inf
    best <- NULL
    
    for (l in 1:nrow(L.edges)) {
      k <- L.edges[l,]
      
      p <- as.numeric(c(k[1], k[2]))
      q <- as.numeric(c(k[3], k[4]))
      d.L <- q - p
      d.L <- d.L / sqrt(sum(d.L^2))
      
      nL <- c(-d.L[2], d.L[1])
      vals <- df.hull.null$z.x * nL[1] + df.hull.null$z.y * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$z.x * nR[1] + df.hull.null$z.y * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
        
        M <- cbind(d.L, -d.R) # calculate the apex
        if (abs(det(M)) < 1e-12) next  # skip near-parallel
        
        sol <- solve(M, x0.R - x0.L)
        apex <- x0.L + sol[1] * d.L
        
        apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
        if (apex.h <= cE) next
        
        M1 <- cbind(AB, -d.L)
        sol1 <- solve(M1, x0.L - A)
        base.L <- A + sol1[1] * AB
        
        M2 <- cbind(AB, -d.R)
        sol2 <- solve(M2, x0.R - A)
        base.R <- A + sol2[1] * AB
        
        v1 <- base.L - apex
        v2 <- base.R - apex
        area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
        
        if (area < best.area) {
          best.area <- area
          best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
        }
      }
    }
    
    if (is.null(best) || !is.finite(best.area)) next
    
    tri.null <- rbind(tri.null,data.frame(
      anchor.edge = i,                                                            # Which segment is E (leading vertex)
      
      apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
      baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
      baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
      
      area = best.area,                                                           # Triangle area
      poly.area = abs(polyarea(df.hull.full.null$z.x, df.hull.full.null$z.y))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri$poly.area[1]) # 0.642

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.690

###### Pareto fronts ######

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
fit <- lm(z.y~z.x, data = par.res.1)

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

NP.scam <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Competitive ability (1/P*)", 
       y = "Competitive ability (1/N*)", 
       title = "F — Nitrogen ~ Phosphorous") +  # labels
  
  ylim(-0.10, 55) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NP.scam  # Display the plot

p.np3 <- p.np  + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) # Adding scam PF fits
p.np3

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

mean(null.df$a.emp.n >= a.emp) # p-value 0.414

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # -0.01522, p 0.88197
summary(q75, se = "boot", R = 1000) # -0.03540, p 0.87148
summary(q90, se = "boot", R = 1000) # 0.23716, p 0.35713

NP.qr <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Competitive ability (1/P*)", 
       y = "Competitive ability (1/N*)", 
       title = "F — Nitrogen ~ Phosphorous") +  # labels
  
  ylim(-0.10, 55) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NP.qr # Display the plot

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%                           # smallest → largest
  slice((floor(0.6667 * n()) + 1):n()) %>%          # keep *second* half
  select(-distance)

NP.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Competitive ability (1/P*)", 
       y = "Competitive ability (1/N*)", 
       title = "F — Nitrogen ~ Phosphorous") +  # labels
  
  ylim(-0.10, 55) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NP.scam2  # Display the plot

null.df3 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df3 <- rbind(null.df3, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.044

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -0.18707, p 0.14409
summary(q75, se = "boot", R = 1000) # -0.22174, p 0.30668
summary(q90, se = "boot", R = 1000) # 0.10430, p 0.71859

NP.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Competitive ability (1/P*)", 
       y = "Competitive ability (1/N*)", 
       title = "F — Nitrogen ~ Phosphorous") +  # labels
  
  ylim(-0.10, 55) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NP.qr2 # Display the plot

# Nitrogen v Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.N,
    z.x = T.br
  ) %>% # Specify the x and y variables
  filter(!is.na(z.y)) %>% 
  filter(!is.na(z.x)) # Only 23 observations

plot(z.y ~ z.x, data= df.filt) # look at the data 

###### Shoval et al triangular analysis ######

hull.id <- chull(df.filt$z.x, df.filt$z.y)

df.hull.full <- df.filt[hull.id,] %>% 
  bind_rows(df.filt[hull.id[1],]) # Close the polygon

ggplot(df.filt, aes(z.x, z.y)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(z.x, z.y, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Find the smallest triangle ######

df.hull <- df.filt[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri <- data.frame(
  anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
  
  apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
  baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
  baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
  
  area = numeric(),                                                             # Triangle area
  
  poly.area = numeric(),                                                        # Polygon size for comparison
  
  stringsAsFactors = FALSE
)


for (i in 1:n){
  
  a <- i # select the vertices to cut open the polygon on
  b <- ifelse(i<n, (i + 1), 1)
  
  idx <- a
  k <- a
  
  while (k != b) {
    k <- ifelse(k == 1, n, k - 1)
    idx <- c(idx, k)
  }
  
  df.cut <- df.hull[idx, ] # OK our loop has succesfully cut open the polygon along the edge between 1 and 2!
  
  # Step 2 - find the antipode E' of E
  
  A <- as.numeric(df.cut[1,c("z.x", "z.y")])                      # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("z.x", "z.y")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x", "z.y")]) # Points to test (calculating distance from E)
  
  AB <- B-A 
  num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
  dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
  
  k <- which.max(dist) # Identify antipode (max distance vertex)
  
  E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
  
  # Split the polygon into left and right sides
  left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
  right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
  
  # Fit data frame of these as edge sequences
  L.edges <- data.frame(
    x1 = left.verts$z.x[-nrow(left.verts)],
    y1 = left.verts$z.y[-nrow(left.verts)],
    x2 = left.verts$z.x[-1],
    y2 = left.verts$z.y[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$z.x[-nrow(right.verts)],
    y1 = right.verts$z.y[-nrow(right.verts)],
    x2 = right.verts$z.x[-1],
    y2 = right.verts$z.y[-1]
  )
  
  t.hat <- AB/ sqrt(sum(AB^2))
  n.hat <- c(-t.hat[2], t.hat[1])
  
  cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
  
  best.area <- Inf
  best <- NULL
  
  for (l in 1:nrow(L.edges)) {
    k <- L.edges[l,]
    
    p <- as.numeric(c(k[1], k[2]))
    q <- as.numeric(c(k[3], k[4]))
    d.L <- q - p
    d.L <- d.L / sqrt(sum(d.L^2))
    
    nL <- c(-d.L[2], d.L[1])
    vals <- df.hull$z.x * nL[1] + df.hull$z.y * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$z.x[idx], df.hull$z.y[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$z.x * nR[1] + df.hull$z.y * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$z.x[idx], df.hull$z.y[idx])
      
      M <- cbind(d.L, -d.R) # calculate the apex
      if (abs(det(M)) < 1e-12) next  # skip near-parallel
      
      sol <- solve(M, x0.R - x0.L)
      apex <- x0.L + sol[1] * d.L
      
      apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
      if (apex.h <= cE) next
      
      M1 <- cbind(AB, -d.L)
      sol1 <- solve(M1, x0.L - A)
      base.L <- A + sol1[1] * AB
      
      M2 <- cbind(AB, -d.R)
      sol2 <- solve(M2, x0.R - A)
      base.R <- A + sol2[1] * AB
      
      v1 <- base.L - apex
      v2 <- base.R - apex
      area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
      
      if (area < best.area) {
        best.area <- area
        best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
      }
    }
  }
  
  tri <- rbind(tri,data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$z.x, df.hull.full$z.y))               # Polygon area
  ))
  
}

best.tri <- tri %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly <- data.frame(
  x = c(
    best.tri$apex.x,
    best.tri$baseL.x,
    best.tri$baseR.x,
    best.tri$apex.x   # close polygon
  ),
  y = c(
    best.tri$apex.y,
    best.tri$baseL.y,
    best.tri$baseR.y,
    best.tri$apex.y
  )
)

p.nt <- ggplot(df.filt, aes(z.x, z.y)) +
  geom_point(size = 2) +
  
  labs(x = "Thermal breadth (°C)", 
       y = "Competitive ability (1/N*)", 
       title = "G — Nitrogen ~ Temperature") +  # labels
  
  geom_polygon(
    data = tri.poly,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  ylim(-1, 45) +
  # xlim(0, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.nt

tri

best.tri$area/best.tri$poly.area # 1.011667

###### Randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.filt %>%          # scramble the data
    
    mutate(
      z.x = sample(df.filt$z.x, replace = FALSE),             # Randomly assign rmax
      z.y = sample(df.filt$z.y, replace = FALSE)              # Separately reassign z.y
    )
  
  hull.id.null <- chull(df.null$z.x, df.null$z.y)    # Get the polygon
  
  df.hull.full.null <- df.null[hull.id.null,] %>% 
    bind_rows(df.null[hull.id.null[1],]) # Close the polygon
  
  df.hull.null <- df.null[hull.id.null,]  # Don't close the polygon
  
  n <- nrow(df.hull.null) # need the number of vertices in the convex hull dataset
  
  df.hull.null <- df.hull.null[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 
  
  tri.null <- data.frame(
    anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
    
    apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
    baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
    baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
    
    area = numeric(),                                                             # Triangle area
    
    poly.area = numeric(),                                                        # Polygon size for z.yarison
    
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n){
    
    a <- i # select the vertices to cut open the polygon on
    b <- ifelse(i<n, (i + 1), 1)
    
    idx <- a
    k <- a
    
    while (k != b) {
      k <- ifelse(k == 1, n, k - 1)
      idx <- c(idx, k)
    }
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("z.x", "z.y")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("z.x","z.y")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x","z.y")]) # Points to test (calculating distance from E)
    
    AB <- B-A 
    num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
    dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
    
    k <- which.max(dist) # Identify antipode (max distance vertex)
    
    E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
    
    # Split the polygon into left and right sides
    left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
    right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
    
    # Fit data frame of these as edge sequences
    L.edges <- data.frame(
      x1 = left.verts$z.x[-nrow(left.verts)],
      y1 = left.verts$z.y[-nrow(left.verts)],
      x2 = left.verts$z.x[-1],
      y2 = left.verts$z.y[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$z.x[-nrow(right.verts)],
      y1 = right.verts$z.y[-nrow(right.verts)],
      x2 = right.verts$z.x[-1],
      y2 = right.verts$z.y[-1]
    )
    
    t.hat <- AB/ sqrt(sum(AB^2))
    n.hat <- c(-t.hat[2], t.hat[1])
    
    cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
    
    best.area <- Inf
    best <- NULL
    
    for (l in 1:nrow(L.edges)) {
      k <- L.edges[l,]
      
      p <- as.numeric(c(k[1], k[2]))
      q <- as.numeric(c(k[3], k[4]))
      d.L <- q - p
      d.L <- d.L / sqrt(sum(d.L^2))
      
      nL <- c(-d.L[2], d.L[1])
      vals <- df.hull.null$z.x * nL[1] + df.hull.null$z.y * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$z.x * nR[1] + df.hull.null$z.y * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
        
        M <- cbind(d.L, -d.R) # calculate the apex
        if (abs(det(M)) < 1e-12) next  # skip near-parallel
        
        sol <- solve(M, x0.R - x0.L)
        apex <- x0.L + sol[1] * d.L
        
        apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
        if (apex.h <= cE) next
        
        M1 <- cbind(AB, -d.L)
        sol1 <- solve(M1, x0.L - A)
        base.L <- A + sol1[1] * AB
        
        M2 <- cbind(AB, -d.R)
        sol2 <- solve(M2, x0.R - A)
        base.R <- A + sol2[1] * AB
        
        v1 <- base.L - apex
        v2 <- base.R - apex
        area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
        
        if (area < best.area) {
          best.area <- area
          best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
        }
      }
    }
    
    if (is.null(best) || !is.finite(best.area)) next
    
    tri.null <- rbind(tri.null,data.frame(
      anchor.edge = i,                                                            # Which segment is E (leading vertex)
      
      apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
      baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
      baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
      
      area = best.area,                                                           # Triangle area
      poly.area = abs(polyarea(df.hull.full.null$z.x, df.hull.full.null$z.y))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri$poly.area[1]) # 0.043

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.032

###### Pareto fronts ######

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

NT.scam <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Thermal breadth (°C)", 
       y = "Competitive ability (1/N*)", 
       title = "G — Nitrogen ~ Temperature") +  # labels
  
  ylim(-0.10, 45) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NT.scam  # Display the plot

p.nt3 <- p.nt  + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) # Adding scam PF fits
p.nt3

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

mean(null.df$a.emp.n >= a.emp) # p-value 0.064

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # -0.27909, p 0.72781
summary(q75, se = "boot", R = 1000) # -1.32741, p 0.28479
summary(q90, se = "boot", R = 1000) # -2.63843, p 0.15974

NT.qr <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Thermal breadth (°C)", 
       y = "Competitive ability (1/N*)", 
       title = "G — Nitrogen ~ Temperature") +  # labels
  
  ylim(-0.10, 45) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NT.qr # Display the plot

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%                           # smallest → largest
  slice((floor(0.6667 * n()) + 1):n()) %>%          # keep *second* half
  select(-distance)

NT.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Thermal breadth (°C)", 
       y = "Competitive ability (1/N*)", 
       title = "G — Nitrogen ~ Temperature") +  # labels
  
  ylim(-0.10, 45) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NT.scam2  # Display the plot

null.df3 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df3 <- rbind(null.df3, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.001

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -2.27203, p 0.20084
summary(q75, se = "boot", R = 1000) # -2.63843, p 0.10835
summary(q90, se = "boot", R = 1000) # -3.28776, p 0.04782

NT.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Thermal breadth (°C)", 
       y = "Competitive ability (1/N*)", 
       title = "G — Nitrogen ~ Temperature") +  # labels
  
  ylim(-0.10, 45) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

NT.qr2 # Display the plot

# Phosphorous v Phosphorous -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.P,
    z.x = r.max.P
  ) %>% # Specify the x and y variables
  filter(!is.na(z.y)) %>% 
  filter(!is.na(z.x)) # Only 31 observations

plot(z.y ~ z.x, data= df.filt) # look at the data 

###### Shoval et al triangular analysis ######

hull.id <- chull(df.filt$z.x, df.filt$z.y)

df.hull.full <- df.filt[hull.id,] %>% 
  bind_rows(df.filt[hull.id[1],]) # Close the polygon

ggplot(df.filt, aes(z.x, z.y)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(z.x, z.y, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Find the smallest triangle ######

df.hull <- df.filt[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri <- data.frame(
  anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
  
  apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
  baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
  baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
  
  area = numeric(),                                                             # Triangle area
  
  poly.area = numeric(),                                                        # Polygon size for comparison
  
  stringsAsFactors = FALSE
)


for (i in 1:n){
  
  a <- i # select the vertices to cut open the polygon on
  b <- ifelse(i<n, (i + 1), 1)
  
  idx <- a
  k <- a
  
  while (k != b) {
    k <- ifelse(k == 1, n, k - 1)
    idx <- c(idx, k)
  }
  
  df.cut <- df.hull[idx, ] # OK our loop has succesfully cut open the polygon along the edge between 1 and 2!
  
  # Step 2 - find the antipode E' of E
  
  A <- as.numeric(df.cut[1,c("z.x", "z.y")])                      # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("z.x", "z.y")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x", "z.y")]) # Points to test (calculating distance from E)
  
  AB <- B-A 
  num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
  dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
  
  k <- which.max(dist) # Identify antipode (max distance vertex)
  
  E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
  
  # Split the polygon into left and right sides
  left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
  right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
  
  # Fit data frame of these as edge sequences
  L.edges <- data.frame(
    x1 = left.verts$z.x[-nrow(left.verts)],
    y1 = left.verts$z.y[-nrow(left.verts)],
    x2 = left.verts$z.x[-1],
    y2 = left.verts$z.y[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$z.x[-nrow(right.verts)],
    y1 = right.verts$z.y[-nrow(right.verts)],
    x2 = right.verts$z.x[-1],
    y2 = right.verts$z.y[-1]
  )
  
  t.hat <- AB/ sqrt(sum(AB^2))
  n.hat <- c(-t.hat[2], t.hat[1])
  
  cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
  
  best.area <- Inf
  best <- NULL
  
  for (l in 1:nrow(L.edges)) {
    k <- L.edges[l,]
    
    p <- as.numeric(c(k[1], k[2]))
    q <- as.numeric(c(k[3], k[4]))
    d.L <- q - p
    d.L <- d.L / sqrt(sum(d.L^2))
    
    nL <- c(-d.L[2], d.L[1])
    vals <- df.hull$z.x * nL[1] + df.hull$z.y * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$z.x[idx], df.hull$z.y[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$z.x * nR[1] + df.hull$z.y * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$z.x[idx], df.hull$z.y[idx])
      
      M <- cbind(d.L, -d.R) # calculate the apex
      if (abs(det(M)) < 1e-12) next  # skip near-parallel
      
      sol <- solve(M, x0.R - x0.L)
      apex <- x0.L + sol[1] * d.L
      
      apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
      if (apex.h <= cE) next
      
      M1 <- cbind(AB, -d.L)
      sol1 <- solve(M1, x0.L - A)
      base.L <- A + sol1[1] * AB
      
      M2 <- cbind(AB, -d.R)
      sol2 <- solve(M2, x0.R - A)
      base.R <- A + sol2[1] * AB
      
      v1 <- base.L - apex
      v2 <- base.R - apex
      area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
      
      if (area < best.area) {
        best.area <- area
        best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
      }
    }
  }
  
  tri <- rbind(tri,data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$z.x, df.hull.full$z.y))               # Polygon area
  ))
  
}

best.tri <- tri %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly <- data.frame(
  x = c(
    best.tri$apex.x,
    best.tri$baseL.x,
    best.tri$baseR.x,
    best.tri$apex.x   # close polygon
  ),
  y = c(
    best.tri$apex.y,
    best.tri$baseL.y,
    best.tri$baseR.y,
    best.tri$apex.y
  )
)

p.p <- ggplot(df.filt, aes(z.x, z.y)) +
  geom_point(size = 2) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/P*)", 
       title = "H — Phosphorous") +  # labels
  
  geom_polygon(
    data = tri.poly,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  ylim(-0.1, 350) +
  # xlim(0, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.p

tri

best.tri$area/best.tri$poly.area # 1.118166

###### Randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.filt %>%          # scramble the data
    
    mutate(
      z.x = sample(df.filt$z.x, replace = FALSE),             # Randomly assign rmax
      z.y = sample(df.filt$z.y, replace = FALSE)              # Separately reassign z.y
    )
  
  hull.id.null <- chull(df.null$z.x, df.null$z.y)    # Get the polygon
  
  df.hull.full.null <- df.null[hull.id.null,] %>% 
    bind_rows(df.null[hull.id.null[1],]) # Close the polygon
  
  df.hull.null <- df.null[hull.id.null,]  # Don't close the polygon
  
  n <- nrow(df.hull.null) # need the number of vertices in the convex hull dataset
  
  df.hull.null <- df.hull.null[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 
  
  tri.null <- data.frame(
    anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
    
    apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
    baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
    baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
    
    area = numeric(),                                                             # Triangle area
    
    poly.area = numeric(),                                                        # Polygon size for z.yarison
    
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n){
    
    a <- i # select the vertices to cut open the polygon on
    b <- ifelse(i<n, (i + 1), 1)
    
    idx <- a
    k <- a
    
    while (k != b) {
      k <- ifelse(k == 1, n, k - 1)
      idx <- c(idx, k)
    }
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("z.x", "z.y")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("z.x","z.y")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x","z.y")]) # Points to test (calculating distance from E)
    
    AB <- B-A 
    num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
    dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
    
    k <- which.max(dist) # Identify antipode (max distance vertex)
    
    E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
    
    # Split the polygon into left and right sides
    left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
    right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
    
    # Fit data frame of these as edge sequences
    L.edges <- data.frame(
      x1 = left.verts$z.x[-nrow(left.verts)],
      y1 = left.verts$z.y[-nrow(left.verts)],
      x2 = left.verts$z.x[-1],
      y2 = left.verts$z.y[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$z.x[-nrow(right.verts)],
      y1 = right.verts$z.y[-nrow(right.verts)],
      x2 = right.verts$z.x[-1],
      y2 = right.verts$z.y[-1]
    )
    
    t.hat <- AB/ sqrt(sum(AB^2))
    n.hat <- c(-t.hat[2], t.hat[1])
    
    cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
    
    best.area <- Inf
    best <- NULL
    
    for (l in 1:nrow(L.edges)) {
      k <- L.edges[l,]
      
      p <- as.numeric(c(k[1], k[2]))
      q <- as.numeric(c(k[3], k[4]))
      d.L <- q - p
      d.L <- d.L / sqrt(sum(d.L^2))
      
      nL <- c(-d.L[2], d.L[1])
      vals <- df.hull.null$z.x * nL[1] + df.hull.null$z.y * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$z.x * nR[1] + df.hull.null$z.y * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
        
        M <- cbind(d.L, -d.R) # calculate the apex
        if (abs(det(M)) < 1e-12) next  # skip near-parallel
        
        sol <- solve(M, x0.R - x0.L)
        apex <- x0.L + sol[1] * d.L
        
        apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
        if (apex.h <= cE) next
        
        M1 <- cbind(AB, -d.L)
        sol1 <- solve(M1, x0.L - A)
        base.L <- A + sol1[1] * AB
        
        M2 <- cbind(AB, -d.R)
        sol2 <- solve(M2, x0.R - A)
        base.R <- A + sol2[1] * AB
        
        v1 <- base.L - apex
        v2 <- base.R - apex
        area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
        
        if (area < best.area) {
          best.area <- area
          best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
        }
      }
    }
    
    if (is.null(best) || !is.finite(best.area)) next
    
    tri.null <- rbind(tri.null,data.frame(
      anchor.edge = i,                                                            # Which segment is E (leading vertex)
      
      apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
      baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
      baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
      
      area = best.area,                                                           # Triangle area
      poly.area = abs(polyarea(df.hull.full.null$z.x, df.hull.full.null$z.y))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri$poly.area[1]) # 0.030

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.077

###### Pareto fronts ######

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

P.scam <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/P*)", 
       title = "H — Phosphorous") +  # labels
  
  ylim(-0.10, 300) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

P.scam  # Display the plot

p.p3 <- p.p  + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) # Adding scam PF fits
p.p3

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

mean(null.df$a.emp.n >= a.emp) # p-value 0.059

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # -26.45475, p 0.30621
summary(q75, se = "boot", R = 1000) # -99.00829, p 0.03322
summary(q90, se = "boot", R = 1000) # -26.45475, p 0.28363

P.qr <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/P*)", 
       title = "H — Phosphorous") +  # labels
  
  ylim(-0.10, 300) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

P.qr # Display the plot

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%                           # smallest → largest
  slice((floor(0.6667 * n()) + 1):n()) %>%          # keep *second* half
  select(-distance)

P.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/P*)", 
       title = "H — Phosphorous") +  # labels
  
  ylim(-0.10, 300) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

P.scam2  # Display the plot

null.df3 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df3 <- rbind(null.df3, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.001

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -173.47772, p 0.00005
summary(q75, se = "boot", R = 1000) # -192.17816, p 0.00028
summary(q90, se = "boot", R = 1000) # -204.68265, p 0.00001

P.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Competitive ability (1/P*)", 
       title = "H — Phosphorous") +  # labels
  
  ylim(-0.10, 300) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

P.qr2 # Display the plot

# Phosphorous v Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = comp.P,
    z.x = T.br
  ) %>% # Specify the x and y variables
  filter(!is.na(z.y)) %>% 
  filter(!is.na(z.x)) # Only 23 observations

plot(z.y ~ z.x, data= df.filt) # look at the data 
# Not going to work. Point in the upper right

###### Shoval et al triangular analysis ######

hull.id <- chull(df.filt$z.x, df.filt$z.y)

df.hull.full <- df.filt[hull.id,] %>% 
  bind_rows(df.filt[hull.id[1],]) # Close the polygon

ggplot(df.filt, aes(z.x, z.y)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(z.x, z.y, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Find the smallest triangle ######

df.hull <- df.filt[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri <- data.frame(
  anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
  
  apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
  baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
  baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
  
  area = numeric(),                                                             # Triangle area
  
  poly.area = numeric(),                                                        # Polygon size for comparison
  
  stringsAsFactors = FALSE
)


for (i in 1:n){
  
  a <- i # select the vertices to cut open the polygon on
  b <- ifelse(i<n, (i + 1), 1)
  
  idx <- a
  k <- a
  
  while (k != b) {
    k <- ifelse(k == 1, n, k - 1)
    idx <- c(idx, k)
  }
  
  df.cut <- df.hull[idx, ] # OK our loop has succesfully cut open the polygon along the edge between 1 and 2!
  
  # Step 2 - find the antipode E' of E
  
  A <- as.numeric(df.cut[1,c("z.x", "z.y")])                      # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("z.x", "z.y")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x", "z.y")]) # Points to test (calculating distance from E)
  
  AB <- B-A 
  num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
  dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
  
  k <- which.max(dist) # Identify antipode (max distance vertex)
  
  E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
  
  # Split the polygon into left and right sides
  left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
  right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
  
  # Fit data frame of these as edge sequences
  L.edges <- data.frame(
    x1 = left.verts$z.x[-nrow(left.verts)],
    y1 = left.verts$z.y[-nrow(left.verts)],
    x2 = left.verts$z.x[-1],
    y2 = left.verts$z.y[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$z.x[-nrow(right.verts)],
    y1 = right.verts$z.y[-nrow(right.verts)],
    x2 = right.verts$z.x[-1],
    y2 = right.verts$z.y[-1]
  )
  
  t.hat <- AB/ sqrt(sum(AB^2))
  n.hat <- c(-t.hat[2], t.hat[1])
  
  cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
  
  best.area <- Inf
  best <- NULL
  
  for (l in 1:nrow(L.edges)) {
    k <- L.edges[l,]
    
    p <- as.numeric(c(k[1], k[2]))
    q <- as.numeric(c(k[3], k[4]))
    d.L <- q - p
    d.L <- d.L / sqrt(sum(d.L^2))
    
    nL <- c(-d.L[2], d.L[1])
    vals <- df.hull$z.x * nL[1] + df.hull$z.y * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$z.x[idx], df.hull$z.y[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$z.x * nR[1] + df.hull$z.y * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$z.x[idx], df.hull$z.y[idx])
      
      M <- cbind(d.L, -d.R) # calculate the apex
      if (abs(det(M)) < 1e-12) next  # skip near-parallel
      
      sol <- solve(M, x0.R - x0.L)
      apex <- x0.L + sol[1] * d.L
      
      apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
      if (apex.h <= cE) next
      
      M1 <- cbind(AB, -d.L)
      sol1 <- solve(M1, x0.L - A)
      base.L <- A + sol1[1] * AB
      
      M2 <- cbind(AB, -d.R)
      sol2 <- solve(M2, x0.R - A)
      base.R <- A + sol2[1] * AB
      
      v1 <- base.L - apex
      v2 <- base.R - apex
      area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
      
      if (area < best.area) {
        best.area <- area
        best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
      }
    }
  }
  
  tri <- rbind(tri,data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$z.x, df.hull.full$z.y))               # Polygon area
  ))
  
}

best.tri <- tri %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly <- data.frame(
  x = c(
    best.tri$apex.x,
    best.tri$baseL.x,
    best.tri$baseR.x,
    best.tri$apex.x   # close polygon
  ),
  y = c(
    best.tri$apex.y,
    best.tri$baseL.y,
    best.tri$baseR.y,
    best.tri$apex.y
  )
)

p.pt <- ggplot(df.filt, aes(z.x, z.y)) +
  geom_point(size = 2) +
  
  labs(x = "Thermal breadth (°C)", 
       y = "Competitive ability (1/P*)", 
       title = "I — Phosphorous ~ Temperature") +  # labels
  
  geom_polygon(
    data = tri.poly,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  ylim(-0.1, 900) +
  # xlim(0, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.pt

tri

best.tri$area/best.tri$poly.area # 2.861394

###### Randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.filt %>%          # scramble the data
    
    mutate(
      z.x = sample(df.filt$z.x, replace = FALSE),             # Randomly assign rmax
      z.y = sample(df.filt$z.y, replace = FALSE)              # Separately reassign z.y
    )
  
  hull.id.null <- chull(df.null$z.x, df.null$z.y)    # Get the polygon
  
  df.hull.full.null <- df.null[hull.id.null,] %>% 
    bind_rows(df.null[hull.id.null[1],]) # Close the polygon
  
  df.hull.null <- df.null[hull.id.null,]  # Don't close the polygon
  
  n <- nrow(df.hull.null) # need the number of vertices in the convex hull dataset
  
  df.hull.null <- df.hull.null[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 
  
  tri.null <- data.frame(
    anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
    
    apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
    baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
    baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
    
    area = numeric(),                                                             # Triangle area
    
    poly.area = numeric(),                                                        # Polygon size for z.yarison
    
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n){
    
    a <- i # select the vertices to cut open the polygon on
    b <- ifelse(i<n, (i + 1), 1)
    
    idx <- a
    k <- a
    
    while (k != b) {
      k <- ifelse(k == 1, n, k - 1)
      idx <- c(idx, k)
    }
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("z.x", "z.y")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("z.x","z.y")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x","z.y")]) # Points to test (calculating distance from E)
    
    AB <- B-A 
    num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
    dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
    
    k <- which.max(dist) # Identify antipode (max distance vertex)
    
    E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
    
    # Split the polygon into left and right sides
    left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
    right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
    
    # Fit data frame of these as edge sequences
    L.edges <- data.frame(
      x1 = left.verts$z.x[-nrow(left.verts)],
      y1 = left.verts$z.y[-nrow(left.verts)],
      x2 = left.verts$z.x[-1],
      y2 = left.verts$z.y[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$z.x[-nrow(right.verts)],
      y1 = right.verts$z.y[-nrow(right.verts)],
      x2 = right.verts$z.x[-1],
      y2 = right.verts$z.y[-1]
    )
    
    t.hat <- AB/ sqrt(sum(AB^2))
    n.hat <- c(-t.hat[2], t.hat[1])
    
    cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
    
    best.area <- Inf
    best <- NULL
    
    for (l in 1:nrow(L.edges)) {
      k <- L.edges[l,]
      
      p <- as.numeric(c(k[1], k[2]))
      q <- as.numeric(c(k[3], k[4]))
      d.L <- q - p
      d.L <- d.L / sqrt(sum(d.L^2))
      
      nL <- c(-d.L[2], d.L[1])
      vals <- df.hull.null$z.x * nL[1] + df.hull.null$z.y * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$z.x * nR[1] + df.hull.null$z.y * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
        
        M <- cbind(d.L, -d.R) # calculate the apex
        if (abs(det(M)) < 1e-12) next  # skip near-parallel
        
        sol <- solve(M, x0.R - x0.L)
        apex <- x0.L + sol[1] * d.L
        
        apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
        if (apex.h <= cE) next
        
        M1 <- cbind(AB, -d.L)
        sol1 <- solve(M1, x0.L - A)
        base.L <- A + sol1[1] * AB
        
        M2 <- cbind(AB, -d.R)
        sol2 <- solve(M2, x0.R - A)
        base.R <- A + sol2[1] * AB
        
        v1 <- base.L - apex
        v2 <- base.R - apex
        area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
        
        if (area < best.area) {
          best.area <- area
          best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
        }
      }
    }
    
    if (is.null(best) || !is.finite(best.area)) next
    
    tri.null <- rbind(tri.null,data.frame(
      anchor.edge = i,                                                            # Which segment is E (leading vertex)
      
      apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
      baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
      baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
      
      area = best.area,                                                           # Triangle area
      poly.area = abs(polyarea(df.hull.full.null$z.x, df.hull.full.null$z.y))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri$poly.area[1]) # 0.927

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.981

###### Pareto fronts ######

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
fit2 <- lm(z.y ~ z.x, data = par.res.2)

fit <- lm(z.y~z.x, data = par.res.1)

pred.curve.2 <- data.frame( # predicted data frame
  z.x = x.vals,
  z.y = predict(fit2, newdata = data.frame(z.x = x.vals))
)

PT.scam <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  #geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Thermal breadth (°C)", 
       y = "Competitive ability (1/P*)", 
       title = "I — Phosphorous ~ Temperature") +  # labels
  
  ylim(-0.10, 200) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

PT.scam  # Display the plot

p.pt3 <- p.pt # Adding scam PF fits
p.pt3

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

mean(null.df$a.emp.n >= a.emp) # p-value 0.892

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 1.19109, p 0.78554
summary(q75, se = "boot", R = 1000) # 9.12497, p 0.06704
summary(q90, se = "boot", R = 1000) # 5.98178, p 0.54264

PT.qr <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Thermal breadth (°C)", 
       y = "Competitive ability (1/P*)", 
       title = "I — Phosphorous ~ Temperature") +  # labels
  
  ylim(-0.10, 200) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

PT.qr # Display the plot

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%                           # smallest → largest
  slice((floor(0.6667 * n()) + 1):n()) %>%          # keep *second* half
  select(-distance)

PT.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  #geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Thermal breadth (°C)", 
       y = "Competitive ability (1/P*)", 
       title = "I — Phosphorous ~ Temperature") +  # labels
  
  ylim(-0.10, 200) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

PT.scam2  # Display the plot

null.df3 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df3 <- rbind(null.df3, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.688

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -6.17759, p 0.71804
summary(q75, se = "boot", R = 1000) # 4.33302, p 0.73040
summary(q90, se = "boot", R = 1000) # -1.55442, p 0.88801

PT.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Thermal breadth (°C)", 
       y = "Competitive ability (1/P*)", 
       title = "I — Phosphorous ~ Temperature") +  # labels
  
  ylim(-0.10, 200) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

PT.qr2 # Display the plot

# Temperature v Temperature -------------------------------------------------------------------

df.filt <- df %>% 
  mutate(
    z.y = T.br,
    z.x = r.max.T
  ) %>% # Specify the x and y variables
  filter(!is.na(z.y)) %>% 
  filter(!is.na(z.x)) # Only 31 observations

plot(z.y ~ z.x, data= df.filt) # look at the data 

###### Shoval et al triangular analysis ######

hull.id <- chull(df.filt$z.x, df.filt$z.y)

df.hull.full <- df.filt[hull.id,] %>% 
  bind_rows(df.filt[hull.id[1],]) # Close the polygon

ggplot(df.filt, aes(z.x, z.y)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(z.x, z.y, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Find the smallest triangle ######

df.hull <- df.filt[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri <- data.frame(
  anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
  
  apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
  baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
  baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
  
  area = numeric(),                                                             # Triangle area
  
  poly.area = numeric(),                                                        # Polygon size for comparison
  
  stringsAsFactors = FALSE
)


for (i in 1:n){
  
  a <- i # select the vertices to cut open the polygon on
  b <- ifelse(i<n, (i + 1), 1)
  
  idx <- a
  k <- a
  
  while (k != b) {
    k <- ifelse(k == 1, n, k - 1)
    idx <- c(idx, k)
  }
  
  df.cut <- df.hull[idx, ] # OK our loop has succesfully cut open the polygon along the edge between 1 and 2!
  
  # Step 2 - find the antipode E' of E
  
  A <- as.numeric(df.cut[1,c("z.x", "z.y")])                      # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("z.x", "z.y")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x", "z.y")]) # Points to test (calculating distance from E)
  
  AB <- B-A 
  num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
  dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
  
  k <- which.max(dist) # Identify antipode (max distance vertex)
  
  E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
  
  # Split the polygon into left and right sides
  left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
  right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
  
  # Fit data frame of these as edge sequences
  L.edges <- data.frame(
    x1 = left.verts$z.x[-nrow(left.verts)],
    y1 = left.verts$z.y[-nrow(left.verts)],
    x2 = left.verts$z.x[-1],
    y2 = left.verts$z.y[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$z.x[-nrow(right.verts)],
    y1 = right.verts$z.y[-nrow(right.verts)],
    x2 = right.verts$z.x[-1],
    y2 = right.verts$z.y[-1]
  )
  
  t.hat <- AB/ sqrt(sum(AB^2))
  n.hat <- c(-t.hat[2], t.hat[1])
  
  cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
  
  best.area <- Inf
  best <- NULL
  
  for (l in 1:nrow(L.edges)) {
    k <- L.edges[l,]
    
    p <- as.numeric(c(k[1], k[2]))
    q <- as.numeric(c(k[3], k[4]))
    d.L <- q - p
    d.L <- d.L / sqrt(sum(d.L^2))
    
    nL <- c(-d.L[2], d.L[1])
    vals <- df.hull$z.x * nL[1] + df.hull$z.y * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$z.x[idx], df.hull$z.y[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$z.x * nR[1] + df.hull$z.y * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$z.x[idx], df.hull$z.y[idx])
      
      M <- cbind(d.L, -d.R) # calculate the apex
      if (abs(det(M)) < 1e-12) next  # skip near-parallel
      
      sol <- solve(M, x0.R - x0.L)
      apex <- x0.L + sol[1] * d.L
      
      apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
      if (apex.h <= cE) next
      
      M1 <- cbind(AB, -d.L)
      sol1 <- solve(M1, x0.L - A)
      base.L <- A + sol1[1] * AB
      
      M2 <- cbind(AB, -d.R)
      sol2 <- solve(M2, x0.R - A)
      base.R <- A + sol2[1] * AB
      
      v1 <- base.L - apex
      v2 <- base.R - apex
      area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
      
      if (area < best.area) {
        best.area <- area
        best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
      }
    }
  }
  
  tri <- rbind(tri,data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$z.x, df.hull.full$z.y))               # Polygon area
  ))
  
}

best.tri <- tri %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly <- data.frame(
  x = c(
    best.tri$apex.x,
    best.tri$baseL.x,
    best.tri$baseR.x,
    best.tri$apex.x   # close polygon
  ),
  y = c(
    best.tri$apex.y,
    best.tri$baseL.y,
    best.tri$baseR.y,
    best.tri$apex.y
  )
)

p.t <- ggplot(df.filt, aes(z.x, z.y)) +
  geom_point(size = 2) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Thermal breadth (°C)", 
       title = "J — Temperature") +  # labels
  
  geom_polygon(
    data = tri.poly,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  ylim(0, 37.5) +
  # xlim(0, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.t

tri

best.tri$area/best.tri$poly.area # 1.227325

###### Randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.filt %>%          # scramble the data
    
    mutate(
      z.x = sample(df.filt$z.x, replace = FALSE),             # Randomly assign rmax
      z.y = sample(df.filt$z.y, replace = FALSE)              # Separately reassign z.y
    )
  
  hull.id.null <- chull(df.null$z.x, df.null$z.y)    # Get the polygon
  
  df.hull.full.null <- df.null[hull.id.null,] %>% 
    bind_rows(df.null[hull.id.null[1],]) # Close the polygon
  
  df.hull.null <- df.null[hull.id.null,]  # Don't close the polygon
  
  n <- nrow(df.hull.null) # need the number of vertices in the convex hull dataset
  
  df.hull.null <- df.hull.null[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 
  
  tri.null <- data.frame(
    anchor.edge = integer(),                                                      # Which segment is E (leading vertex)
    
    apex.x = numeric(), apex.y = numeric(),                                       # Apex coordinates
    baseL.x = numeric(), baseL.y = numeric(),                                     # Left base coordinates
    baseR.x = numeric(), baseR.y = numeric(),                                     # Right base coordinates
    
    area = numeric(),                                                             # Triangle area
    
    poly.area = numeric(),                                                        # Polygon size for z.yarison
    
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n){
    
    a <- i # select the vertices to cut open the polygon on
    b <- ifelse(i<n, (i + 1), 1)
    
    idx <- a
    k <- a
    
    while (k != b) {
      k <- ifelse(k == 1, n, k - 1)
      idx <- c(idx, k)
    }
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("z.x", "z.y")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("z.x","z.y")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("z.x","z.y")]) # Points to test (calculating distance from E)
    
    AB <- B-A 
    num <- AB[1] * (P[,2] - A[2]) - AB[2] * (P[,1] - A[1])   # 2D cross product scalar
    dist <- abs(num) / sqrt(sum(AB^2))                       # perpendicular distance
    
    k <- which.max(dist) # Identify antipode (max distance vertex)
    
    E.vertex <- df.cut[1 + k, ]   # +1 because P starts at row 2 of df.cut12
    
    # Split the polygon into left and right sides
    left.verts  <- df.cut[1:(k+1), ]                         # v1 ... E'
    right.verts  <- df.cut[(k+1):nrow(df.cut), ]             # E' ... v2
    
    # Fit data frame of these as edge sequences
    L.edges <- data.frame(
      x1 = left.verts$z.x[-nrow(left.verts)],
      y1 = left.verts$z.y[-nrow(left.verts)],
      x2 = left.verts$z.x[-1],
      y2 = left.verts$z.y[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$z.x[-nrow(right.verts)],
      y1 = right.verts$z.y[-nrow(right.verts)],
      x2 = right.verts$z.x[-1],
      y2 = right.verts$z.y[-1]
    )
    
    t.hat <- AB/ sqrt(sum(AB^2))
    n.hat <- c(-t.hat[2], t.hat[1])
    
    cE <- A[1]*n.hat[1] + A[2]*n.hat[2]
    
    best.area <- Inf
    best <- NULL
    
    for (l in 1:nrow(L.edges)) {
      k <- L.edges[l,]
      
      p <- as.numeric(c(k[1], k[2]))
      q <- as.numeric(c(k[3], k[4]))
      d.L <- q - p
      d.L <- d.L / sqrt(sum(d.L^2))
      
      nL <- c(-d.L[2], d.L[1])
      vals <- df.hull.null$z.x * nL[1] + df.hull.null$z.y * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$z.x * nR[1] + df.hull.null$z.y * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$z.x[idx], df.hull.null$z.y[idx])
        
        M <- cbind(d.L, -d.R) # calculate the apex
        if (abs(det(M)) < 1e-12) next  # skip near-parallel
        
        sol <- solve(M, x0.R - x0.L)
        apex <- x0.L + sol[1] * d.L
        
        apex.h <- apex[1]*n.hat[1] + apex[2]*n.hat[2] # must be on polygon side of base
        if (apex.h <= cE) next
        
        M1 <- cbind(AB, -d.L)
        sol1 <- solve(M1, x0.L - A)
        base.L <- A + sol1[1] * AB
        
        M2 <- cbind(AB, -d.R)
        sol2 <- solve(M2, x0.R - A)
        base.R <- A + sol2[1] * AB
        
        v1 <- base.L - apex
        v2 <- base.R - apex
        area <- 0.5 * abs(v1[1]*v2[2] - v1[2]*v2[1])
        
        if (area < best.area) {
          best.area <- area
          best <- list(l = l, r = r, apex = apex, base.L = base.L, base.R = base.R)
        }
      }
    }
    
    if (is.null(best) || !is.finite(best.area)) next
    
    tri.null <- rbind(tri.null,data.frame(
      anchor.edge = i,                                                            # Which segment is E (leading vertex)
      
      apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
      baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
      baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
      
      area = best.area,                                                           # Triangle area
      poly.area = abs(polyarea(df.hull.full.null$z.x, df.hull.full.null$z.y))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri$poly.area[1]) # 0.583

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.396

###### Pareto fronts ######

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

T.scam <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Thermal breadth (°C)", 
       title = "J — Temperature") +  # labels
  
  ylim(-0.10, 37.5) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.scam  # Display the plot

p.t3 <- p.t  + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) # Adding scam PF fits
p.t3

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

mean(null.df$a.emp.n >= a.emp) # p-value 0.479

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 1.93522, p 0.03919
summary(q75, se = "boot", R = 1000) # 1.96605, p 0.10820
summary(q90, se = "boot", R = 1000) # 2.55748, p 0.21289

T.qr <- ggplot(df.filt, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Thermal breadth (°C)", 
       title = "J — Temperature") +  # labels
  
  ylim(-0.10, 37.5) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.qr # Display the plot

###### top 33% of data ######

df.filt3 <- df.filt %>%
  arrange(distance) %>%                           # smallest → largest
  slice((floor(0.6667 * n()) + 1):n()) %>%          # keep *second* half
  select(-distance)

T.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Thermal breadth (°C)", 
       title = "J — Temperature") +  # labels
  
  ylim(-0.10, 37.5) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.scam2  # Display the plot

null.df3 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df <- df.filt3 %>%
    
    mutate(
      z.x.sim = sample(df.filt3$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt3$z.y, replace = FALSE),     # Separately reassign y
    )
  
  par.res.n <- par_frt(shuffled.df, xvar = "z.x.sim", yvar = "z.y.sim") #  Pareto front on shuffled data
  
  par.res.n <- par.res.n %>% # Arrange the set of Pareto-optimal points by increasing values of z.x
    arrange(z.x.sim)
  
  x.max.n <- max(shuffled.df$z.x.sim) # Extract the max values for x and y
  y.max.n <- max(shuffled.df$z.y.sim)
  
  poly.n <- par.res.n[, c("z.x.sim", "z.y.sim")] %>% # Create a dataframe with the pareto-optimal points and the maximum value 
    add_row(z.x.sim = x.max.n, z.y.sim = y.max.n)
  
  a.emp.n <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim) # Calculate the area enclosed by these vertices
  
  null.df3 <- rbind(null.df3, data.frame(  # Save the data
    a.emp.n = a.emp.n,                   # Area above the curve 
    n.PF = nrow(par.res.n)              # Number of data points in the PF
  ))
  
}

mean(null.df3$a.emp.n >= a.emp) # p-value 0.021

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -2.76484, p 0.05569
summary(q75, se = "boot", R = 1000) # -3.69575, p 0.03038
summary(q90, se = "boot", R = 1000) # -3.15868, p 0.06728

T.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)", 
       y = "Thermal breadth (°C)", 
       title = "J — Temperature") +  # labels
  
  ylim(-0.10, 37.5) +
  #xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.qr2 # Display the plot

# Assemble the plots ------------------------------------------------------

meta_cross_PF <- plot_grid(p.l3, p.ln3, p.lp3, p.lt3,
                           NULL, p.n3, p.np3, p.nt3,
                           NULL, NULL, p.p3, p.pt3,
                           NULL, NULL, NULL, p.t3,
                     ncol = 4,
                     align = "hv",
                     axis = "tblr")

ggsave("figures/110_fig_extra_synthesis_crosses.jpeg", meta_cross_PF, width = 12, height = 12)
