# Jason R Laurich

# January 13th, 2025

# I am going to analyze interspecific pareto fronts here, working with the temperature, phosphorous, light, and nitrogen data 
# I have compiled. 

# Going to base my analysis on shoval 2012 (science)

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

# Light -------------------------------------------------------------------

# Edwards 2016

df.l.e <- read.csv("data-processed/502a_Edwards_2016_light_monods.csv")
head(df.l.e)

df.l.e$dataset <- "Edwards et al., 2016"

# Lewington-Pearce 2019

df.l.l <- read.csv("data-processed/502c_Lewington_2019_light_monods.csv")
head(df.l.l)

df.l.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.l.lv <- read.csv("data-processed/502e_Levasseur_2025_light_monods.csv")
head(df.l.lv)

df.l.lv$dataset <- "Levasseur et al., 2025"

# Narwani 2015

df.l.n <- read.csv("data-processed/13_Narwani_2015_summary.csv")
head(df.l.n)

df.l.n$dataset <- "Narwani et al., 2015"

df.l.n <- df.l.n %>% 
  filter(!is.na(Istar))

df.l <- bind_rows(
  
  df.l.e %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.l.l %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.l.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.l.n %>% 
    transmute(Species = Species.name,
              dataset = dataset,
              r.max = umax.light,
              R = Istar,
              comp = 1/Istar)
  
)

head(df.l) # 131 entries!

df.l <- df.l %>% 
  filter(R > 0) # 130 viable entries

###### Shoval analysis ######

plot(df.l$comp ~ df.l$r.max) # Some really wonky (high) estimates of 1/R* that suggest error (filter at > 3.5)

df.l <- df.l %>% 
  filter(comp < 3.5,
         r.max < 2.5,
         !is.na(comp))

write.csv(df.l, "data-processed/504a_light_metadata.csv") # Summary

hull.id <- chull(df.l$r.max, df.l$comp)

df.hull.full <- df.l[hull.id,] %>% 
  bind_rows(df.l[hull.id[1],]) # Close the polygon

ggplot(df.l, aes(r.max, comp)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(r.max, comp, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Clean code to run in loop ######

df.hull <- df.l[hull.id,]  # Don't close the polygon

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
  
  A <- as.numeric(df.cut[1,c("r.max", "comp")])                     # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("r.max","comp")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("r.max","comp")]) # Points to test (calculating distance from E)
  
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
    x1 = left.verts$r.max[-nrow(left.verts)],
    y1 = left.verts$comp[-nrow(left.verts)],
    x2 = left.verts$r.max[-1],
    y2 = left.verts$comp[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$r.max[-nrow(right.verts)],
    y1 = right.verts$comp[-nrow(right.verts)],
    x2 = right.verts$r.max[-1],
    y2 = right.verts$comp[-1]
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
    vals <- df.hull$r.max * nL[1] + df.hull$comp * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$r.max[idx], df.hull$comp[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$r.max * nR[1] + df.hull$comp * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$r.max[idx], df.hull$comp[idx])
      
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
    poly.area = abs(polyarea(df.hull.full$r.max, df.hull.full$comp))            # Polygon area
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

p.l <- ggplot(df.l, aes(r.max, comp, colour = dataset)) +
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
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
  ylim(-0.10, 3.25) +
  xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.l

tri

best.tri$area/best.tri$poly.area # 1.09819

###### Null randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.l %>%          # scramble the data
    
    mutate(
      r.max = sample(df.l$r.max, replace = FALSE),           # Randomly assign rmax
      comp = sample(df.l$comp, replace = FALSE)              # Separately reassign comp
    )
  
  hull.id.null <- chull(df.null$r.max, df.null$comp)    # Get the polygon
  
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
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("r.max", "comp")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("r.max","comp")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("r.max","comp")]) # Points to test (calculating distance from E)
    
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
      x1 = left.verts$r.max[-nrow(left.verts)],
      y1 = left.verts$comp[-nrow(left.verts)],
      x2 = left.verts$r.max[-1],
      y2 = left.verts$comp[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$r.max[-nrow(right.verts)],
      y1 = right.verts$comp[-nrow(right.verts)],
      x2 = right.verts$r.max[-1],
      y2 = right.verts$comp[-1]
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
      vals <- df.hull.null$r.max * nL[1] + df.hull.null$comp * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$r.max[idx], df.hull.null$comp[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$r.max * nR[1] + df.hull.null$comp * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$r.max[idx], df.hull.null$comp[idx])
        
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
      poly.area = abs(polyarea(df.hull.full.null$r.max, df.hull.full.null$comp))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri$poly.area[1])
# For the comp <10 df: 0.418
# For the heavily restricted df: 0.003

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) 
# For the comp <10 df: 0.436
# For the heavily restricted df: 0.010

# Nitrogen -------------------------------------------------------------------

# Edwards 2015

df.n.e <- read.csv("data-processed/14_Edwards_2015_summary.csv")
head(df.n.e)

df.n.e <- df.n.e %>% 
  filter(!is.na(k_nit_m),
         !is.na(mu_nit)) %>% 
  mutate(R = 0.1*k_nit_m/(mu_nit - 0.1))

df.n.e$dataset <- "Edwards et al., 2015"

# Lewington-Pearce 2019

df.n.l <- read.csv("data-processed/503a_Lewington_2019_nit_monods.csv")
head(df.n.l)

df.n.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.n.lv <- read.csv("data-processed/503c_Levasseur_2025_nit_monods.csv")
head(df.n.lv)

df.n.lv$dataset <- "Levasseur et al., 2025"

df.n <- bind_rows(
  
  df.n.e %>% 
    transmute(Species = species,
              dataset = dataset,
              r.max = mu_nit,
              R = R,
              comp = 1/R),
  
  df.n.l %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.n.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post)
  
)

head(df.n) # 58 entries!

df.n <- df.n %>% 
  filter(R > 0) # 57 viable entries

###### Shoval analysis ######

plot(df.n$comp ~ df.n$r.max) # Some really wonky (high) estimates of 1/R* that suggest error or nitrogen fixation?

df.n <- df.n %>% 
  filter(comp < 60,
         !is.na(r.max),
         !is.na(comp))

write.csv(df.n, "data-processed/504b_nit_metadata.csv") # Summary

hull.id <- chull(df.n$r.max, df.n$comp)

df.hull.full <- df.n[hull.id,] %>% 
  bind_rows(df.n[hull.id[1],]) # Close the polygon

ggplot(df.n, aes(r.max, comp)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(r.max, comp, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Clean code to run in loop ######

df.hull <- df.n[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri.n <- data.frame(
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
  
  A <- as.numeric(df.cut[1,c("r.max", "comp")])                     # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("r.max","comp")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("r.max","comp")]) # Points to test (calculating distance from E)
  
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
    x1 = left.verts$r.max[-nrow(left.verts)],
    y1 = left.verts$comp[-nrow(left.verts)],
    x2 = left.verts$r.max[-1],
    y2 = left.verts$comp[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$r.max[-nrow(right.verts)],
    y1 = right.verts$comp[-nrow(right.verts)],
    x2 = right.verts$r.max[-1],
    y2 = right.verts$comp[-1]
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
    vals <- df.hull$r.max * nL[1] + df.hull$comp * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$r.max[idx], df.hull$comp[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$r.max * nR[1] + df.hull$comp * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$r.max[idx], df.hull$comp[idx])
      
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
  
  tri.n <- rbind(tri.n, data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$r.max, df.hull.full$comp))            # Polygon area
  ))
  
}

best.tri.n <- tri.n %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly.n <- data.frame(
  x = c(
    best.tri.n$apex.x,
    best.tri.n$baseL.x,
    best.tri.n$baseR.x,
    best.tri.n$apex.x   # close polygon
  ),
  y = c(
    best.tri.n$apex.y,
    best.tri.n$baseL.y,
    best.tri.n$baseR.y,
    best.tri.n$apex.y
  )
)

p.n <- ggplot(df.n, aes(r.max, comp, colour = dataset)) +
  geom_point(size = 2) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/N*)", 
       title = "B — Nitrogen") +  # labels
  
  geom_polygon(
    data = tri.poly.n,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2015" = "chocolate1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1")
  ) +
  
  ylim(0,52) +
  xlim(0,4.5) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.n

tri.n

best.tri.n$area/best.tri.n$poly.area # 1.358692

###### Null randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.n %>%          # scramble the data
    
    mutate(
      r.max = sample(df.n$r.max, replace = FALSE),           # Randomly assign rmax
      comp = sample(df.n$comp, replace = FALSE)              # Separately reassign comp
    )
  
  hull.id.null <- chull(df.null$r.max, df.null$comp)    # Get the polygon
  
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
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("r.max", "comp")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("r.max","comp")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("r.max","comp")]) # Points to test (calculating distance from E)
    
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
      x1 = left.verts$r.max[-nrow(left.verts)],
      y1 = left.verts$comp[-nrow(left.verts)],
      x2 = left.verts$r.max[-1],
      y2 = left.verts$comp[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$r.max[-nrow(right.verts)],
      y1 = right.verts$comp[-nrow(right.verts)],
      x2 = right.verts$r.max[-1],
      y2 = right.verts$comp[-1]
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
      vals <- df.hull.null$r.max * nL[1] + df.hull.null$comp * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$r.max[idx], df.hull.null$comp[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$r.max * nR[1] + df.hull.null$comp * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$r.max[idx], df.hull.null$comp[idx])
        
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
      poly.area = abs(polyarea(df.hull.full.null$r.max, df.hull.full.null$comp))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri.n$poly.area[1])
# For the comp <10 df: 
# For the heavily restricted df: 0.143

mean(null.df$best.tri.area/null.df$poly.area <= min(tri.n$area)/tri.n$poly.area[1]) 
# For the comp <10 df: 
# For the heavily restricted df: 0.595

# Phosphorous -------------------------------------------------------------------

# Bestion 2018

df.p.b <- read.csv("data-processed/503a_Bestion_2018_phos_monods.csv")
head(df.p.b)

df.p.b$dataset <- "Bestion et al., 2018"

# Edwards 2015

df.p.e <- read.csv("data-processed/14_Edwards_2015_summary.csv")
head(df.p.e)

df.p.e <- df.p.e %>% 
  filter(!is.na(k_p_m),
         !is.na(mu_p)) %>% 
  mutate(R = 0.1*k_p_m/(mu_p - 0.1))

df.p.e$dataset <- "Edwards et al., 2015"

# Levasseur 2025

df.p.lv <- read.csv("data-processed/503c_Levasseur_2025_phos_monods.csv")
head(df.p.lv)

df.p.lv$dataset <- "Levasseur et al., 2025"

# Narwani 2015

df.p.n <- read.csv("data-processed/13_Narwani_2015_summary.csv")
head(df.p.n)

df.p.n$dataset <- "Narwani et al., 2015"

df.p <- bind_rows(
  
  df.p.b %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.p.e %>% 
    transmute(Species = species,
              dataset = dataset,
              r.max = mu_p,
              R = R,
              comp = 1/R),
  
  df.p.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max.post,
              R = R.post,
              comp = 1/R.post),
  
  df.p.n %>% 
    transmute(Species = Species.name,
              dataset = dataset,
              r.max = umax.nitrate,
              R = Nstar,
              comp = 1/Nstar)
  
)

head(df.p) # 156 entries!

df.p <- df.p %>% 
  filter(R > 0) # 124 viable entries (NAs in the narwani data set)

###### Shoval analysis ######

plot(df.p$comp ~ df.p$r.max) # Some really wonky (high) estimates of 1/R* that suggest error or nitrogen fixation?

df.p <- df.p %>% 
  filter(comp < 300,
         r.max < 2.5,
         !is.na(r.max),
         !is.na(comp))

write.csv(df.p, "data-processed/504c_phos_metadata.csv") # Summary

hull.id <- chull(df.p$r.max, df.p$comp)

df.hull.full <- df.p[hull.id,] %>% 
  bind_rows(df.p[hull.id[1],]) # Close the polygon

ggplot(df.p, aes(r.max, comp)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(r.max, comp, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Clean code to run in loop ######

df.hull <- df.p[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri.p <- data.frame(
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
  
  A <- as.numeric(df.cut[1,c("r.max", "comp")])                     # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("r.max","comp")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("r.max","comp")]) # Points to test (calculating distance from E)
  
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
    x1 = left.verts$r.max[-nrow(left.verts)],
    y1 = left.verts$comp[-nrow(left.verts)],
    x2 = left.verts$r.max[-1],
    y2 = left.verts$comp[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$r.max[-nrow(right.verts)],
    y1 = right.verts$comp[-nrow(right.verts)],
    x2 = right.verts$r.max[-1],
    y2 = right.verts$comp[-1]
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
    vals <- df.hull$r.max * nL[1] + df.hull$comp * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$r.max[idx], df.hull$comp[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$r.max * nR[1] + df.hull$comp * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$r.max[idx], df.hull$comp[idx])
      
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
  
  tri.p <- rbind(tri.p, data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$r.max, df.hull.full$comp))            # Polygon area
  ))
  
}

best.tri.p <- tri.p %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly.p <- data.frame(
  x = c(
    best.tri.p$apex.x,
    best.tri.p$baseL.x,
    best.tri.p$baseR.x,
    best.tri.p$apex.x   # close polygon
  ),
  y = c(
    best.tri.p$apex.y,
    best.tri.p$baseL.y,
    best.tri.p$baseR.y,
    best.tri.p$apex.y
  )
)

p.p <- ggplot(df.p, aes(r.max, comp, colour = dataset)) +
  geom_point(size = 2) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       title = "C — Phosphorous") +  # labels
  
  geom_polygon(
    data = tri.poly.p,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2015" = "chocolate1",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
  ylim(-0.1,340) +
  xlim(0,2.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.p

tri.p

best.tri.p$area/best.tri.p$poly.area # 1.067351

###### Null randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.p %>%          # scramble the data
    
    mutate(
      r.max = sample(df.p$r.max, replace = FALSE),           # Randomly assign rmax
      comp = sample(df.p$comp, replace = FALSE)              # Separately reassign comp
    )
  
  hull.id.null <- chull(df.null$r.max, df.null$comp)    # Get the polygon
  
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
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("r.max", "comp")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("r.max","comp")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("r.max","comp")]) # Points to test (calculating distance from E)
    
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
      x1 = left.verts$r.max[-nrow(left.verts)],
      y1 = left.verts$comp[-nrow(left.verts)],
      x2 = left.verts$r.max[-1],
      y2 = left.verts$comp[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$r.max[-nrow(right.verts)],
      y1 = right.verts$comp[-nrow(right.verts)],
      x2 = right.verts$r.max[-1],
      y2 = right.verts$comp[-1]
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
      vals <- df.hull.null$r.max * nL[1] + df.hull.null$comp * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$r.max[idx], df.hull.null$comp[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$r.max * nR[1] + df.hull.null$comp * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$r.max[idx], df.hull.null$comp[idx])
        
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
      poly.area = abs(polyarea(df.hull.full.null$r.max, df.hull.full.null$comp))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri.p$poly.area[1])
# The restricted one is the reasonable one here: 0.001

mean(null.df$best.tri.area/null.df$poly.area <= min(tri.p$area)/tri.p$poly.area[1]) 
# The restricted one is the reasonable one here: 0.003

# Temperature -------------------------------------------------------------

# Bestion 2018

df.t.b <- read.csv("data-processed/501c_Bestion_2018_TPCs.csv")
head(df.t.b)

df.t.b$dataset <- "Bestion et al., 2018"

# Edwards 2016

df.t.e <- read.csv("data-processed/501g_Edwards_2016_TPCs.csv")
head(df.t.e)

df.t.e$dataset <- "Edwards et al., 2016"

# Lewington-Pearce 2019

df.t.l <- read.csv("data-processed/501e_Lewington_2019_TPCs.csv")
head(df.t.l)

df.t.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.t.lv <- read.csv("data-processed/501i_Levasseur_2025_TPCs.csv")
head(df.t.lv)

df.t.lv$dataset <- "Levasseur et al., 2025"

# Thomas 2012

df.t.t <- read.csv("data-processed/501a_Thomas_2012_TPCs.csv")
head(df.t.t)

df.t.t$dataset <- "Thomas et al., 2012"

df.t <- bind_rows(
  
  df.t.b %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.e %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.l %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.lv %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
  df.t.t %>% 
    transmute(Species = Sp.name,
              dataset = dataset,
              r.max = r.max,
              T.br = T.br.max - (pmax(-1.8, T.br.min))),
  
)

head(df.t) # 229 entries!

###### Shoval analysis ######

plot(df.t$T.br ~ df.t$r.max) # Some really wonky (high) estimates of 1/R* that suggest error or nitrogen fixation?

df.t <- df.t %>% 
  filter(r.max < 3.5,
         T.br < 36,
         !is.na(r.max),
         !is.na(T.br))

write.csv(df.t, "data-processed/504d_temp_metadata.csv") # Summary

hull.id <- chull(df.t$r.max, df.t$T.br)

df.hull.full <- df.t[hull.id,] %>% 
  bind_rows(df.t[hull.id[1],]) # Close the polygon

ggplot(df.t, aes(r.max, T.br)) +
  geom_point() +
  geom_polygon(
    data = df.hull.full,
    aes(r.max, T.br, group = 1),
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  
  theme_classic()

###### Clean code to run in loop ######

df.hull <- df.t[hull.id,]  # Don't close the polygon

n <- nrow(df.hull) # need the number of vertices in the convex hull dataset

df.hull <- df.hull[n:1, ] # we're going to reverse the order of the polygon so that it runs CCW. 

tri.t <- data.frame(
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
  
  A <- as.numeric(df.cut[1,c("r.max", "T.br")])                     # Edge vertex a
  B <- as.numeric(df.cut[nrow(df.cut), c("r.max","T.br")])          # Edge vertex b
  
  P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("r.max","T.br")]) # Points to test (calculating distance from E)
  
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
    x1 = left.verts$r.max[-nrow(left.verts)],
    y1 = left.verts$T.br[-nrow(left.verts)],
    x2 = left.verts$r.max[-1],
    y2 = left.verts$T.br[-1]
  )
  
  R.edges <- data.frame(
    x1 = right.verts$r.max[-nrow(right.verts)],
    y1 = right.verts$T.br[-nrow(right.verts)],
    x2 = right.verts$r.max[-1],
    y2 = right.verts$T.br[-1]
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
    vals <- df.hull$r.max * nL[1] + df.hull$T.br * nL[2]
    idx <- which.max(vals)
    x0.L <- c(df.hull$r.max[idx], df.hull$T.br[idx])
    
    for (r in 1:nrow(R.edges)) {
      m <- R.edges[r,]
      
      p <- as.numeric(c(m[1], m[2]))
      q <- as.numeric(c(m[3], m[4]))
      d.R <- q - p
      d.R <- d.R / sqrt(sum(d.R^2))
      
      nR <- c(-d.R[2], d.R[1])
      vals <- df.hull$r.max * nR[1] + df.hull$T.br * nR[2]
      idx <- which.max(vals)
      x0.R <- c(df.hull$r.max[idx], df.hull$T.br[idx])
      
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
  
  tri.t <- rbind(tri.t, data.frame(
    anchor.edge = i,                                                            # Which segment is E (leading vertex)
    
    apex.x = best$apex[1], apex.y = best$apex[2],                               # Apex coordinates
    baseL.x = best$base.L[1], baseL.y = best$base.L[2],                         # Left base coordinates
    baseR.x = best$base.R[1], baseR.y = best$base.R[2],                         # Right base coordinates
    
    area = best.area,                                                           # Triangle area
    poly.area = abs(polyarea(df.hull.full$r.max, df.hull.full$T.br))            # Polygon area
  ))
  
}

best.tri.t <- tri.t %>% 
  filter(area == min(area)) %>% 
  filter(anchor.edge == min(anchor.edge))


tri.poly.t <- data.frame(
  x = c(
    best.tri.t$apex.x,
    best.tri.t$baseL.x,
    best.tri.t$baseR.x,
    best.tri.t$apex.x   # close polygon
  ),
  y = c(
    best.tri.t$apex.y,
    best.tri.t$baseL.y,
    best.tri.t$baseR.y,
    best.tri.t$apex.y
  )
)

p.t <- ggplot(df.t, aes(r.max, T.br, colour = dataset)) +
  geom_point(size = 2) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       title = "D — Temperature") +  # labels
  
  geom_polygon(
    data = tri.poly.t,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Thomas et al., 2012" = "magenta2")
  ) +
  
  ylim(-1, 47) +
  xlim(-0.05,3.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.t

tri.t

best.tri.t$area/best.tri.t$poly.area # 1.216089

###### Null randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.t %>%          # scramble the data
    
    mutate(
      r.max = sample(df.t$r.max, replace = FALSE),           # Randomly assign rmax
      T.br = sample(df.t$T.br, replace = FALSE)              # Separately reassign T.br
    )
  
  hull.id.null <- chull(df.null$r.max, df.null$T.br)    # Get the polygon
  
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
    
    df.cut <- df.hull.null[idx, ] # OK cut open the polygon
    
    A <- as.numeric(df.cut[1,c("r.max", "T.br")])                     # Edge vertex a
    B <- as.numeric(df.cut[nrow(df.cut), c("r.max","T.br")])          # Edge vertex b
    
    P <- as.data.frame(df.cut[2:(nrow(df.cut)-1), c("r.max","T.br")]) # Points to test (calculating distance from E)
    
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
      x1 = left.verts$r.max[-nrow(left.verts)],
      y1 = left.verts$T.br[-nrow(left.verts)],
      x2 = left.verts$r.max[-1],
      y2 = left.verts$T.br[-1]
    )
    
    R.edges <- data.frame(
      x1 = right.verts$r.max[-nrow(right.verts)],
      y1 = right.verts$T.br[-nrow(right.verts)],
      x2 = right.verts$r.max[-1],
      y2 = right.verts$T.br[-1]
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
      vals <- df.hull.null$r.max * nL[1] + df.hull.null$T.br * nL[2]
      idx <- which.max(vals)
      x0.L <- c(df.hull.null$r.max[idx], df.hull.null$T.br[idx])
      
      for (r in 1:nrow(R.edges)) {
        m <- R.edges[r,]
        
        p <- as.numeric(c(m[1], m[2]))
        q <- as.numeric(c(m[3], m[4]))
        d.R <- q - p
        d.R <- d.R / sqrt(sum(d.R^2))
        
        nR <- c(-d.R[2], d.R[1])
        vals <- df.hull.null$r.max * nR[1] + df.hull.null$T.br * nR[2]
        idx <- which.max(vals)
        x0.R <- c(df.hull.null$r.max[idx], df.hull.null$T.br[idx])
        
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
      poly.area = abs(polyarea(df.hull.full.null$r.max, df.hull.full.null$T.br))  # Polygon area
    ))
    
  }
  
  null.df <- rbind(null.df, data.frame(               # Null model results
    poly.area = tri.null$poly.area[1],                # Polygon area
    best.tri.area = min(tri.null$area),            # Area of the smallest triangle
    stringsAsFactors = FALSE            
  ))
  
}

mean(null.df$poly.area <= tri.t$poly.area[1])
# The restricted one is the reasonable one here: 0.181

mean(null.df$best.tri.area/null.df$poly.area <= min(tri.t$area)/tri.t$poly.area[1]) 
# The restricted one is the reasonable one here: 0.125

# Assemble the Shoval plots -----------------------------------------------

df.t.leg <- df.t %>% 
  rename(R = T.br) %>% 
  mutate(comp = 1)

df.leg <- rbind(df.l, df.n, df.p, df.t.leg)

legend_plot <- ggplot(df.leg, aes(x = r.max, y = comp)) +
  geom_point(aes(colour = dataset), size = 3, stroke = 1.5) +
  
  scale_color_manual(name = "Dataset",
    values = c("Bestion et al., 2018" = "firebrick2",
               "Edwards et al., 2015"  = "chocolate1",
               "Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3",
               "Thomas et al., 2012" = "magenta2")
  ) +
  
  theme_void() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "lines")
  )


legend_only <- get_legend(legend_plot)

shoval_plots <- plot_grid(p.l, p.n, p.p, p.t, legend_only,
                        ncol = 2,
                        align = "hv")

ggsave("figures/109_fig_3_shoval_meta.jpeg", shoval_plots, width = 8, height = 12)

###### Now with my data added ######

df <- read.csv("data-processed/304_summary_table_final.csv") # Summary file
head(df)

p.l2 <- ggplot(df.l, aes(r.max, comp, colour = dataset)) +
  geom_point(size = 2) +
  
  geom_point(
    data = df,
    aes(x = I.µ.max, y = I.comp),
    size = 2,
    colour = "black",
    inherit.aes = FALSE
  ) +
  
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
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
  ylim(-0.10, 3.25) +
  xlim(-0.25, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.l2

p.n2 <- ggplot(df.n, aes(r.max, comp, colour = dataset)) +
  geom_point(size = 2) +
  
  geom_point(
    data = df,
    aes(x = N.µ.max, y = N.comp),
    size = 2,
    colour = "black",
    inherit.aes = FALSE
  ) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/N*)", 
       title = "B — Nitrogen") +  # labels
  
  geom_polygon(
    data = tri.poly.n,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2015" = "chocolate1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1")
  ) +
  
  ylim(0,52) +
  xlim(0,4.5) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.n2

p.p2 <- ggplot(df.p, aes(r.max, comp, colour = dataset)) +
  geom_point(size = 2) +
  
  geom_point(
    data = df,
    aes(x = P.µ.max, y = P.comp),
    size = 2,
    colour = "black",
    inherit.aes = FALSE
  ) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       title = "C — Phosphorous") +  # labels
  
  geom_polygon(
    data = tri.poly.p,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2015" = "chocolate1",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
  ylim(-0.1,340) +
  xlim(0,2.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.p2

p.t2 <- ggplot(df.t, aes(r.max, T.br, colour = dataset)) +
  geom_point(size = 2) +
  
  geom_point(
    data = df,
    aes(x = T.µ.max, y = T.br),
    size = 2,
    colour = "black",
    inherit.aes = FALSE
  ) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       title = "D — Temperature") +  # labels
  
  geom_polygon(
    data = tri.poly.t,
    aes(x, y),
    fill = "grey60",
    alpha = 0.3,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Thomas et al., 2012" = "magenta2")
  ) +
  
  ylim(-1, 47) +
  xlim(-0.05,5) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  ) 

p.t2

df.leg <- df %>% 
  transmute(Species = "na",
            dataset = "Laurich et al., 2026", 
            r.max = 1,
            R = 1,
            comp = 1)

df.leg2 <- rbind(df.l, df.n, df.p, df.t.leg, df.leg)

legend_plot2 <- ggplot(df.leg2, aes(x = r.max, y = comp)) +
  geom_point(aes(colour = dataset), size = 3, stroke = 1.5) +
  
  scale_color_manual(name = "Dataset",
                     values = c("Bestion et al., 2018" = "firebrick2",
                                "Edwards et al., 2015"  = "chocolate1",
                                "Edwards et al., 2016" = "goldenrod1",
                                "Lewington-Pearce et al., 2019" = "forestgreen",
                                "Levasseur et al., 2025" = "royalblue1",
                                "Narwani et al., 2015" = "darkorchid3",
                                "Thomas et al., 2012" = "magenta2",
                                "Laurich et al., 2026" = "black")
  ) +
  
  theme_void() +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1.2, "lines")
  )

legend_only2 <- get_legend(legend_plot2)

shoval_plots2 <- plot_grid(p.l2, p.n2, p.p2, p.t2, legend_only2,
                          ncol = 2,
                          align = "hv")

ggsave("figures/109a_fig_3_shoval_meta_ourdata.jpeg", shoval_plots2, width = 8, height = 12)

# Our Pareto front analysis -----------------------------------------------

# Light -------------------------------------------------------------------

df.filt <- df.l %>% 
  mutate(
    z.y = comp,
    z.x = r.max
  ) # Specify the x and y variables and their 95% CIs

plot(z.y ~ z.x, data= df.filt) # do we need to exclude any points?

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

I.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
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

mean(null.df$a.emp.n >= a.emp) # p-value 0.261

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.1 %>%
    
    mutate(
      z.x.sim = sample(par.res.1$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(par.res.1$z.y, replace = FALSE)      # Separately reassign y
      
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

mean(null.df2$a.emp.n >= a.emp) # p-value 0!!!

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 0.07434, p 0.51071
summary(q75, se = "boot", R = 1000) # -0.03627, p 0.79385
summary(q90, se = "boot", R = 1000) # 0.24586, p 0.42630 

I.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
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

I.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
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

summary(q50, se = "boot", R = 1000) # -1.13080, p 0.00009
summary(q75, se = "boot", R = 1000) # -1.11815, p 0.00153
summary(q90, se = "boot", R = 1000) # -1.41332, p 0.00027

I.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
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

# Nitrogen -------------------------------------------------------------------

df.filt <- df.n %>% 
  mutate(
    z.y = comp,
    z.x = r.max
  ) # Specify the x and y variables and their 95% CIs

plot(z.y ~ z.x, data= df.filt) # do we need to exclude any points?

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

# fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF
fit <- lm(z.y ~ z.x, data = par.res.1) # Fit a scam to the PF

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

N.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/N*)", 
       title = "B — Nitrogen") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2015" = "chocolate1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1")
  ) +
  
  ylim(0,52) +
  xlim(0, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

N.scam  # Display the plot

p.n3 <- p.n + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE)  # Adding scam PF fits
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

mean(null.df$a.emp.n >= a.emp) # p-value 0.789

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.1 %>%
    
    mutate(
      z.x.sim = sample(par.res.1$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(par.res.1$z.y, replace = FALSE)      # Separately reassign y
      
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

mean(null.df2$a.emp.n >= a.emp) # p-value 0.165

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 5.25457, p 0.18391
summary(q75, se = "boot", R = 1000) # 12.10420, p 0.00495
summary(q90, se = "boot", R = 1000) # 26.00688, p 0.00542 

N.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/N*)", 
       title = "B — Nitrogen") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2015" = "chocolate1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1")
  ) +
  
  ylim(0,52) +
  xlim(0, 3) +
  
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

N.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/N*)", 
       title = "B — Nitrogen") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2015" = "chocolate1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1")
  ) +
  
  ylim(0,52) +
  xlim(0, 3) +
  
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

mean(null.df3$a.emp.n >= a.emp) # p-value 0.236

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -8.19660, p 0.44780
summary(q75, se = "boot", R = 1000) # -10.59224, p 0.43187
summary(q90, se = "boot", R = 1000) # -13.60864, p 0.28216

N.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/N*)", 
       title = "B — Nitrogen") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2015" = "chocolate1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1")
  ) +
  
  ylim(0,52) +
  xlim(0, 3) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

N.qr2 # Display the plot

# Phosphorous -------------------------------------------------------------------

df.filt <- df.p %>% 
  mutate(
    z.y = comp,
    z.x = r.max
  ) # Specify the x and y variables and their 95% CIs

plot(z.y ~ z.x, data= df.filt) # do we need to exclude any points?

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

P.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       title = "C — Phosphorous") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2015" = "chocolate1",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
  ylim(-0.1,340) +
  xlim(0,2.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

P.scam  # Display the plot

p.p3 <- p.p + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE)  # Adding scam PF fits
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

mean(null.df$a.emp.n >= a.emp) # p-value 0.004

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.1 %>%
    
    mutate(
      z.x.sim = sample(par.res.1$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(par.res.1$z.y, replace = FALSE)      # Separately reassign y
      
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

mean(null.df2$a.emp.n >= a.emp) # p-value 0

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # -18.82682, p 0.18324
summary(q75, se = "boot", R = 1000) # -77.86779, p 0.00031
summary(q90, se = "boot", R = 1000) # -114.21657, p 0.00002 

P.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       title = "C — Phosphorous") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2015" = "chocolate1",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
  ylim(-0.1,340) +
  xlim(0,2.6) +
  
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

P.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       title = "C — Phosphorous") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2015" = "chocolate1",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
  ylim(-0.1,340) +
  xlim(0,2.6) +
  
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

mean(null.df3$a.emp.n >= a.emp) # p-value 0

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -178.60295, p 0.00017
summary(q75, se = "boot", R = 1000) # -173.05897, p 0.00000
summary(q90, se = "boot", R = 1000) # -149.80227, p 0.00000

P.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/P*)", 
       title = "C — Phosphorous") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2015" = "chocolate1",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
  ylim(-0.1,340) +
  xlim(0,2.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

P.qr2 # Display the plot

# Temperature -------------------------------------------------------------------

df.filt <- df.t %>% 
  mutate(
    z.y = T.br,
    z.x = r.max
  ) # Specify the x and y variables and their 95% CIs

plot(z.y ~ z.x, data= df.filt) # do we need to exclude any points?

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

T.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       title = "D — Temperature") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Thomas et al., 2012" = "magenta2")
  ) +
  
  ylim(-1, 40) +
  xlim(-0.05,3.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.scam  # Display the plot

p.t3 <- p.t + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE)  # Adding scam PF fits
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

mean(null.df$a.emp.n >= a.emp) # p-value 0.133

# Randomize across optimal and near-optimal points

null.df2 <- data.frame(       # Null model results
  a.emp.n = numeric(),       # Area above the Pareto front (polygon with xmax, ymax) 
  n.PF = numeric(),          # Extract the number of Pareto front points
  stringsAsFactors = FALSE            
)

for (i in 1:1000){
  
  shuffled.df2 <- par.res.1 %>%
    
    mutate(
      z.x.sim = sample(par.res.1$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(par.res.1$z.y, replace = FALSE)      # Separately reassign y
      
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

mean(null.df2$a.emp.n >= a.emp) # p-value 0

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 1.74209, p 0.00004
summary(q75, se = "boot", R = 1000) # 1.07406, p 0.26978
summary(q90, se = "boot", R = 1000) # -0.24125, p 0.90881 

T.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       title = "D — Temperature") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Thomas et al., 2012" = "magenta2")
  ) +
  
  ylim(-1, 40) +
  xlim(-0.05,3.6) +
  
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

T.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       title = "D — Temperature") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Thomas et al., 2012" = "magenta2")
  ) +
  
  ylim(-1, 40) +
  xlim(-0.05,3.6) +
  
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

mean(null.df3$a.emp.n >= a.emp) # p-value 0.3005

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -4.01353, p 0.00034
summary(q75, se = "boot", R = 1000) # -5.21847, p 0.00000
summary(q90, se = "boot", R = 1000) # -4.66452, p 0.00840

T.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, color = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Thermal breadth (°C)", 
       title = "D — Temperature") +  # labels
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Bestion et al., 2018" = "firebrick2", 
               "Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Thomas et al., 2012" = "magenta2")
  ) +
  
  ylim(-1, 40) +
  xlim(-0.05,3.6) +
  
  theme_classic() +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 12, face = "bold"),  
    axis.text = element_text(size = 10, face ="plain"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.03)# theme stuff
  )

T.qr2 # Display the plot

# Assemble and save the figures -------------------------------------------

meta_PF <- plot_grid(I.scam, N.scam, P.scam, T.scam, legend_only,
                          ncol = 2,
                          align = "hv")

ggsave("figures/109c_fig_Sx_meta_parfronts.jpeg", meta_PF, width = 8, height = 12)

meta_PF2 <- plot_grid(I.scam2, N.scam2, P.scam2, T.scam2, legend_only,
                     ncol = 2,
                     align = "hv")

ggsave("figures/109d_fig_Sx_meta_parfronts_33pruned.jpeg", meta_PF2, width = 8, height = 12)

meta_QR <- plot_grid(I.qr, N.qr, P.qr, T.qr, legend_only,
                     ncol = 2,
                     align = "hv")

ggsave("figures/109d_fig_Sx_meta_quantregs.jpeg", meta_QR, width = 8, height = 12)

meta_QR2 <- plot_grid(I.qr2, N.qr2, P.qr2, T.qr2, legend_only,
                     ncol = 2,
                     align = "hv")

ggsave("figures/109e_fig_Sx_meta_quantregs_33pruned.jpeg", meta_QR2, width = 8, height = 12)

# Create panel with both Pareto fronts and triangles ----------------------

meta_both <- plot_grid(p.l3, p.n3, p.p3, p.t3, legend_only,
                     ncol = 2,
                     align = "hv")

ggsave("figures/109f_fig_Sx_meta_both_PF_shoval.jpeg", meta_both, width = 8, height = 12)
