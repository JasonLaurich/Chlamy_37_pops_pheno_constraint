# Jason R Laurich

# February 11th, 2026

# Here I am going to import the various datasets that contain synthesis data for phytoplankton growth across
  # light, nitrogen, phosphorous, and temperature gradients.

# I will then apply analyses based on Shoval 2012 (Science) to calculate the triangularity of our data and compare this 
  # estimate to null randomizations. I will further integrate the Pareto front analyses I have been employing elsewhere
  # in this data analysis pipeline (scripts 08/09)

# Inputs: 47_Edwards_2016_light_monods.csv, 50_Lewington_2019_light_monods.csv, 53_Levasseur_2025_light_monods.csv, 66_Narwani_2015_summary.csv,
  # 68_Edwards_2015_summary.csv, 56_Lewington_2019_nit_monods.csv, 59_Levasseur_2025_nit_monods.csv, 61_Bestion_2018_phos_monods.csv,
  # 64_Levasseur_2025_phos_monods.csv, 32_Bestion_2018_TPCs.csv, 39_Edwards_2016_TPCs.csv, 36_Lewington_2019_TPCs.csv, 45_Levasseur_2025_TPCs.csv,
  # 29_Thomas_2012_TPCs.csv, 67a_light_metadata_sp.csv, 69_nit_metadata_sp.csv, 70_phos_metadata_sp.csv, 71_temp_metadata_sp.csv
# Outputs: in processed-data : 67_light_metadata.csv, 69_nit_metadata.csv, 70_phos_metadata.csv, 71_temp_metadata.csv, 72_synthesis_metadata_sp_summary.csv

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

# Load and compile the data sets -------------------------------------------------------------------

###### Light synthesis data ######

# Edwards 2016

df.l.e <- read.csv("processed-data/47_Edwards_2016_light_monods.csv")
head(df.l.e)

df.l.e$dataset <- "Edwards et al., 2016"

# Lewington-Pearce 2019

df.l.l <- read.csv("processed-data/50_Lewington_2019_light_monods.csv")
head(df.l.l)

df.l.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.l.lv <- read.csv("processed-data/53_Levasseur_2025_light_monods.csv")
head(df.l.lv)

df.l.lv$dataset <- "Levasseur et al., 2025"

# Narwani 2015

df.l.n <- read.csv("processed-data/66_Narwani_2015_summary.csv")
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

plot(df.l$comp ~ df.l$r.max) # Some really wonky (high) estimates of 1/R* that suggest error (filter at > 3.5)

df.l <- df.l %>% 
  filter(comp < 3.5,
         r.max < 2.5,
         !is.na(comp))

write.csv(df.l, "processed-data/67_light_metadata.csv") # Summary

###### Nitrogen synthesis data ######

# Edwards 2015

df.n.e <- read.csv("processed-data/68_Edwards_2015_summary.csv")
head(df.n.e)

df.n.e <- df.n.e %>% 
  filter(!is.na(k_nit_m),
         !is.na(mu_nit)) %>% 
  mutate(R = 0.1*k_nit_m/(mu_nit - 0.1))

df.n.e$dataset <- "Edwards et al., 2015"

# Lewington-Pearce 2019

df.n.l <- read.csv("processed-data/56_Lewington_2019_nit_monods.csv")
head(df.n.l)

df.n.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.n.lv <- read.csv("processed-data/59_Levasseur_2025_nit_monods.csv")
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

plot(df.n$comp ~ df.n$r.max) # Some really wonky (high) estimates of 1/R* that suggest error or nitrogen fixation?

df.n <- df.n %>% 
  filter(comp < 60,
         !is.na(r.max),
         !is.na(comp))

write.csv(df.n, "processed-data/69_nit_metadata.csv") # Summary

###### Phosphorous synthesis data ######

# Bestion 2018

df.p.b <- read.csv("processed-data/61_Bestion_2018_phos_monods.csv")
head(df.p.b)

df.p.b$dataset <- "Bestion et al., 2018"

# Edwards 2015

df.p.e <- read.csv("processed-data/68_Edwards_2015_summary.csv")
head(df.p.e)

df.p.e <- df.p.e %>% 
  filter(!is.na(k_p_m),
         !is.na(mu_p)) %>% 
  mutate(R = 0.1*k_p_m/(mu_p - 0.1))

df.p.e$dataset <- "Edwards et al., 2015"

# Levasseur 2025

df.p.lv <- read.csv("processed-data/64_Levasseur_2025_phos_monods.csv")
head(df.p.lv)

df.p.lv$dataset <- "Levasseur et al., 2025"

# Narwani 2015

df.p.n <- read.csv("processed-data/66_Narwani_2015_summary.csv")
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

plot(df.p$comp ~ df.p$r.max) # Some really wonky (high) estimates of 1/R* that suggest error or nitrogen fixation?

df.p <- df.p %>% 
  filter(comp < 300,
         r.max < 2.5,
         !is.na(r.max),
         !is.na(comp))

write.csv(df.p, "processed-data/70_phos_metadata.csv") # Summary

###### Temperature synthesis data ######

# Bestion 2018

df.t.b <- read.csv("processed-data/32_Bestion_2018_TPCs.csv")
head(df.t.b)

df.t.b$dataset <- "Bestion et al., 2018"

# Edwards 2016

df.t.e <- read.csv("processed-data/39_Edwards_2016_TPCs.csv")
head(df.t.e)

df.t.e$dataset <- "Edwards et al., 2016"

# Lewington-Pearce 2019

df.t.l <- read.csv("processed-data/36_Lewington_2019_TPCs.csv")
head(df.t.l)

df.t.l$dataset <- "Lewington-Pearce et al., 2019"

# Levasseur 2025

df.t.lv <- read.csv("processed-data/45_Levasseur_2025_TPCs.csv")
head(df.t.lv)

df.t.lv$dataset <- "Levasseur et al., 2025"

# Thomas 2012

df.t.t <- read.csv("processed-data/29_Thomas_2012_TPCs.csv")
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

plot(df.t$T.br ~ df.t$r.max)

df.t <- df.t %>% 
  filter(r.max < 3.5,
         T.br < 36,
         !is.na(r.max),
         !is.na(T.br))

write.csv(df.t, "processed-data/71_temp_metadata.csv") # Summary

###### Get sp-level summary data (means) ######

# For comparing across gradients, where we cannot be confident that observations on the x and y axis were obtained from the same study.

df.l <- read.csv("processed-data/67a_light_metadata_sp.csv")
head(df.l)

length(unique(df.l$Sp.name)) # 99 species
sort(unique(df.l$Sp.name))

# Nitrogen

df.n <- read.csv("processed-data/69a_nit_metadata_sp.csv")
head(df.n)

length(unique(df.n$Sp.name)) #30 species
sort(unique(df.n$Sp.name))

# Phosphorous 

df.p <- read.csv("processed-data/70a_phos_metadata_sp.csv")
head(df.p)

length(unique(df.p$Sp.name)) # 57 species
sort(unique(df.p$Sp.name))

# Temperature

df.t <- read.csv("processed-data/71a_temp_metadata_sp.csv")
head(df.t)

length(unique(df.t$Sp.name)) # 117 species
sort(unique(df.t$Sp.name))

# OK so we have r.max and either T.br or comp for all. We need to combine these into a single data frame. 
# Combine the data frames

df.t.1 <- df.t %>%
  group_by(Sp.name) %>%
  filter(Sp.name != "NA") %>% 
  summarise(
    r.max.T = mean(r.max, na.rm = TRUE),
    T.br    = mean(T.br,  na.rm = TRUE),
    .groups = "drop"
  )

df.n.1 <- df.n %>%
  group_by(Sp.name) %>%
  filter(Sp.name != "NA") %>% 
  summarise(
    r.max.N = mean(r.max, na.rm = TRUE),
    comp.N  = mean(comp,  na.rm = TRUE),
    .groups = "drop"
  )

df.p.1 <- df.p %>%
  group_by(Sp.name) %>%
  filter(Sp.name != "NA") %>% 
  summarise(
    r.max.P = mean(r.max, na.rm = TRUE),
    comp.P  = mean(comp,  na.rm = TRUE),
    .groups = "drop"
  )

df.l.1 <- df.l %>%
  group_by(Sp.name) %>%
  filter(Sp.name != "NA") %>% 
  summarise(
    r.max.L = mean(r.max, na.rm = TRUE),
    comp.L  = mean(comp,  na.rm = TRUE),
    .groups = "drop"
  )

df <- df.t.1 %>%
  full_join(df.n.1, by = "Sp.name") %>%
  full_join(df.p.1, by = "Sp.name") %>%
  full_join(df.l.1, by = "Sp.name")

write.csv(df, "processed-data/72_synthesis_metadata_sp_summary.csv") # Summary

# OK so now we will move on to performing the relevant analyses.

# Light v Light -------------------------------------------------------------------

df.filt <- df.l %>% 
  mutate(
    z.y = comp,
    z.x = r.max
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

p.l <- ggplot(df.filt, aes(z.x, z.y, colour = dataset)) +
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
  
  annotate(
    "text",
    x = -0.3 + 0.85 * (3 - -0.3),
    y = -0.1 + 0.95 * (3.25 - -0.1),
    label = "Tri\nPF\nQR",
    hjust = 0,
    vjust = 1,
    size = 3.5,
    fontface = "bold"
  ) +
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
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

mean(null.df$poly.area <= tri$poly.area[1]) # 0.003

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.010

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

I.scam <- ggplot(df.filt, aes(x = z.x, y = z.y, colour = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  ylim(-0.10, 3.25) +
  xlim(-0.25, 3) +
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
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

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt)  

summary(q50, se = "boot", R = 1000) # 0.07434, p 0.51071
summary(q75, se = "boot", R = 1000) # -0.03627, p 0.79385
summary(q90, se = "boot", R = 1000) # 0.24586, p 0.42630 

I.qr <- ggplot(df.filt, aes(x = z.x, y = z.y, colour = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  ylim(-0.10, 3.25) +
  xlim(-0.25, 3) +
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
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

I.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, colour = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
  geom_line(data = pred.curve.2, aes(x = z.x, y = z.y), color = "black", size = 1.1, linetype = "dashed", inherit.aes = FALSE) +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  ylim(-0.10, 3.25) +
  xlim(-0.25, 3) +
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
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

I.qr2 <- ggplot(df.filt3, aes(x = z.x, y = z.y, colour = dataset)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_abline(intercept = coef(q90)[1], slope = coef(q90)[2], lwd = 1.1) +
  geom_abline(intercept = coef(q75)[1], slope = coef(q75)[2], lwd = 1.1, linetype = "dashed") +
  geom_abline(intercept = coef(q50)[1], slope = coef(q50)[2], lwd = 1.1, linetype = "dotted") +
  
  labs(x = "Maximum exponential growth rate (µ max)",    
       y = "Competitive ability (1/I*)", 
       title = "A — Light") +  # labels
  
  ylim(-0.10, 3.25) +
  xlim(-0.25, 3) +
  
  scale_color_manual(
    name = "Data set",  # Update the legend title
    values = c("Edwards et al., 2016" = "goldenrod1",
               "Lewington-Pearce et al., 2019" = "forestgreen",
               "Levasseur et al., 2025" = "royalblue1",
               "Narwani et al., 2015" = "darkorchid3")
  ) +
  
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

df.filt.x <- df.filt %>% 
  filter(z.x < 45)

###### Shoval et al triangular analysis ######

hull.id <- chull(df.filt.x$z.x, df.filt.x$z.y)

df.hull.full <- df.filt.x[hull.id,] %>% 
  bind_rows(df.filt.x[hull.id[1],]) # Close the polygon

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

df.hull <- df.filt.x[hull.id,]  # Don't close the polygon

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

best.tri$area/best.tri$poly.area # 1.012402

###### Randomization ######

null.df <- data.frame(                           # Null model results
  poly.area = numeric(),                         # Polygon area
  best.tri.area = numeric(),                     # Area of the smallest triangle
  stringsAsFactors = FALSE            
)

for (z in 1:1000){
  
  df.null <- df.filt.x %>%          # scramble the data
    
    mutate(
      z.x = sample(df.filt.x$z.x, replace = FALSE),             # Randomly assign rmax
      z.y = sample(df.filt.x$z.y, replace = FALSE)              # Separately reassign z.y
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

mean(null.df$poly.area <= tri$poly.area[1]) # 0.136

mean(null.df$best.tri.area/null.df$poly.area <= min(tri$area)/tri$poly.area[1]) # 0.027

###### Pareto fronts ######

df.filt.x <- df.filt.x %>% # Scale for euclidean distance calculation and point exclusion determination
  mutate(
    z.y2 = scale(z.y)[, 1],
    z.x2 = scale(z.x)[, 1]
  )  

x.ref <- min(df.filt.x$z.x2, na.rm = TRUE) # Min x
y.ref <- min(df.filt.x$z.y2, na.rm = TRUE) # Min y

df.filt.x <- df.filt.x %>% # Calculate Euclidean distance from min
  mutate(
    distance = sqrt((z.x2 - x.ref)^2 + (z.y2 - y.ref)^2),
    dist.sc = scale(distance)
  ) %>%
  arrange(distance) #  Distance for point exclusion and 75th quantile calculation. 

df.filt.x%>%
  { 
    bind_rows(
      arrange(., desc(z.y)) %>% slice_head(n = 3),
      arrange(., z.y) %>% slice_head(n = 3),
      arrange(., desc(z.x)) %>% slice_head(n = 3),
      arrange(., z.x) %>% slice_head(n = 3)
    )
  } %>%
  print() # Display the outliers

par.res.1 <- par_frt(df.filt.x, xvar = "z.x", yvar = "z.y") # Get the raw Pareto Front. Considering remaining extreme points as possible escapees

fit <- scam(z.y ~ s(z.x, bs = "mpd", k = min(6,nrow(par.res.1))), data = par.res.1) # Fit a scam to the PF

x.vals <- seq(min(df.filt.x$z.x), max(df.filt.x$z.x), length.out = 100) # Generate an x sequence for plotting

pred.curve.1 <- data.frame( # Get the corresponding y values
  z.x = x.vals,
  z.y = predict(fit, newdata = data.frame(z.x = x.vals))
)

df.filt2 <- df.filt.x %>% # Filter out based on Euclidean distance from min
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
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
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

p.ln3 <- p.ln  + geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE)  # Adding scam PF fits
p.ln3

###### Polygonal empty space analysis (Li et al 2019) ######

x.max <- max(df.filt.x$z.x) # Extract the max values for x and y.
y.max <- max(df.filt.x$z.y)

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
  
  shuffled.df <- df.filt.x %>%
    
    mutate(
      z.x.sim = sample(df.filt.x$z.x, replace = FALSE),     # Randomly assign x
      
      z.y.sim = sample(df.filt.x$z.y, replace = FALSE),     # Separately reassign y
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

mean(null.df$a.emp.n >= a.emp) # p-value 0.351

###### Quantile regression ######

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt.x) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt.x) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt.x)  

summary(q50, se = "boot", R = 1000) # 0.00783, p 0.48316 
summary(q75, se = "boot", R = 1000) # -0.00073, p 0.97499
summary(q90, se = "boot", R = 1000) # -0.03021, p 0.46629

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

df.filt3 <- df.filt.x %>%
  arrange(distance) %>%                           # smallest → largest
  slice((floor(0.6667 * n()) + 1):n()) %>%          # keep *second* half
  select(-distance)

LN.scam2 <- ggplot(df.filt3, aes(x = z.x, y = z.y)) +  # We'll lay out the PFs onto our raw data
  geom_point(size = 2) +  # Scatter plot of raw data
  
  geom_line(data = pred.curve.1, aes(x = z.x, y = z.y), color = "black", size = 1.1, inherit.aes = FALSE) +  # Adding scam PF fits
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

mean(null.df3$a.emp.n >= a.emp) # p-value 0.002

q50  <- rq(z.y ~ z.x, tau = 0.50, data = df.filt3) 
q75  <- rq(z.y ~ z.x, tau = 0.75, data = df.filt3) 
q90  <- rq(z.y ~ z.x, tau = 0.90, data = df.filt3)  

summary(q50, se = "boot", R = 1000) # -0.02915, p 0.29413
summary(q75, se = "boot", R = 1000) # -0.03021, p 0.19789
summary(q90, se = "boot", R = 1000) # -0.03298, p 0.19928

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
