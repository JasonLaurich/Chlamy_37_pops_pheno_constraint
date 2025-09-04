# Jason R Laurich

# September 3rd, 2025

# Loading all of the disparate files that contain relevant ecological data and combining into a single massive data frame.
# Am also going to (1) re-estimate thermal breadth at relevant levels (0.56 for our study, 0.1 for comparison with other datasets) and 
# (2) calculate the 95% HDPI intervals for all estimates using the posterior distributions from our JAGS models!

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(Deriv)

lactin2 <- function(temp, cf.a, cf.b, cf.delta_t, cf.tmax) {
  exp(cf.a * temp) - exp(cf.a * cf.tmax - ((cf.tmax - temp) / cf.delta_t)) + cf.b
}

find_roots <- function(cf.a, cf.b, cf.delta_t, cf.tmax) {
  fn <- function(temp) lactin2(temp, cf.a, cf.b, cf.delta_t, cf.tmax) - target
  
  root1 <- tryCatch(uniroot(fn, lower = 5, upper = cf.tmax - 10)$root, error = function(e) NA)
  root2 <- tryCatch(uniroot(fn, lower = cf.tmax - 10, upper = 45)$root, error = function(e) NA)
  
  c(root1, root2)
}

lactin2_deriv <- function(temp, cf.a, cf.b, cf.tmax, cf.delta_t) {
  rho <- cf.a
  T_max <- cf.tmax
  delta_T <- cf.delta_t
  
  term1 <- rho * exp(rho * temp)
  term2 <- (1 / delta_T) * exp(rho * T_max - (T_max - temp) / delta_T)
  
  return(term1 - term2)
} # Derivative of the Lactin II function

lactin2_0.56 <- function(temp, cf.a, cf.b, cf.tmax, cf.delta_t) {
  exp(cf.a * temp) - exp(cf.a * cf.tmax - (cf.tmax - temp) / cf.delta_t) + cf.b - 0.56
} # OK we're going to modify the function to calculate T_breadth, based on a modified lactin.


# Load & examine data -----------------------------------------------------

df.hist <- read.csv("data-processed/16_rfus_time_summary.csv") # This file has the history (ancestry, evolutionary treatment) for each well ID.
head(df.hist)

df.hist$anc <- as.factor(df.hist$ancestor_id) # Will port these into combined dataset, sorting by population. 
df.hist$evol <- as.factor(df.hist$treatment)
df.hist$Pop.fac <- as.factor(df.hist$population)

df.tpc <- read.csv("data-processed/05_TPCs.csv") # TPC summary data. We'll use Lactin II model fits. 
head(df.tpc)

df <- df.tpc[,c(2,5:9)] # We'll work with the analytically determined parameters (Deriv package)
head(df)

df.hist.agg <- df.hist %>% # Aggregate the history data for merging.
  group_by(Pop.fac) %>%
  summarise(
    anc = first(anc),
    evol = first(evol),
    .groups = "drop" 
  )

df <- merge(df.hist.agg, df, by = "Pop.fac", all.x = TRUE)
head(df) 
str(df)

df$Pop.fac <- as.factor(df$Pop.fac) # Looks great! Now we layer in other summary metrics.
names(df) <- c("Pop.fac", "anc", "evol", "T.min", "T.max", "T.br", "Topt", "r.max_T") # rename a few variables. 

# Bring in light data

df.I <- read.csv("data-processed/06b_Monod_light_estimates.csv") # I* data 
head(df.I)
str(df.I)

df.I.par <- df.I[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.I.par) <- c("Pop.fac", "K.s_I", "r.max_I", "I.jag", "I.mth") # rename
df.I.par$I.comp <- 1/df.I.par$I.mth

df <- merge(df, df.I.par, by = "Pop.fac", all.x = TRUE)

# Bring in nitrogen data

df.N <- read.csv("data-processed/07b_Monod_nitrogen_estimates.csv") # N* data 
head(df.N)
str(df.N)

df.N.par <- df.N[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.N.par) <- c("Pop.fac", "K.s_N", "r.max_N", "N.jag", "N.mth") # rename
df.N.par$N.comp <- 1/df.N.par$N.mth

df <- merge(df, df.N.par, by = "Pop.fac", all.x = TRUE)

# Bring in phosphorous data

df.P <- read.csv("data-processed/08b_Monod_phosphorous_estimates.csv") # P* data 
head(df.P)
str(df.P)

df.P.par <- df.P[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.P.par) <- c("Pop.fac", "K.s_P", "r.max_P", "P.jag", "P.mth") # rename
df.P.par$P.comp <- 1/df.P.par$P.mth

df <- merge(df, df.P.par, by = "Pop.fac", all.x = TRUE)

# Bring in salt data

df.S <- read.csv("data-Processed/09b_salt_tolerance_estimates.csv") # S* data 
head(df.S)
str(df.S)

df.S$Pop.fac <- gsub("Anc ", "anc", df.S$Pop.fac) # ancestors are written differently here

df.S.par <- df.S[, c("Pop.fac", "r.max", "c.mod", "c.pred")] # need these
names(df.S.par) <- c("Pop.fac", "r.max_S", "S.c.mod", "S.c.pred") # rename

df <- merge(df, df.S.par, by = "Pop.fac", all.x = TRUE)

# Bring in stoichiometry data

df.stoich <- read.csv("data-processed/17_stoich_data.csv") # Stoichiometry data for N and P
head(df.stoich)
str(df.stoich)

levels(as.factor(df.stoich$Name)) # OK these are not the same numbers as are present in our other datasets
length(unique(df.stoich$Name)) # But the total number is the same

df.stoich %>% # Let's look at the mean values here (there are 2 data points for each population)
  group_by(Name) %>% 
  summarize(mean.N.µg.l = mean(N.µg.l), mean.P.µg.l = mean(P.µg.l)) %>% 
  print(n=37)

df.stoich.sum <- df.stoich %>%    # Let's just save the means
  group_by(Name) %>% 
  summarize(mean.N.µg.l = mean(N.µg.l), mean.P.µg.l = mean(P.µg.l)) %>% 
  print(n=37)

# Let's bring in an identity mapping file. 

df.id <- read.csv("data-processed/18_id_mapping.csv") # File containing the identity assignment

df.id

df.stoich.sum <- df.stoich.sum %>% # Join this with the df.id matching schema
  left_join(df.id, by = c("Name" = "sample"))

head(df.stoich.sum) # Correct population assignment in the population column

# Finally the pigment data

df.pig <- read.csv("data-processed/19_pigment_data.csv") # Pigment data

head(df.pig)
str(df.pig)

levels(as.factor(df.pig$HPLC.Nummer)) # These are already numbered 1-37, ie. the missing populations do not feature here. 

# Calculate & store new metrics ---------------------------------------------------

# Recalculate Tbr at 0.56

df.tpc.fits <- read.csv('data-processed/05a_TPC_fits.csv') # Load the data with the TPC shape parameters
head(df.tpc.fits)

df.tpc.fits <- df.tpc.fits %>% # Pivot to wide format and remove the extra population
  filter(Model == 'Lactin 2', Pop.fac != 'cc1629') %>% 
  select(Pop.fac, Parameter, mean) %>%
  pivot_wider(names_from = Parameter, values_from = mean)

head(df.tpc.fits)

target <- 0.56

df.tpc.roots <- df.tpc.fits %>%
  mutate(roots = pmap(list(cf.a, cf.b, cf.delta_t, cf.tmax), find_roots)) %>%
  transmute(
    Pop.fac,
    T.min.0.56 = map_dbl(roots, 1),
    T.max.0.56 = map_dbl(roots, 2)
  )

head(df.tpc.roots)

# OK now we need to load all of the TPC objects and calculate the 95% HPD around my variables of interest (T.br_0.56, µ_max, T.br_0)

hpd.temp.df <- data.frame(             # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
  Pop.fac = character(),               # Population
  Tmin.l.hpdi = numeric(),             # Tmin, lower 95% HDPI 
  Tmin.u.hpdi = numeric(),             # Tmin, upper 95% HDPI 
  Tmax.l.hpdi = numeric(),             # Tmax, lower 95% HDPI 
  Tmax.u.hpdi = numeric(),             # Tmax, upper 95% HDPI
  µmax.l.hpdi = numeric(),             # µmax, lower 95% HDPI
  µmax.u.hpdi = numeric(),             # µmax, upper 95% HDPI 
  Tmin.0.56.l.hpdi = numeric(),        # Tmin at 0.56, lower 95% HDPI 
  Tmin.0.56.u.hpdi = numeric(),        # Tmin at 0.56, upper 95% HDPI 
  Tmax.0.56.l.hpdi = numeric(),        # Tmax at 0.56, lower 95% HDPI 
  Tmax.0.56.u.hpdi = numeric()         # Tmax at 0.56, upper 95% HDPI
)

for (i in c(1:36, 38)){                                                    # Temperature R2jags (model fits), excluding 37 which is population cc1629
  load(paste0("R2jags-objects/pop_", i, "_lactin.RData"))                  # load the models
  posts <- as.data.frame(
    lac_jag$BUGSoutput$sims.list[c("cf.a","cf.b","cf.delta_t","cf.tmax")]  # Extract the posterior estimates for variables of interest (6000 total)
  )
  
  hpd.df <- data.frame(        # A dataframe to store the 6000 estimates for each parameter in! Reset for each population
    Tmin = numeric(),          # Tmin
    Tmax = numeric(),          # Tmax
    µmax = numeric(),          # µmax
    Tmin.0.56 = numeric(),     # Tmin at 0.56
    Tmax.0.56 = numeric()      # Tmax at 0.56
  )
  
  for (p in 1:6000){                                                     # For each posterior estimate
    
    Topt <- uniroot(                                                     # Find Topt
      function(temp) lactin2_deriv(temp, cf.a = posts$cf.a[p],
                                   cf.b = posts$cf.b[p], 
                                   cf.tmax = posts$cf.tmax[p], 
                                   cf.delta_t = posts$cf.delta_t[p]),
      interval = c(10, 45)
    )$root                              
    
    Tmin <- uniroot(lactin2, interval = c(0, T_opt),                      # Find Tmin
                    cf.a = posts$cf.a[p], cf.b = posts$cf.b[p], 
                    cf.tmax = posts$cf.tmax[p], 
                    cf.delta_t = posts$cf.delta_t[p])$root
    
    Tmax <- uniroot(lactin2, interval = c(T_opt,45),                      # Find Tmax
                    cf.a = posts$cf.a[p], cf.b = posts$cf.b[p], 
                    cf.tmax = posts$cf.tmax[p], 
                    cf.delta_t = posts$cf.delta_t[p])$root
    
    µmax <- lactin2(temp=Topt, cf.a = posts$cf.a[p],                      # Find µmax
                    cf.b = posts$cf.b[p], cf.tmax = posts$cf.tmax[p], 
                    cf.delta_t = posts$cf.delta_t[p])
    
    Tmin.0.56 <- uniroot(lactin2_0.56, interval = c(Tmin, Topt),          # Find Tmin at 0.56
                         cf.a = posts$cf.a[p], cf.b = posts$cf.b[p], 
                         cf.tmax = posts$cf.tmax[p], 
                         cf.delta_t = posts$cf.delta_t[p])$root
    
    Tmax.0.56 <- uniroot(lactin2_0.56, interval = c(Topt, Tmax),          # Find Tmax at 0.56
                         cf.a = posts$cf.a[p], cf.b = posts$cf.b[p], 
                         cf.tmax = posts$cf.tmax[p], 
                         cf.delta_t = posts$cf.delta_t[p])$root
    
    hpd.df <- rbind(hpd.df, data.frame(        # A dataframe to store the 6000 estimates for each parameter in! Reset for each population
      Tmin = Tmin,                             # Tmin
      Tmax = Tmax,                             # Tmax
      µmax = µmax,                             # µmax
      Tmin.0.56 = Tmin.0.56,                   # Tmin at 0.56
      Tmax.0.56 = Tmax.0.56                    # Tmax at 0.56
    ))
    
  }
  
  hpd.temp.df <- rbind(hpd.temp.df, data.frame(                    # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
    Pop.fac = df.tpc$Pop.fac (which matches i in the Pop.num column),               # Population
    Tmin.l.hpdi = the 2.5% lower quantile from hpd.df$Tmin,             # Tmin, lower 95% HDPI 
    Tmin.u.hpdi = numeric(),             # Tmin, upper 95% HDPI 
    Tmax.l.hpdi = numeric(),             # Tmax, lower 95% HDPI 
    Tmax.u.hpdi = numeric(),             # Tmax, upper 95% HDPI
    µmax.l.hpdi = numeric(),             # µmax, lower 95% HDPI
    µmax.u.hpdi = numeric(),             # µmax, upper 95% HDPI 
    Tmin.0.56.l.hpdi = numeric(),        # Tmin at 0.56, lower 95% HDPI 
    Tmin.0.56.u.hpdi = numeric(),        # Tmin at 0.56, upper 95% HDPI 
    Tmax.0.56.l.hpdi = numeric(),        # Tmax at 0.56, lower 95% HDPI 
    Tmax.0.56.u.hpdi = numeric()         # Tmax at 0.56, upper 95% HDPI
  )
  
}

# Organize & compile into single file -------------------------------------


write.csv(df, "data-processed/14_summary_metric_table.csv") # Save summary table

