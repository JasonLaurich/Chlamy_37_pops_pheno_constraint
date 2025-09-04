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

df.I <- read.csv("data-processed/06b_Monod_light_estimates.csv") # Light data 
head(df.I)
str(df.I)

df.N <- read.csv("data-processed/07b_Monod_nitrogen_estimates.csv") # Nitrogen data 
head(df.N)
str(df.N)

df <- merge(df, df.N.par, by = "Pop.fac", all.x = TRUE)

df.P <- read.csv("data-processed/08b_Monod_phosphorous_estimates.csv") # Phosphorous data 
head(df.P)
str(df.P)

df.S <- read.csv("data-Processed/09b_salt_tolerance_estimates.csv") # Salt data 
head(df.S)
str(df.S)

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

df.stoich.sum <- as.data.frame(df.stoich %>%    # Let's just save the means
  group_by(Name) %>% 
  summarize(mean.N.µg.l = mean(N.µg.l), mean.P.µg.l = mean(P.µg.l)) %>% 
  print(n=37))

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

df.pig<- df.pig %>%
  mutate(population = df.id$population)

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

###### Temperature HDPIs ######
# OK now we need to load all of the TPC objects and calculate the 95% HPD around my variables of interest (T.br_0.56, µ_max, T.br_0)

hpd.temp.df <- data.frame(             # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
  Pop.fac = character(),               # Population
  Tbr.l.hpdi = numeric(),              # Tbr, lower 95% HDPI 
  Tbr.u.hpdi = numeric(),              # Tbr, upper 95% HDPI              
  µmax.l.hpdi = numeric(),             # µmax, lower 95% HDPI
  µmax.u.hpdi = numeric(),             # µmax, upper 95% HDPI 
  Tbr.0.56.l.hpdi = numeric(),         # Tbr at 0.56, lower 95% HDPI 
  Tbr.0.56.u.hpdi = numeric()          # Tbr at 0.56, upper 95% HDPI
)

for (i in c(1:36, 38)){                                                    # Temperature R2jags (model fits), excluding 37 which is population cc1629
  load(paste0("R2jags-objects/pop_", i, "_lactin.RData"))                  # load the models
  posts <- as.data.frame(
    lac_jag$BUGSoutput$sims.list[c("cf.a","cf.b","cf.delta_t","cf.tmax")]  # Extract the posterior estimates for variables of interest (6000 total)
  )
  
  hpd.df <- data.frame(        # A dataframe to store the 6000 estimates for each parameter in! Reset for each population
    Tbr = numeric(),           # Tbr
    µmax = numeric(),          # µmax
    Tbr.0.56 = numeric()       # Tbr at 0.56
  )
  
  for (p in 1:6000){                                                     # For each posterior estimate
    
    Topt <- uniroot(                                                     # Find Topt
      function(temp) lactin2_deriv(temp, cf.a = posts$cf.a[p],
                                   cf.b = posts$cf.b[p], 
                                   cf.tmax = posts$cf.tmax[p], 
                                   cf.delta_t = posts$cf.delta_t[p]),
      interval = c(10, 45)
    )$root                              
    
    Tmin <- uniroot(lactin2, interval = c(-200, Topt),                      # Find Tmin
                    cf.a = posts$cf.a[p], cf.b = posts$cf.b[p], 
                    cf.tmax = posts$cf.tmax[p], 
                    cf.delta_t = posts$cf.delta_t[p])$root
    
    Tmax <- uniroot(lactin2, interval = c(Topt,65),                      # Find Tmax
                    cf.a = posts$cf.a[p], cf.b = posts$cf.b[p], 
                    cf.tmax = posts$cf.tmax[p], 
                    cf.delta_t = posts$cf.delta_t[p])$root
    
    Tbr <- Tmax - Tmin                                                    # Calculate Tbr
    
    µmax <- lactin2(temp=Topt, cf.a = posts$cf.a[p],                      # Find µmax
                    cf.b = posts$cf.b[p], cf.tmax = posts$cf.tmax[p], 
                    cf.delta_t = posts$cf.delta_t[p])
    
    Tmin.0.56 <- uniroot(lactin2_0.56, interval = c(-50, Topt),          # Find Tmin at 0.56
                         cf.a = posts$cf.a[p], cf.b = posts$cf.b[p], 
                         cf.tmax = posts$cf.tmax[p], 
                         cf.delta_t = posts$cf.delta_t[p])$root
    
    Tmax.0.56 <- uniroot(lactin2_0.56, interval = c(Topt, 65),          # Find Tmax at 0.56
                         cf.a = posts$cf.a[p], cf.b = posts$cf.b[p], 
                         cf.tmax = posts$cf.tmax[p], 
                         cf.delta_t = posts$cf.delta_t[p])$root
    
    Tbr.0.56 <- Tmax.0.56 - Tmin.0.56                                     # Calculate Tbr at 0.56
    
    hpd.df <- rbind(hpd.df, data.frame(        # A dataframe to store the 6000 estimates for each parameter in! Reset for each population
      Tbr = Tbr,                               # Tbr
      µmax = µmax,                             # µmax
      Tbr.0.56 = Tbr.0.56                      # Tbr at 0.56
    ))
    
    print(p)
    
  }
  
  hpd.temp.df <- rbind(hpd.temp.df, data.frame(                          # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
    Pop.fac = df.tpc %>% filter(Pop.num == i) %>% 
      pull(Pop.fac) %>% dplyr::first(),                                  # Population as a factor, pulled from the df.tpc summary table
    Tbr.l.hpdi = quantile(hpd.df$Tbr, probs = 0.025),                    # Tbr, lower 95% HDPI 
    Tbr.u.hpdi = quantile(hpd.df$Tbr, probs = 0.975),                    # Tbr, upper 95% HDPI 
    µmax.l.hpdi = quantile(hpd.df$µmax, probs = 0.025),                  # µmax, lower 95% HDPI
    µmax.u.hpdi = quantile(hpd.df$µmax, probs = 0.975),                  # µmax, upper 95% HDPI 
    Tbr.0.56.l.hpdi = quantile(hpd.df$Tbr.0.56, probs = 0.025),          # Tbr at 0.56, lower 95% HDPI 
    Tbr.0.56.u.hpdi = quantile(hpd.df$Tbr.0.56, probs = 0.975)           # Tbr at 0.56, upper 95% HDPI 
  ))
  
  print(i)
  
}

hpd.temp.df

###### Light Monod HDPIs ######
# Load the light Monods and estimate 95% HPDIs for I.ks, I.comp and µmax
# I want to calculate I.comp at each posterior after adjusting for conversion to photons

hpd.light.df <- data.frame(              # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
  Pop.fac = character(),                 # Population
  I_Ks.l.hpdi = numeric(),               # IKs, lower 95% HDPI 
  I_Ks.u.hpdi = numeric(),               # IKs, upper 95% HDPI
  Icomp.0.1.l.hpdi = numeric(),          # Icomp.0.1, lower 95% HDPI 
  Icomp.0.1.u.hpdi = numeric(),          # Icomp.0.1, upper 95% HDPI  
  Icomp.0.56.l.hpdi = numeric(),         # Icomp.0.56, lower 95% HDPI 
  Icomp.0.56.u.hpdi = numeric(),         # Icomp.0.56, upper 95% HDPI              
  I_µmax.l.hpdi = numeric(),             # µmax, lower 95% HDPI
  I_µmax.u.hpdi = numeric()              # µmax, upper 95% HDPI 
)

for (i in c(1:37)){                                                    # Light R2jags (model fits)
  load(paste0("R2jags-objects/pop_", i, "_light_monod.RData"))         # load the models
  
  I_Ks.l.hpdi <- monod_jag$BUGSoutput$summary[1,3]*2.5                 # Convert %s to photons
  I_Ks.u.hpdi <- monod_jag$BUGSoutput$summary[1,7]*2.5
  I_µmax.l.hpdi <- monod_jag$BUGSoutput$summary[3,3]
  I_µmax.u.hpdi <- monod_jag$BUGSoutput$summary[3,7]
  
  posts <- as.data.frame(
    monod_jag$BUGSoutput$sims.list[c("K_s","r_max")]  # Extract the posterior estimates for variables of interest (6000 total)
  )
  
  hpd.df <- data.frame(        # A dataframe to store the 6000 estimates for each parameter in! Reset for each population
    I.comp.0.1 = numeric(),    # I.comp at 0.1
    I.comp.0.56 = numeric()    # I.comp at 0.56
  )
  
  for (p in 1:6000){                                                           # For each posterior estimate
  
    I.comp.0.1 <- 1 / (0.1 * posts$K_s[p] * 2.5 / (posts$r_max[p] - 0.1))      # Converting Ks to photons from % (for consistency with other studies)
    
    I.comp.0.56 <- 1 / (0.56 * posts$K_s[p] * 2.5 / (posts$r_max[p] - 0.56))   # Converting Ks to photons from % (for consistency with other studies)
    
    hpd.df <- rbind(hpd.df, data.frame(        # A dataframe to store the 6000 estimates for each parameter in! Reset for each population
      I.comp.0.1 = I.comp.0.1,                 # I.comp at 0.1             
      I.comp.0.56 = I.comp.0.56                # I.comp at 0.56
    ))
    
    print(p)
    
  }

  hpd.light.df <- rbind(hpd.light.df, data.frame(                           # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
    Pop.fac = df.tpc %>% filter(Pop.num == i) %>% 
      pull(Pop.fac) %>% dplyr::first(),                                     # Population as a factor, pulled from the df.tpc summary table
    I_Ks.l.hpdi = I_Ks.l.hpdi,                                              # IKs, lower 95% HDPI 
    I_Ks.u.hpdi = I_Ks.u.hpdi,                                              # IKs, upper 95% HDPI
    Icomp.0.1.l.hpdi = quantile(hpd.df$I.comp.0.1, probs = 0.025),          # Icomp.0.1, lower 95% HDPI 
    Icomp.0.1.u.hpdi = quantile(hpd.df$I.comp.0.1, probs = 0.975),          # Icomp.0.1, upper 95% HDPI  
    Icomp.0.56.l.hpdi = quantile(hpd.df$I.comp.0.56, probs = 0.025),        # Icomp.0.56, lower 95% HDPI 
    Icomp.0.56.u.hpdi = quantile(hpd.df$I.comp.0.56, probs = 0.975),        # Icomp.0.56, upper 95% HDPI 
    I_µmax.l.hpdi = I_µmax.l.hpdi,                                          # µmax, lower 95% HDPI
    I_µmax.u.hpdi = I_µmax.u.hpdi                                           # µmax, upper 95% HDPI 
  ))
  
  print(i)
  
}

hpd.light.df

hpd.light.df$Pop.fac[37] <- "cc1690"  # Fix the name of the last entry (the raw data for light, nitrogen, phosphorous, and salt did not include cc1629)

###### Nitrogen Monod HDPIs ######
# Load the nitrogen Monods and extract 95% HPDIs for N.Ks, N.comp and µmax
# Calculate N.comp at 0.1 and 0.56

hpd.nit.df <- data.frame(                # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
  Pop.fac = character(),                 # Population
  N_Ks.l.hpdi = numeric(),               # NKs, lower 95% HDPI 
  N_Ks.u.hpdi = numeric(),               # NKs, upper 95% HDPI
  Ncomp.0.1.l.hpdi = numeric(),          # Ncomp.0.1, lower 95% HDPI 
  Ncomp.0.1.u.hpdi = numeric(),          # Ncomp.0.1, upper 95% HDPI  
  Ncomp.0.56.l.hpdi = numeric(),         # Ncomp.0.56, lower 95% HDPI 
  Ncomp.0.56.u.hpdi = numeric(),         # Ncomp.0.56, upper 95% HDPI               
  N_µmax.l.hpdi = numeric(),             # µmax, lower 95% HDPI
  N_µmax.u.hpdi = numeric()              # µmax, upper 95% HDPI 
)

for (i in c(1:37)){                                                    # Nitrogen R2jags (model fits)
  load(paste0("R2jags-objects/pop_", i, "_nitrogen_monod.RData"))      # load the models
  
  N_Ks.l.hpdi <- monod_jag$BUGSoutput$summary[1,3]
  N_Ks.u.hpdi <- monod_jag$BUGSoutput$summary[1,7]
  N_µmax.l.hpdi <- monod_jag$BUGSoutput$summary[3,3]
  N_µmax.u.hpdi <- monod_jag$BUGSoutput$summary[3,7]
  
  posts <- as.data.frame(
    monod_jag$BUGSoutput$sims.list[c("K_s","r_max")]  # Extract the posterior estimates for variables of interest (6000 total)
  )
  
  hpd.df <- data.frame(        # A dataframe to store the 6000 estimates for each parameter in! Reset for each population
    N.comp.0.1 = numeric(),    # N.comp at 0.1
    N.comp.0.56 = numeric()    # N.comp at 0.56
  )
  
  for (p in 1:6000){                                                           # For each posterior estimate
    
    N.comp.0.1 <- 1 / (0.1 * posts$K_s[p] / (posts$r_max[p] - 0.1))      
    
    N.comp.0.56 <- 1 / (0.56 * posts$K_s[p] / (posts$r_max[p] - 0.56))   
    
    hpd.df <- rbind(hpd.df, data.frame(        # A dataframe to store the 6000 estimates for each parameter in! Reset for each population
      N.comp.0.1 = N.comp.0.1,                 # N.comp at 0.1             
      N.comp.0.56 = N.comp.0.56                # N.comp at 0.56
    ))
    
    print(p)
    
  }
  
  hpd.nit.df <- rbind(hpd.nit.df, data.frame(                               # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
    Pop.fac = df.tpc %>% filter(Pop.num == i) %>% 
      pull(Pop.fac) %>% dplyr::first(),                                     # Population as a factor, pulled from the df.tpc summary table
    N_Ks.l.hpdi = N_Ks.l.hpdi,                                              # NKs, lower 95% HDPI 
    N_Ks.u.hpdi = N_Ks.u.hpdi,                                              # NKs, upper 95% HDPI
    Ncomp.0.1.l.hpdi = quantile(hpd.df$N.comp.0.1, probs = 0.025),          # Ncomp.0.1, lower 95% HDPI 
    Ncomp.0.1.u.hpdi = quantile(hpd.df$N.comp.0.1, probs = 0.975),          # Ncomp.0.1, upper 95% HDPI  
    Ncomp.0.56.l.hpdi = quantile(hpd.df$N.comp.0.56, probs = 0.025),        # Ncomp.0.56, lower 95% HDPI 
    Ncomp.0.56.u.hpdi = quantile(hpd.df$N.comp.0.56, probs = 0.975),        # Ncomp.0.56, upper 95% HDPI 
    N_µmax.l.hpdi = N_µmax.l.hpdi,                                          # µmax, lower 95% HDPI
    N_µmax.u.hpdi = N_µmax.u.hpdi                                           # µmax, upper 95% HDPI 
  ))
  
  print(i)
  
}

hpd.nit.df

hpd.nit.df$Pop.fac[37] <- "cc1690"

###### Phosphorous Monod HDPIs ######
# Load the phosphorous Monods and extract 95% HPDIs for P.Ks, P.comp and µmax
# Calculate P.comp at 0.1 and 0.56

hpd.phos.df <- data.frame(               # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
  Pop.fac = character(),                 # Population
  P_Ks.l.hpdi = numeric(),               # PKs, lower 95% HDPI 
  P_Ks.u.hpdi = numeric(),               # PKs, upper 95% HDPI
  Pcomp.0.1.l.hpdi = numeric(),          # Pcomp.0.1, lower 95% HDPI 
  Pcomp.0.1.u.hpdi = numeric(),          # Pcomp.0.1, upper 95% HDPI  
  Pcomp.0.56.l.hpdi = numeric(),         # Pcomp.0.56, lower 95% HDPI 
  Pcomp.0.56.u.hpdi = numeric(),         # Pcomp.0.56, upper 95% HDPI               
  P_µmax.l.hpdi = numeric(),             # µmax, lower 95% HDPI
  P_µmax.u.hpdi = numeric()              # µmax, upper 95% HDPI 
)

for (i in c(1:37)){                                                    # Phosphorous R2jags (model fits)
  load(paste0("R2jags-objects/pop_", i, "_phosphorous_monod.RData"))   # load the models
  
  P_Ks.l.hpdi <- monod_jag$BUGSoutput$summary[1,3]
  P_Ks.u.hpdi <- monod_jag$BUGSoutput$summary[1,7]
  P_µmax.l.hpdi <- monod_jag$BUGSoutput$summary[3,3]
  P_µmax.u.hpdi <- monod_jag$BUGSoutput$summary[3,7]
  
  posts <- as.data.frame(
    monod_jag$BUGSoutput$sims.list[c("K_s","r_max")]  # Extract the posterior estimates for variables of interest (6000 total)
  )
  
  hpd.df <- data.frame(        # A dataframe to store the 6000 estimates for each parameter in! Reset for each population
    P.comp.0.1 = numeric(),    # P.comp at 0.1
    P.comp.0.56 = numeric()    # P.comp at 0.56
  )
  
  for (p in 1:6000){                                                           # For each posterior estimate
    
    P.comp.0.1 <- 1 / (0.1 * posts$K_s[p] / (posts$r_max[p] - 0.1))      
    
    P.comp.0.56 <- 1 / (0.56 * posts$K_s[p] / (posts$r_max[p] - 0.56))   
    
    hpd.df <- rbind(hpd.df, data.frame(        # A dataframe to store the 6000 estimates for each parameter in! Reset for each population
      P.comp.0.1 = P.comp.0.1,                 # P.comp at 0.1             
      P.comp.0.56 = P.comp.0.56                # P.comp at 0.56
    ))
    
    print(p)
    
  }
  
  hpd.phos.df <- rbind(hpd.phos.df, data.frame(                             # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
    Pop.fac = df.tpc %>% filter(Pop.num == i) %>% 
      pull(Pop.fac) %>% dplyr::first(),                                     # Population as a factor, pulled from the df.tpc summary table
    P_Ks.l.hpdi = P_Ks.l.hpdi,                                              # NKs, lower 95% HDPI 
    P_Ks.u.hpdi = P_Ks.u.hpdi,                                              # NKs, upper 95% HDPI
    Pcomp.0.1.l.hpdi = quantile(hpd.df$P.comp.0.1, probs = 0.025),          # Pcomp.0.1, lower 95% HDPI 
    Pcomp.0.1.u.hpdi = quantile(hpd.df$P.comp.0.1, probs = 0.975),          # Pcomp.0.1, upper 95% HDPI  
    Pcomp.0.56.l.hpdi = quantile(hpd.df$P.comp.0.56, probs = 0.025),        # Pcomp.0.56, lower 95% HDPI 
    Pcomp.0.56.u.hpdi = quantile(hpd.df$P.comp.0.56, probs = 0.975),        # Pcomp.0.56, upper 95% HDPI 
    P_µmax.l.hpdi = P_µmax.l.hpdi,                                          # µmax, lower 95% HDPI
    P_µmax.u.hpdi = P_µmax.u.hpdi                                           # µmax, upper 95% HDPI 
  ))
  
  print(i)
  
}

hpd.phos.df

hpd.phos.df$Pop.fac[37] <- "cc1690"

###### Salt tolerance HDPIs ######
# Load the salt tolerance curves and extract 95% HPDIs for c and µmax.
# Again, these are terms directly estimated by the models, so I can just extract and tabulate them.

hpd.salt.df <- data.frame(                # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
  Pop.fac = character(),                  # Population
  S.c.l.hpdi = numeric(),                 # c, lower 95% HDPI 
  S.c.u.hpdi = numeric(),                 # c, upper 95% HDPI
  S_µmax.l.hpdi = numeric(),              # µmax, lower 95% HDPI
  S_µmax.u.hpdi = numeric()               # µmax, upper 95% HDPI 
)

for (i in c(1:37)){                                                    # Salt R2jags (model fits)
  load(paste0("R2jags-objects/pop_", i, "_salt_tolerance.RData"))      # load the models
  
  S.c.l.hpdi <- monod_jag$BUGSoutput$summary[3,3] 
  S.c.u.hpdi <- monod_jag$BUGSoutput$summary[3,7]
  S_µmax.l.hpdi <- monod_jag$BUGSoutput$summary[1,3]
  S_µmax.u.hpdi <- monod_jag$BUGSoutput$summary[1,7]
  
  hpd.salt.df <- rbind(hpd.salt.df, data.frame(         # A dataframe to store the summary data (highest posterior density intervals, HDPIs) for each population
    Pop.fac = df.tpc %>% filter(Pop.num == i) %>% 
      pull(Pop.fac) %>% dplyr::first(),               # Population as a factor, pulled from the df.tpc summary table
    S.c.l.hpdi = S.c.l.hpdi,                          # c, lower 95% HDPI 
    S.c.u.hpdi = S.c.u.hpdi,                          # c, upper 95% HDPI
    S_µmax.l.hpdi = S_µmax.l.hpdi,                    # µmax, lower 95% HDPI
    S_µmax.u.hpdi = S_µmax.u.hpdi                     # µmax, upper 95% HDPI 
  ))
  
  print(i)
  
}

hpd.salt.df

hpd.salt.df$Pop.fac[37] <- "cc1690"

# Organize & compile into single file -------------------------------------

df.summ <- df.hist.agg %>%                                                      # Make a beautiful, comprehensive final data frame
  filter(Pop.fac != 'cc1629') %>% 
  transmute(Pop.fac = Pop.fac,
            Anc = anc,
            Evol = evol) %>% 
  
  left_join(                                                               
    df.tpc %>%                                                                  # Add in µmax and thermal breadth at 0
      filter(Model == 'Lactin 2', Pop.fac != 'cc1629') %>%
      mutate(T.br.0 = T.max - T.min) %>% 
      transmute(T.µ.max = r.max,
                T.br.0 = T.br.0,
                Pop.fac = Pop.fac), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(                                  
    df.tpc.roots %>%                                                            # Add in thermal breadth at 0.56
      mutate(T.br.0.56 = T.max.0.56 - T.min.0.56) %>% 
      select(T.br.0.56, Pop.fac), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(
    hpd.temp.df %>%                                                             # Add in 95% HDPIs for thermal traits
      transmute(Pop.fac = Pop.fac,
                T.µ.max.L = µmax.l.hpdi,
                T.µ.max.U = µmax.u.hpdi,
                T.br.0.L = Tbr.l.hpdi,
                T.br.0.U = Tbr.u.hpdi,
                T.br.0.56.L = Tbr.0.56.l.hpdi,
                T.br.0.56.U = Tbr.0.56.u.hpdi), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  
    left_join(
      df.I %>%                                                                  # Add in light competition traits
        mutate(I.star.0.56 = R.mth*2.5,                                         # Need to convert light %s to µmol photons (consistent with other interspecific datasets)
               I.comp.0.56 = 1/I.star.0.56,
               I.Ks = K.s*2.5,
               I.comp.0.1 = 1/(0.1*I.Ks/(r.max - 0.1))) %>%                     # Calculate I.comp at 0.10
        transmute(Pop.fac = Pop.fac,
                  I.µ.max = r.max,
                  I.comp.0.1 = I.comp.0.1,
                  I.comp.0.56 = I.comp.0.56), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(
    hpd.light.df %>%                                                            # Add in 95% HDPIs for light traits
      transmute(Pop.fac = Pop.fac,
                I.µ.max.L = I_µmax.l.hpdi,
                I.µ.max.U = I_µmax.u.hpdi,
                I.comp.0.1.L = Icomp.0.1.l.hpdi,
                I.comp.0.1.U = Icomp.0.1.u.hpdi,
                I.comp.0.56.L = Icomp.0.56.l.hpdi,
                I.comp.0.56.U = Icomp.0.56.u.hpdi), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(
    df.N %>%                                                                    # Add in nitrogen competition traits
              mutate(N.comp.0.56 = 1/R.mth,
                     N.comp.0.1 = 1/(0.1*K.s/(r.max - 0.1))) %>%               
              transmute(Pop.fac = Pop.fac,
                        N.µ.max = r.max,
                        N.comp.0.1 = N.comp.0.1,
                        N.comp.0.56 = N.comp.0.56), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(
    hpd.nit.df %>%                                                              # Add in 95% HDPIs for nitrogen traits
      transmute(Pop.fac = Pop.fac,
                N.µ.max.L = N_µmax.l.hpdi,
                N.µ.max.U = N_µmax.u.hpdi,
                N.comp.0.1.L = Ncomp.0.1.l.hpdi,
                N.comp.0.1.U = Ncomp.0.1.u.hpdi,
                N.comp.0.56.L = Ncomp.0.56.l.hpdi,
                N.comp.0.56.U = Ncomp.0.56.u.hpdi), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(
    df.P %>%                                                                    # Add in phosphorous competition traits
      mutate(P.comp.0.56 = 1/R.mth,
             P.comp.0.1 = 1/(0.1*K.s/(r.max - 0.1))) %>%               
      transmute(Pop.fac = Pop.fac,
                P.µ.max = r.max,
                P.comp.0.1 = P.comp.0.1,
                P.comp.0.56 = P.comp.0.56), by = c("Pop.fac" = "Pop.fac")) %>%
    
  left_join(
    hpd.phos.df %>%                                                             # Add in 95% HDPIs for phosphorous traits
      transmute(Pop.fac = Pop.fac,
                P.µ.max.L = P_µmax.l.hpdi,
                P.µ.max.U = P_µmax.u.hpdi,
                P.comp.0.1.L = Pcomp.0.1.l.hpdi,
                P.comp.0.1.U = Pcomp.0.1.u.hpdi,
                P.comp.0.56.L = Pcomp.0.56.l.hpdi,
                P.comp.0.56.U = Pcomp.0.56.u.hpdi), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(
    df.S %>%                                                                    # Add in salt competition traits
      transmute(Pop.fac = Pop.fac,
                S.µ.max = r.max,
                S.c = c.mod), by = c("Pop.fac" = "Pop.fac")) %>%  
  
  left_join(
    hpd.salt.df %>%                                                             # Add in 95% HDPIs for salt traits
      transmute(Pop.fac = Pop.fac,
                S.µ.max.L = S_µmax.l.hpdi,
                S.µ.max.U = S_µmax.u.hpdi,
                S.c.L = S.c.l.hpdi,
                S.c.U = S.c.u.hpdi), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(
    df.stoich.sum %>%                                                           # Add in stoichiometry data
      transmute(Pop.fac = population,
                mean.N.µg.l = mean.N.µg.l,
                mean.P.µg.l = mean.P.µg.l), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(
    df.pig %>%                                                           # Finally, add in pigment data!
      transmute(Pop.fac = population,
                chl.a = chl.a,
                chl.b = chl.b,
                luthein = luthein), by = c("Pop.fac" = "Pop.fac")) 
        
View(df.summ)

write.csv(df.summ, "data-processed/20_summary_table.csv") # Save summary table!
