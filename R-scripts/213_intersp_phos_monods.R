# Jason R Laurich

# January 12th, 2025

# I am going to fit all of phosphorous Monod curves for all of the relevant datasets: 
# That's the Bestion 2018, Narwani 2015, Edwards 2015, and Levasseur 2025 datasets
# Edwards 2015 and Narwani 2015 only report estimates, not raw data. 


# Where possible (ie. phos and T co-vary in the data — Bestion, and Levasseur),
# I will fit Monod curves at the temperature closest to a species' Topt (as determined by TPC fitting)
# If I was unable to fit a TPC to the data (or if it was removed from consideration) I will fit Monod functions
# at the temperature where the highest growth occurs in the raw data. I will record both numbers just for fun!

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(dplyr)
library(Deriv)
library(car)
library(minpack.lm)
library(cowplot)

library(R2jags)
library(mcmcplots)
library(bayestestR)

# And R2jags settings
inits.monod <- function() { # Set the initial values for our Monod curve
  list(
    r_max = runif(1, 0.1, 5), # Initial guess for r_max
    K_s = runif(1, 0.1, 5),   # Initial guess for K_s
    sigma = runif(1, 0.1, 1)  # Initial guess for error
  )
}

parameters.monod <- c("r_max", "K_s", "sigma", "r_pred_new") # Save these

# Let's do larger models for the final things (10 times larger)
ni.fit <- 330000    # iterations / chain
nb.fit <- 30000     # burn in periods for each chain
nt.fit <- 300       # thinning interval : (330,000 - 30,000) / 300 = 1000 posterior estimates / chain
nc.fit <- 3         # number of chains, total of 3,000 estimates for each model. 

# Bestion 2018 ------------------------------------------------------------

###### Load the data ######

df.b.raw <- read.csv("data-processed/11_Bestion_2018_raw_data.csv") # Raw growth data
head(df.b.raw)

df.b <-df.b.raw %>% 
  rename(temp = Temperature_c,
         id.number = SpeciesNb,
         phos = Phosphate_c,
         Species.name = SpeciesName)

length(unique(df.b$id.number))

df.b.t <- read.csv("data-processed/501c_Bestion_2018_TPCs.csv") # Raw data file
head(df.b.t)

###### Model fitting ######

bestion.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  
  K.s.mod = numeric(),      # Half saturation constant (model output)
  K.s.post = numeric(),     # Half saturation constant (posterior median)
  K.s.min = numeric(),      # Half saturation constant (lower HDPI)
  K.s.max = numeric(),      # Half saturation constant (upper HDPI)
  K.s.na = numeric(),       # % NA returns
  
  r.max.mod = numeric(),    # Maximum growth rate (model output)
  r.max.post = numeric(),   # Maximum growth rate (posterior median)
  r.max.min = numeric(),    # Maximum growth rate (lower HDPI)
  r.max.max = numeric(),    # Maximum growth rate (upper HDPI)
  r.max.na = numeric(),     # % NA returns
  
  R.post = numeric(),       # R* (m = 0.1) (posterior median)
  R.min = numeric(),        # R* (m = 0.1) (lower HDPI)
  R.max = numeric(),        # R* (m = 0.1) (upper HDPI)
  R.na = numeric(),         # % NA returns
  
  T.opt.TPC = numeric(),    # Topt (if there is a fitted TPC) — temp closest to this
  T.max.raw = numeric(),    # ~ Topt (temp in the raw data with highest growth rate)
  
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~3000)
  stringsAsFactors = FALSE            
)

n <-0 # progression tracker

for (i in unique(df.b$id.number[df.b$id.number >= 1])) { # for each species ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.b %>% 
    filter(id.number == i) %>% 
    filter(!is.na(mu))
  
  t.max <- df.i %>%                             # Find the temperature with the highest growth data
    filter(mu == max(mu, na.rm = TRUE)) %>% 
    pull(temp)
  
  if (i %in% df.b.t$Sp.id) {                    # pull T.opt from external dataset
    t.opt <- df.b.t %>%
      filter(Sp.id == i) %>%
      pull(T.opt) 
    fit <- "y"
  } else {                                      # fallback: Topt from raw data
    t.opt <- t.max
    fit <- "n"
  }
  
  t.sel <- df.i$temp[which.min(abs(df.i$temp - t.opt))] # Pick the temperature closest to Topt
  
  df.i <- df.i %>% 
    filter(temp == t.sel)
  
  trait <- df.i$mu    # format the data for jags
  N.obs <- length(trait)
  
  phos <- df.i$phos
  
  S.pred <- seq(0, 50, 0.05) # Phosphorous gradient we're interested in - upped the granularity here
  N.S.pred <-length(S.pred) # We'll reset this internally since the gradient varies substantially
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = phos, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod.jag <- jags( # Run the phos Monod function. 
    data = jag.data,
    inits = inits.monod,
    parameters.to.save = parameters.monod,
    model.file = "monod.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  post <- as.data.frame(monod.jag$BUGSoutput$sims.matrix) # The posteriors
  post$R <- 0.1*post$K_s/(post$r_max - 0.1)
  
  bestion.summ.df <- rbind(bestion.summ.df, data.frame(                     # Add data
    Sp.id = df.i$id.number[1],                                                  # Species id
    Sp.name = df.i$Species.name[1],                                             # Species name
    
    K.s.mod = monod.jag$BUGSoutput$summary[1,1],                                # Half saturation constant (model output)
    K.s.post = median(post$K_s, na.rm = T),                                     # Half saturation constant (posterior median)
    K.s.min = hdi(post$K_s, ci = 0.95)$CI_low,                                  # Half saturation constant (lower HDPI)
    K.s.max = hdi(post$K_s, ci = 0.95)$CI_high,                                 # Half saturation constant (upper HDPI)
    K.s.na = mean(is.na(post$K_s)),                                             # % NA returns
    
    r.max.mod = monod.jag$BUGSoutput$summary[3,1],                              # Maximum growth rate (model output)
    r.max.post = median(post$r_max, na.rm = T),                                 # Maximum growth rate (posterior median)
    r.max.min = hdi(post$r_max, ci = 0.95)$CI_low,                              # Maximum growth rate (lower HDPI)
    r.max.max = hdi(post$r_max, ci = 0.95)$CI_high,                             # Maximum growth rate (upper HDPI)
    r.max.na = mean(is.na(post$r_max)),                                         # % NA returns
    
    R.post = median(post$R, na.rm = T),                                         # R* (m = 0.1) (posterior median)
    R.min = hdi(post$R, ci = 0.95)$CI_low,                                      # R* (m = 0.1) (lower HDPI)
    R.max = hdi(post$R, ci = 0.95)$CI_high,                                     # R* (m = 0.1) (upper HDPI)
    R.na = mean(is.na(post$R)),                                                 # % NA returns
    
    T.opt.TPC = ifelse(fit == "y", t.opt, NA_real_),                            # Topt (if there is a fitted TPC) — temp closest to this
    T.max.raw = t.max                                                           # ~ Topt (temp in the raw data with highest growth rate)
  ))
  
  phos_sum <- monod.jag$BUGSoutput$summary[c(1:3, (50/0.05) + 5),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = df.i$id.number[1],                 # Species #
      Sp.name = df.i$Species.name[1],            # Species name      
      Parameter = rownames(phos_sum)[j],         # Model parameter (e.g. K_s, r_max, etc.)
      mean = phos_sum[j,1],                      # Posterior mean
      Rhat = phos_sum[j,8],                      # Rhat values
      n.eff = phos_sum[j,9]                      # Sample size estimates (should be ~3000)
    ))
    
  }
  
  print(paste("Done", n, "of ", length(unique(df.b$id.number))))
  
}

write.csv(bestion.summ.df, "data-processed/503a_Bestion_2018_phos_monods.csv") # Bestion_2018 summary table
write.csv(fit.df, "data-processed/503b_Bestion_2018_phos_monods_fits.csv") # Save model fit summary table

# Levasseur 2025 ------------------------------------------------------------

###### Load the data ######

df.lv.raw <- read.csv("data-processed/406_Levasseur_2025_µ_estimates_phosphorous.csv") # Raw data file

head(df.lv.raw)

df.lv <-df.lv.raw %>% 
  rename(mu = r.exp,
         Species.name = Sp.id) %>% 
  mutate(id.number = as.integer(factor(df.lv.raw$Sp.id)))

length(unique(df.lv$id.number))

df.lv.t <- read.csv("data-processed/501i_Levasseur_2025_TPCs.csv") # Raw data file

###### Model fitting ######

levasseur.summ.df <- data.frame(   # We'll create a dataframe to store the data as we fit models.
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  
  K.s.mod = numeric(),      # Half saturation constant (model output)
  K.s.post = numeric(),     # Half saturation constant (posterior median)
  K.s.min = numeric(),      # Half saturation constant (lower HDPI)
  K.s.max = numeric(),      # Half saturation constant (upper HDPI)
  K.s.na = numeric(),       # % NA returns
  
  r.max.mod = numeric(),    # Maximum growth rate (model output)
  r.max.post = numeric(),   # Maximum growth rate (posterior median)
  r.max.min = numeric(),    # Maximum growth rate (lower HDPI)
  r.max.max = numeric(),    # Maximum growth rate (upper HDPI)
  r.max.na = numeric(),     # % NA returns
  
  R.post = numeric(),       # R* (m = 0.1) (posterior median)
  R.min = numeric(),        # R* (m = 0.1) (lower HDPI)
  R.max = numeric(),        # R* (m = 0.1) (upper HDPI)
  R.na = numeric(),         # % NA returns
  
  T.opt.TPC = numeric(),    # Topt (if there is a fitted TPC) — temp closest to this
  T.max.raw = numeric(),    # ~ Topt (temp in the raw data with highest growth rate)
  
  stringsAsFactors = FALSE  # Avoid factor conversion
)

fit.df <- data.frame(       # Save model fit estimates for examination
  Sp.id = numeric(),        # Species id
  Sp.name = character(),    # Species name
  Parameter = character(),  # Model parameter (e.g. cf.a, cf.tmax, etc.)
  mean = numeric(),         # Posterior mean
  Rhat = numeric(),         # Rhat values
  n.eff = numeric(),        # Sample size estimates (should be ~3000)
  stringsAsFactors = FALSE            
)

n <-0 # progression tracker

for (i in unique(df.lv$id.number[df.lv$id.number >= 1])) { # for each species ID. Can adjust starting point if running in chunks
  
  n <- n + 1
  
  df.i <- df.lv %>% 
    filter(id.number == i) %>% 
    filter(!is.na(mu))
  
  t.max <- df.i %>%                             # Find the temperature with the highest growth data
    filter(mu == max(mu, na.rm = TRUE)) %>% 
    pull(temp)
  
  if (i %in% df.lv.t$Sp.id) {                    # pull T.opt from external dataset
    t.opt <- df.lv.t %>%
      filter(Sp.id == i) %>%
      pull(T.opt) 
    fit <- "y"
  } else {                                      # fallback: Topt from raw data
    t.opt <- t.max
    fit <- "n"
  }
  
  t.sel <- df.i$temp[which.min(abs(df.i$temp - t.opt))] # Pick the temperature closest to Topt
  
  df.i <- df.i %>% 
    filter(temp == t.sel)
  
  trait <- df.i$mu    # format the data for jags
  N.obs <- length(trait)
  
  phos <- df.i$phos
  
  S.pred <- seq(0, 50, 0.05) # Phosphorous gradient we're interested in - upped the granularity here
  N.S.pred <-length(S.pred) # We'll reset this internally since the gradient varies substantially
  
  jag.data <- list(trait = trait, N.obs = N.obs, S = phos, S.pred = S.pred, N.S.pred = N.S.pred)
  
  monod.jag <- jags( # Run the phos Monod function. 
    data = jag.data,
    inits = inits.monod,
    parameters.to.save = parameters.monod,
    model.file = "monod.txt",
    n.thin = nt.fit,
    n.chains = nc.fit,
    n.burnin = nb.fit,
    n.iter = ni.fit,
    DIC = TRUE,
    working.directory = getwd()
  )
  
  post <- as.data.frame(monod.jag$BUGSoutput$sims.matrix) # The posteriors
  post$R <- 0.1*post$K_s/(post$r_max - 0.1)
  
  levasseur.summ.df <- rbind(levasseur.summ.df, data.frame(                     # Add data
    Sp.id = df.i$id.number[1],                                                  # Species id
    Sp.name = df.i$Species.name[1],                                             # Species name
    
    K.s.mod = monod.jag$BUGSoutput$summary[1,1],                                # Half saturation constant (model output)
    K.s.post = median(post$K_s, na.rm = T),                                     # Half saturation constant (posterior median)
    K.s.min = hdi(post$K_s, ci = 0.95)$CI_low,                                  # Half saturation constant (lower HDPI)
    K.s.max = hdi(post$K_s, ci = 0.95)$CI_high,                                 # Half saturation constant (upper HDPI)
    K.s.na = mean(is.na(post$K_s)),                                             # % NA returns
    
    r.max.mod = monod.jag$BUGSoutput$summary[3,1],                              # Maximum growth rate (model output)
    r.max.post = median(post$r_max, na.rm = T),                                 # Maximum growth rate (posterior median)
    r.max.min = hdi(post$r_max, ci = 0.95)$CI_low,                              # Maximum growth rate (lower HDPI)
    r.max.max = hdi(post$r_max, ci = 0.95)$CI_high,                             # Maximum growth rate (upper HDPI)
    r.max.na = mean(is.na(post$r_max)),                                         # % NA returns
    
    R.post = median(post$R, na.rm = T),                                         # R* (m = 0.1) (posterior median)
    R.min = hdi(post$R, ci = 0.95)$CI_low,                                      # R* (m = 0.1) (lower HDPI)
    R.max = hdi(post$R, ci = 0.95)$CI_high,                                     # R* (m = 0.1) (upper HDPI)
    R.na = mean(is.na(post$R)),                                                 # % NA returns
    
    T.opt.TPC = ifelse(fit == "y", t.opt, NA_real_),                            # Topt (if there is a fitted TPC) — temp closest to this
    T.max.raw = t.max                                                           # ~ Topt (temp in the raw data with highest growth rate)
  ))
  
  phos_sum <- monod.jag$BUGSoutput$summary[c(1:3, (50/0.05) + 5),] # Have to create a new frame for summaries (not listed 1 to 6)
  
  for (j in 1:4){
    fit.df <- rbind(fit.df, data.frame(          # Model performance data
      Sp.id = df.i$id.number[1],                 # Species #
      Sp.name = df.i$Species.name[1],            # Species name      
      Parameter = rownames(phos_sum)[j],         # Model parameter (e.g. K_s, r_max, etc.)
      mean = phos_sum[j,1],                      # Posterior mean
      Rhat = phos_sum[j,8],                      # Rhat values
      n.eff = phos_sum[j,9]                      # Sample size estimates (should be ~3000)
    ))
    
  }
  
  print(paste("Done", n, "of ", length(unique(df.lv$id.number))))
  
}

write.csv(levasseur.summ.df, "data-processed/503c_Levasseur_2025_phos_monods.csv") # Levasseur_2025 summary table
write.csv(fit.df, "data-processed/503d_Levasseur_2025_phos_monods_fits.csv") # Save model fit summary table

# Plotting & troubleshooting ----------------------------------------------

p.b <- ggplot(data=bestion.summ.df, aes(x = T.opt.TPC, y = T.max.raw)) +
  geom_point() +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  #ylim(10,35) +
  #xlim(10,35) +
  labs(
    x = "Topt (from TPC)",
    y = "T at max growth (raw data)",
    title = "Bestion 2018 - 5/6 TPCs fit"
  )

p.lv <- ggplot(data=levasseur.summ.df, aes(x = T.opt.TPC, y = T.max.raw)) +
  geom_point() +
  theme_classic() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  #ylim(5,35) +
  #xlim(5,35) +
  labs(
    x = "Topt (from TPC)",
    y = "T at max growth (raw data)",
    title = "Levasseur 2025 - 20/20 TPCs fit"
  )

phos.ts <- plot_grid(p.b, p.lv,
                    nrow = 1,
                    align = "hv")
