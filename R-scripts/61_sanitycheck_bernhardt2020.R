# Jason R Laurich

# November 2nd, 2025

# Moving towards the end of this project — submission soon hopefully! I want to make sure our data makes sense in the light
# of previous papers that used the same data (metanalyses and Joey's 2020 paper)!

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(cowplot)
library(nls.multstart)

# Load & examine the data -------------------------------------------------

df.hpd <- read.csv("data-processed/20_summary_table.csv") # Summary file, with the confidence intervals
str(df.hpd)
head(df.hpd)

df <- read.csv("data-processed/53_summary_table.csv") # Summary file, no hpds. Re-estimated for each replicate?
str(df)
head(df)

# (a) Evolutionary changes in R* and salt tolerance -----------------------

# We need to calculate changes in values, and then plot them. Also determine whether or not estimate diverge significantly from ancestral counterparts.

df.sum <- df.hpd %>%
  mutate(
    I.star = 1 / I.comp.0.56,                            # calculate I* for each population
    I.anc = 1 / I.comp.0.56[ match(Anc, Pop.fac) ],      # Extract I* of the ancestor
    d.I.star = I.star - I.anc,                           # Calculate difference for the mean.
    d.I.star.L = (1 / I.comp.0.56.L) - I.anc,            # Calculate CIs for the difference
    d.I.star.U = (1 / I.comp.0.56.U) - I.anc,
    
    N.star = 1 / N.comp.0.56,                            # same for N*
    N.anc = 1 / N.comp.0.56[ match(Anc, Pop.fac) ],      
    d.N.star = N.star - N.anc,                          
    d.N.star.L = (1 / N.comp.0.56.L) - N.anc,            
    d.N.star.U = (1 / N.comp.0.56.U) - N.anc,
    
    P.star = 1 / P.comp.0.56,                            # and P*
    P.anc = 1 / P.comp.0.56[ match(Anc, Pop.fac) ],      
    d.P.star = P.star - P.anc,                          
    d.P.star.L = (1 / P.comp.0.56.L) - P.anc,            
    d.P.star.U = (1 / P.comp.0.56.U) - P.anc,
    
    S.c.anc = S.c[ match(Anc, Pop.fac) ],                # and salt tolerance (c) 
    d.S.c = S.c - S.c.anc,                          
    d.S.c.L = S.c.L - S.c.anc,            
    d.S.c.U = S.c.U - S.c.anc
  ) %>%
  select(Pop.fac, Anc, Evol, I.star, d.I.star, d.I.star.L, d.I.star.U, I.µ.max, 
         N.star, d.N.star, d.N.star.L, d.N.star.U, N.µ.max,
         P.star, d.P.star, d.P.star.L, d.P.star.U, P.µ.max,
         S.c, d.S.c, d.S.c.L, d.S.c.U, S.µ.max)

# Recreate Figure 2 - one panel at a time

df.sum <- df.sum %>%
  mutate(Evol = if_else(Evol == "none", "A", Evol)) %>% 
  mutate(Evol = factor(Evol, levels = c("A","C","L","N","P","B","S","BS")))

# Panel A - change in P* ~ evol condition. 

p.P <- ggplot(df.sum, aes(x = Evol, y = d.P.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +

  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * P^"*" ~ "(µM P)")  # or I^"*", N^"*", etc.
  ) +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.P

# Panel B - change in N* ~ evol condition. 

p.N <- ggplot(df.sum, aes(x = Evol, y = d.N.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * N^"*" ~ "(µM N)")  # or I^"*", N^"*", etc.
  ) +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.N

# Panel C - change in I* ~ evol condition. 

p.I <- ggplot(df.sum, aes(x = Evol, y = d.I.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * I^"*" ~ "(µmol m-2s-1)")  # or I^"*", N^"*", etc.
  ) +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.I

# Panel D - change in N* ~ evol condition. 

p.S <- ggplot(df.sum, aes(x = Evol, y = d.S.c, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in salt tolerance" ~ "(g/L)")  # or I^"*", N^"*", etc.
  ) +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.S

fig <- plot_grid(p.P, p.N, p.I, p.S, nrow = 2, align='hv', rel_widths = c(1,1))
fig

# Let's look at the growth rates. Is this explaining it? ------------------

df.joey.p <- read.csv("data-processed/100_P_growth_Joey.csv")  # Joey's estimates of µ
head(df.joey.p)

df.joey.p$phos.lvl <- df.joey.p$phosphate_concentration 
df.joey.p$well.ID <- df.joey.p$well_plate

df.jason.p <- read.csv("data-processed/08a_µ_estimates_phosphorous.csv")  # My estimates of µ
head(df.jason.p)

df.growth <- df.jason.p %>% 
  left_join(df.joey.p %>% 
              select(population, well.ID, phos.lvl, estimate),
            by = c("population", "well.ID", "phos.lvl"))


df.sum$population <- df.sum$Pop.fac

df.growth2 <- df.growth %>% 
  left_join(df.sum %>% 
              select(Evol, population), by = "population")

reg <- ggplot(df.growth2, aes(x= r.exp, y =estimate, colour = phos.lvl)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

reg

# Let's look at N?

df.joey.n <- read.csv("data-processed/101_N_growth_Joey.csv")  # Joey's estimates of µ
head(df.joey.n)

df.joey.n$nitrate.lvl <- df.joey.n$nitrate_concentration 
df.joey.n$well.ID <- df.joey.n$well_plate

df.jason.n <- read.csv("data-processed/07a_µ_estimates_nitrogen_old.csv")  # My estimates of µ
head(df.jason.n)

df.growth.n <- df.jason.n %>% 
  left_join(df.joey.n %>% 
              select(population, well.ID, nitrate.lvl, estimate),
            by = c("population", "well.ID", "nitrate.lvl"))

reg2 <- ggplot(df.growth.n, aes(x= r.exp, y =estimate, colour = nitrate.lvl)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

reg2

nold <- lm(estimate ~ r.exp, data= df.growth.n)
summary(nold)

# Let's look at L?

df.joey.l <- read.csv("data-processed/102_L_growth_Joey.csv")  # Joey's estimates of µ
head(df.joey.l)

df.joey.l$percentage <- df.joey.l$light/2.5
df.joey.l$well.ID <- df.joey.l$well_plate

df.jason.l <- read.csv("data-processed/06a_µ_estimates_light.csv")  # My estimates of µ
head(df.jason.l)

df.growth.l <- df.jason.l %>% 
  left_join(df.joey.l %>% 
              select(population, well.ID, percentage, estimate),
            by = c("population", "well.ID", "nitrate.lvl"))

reg2 <- ggplot(df.growth.n, aes(x= r.exp, y =estimate, colour = nitrate.lvl)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

reg2

nold <- lm(estimate ~ r.exp, data= df.growth.n)
summary(nold)

# (b) The structure of trait variance: observed trade-offs ----------------

# Need to create Figure 4 a-c

# relationships between µmax and R*s

# Panel A - Phosphorous

p.P2 <- ggplot(df.sum, aes(x = P.star, y = P.µ.max)) +
  
  geom_point() +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.P2  

p.N2 <- ggplot(df.sum, aes(x = N.star, y = N.µ.max)) +
  
  geom_point() +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.N2  

p.I2 <- ggplot(df.sum, aes(x = I.star, y = I.µ.max)) +
  
  geom_point() +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.I2 

fig2 <- plot_grid(p.P2, p.N2, p.I2, nrow = 1, align='hv', rel_widths = c(1,1))
fig2

# (c) Correlations in changes across traits -------------------------------

# Panel a - d.P* v d.I*

p.PI <- ggplot(df.sum, aes(x = d.I.star, y = d.P.star)) +
  
  geom_point() +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.PI 

# Panel b - d.N* ~ d. I*

p.NI <- ggplot(df.sum, aes(x = d.I.star, y = d.N.star)) +
  
  geom_point() +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.NI 

# Panel b - d.N* ~ d. I*

p.PN <- ggplot(df.sum, aes(x = d.N.star, y = d.P.star)) +
  
  geom_point() +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.PN

fig3 <- plot_grid(p.PI, p.NI, p.PN, nrow = 1, align='hv', rel_widths = c(1,1))
fig3

# OK let's dive deeper - refitting µmax without my fancy cutoff approach --------

# Will this now faithfully reproduce Joey's results?

# Start with P. 

# Upload & examine data ---------------------------------------------------

df <- read.csv("data-processed/08_phosphorous_rfus_time.csv")
head(df) # RFU is density, days is time, phosphate_level is a factor (1 to 10). phosphate_concentration is what we want
str(df)

df<-df[,-c(1,2,4:6,8,10,11,13,14,16,17)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df$phos.conc <- as.numeric(df$phosphate_concentration)
df$phos_level <- factor(df$phosphate_level, levels = sort(unique(df$phosphate_level)), ordered = TRUE) # Keep the numerical sorting.
df$well_plate <- as.factor(df$well_plate)

df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac) # I don't recognize the COMBO group, I'm guessing this is a control?

df.exp <- subset(df, df$pop.fac != "COMBO") 
df.exp$well.ID<-as.factor(df.exp$well_plate)

N0.df <- df.exp %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df.exp <- df.exp %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.exp <- split(df.exp, df.exp$pop.num)  # Each element is a data frame for one population in df.exp

# Estimate µ --------------------------------------------------------------

df.r.exp <- data.frame( # Initializing a dataframe to store the results for each well, pop, and phosphorous level
  population = character(),
  population.number = numeric(),
  phos.lvl = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

phos <- as.vector(as.numeric(as.character(unique(df.exp$phos.conc)))) # for looping through nitrate levels
ord.phos<- sort(phos)

for (i in 1:length(mat.exp)){ # Looping through all of the populations
  
  for (t in ord.phos){ # and all of the phosphate levels
    
    df.it <- subset(mat.exp[[i]], phos.conc==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$days), ] # Just in case the phosphate data is also disordered. 
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
      s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
      
      df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
      # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl.th$pop.fac[1],          
          population.number = df.it.wl$pop.num[1],      
          nitrate.lvl = df.it.wl$nitrate.conc[1],        
          well.ID = df.it.wl$well.ID[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl$pop.fac[1],          # Population as factor
          population.number = df.it.wl$pop.num[1],   # Numeric population number
          phos.lvl = df.it.wl$phos.conc[1],          # Phosphorous level
          well.ID = df.it.wl$well.ID[1],             # Well ID (assuming one well ID per subset)
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
        
      }
      
    }
    
  }
  
}

# Save the files, with different subsets of data considered (e.g. the s parameter)
# write.csv(df.r.exp, "data-processed/200_µ_estimates_phosphorous.3-7.csv") # let's save the file.
# write.csv(df.r.exp, "data-processed/200_µ_estimates_phosphorous.2-7.csv") # let's save the file.

df.growth.p37 <- df.r.exp %>% 
  left_join(df.joey.p %>% 
              select(population, well.ID, phos.lvl, estimate),
            by = c("population", "well.ID", "phos.lvl"))

reg.p37 <- ggplot(df.growth.p37, aes(x= r.exp, y =estimate, colour = phos.lvl)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

reg.p37

df.growth.p27 <- df.r.exp %>% 
  left_join(df.joey.p %>% 
              select(population, well.ID, phos.lvl, estimate),
            by = c("population", "well.ID", "phos.lvl"))

reg.p27 <- ggplot(df.growth.p27, aes(x= r.exp, y =estimate, colour = phos.lvl)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

reg.p27

# Now let's move on to nitrogen -------------------------------------------

df <- read.csv("data-processed/07_nitrogen_rfus_time.csv")
head(df) # RFU is density, days is time, nitrate_level is a factor (1 to 10). nitrate_concentration is what we want
str(df)

df<-df[,-c(1,2,4,5,8:10,12:14)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df$nitrate.conc <- as.numeric(df$nitrate_concentration)
df$nitrate_level <- factor(df$nitrate_level, levels = sort(unique(df$nitrate_level)), ordered = TRUE) # Keep the numerical sorting.
df$well_plate <- as.factor(df$well_plate)

df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac) # Removing the COMBO treatment, which is simply a control. 

df.exp <- subset(df, df$pop.fac != "COMBO")
df.exp$well.ID<-as.factor(df.exp$well_plate)

N0.df <- df.exp %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df.exp <- df.exp %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.exp <- split(df.exp, df.exp$pop.num)  # Each element is a data frame for one population in df.exp

# OK let's fit µ ----------------------------------------------------------

df.r.exp <- data.frame( # Initializing a dataframe to store the results for each well, pop, and nitrogen level
  population = character(),
  population.number = numeric(),
  nitrate.lvl = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

nitrate <- as.vector(as.numeric(as.character(unique(df.exp$nitrate.conc)))) # for looping through nitrate levels
ord.nit<- sort(nitrate)

for (i in 1:length(mat.exp)){ # Looping through all of the populations
  
  for (t in ord.nit){ # and all of the nitrate levels
    
    df.it <- subset(mat.exp[[i]], nitrate.conc==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$days), ] # For whatever reason, the nitrogen data is not ordered properly by date in some cases.
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
      s <- which.max(ln.slopes[3:length(ln.slopes)])  # We need at least 3 data points
      
      df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s + 2], ] # Get the thresholded data according to our sliding window approach
      # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl.th$pop.fac[1],          
          population.number = df.it.wl$pop.num[1],      
          nitrate.lvl = df.it.wl$nitrate.conc[1],        
          well.ID = df.it.wl$well.ID[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl$pop.fac[1],          # Population as factor
          population.number = df.it.wl$pop.num[1],   # Numeric population number
          nitrate.lvl = df.it.wl$nitrate.conc[1],    # Nitrogen level
          well.ID = df.it.wl$well.ID[1],             # Well ID (assuming one well ID per subset)
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
        
      }
      
    }
    
  }
  
}

# write.csv(df.r.exp, "data-processed/201_µ_estimates_nitrogen_new.csv") # let's save the file.
write.csv(df.r.exp, "data-processed/201_µ_estimates_nitrogen_new_3-7.csv") # let's save the file.

df.growth.n27 <- df.r.exp %>% 
  left_join(df.joey.n %>% 
              select(population, well.ID, nitrate.lvl, estimate),
            by = c("population", "well.ID", "nitrate.lvl"))

reg.n27 <- ggplot(df.growth.n27, aes(x= r.exp, y =estimate, colour = nitrate.lvl)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

reg.n27

n27 <- lm(estimate ~ r.exp, data = df.growth.n27)
summary(n27)

df.growth.n37 <- df.r.exp %>% 
  left_join(df.joey.n %>% 
              select(population, well.ID, nitrate.lvl, estimate),
            by = c("population", "well.ID", "nitrate.lvl"))

reg.n37 <- ggplot(df.growth.n37, aes(x= r.exp, y =estimate, colour = nitrate.lvl)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

reg.n37

n37 <- lm(estimate ~ r.exp, data = df.growth.n37)
summary(n37)

# Now let's move on to light -------------------------------------------

df <- read.csv("data-processed/06_light_rfus_time.csv")
head(df) # RFU is density, days is time, light_level is a factor (1 to 10). Percentage is also a measurement of light I think?
str(df)

df<-df[,-c(1,2,4,5,8,10,12,13,15,16,17)]

head(df)

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df$percentage <- as.numeric(df$percentage) # From an examination of the csv, the light level 2 corresponds to "0.5-0.7"
df$percentage[is.na(df$percentage)] <- 0.6 # For now, let's set this to 0.6, but I need to talk with Joey about this. 
df$light_level <- factor(df$light_level, levels = sort(unique(df$light_level)), ordered = TRUE) # Keep the numerical sorting.
df$well_plate <- as.factor(df$well_plate)

# df$days <- df$days + 0.001 # Can't have 0s
df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac) # The COMBO group is a control, not relevant to this experiment. 

df.exp <- subset(df, df$pop.fac != "COMBO") 
df.exp$well.ID<-as.factor(df.exp$well_plate)

N0.df <- df.exp %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df.exp <- df.exp %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.exp <- split(df.exp, df.exp$pop.num)  # Each element is a data frame for one population in df.exp

# OK let's fit µ ----------------------------------------------------------

df.r.exp <- data.frame( # Initializing a dataframe to store the results for each well, pop, and nitrogen level
  population = character(),
  population.number = numeric(),
  nitrate.lvl = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

light <- as.vector(as.numeric(as.character(unique(df.exp$percentage)))) # for looping through light levels
ord.light<- sort(light)

for (i in 1:length(mat.exp)){ # Looping through all of the populations
  
  for (t in ord.light){ # and all of the light levels
    
    df.it <- subset(mat.exp[[i]], percentage==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$days), ] # For whatever reason, the nitrogen data is not ordered properly by date in some cases.
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
      s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
      
      df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
      # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl.th$pop.fac[1],          
          population.number = df.it.wl$pop.num[1],      
          percentage = df.it.wl$percentage[1],        
          well.ID = df.it.wl$well.ID[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl$pop.fac[1],          # Population as factor
          population.number = df.it.wl$pop.num[1],   # Numeric population number
          percentage = df.it.wl$percentage[1],       # Light level
          well.ID = df.it.wl$well.ID[1],             # Well ID (assuming one well ID per subset)
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
        
      }
      
    }
    
  }
  
}

write.csv(df.r.exp, "data-processed/202_µ_estimates_light_new.csv") # let's save the file.
# write.csv(df.r.exp, "data-processed/201_µ_estimates_light_new_3-7.csv") # let's save the file.

df.growth.l27 <- df.r.exp %>% 
  left_join(df.joey.l %>% 
              select(population, well.ID, percentage, estimate),
            by = c("population", "well.ID", "percentage"))

reg.l27 <- ggplot(df.growth.l27, aes(x= r.exp, y =estimate, colour = percentage)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

reg.l27

l27 <- lm(estimate ~ r.exp, data = df.growth.l27)
summary(l27)

# Now let's finish with salt -------------------------------------------

df <- read.csv("data-processed/09_salt_rfus_time.csv")
head(df) # RFU is density, time is time, treatment is the salt level s0 to S9.
str(df)

df<-df[,-c(2:4,8:9)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df$salt_level <- factor(df$treatment, levels = sort(unique(df$treatment)), ordered = TRUE) # Keep the numerical sorting.
df$salt <- as.numeric(df$salt_level) # These are just numbers for now reflecting the ordered factor.
df$salt <- df$salt - 1 # Convert this to proper levels (0-9)
df$days <- df$time

df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac)

df.exp <- df
df.exp$well.ID<-as.factor(df.exp$unique_id)

N0.df <- df.exp %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(time)]) %>%
  ungroup()

df.exp <- df.exp %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.exp <- split(df.exp, df.exp$pop.num)  # Each element is a data frame for one population in df.exp

# OK let's fit µ ----------------------------------------------------------

df.r.exp <- data.frame( # Initializing a dataframe to store the results for each well, pop, and nitrogen level
  population = character(),
  population.number = numeric(),
  nitrate.lvl = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

salt <- as.vector(as.numeric(as.character(unique(df.exp$salt)))) # for looping through nitrate levels
ord.salt<- sort(salt)

for (i in 1:length(mat.exp)){ # Looping through all of the populations
  
  for (t in ord.salt){ # and all of the salt levels
    
    df.it <- subset(mat.exp[[i]], salt==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$days), ] # For whatever reason, the nitrogen data is not ordered properly by date in some cases.
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
      s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
      
      df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
      # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl.th$pop.fac[1],          
          population.number = df.it.wl$pop.num[1],      
          salt = df.it.wl$salt[1],        
          well.ID = df.it.wl$well.ID[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl$pop.fac[1],          # Population as factor
          population.number = df.it.wl$pop.num[1],   # Numeric population number
          salt = df.it.wl$salt[1],                   # Salt level
          well.ID = df.it.wl$well.ID[1],             # Well ID (assuming one well ID per subset)
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
        
      }
      
    }
    
  }
  
}

write.csv(df.r.exp, "data-processed/203_µ_estimates_salt_new.csv") # let's save the file.
# write.csv(df.r.exp, "data-processed/201_µ_estimates_light_new_3-7.csv") # let's save the file.

df.growth.s27 <- df.r.exp %>% 
  left_join(df.joey.s %>% 
              select(population, well.ID, percentage, estimate),
            by = c("population", "well.ID", "percentage"))

reg.l27 <- ggplot(df.growth.l27, aes(x= r.exp, y =estimate, colour = percentage)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)

reg.l27

l27 <- lm(estimate ~ r.exp, data = df.growth.l27)
summary(l27)

# OK let's recreate Fig 2 with the new Monod data -------------------------

df.p <- read.csv("data-processed/210_Monod_phosphorous_pops_estimates.csv") # Summary file, with the confidence intervals
str(df.p)
head(df.p)

df.p <- df.p %>% 
  left_join(df.hpd %>% 
              select(Anc, Evol, Pop.fac), by = "Pop.fac")

df.p <- df.p %>%
  mutate(Evol = if_else(Evol == "none", "A", Evol)) %>% 
  mutate(Evol = factor(Evol, levels = c("A","C","L","N","P","B","S","BS")))

df.p <- df.p %>%
  mutate(
    P.star = R.mth,                                      # P* for each population
    P.anc = P.star[ match(Anc, Pop.fac) ],               # Extract P* of the ancestor
    d.P.star = P.star - P.anc                         # Calculate difference for the mean.
  ) %>%
  select(Pop.fac, Anc, Evol, P.star, d.P.star, r.max)

# Panel A - change in P* ~ evol condition. 

p.P.new <- ggplot(df.p, aes(x = Evol, y = d.P.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * P^"*" ~ "(µM P)")  # or I^"*", N^"*", etc.
  ) +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.P.new

# OK, Nitrogen
 
df.n <- read.csv("data-processed/212_Monod_nitrogen_pops_estimates.csv") # Summary file
str(df.n)
head(df.n)

df.n <- df.n %>% 
  left_join(df.hpd %>% 
              select(Anc, Evol, Pop.fac), by = "Pop.fac")

df.n <- df.n %>%
  mutate(Evol = if_else(Evol == "none", "A", Evol)) %>% 
  mutate(Evol = factor(Evol, levels = c("A","C","L","N","P","B","S","BS")))

df.n <- df.n %>%
  mutate(
    N.star = R.mth,                                      # N* for each population
    N.anc = N.star[ match(Anc, Pop.fac) ],               # Extract N* of the ancestor
    d.N.star = N.star - N.anc                            # Calculate difference for the mean.
  ) %>%
  select(Pop.fac, Anc, Evol, N.star, d.N.star, r.max)

# Panel B - change in N* ~ evol condition. 

p.N.new <- ggplot(df.n, aes(x = Evol, y = d.N.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * P^"*" ~ "(µM P)")  # or I^"*", N^"*", etc.
  ) +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.N.new

# Now, light

df.l <- read.csv("data-processed/214_Monod_light_pops_estimates.csv") # Summary file
str(df.l)
head(df.l)

df.l <- df.l %>% 
  left_join(df.hpd %>% 
              select(Anc, Evol, Pop.fac), by = "Pop.fac")

df.l <- df.l %>%
  mutate(Evol = if_else(Evol == "none", "A", Evol)) %>% 
  mutate(Evol = factor(Evol, levels = c("A","C","L","N","P","B","S","BS")))

df.l <- df.l %>%
  mutate(
    L.star = R.mth,                                      # L* for each population
    L.anc = L.star[ match(Anc, Pop.fac) ],               # Extract L* of the ancestor
    d.L.star = L.star - L.anc                            # Calculate difference for the mean.
  ) %>%
  select(Pop.fac, Anc, Evol, L.star, d.L.star, r.max)

# Panel C - change in I* ~ evol condition. 

p.L.new <- ggplot(df.l, aes(x = Evol, y = d.L.star, colour = Evol)) +
  
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey30") +
  
  geom_point(position = position_jitter(width = 0.15, height = 0),
             size = 2.5) +
  
  scale_colour_manual(values = c(
    A  = "black",
    C  = "forestgreen",
    L  = "gold",
    N  = "magenta3",
    P  = "firebrick",
    B  = "darkorange",
    S  = "deepskyblue1",
    BS = "blue"
  ), guide = "none") +   # drop colour legend if you don’t want it
  
  labs(
    x = NULL,
    y = expression("change in " * P^"*" ~ "(µM P)")  # or I^"*", N^"*", etc.
  ) +
  
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x  = element_text(size = 10),
    axis.text.y  = element_text(size = 10)
  )

p.L.new

fig.new <- plot_grid(p.P.new, p.N.new, p.I.new, nrow = 2, align='hv', rel_widths = c(1,1))
fig.new

# OK let's try figure 4 again ---------------------------------------------

p.P2.new <- ggplot(df.sum, aes(x = P.star, y = P.µ.max)) +
  
  geom_point() +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.P2.new  

p.N2.new <- ggplot(df.n, aes(x = N.star, y = r.max)) +
  
  geom_point() +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.N2.new  

p.I2 <- ggplot(df.sum, aes(x = I.star, y = I.µ.max)) +
  
  geom_point() +
  
  geom_smooth(method = 'lm') +
  
  theme_classic()

p.I2 

fig2 <- plot_grid(p.P2, p.N2, p.I2, nrow = 1, align='hv', rel_widths = c(1,1))
fig2

# Second last: temperature -------------------------------------------

df <- read.csv("data-processed/00_temp_rfus_time.csv")
head(df) #RFU is density, days is time, temperature is temperature
str(df)

df<-df[,-c(1,2,4,5,6,8,13)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df<-subset(df,df$temperature!=20) # 20C only includes 'single' measurements, which we will be excluding.

# df$days <- df$days + 0.001 # Can't have 0s
df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac) # The cc1629 group is not relevant to the rest of the experimental data we are looking at, but I will keep it for now. 

df.rep <- subset(df, df$plate_type == "repeat") # We're going to work only with repeat data (this is most of the well_plates) which I'll use as replicates. 
df.rep$well.ID<-as.factor(df.rep$well_plate)

N0.df <- df.rep %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df.rep <- df.rep %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.rep <- split(df.rep, df.rep$pop.num)  # Each element is a data frame for one population in df.rep contained within a matrix

# OK let's fit µ ----------------------------------------------------------

df.r.exp <- data.frame( # Initializing a dataframe to store the results for each well, pop, and temp level
  population = character(),
  population.number = numeric(),
  nitrate.lvl = numeric(),
  well.ID = character(),
  r.exp = numeric()
)

tmp <- as.vector(as.numeric(as.character(unique(df.rep$temperature)))) # for looping through temperatures
ord.tmp<- sort(tmp)

for (i in 1:length(mat.rep)){ # Looping through all of the populations
  
  for (t in ord.tmp){ # and all of the temp levels
    
    df.it <- subset(mat.rep[[i]], temperature==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$days), ] # For whatever reason, the nitrogen data is not ordered properly by date in some cases.
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x N x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
      s <- which.max(ln.slopes[2:length(ln.slopes)])  # We need at least 3 data points
      
      df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s + 1], ] # Get the thresholded data according to our sliding window approach
      # The + 1 here just corrects for the labelling mismatch (e.g. in the above line, s will return 1 when the 2nd slope is the highest)
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      if (is.null(r_exp)){
        
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl.th$pop.fac[1],          
          population.number = df.it.wl$pop.num[1],      
          temperature = df.it.wl$temperature[1],        
          well.ID = df.it.wl$well.ID[1],                
          r.exp = NA        
        ))
        
      }else{
        # Add data to our summary table
        df.r.exp <- rbind(df.r.exp, data.frame(
          population = df.it.wl$pop.fac[1],          # Population as factor
          population.number = df.it.wl$pop.num[1],   # Numeric population number
          temperature = df.it.wl$temperature[1],     # Temp level
          well.ID = df.it.wl$well.ID[1],             # Well ID (assuming one well ID per subset)
          r.exp = summary(r_exp)$parameters[1,1]     # The calculated r.exp value
        ))
        
      }
      
    }
    
  }
  
}

write.csv(df.r.exp, "data-processed/204_µ_estimates_temp_new.csv") # let's save the file.
# write.csv(df.r.exp, "data-processed/201_µ_estimates_light_new_3-7.csv") # let's save the file.
