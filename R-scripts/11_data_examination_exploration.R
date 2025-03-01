# Jason R Laurich
# February 24, 2025

# OK so here I am going upload data for TPC shapes, Light, Nitrogen, and Phosphorous Monod curves, and Salt tolerance data.

# Additionally, I'll look at the data closely and troubleshoot (e.g. looking for outliers, figuring out why µmax is so much higher for TPCs)
# And I'll do a PCA and RDA of these data, generating some plots.

# Originally was going to do analysis here too, but decided to move that over to a separate script (12)

############# Packages ########################

library(dplyr)
library(ggplot2)
library(rPref)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(ggrepel)

############# Upload and organize data #######################

df.hist <- read.csv("data-processed/chlamee-acute-exponential.csv") # This file has the history (ancestry, evolutionary treatment) for each well ID.
head(df.hist)
str(df.hist)
df.hist$anc <- as.factor(df.hist$ancestor_id) # Will port these into combined dataset, remember: we are sorting by population. 
df.hist$evol <- as.factor(df.hist$treatment)  # so we really only need population level - evol - anc mapping, not well-specific.
df.hist$Pop.fac <- as.factor(df.hist$population)

df.tpc <- read.csv("data-processed/09_TPC_shape_values_BayesTPC.csv") # TPC summary data. We'll use Lactin II model fits. 
head(df.tpc)
str(df.tpc)

df.tpc <- df.tpc[df.tpc$Pop.fac != 'cc1629' & df.tpc$Model == 'Lactin 2', ] # Trim off Thomas models and that one weird population.

df <- df.tpc[,c(2,5:9)] # We'll work with the analytically determined parameters (Deriv package)

df.hist.agg <- df.hist %>% # Aggregate the history data for merging.
  group_by(Pop.fac) %>%
  summarise(
    anc = first(anc),
    evol = first(evol),
    .groups = "drop" 
  )

df <- merge(df, df.hist.agg, by = "Pop.fac", all.x = TRUE)
head(df) 
str(df)
df$Pop.fac <- as.factor(df$Pop.fac) # Looks great! Now we layer in other summary metrics.
names(df) <- c("Pop.fac", "T.min", "T.max", "T.br", "Topt", "r.max_T", "anc", "evol") # rename a few variables. 

df.I <- read.csv("data-processed/10b_light_Monod_estimates.csv") # I* data 
head(df.I)
str(df.I)

df.I.par <- df.I[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.I.par) <- c("Pop.fac", "K.s_I", "r.max_I", "I.jag", "I.mth") # rename
df.I.par$I.comp <- 1/df.I.par$I.mth

df <- merge(df, df.I.par, by = "Pop.fac", all.x = TRUE)

df.N <- read.csv("data-processed/11b_nitrogen_Monod_estimates.csv") # N* data 
head(df.N)
str(df.N)

df.N.par <- df.N[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.N.par) <- c("Pop.fac", "K.s_N", "r.max_N", "N.jag", "N.mth") # rename
df.N.par$N.comp <- 1/df.N.par$N.mth

df <- merge(df, df.N.par, by = "Pop.fac", all.x = TRUE)

df.P <- read.csv("data-processed/12b_phosphorous_Monod_estimates.csv") # P* data 
head(df.P)
str(df.P)

df.P.par <- df.P[, c("Pop.fac", "K.s", "r.max", "R.jag", "R.mth")] # need these
names(df.P.par) <- c("Pop.fac", "K.s_P", "r.max_P", "P.jag", "P.mth") # rename
df.P.par$P.comp <- 1/df.P.par$P.mth

df <- merge(df, df.P.par, by = "Pop.fac", all.x = TRUE)

df.S <- read.csv("data-Processed/13b_salt_tolerance_estimates.csv") # S* data 
head(df.S)
str(df.S)

df.S$Pop.fac <- gsub("Anc ", "anc", df.S$Pop.fac) # ancestors are written differently here

df.S.par <- df.S[, c("Pop.fac", "r.max", "c.mod", "c.pred")] # need these
names(df.S.par) <- c("Pop.fac", "r.max_S", "S.c.mod", "S.c.pred") # rename

df <- merge(df, df.S.par, by = "Pop.fac", all.x = TRUE)

write.csv(df, "data-processed/14_summary_metric_table.csv") # Save summary table

############# Let's visualize some trade-offs #######################

df<-read.csv("data-processed/14_summary_metric_table.csv") # To start back up from this point forward

df$shape <- ifelse(df$evol == "none", 22, 16) # I want to add a shape column to the dataframe that I will update
# The idea is to label un-evolved populations with a square, and then later (not for T) relevant experimental evolution nutrient conditions with a star

T_gs <- ggplot(df, aes(x = r.max_T, y = T.br)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Thermal breadth", 
       title = "Thermal performance") +
  scale_shape_manual(values = c(16, 3)) +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

T_gs # Negative relationship, data looks good

# Phosphorous

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "P", 8, 16)) # Ps are now equivalent to 8, for later mapping

P_go <- ggplot(df, aes(x = r.max_P, y = P.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Phosphorous limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

P_go # Negative correlation, one outlier (upper left) to look into.

new_data <- data.frame(r.max_P = seq(min(par.res$r.max_P), max(par.res$r.max_P), length.out = 100))
preds <- as.data.frame(fitted(brm_ns, newdata = new_data, probs = c(0.05, 0.95))) 

# Nitrogen

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "N", 8, 16)) # Ps are now equivalent to 8, for later mapping

N_go <- ggplot(df, aes(x = r.max_N, y = N.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Nitrogen limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

N_go # No evidence of a trade-off really. One data point in the top right? Not much interesting going on here. 

# Light

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "L", 8, 16)) # Ls are now equivalent to 8, for later mapping

L_go <- ggplot(df, aes(x = r.max_I, y = I.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Light limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

L_go # One points seems WAAY off. Look into that. Replot without it? I don't believe that data.

L_go.1 <- ggplot(df[df$I.comp<10,], aes(x = r.max_I, y = I.comp)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Competitive ability (1/R*)", 
       title = "Light limitation") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

L_go.1 # Looks more interesting. A bit of a positive correlation?

# Salt

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol == "S", 8, 16)) # Ls are now equivalent to 8, for later mapping

S_gs <- ggplot(df, aes(x = r.max_S, y = S.c.mod)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Salt tolerance (c)", 
       title = "Salt stress") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

S_gs # Looks negative

unique(df$evol)

df$shape <- ifelse(df$evol == "none", 22, 
                   ifelse(df$evol %in% c("S", "BS"), 8, 16)) # add BS to our list here?

S_gs.1 <- ggplot(df, aes(x = r.max_S, y = S.c.mod)) +  # Remove shape from aes() for regression
  geom_point(aes(shape = as.factor(shape)), size = 2) +  # Keep shape only for points
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1, linetype = "dashed") +  # Single regression
  labs(x = "Maximum exponential growth rate", 
       y = "Salt tolerance (c)", 
       title = "Salt stress") +
  scale_shape_manual(values = c("16" = 16, "22" = 3, "8" = 8)) +  # Assign stars to 8, circles to 16, pluses to 22 +  # Keep custom shapes
  theme_classic() +
  theme(legend.position = "none")  # Remove legend

S_gs.1 # LHmm this is quite compelling - all of the derived populations are on the top right!

################ Troubleshooting & Investigating ################################

plot_grid(T_par, P_par, N_par, S_par, L_par) # Different max growth rates being observed across experiments! TPC > salt > nutrients..

# Is this real? Let's load the raw data and find out.

df.t <- read.csv("data-processed/chlamee-acute-rfu-time.csv")
head(df.t)

df.l <- read.csv("data-processed/10_light_rstar_rfus_time.csv")
head(df.l)

df.n <- read.csv("data-processed/11_nitrate_abundances_processed.csv")
head(df.n)

df.p <- read.csv("data-processed/12_phosphate_rstar_rfus_time.csv")
head(df.p)

df.s <- read.csv("data-processed/13_chlamee_salt_rfus_time.csv")
head(df.s)

# OK let's make a grid plot with the raw data to look at it (< 4 days to capture the exponential growth phase of most samples)

p.t <- ggplot(df.t[df.t$days<4,], aes(x=days, y=RFU)) +
  geom_point() +
  geom_smooth(method='lm', colour='red') +
  ylim(c(0,2500)) + # Some weird outliers, we'll keep this constant
  theme_classic() +
  labs(title = 'Temperature')

p.l <- ggplot(df.l[df.l$days<4,], aes(x=days, y=RFU)) +
  geom_point() +
  geom_smooth(method='lm', colour='red') +
  ylim(c(0,2500)) + # Some weird outliers, we'll keep this constant
  theme_classic() +
  labs(title = 'Light')

p.n <- ggplot(df.n[df.n$days<4,], aes(x=days, y=RFU)) +
  geom_point() +
  geom_smooth(method='lm', colour='red') +
  ylim(c(0,2500)) + # Some weird outliers, we'll keep this constant
  theme_classic() +
  labs(title = 'Nitrogen')

p.p <- ggplot(df.p[df.p$days<4,], aes(x=days, y=RFU)) +
  geom_point() +
  geom_smooth(method='lm', colour='red') +
  ylim(c(0,2500)) + # Some weird outliers, we'll keep this constant
  theme_classic() +
  labs(title = 'Phosphorous')

p.s <- ggplot(df.s[df.s$time<4,], aes(x=time, y=RFU)) +
  geom_point() +
  geom_smooth(method='lm', colour='red') +
  ylim(c(0,2500)) + # Some weird outliers, we'll keep this constant
  theme_classic() +
  labs(title = 'Salt')

plot_grid(p.t, p.l, p.n, p.p, p.s, nrow = 1)

# OK let's look into outliers

#Phosphorous - 1 data point has very high P.comp
P_par # checked data, this is population 35 (pop.fac) — corresponds to pop.num 26

ran <- sample(c(1:25,27:37), 5) # five random populations to compare it to.

for (i in c(26,ran)){
  load(paste0("R2jags-objects/pop_", i, "_phosphorous_monod.RData"))
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:3,2005),]
  df.jags.plot$phos <- seq(0, 50, 0.025)
  assign(paste0("df.P.jags", i), df.jags.plot)
}

P.comp.plot <- ggplot(data = df.P.jags26, aes(x = phos)) +
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  
  geom_line(data = df.P.jags8, aes(x = phos, y = mean), color = "darkorange", size = 1) + 
  
  geom_line(data = df.P.jags9, aes(x = phos, y = mean), color = "forestgreen", size = 1) + 
  
  geom_line(data = df.P.jags18, aes(x = phos, y = mean), color = "goldenrod1", size = 1) + 
  
  geom_line(data = df.P.jags22, aes(x = phos, y = mean), color = "dodgerblue3", size = 1) + 
  
  geom_line(data = df.P.jags23, aes(x = phos, y = mean), color = "black", size = 1) + 
  
  scale_x_continuous(limits = c(0, 50)) + 
  scale_y_continuous(limits = c(0, 2)) + # Customize the axes and labels +
  labs(
    x = "Phosphorous concentration",
    y = "Growth rate",
    title = "P limitation, Monod"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

P.comp.plot # So it does actually seem like the 26th population shoots up, with a significantly higher P.comp/ lower R*

# Is this reflected in the raw data?

df.p.r <- read.csv("data-processed/12a_phosphorous_r_estimates.csv")
head(df.p.r)

df.p.sub <- df.p.r[df.p.r$population.number %in% c(26,ran) & df.p.r$phos.lvl < 5 ,]

ggplot(df.p.sub, aes(x=phos.lvl, r.exp, colour = as.factor(as.character(population.number)))) +
  geom_point() +
  theme_classic() +
  labs(color = "Pop") +  # Change legend title
  scale_color_manual(values = c('darkorange', 'forestgreen', 'goldenrod1', 'dodgerblue3', 'darkorchid', 'black')) +
  labs(title = 'P r.exp comparison') # OK so it seems like our r.exp metrices really are higher for population 26...

## OK so we need to refit the growth curves and plot the data to check this out here..
  
df.p<-df.p[,-c(1,2,4:6,8,10,11,13,14,16,17)]
df.p.sub2 <- df.p[df.p$phosphate_concentration < 5 ,]

df.p.sub2$pop.fac <- as.factor(df.p.sub2$population)
df.p.sub2$pop.num <- as.numeric(df.p.sub2$pop.fac)
df.p.sub2$phos.conc <- as.numeric(df.p.sub2$phosphate_concentration)
df.p.sub2$phos_level <- factor(df.p.sub2$phosphate_level, levels = sort(unique(df.p.sub2$phosphate_level)), ordered = TRUE) # Keep the numerical sorting.
df.p.sub2$well_plate <- as.factor(df.p.sub2$well_plate)

df.p.sub2$logRFU <- log(df.p.sub2$RFU + 0.001)

df.p.sub2$well.ID<-as.factor(df.p.sub2$well_plate)

N0.df <- df.p.sub2 %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df.p.sub2 <- df.p.sub2 %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.exp <- split(df.p.sub2, df.p.sub2$pop.num)  # Each element is a data frame for one population in df.p.sub2

pops<-names(mat.exp) # We're going to plot out the growth curves for 6 populations - all 4 phosphorous levels. 

plot.list <- list() # make a list of the plots

phos <- as.vector(as.numeric(as.character(unique(df.p.sub2$phos.conc)))) # for looping through nitrate levels
ord.phos<- sort(phos)

for (i in  c(26,ran)){ # Looping through chosen populations
  
  for (t in ord.phos){ # and all of the P levels
    
    df.it <- subset(mat.exp[[i]], phos.conc==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    pred.df.mod <- data.frame( # Dataframe to store fitted data for P levels growth curves. Will reset for each level, but that's OK, as we will save the plots first
      population.fac = character(),
      population.num = numeric(),
      phos.conc = numeric(),
      well.ID = character(),
      smt.days = numeric(),
      fit.RFU = numeric()
    )
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      df.it.wl <- df.it.wl[order(df.it.wl$days), ] # Not sure if this is needed for P, but can't hurt.
      
      if (df.it.wl$RFU[2] <  df.it.wl$RFU[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.wl <- df.it.wl[-1,]
        df.it.wl$N0 <- df.it.wl$RFU[1]
      }
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x P x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x P level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 10%
        
        if (ln.slopes[r] != 0 & round(ln.slopes[r], digits = 5) !=0) { # We also need to account for tiny values that are basically 0 (e.g. 5 e-16, but are messing up our loops)
          # percent.chg <- (ln.slopes[s] - ln.slopes[s + 1]) / ln.slopes[s] This was the reason we were getting weird negative results! If there was a stochastic drop.
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] # This should (along with the next line) fix it.
          
          if (percent.chg >= 1.10 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { # Because now, the drop-off ignores transiently negative slopes. 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
      
      df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      smt.days <- seq(min(df.it.wl.th$days, na.rm = TRUE), max(df.it.wl.th$days, na.rm = TRUE), length.out = 500) # Get a smooth distribution of time points
      
      # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
      fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
      pred.df.mod <- rbind(pred.df.mod, data.frame(
        population.fac = df.it.wl$pop.fac[1],
        population.num = i,
        phos.conc = t,
        well.ID = unique(df.it.wl$well.ID),
        smt.days = smt.days,
        fit.RFU = fit.RFU.mod
      ))
      
      # Now we plot it out
      
      p <- ggplot() + 
        geom_point(data=df.it, aes(x=days, y=RFU, colour= well.ID), size = 2) + # Observed data, not specific to individual well.IDs
        geom_line(data = pred.df.mod, aes(x = smt.days, y = fit.RFU, colour = well.ID), size = 1) + # Fitted line
        labs(title = paste("Pop", i, "[P]:", t), x = "Days", y = "RFU") +
        theme_minimal() +
        theme(legend.position = "none")
      
      plot.list[[paste0("Pop", df.it.wl.th$pop.fac[1], " Phosphate ", t)]] <- p
      
    }
  }
}

plot_grid <- arrangeGrob(grobs=plot.list, ncol = 4, nrow = 6)

grid.draw(plot_grid)  # OK this does seem real! The 26th population is reaching much higher densities of RFUs early on for the first 2 P concentrations.

# OK so now let's look at the light outlier, which seems crazy.Population # 31 (fac 8) has an R* of 0.05? Unlikely!

ran2 <- sample(c(1:30,32:37), 5) # five random populations to compare it to.

for (i in c(31,ran2)){
  load(paste0("R2jags-objects/pop_", i, "_light_monod.RData"))
  df.jags <- data.frame(monod_jag$BUGSoutput$summary)
  df.jags.plot <- df.jags[-c(1:3,2005),]
  df.jags.plot$light <- seq(0, 100, 0.05)
  assign(paste0("df.L.jags", i), df.jags.plot)
}

L.comp.plot <- ggplot(data = df.L.jags31, aes(x = light)) +
  geom_line(aes(y = mean), color = "darkorchid", size = 1) + # Add the mean prediction line
  
  geom_line(data = df.L.jags8, aes(x = light, y = mean), color = "darkorange", size = 1) + 
  
  geom_line(data = df.L.jags13, aes(x = light, y = mean), color = "forestgreen", size = 1) + 
  
  geom_line(data = df.L.jags15, aes(x = light, y = mean), color = "goldenrod1", size = 1) + 
  
  geom_line(data = df.L.jags25, aes(x = light, y = mean), color = "dodgerblue3", size = 1) + 
  
  geom_line(data = df.L.jags35, aes(x = light, y = mean), color = "black", size = 1) + 
  
  scale_x_continuous(limits = c(0, 100)) + 
  scale_y_continuous(limits = c(0, 2)) + # Customize the axes and labels +
  labs(
    x = "Light level",
    y = "Growth rate",
    title = "Light limitation, Monod"
  ) +
  theme_classic() +
  geom_hline(yintercept = 0)

L.comp.plot # So it does actually seem like the 26th population shoots up, with a significantly higher P.comp/ lower R*

df.l.sub<-df.l[,-c(1,2,4,5,8,10,12,13,15,16,17)]
df.l.sub$pop.fac <- as.factor(df.l.sub$population)
df.l.sub$pop.num <- as.numeric(df.l.sub$pop.fac)
df.l.sub$percentage <- as.numeric(df.l.sub$percentage) # From an examination of the csv, the light level 2 corresponds to "0.5-0.7"
df.l.sub$percentage[is.na(df.l.sub$percentage)] <- 0.6 # For now, let's set this to 0.6, but I need to talk with Joey about this. 
df.l.sub$light_level <- factor(df.l.sub$light_level, levels = sort(unique(df.l.sub$light_level)), ordered = TRUE) # Keep the numerical sorting.
df.l.sub$well_plate <- as.factor(df.l.sub$well_plate)

# df$days <- df$days + 0.001 # Can't have 0s
df.l.sub$logRFU <- log(df.l.sub$RFU + 0.001)

df.l.sub$well.ID<-as.factor(df.l.sub$well_plate)

N0.df <- df.l.sub %>% # We want to create an additional column that holds N0 data - we will pull this from the RFU column for the first time point for each well.ID
  group_by(well.ID) %>%
  summarize(N0 = RFU[which.min(days)]) %>%
  ungroup()

df.l.sub <- df.l.sub %>% # Recombine this with our dataframe
  left_join(N0.df, by = "well.ID") # Viewed it, looks good.

mat.exp <- split(df.l.sub, df.l.sub$pop.num)  # Each element is a data frame for one population in df.p.sub2

plot.list <- list() # make a list of the plots

light <- as.vector(as.numeric(as.character(unique(df.l.sub$percentage)))) # for looping through light levels
ord.light<- sort(light)

for (i in c(31, ran2)){ # Looping through all of the populations
  
  for (t in ord.light){ # and all of the light levels
    
    df.it <- subset(mat.exp[[i]], percentage==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    pred.df.mod <- data.frame( # Dataframe to store fitted data for light levels growth curves. Will reset for each level, but that's OK, as we will save the plots first
      population.fac = character(),
      population.num = numeric(),
      light.lvl = numeric(),
      well.ID = character(),
      smt.days = numeric(),
      fit.RFU = numeric()
    )
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      if (df.it.wl$RFU[2] <  df.it.wl$RFU[1]){ # Get rid of 1st data points where there is a drop off after the first observation
        df.it.wl <- df.it.wl[-1,]
        df.it.wl$N0 <- df.it.wl$RFU[1]
      }
      
      t.series <- unique(df.it.wl$days) # Re-initialize this internally - we will only save summary data for each unique pop x I x well combo
      t.series <- t.series[-1] # Trim off the first entry to make tracking easier.
      
      ln.slopes <- c() # Re-initialize this too!
      
      for (z in t.series){
        
        df.it.wl.sl <- df.it.wl[df.it.wl$days <= z, ] # Subset the data to exclude time points above our window
        
        ln_slope <- lm(logRFU~days, data = df.it.wl.sl)
        
        ln.slopes <- c(ln.slopes, summary(ln_slope)$coefficients[2,1])
        
      } # So now we have our slopes for each well.ID x Pop x Light level
      
      s <- length(ln.slopes) # Initialize the full length of the dataset, for cases in which the entire period is exponential growth.
      
      for (r in 1:(length(ln.slopes) - 1)) { # Loop through slopes to find when the drop-off exceeds 10%
        
        if (ln.slopes[r] != 0 & round(ln.slopes[r], digits = 5) !=0) { # We also need to account for tiny values that are basically 0 (e.g. 5 e-16, but are messing up our loops)
          # percent.chg <- (ln.slopes[s] - ln.slopes[s + 1]) / ln.slopes[s] This was the reason we were getting weird negative results! If there was a stochastic drop.
          percent.chg <- ln.slopes[r] / ln.slopes[r + 1] # This should (along with the next line) fix it.
          
          if (percent.chg >= 1.10 & ln.slopes[r] > 0 & ln.slopes[r+1] > 0) { # Because now, the drop-off ignores transiently negative slopes. 
            s <- r # If the condition is met, reassign s to the corresponding drop-off point!
            break  # Exit loop when condition is met
          }
        }
      } # Now I have s!
      
      df.it.wl.th <- df.it.wl[df.it.wl$days <= t.series[s], ] # Get the thresholded data according to our sliding window approach
      
      r_exp <- nls_multstart(RFU ~ N0 * exp(r*days),  # Exponential growth model (N0 is in our dataframe)
                             data = df.it.wl.th,
                             start_lower = c(r = -4.5), 
                             start_upper = c(r = 4.5),   
                             iter = 500,
                             supp_errors = 'Y',
                             control = nls.control(maxiter = 200))
      
      smt.days <- seq(min(df.it.wl.th$days, na.rm = TRUE), max(df.it.wl.th$days, na.rm = TRUE), length.out = 500) # Get a smooth distribution of time points
      
      # Predict fitted RFU values for the exponential growth model at each well.ID and store it for plotting
      fit.RFU.mod <- predict(r_exp, newdata = data.frame(days = smt.days))
      pred.df.mod <- rbind(pred.df.mod, data.frame(
        population.fac = df.it.wl$pop.fac[1],
        population.num = i,
        light.lvl = t,
        well.ID = unique(df.it.wl$well.ID),
        smt.days = smt.days,
        fit.RFU = fit.RFU.mod
      ))
      
      # Now we plot it out
      
      p <- ggplot() + 
        geom_point(data=df.it, aes(x=days, y=RFU, colour= well.ID), size = 2) + # Observed data, not specific to individual well.IDs
        geom_line(data = pred.df.mod, aes(x = smt.days, y = fit.RFU, colour = well.ID), size = 1) + # Fitted line
        labs(title = paste("Pop", i, "Light:", t), x = "Days", y = "RFU") +
        theme_minimal() +
        theme(legend.position = "none")
      
      plot.list[[paste0("Pop", df.it.wl.th$pop.fac[1], " Light", t)]] <- p
      
    }
  }
}

plot_grid <- arrangeGrob(grobs=plot.list, ncol = 5, nrow = 6)

grid.draw(plot_grid)   


df.l.sub2 <- df.l.sub[df.l.sub$percentage < 1,]
df.l.sub2 <- droplevels(df.l.sub2)
  
ggplot(df.l.sub2, aes(y=RFU, x=days, colour= as.factor(as.character(pop.num)))) +
  geom_point() +
  facet_wrap(~percentage) +
  theme_classic() # OK what, this population might actually be off the chain? Does it require no light basically?

#################### Conclusions #######################

# 1 Variation in µmax is fine
# 2 the P outlier (high P.comp) is gine
# 3 the I outlier (insanely high I.comp) is not — remove it for now.