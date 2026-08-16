# ============================================================
# 01_JB_v_JL_analysis.R
#
# Systematic comparison of Bernhardt (JB) vs Laurich (JL)
# phosphorus R* estimation pipelines.
#
# ============================================================

# Jason R Laurich

# July 14th, 2026

# Here we are going to investigate discrepancies between my analysis and Joey's original results (ProcB paper)

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(minpack.lm)
library(nlstools)
library(broom)
library(janitor)
library(readxl)
library(cowplot)

m <- 0.56   # mortality / dilution rate (day^-1)

# Load & examine the data -------------------------------------------------

###### JB data ######

df.jb.rfu <- read.csv("JB-v-JL/chlamee-r-star/data-processed/phosphate-rstar-rfus-time.csv") # JB raw RFU timeseries
head(df.jb.rfu)

df.trt.jb <- read_excel("JB-v-JL/chlamee-r-star/data-general/ChlamEE_Treatments_JB.xlsx") %>% clean_names() # JB treatment / ancestor key
head(df.trt.jb)

df.jb.pt <- read.csv("JB-v-JL/chlamee-r-star/data-processed/phosphate-rstars-direct.csv") # JB published P* values

df.jb.boot <- read.csv("JB-v-JL/chlamee-r-star/data-processed/change-pstar-monod-boot-0.56.csv") # boot strapped data

###### JL data ######

df.jl.rfu <- read.csv("processed-data/14_phosphorous_rfus_time.csv") # JL raw RFU timeseries (should be the same as df.jb.rfu)

df.jl.mu <- read.csv("processed-data/15_µ_estimates_phosphorous.csv") # JL per-replicate mu estimates (estimates of growth rates)

df.jl.bayes <- read.csv("processed-data/16_phos_monod_summary.csv") # JL Bayesian Monod fits (per replicate)

df.jl.sum <- read.csv("processed-data/27_summary_table.csv") # Summary file
head(df.jl.sum)

# Section 1: replicate JB results exactly ---------------------------------
# Direct fitting of Monod curves using pooled replicates in nls, bootstrapping.  

df.jb.rfu <- df.jb.rfu %>%
  filter(population != "COMBO",
         !is.na(RFU)) %>%
  mutate(population = as.character(population)) %>%
  arrange(well_plate, days)

df.jb.rfu <- df.jb.rfu %>%  # Extract N0 (first estimate of RFUs) per pop. 
  group_by(well_plate) %>%
  mutate(N0 = first(RFU)) %>%
  ungroup() %>%
  group_by(population) %>%
  mutate(N0_mean = mean(N0)) %>%
  ungroup()

slopes.all <- data.frame() # estimate the log-linear slopes to identify the exponential phase

for (n in 3:7) {
  slopes.n <- df.jb.rfu %>%
    
    group_by(well_plate) %>%
    slice_min(order_by = days, n = n, with_ties = FALSE) %>%
    filter(RFU > 0) %>%
    summarise(slope = coef(lm(log(RFU) ~ days))[["days"]], .groups = "drop") %>%
    mutate(n_pts = n)
  
  slopes.all <- rbind(slopes.all, slopes.n)
}


exp.window <- slopes.all %>% # Pick the window with the steepest slope
  group_by(well_plate) %>%
  slice_max(order_by = slope, n = 1, with_ties = FALSE) %>%
  ungroup()

df.jb.rfu.exp <- df.jb.rfu %>% # Filter to exponential phase and add treatment info
  group_by(well_plate) %>%
  mutate(time_point = row_number()) %>%
  ungroup() %>%
  left_join(exp.window, by = "well_plate") %>%
  filter(time_point <= n_pts) %>%
  filter(well_plate != "D11_30") %>%          # JB outlier removal
  left_join(df.trt.jb, by = "population") %>%
  filter(!is.na(ancestor_id)) %>%
  group_by(population) %>%
  mutate(N0_mean = mean(N0)) %>%              # recompute after filter
  ungroup()

jb.dir.params <- data.frame()

for (pop in unique(df.jb.rfu.exp$population)) {
  
  df.i <- df.jb.rfu.exp %>% 
    filter(population == pop)
  
  mod <- nlsLM(RFU ~ N0_mean * exp((umax * (phosphate_concentration / (ks + phosphate_concentration))) * days),
               data    = df.i,
               start   = c(umax = 1, ks = 1),
               lower   = c(umax = 0, ks = 0),
               upper   = c(umax = 3, ks = 100),
               control = nls.control(maxiter = 1024, minFactor = 1/204800000))
               
  jb.dir.params <- rbind(jb.dir.params, data.frame(population = pop,
                                                   umax = coef(mod)["umax"],
                                                   ks = coef(mod)["ks"]))
  
}

comp <- jb.dir.params %>% # join data for comparison purposes
  left_join(df.jb.pt %>% select(population, ks.jb = ks, umax.jb = umax), 
            by = "population")

cor(comp$umax, comp$umax.jb)^2  # 0.9999978
cor(comp$ks, comp$ks.jb)^2      # 0.9999999

p.1A <- ggplot(comp, aes(x = umax.jb, y = umax)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "JB umax", y = "Replicated umax", 
       title = expression(paste("A — umax,  ", R^2, " = ", 0.9999978)))

p.1A

p.1B <- ggplot(comp, aes(x = ks.jb, y = ks)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "JB ks", y = "Replicated ks", 
       title = expression(paste("B — ks,  ", R^2, " = ", 0.9999999)))

p.1B

plot_grid(p.1A, p.1B)

jb.dir.params$rstar <- jb.dir.params$ks * m / (jb.dir.params$umax - m) # m 0.56

df.jb.pt$rstar.1  <- df.jb.pt$ks * 0.1  / (df.jb.pt$umax - 0.1) # recalculate at 0.1 m
df.jb.pt$rstar.56 <- df.jb.pt$ks * 0.56 / (df.jb.pt$umax - 0.56) # and at 0.56

ggplot(df.jb.pt) + # which one matches figure 4? should be 0.56
  geom_point(aes(x = rstar.1, y = umax, color = "m = 0.1")) +
  geom_point(aes(x = rstar.56, y = umax, color = "m = 0.56")) +
  labs(x = "P* (µM P)", y = expression(mu[max] ~ (day^-1)), color = "") +
  theme_cowplot() # yep its 0.56, so now we can compare our results. 

comp <- jb.dir.params %>%
  left_join(df.jb.pt %>% select(population, rstar.jb = rstar.56), by = "population")

round(cor(comp$rstar, comp$rstar.jb)^2, 6) # 0.999999

ggplot(comp, aes(x = rstar.jb, y = rstar)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "JB P* (µM P)", y = "Replicated P* (µM P)",
       title = expression(paste("P* comparison  ", R^2, " = 0.999999"))) +
  theme_cowplot()

df.trt.jb$treatment[is.na(df.trt.jb$treatment)] <- "none"

jb.dir.params <- jb.dir.params %>%
  left_join(df.trt.jb %>% select(population, ancestor_id, treatment), by = "population")

anc.means <- jb.dir.params %>%  # Ancestral mean P* per lineage
  filter(treatment == "none") %>%
  group_by(ancestor_id) %>%
  summarise(anc.mean = mean(rstar))

jb.dir.params <- jb.dir.params %>% # Change from ancestor
  left_join(anc.means, by = "ancestor_id") %>%
  mutate(change.pstar = rstar - anc.mean)

jb.dir.params <- jb.dir.params %>%
  mutate(trt.plot = ifelse(treatment == "none", "A", treatment)) %>%
  mutate(trt.plot = factor(trt.plot, levels = c("A", "C", "L", "N", "P", "B", "S", "BS")))

cols <- c("A"  = "black",
          "C"  = "#666666",
          "L"  = "#2166AC",
          "N"  = "#74ADD1",
          "P"  = "#4DAC8A",
          "B"  = "#F4A520",
          "S"  = "#F08080",
          "BS" = "#8B1A4A")

ggplot(jb.dir.params, aes(x = trt.plot, y = change.pstar, color = trt.plot)) +
  geom_hline(yintercept = 0) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 5, shape = 16) +
  scale_color_manual(values = cols, guide = "none") +
  coord_cartesian(ylim = c(-1.0, 1.5)) +
  labs(x = "", y = expression("change in P* " ~ (mu*M ~ P))) +
  theme_cowplot()

###### Save the files we need for a Shiny app comparison ######

write.csv(jb.dir.params %>% select(population, umax, ks),
          "JB-v-JL/shiny_mu_compare/data/jb_dir_params.csv", row.names = FALSE)

well_info <- df.jb.rfu %>%
  select(well_plate, phosphate_concentration, population) %>%
  distinct() %>%
  mutate(well = sub("_.*", "", well_plate),
         population = as.character(population))

jb.slopes <- exp.window %>%
  left_join(well_info, by = "well_plate") %>%
  rename(phos = phosphate_concentration) %>%
  group_by(population) %>%
  mutate(rep = as.integer(factor(well, levels = sort(unique(well))))) %>%
  ungroup()

write.csv(jb.slopes, "JB-v-JL/shiny_mu_compare/data/jb_exp_slopes.csv", row.names = FALSE)

# Section 2: first estimating mu (indirect approach) ---------------------------------------------------------------

#OK we will use my µ estimates for P

jl.ind.params <- data.frame()

for (pop in unique(df.jl.mu$population)) {
  
  df.i <- df.jl.mu %>%
    filter(population == pop,
           `µ` > 0)
  
  mod <- nlsLM(`µ` ~ umax * (phos / (ks + phos)),
               data    = df.i,
               start   = c(umax = 1, ks = 1),
               lower   = c(umax = 0, ks = 0),
               upper   = c(umax = 3, ks = 100),
               control = nls.control(maxiter = 1024, minFactor = 1/204800000))
  
  jl.ind.params <- rbind(jl.ind.params, data.frame(population = pop,
                                                   umax = coef(mod)["umax"],
                                                   ks   = coef(mod)["ks"]))
}

comp <- jl.ind.params %>% # join data for comparison purposes
  left_join(df.jb.pt %>% select(population, ks.jb = ks, umax.jb = umax), 
            by = "population")

cor(comp$umax, comp$umax.jb)^2  # 0.06965971
cor(comp$ks, comp$ks.jb)^2      # 0.1733072

p.2A <- ggplot(comp, aes(x = umax.jb, y = umax)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "JB umax", y = "Indirect umax", 
       title = expression(paste("A — umax,  ", R^2, " = ", 0.06965971)))

p.2A

p.2B <- ggplot(comp, aes(x = ks.jb, y = ks)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "JB ks", y = "Indirect ks", 
       title = expression(paste("B — ks,  ", R^2, " = ", 0.1733072)))

p.2B

plot_grid(p.2A, p.2B)

jl.ind.params$rstar <- jl.ind.params$ks * m / (jl.ind.params$umax - m)

comp <- jl.ind.params %>%
  left_join(df.jb.pt %>% select(population, rstar.jb = rstar.56), by = "population")

round(cor(comp$rstar, comp$rstar.jb)^2, 6) # 0.367466

ggplot(comp, aes(x = rstar.jb, y = rstar)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "JB P* (µM P)", y = "Inidrect P* (µM P)",
       title = expression(paste("Indirect P* comparison  ", R^2, " = 0.367466"))) +
  theme_cowplot()

jl.ind.params <- jl.ind.params %>%
  left_join(df.trt.jb %>% select(population, ancestor_id, treatment), by = "population")

anc.means <- jl.ind.params %>%  # Ancestral mean P* per lineage
  filter(treatment == "none") %>%
  group_by(ancestor_id) %>%
  summarise(anc.mean = mean(rstar))

jl.ind.params <- jl.ind.params %>% # Change from ancestor
  left_join(anc.means, by = "ancestor_id") %>%
  mutate(change.pstar = rstar - anc.mean)

jl.ind.params <- jl.ind.params %>%
  mutate(trt.plot = ifelse(treatment == "none", "A", treatment)) %>%
  mutate(trt.plot = factor(trt.plot, levels = c("A", "C", "L", "N", "P", "B", "S", "BS")))

ggplot(jl.ind.params, aes(x = trt.plot, y = change.pstar, color = trt.plot)) +
  geom_hline(yintercept = 0) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 5, shape = 16) +
  scale_color_manual(values = cols, guide = "none") +
  coord_cartesian(ylim = c(-1.0, 1.5)) +
  labs(x = "", y = expression("change in P* " ~ (mu*M ~ P))) +
  theme_cowplot()

# Testing curve plotting --------------------------------------------------

head(df.jb.rfu)

rfu.test <- df.jb.rfu %>%
  filter(population == 1,
         grepl("^B09_", well_plate))

rfu.test <- rfu.test %>% 
  filter(phosphate_concentration == 2)

rfu.test <- rfu.test %>%
  left_join(exp.window %>% select(well_plate, n_pts), by = "well_plate") %>%
  arrange(days) %>%
  mutate(time_point = row_number(),
         in_exp = time_point <= n_pts)

params.1 <- jb.dir.params %>% filter(population == "1")
mu.pred   <- params.1$umax * 2 / (params.1$ks + 2)
rfu.test$RFU_pred <- rfu.test$N0_mean * exp(mu.pred * rfu.test$days)

ggplot(rfu.test, aes(x = days, y = RFU)) +
  geom_point(aes(color = in_exp), size = 3) +
  geom_line(data = rfu.test %>% filter(in_exp), aes(y = RFU_pred)) +
  theme_cowplot()
