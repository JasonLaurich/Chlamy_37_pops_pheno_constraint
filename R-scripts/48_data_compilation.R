# Jason R Laurich

# October 7th, 2025

# We are going to combine all of the data we need at the individual level to generate a massive summary data frame.
# I'm also going to explore collapsing variation in pigments into a single axis using PCA.

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(vegan)  # For PCA

# Load & examine the data -------------------------------------------------

# Temperature

df.t <- read.csv("data-processed/303a_TPC_estimates_final.csv") # This file has the TPC data for each replicate.
head(df.t)

df.t <- df.t[,c(2:4,8,10)] # We'll work with the analytically determined parameters (Deriv package)
head(df.t)

df.t <-df.t %>% # Need to convert 38 to 37 (from the inclusion of cc1629 in this dataset only)
  mutate(Pop.num = replace(Pop.num, Pop.num == 38, 37L),
         unique.id = case_when(
           unique.id == "38.B04" ~ "37.B04",
           unique.id == "38.C05" ~ "37.C05",
           unique.id == "38.E05" ~ "37.E05",
           unique.id == "38.F07" ~ "37.F07",
           TRUE ~ unique.id
         ))

# Light

df.l <- read.csv("data-processed/300a_Monod_light_estimates_new.csv") # This file has the light Monod curve data for each replicate.
head(df.l)

df.l <- df.l %>% 
  mutate(R.mth = R.mth*2.5) %>%    # Convert to µmol/m^2/s PAR
  mutate(I.comp = 1/R.mth)

df.l <- df.l[,c(2:4,6,9)] # We'll work with the analytically determined parameters (Deriv package)
head(df.l)

# Nitrogen

df.n <- read.csv("data-processed/301a_Monod_nit_estimates_new.csv") # This file has the TPC data for each replicate.
head(df.n)

df.n <- df.n %>% 
  mutate(N.comp = 1/R.mth)

df.n <- df.n[,c(2:4,6,9)] # We'll work with the analytically determined parameters (Deriv package)
head(df.n)

# Phosphorous

df.p <- read.csv("data-processed/302a_Monod_phosphorous_estimates_new.csv") # This file has the TPC data for each replicate.
head(df.p)

df.p <- df.p %>% 
  mutate(P.comp = 1/R.mth)

df.p <- df.p[,c(2:4,6,9)] # We'll work with the analytically determined parameters (Deriv package)
head(df.p)

# Combine T, L, N, and P --------------------------------------------------


df.summ <- df.t %>% 
  transmute(Pop.fac = Pop.fac,
            Pop.num = Pop.num,
            unique.id = unique.id,
            T.µ.max = r.max,
            T.br = T.br) %>% 
  
  left_join(df.l %>% 
              transmute(I.µ.max = r.max,
                        I.comp = I.comp,
                        unique.id = unique.id),
            by = c("unique.id" = "unique.id")) %>% 
  
  left_join(df.n %>% 
              transmute(N.µ.max = r.max,
                        N.comp = N.comp,
                        unique.id = unique.id),
            by = c("unique.id" = "unique.id")) %>% 
  
  left_join(df.p %>% 
              transmute(P.µ.max = r.max,
                        P.comp = P.comp,
                        unique.id = unique.id),
            by = c("unique.id" = "unique.id"))
  
head(df.summ)

# OK now we are going to add in other stuff.

# Salt

df.s <- read.csv("data-processed/216_salt_pops_estimates_final.csv") # This file has the TPC data for each replicate.
head(df.s)

# Stoichiometry data

df.stoich <- read.csv("data-processed/17_stoich_data.csv") # Stoichiometry data for N and P
head(df.stoich)

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


df.id <- read.csv("data-processed/18_id_mapping.csv") # Let's bring in an identity mapping file. 

df.id

df.stoich.sum <- df.stoich.sum %>% # Join this with the df.id matching schema
  left_join(df.id, by = c("Name" = "sample"))

head(df.stoich.sum) # Correct population assignment in the population column

# Pigment data

df.pig <- read.csv("data-processed/19_pigment_data.csv") # Pigment data

head(df.pig)

levels(as.factor(df.pig$HPLC.Nummer)) # These are already numbered 1-37, ie. the missing populations do not feature here. 

df.pig<- df.pig %>%
  mutate(population = df.id$population)

# PCA of pigment data — we know these data are tightly correlated. 
df.pca <- df.pig %>% select(chl.a, chl.b, luthein) # Prepare the data: selecting only the relevant columns
df.pca

pca.result <- prcomp(df.pca, center = TRUE, scale. = TRUE) # Perform PCA
pca.result

sdev <- pca.result$sdev  

sdev %>%
  { .^2 / sum(.^2) } %>%
  print() 

df.pca.res <- data.frame(pca.result$x)
df.pca.res <- df.pca.res %>%
  mutate(population = df.id$population)

PCA <- ggplot(df.pca.res, aes(x = PC1, y = PC2)) +  # PCA biplot visualization
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("PC1 (", round(pca.result$sdev[1]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pca.result$sdev[2]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       title = "PCA of pigmentation") 

PCA

# Finally, we need cell size data

df.size <- read.csv("data-processed/52_cell_size.csv") # Biovolume data

head(df.size)

# Actually finally, we need the treatment history data!

df.hist <- read.csv("data-processed/16_rfus_time_summary.csv") # This file has the history (ancestry, evolutionary treatment) for each well ID.
head(df.hist)

df.hist$anc <- as.factor(df.hist$ancestor_id) # Will port these into combined dataset, sorting by population. 
df.hist$evol <- as.factor(df.hist$treatment)
df.hist$Pop.fac <- as.factor(df.hist$population)

df.hist.agg <- df.hist %>% # Aggregate the history data for merging.
  group_by(Pop.fac) %>%
  summarise(
    anc = first(anc),
    evol = first(evol),
    .groups = "drop" 
  )

# Add in the data that appears only once for each population

df.summ <- df.summ %>% 
  left_join(df.s %>% 
              transmute(S.µ.max = r.max,
                        S.c = c.mod,
                        Pop.num = Pop.num),
            by = c("Pop.num" = "Pop.num")) %>% 
  
  left_join(df.stoich.sum %>%                                                           # Add in stoichiometry data
              transmute(Pop.fac = population,
                        mean.N.µg.l = mean.N.µg.l,
                        mean.P.µg.l = mean.P.µg.l), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(df.pig %>%                                                           # Add in pigment data
              transmute(Pop.fac = population,
                        chl.a = chl.a,
                        chl.b = chl.b,
                        luthein = luthein), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(df.pca.res %>%                                                           # Add in pigment  PCA data
              transmute(Pop.fac = population,
                        pig.PC = PC1), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(df.size %>%                                                               # Add in size
              transmute(Pop.fac = population,
                        bio.vol = bio.vol.mean), by = c("Pop.fac" = "Pop.fac")) %>% 
  
  left_join(df.hist.agg %>%                                                               # Add in history data
              transmute(Pop.fac = Pop.fac,
                        Anc = anc,
                        Evol = evol), by = c("Pop.fac" = "Pop.fac"))

write.csv(df.summ, "data-processed/304_summary_table_final.csv") # Save summary table!
