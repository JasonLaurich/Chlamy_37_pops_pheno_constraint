# Jason R Laurich

# February 8th, 2026

# This file will assemble summary data for all C. reinhardtii populations from our experimental evolution project.
# We incorporate niche-determining traits from scripts 01-05.R as well as pigmentation, stoichiometry, and biovolume data.

# Inputs: in processed-data : 04_TPC_summary.csv, 08_light_monod_summary.csv, 12_nit_monod_summary.csv, 16_phos_monod_summary.csv,
  # 20_salt_tol_summary.csv, 22_stoich_data.csv, 23_id_mapping.csv, 24_pigment_data.csv, 25_cell_size.csv, 26_labels_summary.csv
# Outputs: in processed-data : 27_summary_table.csv

# Packages & functions ----------------------------------------------------

library(tidyverse)
library(vegan)  # For PCA

# Load & examine the data -------------------------------------------------

# Temperature

df.t <- read_csv("processed-data/04_TPC_summary.csv") # This file has the TPC data for each replicate.
head(df.t)

df.t <- df.t %>% 
  mutate(T.br = T.br.max - T.br.min) %>% 
  select(population, rep.ID, T.min, T.max, T.opt, r.max, T.br, T.br.min, T.br.max, a, b, tmax, d.t)

# Light

df.l <- read_csv("processed-data/08_light_monod_summary.csv") # This file has the light Monod curve data for each replicate.
head(df.l)

df.l <- df.l %>% 
  select(rep.ID, K.s.post, r.max.post, R.post) %>% 
  mutate(comp = 1/R.post)

# Nitrogen

df.n <- read_csv("processed-data/12_nit_monod_summary.csv") # This file has the nitrogen Monod curve data for each replicate.
head(df.n)

df.n <- df.n %>% 
  select(rep.ID, K.s.post, r.max.post, R.post) %>% 
  mutate(comp = 1/R.post)

# Phosphorous

df.p <- read_csv("processed-data/16_phos_monod_summary.csv") # This file has the phosphorous Monod curve data for each replicate.
head(df.p)

df.p <- df.p %>% 
  select(rep.ID, K.s.post, r.max.post, R.post) %>% 
  mutate(comp = 1/R.post)

# Salt 

df.s <- read_csv("processed-data/20_salt_tol_summary.csv") # This file has the salt tolerance curve data for each population. 
head(df.s)

df.s <- df.s %>% 
  select(population, r.max.post, c.post, b.post) 

df.s <- df.s %>%
  mutate(population = tolower(gsub(" ", "", population))) # Convert Anc x to ancx, so that it matches other data frames. 

# Combine T, L, N, and P --------------------------------------------------


df.summ <- df.t %>% 
  transmute(population = population,
            rep.ID = rep.ID,
            T.min = T.min,
            T.opt = T.opt,
            T.max = T.max,
            T.µ.max = r.max,
            T.br = T.br,
            T.br.min = T.br.min,
            T.br.max = T.br.max,
            T.a = a,
            T.b = b,
            T.tmax = tmax,
            T.d.t = d.t) %>% 
  
  left_join(df.l %>% 
              transmute(I.µ.max = r.max.post,
                        I.comp = comp,
                        I.K.s = K.s.post,
                        rep.ID = rep.ID),
            by = c("rep.ID" = "rep.ID")) %>% 
  
  left_join(df.n %>% 
              transmute(N.µ.max = r.max.post,
                        N.comp = comp,
                        N.K.s = K.s.post,
                        rep.ID = rep.ID),
            by = c("rep.ID" = "rep.ID")) %>% 
  
  left_join(df.p %>% 
              transmute(P.µ.max = r.max.post,
                        P.comp = comp,
                        P.K.s = K.s.post,
                        rep.ID = rep.ID),
            by = c("rep.ID" = "rep.ID"))

head(df.summ)

# OK now we are going to add in other data, where observations exist only at the population level. 

# Stoichiometry data

df.stoich <- read.csv("processed-data/22_stoich_data.csv") # Stoichiometry data for N and P
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

df.id <- read.csv("processed-data/23_id_mapping.csv") # Let's bring in an identity mapping file. 

df.id

df.stoich.sum <- df.stoich.sum %>% # Join this with the df.id matching schema
  left_join(df.id, by = c("Name" = "sample"))

head(df.stoich.sum) # Correct population assignment in the population column

# Pigment data

df.pig <- read.csv("processed-data/24_pigment_data.csv") # Pigment data

head(df.pig)

levels(as.factor(df.pig$HPLC.Nummer)) # These are already numbered 1-37, ie. the missing populations do not feature here.

df.id <- df.id %>%  # But we still need to convert to population labels that are consistent with other data sets.
  
  mutate(id.pig = match(sample, sort(unique(sample))))

df.pig <- df.pig %>% # Join this with the df.id matching schema
  left_join(df.id, by = c("HPLC.Nummer" = "id.pig"))

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

loadings <- as.data.frame(pca.result$rotation) # Extract PCA loadings (rotation matrix)

loadings$PC1 <- loadings$PC1 * max(abs(df.pca.res$PC1)) # Scale loadings to fit within the PCA plot (can adjust scaling factor)
loadings$PC2 <- loadings$PC2 * max(abs(df.pca.res$PC2))

loadings$variable <- rownames(loadings) # Add variable names for annotation

PCA <- ggplot(df.pca.res, aes(x = PC1, y = PC2)) +  # PCA biplot visualization
  geom_point(size = 3) +
  theme_classic() +
  labs(x = paste("PC1 (", round(pca.result$sdev[1]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       y = paste("PC2 (", round(pca.result$sdev[2]^2 / sum(pca.result$sdev^2) * 100, 2), "%)", sep = ""),
       title = "PCA of pigmentation")  +
  
  geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size= 1.1)

PCA

# Cell size (biovolume)

df.size <- read.csv("processed-data/25_cell_size.csv") # Biovolume data

head(df.size)

# Finally, we need the treatment history data!

df.hist <- read.csv("processed-data/26_labels_summary.csv") # This file has the history (ancestry, evolutionary treatment) for each well ID.
head(df.hist)

df.hist$anc <- as.factor(df.hist$ancestor_id) # Will port these into combined dataset, sorting by population. 
df.hist$evol <- as.factor(df.hist$treatment)

df.hist.agg <- df.hist %>% # Aggregate the history data for merging.
  group_by(population) %>%
  summarise(
    anc = first(anc),
    evol = first(evol),
    .groups = "drop" 
  ) %>% 
  filter(population != "cc1629")

# Add in the data that appears only once for each population

df.summ <- df.summ %>% 
  left_join(df.s %>% 
              transmute(S.µ.max = r.max.post,
                        S.c = c.post,
                        S.b = b.post,
                        population = population),
            by = c("population" = "population")) %>% 
  
  left_join(df.stoich.sum %>%                                                   # Add in stoichiometry data
              transmute(population = population,
                        mean.N.µg.l = mean.N.µg.l,
                        mean.P.µg.l = mean.P.µg.l), 
            by = c("population" = "population")) %>% 
  
  left_join(df.pig %>%                                                          # Add in pigment data
              transmute(population = population,
                        chl.a = chl.a,
                        chl.b = chl.b,
                        luthein = luthein), 
            by = c("population" = "population")) %>% 
  
  left_join(df.pca.res %>%                                                      # Add in pigment  PCA data
              transmute(population = population,
                        pig.PC = PC1), 
            by = c("population" = "population")) %>% 
  
  left_join(df.size %>%                                                         # Add in size
              transmute(population = population,
                        bio.vol = bio.vol.mean), 
            by = c("population" = "population")) %>% 
  
  left_join(df.hist.agg %>%                                                     # Add in history data
              transmute(population = population,
                        Anc = anc,
                        Evol = evol), 
            by = c("population" = "population"))

write.csv(df.summ, "processed-data/27_summary_table.csv") # Save summary table!
