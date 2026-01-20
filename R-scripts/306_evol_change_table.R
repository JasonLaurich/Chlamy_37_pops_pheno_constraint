# Jason R Laurich

# January 18th, 2025

# The final analysis — how much each trait shifted from its ancestor coupled with modelling.

# Packages & functions ----------------------------------------------------

# Packages & Functions ----------------------------------------------------

library(tidyverse)
library(emmeans)
library(lme4)
library(performance)

# Load & examine the data -------------------------------------------------

df <- read.csv("data-processed/304_summary_table_final.csv") # Summary file
head(df)

# (a) Evolutionary changes in traits -----------------------

# We need to calculate changes in values

df.sum <- df %>%
  
  group_by(Anc) %>%
  
  mutate(
    I.comp = I.comp,                                     # Comp for I
    I.anc = mean(I.comp[match(Anc, Pop.fac) ]),          # Extract /I* of the ancestor
    d.I.comp = I.comp - I.anc,                           # Calculate difference for the mean.
    I.µ = I.µ.max,                                       # µmax
    I.µ.anc = mean(I.µ.max[match(Anc, Pop.fac) ]),       # Extract ancestral µmax
    d.I.µ = I.µ - I.µ.anc,                               # Calculate the change in µmax
  
    N.comp = N.comp,                                     # same for N
    N.anc = mean(N.comp[match(Anc, Pop.fac) ]),                
    d.N.comp = N.comp - N.anc,                           
    N.µ = N.µ.max,                                       
    N.µ.anc = mean(N.µ.max[match(Anc, Pop.fac) ]),             
    d.N.µ = N.µ - N.µ.anc,                               
    
    P.comp = P.comp,                                     # and P
    P.anc = mean(P.comp[match(Anc, Pop.fac) ]),                
    d.P.comp = P.comp - P.anc,                           
    P.µ = P.µ.max,                                       
    P.µ.anc = mean(P.µ.max[match(Anc, Pop.fac) ]),             
    d.P.µ = P.µ - P.µ.anc,                
    
    S.c = S.c,                                     # and salt tolerance (c)
    S.c.anc = mean(S.c[match(Anc, Pop.fac) ]),                
    d.S.c = S.c - S.c.anc,                           
    S.µ = S.µ.max,                                       
    S.µ.anc = mean(S.µ.max[match(Anc, Pop.fac) ]),             
    d.S.µ = S.µ - S.µ.anc,        
    
    T.br = T.br,                                     # and temperature
    T.br.anc = mean(T.br[match(Anc, Pop.fac) ]),                
    d.T.br = T.br - T.br.anc,                           
    T.µ = T.µ.max,                                       
    T.µ.anc = mean(T.µ.max[match(Anc, Pop.fac) ]),             
    d.T.µ = T.µ - T.µ.anc,        
    
    mean.N.µg.l = mean.N.µg.l,                       # and nitrogen content
    mean.N.µg.l.anc = mean.N.µg.l[match(Anc, Pop.fac) ],                
    d.mean.N.µg.l = mean.N.µg.l - mean.N.µg.l.anc,     
    
    mean.P.µg.l = mean.P.µg.l,                       # and phosphorous content
    mean.P.µg.l.anc = mean(mean.P.µg.l[match(Anc, Pop.fac) ]),                
    d.mean.P.µg.l = mean.P.µg.l - mean.P.µg.l.anc, 
    
    pig.PC = pig.PC,                                 # and pigmentation
    pig.PC.anc = mean(pig.PC[match(Anc, Pop.fac) ]),                
    d.pig.PC = pig.PC - pig.PC.anc, 
    
    bio.vol = bio.vol,                               # and biovolume
    bio.vol.anc = mean(bio.vol[match(Anc, Pop.fac) ]),                
    d.bio.vol = bio.vol - bio.vol.anc
    
  ) %>%
  
  ungroup() %>% 
  
  select(Pop.fac, Anc, Evol, I.comp, d.I.comp, I.µ, d.I.µ, 
         N.comp, d.N.comp, N.µ, d.N.µ,
         P.comp, d.P.comp, P.µ, d.P.µ,
         S.c, d.S.c, S.µ, d.S.µ,
         T.br, d.T.br, T.µ, d.T.µ,
         mean.N.µg.l, d.mean.N.µg.l,
         mean.P.µg.l, d.mean.P.µg.l,
         pig.PC, d.pig.PC,
         bio.vol, d.bio.vol)

###### Modelling ######

df.sum <- df.sum %>% filter(Evol != "none")

df.sum$Evol <- as.factor(df.sum$Evol)
df.sum$Anc <- as.factor(df.sum$Anc)

# Light

mod.I.comp.null <- lmer(d.I.comp ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.I.comp.full <- lmer(d.I.comp ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.I.comp.null, mod.I.comp.full)
performance::r2(mod.I.comp.full)

emm.I.comp <- emmeans(mod.I.comp.full, ~ Evol)
test(emm.I.comp)

mod.I.µ.null <- lmer(d.I.µ ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.I.µ.full <- lmer(d.I.µ ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.I.µ.null, mod.I.µ.full)
performance::r2(mod.I.µ.full)

emm.I.µ <- emmeans(mod.I.µ.full, ~ Evol)
test(emm.I.µ)

# Nitrogen

mod.N.comp.null <- lmer(d.N.comp ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.N.comp.full <- lmer(d.N.comp ~ Evol + (1 | Anc), data = df.sum, REML = FALSE) # Ancestry has no effect!

mod.N.comp.full <- lm(d.N.comp ~ Evol, data = df.sum)
mod.N.comp.null <- lm(d.N.comp ~ 1, data = df.sum)

anova(mod.N.comp.null, mod.N.comp.full)

performance::r2(mod.N.comp.full)

emm.N.comp <- emmeans(mod.N.comp.full, ~ Evol)
test(emm.N.comp)

mod.N.µ.null <- lmer(d.N.µ ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.N.µ.full <- lmer(d.N.µ ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.N.µ.null, mod.N.µ.full)
performance::r2(mod.N.µ.full)

emm.N.µ <- emmeans(mod.N.µ.full, ~ Evol)
test(emm.N.µ)

# Phosphorous

mod.P.comp.null <- lmer(d.P.comp ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.P.comp.full <- lmer(d.P.comp ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.P.comp.null, mod.P.comp.full)
performance::r2(mod.P.comp.full)

emm.P.comp <- emmeans(mod.P.comp.full, ~ Evol)
test(emm.P.comp)

mod.P.µ.null <- lmer(d.P.µ ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.P.µ.full <- lmer(d.P.µ ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.P.µ.null, mod.P.µ.full)
performance::r2(mod.P.µ.full)

emm.P.µ <- emmeans(mod.P.µ.full, ~ Evol)
test(emm.P.µ)

# Salt

mod.S.c.null <- lmer(d.S.c ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.S.c.full <- lmer(d.S.c ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.S.c.null, mod.S.c.full)
performance::r2(mod.S.c.full)

emm.S.c <- emmeans(mod.S.c.full, ~ Evol)
test(emm.S.c)

mod.S.µ.null <- lmer(d.S.µ ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.S.µ.full <- lmer(d.S.µ ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.S.µ.null, mod.S.µ.full)
performance::r2(mod.S.µ.full)

emm.S.µ <- emmeans(mod.S.µ.full, ~ Evol)
test(emm.S.µ)

# Temperature

mod.T.br.null <- lmer(d.T.br ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.T.br.full <- lmer(d.T.br ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.T.br.null, mod.T.br.full)
performance::r2(mod.T.br.full)

emm.T.br <- emmeans(mod.T.br.full, ~ Evol)
test(emm.T.br)

mod.T.µ.null <- lmer(d.T.µ ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.T.µ.full <- lmer(d.T.µ ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.T.µ.null, mod.T.µ.full)
performance::r2(mod.T.µ.full)

emm.T.µ <- emmeans(mod.T.µ.full, ~ Evol)
test(emm.T.µ)

# Stoichiometry

mod.mean.N.µg.l.null <- lmer(d.mean.N.µg.l ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.mean.N.µg.l.full <- lmer(d.mean.N.µg.l ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.mean.N.µg.l.null, mod.mean.N.µg.l.full)
performance::r2(mod.mean.N.µg.l.full)

emm.mean.N.µg.l <- emmeans(mod.mean.N.µg.l.full, ~ Evol)
test(emm.mean.N.µg.l)

mod.mean.P.µg.l.null <- lmer(d.mean.P.µg.l ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.mean.P.µg.l.full <- lmer(d.mean.P.µg.l ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.mean.P.µg.l.null, mod.mean.P.µg.l.full)
performance::r2(mod.mean.P.µg.l.full)

emm.mean.P.µg.l <- emmeans(mod.mean.P.µg.l.full, ~ Evol)
test(emm.mean.P.µg.l)

# Pigmentation

mod.pig.PC.null <- lmer(d.pig.PC ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.pig.PC.full <- lmer(d.pig.PC ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.pig.PC.null, mod.pig.PC.full)
performance::r2(mod.pig.PC.full)

emm.pig.PC <- emmeans(mod.pig.PC.full, ~ Evol)
test(emm.pig.PC)

# Pigmentation

mod.bio.vol.null <- lmer(d.bio.vol ~ (1 | Anc), data = df.sum, REML = FALSE)
mod.bio.vol.full <- lmer(d.bio.vol ~ Evol + (1 | Anc), data = df.sum, REML = FALSE)

anova(mod.bio.vol.null, mod.bio.vol.full)
performance::r2(mod.bio.vol.full)

emm.bio.vol <- emmeans(mod.bio.vol.full, ~ Evol)
test(emm.bio.vol)
