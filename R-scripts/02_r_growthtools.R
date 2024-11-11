# Jason R Laurich
# November 11, 2024

# This script will fit growth rates to Chlamydomonas reinhardtii at different temperatures,
# using the growthTools package (https://github.com/ctkremer/growthTools/tree/master)

############# Packages #####################

#remotes::install_github("ctkremer/mleTools")
#remotes::install_github("ctkremer/growthTools")

library(growthTools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(png)
library(cowplot)

############# Upload and examine data

df <- read.csv("data-processed/chlamee-acute-rfu-time.csv")
head(df) #RFU is density, days is time, temperature is temperature
str(df)

df<-df[,-c(1,2,4,5,6,8,13)]

df$pop.fac <- as.factor(df$population)
df$pop.num <- as.numeric(df$pop.fac)
df<-subset(df,df$temperature!=20)

# df$days <- df$days + 0.001 # Can't have 0s
df$logRFU <- log(df$RFU + 0.001)

levels(df$pop.fac) # I don't recognize the cc1629 group, but we'll keep it for now

df.rep <- subset(df, df$plate_type == "repeat") # We're going to work only with repeat data (this is most of the well_plates) which I'll use as replicates. 
df.rep$well.ID<-as.factor(df.rep$well_plate)

mat.rep <- split(df.rep, df.rep$pop.num)  # Each element is a data frame for one population in df.rep

############# Data exploration #############

# We'll eventually loop this through all of my temperatures and populations, but for now:

df.it <- subset(mat.rep[[3]], temperature==28)
df.it <- droplevels(df.it)

df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == 1)

res<-get.growth.rate(df.it.wl$days, df.it.wl$logRFU, plot.best.Q = T, id = 'Population 11, 28C, Well B07_46')

# Extract some summary statistics
res$best.model # gr.sat, no lag period.
res$best.slope
res$best.model.rsqr
res$best.se

# OK, let's think about modelling this for all temperatures on populations.
# We'll want to record.
  # 1. best.model - are they all the same? Except for T = 40?
  # 2. best.slope and best.se, obviously!
  # 3. We are not going to want to plot anything yet, we'll likely do this for a couple of populations later.

################# Fit models and extract r from all populations and treatments ################

# I have a strong suspicion that this approach won't fit our 40 C data, but I'll try

tmp <- as.vector(as.numeric(as.character(unique(df$temperature)))) # for looping through temperatures
ord.tmp<- sort(tmp)

df.r.gt <- data.frame( # dataframe for storage
  population = character(),
  population.number = numeric(),
  temperature = numeric(),
  well.ID = character(),
  r.gt = numeric(),
  r.best.model = character()
)

# Use the same looping structure as the 01_r_estimation file.

for (i in 1:38){ #population
  for (t in ord.tmp){ # temperature
    
    df.it <- subset(mat.rep[[i]], temperature==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      res<-get.growth.rate(df.it.wl$days, df.it.wl$logRFU, plot.best.Q = F)
      
      # Add data to our summary table
      df.r.gt <- rbind(df.r.gt, data.frame(
        population = df.it.wl$pop.fac[1],             # Population as factor
        population.number = df.it.wl$pop.num[1],      # Numeric population number
        temperature = df.it.wl$temperature[1],        # Temperature
        well.ID = df.it.wl$well.ID[1],                # Well ID (assuming one well ID per subset)
        r.gt =  res$best.slope,                        # The calculated r value (slope) from the best model
        r.best.model = res$best.model                 # The model that fit the data best
      ))
      
    }
  }
}

# OK, that actually worked for 40 C even, and the data looks good. Let's save the file and generate some summary plots to discuss.

write.csv(df.r.gt, "data-processed/03_r_growthtools.csv") # Save r estimates and other info

# randomly select 3 populations
pops<-names(mat.rep)
ran <- sample(pops, 3, replace = F)

plot.list <- list()

for (i in ran){ #population
  
  for (t in ord.tmp){ # temperature
    
    df.it <- subset(mat.rep[[i]], temperature==t)
    df.it <- droplevels(df.it) # Drop unused levels to isolate well replicate IDs at given t
    
    for (w in 1:length(levels(df.it$well.ID))){
      
      df.it.wl <- subset(df.it, as.numeric(df.it$well.ID) == w)
      
      res<-get.growth.rate(df.it.wl$days, df.it.wl$logRFU, plot.best.Q = T, id = paste("Pop:", df.it.wl$pop.fac[1], "Temp:", t, "Well:", w))
      
      p <- recordPlot()
      plot.list[[paste0("Pop", df.it.wl$pop.fac[1], "_temp", t, "_well", w)]] <- p
    
    }
  }
}

output_dir <- "figures"

# Save the first 10 recorded plots from plot.list as PNG files
for (i in 1:72) {
  # Specify the file path
  file_path <- file.path(output_dir, paste0("plot_", i, ".png"))
  
  # Open a PNG device
  png(filename = file_path, width = 800, height = 600)
  
  # Replay the recorded plot to the device
  replayPlot(plot.list[[i]])
  
  # Close the device
  dev.off()
}

# Load saved PNGs into a list of rasterGrob objects for grid arrangement
image_grobs_1 <- list()
for (i in 1:24) {
  # Load the image
  img_path <- file.path(output_dir, paste0("plot_", i, ".png"))
  img <- readPNG(img_path)
  
  # Convert to a rasterGrob
  image_grobs_1[[i]] <- rasterGrob(img, interpolate = TRUE)
}

image_grobs_2 <- list()
for (i in 25:48) {
  # Load the image
  img_path <- file.path(output_dir, paste0("plot_", i, ".png"))
  img <- readPNG(img_path)
  
  # Convert to a rasterGrob
  image_grobs_2[[i]] <- rasterGrob(img, interpolate = TRUE)
}

image_grobs_3 <- list()
for (i in 49:72) {
  # Load the image
  img_path <- file.path(output_dir, paste0("plot_", i, ".png"))
  img <- readPNG(img_path)
  
  # Convert to a rasterGrob
  image_grobs_3[[i]] <- rasterGrob(img, interpolate = TRUE)
}

# Arrange the images in a grid
grid.arrange(grobs = image_grobs_1, ncol = 4, nrow = 6, top = "Population 5")

grid.arrange(grobs = image_grobs_2, ncol = 4, nrow = 6, top = "Population 27")

grid.arrange(grobs = image_grobs_2, ncol = 4, nrow = 6, top = "Population 30")
