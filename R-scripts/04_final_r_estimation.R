# Jason R Laurich
# November 20, 2024

# This script will estimate r one final way (which is the way we intend to calculate estimates used for the fitting of TPCs)
# For all sub-40 C populations, I will fit logged linear models to our data, within a window beginning at our first time point.

# As the width of that window increases, I will calculate the slope of the line, and interpret a decrease in the slope as the
# end of the exponential growth phase.
# I will then fit an exponential growth curve to that period of data for each population and well replicate of Chlamydomonas.

# For 40 C, our data shows a brief growth spike that represents anomalous data (a completion of cell division already initiated)
# that is not biologically meaningful. 
# Thus, for this treatment, I will threshold data to after maximum RFU count, and fit the same exponential growth curve to it.

############# Packages ########################

library(nls.multstart)
library(tidyr)
library(cowplot)
library(ggplot2)

############# Upload and examine data #######################