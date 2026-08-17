# Pareto fronts reveal constraints on the evolution of niche-determining traits in phytoplankton

Code and data repository for Laurich et al., 2026 (PNAS submission, manuscript ID - 2026-07212R)

## Project overview

This repository contains the data, processing code, analytical code, figure generation scripts, and figures for the manuscript "Widespread Pareto fronts
constrain the evolution of ecological limits in phytoplankton". This project investigates whether multivariate trade-offs among ecological traits generate 
Pareto fronts that constrain the evolution of niche limits in phytoplankton.

The study combines experimental evolution of *Chlamydomonas reinhardtii* populations across multiple environmental gradients (including temperature, nutrient
availability, light limitation, and salinity) with comparative analyses of trait data compiled from published phytoplankton datasets. Growth rates, niche 
breadths, and competitive ability were estimated for evolved and ancestral populations of *C. reinhardtii* (and for other species) and used to test for the 
presence of Pareto fronts using several complementary statistical approaches.

The repository includes scripts for processing raw experimental data, estimating growth rates, fitting thermal performance curves, Monod curves, and salt
tolerance curves, fitting and testing the significance of Pareto front models and quantile regressions, and generating the figures used in the manuscript.

All figures and analyses included in the manuscript can be reproduced from the scripts contained in this repository. Note that we do not include here objects
produced using Bayesian analysis in R2jags (TPCs, Monod curves, and salt tolerance curves) themselves due to memory constraints. However, the scripts necessary
to create them and summary tables containing the relevant model outputs (required for analysis downstream) are included. 

## Repository structure

repository/
├── README.md              README
├── .gitignore
├── lactin.txt             
├── monod.light.txt
├── monod.txt
├── salt.tolerance.txt     .txt files for running R2jags models
│
├── scripts-R/             Analytical scripts
│   └── scripts 00-21         Information on purpose of scripts, required files in metadata at the top. Master script (00) that runs the full analysis pipeline is included herein.
│
├── processed-data/        Data files (sourced from other publications, and processed data)
│   ├── README.md             Detailed description of files
│   └── files 01-79
│
├── Shiny-apps/            Shiny apps for confirming model fits to raw data (Chlamydomonas exp evolution and synthesis data)
│   └── folders 01-05
│
├── R2jags-models/         Folder for storing R2 jags models (not uploaded to github due to size constraints)
│
├── figures-main/          Main text figures
│   └── figures 01-05
│
├── figures-supp/          Supporting text figures
│   └── figures 01-08
│
└── figures-misc/          Additional exploratory figures
    └── figure 01

## Reproducing the analysis

To reproduce the full analysis:

1. Clone the repository
2. Install required packages
3. Run the master script. NOTE: currently the scripts that require heavy model-fitting are commented out. Remove the '#'s to run the whole pipeline.
4. OR run the scripts in scripts-R in numerical order (01-21)

## Software requirements

R version 4.4.1+

Required packages include:
bayestestR
boot
car
cowplot
Deriv
emmeans
lme4
minpack.lm
MuMIn
nls.multstart
patchwork
performance
pracma
rTPC
R2jags
scam
smatr
sp
tidyverse
vegan
ggforce

## Author contact information
Jason R. Laurich
University of Guelph
Department of Integrative Biology
jlaurich@uoguelph.ca
jason@laurich.ca