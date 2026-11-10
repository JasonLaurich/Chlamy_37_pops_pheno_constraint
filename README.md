# Pareto fronts reveal constraints on the evolution of ecological niches

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19032401.svg)](https://doi.org/10.5281/zenodo.19032401)

Code and data for Laurich, Narwani & Bernhardt (2026), *Proceedings of the National Academy of Sciences* (in press).

## Project overview

This repository contains the data, processing code, analytical code, figure generation scripts, and figures for the manuscript "Pareto fronts reveal 
constraints on the evolution of ecological niches". This project investigates whether multivariate trade-offs among ecological traits generate 
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

- R version 4.4.1+
- JAGS 4.3.x (required by R2jags; install separately from https://sourceforge.net/projects/mcmc-jags/ before installing the R packages)

### R packages

Exact versions of all packages (including dependencies) are recorded in [`renv.lock`](renv.lock). To recreate the package environment:

    install.packages("renv")
    renv::restore()

Require R packages include:

| Package | Version | 
|---|---|
| bayestestR | 0.19.0 |
| boot | 1.3-32 |
| car | 3.1-5 |
| cowplot | 1.2.0 |
| Deriv | 4.3.5 |
| emmeans | 2.0.4 |
| ggforce | 0.5.0 |
| gridExtra | 2.3.1 |
| lme4 | 2.0-6 |
| minpack.lm | 1.2-4 |
| MuMIn | 1.48.19 |
| nls.multstart | 2.0.0 |
| patchwork | 1.3.2 |
| performance | 0.19.0 |
| pracma | 2.4.6 |
| quantreg | 6.1 |
| R2jags | 0.8-9 |
| rTPC | 1.1.0 |
| scales | 1.4.0 |
| scam | 1.2-23 |
| shiny | 1.14.0 |
| smatr | 3.5-2 |
| sp | 2.2-4 |
| svglite | 2.2.2 |
| tibble | 3.3.1 |
| tidyr | 1.3.2 |
| tidyverse | 2.0.0 |
| vegan | 2.7-6 |

## Author contact information
Jason R. Laurich
University of Guelph
Department of Integrative Biology
jlaurich@uoguelph.ca
jason@laurich.ca