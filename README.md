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

<pre>
Chlamy_37_pops_pheno_constraint/
├── <a href="README.md">README.md</a>
├── <a href="CITATION.cff">CITATION.cff</a>                Citation metadata
├── <a href="LICENSE">LICENSE</a>
├── <a href="Chlamy_37_pops_pheno_constraint.Rproj">Chlamy_37_pops_pheno_constraint.Rproj</a>   Open first (sets file paths)
├── <a href="renv.lock">renv.lock</a>                   Exact R package versions (see renv::restore())
├── <a href="renv">renv/</a>                       renv environment files
├── <a href="lactin.txt">lactin.txt</a>                  JAGS model files used by R2jags
├── <a href="monod.txt">monod.txt</a>
├── <a href="monod.light.txt">monod.light.txt</a>
├── <a href="salt.tolerance.txt">salt.tolerance.txt</a>
│
├── <a href="scripts-R">scripts-R/</a>                  Analysis scripts (00–21; see Script → output map)
│   └── <a href="scripts-R/00_master_script.R">00_master_script.R</a>      Runs the full pipeline
│
├── <a href="processed-data">processed-data/</a>             Data files (01–79): experimental data, model outputs,
│   │                           and data compiled from published studies
│   └── <a href="processed-data/README.md">README.md</a>               Description of each data file
│
├── <a href="Shiny-apps">Shiny-apps/</a>                 Shiny apps for checking model fits against raw data
│   ├── <a href="Shiny-apps/01_Chlamydomonas_exp_evol_model_quality_control">01_Chlamydomonas_exp_evol_model_quality_control/</a>
│   ├── <a href="Shiny-apps/02_synthesis_light_monod_quality_control">02_synthesis_light_monod_quality_control/</a>
│   ├── <a href="Shiny-apps/03_synthesis_nitrogen_monod_quality_control">03_synthesis_nitrogen_monod_quality_control/</a>
│   ├── <a href="Shiny-apps/04_synthesis_phosphorous_monod_quality_control">04_synthesis_phosphorous_monod_quality_control/</a>
│   └── <a href="Shiny-apps/05_synthesis_TPC_quality_control">05_synthesis_TPC_quality_control/</a>
│
├── <a href="figures-main">figures-main/</a>               Main text figures (Figs. 1–5)
├── <a href="figures-supplemental">figures-supplemental/</a>       Supporting figures (Figs. S1–S8)
├── <a href="figures-misc">figures-misc/</a>               Additional figures from reviewer-requested analyses
│
└── R2jags-models/              Created locally when scripts 01–05 are run; R2jags model
                                objects are not included in the repository (file size)
</pre>

## Script → output map

All scripts are in [`scripts-R/`](scripts-R). [`00_master_script.R`](scripts-R/00_master_script.R) runs the full pipeline in order.

### Main text figures

| Figure | Script | Output |
|---|---|---|
| Fig. 1 | [`07_conceptual_figure.R`](scripts-R/07_conceptual_figure.R) | [`01_fig_1_conceptual_figure.jpeg`](figures-main/01_fig_1_conceptual_figure.jpeg) |
| Fig. 2 | [`08_intra-gradient_trade-offs.R`](scripts-R/08_intra-gradient_trade-offs.R) | [`02_fig_2_intra-gradient_tradeoffs.pdf`](figures-main/02_fig_2_intra-gradient_tradeoffs.pdf) |
| Fig. 3 | [`09_inter-gradient_trade-offs.R`](scripts-R/09_inter-gradient_trade-offs.R) | [`03_fig_3_inter-gradient_toffs.pdf`](figures-main/03_fig_3_inter-gradient_toffs.pdf) |
| Fig. 4 | [`18_synthesis_pareto.R`](scripts-R/18_synthesis_pareto.R) | [`04_fig_4_synthesis_comps.pdf`](figures-main/04_fig_4_synthesis_comps.pdf) |
| Fig. 5 | [`19_summary_fig.R`](scripts-R/19_summary_fig.R) | [`05_fig_5_summary.pdf`](figures-main/05_fig_5_summary.pdf) |

### Supporting figures and tables

| Figure/Table | Script | Output |
|---|---|---|
| Fig. S1 | [`20_model_fitting_schema.R`](scripts-R/20_model_fitting_schema.R) | [`01_fig_s1_model_fitting.pdf`](figures-supplemental/01_fig_s1_model_fitting.pdf) |
| Fig. S2 | [`11_sensitivity_analysis.R`](scripts-R/11_sensitivity_analysis.R) | [`02_fig_s2_sensitivity.pdf`](figures-supplemental/02_fig_s2_sensitivity.pdf) |
| Fig. S3 | [`08_intra-gradient_trade-offs.R`](scripts-R/08_intra-gradient_trade-offs.R) | [`03_fig_s3_evol_opt_test.pdf`](figures-supplemental/03_fig_s3_evol_opt_test.pdf) |
| Fig. S4 | [`09_inter-gradient_trade-offs.R`](scripts-R/09_inter-gradient_trade-offs.R) | [`04_fig_s4_dist_to_pf.pdf`](figures-supplemental/04_fig_s4_dist_to_pf.pdf) |
| Fig. S5 | [`13_PCA_RDA.R`](scripts-R/13_PCA_RDA.R) | [`05_fig_s5_RDAs.pdf`](figures-supplemental/05_fig_s5_RDAs.pdf) |
| Fig. S6 | [`13_PCA_RDA.R`](scripts-R/13_PCA_RDA.R) | [`06_fig_s6_PCA.pdf`](figures-supplemental/06_fig_s6_PCA.pdf) |
| Fig. S7 | [`10_secondary_trait_trade-offs.R`](scripts-R/10_secondary_trait_trade-offs.R) | [`07_fig_s7_extra_inter-gradient_toffs.pdf`](figures-supplemental/07_fig_s7_extra_inter-gradient_toffs.pdf) |
| Fig. S8 | [`18_synthesis_pareto.R`](scripts-R/18_synthesis_pareto.R) | [`08_fig_s8_synthesis_with_Laurich_2026_data.pdf`](figures-supplemental/08_fig_s8_synthesis_with_Laurich_2026_data.pdf) |
| Table S1 | [`11_sensitivity_analysis.R`](scripts-R/11_sensitivity_analysis.R) | [`29_sensitivity_grid_table.csv`](processed-data/29_sensitivity_grid_table.csv) |

### Data processing and model fitting (no figures)

| Scripts | Purpose | Key output |
|---|---|---|
| [`01`](scripts-R/01_light.R)–[`05`](scripts-R/05_temperature.R) | Estimate µ and fit Monod, salt tolerance and thermal performance curves (R2jags; computationally heavy) for *C. reinhardtii* | Summary tables in `processed-data/` (files 02–21) |
| [`06_data_compilation.R`](scripts-R/06_data_compilation.R) | Compile all *C. reinhardtii* traits | [`27_summary_table.csv`](processed-data/27_summary_table.csv) (input to scripts 08–13, 19–21) |
| [`17`](scripts-R/17_synthesis_temperature.R), then [`14`](scripts-R/14_synthesis_light.R)–[`16`](scripts-R/16_synthesis_phosphorus.R) | Fit TPCs, then Monod curves, to the phytoplankton synthesis data (run 17 first) | Files 37–67 in `processed-data/` |
| [`12_evol_effects_modelling.R`](scripts-R/12_evol_effects_modelling.R) | Mixed models of evolutionary trait shifts | Statistics reported in text |
| [`08b_figure2_log_intra-gradient_trade-offs.R`](scripts-R/08b_figure2_log_intra-gradient_trade-offs.R), [`21_reviewer3_data_scaling.R`](scripts-R/21_reviewer3_data_scaling.R) | Reviewer-requested checks on log-transformed and rescaled data | [`figures-misc/`](figures-misc) |

## Reproducing the analysis

To reproduce the full analysis:

1. Clone the repository
2. Open Chlamy_37_pops_pheno_constraint.Rproj. All file paths are relative to the project root
3. Install required packages (see below)
4. Run the master script. NOTE: currently the scripts that require heavy model-fitting are commented out. Remove the '#'s to run the whole pipeline.
5. OR run the scripts in scripts-R in numerical order (01-21)

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

# Citation

If you use this code or data, please cite both the paper and the repository:

**Paper:** Laurich, J. R., Narwani, A., & Bernhardt, J. R. (2026). Pareto fronts reveal constraints on the evolution of ecological niches. *Proceedings of the National Academy of Sciences*. (in press; DOI to be added)

**Code and data:** Laurich, J. R., Narwani, A., & Bernhardt, J. R. (2026). Repository for "Pareto fronts reveal constraints on the evolution of ecological niches" (v2.1). Zenodo. https://doi.org/10.5281/zenodo.19032401

## Author contact information
Jason R. Laurich
University of Guelph
Department of Integrative Biology
jlaurich@uoguelph.ca
jason@laurich.ca