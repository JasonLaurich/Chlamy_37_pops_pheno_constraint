# Processed data

This folder contains datasets used in the analyses for Laurich et al. (2026), "Pareto fronts reveal constraints on the evolution
of ecological niches".

The folder includes both raw datasets obtained from previously published studies and derived datasets generated during the analysis 
pipeline. Scripts in scripts-R/ read from and write to this directory.

## File naming convention

Processed data files are numbered according to the order in which they are generated or used in the analysis pipeline.

## File breakdown

### Raw data: *C. reinhardtii experimental evolution

**Column definitions:** each data file has a matching metadata file (e.g. `01_light_rfus_time_metadata.csv`) describing every column and its units.

| File | Description |
|---|---|
| [`01_light_rfus_time.csv`](01_light_rfus_time.csv) | Raw growth time series (RFU) across the light gradient |
| [`05_nitrogen_rfus_time.csv`](05_nitrogen_rfus_time.csv) | Raw growth time series (RFU) across the nitrogen gradient |
| [`09_phosphorus_rfus_time.csv`](09_phosphorus_rfus_time.csv) | Raw growth time series (RFU) across the phosphorus gradient |
| [`13_salt_rfus_time.csv`](13_salt_rfus_time.csv) | Raw growth time series (RFU) across the salt gradient |
| [`17_temperature_rfus_time.csv`](17_temperature_rfus_time.csv) | Raw growth time series (RFU) across the temperature gradient |
| [`22_stoich_data.csv`](22_stoich_data.csv) | Cellular nitrogen and phosphorus content |
| [`23_id_mapping.csv`](23_id_mapping.csv) | Maps population labels in file 22 to those used elsewhere |
| [`24_pigment_data.csv`](24_pigment_data.csv) | Pigment data (chlorophyll *a*, chlorophyll *b*, lutein) |
| [`25_cell_size.csv`](25_cell_size.csv) | Biovolume data |
| [`26_labels_summary.csv`](26_labels_summary.csv) | Treatment information for all populations |

Files 01, 05, 09, 13, 17, 22, and 23 are from Bernhardt et al., 2020:

J. R. Bernhardt, et al., The evolution of competitive ability for essential resources. Philos. Trans. R. Soc. B Biol. Sci. 375, 20190247 (2020).

### Raw data: published pyhtoplankton studies

| File | Source |
|---|---|
| [`30_Lewington-Pearce_2019_raw_data.csv`](30_Lewington-Pearce_2019_raw_data.csv) | Lewington-Pearce et al. (2019) |
| [`31_Edwards_2016_raw_data.csv`](31_Edwards_2016_raw_data.csv) | Edwards et al. (2016) |
| [`32_Levasseur2025_l_rawdata.csv`](32_Levasseur2025_l_rawdata.csv), [`33`](33_Levasseur2025_n_rawdata.csv), [`35`](35_Levasseur2025_p_rawdata.csv) | Levasseur et al. (2025): light, nitrogen and phosphorus |
| [`34_Bestion_2018_raw_data.csv`](34_Bestion_2018_raw_data.csv) | Bestion et al. (2018) |
| [`36_Thomas_2012_raw_data.csv`](36_Thomas_2012_raw_data.csv) | Thomas et al. (2012) |
| [`68_Narwani_2015_summary.csv`](68_Narwani_2015_summary.csv) | Narwani et al. (2015) |
| [`69_Edwards_2015_summary.csv`](69_Edwards_2015_summary.csv) | Edwards et al. (2015) |

#### Data sources

- E. Bestion, B. García-Carreras, C.-E. Schaum, S. Pawar, G. Yvon-Durocher, Metabolic traits predict the effects of warming on phytoplankton competition. Ecol. Lett. 21, 655–664 (2018).
- K. F. Edwards, C. A. Klausmeier, E. Litchman, Nutrient utilization traits of phytoplankton. Ecology 96, 2311–2311 (2015).
- K. F. Edwards, M. K. Thomas, C. A. Klausmeier, E. Litchman, Phytoplankton growth and the interaction of light and temperature: a synthesis at the species and community level. Limnol. Oceanogr. 61, 1232–1244 (2016).
- L. Lewington-Pearce, et al., Temperature-dependence of minimum resource requirements alters competitive hierarchies in phytoplankton. Oikos 128, 1194–1205 (2019).
- S. Levasseur, et al., Warm-loving species perform well under limiting resources: trait combinations for future climate. Glob. Change Biol. 31, e70554 (2025).
- A. Narwani, et al., Common ancestry is a poor predictor of competitive traits in freshwater green algae. PLOS ONE 10, e0137085 (2015).
- M. K. Thomas, C. T. Kremer, C. A. Klausmeier, E. Litchman, A global pattern of thermal adaptation in marine phytoplankton. Science 338, 1085–1088 (2012).

## Derived data

| Files | Description | Script |
|---|---|---|
| 02, 06, 10, 14, 18 | Growth rate (µ) estimates for *C. reinhardtii* across each gradient | [01](../scripts-R/01_light.R)–[05](../scripts-R/05_temperature.R) |
| 03–04, 07–08, 11–12, 15–16, 20–21 | Model parameter summaries and fitted curves (Monod, salt tolerance, TPC) | 01–05 |
| [`19_temperature_nls_rTPC_comparison_AICcs.csv`](19_temperature_nls_rTPC_comparison_AICcs.csv) | AICc comparison of candidate TPC models | [05](../scripts-R/05_temperature.R) |
| [`27_summary_table.csv`](27_summary_table.csv) | Population-level trait summary (main input to downstream analyses) | [06](../scripts-R/06_data_compilation.R) |
| [`28_sensitivity_pvalues.csv`](28_sensitivity_pvalues.csv), [`29_sensitivity_grid_table.csv`](29_sensitivity_grid_table.csv) | Sensitivity analysis results; file 29 = Table S1 | [11](../scripts-R/11_sensitivity_analysis.R) |
| 37–44 | Light: µ estimates, Monod parameters and fits (synthesis data) | [14](../scripts-R/14_synthesis_light.R) |
| 45–50 | Nitrogen: µ estimates, Monod parameters and fits | [15](../scripts-R/15_synthesis_nitrogen.R) |
| 51–55 | Phosphorus: µ estimates, Monod parameters and fits | [16](../scripts-R/16_synthesis_phosphorus.R) |
| 56–67 | Temperature: µ estimates, TPC parameters and fits | [17](../scripts-R/17_synthesis_temperature.R) |
| 70, 72, 74, 76, 78 | Compiled synthesis datasets by niche axis, and species summary | [18](../scripts-R/18_synthesis_pareto.R) |
| 71, 73, 75, 77 | Synthesis datasets with species specified | — |
| [`79_summary_fig_data.csv`](79_summary_fig_data.csv) | P values summarized for Fig. 5 | — |

Derived files can be regenerated by running [`00_master_script.R`](../scripts-R/00_master_script.R). Raw data files cannot be regenerated.