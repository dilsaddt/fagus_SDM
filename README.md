# fagus_project
Spatial distribution of *Fagus orientalis* (Oriental beech) in Eurasia with SDM by using biomod2 package
Model trained with current climate data. After that predictions were done for past projections: Last Glacial Maximum (~11k bp) and Mid-Holocene (~5k bp), and future projections 2050 and 2070. For future projections different IPCC future climate scenarios were applied: RCP 2.6 good scenario, RCP 4.5 moderata scenario, and RCP 8.5 bad scenario.

**bioclim.R**: preparation of bioclimatic variables downloaded from WorldClim for current climate

**new_climate.R**: preparation of bioclimatic variables downloaded from WorldClim for
- current climate
- MIROC Global Climate Model - Last Glacial Maximum
- MIROC Global Climate Model - Mid-Holocene
- MIROC Global Climate Model - 2050 RCP 2.6
- MIROC Global Climate Model - 2050 RCP 4.5
- MIROC Global Climate Model - 2050 RCP 8.5
- MIROC Global Climate Model - 2070 RCP 2.6
- MIROC Global Climate Model - 2070 RCP 4.5
- MIROC Global Climate Model - 2070 RCP 8.5
- CCSM4 Global Climate Model - Last Glacial Maximum
- CCSM4 Global Climate Model - Mid-Holocene
- CCSM4 Global Climate Model - 2050 RCP 2.6
- CCSM4 Global Climate Model - 2050 RCP 4.5
- CCSM4 Global Climate Model - 2050 RCP 8.5
- CCSM4 Global Climate Model - 2070 RCP 2.6
- CCSM4 Global Climate Model - 2070 RCP 4.5
- CCSM4 Global Climate Model - 2070 RCP 8.5

**biomod2_models.R**: original models with biomod2 package using 5 algorithms - GLM, GAM, RF, SRE (BIOCLIM), and MAXENT

**makale_rf_mart19_run.R**: re-run of 2070 MIROC

**biomod2_models_env.RData**: model output files, R environment
