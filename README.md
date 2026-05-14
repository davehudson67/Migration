# Migration

This repository contains R code for analysing associations between migratory status, flyway/hemisphere groupings, extent of occurrence, phylogeny, and IUCN population trend status.

The main response variable is whether a species is classified as declining, derived from IUCN population trend categories. The main modelling approach uses binomial GLMs and phylogenetic logistic regression models fitted with `phyloglm`.

## Main directory files

### `AnalysisCode.R`

Early/main working analysis script. It loads the IUCN dataset, prepares hemisphere and flyway groupings, filters the data to species with usable population trend and EOO values, reads the phylogenetic tree, and fits a series of GLM and phylogenetic logistic regression models. It also explores different model formulations involving migratory status, combined hemisphere/flyway categories, EOO, and phylogenetic signal.

### `FullAnalysis_and_alternatePlots.R`

Most complete all-in-one analysis script. It prepares the data, scales and matches the phylogenetic tree, fits the main GLM and phylogenetic logistic regression models, saves fitted model objects to `Outputs/`, and generates the main summary figures. This is the best candidate for the main reproducible analysis workflow.

### `Combined_HemFly_DensityPlots_EOO&AlphaPlots_SummaryEffectPlots.R`

Plotting script that uses saved model objects from `Outputs/` to generate combined figure panels. It focuses on bootstrap-based visual summaries, including migratory effects, predicted probability of decline by hemisphere/flyway group, migrant–resident differences, EOO effects, and phylogenetic signal.

### `Single_HemFly_DensityPlots.R`

Plotting script for making small, focused density plots for a single selected hemisphere/flyway group. It compares resident and migrant predicted probabilities on the same axis and saves compact transparent-background figures for individual groups.


## Folders

### `Data/`

Contains the input datasets and phylogenetic tree objects used by the analysis scripts.

### `Outputs/`

Stores fitted model objects and generated figures produced by the analysis and plotting scripts.

### `OldCode/`

Contains older or superseded code retained for reference.
