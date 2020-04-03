# Overview

Analysis scripts for processing Abeta deep mutational scanning (DMS) data.

# Required Software

To run the abetadms pipeline you will need the following software and associated packages:

* **[R](https://www.r-project.org/) >=v3.5.2** (Biostrings, caTools, corpcor, cowplot, data.table, fpc, gdata, ggplot2, GGally, hexbin, lemon, optparse, parallel, pdist, plyr, ppcor, raster, reshape2, Rpdb, RColorBrewer)

The following packages are optional:

* **[DiMSum](https://github.com/lehner-lab/DiMSum)** (pipeline for pre-processing deep mutational scanning data i.e. FASTQ to counts)
* **[DMS2structure](https://github.com/lehner-lab/DMS2structure)** (scripts used for epistasis and structure analysis of deep mutational scanning data in Schmiedel & Lehner, bioRxiv 2018)

# Installation and loading

Open R and enter:

```
# Install
if(!require(devtools)) install.packages("devtools")
devtools::install_github("lehner-lab/abetadms")

# Load
library(abetadms)

# Help
?abetadms
```

# Required Data

DiMSum fitness estimates and required miscellaneous files should be downloaded from [here]() to your project directory (see 'base_dir' argument) i.e. where output files should be written, and unzipped.

# Running

There are a number of options available for running the abetadms pipeline depending on user requirements.

* ## Basic (default)

Default pipeline functionality uses DiMSum fitness estimates (see 'Required Data'). Neither **[DiMSum](https://github.com/lehner-lab/DiMSum)** nor **[DMS2structure](https://github.com/lehner-lab/DMS2structure)** packages are required for this default functionality.

* ## Raw read processing

Raw read processing is not handled by the abetadms pipeline. FastQ files from paired-end sequencing of replicate deep mutational scanning (DMS) libraries before ('input') and after selection ('output') were processed using **[DiMSum](https://github.com/lehner-lab/DiMSum)** (manuscript in prep.), an R package that wraps common biological sequence processing tools.

* ## Preprocess fitness

Pipeline stage 1 ('abetadms_preprocess_fitness') reformats DiMSum files and re-estimates fitness of doubles mutants using a bayesian framework ('bayesian_double_fitness = T'). The latter is computationally intensive (~30minutes on 10 cores) and is therefore not run by default.

* ## Epistasis analysis

Pipeline stage 10 ('abetadms_epistasis_analysis') performs epistasis calculations. This stage is computationally intensive (~30minutes on 10 cores) and is therefore not run by default. **Note:** 'Required Data' (see above) already includes precomputed results of the epistasis analysis. However, to force re-execution of this stage set 'rerun_epistasis = T'. Additionally, the correct path to your local copy of the [DMS2structure](https://github.com/lehner-lab/DMS2structure) repository must be specified with 'DMS2structure_path = MY_LOCAL_PATH'.

* ## Structure analyses

Pipeline stages 11 and 12 ('abetadms_secondary_structure_predictions', 'abetadms_AB42_structure_propensities') perform secondary structure predictions and structure propensity calculations for PDB-structure derived
contact matrices respectively. Secondary structure predictions and propensity calculations are computationally intensive and are therefore not re-run by default. **Note:** 'Required Data' (see above) already includes precomputed results of the structure analyses. To force re-execution set 'rerun_structure = T'. Additionally, the correct path to your local copy of the [DMS2structure](https://github.com/lehner-lab/DMS2structure) repository must be specified with 'DMS2structure_path = MY_LOCAL_PATH'.

# Pipeline

The top-level function **abetadms()** is the recommended entry point to the pipeline. See section on "Required Data" above for instructions on how to obtain all required data and miscellaneous files before running the pipeline.

## Stage 1: Preprocess fitness

This stage ('abetadms_preprocess_fitness') reformats DiMSum files and re-estimates fitness of doubles mutants using a bayesian framework ('bayesian_double_fitness = T'). The latter is computationally intensive (~30minutes on 10 cores) and is therefore not run by default.

## Stage 2: Quality control plots

This stage ('abetadms_quality_control') produces quality control plots of fitness estimates.

## Stage 3: Combine fitness estimates

This stage ('abetadms_combine_fitness') performs normalisation of fitness estimates (based on silent mutants), fitness distribution plots and position-wise fitness plots.

## Stage 4: Calculate single and double mutant effects from AA PCA

This stage ('abetadms_aa_properties_mutant_effects') performs principal component analysis (PCA) of a curated collection of numerical indices representing various physicochemical and biochemical properties of amino acid (AA) properties. AA property feature values represent the difference between the WT and mutant PC scores.

## Stage 5: Calculate single and double mutant effects from aggregation tool predictions

This stage ('abetadms_agg_tools_mutant_effects') calculates aggregation / disorder algorithm feature values for single and double mutant variants (similar to stage 4).

## Stage 6: Single mutant heatmaps

This stage ('abetadms_single_mutant_heatmaps') produces single mutant heatmaps of fitness effects.

## Stage 7: Human disease mutations

This stage ('abetadms_human_disease_mutations') tests whether human disease mutations have biased fitness estimates.

## Stage 8: Dot plots showing explained variance of models to predict variant fitness

This stage ('abetadms_fitness_model_summary') produces plots of results from simple linear regression models to predict variant fitness.

## Stage 9: Helix propensity of WT and fitness hotspot line plot

This stage ('abetadms_wt_helix_propensity') plots helix propensity score and mean fitness effect along the length of WT Abeta.

## Stage 10: Epistasis analysis

This stage ('abetadms_epistasis_analysis') performs epistasis calculations. This stage is computationally intensive (~30minutes on 10 cores) and is therefore not run by default. To force re-execution of this stage set 'rerun_epistasis = T'. Additionally, the corect path to your local copy of the [DMS2structure](https://github.com/lehner-lab/DMS2structure) repository must be specified with 'DMS2structure_path = MY_LOCAL_PATH'.

## Stage 11: Secondary structure predictions

This stage ('abetadms_secondary_structure_predictions') performs secondary structure predictions and produces combined summary plots. Secondary structure predictions are computationally intensive and are therefore not re-run by default. To force re-execution of secondary structure predictions set 'rerun_structure = T'. Additionally, the corect path to your local copy of the [DMS2structure](https://github.com/lehner-lab/DMS2structure) repository must be specified with 'DMS2structure_path = MY_LOCAL_PATH'.

## Stage 12: PDB structure propensities

This stage ('abetadms_AB42_structure_propensities') performs structure propensity calculations for PDB-structure derived
contact matrices. Structure propensity calculation are computationally intensive and are therefore not re-run by default. To force re-execution of structure propensity calculations set 'rerun_structure = T'. Additionally, the corect path to your local copy of the [DMS2structure](https://github.com/lehner-lab/DMS2structure) repository must be specified with 'DMS2structure_path = MY_LOCAL_PATH'.

## Stage 13: PWI heatmaps

This stage ('abetadms_PWI_heatmaps') plots pair-wise interaction (PWI) score heatmaps.


