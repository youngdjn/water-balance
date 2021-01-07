# Water balance modeling and analysis scripts and data
Scripts and data supporting XXXXX et al. 2021. Inherent assumptions regarding vegetation physiology affect the utility of climatic water balance for ecological inference. Global Ecology and Biogeography.

This includes scripts for modeling water balance as described in the paper as well as for all downstream analyses (comparison of modeled water balance values between different PET coefficient assumptions; distribution models) and figures

These scripts are intended to be run with the working directory set to the root level of the repo (i.e., the folder this readme is in). Most of the scripts depend on the one prior to it (according to the numbering in the filename) having been run first to generate the necessary data.

## Description of scripts

### scripts/01_run_water_balance.R
Compute modified Dobrowski water balance values across the study landscape. It first calls a separate script (`scripts/01b_prep_water_balance_inputs.R`) to prep the inputs in the correct format (from the original input variable dataset which is a multilayer raster), and to compute inputs for alternative stylized climate regimes. Then it runs the water balance model, depending on functions in `scripts/water_balance_functions`.

This script relies on a raster stack of environmental variables as input to the water balance model. The expected structure of this raster stack is described in the header of this script.

Outputs (modeled water balance values, both tabular and raster format) are stored in `data/wb_output`.

This script (and the scripts it calls) implements a modification of the water balance model published by Dobrowski et al. 2013. The climate velocity of the contiguous United States during the 20th century. Global Change Biology 19:241-251. The modified Dobrowski water balance modeling functions, along with wrapper functions which facilitate their use, are in `scripts/water_balance_functions/dobrowski_wb.R`.



### scripts/02_compare_wb.R
For each climate scenario, compare the water balance outputs based on PET coef = 0.25 to the outputs based on PET coef = 1.00 and produce the mapped water balance figure (Fig. 3) and "heatmap" comparison figure (Fig. 5).

### scripts/03_summarize_gradient.R
Summarize landscape water balance along an elevation gradient and across aspects (Fig. 4).

### scripts/04_map_area.R
Make a map of the study area and a vicinity map (Fig. 2).

### scripts/05_predict_veg_crossval.R
Use computed water balance values to predict vegetation distributions across the study landscape, and compute model fit/validation and agreement for 1000 (or other specified) iterations of different training data partitions. Map predictions (Fig. 6).

### scripts/06_map_global_climate.R
Summarize WorldClim temperature and precip in different Koppen climate zones and build stylized climate regimes.

### scripts/07_pet_scaling_plots.R
Make the hypothetical PET scaling plots (Fig. 1).

## Data

Necessary input data (as described in the comments of each script that requires each data source) for the focal area of analysis are located in the `data` folder. This includes EVEG vegetation distribution data, climate and other abiotic environment data, study area perimeter, and cartographic basemap.

## Modeling water balance for a new area

These scripts can be adapted to model water balance using the Dobrowski method for a new location. There are two options for this.

1. Use `01_run_water_balance.R`, but first prepare an alternative input raster stack (`data/wb_input/wb_input_vars2.grd`) with the variables as described in the header of `01_run_water_balance.R`.

2. Use the function `run_dob_wb()` defined in `scripts/water_balance_functions/dobrowski_wb.R`. This function takes vectors of monthly values of environmental input variables, as described in the header of the function. This wrapper function is set up to model *climatic normal* water balance (that is, you supply it 12 months of climate normal data, and it computes climatic normal water balance, not water balance for specific individual years). This is a wrapper function around the functions published by Dobrowski et al. (2013). The versions of the Dobworski water balance functions included in this repo include the modifications described in the XXXXX et al. (2021) GEB paper that this repository accompanies. 

Both of these options will by default compute water balance two ways, one assuming a PET coefficient of 1.00 and one assuming a PET coefficient of 0.25.
