## Compute Dobrowski water balance values across the study landscape
## But first, call a separate script ("scripts/01b_prep_water_balance_inputs.R") to prep the inputs in the correct format (from the original input variable dataset which is a multilayer raster), and to compute inputs for alternative stylized climate regimes

## Input values are provided as a raster with parameters named as follows (where 'm' is replaced with the month number):
#		tmean.m, tmin.m, tmax.m		monthly mean, min, max temperature (degrees C)
#   td.m   monthly mean dewpoint temperature (degrees C)
#		ppt.m		monthly precipitation (mm)
#		rad.m		monthly fine-scale total solar radiation (W hrs per sq. m), computed following Dobrowski et al. 2013
#   rad.nldas.m   monthly coarse-scale total solar radiation (W hrs per sq. m) from NLDAS-2 dataset
#   wind.m    monthly average wind speed in m/s at 10m above ground
#   elev    elevation (m)

# For an example input raster with properly named layers, see data/wb_input/wb_input_vars2.grd

library(raster)
library(tidyverse)

source("scripts/water_balance_functions/set_wb.R")


## Take input parameter values from a raster and prepare them in the format (data frame) needed to input to water balance model
## Also apply alternative climate scenarios
## This is all done in 01b_prep_water_balance_inputs.R, which loads those input data frames into memory
## It also loads the relevant rasters into memory so they can be used at the bottom of this script as a template for making rasters of modeled water balance
## See that script for a description of the layers needed in that input data raster
source("scripts/01b_prep_water_balance_inputs.R")



## All of these data frames below (wbparams_*) were loaded into memory by the previous line. They are the data frames of water balance inputs for each climate scenario.
# Put them into a list so we can map water balance computation over them
inputs = list(base = wbparams_base,
              z04 = wbparams_z04,
              z07 = wbparams_z07,
              z11 = wbparams_z11,
              z14 = wbparams_z14,
              z29 = wbparams_z29)



# Apply the set water balance function (to compute water balance for every grid cell) to each set of inputs (climate regimes)
# The set_wb function is defined in scripts/water_balance_functions/set_wb.R
wb = map(inputs,set_wb)
# Columns returned ending in cc100 and cc025 indicate water balance values modeled assuming a PET coefficient (or crop coefficient, "cc") of 1.00 and 0.25 respectively.



# Compile and save data frames of inputs and outputs, plus save rasters: one for each climate scenario

for(i in 1:length(inputs)) {
  scenario = names(inputs)[i]
  inputs_foc = inputs[[scenario]]
  outputs_foc = wb[[scenario]]
  in_out_foc = left_join(inputs_foc, outputs_foc,by="ID") %>%
    select(-scenario)
  write.csv(in_out_foc,paste0("data/wb_output/",scenario,".csv"),row.names=FALSE)
  
  
  
  ## Add all the modeled WB values/layers to the WB inputs DF so that outputs are stacked with inputs
  
  ## arrange the wb.nonwb DF with IDs in the same order as the raster cells

  order <- match(raster::values(id.layer),in_out_foc$ID)
  
  d.rast.compare <- in_out_foc[order,]
  
  wb.inputs.foc = wb.inputs
  
  for (j in 1:ncol(d.rast.compare)) {
    col <- d.rast.compare[,j]
    name <- names(d.rast.compare)[j]

    if(!(name %in% names(wb.inputs.foc))) { # if the layer we are adding is not a name in the raster, it was just added at the end (newly computed wb model output) with a genericc name, so rename it
      newlayer = wb.inputs.foc[[1]] * 0
      values(newlayer) = col
      names(newlayer) = name
      wb.inputs.foc = addLayer(wb.inputs.foc,newlayer)
    } else { # it already exists, so replace it
      wb.inputs.foc[[name]] <- col
    }
    
  }
  
  
  writeRaster(wb.inputs.foc,paste0("data/wb_output/rasters/",scenario,".grd"), overwrite=TRUE)
  
  
}