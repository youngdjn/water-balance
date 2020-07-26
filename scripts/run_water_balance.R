library(raster)
library(tidyverse)
library(furrr)
plan(multiprocess,workers=8)

source("scripts/water_balance_models/set_wb.R")

## Load raster grid of input values
wb.inputs <- brick("data/wb_input/wb_input_vars.grd")

## Make a layer that contains cell IDs, for loading values back in later
id.layer <- wb.inputs[[1]]
values(id.layer) <- seq_along(id.layer)
wb.inputs[["ID"]] <- id.layer

## Convert to data frame
wb.inputs.df <- as.data.frame(wb.inputs)
wb.inputs.df <- as(wb.inputs,"SpatialPointsDataFrame")
# project to lat/lon
wb.inputs.df <- spTransform(wb.inputs.df,CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
wb.inputs.df$Lat <- coordinates(wb.inputs.df)[,2]
wb.inputs.df$Lon <- coordinates(wb.inputs.df)[,1]
wb.inputs.df <- as.data.frame(wb.inputs.df)

PET.methods <- "PT"
PET.mods <- c("STD","cc025")
AET.methods <- "Wil150mm"

## Define base parameter DF
wbparams_base <- wb.inputs.df[complete.cases(wb.inputs.df),] %>%
  mutate(scenario = "base")

## ALTERNATE PARAMETER DFS
## Aseasonal temp
tmin_mean = wbparams_base %>% select(starts_with("tmin.")) %>% rowMeans
tmean_mean = wbparams_base %>% select(starts_with("tmean.")) %>% rowMeans
tmax_mean = wbparams_base %>% select(starts_with("tmax.")) %>% rowMeans
td_mean = wbparams_base %>% select(starts_with("td.")) %>% rowMeans
wbparams_consttemp = wbparams_base %>%
  mutate_at(vars(starts_with("tmin.")), ~tmin_mean) %>%
  mutate_at(vars(starts_with("tmean.")), ~tmean_mean) %>%
  mutate_at(vars(starts_with("tmax.")), ~tmax_mean) %>%
  mutate_at(vars(starts_with("td.")), ~td_mean) %>%
  mutate(scenario = "consttemp")

## Aseasonal precip
ppt_mean = wbparams_base %>% select(starts_with("ppt")) %>% rowMeans
wbparams_constppt = wbparams_base %>%
  mutate_at(vars(starts_with("ppt.")), ~ppt_mean) %>%
  mutate(scenario = "constppt")

## Aseasonal temp and precip
wbparams_consttempppt = wbparams_base %>%
  mutate_at(vars(starts_with("ppt.")), ~ppt_mean) %>%
  mutate_at(vars(starts_with("tmin.")), ~tmin_mean) %>%
  mutate_at(vars(starts_with("tmean.")), ~tmean_mean) %>%
  mutate_at(vars(starts_with("tmax.")), ~tmax_mean) %>%
  mutate_at(vars(starts_with("td.")), ~td_mean) %>%
  mutate(scenario = "consttempppt")

## Double precip
wbparams_doubleppt = wbparams_base %>%
  mutate_at(vars(starts_with("ppt.")), ~.*2) %>%
  mutate(scenario = "doubleppt")

## Double aseasonal precip, aseasonal temp
wbparams_consttemppptdoubleppt = wbparams_consttempppt %>%
  mutate_at(vars(starts_with("ppt.")), ~.*2) %>%
  mutate(scenario = "consttemppptdoubleppt")


## LIST OF PARAMETERS DFs
inputs = list(base = wbparams_base,
              consttemp = wbparams_consttemp,
              constppt = wbparams_constppt,
              consttempppt = wbparams_consttempppt,
              doubleppt = wbparams_doubleppt,
              consttemppptdoubleppt = wbparams_consttemppptdoubleppt)

a = future_map_dfr(inputs,set_wb,PET.methods = PET.methods, PET.mods = PET.mods, AET.methods = AET.methods, monthly=FALSE, dobr.wb = TRUE,
               .progress = TRUE,
               .options = future_options(globals = c("PET.PT","Point_WB","monthlengths","SatVapPresSlope","ElevToPress","PET.mod.STD","PET.mod.cc025","AET.Wil150mm","AET.Wil","run_dob_wb","dobr_wb","snowmod","monthlyET0","aetmod"),
                                         packages = "dplyr"))

## Compile all into a long DF (with a column indicating what scenario)

a = set_wb(inputs[[1]],PET.methods = PET.methods, PET.mods = PET.mods, AET.methods = AET.methods, monthly=FALSE, dobr.wb = TRUE)

### pull together the output of the different wb methods and write to a file
write.csv(wb,"data/wb_output/landscape_wb_outputs.csv",row.names=FALSE)


#### Merge WB values with non-wb predictors and FIA responses ####
wb <- read.csv("data/wb_output/landscape_wb_outputs.csv",header=TRUE)
wb.nonwb <- left_join(wb,wb.inputs.df,by="ID")
d = wb.nonwb
write.csv(d,"data/wb_output/landscape_wb_inputs_outputs.csv",row.names=FALSE)

### rasterize
wb.nonwb <- read.csv("data/wb_output/landscape_wb_inputs_outputs.csv",header=TRUE,stringsAsFactors = FALSE)
d <- wb.nonwb


# template raster is wb.inputs, already loaded at beginning of script

id.layer <- wb.inputs[["ID"]]
#values(id.layer) <- seq_along(id.layer) # this is the way we assigned IDs to that grid

## Add all the modeled WB values/layers to the WB inputs DF so that outputs are stacked with inputs
## arrange the wb.nonwb DF with IDs in the same order as the raster cells
order <- match(raster::values(id.layer),d$ID)
d.rast.compare <- d[order,]

for (i in 1:ncol(d.rast.compare)) {
  col <- d.rast.compare[,i]
  name <- names(d.rast.compare)[i]
  
  wb.inputs[[name]] <- col
  if(!(name %in% names(wb.inputs))) { # if the layer we just added is not a name in the raster, it was just added at the end with a genericc name, so rename it
    names(wb.inputs)[length(names(wb.inputs))] <- name #set the last layer (the one just created) to the correct name
  }
    
}

d.rast <- wb.inputs

writeRaster(d.rast,"data/wb_output/rasters/landscape_wb_inputs_outputs.grd", overwrite=TRUE)


# ## output the AET and CWD rasters of the different crop coefficients
# writeRaster(d.rast$AET.Dobr.cc025,"data/wb_output/rasters/aet_025.tif",overwrite=TRUE)
# writeRaster(d.rast$AET.Dobr.cc100,"data/wb_output/rasters/aet_100.tif",overwrite=TRUE)
# writeRaster(d.rast$Deficit.Dobr.cc025,"data/wb_output/rasters/cwd_025.tif",overwrite=TRUE)
# writeRaster(d.rast$Deficit.Dobr.cc100,"data/wb_output/rasters/cwd_100.tif",overwrite=TRUE)
# 


