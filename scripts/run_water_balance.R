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
ppt_mean = wbparams_base %>% select(starts_with("ppt.")) %>% rowMeans
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

wb = future_map(inputs,set_wb,PET.methods = PET.methods, PET.mods = PET.mods, AET.methods = AET.methods, monthly=FALSE, dobr.wb = TRUE,
               .progress = TRUE,
               .options = future_options(globals = c("PET.PT","Point_WB","monthlengths","SatVapPresSlope","ElevToPress","PET.mod.STD","PET.mod.cc025","AET.Wil150mm","AET.Wil","run_dob_wb","dobr_wb","snowmod","monthlyET0","aetmod"),
                                         packages = "dplyr"))

# Compile and save Data frames of inputs and outputs , plus rasters: one for each scenario

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

    if(!(name %in% names(wb.inputs.foc))) { # if the layer we are adding is not a name in the raster, it was just added at the end with a genericc name, so rename it
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



