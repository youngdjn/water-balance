library(raster)
library(tidyverse)

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
PET.mods <- c("STD","cc050","cc025")
AET.methods <- "Wil150mm"

out.name <- "data/wb_output/landscape_wb"

plots.non.wb.df.complete <- plots.non.wb.df[complete.cases(plots.non.wb.df),]

wb <- set_wb(plots.non.wb.df.complete, out.name=out.name, PET.methods = PET.methods, PET.mods = PET.mods, AET.methods = AET.methods, monthly=FALSE, ret.soil.water=TRUE, dobr.wb = TRUE)

