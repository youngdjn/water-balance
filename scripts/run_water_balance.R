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

wb.inputs.df.complete <- wb.inputs.df[complete.cases(wb.inputs.df),]

wb <- set_wb(wb.inputs.df.complete, out.name=out.name, PET.methods = PET.methods, PET.mods = PET.mods, AET.methods = AET.methods, monthly=FALSE, ret.soil.water=TRUE, dobr.wb = TRUE)












### pull together the output of the different wb methods and write to a file
# wb <- merge(wb.norad,wb.rad,by="ID")
write.csv(wb,paste0(out.name,"wb_output.csv"),row.names=FALSE)


#### Merge WB values with non-wb predictors and FIA responses ####
wb <- read.csv(paste0(out.name,"wb_output.csv"),header=TRUE)
wb.nonwb <- left_join(wb,wb.inputs.df,by="ID")

# ## Get the dobrowski and dingman columns
# dobr.cols <- grep("Dobr|alt",names(wb.nonwb))
# id.col <- grep("ID",names(wb.nonwb))
# keep.cols <- c(id.col,dobr.cols)
# d <- wb.nonwb[,keep.cols]
d = wb.nonwb

write.csv(d,paste0(out.name,"wb_nonwb_output.csv"),row.names=FALSE)

### rasterize

wb.nonwb <- read.csv(paste0(out.name,"wb_nonwb_output.csv"),header=TRUE,stringsAsFactors = FALSE)
d <- wb.nonwb


######!!!!! Need to fix something here with layer naming in output

# load the template raster (wb.inputs, already loaded at beginning of script)

id.layer <- wb.inputs[["ID"]]
#values(id.layer) <- seq_along(id.layer) # this is the way we assigned IDs to that grid

## arrange the wb.nonwb DF with IDs in the same order as the raster cells

order <- match(values(id.layer),d$ID)
d.rast.compare <- d[order,]

for (i in 1:ncol(d.rast.compare)) {
  col <- d.rast.compare[,i]
  name <- names(d.rast.compare)[i]
  
  wb.inputs[[name]] <- col
  names(wb.inputs)[length(names(wb.inputs))] <- name #set the last layer (the one just created) to the correct name
  
}



d.rast <- wb.inputs

## output the AET and CWD rasters of the different crop coefficients
writeRaster(d.rast$AET.Dobr.cc025,"data/wb_output/rasters/aet_025.tif",overwrite=TRUE)
writeRaster(d.rast$AET.Dobr.cc100,"data/wb_output/rasters/aet_100.tif",overwrite=TRUE)
writeRaster(d.rast$Deficit.Dobr.cc025,"data/wb_output/rasters/cwd_025.tif",overwrite=TRUE)
writeRaster(d.rast$Deficit.Dobr.cc100,"data/wb_output/rasters/cwd_100.tif",overwrite=TRUE)

# writeRaster(d.rast$AET.altwb.cc025,"data/wb_output/rasters/aet2_025.tif",overwrite=TRUE)
# writeRaster(d.rast$AET.altwb.cc100,"data/wb_output/rasters/aet2_100.tif",overwrite=TRUE)
# writeRaster(d.rast$Deficit.altwb.cc025,"data/wb_output/rasters/cwd2_025.tif",overwrite=TRUE)
# writeRaster(d.rast$Deficit.altwb.cc100,"data/wb_output/rasters/cwd2_100.tif",overwrite=TRUE)







d.rast$AET.Dobr.prop.025.100 <- (d.rast$AET.Dobr.cc025 -  d.rast$AET.Dobr.cc100) / d.rast$AET.Dobr.cc100
d.rast$AET.Dobr.prop.050.100 <- d.rast$AET.Dobr.cc050 / d.rast$AET.Dobr.cc100

plot(d.rast$AET.Dobr.prop.025.100)


res(d.rast)

writeRaster(d.rast$AET.Dobr.cc025,"aet_025.tif",overwrite=TRUE)
writeRaster(d.rast$AET.Dobr.cc100,"aet_100.tif",overwrite=TRUE)
writeRaster(d.rast$AET.Dobr.prop.025.100,"aet_prop.tif",overwrite=TRUE)

write.csv(d,"data/FIA/Dob_WB_simp.csv",row.names=FALSE)


