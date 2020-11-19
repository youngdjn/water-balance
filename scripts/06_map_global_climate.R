library(tidyverse)
library(sf)
library(raster)


## load precip
precip_files = list.files("/home/derek/gis/Climate/worldclim/wc2.1_5m_prec/",pattern="tif",full.names = TRUE)
precip = stack(precip_files)

## load temp
temp_files = list.files("/home/derek/gis/Climate/worldclim/wc2.1_5m_tavg/",pattern="tif",full.names = TRUE)
temp = stack(temp_files)

## mean and sd precip
ppt_mean = mean(precip)
ppt_sd = stackApply(precip,indices=1,fun=sd)
ppt_cv = ppt_sd/ppt_mean

# make a ppt_cv mask to exclude areas with very low precip
low = quantile(ppt_mean,0.01)
notlow_ppt = ppt_mean > low
notlow_ppt[notlow_ppt == 0] = NA
ppt_cv = ppt_cv * notlow_ppt

writeRaster(ppt_mean,"/home/derek/gis/Climate/worldclim/summarized/ppt_mean.tif")
writeRaster(ppt_sd,"/home/derek/gis/Climate/worldclim/summarized/ppt_sd.tif")

## mean and sd temp
temp_mean = mean(temp)
temp_sd = stackApply(temp,indices=1,fun=sd)
temp_sdk = stackApply(temp+273.15,indices=1,fun=sd)
temp_cv = temp_sdk/(temp_mean+273.15)

writeRaster(temp_mean,"/home/derek/gis/Climate/worldclim/summarized/temp_mean.tif")
writeRaster(temp_sd,"/home/derek/gis/Climate/worldclim/summarized/temp_sd.tif")

### Stack all
clim = stack(ppt_mean,ppt_sd,ppt_cv,temp_mean,temp_sd,temp_cv)
names(clim) = c("ppt_mean","ppt_sd","ppt_cv","temp_mean","temp_sd","temp_cv")


#### Filtering ####

## How much flexibility to grant in each variable?
# 10% of the precip range
ppt_flex = (cellStats(ppt_mean,stat="max") - cellStats(ppt_mean,stat="min"))*0.1

# 10% of the temperature range
temp_flex = (cellStats(temp_mean,stat="max") - cellStats(temp_mean,stat="min"))*0.1


## Areas with CA precip&temp but low seasonality
ppt_match = (ppt_mean > (120-ppt_flex)) & (ppt_mean < (120+ppt_flex))
temp_match = (temp_mean > (10-temp_flex)) & (temp_mean < (10+temp_flex))
mean_match = ppt_match & temp_match
match = mean_match & ppt_cv < 0.2


writeRaster(match,"/home/derek/gis/Climate/worldclim/summarized/testing_match_temporary.tif")













