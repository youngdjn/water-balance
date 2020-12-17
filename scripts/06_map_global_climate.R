library(tidyverse)
library(sf)
library(raster)


## load precip
precip_files = list.files("/home/derek/gis/Climate/worldclim/wc2.1_5m_prec/",pattern="tif",full.names = TRUE)
precip = stack(precip_files)

## load temp
tmin_files = list.files("/home/derek/gis/Climate/worldclim/wc2.1_5m_tmin/",pattern="tif",full.names = TRUE)
tmin = stack(tmin_files)

tmax_files = list.files("/home/derek/gis/Climate/worldclim/wc2.1_5m_tmax/",pattern="tif",full.names = TRUE)
tmax = stack(tmax_files)


## mean and sd precip
ppt_mean = mean(precip)
# ppt_sd = stackApply(precip,indices=1,fun=sd)
# ppt_cv = ppt_sd/ppt_mean

# # make a ppt_cv mask to exclude areas with very low precip
# low = quantile(ppt_mean,0.01)
# notlow_ppt = ppt_mean > low
# notlow_ppt[notlow_ppt == 0] = NA
# ppt_cv = ppt_cv * notlow_ppt

# writeRaster(ppt_mean,"/home/derek/gis/Climate/worldclim/summarized/ppt_mean.tif")
# writeRaster(ppt_sd,"/home/derek/gis/Climate/worldclim/summarized/ppt_sd.tif")

## mean and sd temp
tmin_mean = mean(tmin)
# tmin_meank = mean(tmin+273.15)
# tmin_sd = stackApply(tmin,indices=1,fun=sd)
# tmin_sdk = stackApply(tmin+273.15,indices=1,fun=sd)
# tmin_cv = tmin_sdk/(tmin_mean+273.15)

tmax_mean = mean(tmax)
# tmax_meank = mean(tmax+273.15)
# tmax_sd = stackApply(tmax,indices=1,fun=sd)
# tmax_sdk = stackApply(tmax+273.15,indices=1,fun=sd)
# tmax_cv = tmax_sdk/(tmax_mean+273.15)

# writeRaster(temp_mean,"/home/derek/gis/Climate/worldclim/summarized/temp_mean.tif")
# writeRaster(temp_sd,"/home/derek/gis/Climate/worldclim/summarized/temp_sd.tif")

### Stack all
clim = stack(ppt_mean,tmin_mean,tmax_mean,tmin[[1]],tmin[[7]],tmax[[1]],tmax[[7]],precip[[1]],precip[[7]])
names(clim) = c("ppt_mean","tmin_mean","tmax_mean","tmin_jan","tmin_jul","tmax_jan","tmax_jul","precip_jan","precip_jul")


## workspace saved here

#### Intersect with climate classifications ####

### For each zone get:
## the average intra-annual SD in precip, temp
## the average annual precip, temp

## crop to 33-45 deg lat N and S
coords = matrix(c(-179,30,
                  -179,45,
                  179,45,
                  179,30,
                  -179,30),
                ncol=2, byrow=TRUE)
northern = st_polygon(list(coords)) %>% st_sfc(crs=4326)

coords = matrix(c(-179,-30,
                  -179,-45,
                  179,-45,
                  179,-30,
                  -179,-30),
                ncol=2, byrow=TRUE)
southern = st_polygon(list(coords)) %>% st_sfc(crs=4326)

both = st_union(northern,southern)



## Load classifications
class = stack("/home/derek/repos/water-balance/data/koppen_clasifications/Beck_KG_V1_present_0p5.tif")
class_n = crop(class,northern %>% as("Spatial"))
class = crop(class,both %>% as("Spatial")) %>% mask(both %>% as("Spatial"))
class_poly = rasterToPolygons(class) %>% as("sf")
class_poly_n = rasterToPolygons(class_n) %>% as("sf")


### Once for whole globe

zone = class_poly[class_poly$layer=="4",] %>% as("sf") %>% st_union
summ_z4 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
summ_z4$zone = "4"

zone = class_poly[class_poly$layer=="14",] %>% as("sf") %>% st_union
summ_z14 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
summ_z14$zone = "14"

zone = class_poly[class_poly$layer=="7",] %>% as("sf") %>% st_union
summ_z7 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
summ_z7$zone = "7"

zone = class_poly[class_poly$layer=="29",] %>% as("sf") %>% st_union
summ_z29 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
summ_z29$zone = "29"

zone = class_poly[class_poly$layer=="9",] %>% as("sf") %>% st_union
summ_z9 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
summ_z9$zone = "9"

zone = class_poly[class_poly$layer=="18",] %>% as("sf") %>% st_union
summ_z18 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
summ_z18$zone = "18"

summ = bind_rows(summ_z4,summ_z7, summ_z14, summ_z29, summ_z9, summ_z18)


### Repeat for northern hemisphere only

zone = class_poly[class_poly_n$layer=="4",] %>% as("sf") %>% st_union
summ_z4 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
summ_z4$zone = "4N"

# zone = class_poly[class_poly_n$layer=="5",] %>% as("sf") %>% st_union
# summ_z5 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
# summ_z5$zone = "5N"

zone = class_poly[class_poly_n$layer=="14",] %>% as("sf") %>% st_union
summ_z14 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
summ_z14$zone = "14N"

zone = class_poly[class_poly_n$layer=="7",] %>% as("sf") %>% st_union
summ_z7 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
summ_z7$zone = "7N"

zone = class_poly[class_poly_n$layer=="29",] %>% as("sf") %>% st_union
summ_z29 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
summ_z29$zone = "29N"

summ_n = bind_rows(summ_z4,summ_z7, summ_z14, summ_z29)

summ_all = bind_rows(summ,summ_n)
summ_all


#### Functions to make climate ####

## Stretch a sin wave over 12 months with specified low and high
clim_cycle = function(jan, jul) {
  full_cycle_rad = 2*pi
  x = seq(from=0,to=full_cycle_rad,length.out=13)[1:12]
  #if(jan_high == FALSE) {x = x+pi}
  cycle = cos(x)
  ## scale the cycle to the min and max
  # make it a unit cycle
  cycle = (cycle+1)/2
  amplitude = jan-jul
  cycle = cycle*amplitude + jul
}


#### Make zone 4, desert climate ####

z4_tmin = clim_cycle(3.38,24.12)
z4_tmax = clim_cycle(15.58, 39.00)
z4_ppt = clim_cycle(18.29, 6.43)

dput(z4_tmin)
dput(z4_ppt)


#### Zone 7

z_tmin = clim_cycle(-8.52,13.69)
z_tmax = clim_cycle(3.35, 28.82)
z_ppt = clim_cycle(24.67, 35.32)

dput(z_tmin)
dput(z_tmax)
dput(z_ppt)


#### Zone 14

z_tmin = clim_cycle(-0.4639839,20.207419)
z_tmax = clim_cycle(10.0245175, 31.22183)
z_ppt = clim_cycle(81.72432, 125.640048)

dput(z_tmin)
dput(z_tmax)
dput(z_ppt)



#### Zone 29

z_tmin = clim_cycle(-22.0210038,2.081317)
z_tmax = clim_cycle(-7.6386230, 13.69826)
z_ppt = clim_cycle(11.68857, 84.840125)

dput(z_tmin)
dput(z_tmax)
dput(z_ppt)




#### Focal zone (Sierra Nevada)







# 
# #### Filtering ####
# 
# ## How much flexibility to grant in each variable?
# # 10% of the precip range
# ppt_flex = (cellStats(ppt_mean,stat="max") - cellStats(ppt_mean,stat="min"))*0.1
# 
# # 10% of the temperature range
# temp_flex = (cellStats(temp_mean,stat="max") - cellStats(temp_mean,stat="min"))*0.1
# 
# 
# ## Areas with CA precip&temp but low seasonality
# ppt_match = (ppt_mean > (120-ppt_flex)) & (ppt_mean < (120+ppt_flex))
# temp_match = (temp_mean > (10-temp_flex)) & (temp_mean < (10+temp_flex))
# mean_match = ppt_match & temp_match
# match = mean_match & ppt_cv < 0.2
# 
# 
# writeRaster(match,"/home/derek/gis/Climate/worldclim/summarized/testing_match_temporary.tif")
# 
# 
# 
# 









