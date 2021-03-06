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

zone = class_poly[class_poly_n$layer=="11",] %>% as("sf") %>% st_union
summ_z11 = extract(clim,zone %>% as("Spatial"),fun=mean, na.rm=TRUE) %>% as.data.frame
summ_z11$zone = "11N"

summ_n = bind_rows(summ_z4,summ_z7, summ_z14, summ_z29, summ_z11)

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

z7_tmin = clim_cycle(-8.52,13.69)
z7_tmax = clim_cycle(3.35, 28.82)
z7_ppt = clim_cycle(24.67, 35.32)

dput(z7_tmin)
dput(z7_tmax)
dput(z7_ppt)


#### Zone 14

z14_tmin = clim_cycle(-0.4639839,20.207419)
z14_tmax = clim_cycle(10.0245175, 31.22183)
z14_ppt = clim_cycle(81.72432, 125.640048)

dput(z14_tmin)
dput(z14_tmax)
dput(z14_ppt)



#### Zone 29

z29_tmin = clim_cycle(-22.0210038,2.081317)
z29_tmax = clim_cycle(-7.6386230, 13.69826)
z29_ppt = clim_cycle(11.68857, 84.840125)

dput(z29_tmin)
dput(z29_tmax)
dput(z29_ppt)


#### Zone 11

z11_tmin = clim_cycle(-0.3293741,21.581370)
z11_tmax = clim_cycle(8.4688886, 29.70448)
z11_ppt = clim_cycle(25.15121, 160.344716)

dput(z11_tmin)
dput(z11_tmax)
dput(z11_ppt)




#### Focal zone (Sierra Nevada)

zfoc_tmin = structure(list(tmin.01 = -2.41072513124655, tmin.02 = -2.09452284281255, 
                           tmin.03 = -0.928372148313223, tmin.04 = 1.00936144375414, 
                           tmin.05 = 4.94143029052125, tmin.06 = 9.10159145224969, tmin.07 = 12.8864247871696, 
                           tmin.08 = 12.4547131428914, tmin.09 = 9.62951810396863, tmin.10 = 5.07600902808213, 
                           tmin.11 = 0.399353332989514, tmin.12 = -2.39636750663098), row.names = c(NA, 
                                                                                                    -1L), class = "data.frame") %>% unname %>% as.numeric

zfoc_tmax = structure(list(tmax.01 = 8.3470455881612, tmax.02 = 9.07445750112484, 
                           tmax.03 = 11.2493209703664, tmax.04 = 13.9233650361225, tmax.05 = 18.6152324031166, 
                           tmax.06 = 23.4542322550366, tmax.07 = 27.6599998407719, tmax.08 = 27.3078510885403, 
                           tmax.09 = 24.0352129744094, tmax.10 = 18.499809379395, tmax.11 = 11.9916004397836, 
                           tmax.12 = 8.10329069557953), row.names = c(NA, -1L), class = "data.frame") %>% unname() %>% as.numeric

zfoc_ppt = structure(list(ppt.01 = 199.006025727205, ppt.02 = 179.143073337826, 
                          ppt.03 = 163.420968920676, ppt.04 = 84.038028271778, ppt.05 = 49.743536874516, 
                          ppt.06 = 14.5905457488015, ppt.07 = 7.02463945255649, ppt.08 = 4.75581887294955, 
                          ppt.09 = 18.7976216164704, ppt.10 = 61.5426818049231, ppt.11 = 115.551510455076, 
                          ppt.12 = 172.168588665293), row.names = c(NA, -1L), class = "data.frame") %>% unname() %>% as.numeric



#### Compile into a data frame for plotting

clim = rbind(z4_tmin,
             z4_tmax,
             z4_ppt,
             z7_tmin,
             z7_tmax,
             z7_ppt,
             z11_tmin,
             z11_tmax,
             z11_ppt,
             z14_tmin,
             z14_tmax,
             z14_ppt,
             z29_tmin,
             z29_tmax,
             z29_ppt,
             zfoc_tmin,
             zfoc_tmax,
             zfoc_ppt) %>% as.data.frame

colnames(clim) = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

clim$var = c("tmin","tmax","ppt","tmin","tmax","ppt","tmin","tmax","ppt","tmin","tmax","ppt","tmin","tmax","ppt","tmin","tmax","ppt")
clim$zone = c(4,4,4,7,7,7,11,11,11,14,14,14,29,29,29,"Focal","Focal","Focal")
clim$zone = paste0("Zone ",clim$zone)
clim[clim$zone == "Zone Focal","zone"] = "Tuolumne"
clim$zone = factor(clim$zone,levels=c("Tuolumne",paste0("Zone ",c(4,7,11,14,29))))

clim = pivot_longer(clim,cols=c(-var,-zone),names_to="month") %>%
  pivot_wider(names_from=var,values_from=value)

clim = clim %>%
  mutate(month = factor(month,levels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")))



#### Plot them

g = ggplot(clim,aes(x=month,y=ppt, group=1)) +
  geom_bar(stat="identity") +
  geom_line(aes(x=month,y=(tmin+25)*3),stat="identity") +
  geom_line(aes(x=month,y=(tmax+25)*3),stat="identity") +
  facet_wrap(vars(zone)) +
  scale_y_continuous("Precipitation (mm)",
                     sec.axis = sec_axis(~(./3)-25, name="Temperature (°C)")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

png(paste0("figures/climatographs.png"), width = 1400,height = 1000, res=250)
g
dev.off()








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









