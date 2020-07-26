library(raster)
library(sf)
library(tidyverse)
library(scales)
library(viridis)

d_rast = brick("data/wb_output/rasters/landscape_wb_inputs_outputs.grd")

d_orig = d_rast %>% as("SpatialPointsDataFrame") %>% as("sf") %>% st_transform(3310)

d_orig$x = st_coordinates(d_orig)[,1]
d_orig$y = st_coordinates(d_orig)[,2]

## Compute differences between cc levels

trunc5 = function(x) {
  x = ifelse(x > 5, NA, x)
  return (x)
}

d = d_orig %>%
  rename(dob_aet100 = AET.Dobr.cc100,
         dob_aet025 = AET.Dobr.cc025,
         dob_cwd100 = Deficit.Dobr.cc100,
         dob_cwd025 = Deficit.Dobr.cc025,
         wil_aet100 = AET.PT.STD.Wil150mm,
         wil_aet025 = AET.PT.cc025.Wil150mm,
         wil_cwd100 = Deficit.PT.STD.Wil150mm,
         wil_cwd025 = Deficit.PT.cc025.Wil150mm) %>%
  # rescale AET and CWD so that they are stretched from 0 to 1
  mutate_at(vars(starts_with("dob_"), starts_with("wil_")), rescale, to=c(0,1)) %>%
  mutate(diff_dob_aet = dob_aet100 - dob_aet025,
         diff_dob_cwd = dob_cwd100 - dob_cwd025,
         diff_wil_aet = wil_aet100 - wil_aet025,
         diff_wil_cwd = wil_cwd100 - wil_cwd025)


## bin it myself and require a minimum number of rows per tile

elev_breaks = seq(min(d$elev),max(d$elev),length.out=6)
elev_mids = elev_breaks[1:(length(elev_breaks)-1)] + 0.5*(elev_breaks[2]-elev_breaks[1])
rad_breaks = seq(min(d$rad.03),max(d$rad.03),length.out=6)
rad_mids = rad_breaks[1:(length(rad_breaks)-1)] + 0.5*(rad_breaks[2]-rad_breaks[1])

d_plot = d

st_geometry(d_plot) = NULL

d_plot = d_plot %>%
  mutate(elev_cat = cut(elev,elev_breaks,labels=elev_mids),
         rad_cat = cut(rad.03,rad_breaks,labels=rad_mids)) %>%
  filter(!is.na(elev_cat)) %>%
  group_by(elev_cat,rad_cat)

d_plot1 = d_plot %>%
  summarize(n = n())
  
d_plot = d_plot %>%
  summarize_at(vars(contains("dob_"),contains("wil_")), mean)   %>%
  left_join(d_plot1) %>%
  filter(n > 50) %>%
  ungroup() %>%
  mutate(rad_cat = rad_cat %>% as.character %>% as.numeric,
         elev_cat = elev_cat %>% as.character %>% as.numeric)



## For Dobr100, plot AET and CWD diff as a function of elev and rad

ggplot(d_plot,aes(x=elev_cat,y=rad_cat,fill=diff_wil_aet)) +
  geom_tile() +
  scale_fill_viridis()


## Map the difference between scaled layers
ggplot(d,aes(x=x,y=y,fill=diff_dob_aet)) +
  geom_raster() +
  scale_fill_viridis() +
  coord_equal()
