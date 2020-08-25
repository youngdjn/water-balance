library(raster)
library(sf)
library(tidyverse)
library(ggthemes)
library(gridExtra)

d = brick("data/wb_output/rasters/base.grd")

d = d %>% as("SpatialPointsDataFrame") %>% as("sf") %>% st_transform(3310)

d$x = st_coordinates(d)[,1]
d$y = st_coordinates(d)[,2]

## Load study area mask
tuol_mask <- st_read("data/tuolomne_mask/tuolomne_grid_mask.shp")

## Elevation panel
p_elev <- ggplot(d) +
  geom_raster(aes(x=x,y=y,fill=elev)) +
  #coord_equal() +
  theme_void() +
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
  #theme(strip.text=element_text(size=12),strip.background=element_blank()) +
  #theme(panel.spacing=unit(1,"lines")) +
  #theme(legend.text=element_text(size=9)) +
  labs(title="   b) Elevation") +
  #theme(plot.title=element_text(hjust=0.5,size=12)) +
  scale_fill_continuous(guide=FALSE, low = "black",high="white") +
  geom_sf(data=tuol_mask, fill=NA,color="red")

## Radiation panel
p_rad <- ggplot(d) +
  geom_raster(aes(x=x,y=y,fill=rad.03)) +
  #coord_equal() +
  theme_void() +
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
  #theme(strip.text=element_text(size=12),strip.background=element_blank()) +
  #theme(panel.spacing=unit(1,"lines")) +
  #theme(legend.text=element_text(size=9)) +
  labs(title="   c) Solar radiation") +
  #theme(plot.title=element_text(hjust=0.5,size=12)) +
  scale_fill_continuous(guide=FALSE, low = "black",high="white") +
  geom_sf(data=tuol_mask, fill=NA,color="red")


## California panel: just a blank plot with a title
p_a = ggplot() +
  labs(title = "a) Context map") +
  theme_void()

layout = rbind(c(NA,NA,2),
               c(NA,1,2),
               c(NA,1,3))

png("figures/context_map.png", width = 1600, height = 800, res = 200)
grid.arrange(p_a, p_elev, p_rad, layout_matrix = layout, heights = c(.05,1,1.05), widths = c(.1,1,2))
dev.off()
