# Make a map of the study area and a vicinity map (Fig. 2)

library(raster)
library(sf)
library(tidyverse)
library(ggthemes)
library(gridExtra)
library(ggnewscale)
library(viridis)
library(RStoolbox)
library(cowplot)

d_rast = brick("data/wb_output/rasters/base.grd")

d = d_rast %>% as("SpatialPointsDataFrame") %>% as("sf") %>% st_transform(3310)

d$x = st_coordinates(d)[,1]
d$y = st_coordinates(d)[,2]

## Load study area mask
tuol_mask <- st_read("data/tuolomne_mask/tuolomne_grid_mask.shp")


### Make contour lines
elev = d_rast[["elev"]]
contour = rasterToContour(elev,levels=c(500,1000,1500,2000,2500,3000,3500))


## Combo panel (elevation using color and radiation using grayscale)
p_comb = ggplot(d) +
  geom_raster(aes(x=x,y=y,fill=rad.03/1000), alpha=1) +
  scale_fill_continuous(low = "black",high="white", name = bquote('Insolation'~(kWh~m^-2~d^-1)), guide = guide_colorbar(title.position = "top", title.hjust = 0.5)) +

  new_scale("fill") +
  geom_raster(aes(x=x,y=y,fill=elev), alpha=0.4) +
  scale_fill_viridis(name = bquote('Elevation'~(m)^phantom(2)), guide = guide_colorbar(title.position = "top", title.hjust = 0.5), direction=-1) +

  theme_bw() +
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"),
        legend.position = "bottom", legend.title = element_text(size=10),
        legend.key.width = unit(0.8,"cm")) +
  geom_sf(data=tuol_mask, fill=NA,color="red") +
  geom_sf(data=contour %>% as("sf"), color="black",size=0.2,alpha=1) +
  coord_sf() +
  scale_x_continuous(expand=c(0,0), n.breaks = 2, breaks = waiver(), name = NULL) +
  scale_y_continuous(expand=c(0,0), n.breaks = 2, breaks = waiver(), name = NULL)



#### Make a california (vicinity map) panel

basemap = stack("data/California_TIFF_Basemap.tiff")

### get bbox as polygon
# points
points = st_as_sf(d,coords=c("x","y"),crs=3310)
bbox = st_bbox(points) %>% st_as_sfc

e = extent(basemap)

a = ggplot(data=bbox) +
  ggRGB(img = basemap,r=1,g=2,b=3,ggLayer=TRUE) +
  geom_sf(fill=NA,color="red") +
  xlim(e[1:2]) +
  ylim(c(-950000,450000)) +
  theme_void()

## Write the full two-panel figure
png("figures/context_map_comb.png", width = 2*1400, height = 2*650, res = 2*150)
plot_grid(a,p_comb, rel_widths=c(1,3.4/1.2),labels=c("a) Vicinity map","b) Study area"),hjust=c(0,-0.5),label_fontface="plain")
dev.off()
