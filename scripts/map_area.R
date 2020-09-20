library(raster)
library(sf)
library(tidyverse)
library(ggthemes)
library(gridExtra)
library(ggnewscale)
library(viridis)
library(RStoolbox)

d = brick("data/wb_output/rasters/base.grd")

d = d %>% as("SpatialPointsDataFrame") %>% as("sf") %>% st_transform(3310)

d$x = st_coordinates(d)[,1]
d$y = st_coordinates(d)[,2]

## Load study area mask
tuol_mask <- st_read("data/tuolomne_mask/tuolomne_grid_mask.shp")

# ## Elevation panel
# p_elev <- ggplot(d) +
#   geom_raster(aes(x=x,y=y,fill=elev)) +
#   #coord_equal() +
#   theme_void() +
#   theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
#   #theme(strip.text=element_text(size=12),strip.background=element_blank()) +
#   #theme(panel.spacing=unit(1,"lines")) +
#   #theme(legend.text=element_text(size=9)) +
#   labs(title="   b) Elevation") +
#   #theme(plot.title=element_text(hjust=0.5,size=12)) +
#   scale_fill_continuous(guide=FALSE, low = "black",high="white") +
#   geom_sf(data=tuol_mask, fill=NA,color="red")
# 
# ## Radiation panel
# p_rad <- ggplot(d) +
#   geom_raster(aes(x=x,y=y,fill=rad.03)) +
#   #coord_equal() +
#   theme_void() +
#   theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
#   #theme(strip.text=element_text(size=12),strip.background=element_blank()) +
#   #theme(panel.spacing=unit(1,"lines")) +
#   #theme(legend.text=element_text(size=9)) +
#   labs(title="   c) Solar radiation") +
#   #theme(plot.title=element_text(hjust=0.5,size=12)) +
#   scale_fill_continuous(guide=FALSE, low = "black",high="white") +
#   geom_sf(data=tuol_mask, fill=NA,color="red")


## Combo panel
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
  #theme(strip.text=element_text(size=12),strip.background=element_blank()) +
  #theme(panel.spacing=unit(1,"lines")) +
  #theme(legend.text=element_text(size=9)) +
  #labs(title="b) Study area") +
  #theme(plot.title=element_text(hjust=0.5,size=12)) +
  geom_sf(data=tuol_mask, fill=NA,color="red") +
  coord_sf() +
  scale_x_continuous(expand=c(0,0), n.breaks = 2, breaks = waiver(), name = NULL) +
  scale_y_continuous(expand=c(0,0), n.breaks = 2, breaks = waiver(), name = NULL)



#### Make a california panel

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
  #labs(title="a) Vicinity map")
  



## California panel: just a blank plot with a title
# p_a = ggplot() +
#   labs(title = "a) Vicinity map") +
#   theme_void()
# 
# layout = rbind(c(NA,NA,2),
#                c(NA,1,2),
#                c(NA,1,3))
# 
# png("figures/context_map.png", width = 1600, height = 800, res = 200)
# grid.arrange(p_a, p_elev, p_rad, layout_matrix = layout, heights = c(.05,1,1.05), widths = c(.1,1,2))
# dev.off()


# 
# ## Alternative for single (combined) study area panel
# png("figures/context_map_comb.png", width = 1200, height = 540, res = 150)
# grid.arrange(a, p_comb, widths = c(1.2*1,3.4),ncol=2)
# dev.off()

png("figures/context_map_comb.png", width = 2*1400, height = 2*650, res = 2*150)
plot_grid(a,p_comb, rel_widths=c(1,3.4/1.2),labels=c("a) Vicinity map","b) Study area"),hjust=c(0,-0.5),label_fontface="plain")
dev.off()
