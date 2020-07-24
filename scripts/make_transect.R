library(raster)
library(sf)
library(tidyverse)

d = brick("data/wb_output/rasters/landscape_wb_inputs_outputs.grd")

rad = d$rad.03
elev = d$elev

### Splitting the elevation into 10 bins, getting the N and S extremes for each
aspect_masks = list()

for(bin in 0:9) {
  elev_low = quantile(elev,bin*.1)
  elev_high = quantile(elev,bin*.1+0.1)
  
  cells_in_elev = elev > elev_low & elev < elev_high
  rad_low = quantile(rad[cells_in_elev],0.10)
  rad_high = quantile(rad[cells_in_elev],0.90)
  
  s_facing = rad > rad_high & cells_in_elev
  n_facing = rad < rad_low & cells_in_elev
  
  aspect_mask = raster(elev)
  values(aspect_mask) = 0
  aspect_mask[s_facing] = 2
  aspect_mask[n_facing] = 1
  aspect_mask[!s_facing & !n_facing] = NA
  
  aspect_masks[[as.character(bin)]] = aspect_mask
}

aspect_mask = merge(aspect_masks[[1]],aspect_masks[[2]],aspect_masks[[3]],aspect_masks[[4]],aspect_masks[[5]],
                    aspect_masks[[6]],aspect_masks[[7]],aspect_masks[[8]],aspect_masks[[9]],aspect_masks[[10]])

#append it to the stack of all values
d = stack(d,aspect_mask)
names(d)[nlayers(d)] = "aspect_mask"

### at each cell, get WB values and N-S values
# turn stack to data frame
d_df = as.data.frame(d) %>%
  filter(!is.na(aspect_mask) & aspect_mask != 0) %>% # keep only N and S
  mutate(radiation = recode(aspect_mask,"2" = "high", "1" = "low"))

### within N and S, compute mean, lwr, and upr for each elevation band

## make elevation bands
breaks = seq(400,3100,by=200)
labels = breaks[1:(length(breaks)-1)] + 100
d_df = d_df %>%
  mutate(elev_bins = cut(elev,breaks =breaks , labels = labels)) %>%
  filter(!is.na(elev_bins))

lwr = function(x) {
  quantile(x,0.1)
}

upr = function(x) {
  quantile(x,0.9)
}

d_df_summ = d_df %>%
  group_by(radiation,elev_bins) %>%
  summarize_at(vars(contains(c("AET","PET","Deficit"))), funs(mean,lwr,upr)) %>%
  pivot_longer(cols = c(-elev_bins,-radiation)) %>%
  separate(name,into=c("wb_metric","stat"),sep="_") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  ungroup()


### Make plots ###

d_plot = d_df_summ %>%
  filter(wb_metric %in% c("AET.Dobr.cc025",
                          "AET.Dobr.cc100")) %>%
  mutate(wb_metric = recode(wb_metric,"AET.Dobr.cc100" = "1","AET.Dobr.cc025" =  "0.25")) %>%
  mutate(wb_metric = factor(wb_metric,levels=c("1","0.25"))) %>%
  mutate(radiation = recode(radiation,"high" = "South","low" =  "North") %>% as.factor) %>%
  mutate(radiation = factor(radiation,levels=c("South","North"))) %>%
  mutate(elev = as.numeric(as.character(elev_bins)))

p = ggplot(d_plot,aes(x=elev,y=mean,color=wb_metric,linetype=radiation, group=interaction(wb_metric,radiation))) +
  geom_line(size=1) +
  scale_color_manual(name="PET coefficient", values = c("slateblue","darkorange")) +
  scale_linetype_discrete(name="Aspect") +
  theme_bw(12) +
  labs(x = "Elevation (m)",
       y = "Annual AET (mm)")

png("figures/transect_gradient.png", width=1100,height=800,res=200)
p
dev.off()

ggplot(d_df,aes(x=elev,y=AET.Dobr.cc100, color=aspect_mask)) +
  geom_point()




