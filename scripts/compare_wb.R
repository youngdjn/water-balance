library(raster)
library(sf)
library(tidyverse)
library(scales)
library(viridis)
library(ggthemes)
library(gridExtra)

d_base = brick("data/wb_output/rasters/base.grd")
d_constppt = brick("data/wb_output/rasters/constppt.grd")
d_consttemp = brick("data/wb_output/rasters/consttemp.grd")
d_consttempppt = brick("data/wb_output/rasters/consttempppt.grd")
d_doubleppt = brick("data/wb_output/rasters/doubleppt.grd")
d_consttemppptdoubleppt = brick("data/wb_output/rasters/consttemppptdoubleppt.grd")

d_base = d_base %>% as("SpatialPointsDataFrame") %>% as("sf") %>% st_transform(3310) %>%
  mutate(scenario = "base")
d_constppt = d_constppt %>% as("SpatialPointsDataFrame") %>% as("sf") %>% st_transform(3310) %>%
  mutate(scenario = "constppt")
d_consttemp = d_consttemp %>% as("SpatialPointsDataFrame") %>% as("sf") %>% st_transform(3310) %>%
  mutate(scenario = "consttemp")
d_consttempppt = d_consttempppt %>% as("SpatialPointsDataFrame") %>% as("sf") %>% st_transform(3310) %>%
  mutate(scenario = "consttempppt")
d_doubleppt = d_doubleppt %>% as("SpatialPointsDataFrame") %>% as("sf") %>% st_transform(3310) %>%
  mutate(scenario = "doubleppt")
d_consttemppptdoubleppt = d_consttemppptdoubleppt %>% as("SpatialPointsDataFrame") %>% as("sf") %>% st_transform(3310) %>%
  mutate(scenario = "consttemppptdoubleppt")

d_orig = rbind(d_base,d_constppt,d_consttemp,d_consttempppt,d_doubleppt,d_consttemppptdoubleppt)

d_orig$x = st_coordinates(d_orig)[,1]
d_orig$y = st_coordinates(d_orig)[,2]

## Compute differences between cc levels

d = d_orig %>%
  rename(dob_aet100 = AET.Dobr.cc100,
         dob_aet025 = AET.Dobr.cc025,
         dob_cwd100 = Deficit.Dobr.cc100,
         dob_cwd025 = Deficit.Dobr.cc025,
         wil_aet100 = AET.PT.STD.Wil150mm,
         wil_aet025 = AET.PT.cc025.Wil150mm,
         wil_cwd100 = Deficit.PT.STD.Wil150mm,
         wil_cwd025 = Deficit.PT.cc025.Wil150mm) %>%
  # rescale AET and CWD (and solar exposure) so that they are stretched from 0 to 1
  mutate_at(vars(starts_with("dob_"), starts_with("wil_")), funs(scaled = scale)) %>%
  mutate(rad.03_scaled=rescale(rad.03, to=c(0,1))) %>%
  mutate(diff_dob_aet = dob_aet100_scaled - dob_aet025_scaled,
         diff_dob_cwd = dob_cwd100_scaled - dob_cwd025_scaled,
         diff_wil_aet = wil_aet100_scaled - wil_aet025_scaled,
         diff_wil_cwd = wil_cwd100_scaled - wil_cwd025_scaled) %>%
  # remove model cells that didn't run
  filter(wil_aet100 > -1 & wil_aet025 > -1)







### Make a 6-panel plot of AET and CWD for each method (and each scenario)



mapfun <- function(d,column,title,scale.limits=NULL, colorscale, colorscale_dir) {
  
  if(colorscale == "magma") {
    begin = 0.1
    end = 0.9
  } else {
    begin = 0
    end = 1
  }
  
  p <- ggplot(d) +
    geom_raster(aes(x=x,y=y,fill=!!sym(column))) +
    coord_equal() +
    theme_map() +
    theme(legend.position=c(.008,.35)) +
    theme(legend.title=element_blank()) +
    theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
    theme(strip.text=element_text(size=12),strip.background=element_blank()) +
    theme(panel.spacing=unit(1,"lines")) +
    theme(panel.border=element_rect(fill=NA)) +
    theme(legend.text=element_text(size=9)) +
    labs(title=title) +
    theme(plot.title=element_text(size=14)) +
    scale_fill_viridis(na.value="white",option=colorscale, direction = colorscale_dir,limits=scale.limits, begin = begin, end = end) +
    scale_x_continuous(limit=c(-55000,60000))
  
  return(p)
  
}

scenarios = unique(d$scenario)

for(scen in scenarios) {

  d_plot = d %>%
    filter(scenario == scen)
    # remove model cells that didn't run
  
  # Dobrowski model plots
  a <- ggplotGrob(mapfun(d_plot,column = "dob_aet025", colorscale = "viridis", colorscale_dir = -1, title = "c) AET (PET coefficient = 0.25)"))
  b <- ggplotGrob(mapfun(d_plot,column = "dob_aet100", colorscale = "viridis", colorscale_dir = -1, title = "a) AET (PET coefficient = 1.00)"))
  c <- ggplotGrob(mapfun(d_plot,column = "diff_dob_aet", colorscale = "magma", colorscale_dir = 1, title = "e) AET difference"))
  
  e <- ggplotGrob(mapfun(d_plot,column = "dob_cwd025", colorscale = "viridis", colorscale_dir = 1, title = "d) CWD (PET coefficient = 0.25)"))
  f <- ggplotGrob(mapfun(d_plot,column = "dob_cwd100", colorscale = "viridis", colorscale_dir = 1, title = "b) CWD (PET coefficient = 1.00)"))
  g <- ggplotGrob(mapfun(d_plot,column = "diff_dob_cwd", colorscale = "magma", colorscale_dir = 1, title = "f) CWD difference"))
  
  png(paste0("figures/map_wb/dob_",scen,".png"), width = 3300,height = 2000, res=250)
  grid.arrange(b,f,a,e,c,g,ncol=2)
  dev.off()
  
  # Wilmott model plots
  a <- ggplotGrob(mapfun(d_plot,column = "wil_aet025", colorscale = "viridis", colorscale_dir = -1, title = "c) AET (PET coefficient = 0.25)"))
  b <- ggplotGrob(mapfun(d_plot,column = "wil_aet100", colorscale = "viridis", colorscale_dir = -1, title = "a) AET (PET coefficient = 1.00)"))
  c <- ggplotGrob(mapfun(d_plot,column = "diff_wil_aet", colorscale = "magma", colorscale_dir = 1, title = "e) AET difference"))
  
  e <- ggplotGrob(mapfun(d_plot,column = "wil_cwd025", colorscale = "viridis", colorscale_dir = 1, title = "d) CWD (PET coefficient = 0.25)"))
  f <- ggplotGrob(mapfun(d_plot,column = "wil_cwd100", colorscale = "viridis", colorscale_dir = 1, title = "b) CWD (PET coefficient = 1.00)"))
  g <- ggplotGrob(mapfun(d_plot,column = "diff_wil_cwd", colorscale = "magma", colorscale_dir = 1, title = "f) CWD difference"))
  
  png(paste0("figures/map_wb/wil_",scen,".png"), width = 3300,height = 2000, res=250)
  grid.arrange(b,f,a,e,c,g,ncol=2)
  dev.off()
}









### Plot heatmap of differences by elev and rad

heatmapfun <- function(d,columns,scale.limits=NULL, scenario) {
  
  d = d %>%
    select(elev_cat,rad_cat,!!sym(columns[1]),!!sym(columns[2])) %>%
    pivot_longer(cols=c(!!sym(columns[1]),!!sym(columns[2])),
                 names_to = "metric",
                 values_to = "value") %>%
    mutate(metric = ifelse(str_ends(metric,"_aet"),"AET","CWD"))
  
  if(scenario != "") {
    d = d %>%
      mutate(metric = paste0(scenario,": ",metric))
  }
  
  ggplot(d,aes(x=elev_cat,y=rad_cat,fill=value)) +
    geom_tile() +
    scale_fill_viridis(option="magma", begin = 0.1, end = 0.9) +
    labs(x = "Elevation (m)",
         y = "Relative solar exposure") +
    theme_classic(16) +
    theme(legend.title=element_blank()) +
    facet_wrap("metric") +
    theme(strip.text=element_text(size=12),strip.background=element_blank())
}

## bin it myself and require a minimum number of pixels per tile

elev_breaks = seq(min(d$elev),max(d$elev),length.out=6)
elev_breaks = seq(from=0,to=4000,length.out = 6) # manual
elev_mids = elev_breaks[1:(length(elev_breaks)-1)] + 0.5*(elev_breaks[2]-elev_breaks[1])
rad_breaks = seq(min(d$rad.03_scaled),max(d$rad.03_scaled),length.out=6)
rad_mids = rad_breaks[1:(length(rad_breaks)-1)] + 0.5*(rad_breaks[2]-rad_breaks[1])


prep_plot_scenario = function(d,scen) {

  d_plot = d %>%
    filter(scenario == scen)
  
  st_geometry(d_plot) = NULL
  
  d_plot = d_plot %>%
    mutate(elev_cat = cut(elev,elev_breaks,labels=elev_mids),
           rad_cat = cut(rad.03_scaled,rad_breaks,labels=rad_mids)) %>%
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
  
  return(d_plot)
}

d_plot = prep_plot_scenario(d,scen="base")
a <- ggplotGrob(heatmapfun(d_plot,c("diff_dob_aet","diff_dob_cwd"), scenario = "Base"))

d_plot = prep_plot_scenario(d,scen="doubleppt")
b <- ggplotGrob(heatmapfun(d_plot,c("diff_dob_aet","diff_dob_cwd"), scenario = "2x ppt"))

d_plot = prep_plot_scenario(d,scen="constppt")
c <- ggplotGrob(heatmapfun(d_plot,c("diff_dob_aet","diff_dob_cwd"), scenario = "Flat ppt"))

d_plot = prep_plot_scenario(d,scen="consttemp")
e <- ggplotGrob(heatmapfun(d_plot,c("diff_dob_aet","diff_dob_cwd"), scenario = "Flat temp"))

d_plot = prep_plot_scenario(d,scen="consttempppt")
f <- ggplotGrob(heatmapfun(d_plot,c("diff_dob_aet","diff_dob_cwd"), scenario = "Flat temp, flat ppt"))

d_plot = prep_plot_scenario(d,scen="consttemppptdoubleppt")
g <- ggplotGrob(heatmapfun(d_plot,c("diff_dob_aet","diff_dob_cwd"), scenario = "Flat temp, 2x flat ppt"))

png(paste0("figures/heatmap_diff/dob_all.png"), width = 1800,height = 5000, res=250)
grid.arrange(a,b,c,e,f,g,ncol=1)
dev.off()

### Make a base Dobrowski and Wilmott figures

d_plot = prep_plot_scenario(d,scen="base")
g <- ggplotGrob(heatmapfun(d_plot,c("diff_dob_aet","diff_dob_cwd"), scenario = ""))

png(paste0("figures/heatmap_diff/dob_base.png"), width = 1800,height = 1000, res=250)
plot(g)
dev.off()


d_plot = prep_plot_scenario(d,scen="base")
g <- ggplotGrob(heatmapfun(d_plot,c("diff_wil_aet","diff_wil_cwd"), scenario = ""))

png(paste0("figures/heatmap_diff/wil_base.png"), width = 1800,height = 1000, res=250)
plot(g)
dev.off()


