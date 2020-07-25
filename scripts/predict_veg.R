library(raster)
library(sf)
library(tidyverse)
library(viridis)
library(ggthemes)
library(ggspatial)
library(mgcv)

d_rast = brick("data/wb_output/rasters/landscape_wb_inputs_outputs.grd")
cv.type <- raster("data/calveg_clipped/calveg_clipped_raster_whrtype.tif")
cv_crosswalk <- read.csv("data/calveg_clipped/type_conversion_table.csv",header=TRUE,stringsAsFactors=FALSE)
tuol_mask <- st_read("data/tuolomne_mask/tuolomne_grid_mask.shp")

### Get the env data synced to the calveg data (which first needs to be aggregated to the scale of analysis -- 240 m)
set.seed(1)
cv_project = projectRaster(cv.type,to = d_rast, alignOnly=TRUE, method="ngb")
cv_agg = aggregate(cv_project,fact = 240/30, fun=modal)
d_rast = projectRaster(d_rast,cv_agg)

cv = crop(cv_agg,d_rast)
d_rast = crop(d,cv_agg)

d_mod = stack(d_rast,cv)

d_mod = mask(d_mod,tuol_mask %>% st_transform(crs(d_rast)))
d_mod = mask(d_mod,tuol_mask %>% st_transform(crs(d_rast)))

# index_layer = d_mod[[1]]*0
# names(index_layer) = "index_mod"
# values(index_layer) = 1:ncell(index_layer)
# d_mod = stack(d_mod,index_layer)

### alternate approach to bboxes
grid = st_make_grid(tuol_mask %>% st_transform(3310), cellsize=10000,what="corners",square = FALSE)
plots = grid %>% st_buffer(1000) %>% st_as_sf
bboxes = plots
bboxes$training = 1:nrow(bboxes)


## get bbox overlap yes/no
bbox_overlap = rasterize(bboxes,field="training",d_mod[[1]])
d_mod = stack(d_mod,bbox_overlap)
names(d_mod)[length(names(d_mod))] = "training"

# get whether in tuolomne buffer
tuol_mask$tuol_mask = 1
tuol_overlap = rasterize(tuol_mask %>% st_transform(crs(d_mod)),field="tuol_mask",d_mod[[1]])
d_mod = stack(d_mod,tuol_overlap)
names(d_mod)[length(names(d_mod))] = "tuol_overlap"

## Convert raster to sf
d_mod_sf = d_mod %>% as("SpatialPointsDataFrame") %>% as("sf") #     as.data.frame(d_mod)
d_mod_sf$index_df = 1:nrow(d_mod_sf)

d_mod_sf = d_mod_sf %>%
  rename(cv = "calveg_clipped_raster_whrtype") %>%
  filter(!is.na(tuol_overlap))

d = d_mod_sf %>% st_transform(3310)
d$x = st_coordinates(d)[,1]
d$y = st_coordinates(d)[,2]


### Compute annual temp and precip values
d = d %>%
  mutate(tmean_annual = tmean.01+tmean.02+tmean.03+tmean.04+tmean.05+tmean.06+tmean.07+tmean.08+tmean.09+tmean.10+tmean.11+tmean.12) %>%
  mutate(ppt_annual = ppt.01+ppt.02+ppt.03+ppt.04+ppt.05+ppt.06+ppt.07+ppt.08+ppt.09+ppt.10+ppt.11+ppt.12)

### remove vals we don't need, standardize others
d = d %>%
  dplyr::select(rad.03,ID,starts_with("AET"),starts_with("Deficit"),cv,training,tuol_overlap,tmean_annual,ppt_annual, x, y) %>%
  rename(aet_wil100 = AET.PT.STD.Wil150mm,
         aet_wil025 = AET.PT.cc025.Wil150mm,
         aet_dob100 = AET.Dobr.cc100,
         aet_dob025 = AET.Dobr.cc025,
         cwd_wil100 = Deficit.PT.STD.Wil150mm,
         cwd_wil025 = Deficit.PT.cc025.Wil150mm,
         cwd_dob100 = Deficit.Dobr.cc100,
         cwd_dob025 = Deficit.Dobr.cc025,
         aet_tempppt = ppt_annual,
         cwd_tempppt = tmean_annual) %>%
  dplyr::select(-rad.03,-AET.PT.cc050.Wil150mm,-Deficit.PT.cc050.Wil150mm) %>%
  mutate(ID = 1:nrow(d))



### Get Calveg WHR types (text), and keep top 10
d$cv_text = plyr::mapvalues(d$cv,cv_crosswalk$code,as.character(cv_crosswalk$val))
# 
# ## abundance of diff whr types in study area
# t <- table(d$cv_text)
# t2 <- t[order(t,decreasing=TRUE)]
# t2
# plot(t2)
# 
# #remove lake
# types.drop <- "LAC"
# t3 <- t2[!(names(t2) %in% types.drop)]
# 
# #take top 10 types
# top.whrtypes <- names(t3[1:10])

#convert to numeric
top.whrtypes.numeric <- plyr::mapvalues(top.whrtypes,cv_crosswalk$val,cv_crosswalk$code)


## Calc type presences
d = d %>%
  mutate(mhw = cv_text == "MHW",
         mch = cv_text == "MCH",
         smc = cv_text == "SMC")



## Separate out the training dataset
d_train <- d[!is.na(d$training),]
d_train <- d_train[d_train$cv_text %in% top.whrtypes,]






#### Just fit a few veg types
m_dob025_mhw <- gam(mhw ~ s(aet_dob025, k=3) + s(cwd_dob025, k=3), data=d_train, family="binomial")
m_dob100_mhw <- gam(mhw ~ s(aet_dob100, k=3) + s(cwd_dob100, k=3), data=d_train, family="binomial")
m_wil100_mhw <- gam(mhw ~ s(aet_wil100, k=3) + s(cwd_wil100, k=3), data=d_train, family="binomial")
m_wil025_mhw <- gam(mhw ~ s(aet_wil025, k=3) + s(cwd_wil025, k=3), data=d_train, family="binomial")
m_tpp_mhw = gam(mhw ~ s(aet_tempppt, k=3) + s(cwd_tempppt, k=3), data=d_train, family="binomial")

m_dob025_mch <- gam(mch ~ s(aet_dob025, k=3) + s(cwd_dob025, k=3), data=d_train, family="binomial")
m_dob100_mch <- gam(mch ~ s(aet_dob100, k=3) + s(cwd_dob100, k=3), data=d_train, family="binomial")
m_wil100_mch <- gam(mch ~ s(aet_wil100, k=3) + s(cwd_wil100, k=3), data=d_train, family="binomial")
m_wil025_mch <- gam(mch ~ s(aet_wil025, k=3) + s(cwd_wil025, k=3), data=d_train, family="binomial")
m_tpp_mch = gam(mch ~ s(aet_tempppt, k=3) + s(cwd_tempppt, k=3), data=d_train, family="binomial")

m_dob025_smc <- gam(smc ~ s(aet_dob025, k=3) + s(cwd_dob025, k=3), data=d_train, family="binomial")
m_dob100_smc <- gam(smc ~ s(aet_dob100, k=3) + s(cwd_dob100, k=3), data=d_train, family="binomial")
m_wil100_smc <- gam(smc ~ s(aet_wil100, k=3) + s(cwd_wil100, k=3), data=d_train, family="binomial")
m_wil025_smc <- gam(smc ~ s(aet_wil025, k=3) + s(cwd_wil025, k=3), data=d_train, family="binomial")
m_tpp_smc = gam(smc ~ s(aet_tempppt, k=3) + s(cwd_tempppt, k=3), data=d_train, family="binomial")





## predict to remaining dataset
d$p_dob025_mhw = predict(m_dob025_mhw, newdata=d, type="response")
d$p_dob100_mhw = predict(m_dob100_mhw, newdata=d, type = "response")
d$p_wil025_mhw = predict(m_wil025_mhw, newdata=d, type="response")
d$p_wil100_mhw = predict(m_wil100_mhw, newdata=d, type = "response")
d$p_tempppt_mhw = predict(m_tpp_mhw, newdata=d, type = "response")

d$p_dob025_mch = predict(m_dob025_mch, newdata=d, type="response")
d$p_dob100_mch = predict(m_dob100_mch, newdata=d, type = "response")
d$p_wil025_mch = predict(m_wil025_mch, newdata=d, type="response")
d$p_wil100_mch = predict(m_wil100_mch, newdata=d, type = "response")
d$p_tempppt_mch = predict(m_tpp_mch, newdata=d, type = "response")

d$p_dob025_smc = predict(m_dob025_smc, newdata=d, type="response")
d$p_dob100_smc = predict(m_dob100_smc, newdata=d, type = "response")
d$p_wil025_smc = predict(m_wil025_smc, newdata=d, type="response")
d$p_wil100_smc = predict(m_wil100_smc, newdata=d, type = "response")
d$p_tempppt_smc = predict(m_tpp_smc, newdata=d, type = "response")


# what prop of pixels are mhw?
mhw_prop = sum(d$mhw,na.rm=TRUE)/sum(!is.na(d$mhw))
# what predicted prob would give us that proportion of pixels?
mhw_thresh_dob025 = quantile(d$p_dob025_mhw,1-mhw_prop)
mhw_thresh_dob100 = quantile(d$p_dob100_mhw,1-mhw_prop)

# what prop of pixels are smc?
smc_prop = sum(d$smc,na.rm=TRUE)/sum(!is.na(d$smc))
# what predicted prob would give us that proportion of pixels?
smc_thresh_dob025 = quantile(d$p_dob025_smc,1-smc_prop)
smc_thresh_dob100 = quantile(d$p_dob100_smc,1-smc_prop)

# what prop of pixels are mch?
mch_prop = sum(d$mch,na.rm=TRUE)/sum(!is.na(d$mch))
# what predicted prob would give us that proportion of pixels?
mch_thresh_dob025 = quantile(d$p_dob025_mch,1-mch_prop)
mch_thresh_dob100 = quantile(d$p_dob100_mch,1-mch_prop)

d = d %>%
  mutate(p_dob025_mhw_pres = p_dob025_mhw > mhw_thresh_dob025) %>%
  mutate(p_dob100_mhw_pres = p_dob100_mhw > mhw_thresh_dob100) %>%
  mutate(p_dob025_mch_pres = p_dob025_mch > mch_thresh_dob025) %>%
  mutate(p_dob100_mch_pres = p_dob100_mch > mch_thresh_dob100) %>%
  mutate(p_dob025_smc_pres = p_dob025_smc > smc_thresh_dob025) %>%
  mutate(p_dob100_smc_pres = p_dob100_smc > smc_thresh_dob100)
  












#### Plotting ####

ggplot(d) +
  geom_tile(aes(fill=p_dob100_mch_pres, x=x,y=y)) +
  #coord_equal() +
  theme_map() +
  geom_sf(data=bboxes,fill=NA,color="red",size=1) +
  #scale_fill_viridis(na.value="white",discrete=TRUE,option="inferno",name="Vegetation type") +
  theme(legend.position="bottom") +
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
  theme(strip.text=element_text(size=12),strip.background=element_blank()) +
  theme(panel.spacing=unit(1,"lines")) +
  theme(panel.border=element_rect(fill=NA)) +
  theme(legend.text=element_text(size=9),legend.title=element_text(size=11))
  


## plot of the AET and CWD space of the training units and in general

d_plot = d %>%
  filter(!is.na(tuol_overlap)) %>%
  mutate(training_bool = ifelse(!is.na(training),"Yes","No")) %>%
  arrange(training_bool)

ggplot(d, aes(x = AET.Dobr.cc100, y = Deficit.Dobr.cc100, color=training)) +
  geom_point(size=0.5)




#### Evaluate/explain changes in veg predictions ####

## for each CWB method and veg type, calculate 025, both, and 100

d_eval = d
st_geometry(d_eval) = NULL

d_eval_long = d_eval %>%
  pivot_longer(cols=c(starts_with("aet_"),starts_with("cwd_"), starts_with("p_"), )) %>%
  # get the second part of the name (it's the climate method)
  mutate(clim_metric = str_split(name,"_") %>% map(2)) %>%
  mutate(valtype_veg = str_c(str_split(name,"_") %>% map(1),str_split(name,"_") %>% map(3),str_split(name,"_") %>% map(4), sep= "_")) %>%
  mutate(valtype_veg = str_replace(valtype_veg,"_NULL","")) %>%
  mutate(valtype_veg = str_replace(valtype_veg,"_NULL","")) %>%
  #arrange(ID, name, valtype_veg)
  select(-name) %>%
  pivot_wider(names_from = valtype_veg, values_from = value) %>%
  rename(obs_mhw = mhw,
         obs_smc = smc,
         obs_mch = mch) %>%
  pivot_longer(cols = c(starts_with("p_"), starts_with("obs_"))) %>%
  mutate(presence_metric = case_when(str_ends(name,"_pres") ~ "presab",
                                     str_starts(name,"obs_") ~ "obs",
                                      TRUE ~ "prob")) %>%
  mutate(vegtype = str_split(name,"_") %>% map(2)) %>%
  select(-name) %>%
  pivot_wider(names_from = presence_metric, values_from = value)

## Put the predictions resulting from different CC assumptions side by side
d_eval_presab = d_eval_long %>%
  filter(clim_metric != "tempppt") %>%
  mutate(cc = str_sub(clim_metric,-3,-1)) %>%
  mutate(clim_metric = str_sub(clim_metric,1,3)) %>%
  select(-aet,-cwd) %>%
  pivot_wider(names_from = cc, values_from = c(prob,presab))

## Compute a multi-level venn change classification

d_eval_presab = d_eval_presab %>%
  mutate(presab_summ = case_when(presab_100 & !presab_025 ~ "high",
                                 presab_025 & !presab_100 ~ "low",
                                 presab_100 & presab_025 ~ "both",
                                 TRUE ~ "neither"))



## quick test plot

d_plot = d_eval_presab %>%
  filter(clim_metric == "dob",
         vegtype == "mhw")


ggplot(d_plot,aes(x=x,y=y,fill=presab_summ)) +
  geom_tile() +
  coord_equal()













npixels = d %>%
  nrow()

npixels_outside_training = d %>%
  filter(is.na(training)) %>%
  nrow()

npixels_outside_training_changed_dob = d %>%
  filter(is.na(training)) %>%
  filter(dob_change == TRUE) %>%
  nrow()

npixels_outside_training_changed_wil = d %>%
  filter(!is.na(p_wil025_name)) %>% # this is redundant
  filter(is.na(training)) %>%
  filter(wil_change == TRUE) %>%
  nrow()

npixels_outside_training_changed_dob / npixels_outside_training

npixels_outside_training_changed_wil / npixels_outside_training

