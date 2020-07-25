library(raster)
library(sf)
library(tidyverse)

d = brick("data/wb_output/rasters/landscape_wb_inputs_outputs.grd")
cv.type <- raster("data/calveg_clipped/calveg_clipped_raster_whrtype.tif")
cv_crosswalk <- read.csv("data/calveg_clipped/type_conversion_table.csv",header=TRUE,stringsAsFactors=FALSE)
tuol_mask <- st_read("data/tuolomne_mask/tuolomne_grid_mask.shp")


### Get the env data synced to the calveg data (which first needs to be aggregated to the scale of analysis -- 240 m)
set.seed(1)
cv_project = projectRaster(cv.type,to = d, alignOnly=TRUE, method="ngb")
cv_agg = aggregate(cv_project,fact = 240/30, fun=modal)
d = projectRaster(d,cv_agg)

cv = crop(cv_agg,d)
d = crop(d,cv_agg)

d_mod = stack(d,cv)

d_mod = mask(d_mod,tuol_mask %>% st_transform(crs(d)))
d_mod = mask(d_mod,tuol_mask %>% st_transform(crs(d)))


index_layer = d_mod[[1]]*0
names(index_layer) = "index_mod"
values(index_layer) = 1:ncell(index_layer)
d_mod = stack(d_mod,index_layer)

### Create sampling points
bbox_fun <- function(x1,y1,size) { # give it the top-left corner
  x2 <- x1+size
  y2 <- y1-size
  bbox <- st_polygon(list(cbind(c(x1,x2,x2,x1,x1),c(y1,y1,y2,y2,y1))))
}

bbox1 <- bbox_fun(-32500,-12500,10000)
bbox2 <- bbox_fun(2500,-10000,10000)
bbox3 <- bbox_fun(36000,-4500,10000)
bboxes <- st_sfc(bbox1,bbox2,bbox3,crs=3310) %>% st_as_sf


### alternate approach to bboxes
grid = st_make_grid(tuol_mask %>% st_transform(3310), cellsize=10000,what="corners",square = FALSE)
plots = grid %>% st_buffer(1000) %>% st_as_sf
bboxes = plots


# ## second four instead of 3
# bbox1 <- bbox_fun(-34500,-15000,10000)
# bbox2 <- bbox_fun(-9000,-13000,10000)
# bbox3 <- bbox_fun(20000,-5000,10000)
# bbox4 <- bbox_fun(46000,-1500,10000)
# bboxes <- st_sfc(bbox1,bbox2,bbox3,bbox4,crs=3310) %>% st_sf()


bboxes$training = 1:nrow(bboxes)
tuol_mask$tuol_mask = 1




## get bbox overlap yes/no
bbox_overlap = rasterize(bboxes,field="training",d_mod[[1]])
d_mod = stack(d_mod,bbox_overlap)
names(d_mod)[length(names(d_mod))] = "training"

# get whether in tuolomne buffer
tuol_overlap = rasterize(tuol_mask %>% st_transform(crs(d_mod)),field="tuol_mask",d_mod[[1]])
d_mod = stack(d_mod,tuol_overlap)
names(d_mod)[length(names(d_mod))] = "tuol_overlap"


## Convert raster to sf
d_mod_sf = d_mod %>% as("SpatialPointsDataFrame") %>% as("sf") #     as.data.frame(d_mod)
d_mod_sf$index_df = 1:nrow(d_mod_sf)

d_mod_sf = d_mod_sf %>%
  rename(cv = "calveg_clipped_raster_whrtype")
  filter(!is.na(AET.PT.STD.Wil150mm) & !is.na(cv))

d = d_mod_sf

### Compute annual temp and precip values
d = d %>%
  mutate(tmean_annual = tmean.01+tmean.02+tmean.03+tmean.04+tmean.05+tmean.06+tmean.07+tmean.08+tmean.09+tmean.10+tmean.11+tmean.12) %>%
  mutate(ppt_annual = ppt.01+ppt.02+ppt.03+ppt.04+ppt.05+ppt.06+ppt.07+ppt.08+ppt.09+ppt.10+ppt.11+ppt.12)


### Get Calveg WHR types, and keep top 10
d$cv_text = plyr::mapvalues(d$cv,cv_crosswalk$code,as.character(cv_crosswalk$val))

## abundance of diff whr types in study area
t <- table(d$cv_text)
t2 <- t[order(t,decreasing=TRUE)]
t2
plot(t2)

#remove lake
types.drop <- "LAC"
t3 <- t2[!(names(t2) %in% types.drop)]

#take top 10 types
top.whrtypes <- names(t3[1:10])

#convert to numeric
top.whrtypes.numeric <- plyr::mapvalues(top.whrtypes,cv_crosswalk$val,cv_crosswalk$code)

## thin to 1000 training points  ##!! need to remove?
d_train <- d[!is.na(d$training),]
d_train <- d_train[d_train$cv_text %in% top.whrtypes,]
thin_to <- 1000
full <- nrow(d_train)
thin_factor <- round(full/thin_to)
d_train <- d_train[seq(1,nrow(d_train),by=thin_factor),]


#### Fit multinomial regression models ####
library(nnet)

m.1 <- nnet::multinom(cv_text~(AET.Dobr.cc025+Deficit.Dobr.cc025),data=d_train,maxit=10000)
m_dob025 <- nnet::multinom(cv_text~(scale(AET.Dobr.cc025)+scale(Deficit.Dobr.cc025)+I(scale(AET.Dobr.cc025)^2)+I(scale(Deficit.Dobr.cc025)^2)),data=d_train,maxit=10000)
m_dob100 <- nnet::multinom(cv_text~(scale(AET.Dobr.cc100)+scale(Deficit.Dobr.cc100)+I(scale(AET.Dobr.cc100)^2)+I(scale(Deficit.Dobr.cc100)^2)),data=d_train,maxit=10000)
m_wil025 <- nnet::multinom(cv_text~(scale(AET.PT.cc025.Wil150mm)+scale(Deficit.PT.cc025.Wil150mm)+I(scale(AET.PT.cc025.Wil150mm)^2)+I(scale(Deficit.PT.cc025.Wil150mm)^2)),data=d_train,maxit=10000)
m_wil100 <- nnet::multinom(cv_text~(scale(AET.PT.STD.Wil150mm)+scale(Deficit.PT.STD.Wil150mm)+I(scale(AET.PT.STD.Wil150mm)^2)+I(scale(Deficit.PT.STD.Wil150mm)^2)),data=d_train,maxit=10000)
m_tpp <- nnet::multinom(cv_text~((scale(tmean_annual))+(scale(ppt_annual))),data=d_train,maxit=10000)
m_ttp <- nnet::multinom(cv_text~(scale(tmean_annual)*scale(ppt_annual)),data=d_train,maxit=10000)
m_tpp2 <- nnet::multinom(cv_text~scale(tmean_annual)+scale(ppt_annual) + I(scale(tmean_annual)^2)+ I(scale(ppt_annual)^2),data=d_train,maxit=10000)
m_ttp2 <- nnet::multinom(cv_text~(scale(tmean_annual)*scale(ppt_annual))+I(scale(tmean_annual)^2)*I(scale(ppt_annual)^2),data=d_train,maxit=10000)


#### Just predict montane hardwood
m_dob025 <- nnet::multinom(cv_text~(scale(AET.Dobr.cc025)+scale(Deficit.Dobr.cc025)+I(scale(AET.Dobr.cc025)^2)+I(scale(Deficit.Dobr.cc025)^2)),data=d_train,maxit=10000)






# d_train$p_dob025 = predict(m_dob025)
# d_train$p_dob100 = predict(m_dob100)
# d_train$p_wil025 = predict(m_wil025)
# d_train$p_wil100 = predict(m_wil100)

## predict to remaining dataset
d$p_dob025 = predict(m_dob025,newdat=d)
d$p_dob100 = predict(m_dob100,newdat=d)
d$p_wil025 = predict(m_wil025,newdat=d)
d$p_wil100 = predict(m_wil100,newdat=d)
d$p_tpp = predict(m_tpp,newdat=d)
d$p_ttp = predict(m_ttp,newdat=d)
d$p_tpp2 = predict(m_tpp2,newdat=d)
d$p_ttp2 = predict(m_ttp2,newdat=d)

## whrtype back to numeric code
d$p_dob025_num <- plyr::mapvalues(d$p_dob025,cv_crosswalk$val,as.numeric(cv_crosswalk$code)) %>% as.character %>% as.numeric
d$p_dob100_num <- plyr::mapvalues(d$p_dob100,cv_crosswalk$val,as.numeric(cv_crosswalk$code)) %>% as.character %>% as.numeric
d$p_wil025_num <- plyr::mapvalues(d$p_wil025,cv_crosswalk$val,as.numeric(cv_crosswalk$code)) %>% as.character %>% as.numeric
d$p_wil100_num <- plyr::mapvalues(d$p_wil100,cv_crosswalk$val,as.numeric(cv_crosswalk$code)) %>% as.character %>% as.numeric
d$p_tpp_num <- plyr::mapvalues(d$p_tpp,cv_crosswalk$val,as.numeric(cv_crosswalk$code)) %>% as.character %>% as.numeric
d$p_ttp_num <- plyr::mapvalues(d$p_ttp,cv_crosswalk$val,as.numeric(cv_crosswalk$code)) %>% as.character %>% as.numeric
d$p_tpp2_num <- plyr::mapvalues(d$p_tpp2,cv_crosswalk$val,as.numeric(cv_crosswalk$code)) %>% as.character %>% as.numeric
d$p_ttp2_num <- plyr::mapvalues(d$p_ttp2,cv_crosswalk$val,as.numeric(cv_crosswalk$code)) %>% as.character %>% as.numeric


## to spatial points for plotting
# pred_df = as(d_mod,"SpatialPointsDataFrame") %>% as("sf")
# 
# pred_df[d$index,"dob025"] <- d$p_dob025
# pred_df[d$index,"dob100"] <- d$p_dob100
# pred_df[d$index,"wil025"] <- d$p_dob025
# pred_df[d$index,"wil100"] <- d$p_dob100
# 
# pred_df[d$index,"observed"] <- d$cv_text
# pred_df[d$index,"observed"] = ifelse((pred_df[d$index,]$observed %in% top.whrtypes),  pred_df[d$index,]$observed, NA)
# 
# pred_df$x = st_coordinates(pred_df)[,1]
# pred_df$y = st_coordinates(pred_df)[,2]

## turn d into a raster (for purposes of syncing everything)

p_dob025_num = d_mod[[1]] * 0
values(p_dob025_num) = d$p_dob025_num
names(p_dob025_num) = "p_dob025"

p_dob100_num = d_mod[[1]] * 0
values(p_dob100_num) = d$p_dob100_num
names(p_dob100_num) = "p_dob100"

p_wil025_num = d_mod[[1]] * 0
values(p_wil025_num) = d$p_wil025_num
names(p_wil025_num) = "p_wil025"

p_wil100_num = d_mod[[1]] * 0
values(p_wil100_num) = d$p_wil100_num
names(p_wil100_num) = "p_wil100"

p_tpp_num = d_mod[[1]] * 0
values(p_tpp_num) = d$p_tpp_num
names(p_tpp_num) = "p_tpp"

p_ttp_num = d_mod[[1]] * 0
values(p_ttp_num) = d$p_ttp_num
names(p_ttp_num) = "p_ttp"

p_tpp2_num = d_mod[[1]] * 0
values(p_tpp2_num) = d$p_tpp2_num
names(p_tpp2_num) = "p_tpp2"

p_ttp2_num = d_mod[[1]] * 0
values(p_ttp2_num) = d$p_ttp2_num
names(p_ttp2_num) = "p_ttp2"

d_fullstack = stack(d_mod,p_dob025_num,p_dob100_num,p_wil025_num,p_wil100_num, p_tpp_num, p_ttp_num, p_ttp2_num, p_tpp2_num)

### Back to sf points for plotting

d_fullstack_sf = d_fullstack %>% as("SpatialPointsDataFrame") %>% as("sf")
d_fullstack_sf$x = st_coordinates(d_fullstack_sf)[,1]
d_fullstack_sf$y = st_coordinates(d_fullstack_sf)[,2]

d_fullstack_sf = d_fullstack_sf %>%
  mutate_at(vars(starts_with("p_")),  funs(name=plyr::mapvalues(.,cv_crosswalk$code,cv_crosswalk$val,) )      ) %>%
  mutate(calveg_clipped_raster_whrtype_name = plyr::mapvalues(calveg_clipped_raster_whrtype,cv_crosswalk$code,cv_crosswalk$val)) %>%
  mutate(calveg_clipped_raster_whrtype_name_top = ifelse(calveg_clipped_raster_whrtype_name %in% top.whrtypes, calveg_clipped_raster_whrtype_name, NA)) %>%
  filter(!is.na(p_wil100))

train_region = st_sf(bboxes)

library(viridis)
library(ggthemes)
library(ggspatial)
ggplot(d_fullstack_sf) +
  geom_tile(aes(fill=p_dob100_name, x=x,y=y)) +
  #coord_equal() +
  theme_map() +
  geom_sf(data=train_region,fill=NA,color="red",size=1) +
  scale_fill_viridis(na.value="white",discrete=TRUE,option="inferno",name="Vegetation type") +
  theme(legend.position="bottom") +
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
  theme(strip.text=element_text(size=12),strip.background=element_blank()) +
  theme(panel.spacing=unit(1,"lines")) +
  theme(panel.border=element_rect(fill=NA)) +
  theme(legend.text=element_text(size=9),legend.title=element_text(size=11))
  


## plot of the AET and CWD space of the training units and in general

d_plot = d_fullstack_sf %>%
  filter(!is.na(tuol_overlap)) %>%
  mutate(training_bool = ifelse(!is.na(training),"Yes","No")) %>%
  arrange(training_bool)

ggplot(d_plot, aes(x = AET.Dobr.cc100, y = Deficit.Dobr.cc100, color=training)) +
  geom_point(size=0.1)




#### Evaluate/explain changes in veg predictions ####

d_fullstack_sf = d_fullstack_sf %>%
  mutate(dob_change = p_dob025 != p_dob100,
         wil_change = p_wil025 != p_wil100)

npixels = d_fullstack_sf %>%
  filter(!is.na(p_dob025_name)) %>% # this is redundant
  nrow()

npixels_outside_training = d_fullstack_sf %>%
  filter(!is.na(p_dob025_name)) %>% # this is redundant
  filter(is.na(training)) %>%
  nrow()

npixels_outside_training_changed_dob = d_fullstack_sf %>%
  filter(!is.na(p_dob025_name)) %>% # this is redundant
  filter(is.na(training)) %>%
  filter(dob_change == TRUE) %>%
  nrow()

npixels_outside_training_changed_wil = d_fullstack_sf %>%
  filter(!is.na(p_wil025_name)) %>% # this is redundant
  filter(is.na(training)) %>%
  filter(wil_change == TRUE) %>%
  nrow()

npixels_outside_training_changed_dob / npixels_outside_training

npixels_outside_training_changed_wil / npixels_outside_training

