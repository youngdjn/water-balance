# Use computed water balance values to predict vegetation distributions across the study landscape, and compute model fit/validation and agreement for 1000 iterations of different training data partitions
# Note that throughout this script, variables and columns that end in "025" and "100" refer to PET coefficient 0.25 and 1.00

### Specify how many cross-validation iterations to run for. If 1, no cross-validation, just a single fit. Need to set to 1 to make Fig. 6
cross_validation_iterations = 1


library(raster)
library(sf)
library(tidyverse)
library(viridis)
library(ggthemes)
library(ggspatial)
library(mgcv)
library(pROC)
library(irr)

## Convenience functions
#true skill statistic
tss <- function(pred, obs){
  tpr = sum(pred & obs) / (sum(pred&obs) + sum(!pred&obs))
  tnr = sum(!pred & !obs) / (sum(!pred & !obs) + sum(pred&!obs))
  return(tpr+tnr-1)
}

# Load base scenario water balance
d_rast = brick("data/wb_output/rasters/base.grd")

# Load calveg eveg veg distribution data (in this script cv refers to "calveg")
cv.type <- raster("data/calveg_clipped/calveg_clipped_raster_whrtype.tif")

# Crosswalk to go from Calveg numbers to WHR veg types
cv_crosswalk <- read.csv("data/calveg_clipped/type_conversion_table.csv",header=TRUE,stringsAsFactors=FALSE)

# Study area mask
tuol_mask <- st_read("data/tuolomne_mask/tuolomne_grid_mask.shp")


##### Get the env data synced to the calveg data (which first needs to be aggregated to the scale of analysis -- 240 m) #####
set.seed(1)
cv_project = projectRaster(cv.type,to = d_rast, alignOnly=TRUE, method="ngb")
cv_agg = aggregate(cv_project,fact = 240/30, fun=modal)
d_rast = projectRaster(d_rast,cv_project)
d_rast = aggregate(d_rast,fact = 240/30, fun=mean)

cv = crop(cv_agg,d_rast)
d_rast = crop(d_rast,cv_agg)

d_mod = stack(d_rast,cv)

d_mod = mask(d_mod,tuol_mask %>% st_transform(crs(d_rast)))
d_mod = mask(d_mod,tuol_mask %>% st_transform(crs(d_rast)))




### Prepare to do 1000 iterations of fitting distribution models (each with a different offset to the training data locations)

# rename the data so we can modify d_mod during each iteration but then retrieve the original data frame at the start of the next iteration
d_mod_orig = d_mod

# data frame to store model fit stats
veg_fits_main_overall = data.frame()
d_eval_presab = data.frame()

for(k in 1:cross_validation_iterations) { # takes ~8 hours on a recent laptop

            # make a grid of training data locations (call them `bboxes` due to a vestige of an older approach)
            if(cross_validation_iterations == 1) { # not doing cross-validation, so no need to randomly offset training points
              grid = st_make_grid(tuol_mask %>% st_transform(3310), cellsize=10000,what="corners",square = FALSE)
            } else {
              # random offsets to use for the grid of training data locations. a different set of offsets for each iteration
              offsets = runif(2, min=0,max=10000)
              grid = st_make_grid(tuol_mask %>% st_transform(3310), cellsize=10000,what="corners",square = FALSE, offset = offsets)
            }
  
            plots = grid %>% st_buffer(2000) %>% st_as_sf
            bboxes = plots
            
            #clip them to the focal area
            bboxes = st_intersection(bboxes,tuol_mask %>% st_transform(st_crs(bboxes)))
            
            bboxes$training = 1:nrow(bboxes)
            
            
            ## get whether raster cells overlap the trainin data locations
            bbox_overlap = rasterize(bboxes,field="training",d_mod[[1]])
            d_mod = stack(d_mod_orig,bbox_overlap)
            names(d_mod)[length(names(d_mod))] = "training"
            
            # get whether raster cells are in tuolomne buffer
            tuol_mask$tuol_mask = 1
            tuol_overlap = rasterize(tuol_mask %>% st_transform(crs(d_mod)),field="tuol_mask",d_mod[[1]])
            d_mod = stack(d_mod,tuol_overlap)
            names(d_mod)[length(names(d_mod))] = "tuol_overlap"
            
            ## Convert raster to sf points, only points within study area
            d_mod_sf = d_mod %>% as("SpatialPointsDataFrame") %>% as("sf") #     as.data.frame(d_mod)
            d_mod_sf$index_df = 1:nrow(d_mod_sf)
            
            d_mod_sf = d_mod_sf %>%
              rename(cv = "calveg_clipped_raster_whrtype") %>%
              filter(!is.na(tuol_overlap))
            
            # get coords
            d = d_mod_sf %>% st_transform(3310)
            d$x = st_coordinates(d)[,1]
            d$y = st_coordinates(d)[,2]
            
            
            ### Compute annual temp and precip values (for alternative model that uses temp and precip predictors)
            d = d %>%
              mutate(tmean_annual = tmean.01+tmean.02+tmean.03+tmean.04+tmean.05+tmean.06+tmean.07+tmean.08+tmean.09+tmean.10+tmean.11+tmean.12) %>%
              mutate(ppt_annual = ppt.01+ppt.02+ppt.03+ppt.04+ppt.05+ppt.06+ppt.07+ppt.08+ppt.09+ppt.10+ppt.11+ppt.12)
            
            ### remove vals we don't need, rename the ones we keep. 'cv' refers to calveg eveg veg type data
            d = d %>%
              dplyr::select(rad.03,ID,starts_with("AET"),starts_with("Deficit"),cv,training,tuol_overlap,tmean_annual,ppt_annual, x, y) %>%
              rename(aet_dob100 = AET.Dobr.cc100,
                     aet_dob025 = AET.Dobr.cc025,
                     cwd_dob100 = Deficit.Dobr.cc100,
                     cwd_dob025 = Deficit.Dobr.cc025,
                     aet_tempppt = ppt_annual,
                     cwd_tempppt = tmean_annual) %>%
              mutate(ID = 1:nrow(d))

            
            ### Get Calveg WHR types (text description) from the crosswalk
            d$cv_text = plyr::mapvalues(d$cv,cv_crosswalk$code,as.character(cv_crosswalk$val))
            
            
            ### Calc type presences T/F
            d = d %>%
              mutate(mhw = cv_text == "MHW",
                     mch = cv_text == "MCH",
                     smc = cv_text == "SMC")
            
      
            ## Separate out the training dataset
            d_train <- d[!is.na(d$training),]

            #### Fit distrib models for veg types for the two PET coefficient assumptions, and for temperature-precip, and for temp-ppt-rad and interactions
            m_dob025_mhw <- gam(mhw ~ s((aet_dob025), k=3) + s((cwd_dob025), k=3), data=d_train, family="binomial")
            m_dob100_mhw <- gam(mhw ~ s((aet_dob100), k=3) + s((cwd_dob100), k=3), data=d_train, family="binomial")
            m_tpp_mhw = gam(mhw ~ s((aet_tempppt), k=3) + s((cwd_tempppt), k=3), data=d_train, family="binomial")
            m_tppr_mhw = gam(mhw ~ te((aet_tempppt),(cwd_tempppt),(rad.03), k=3), data=d_train, family="binomial")
            
            m_dob025_mch <- gam(mch ~ s(aet_dob025, k=3) + s(cwd_dob025, k=3), data=d_train, family="binomial")
            m_dob100_mch <- gam(mch ~ s(aet_dob100, k=3) + s(cwd_dob100, k=3), data=d_train, family="binomial")
            m_tpp_mch = gam(mch ~ s(aet_tempppt, k=3) + s(cwd_tempppt, k=3), data=d_train, family="binomial")
            m_tppr_mch = gam(mch ~ te(aet_tempppt,cwd_tempppt,rad.03, k=3), data=d_train, family="binomial")
            
            m_dob025_smc <- gam(smc ~ s(aet_dob025, k=3) + s(cwd_dob025, k=3), data=d_train, family="binomial")
            m_dob100_smc <- gam(smc ~ s(aet_dob100, k=3) + s(cwd_dob100, k=3), data=d_train, family="binomial")
            m_tpp_smc = gam(smc ~ s(aet_tempppt, k=3) + s(cwd_tempppt, k=3), data=d_train, family="binomial")
            m_tppr_smc = gam(smc ~ te(aet_tempppt,cwd_tempppt,rad.03, k=3), data=d_train, family="binomial")
            

            ### predict to landscape
            d$p_dob025_mhw = predict(m_dob025_mhw, newdata=d, type="response")
            d$p_dob100_mhw = predict(m_dob100_mhw, newdata=d, type = "response")
            d$p_tempppt_mhw = predict(m_tpp_mhw, newdata=d, type = "response")
            d$p_temppptr_mhw = predict(m_tppr_mhw, newdata=d, type = "response")
            
            d$p_dob025_mch = predict(m_dob025_mch, newdata=d, type="response")
            d$p_dob100_mch = predict(m_dob100_mch, newdata=d, type = "response")
            d$p_tempppt_mch = predict(m_tpp_mch, newdata=d, type = "response")
            d$p_temppptr_mch = predict(m_tppr_mch, newdata=d, type = "response")
            
            d$p_dob025_smc = predict(m_dob025_smc, newdata=d, type="response")
            d$p_dob100_smc = predict(m_dob100_smc, newdata=d, type = "response")
            d$p_tempppt_smc = predict(m_tpp_smc, newdata=d, type = "response")
            d$p_temppptr_smc = predict(m_tppr_smc, newdata=d, type = "response")
            
            ### For each veg type, calc a predicted pres/ab such that the proportion of predicted presence equals the observed proportion of presence
            
            # what prop of pixels are mhw?
            mhw_prop = sum(d$mhw,na.rm=TRUE)/sum(!is.na(d$mhw))
            # what predicted prob would give us that proportion of pixels?
            mhw_thresh_dob025 = quantile(d$p_dob025_mhw,1-mhw_prop)
            mhw_thresh_dob100 = quantile(d$p_dob100_mhw,1-mhw_prop)
            mhw_thresh_tempppt = quantile(d$p_tempppt_mhw,1-mhw_prop)
            mhw_thresh_temppptr = quantile(d$p_temppptr_mhw,1-mhw_prop)
            
            # what prop of pixels are smc?
            smc_prop = sum(d$smc,na.rm=TRUE)/sum(!is.na(d$smc))
            # what predicted prob would give us that proportion of pixels?
            smc_thresh_dob025 = quantile(d$p_dob025_smc,1-smc_prop)
            smc_thresh_dob100 = quantile(d$p_dob100_smc,1-smc_prop)
            smc_thresh_tempppt = quantile(d$p_tempppt_smc,1-smc_prop)
            smc_thresh_temppptr = quantile(d$p_temppptr_smc,1-smc_prop)
            
            # what prop of pixels are mch?
            mch_prop = sum(d$mch,na.rm=TRUE)/sum(!is.na(d$mch))
            # what predicted prob would give us that proportion of pixels?
            mch_thresh_dob025 = quantile(d$p_dob025_mch,1-mch_prop)
            mch_thresh_dob100 = quantile(d$p_dob100_mch,1-mch_prop)
            mch_thresh_tempppt = quantile(d$p_tempppt_mch,1-mch_prop)
            mch_thresh_temppptr = quantile(d$p_temppptr_mch,1-mch_prop)
            
            ## Optionally manually override the thresholds for a threshold sensitivity analysis
            # mhw_thresh_dob025 = mhw_thresh_dob100 = mhw_thresh_tempppt = mhw_thresh_temppptr = smc_thresh_dob025 = smc_thresh_dob100 = smc_thresh_tempppt =
            #   smc_thresh_temppptr = mch_thresh_dob025 = mch_thresh_dob100 = mch_thresh_tempppt = mch_thresh_temppptr = .30
            
            d = d %>%
              mutate(p_dob025_mhw_pres = p_dob025_mhw > mhw_thresh_dob025) %>%
              mutate(p_dob100_mhw_pres = p_dob100_mhw > mhw_thresh_dob100) %>%
              mutate(p_dob025_mch_pres = p_dob025_mch > mch_thresh_dob025) %>%
              mutate(p_dob100_mch_pres = p_dob100_mch > mch_thresh_dob100) %>%
              mutate(p_dob025_smc_pres = p_dob025_smc > smc_thresh_dob025) %>%
              mutate(p_dob100_smc_pres = p_dob100_smc > smc_thresh_dob100) %>%
              mutate(p_tempppt_smc_pres = p_tempppt_smc > smc_thresh_tempppt) %>%
              mutate(p_tempppt_mhw_pres = p_tempppt_mhw > mhw_thresh_tempppt) %>%
              mutate(p_tempppt_mch_pres = p_tempppt_mch > mch_thresh_tempppt) %>%
              mutate(p_temppptr_smc_pres = p_temppptr_smc > smc_thresh_temppptr) %>%
              mutate(p_temppptr_mhw_pres = p_temppptr_mhw > mhw_thresh_temppptr) %>%
              mutate(p_temppptr_mch_pres = p_temppptr_mch > mch_thresh_temppptr)

            # #### Plotting of the climate space sampled (Fig. S2) ####
            # 
            # ggplot(d) +
            #   geom_tile(aes(fill=p_dob100_mhw_pres, x=x,y=y)) +
            #   #coord_equal() +
            #   theme_map() +
            #   geom_sf(data=bboxes,fill=NA,color="red",size=1) +
            #   #scale_fill_viridis(na.value="white",discrete=TRUE,option="inferno",name="Vegetation type") +
            #   theme(legend.position="bottom") +
            #   theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
            #   theme(strip.text=element_text(size=12),strip.background=element_blank()) +
            #   theme(panel.spacing=unit(1,"lines")) +
            #   theme(panel.border=element_rect(fill=NA)) +
            #   theme(legend.text=element_text(size=9),legend.title=element_text(size=11))
            #   
            # 
            # 
            # ## plot of the AET and CWD space of the training units and in general
            # 
            # d_plot = d %>%
            #   filter(!is.na(tuol_overlap)) %>%
            #   mutate(training_bool = ifelse(!is.na(training),"Training","Validation")) %>%
            #   mutate(training_bool = factor(training_bool, levels = c("Validation","Training"))) %>%
            #   arrange(training_bool)
            # 
            # p = ggplot(d_plot, aes(x = aet_dob100, y = cwd_dob100, color=training_bool)) +
            #   geom_point(size=0.5) +
            #   scale_color_manual(values = c("blue","cornflowerblue"),name = "Pixel type") +
            #   labs(x = "AET (mm)", y = "CWD (mm)") +
            #   theme_classic(16) +
            #   guides(colour = guide_legend(override.aes = list(size=2)))
            #   
            # # png("figures/climate_space_sampled.png", width = 1300, height = 1000, res = 200)
            # p
            # dev.off()
            

            #### Evaluate/explain changes in veg predictions ####
            
            ## for each veg type, determine whether presence predicted by PET coefficient 0.25, both, or 1.00
            
            d_eval = d
            st_geometry(d_eval) = NULL
            
            d_eval_long = d_eval %>%
              mutate(aet_temppptr = aet_tempppt,
                     cwd_temppptr = cwd_tempppt) %>%
              select(-rad.03) %>%
              pivot_longer(cols=c(starts_with("aet_"),starts_with("cwd_"), starts_with("p_"), )) %>%
              # get the second part of the name (it's the climate method)
              mutate(clim_metric = str_split(name,"_") %>% map(2)) %>%
              mutate(valtype_veg = str_c(str_split(name,"_") %>% map(1),str_split(name,"_") %>% map(3),str_split(name,"_") %>% map(4), sep= "_")) %>%
              mutate(valtype_veg = str_replace(valtype_veg,"_NULL","")) %>%
              mutate(valtype_veg = str_replace(valtype_veg,"_NULL","")) %>%
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
            
            ## Put the predictions resulting from different PET coef assumptions side by side
            d_eval_presab = d_eval_long %>%
              filter(!(clim_metric %in% c("tempppt","temppptr"))) %>%
              mutate(cc = str_sub(clim_metric,-3,-1)) %>%
              mutate(clim_metric = str_sub(clim_metric,1,3)) %>%
              select(-aet,-cwd) %>%
              pivot_wider(names_from = cc, values_from = c(prob,presab))
            
            ## Compute a multi-level venn change classification for predicted presence with high PET coef, low, PET coef, both, or neither
            d_eval_presab = d_eval_presab %>%
              mutate(presab_summ = case_when(presab_100 & !presab_025 ~ "high",
                                             presab_025 & !presab_100 ~ "low",
                                             presab_100 & presab_025 ~ "both",
                                             TRUE ~ "neither"))
            
            ### Make a ppt predicted presab file to join to the wb one for computing agreement between the two
            d_eval_presab_ppt = d_eval_long %>%
              filter(clim_metric == "tempppt") %>%
              select(ID,vegtype,presab_ppt = presab) %>%
              mutate(vegtype = as.character(vegtype))
            
            ### And same for ppt radiation
            d_eval_presab_pptr = d_eval_long %>%
              filter(clim_metric == "temppptr") %>%
              select(ID,vegtype,presab_pptr = presab) %>%
              mutate(vegtype = as.character(vegtype))
            
            # Compare the predictions and summarize agreement stats
            d_eval_presab_summ = d_eval_presab %>%
              mutate(vegtype = as.character(vegtype)) %>%
              left_join(d_eval_presab_ppt,by=c("ID","vegtype")) %>%
              left_join(d_eval_presab_pptr,by=c("ID","vegtype")) %>%
              mutate_at(vars(obs,presab_025,presab_100,presab_ppt),as.logical) %>%
              group_by(clim_metric,vegtype) %>%
              summarize(prop_same = sum(presab_summ == "both")/sum(presab_summ != "neither"),
                        tss = tss(presab_025[is.na(training)],presab_100[is.na(training)]),
                        auc_100_val = auc(obs[is.na(training)], prob_100[is.na(training)]) %>% as.numeric(),
                        auc_025_val = auc(obs[is.na(training)], prob_025[is.na(training)]) %>% as.numeric(),
                        tss_ppt_025 = tss(presab_025[is.na(training)],presab_ppt[is.na(training)]),
                        tss_ppt_100 = tss(presab_100[is.na(training)],presab_ppt[is.na(training)]),
                        tss_pptr_025 = tss(presab_025[is.na(training)],presab_pptr[is.na(training)]),
                        tss_pptr_100 = tss(presab_100[is.na(training)],presab_pptr[is.na(training)]),
                        prop_present_100 = sum(presab_summ %in% c("high","both"))/sum(!is.na(presab_summ)),
                        prop_present_obs = sum(obs,na.rm=TRUE)/sum(!is.na(obs))) %>%
              select(clim_metric, vegtype, auc_025_val, auc_100_val, tss, tss_ppt_025, tss_ppt_100, tss_pptr_025, tss_pptr_100, prop_present_100, prop_present_obs) # removed the following for simplicity: , auc_025_train, auc_100_train, auc_025_valid, auc_100_valid, f_stat, kappa, 
            
            
            ## get AUCs for ppt models
            d_eval_ppt = d_eval_long %>%
              filter(clim_metric %in% c("tempppt","temppptr")) %>%
              mutate(vegtype = as.character(vegtype)) %>%
              mutate_at(vars(obs,presab),as.logical) %>%
              group_by(clim_metric,vegtype) %>%
              # what fraction of the cells in that had presence in either were the same in both?
              summarize(auc_ppt = auc(obs[is.na(training)], prob[is.na(training)]) %>% as.numeric) %>%
              mutate(clim_metric = str_sub(clim_metric,-3,-1)) %>%
              pivot_wider(names_from = clim_metric, values_from = auc_ppt) %>%
              rename(auc_ppt = "ppt",
                     auc_ptr = "ptr")
            
            
            ## Table 2 and 3 for main ms: filter to relevant fields, join in PPT AUCs, and display nicely
            
            veg_fits_main = d_eval_presab_summ %>%
              ungroup() %>%
              filter(clim_metric == "dob") %>%
              select(-clim_metric) %>%
              left_join(d_eval_ppt) %>%
              mutate_if(is.numeric,round,digits=2)
            
            veg_fits_main$fold = k
            
            veg_fits_main_overall = bind_rows(veg_fits_main_overall, veg_fits_main)
}

write_csv(veg_fits_main_overall, "tables/veg_fits_main_crossval.csv")

### Summarize the cross-validation veg fits (compute mean and SD across all iterations)

crossval = read_csv("tables/veg_fits_main_crossval.csv")

## get mean and sd for each metric
k_summ = crossval %>%
  select(-fold) %>%
  group_by(vegtype) %>%
  summarize(across(everything(),list(mean = mean, sd = sd), .names="{.col}.{.fn}")) %>%
  ungroup() %>%
  pivot_longer(cols= -vegtype, names_to = "var", values_to="val") %>%
  mutate(metric = str_split(var,fixed(".")) %>% map(1)) %>%
  mutate(metric_type = str_split(var,fixed(".")) %>% map(2)) %>%
  select(-var) %>%
  pivot_wider(names_from = metric_type, values_from = "val")

## make pub-friendly

nicenum = function(x) {
  x %>% round(2) %>% as.string() %>% str_pad(width = 4, side="right",pad=0)
  
}

k_pub = k_summ %>%
  mutate(across(starts_with(c("mean","sd")), ~round(.x,2))) %>%
  mutate(mean_sd = paste0(mean," (",sd,")")) %>%
  select(-mean,-sd) %>%
  pivot_wider(names_from = metric, values_from = mean_sd)

table2 = k_pub %>%
  select(vegtype, auc_100_val, auc_025_val, tss)

write_csv(table2,"tables/Table2_veg_fits_main_crossval_summ.csv")

table3 = k_pub %>%
  select(vegtype, auc_ppt, auc_ptr, tss_ppt_100, tss_ppt_025, tss_pptr_100, tss_pptr_025)

write_csv(table3,"tables/Table3_veg_fits_ppt_crossval_summ.csv")





##### Plot maps of predictions ####

if(cross_validation_iterations == 1) {

  
  #### Presence-absence predictions (and actual vegetation data) maps (Fig. 6)
  
  
  d_plot_pre = d_eval_presab %>%
    mutate(presab_summ = recode(presab_summ, high = "1.00 only",
                                low = "0.25 only",
                                both = "Both",
                                neither = "Neither")) %>%
    mutate(presab_summ = factor(presab_summ, levels=c("0.25 only","Both","1.00 only","Neither")))
  
  
  map_preds = function(d_plot, title, legend) {
    
    p <- ggplot(d_plot) +
      geom_raster(aes(x=x,y=y,fill=presab_summ)) +
      coord_equal() +
      theme_map() +
      theme(legend.position="right") +
      #theme(legend.title=element_blank()) +
      theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
      theme(strip.text=element_text(size=12),strip.background=element_blank()) +
      theme(panel.spacing=unit(1,"lines")) +
      theme(panel.border=element_rect(fill=NA)) +
      theme(legend.text=element_text(size=8)) +
      theme(plot.title = element_text(size=9)) +
      labs(title=title) +
      #theme(plot.title=element_text(hjust=0.5,size=12)) +
      #scale_x_continuous(limit=c(-55000,60000)) +
      scale_fill_manual(values = c("#264653", "#e9c46a", "#e76f51", "grey90"), name = "Presence predicted")
    
    if(!legend) {
      p = p + theme(legend.position = "none")
    }
    
    return(p) 
    
  }
  
  d_plot = d_plot_pre %>%
    filter(clim_metric == "dob",
           vegtype == "mch")
  p2 = map_preds(d_plot, title = "b) Predictions: montane chaparral", legend=FALSE)
  
  d_plot = d_plot_pre %>%
    filter(clim_metric == "dob",
           vegtype == "mhw")
  p3 = map_preds(d_plot, title = "c) Predictions: montane hardwoods", legend=TRUE)
  
  d_plot = d_plot_pre %>%
    filter(clim_metric == "dob",
           vegtype == "smc")
  p4 = map_preds(d_plot, title = "d) Predictions: Sierra mixed-conifer", legend=FALSE)
  
  
  
  #### Map of the 3 focal veg types
  
  d_plot = d %>%
    filter(cv_text %in% c("MCH","MHW","SMC")) %>%
    mutate(cv_text = recode(cv_text,SMC = "Sierra\nmixed-conifer",MHW = "Montane\nhardwoods",MCH = "Montane\nchaparral")) %>%
    mutate(cv_text = factor(cv_text,levels=c("Montane\nchaparral","Montane\nhardwoods","Sierra\nmixed-conifer")))
  
  
  p1 = ggplot(d_plot) +
    geom_sf(data = tuol_mask, fill="grey90", color=NA) +
    geom_tile(aes(fill=cv_text, x=x,y=y)) +
    #coord_equal() +
    theme_map() +
    geom_sf(data=bboxes,fill=NA,color="grey20",size=0.25) +
    #geom_sf(data = tuol_mask, fill=NA, color="grey20") +
    #scale_fill_viridis(na.value="white",discrete=TRUE,option="inferno",name="Vegetation type") +
    theme(legend.position="right") +
    theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
    theme(panel.border=element_rect(fill=NA)) +
    theme(strip.text=element_text(size=12),strip.background=element_blank()) +
    theme(panel.spacing=unit(1,"lines")) +
    #theme(panel.border=element_rect(fill=NA, color=NA)) +
    theme(legend.text=element_text(size=9),legend.title=element_text(size=9)) +
    scale_fill_manual(values=c("#22BBD3","#5C6070", "#C1495D","grey50"), name = "Vegetation type") +
    labs(title="a) Actual vegetation distribution") +
    theme(legend.text = element_text(margin = margin(t = 2,b=2, unit = "pt"),size=8)) +
    theme(plot.title = element_text(size=9))
  
  
  library(gridExtra)
  
  layout = rbind(c(1,1,1),
                 c(NA,2,NA),
                 c(NA,3,3),
                 c(NA,4,NA))
  
  png("figures/veg_pred/veg_pred_main.png",res=200, width = 1000, height = 1000)
  grid.arrange(p1,p2,p3,p4,
               ncol = 3,
               heights = c(1,1,1,1),
               widths = c(0.03,1,0.375),
               layout_matrix = layout)
  dev.off()
  
  
  
  
  
  
  #### Probability-based maps (Fig. S11)
  
  
  ### General data prep
  
  d_plot_pre = d_eval_presab %>%
    mutate(presab_summ = recode(presab_summ, high = "1.00 only",
                                low = "0.25 only",
                                both = "Both",
                                neither = "Neither")) %>%
    mutate(presab_summ = factor(presab_summ, levels=c("0.25 only","Both","1.00 only","Neither"))) %>%
    filter(clim_metric == "dob")
  
  ## pull in precip
  d_ppt = d_eval_long %>%
    filter(clim_metric %in% c("tempppt","temppptr")) %>%
    rename(cc_assumption = "clim_metric",
           prob_presence = "prob") %>%
    mutate(cc_assumption = as.character(cc_assumption))
  
  d_plot_pre = d_plot_pre %>%
    pivot_longer(cols=c(prob_100,prob_025), names_to = "cc_assumption", values_to = "prob_presence" ) %>%
    bind_rows(d_ppt) %>%
    mutate(cc_assumption = recode(cc_assumption, "prob_025" = "Water balance: PET coefficient = 0.25",
                                  "prob_100" = "Water balance: PET coefficient = 1.00",
                                  "tempppt" = "Temperature and precipitation",
                                  "temppptr" = "Temperature, precip., solar rad., and interactions")) %>%
    mutate(cc_assumption = factor(cc_assumption,levels = c("Water balance: PET coefficient = 0.25", "Water balance: PET coefficient = 1.00", "Temperature and precipitation", "Temperature, precip., solar rad., and interactions")))
  
  
  ## Plot for montane chaparral
  d_plot = d_plot_pre %>%
    filter(vegtype == "mch")
  
  p_mch <- ggplot(d_plot) +
    geom_raster(aes(x=x,y=y,fill=prob_presence)) +
    coord_equal() +
    theme_map() +
    theme(legend.position="right") +
    #theme(legend.title=element_blank()) +
    theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
    theme(strip.text=element_text(size=11),strip.background=element_blank()) +
    theme(panel.spacing=unit(1,"lines")) +
    theme(panel.border=element_rect(fill=NA)) +
    theme(legend.text=element_text(size=8)) +
    theme(plot.title=element_text(size=16,hjust = 0.5)) +
    labs(title="Montane chaparral") +
    #theme(plot.title=element_text(hjust=0.5,size=12)) +
    #scale_x_continuous(limit=c(-55000,60000)) +
    #scale_fill_manual(values = c("#264653", "#e9c46a", "#e76f51", "grey90"), name = "Presence predicted") +
    scale_fill_viridis(limits=c(0.15,0.67), name="Prob.\n presence") +
    facet_wrap(~cc_assumption)
  
  
  ## Plot for montane hardwoods
  d_plot = d_plot_pre %>%
    filter(vegtype == "mhw")
  
  p_mhw <- ggplot(d_plot) +
    geom_raster(aes(x=x,y=y,fill=prob_presence)) +
    coord_equal() +
    theme_map() +
    theme(legend.position="right") +
    #theme(legend.title=element_blank()) +
    theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
    theme(strip.text=element_text(size=11),strip.background=element_blank()) +
    theme(panel.spacing=unit(1,"lines")) +
    theme(panel.border=element_rect(fill=NA)) +
    theme(legend.text=element_text(size=8)) +
    #theme(plot.title = element_text(size=9)) +
    labs(title="Montane hardwoods") +
    theme(plot.title=element_text(size=16,hjust = 0.5)) +
    #theme(plot.title=element_text(hjust=0.5,size=12)) +
    #scale_x_continuous(limit=c(-55000,60000)) +
    #scale_fill_manual(values = c("#264653", "#e9c46a", "#e76f51", "grey90"), name = "Presence predicted") +
    scale_fill_viridis(limits=c(0.15,0.99), name="Prob.\n presence") +
    facet_wrap(~cc_assumption)
  
  
  ## Plot for sierra mixed-conifer
  d_plot = d_plot_pre %>%
    filter(vegtype == "smc")
  
  p_smc <- ggplot(d_plot) +
    geom_raster(aes(x=x,y=y,fill=prob_presence)) +
    coord_equal() +
    theme_map() +
    theme(legend.position="right") +
    #theme(legend.title=element_blank()) +
    theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white")) +
    theme(strip.text=element_text(size=11),strip.background=element_blank()) +
    theme(panel.spacing=unit(1,"lines")) +
    theme(panel.border=element_rect(fill=NA)) +
    theme(legend.text=element_text(size=8)) +
    theme(plot.title=element_text(size=16,hjust = 0.5)) +
    labs(title="Sierra mixed-conifer") +
    #theme(plot.title=element_text(hjust=0.5,size=12)) +
    #scale_x_continuous(limit=c(-55000,60000)) +
    #scale_fill_manual(values = c("#264653", "#e9c46a", "#e76f51", "grey90"), name = "Presence predicted") +
    scale_fill_viridis(limits=c(0.15,0.93), name="Prob.\n presence") +
    facet_wrap(~cc_assumption)
  
  png("figures/veg_pred/veg_pred_prob.png",res=150, width = 1200, height = 1500)
  grid.arrange(p_mch, p_mhw, p_smc,
               ncol = 1)
  dev.off()

}


