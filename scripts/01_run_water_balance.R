## Compute Dobrowski and Thornthwaite water balance values across the study landscape
## Input values are provided as a raster with parameters named as expected by (and specified in) "water_balance_models/point_wb.R" 

library(raster)
library(tidyverse)
library(furrr)
plan(multiprocess,workers=8)

source("scripts/water_balance_models/set_wb.R")

## Load raster grid of input values
wb.inputs <- brick("data/wb_input/wb_input_vars.grd")

## Make a layer that contains cell IDs, for loading values back in later
id.layer <- wb.inputs[[1]]
values(id.layer) <- seq_along(id.layer)
wb.inputs[["ID"]] <- id.layer

## Convert to data frame
wb.inputs.df <- as.data.frame(wb.inputs)
wb.inputs.df <- as(wb.inputs,"SpatialPointsDataFrame")
# project to lat/lon
wb.inputs.df <- spTransform(wb.inputs.df,CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs"))
wb.inputs.df$Lat <- coordinates(wb.inputs.df)[,2]
wb.inputs.df$Lon <- coordinates(wb.inputs.df)[,1]
wb.inputs.df <- as.data.frame(wb.inputs.df)

PET.methods <- "PT"
PET.mods <- c("STD","cc025")
AET.methods <- "Wil150mm"

## Define base parameter DF
wbparams_base <- wb.inputs.df[complete.cases(wb.inputs.df),] %>%
  mutate(scenario = "base")



#### Summarize site climate across all sites

wbparams_summ = wbparams_base %>%
  select(starts_with("ppt."),starts_with("tmin."),starts_with("tmax.")) %>%
  summarize_all(mean)

wbparams_summ %>% select(starts_with("ppt.")) %>% as.vector %>% dput
wbparams_summ %>% select(starts_with("tmin.")) %>% as.vector %>% dput
wbparams_summ %>% select(starts_with("tmax.")) %>% as.vector %>% dput




## ALTERNATE PARAMETER DFS

## Get the mean temp and ppt across our dataset
tmin01_mean = wbparams_base %>% select(tmin.01) %>% pull %>% mean
tmin07_mean = wbparams_base %>% select(tmin.07) %>% pull %>% mean

tmax01_mean = wbparams_base %>% select(tmax.01) %>% pull %>% mean
tmax07_mean = wbparams_base %>% select(tmax.07) %>% pull %>% mean

ppt01_mean = wbparams_base %>% select(ppt.01) %>% pull %>% mean
ppt07_mean = wbparams_base %>% select(ppt.07) %>% pull %>% mean


#### Get temp absolute deviation from overall mean for each month, and ppt relative deviation

wbparams_rel = wbparams_base %>%
  mutate(tot_ppt = ppt.01+ppt.02+ppt.03+ppt.04+ppt.05+ppt.06+ppt.07+ppt.08+ppt.09+ppt.10+ppt.11+ppt.12,
         mean_temp = (tmin.01+tmin.02+tmin.03+tmin.04+tmin.05+tmin.06+tmin.07+tmin.08+tmin.09+tmin.10+tmin.11+tmin.12+tmax.01+tmax.02+tmax.03+tmax.04+tmax.05+tmax.06+tmax.07+tmax.08+tmax.09+tmax.10+tmax.11+tmax.12)/24)

overall_mean_annual_ppt = mean(wbparams_rel$tot_ppt)
overall_mean_annual_temp = mean(wbparams_rel$mean_temp)

## temp and ppt offset is the same across months for each plot
wbparams_rel = wbparams_rel %>%
   mutate(across(starts_with("ppt."),~tot_ppt/overall_mean_annual_ppt,.names="dev_{.col}"),
          across(starts_with("tmin."),~mean_temp-overall_mean_annual_temp,.names="dev_{.col}"),
          across(starts_with("tmax."),~mean_temp-overall_mean_annual_temp,.names="dev_{.col}"))
 
adjust_clim = function(df, base_name,dev_name,monthly_mean, multiply=FALSE) {
  monthly_mean_rep = replicate(nrow(df),monthly_mean) %>% t
  base_col_range = which(names(df) %>% str_starts(fixed(base_name)))
  dev_col_range = which(names(df) %>% str_starts(fixed(dev_name)))
  if(multiply) {
    df[,base_col_range] = df[,dev_col_range] * monthly_mean_rep
  } else {
    df[,base_col_range] = df[,dev_col_range] + monthly_mean_rep
  }
  return(df)
}

compute_tmean = function(df) {
  tmean_cols = which(names(df) %>% str_starts(fixed("tmean.")))
  tmin_cols = which(names(df) %>% str_starts(fixed("tmin.")))
  tmax_cols = which(names(df) %>% str_starts(fixed("tmax.")))
  df[,tmean_cols] = (df[,tmin_cols] + df[,tmax_cols])/2
  return(df)  
}


#### To shift precipitation:
## Get each site's tot annual precip
## Get each site's deviation from tot annual across all sites
## Apply this same deviation to all months




#### Zone 4 (hot desert)

tmin_monthly_mean = c(3.38, 4.76931656275537, 8.565, 13.75, 18.935, 22.7306834372446, 
                      24.12, 22.7306834372446, 18.935, 13.75, 8.56500000000001, 4.76931656275537
)
tmax_monthly_mean = c(15.58, 17.1488425216842, 21.435, 27.29, 33.145, 37.4311574783158, 
                      39, 37.4311574783158, 33.145, 27.29, 21.435, 17.1488425216842
)
ppt_monthly_mean = c(18.29, 17.4955306444417, 15.325, 12.36, 9.395, 7.22446935555828, 
                     6.43, 7.22446935555828, 9.395, 12.36, 15.325, 17.4955306444417
)

wbparams_foc = wbparams_rel
wbparams_foc = adjust_clim(wbparams_foc, base_name = "tmax.", dev_name = "dev_tmax.", monthly_mean = tmax_monthly_mean)
wbparams_foc = adjust_clim(wbparams_foc, base_name = "tmin.", dev_name = "dev_tmin.", monthly_mean = tmin_monthly_mean)
wbparams_foc = adjust_clim(wbparams_foc, base_name = "ppt.", dev_name = "dev_ppt.", monthly_mean = ppt_monthly_mean, multiply=TRUE)
wbparams_z04 = wbparams_foc %>% compute_tmean() %>% select(-starts_with("dev_"),-tot_ppt,-mean_temp) %>%
  mutate(scenario = "Zone_04")


#### Zone 7

tmin_monthly_mean = c(-8.52, -7.0322121090262, -2.9675, 2.585, 8.1375, 12.2022121090262, 
                      13.69, 12.2022121090262, 8.1375, 2.585, -2.96749999999999, -7.03221210902619
)
tmax_monthly_mean = c(3.35, 5.05616648280517, 9.7175, 16.085, 22.4525, 27.1138335171948, 
                      28.82, 27.1138335171948, 22.4525, 16.085, 9.71750000000001, 5.05616648280518
)
ppt_monthly_mean = c(24.67, 25.3834147248479, 27.3325, 29.995, 32.6575, 34.6065852751521, 
                     35.32, 34.6065852751521, 32.6575, 29.995, 27.3325, 25.3834147248479
)

wbparams_foc = wbparams_rel
wbparams_foc = adjust_clim(wbparams_foc, base_name = "tmax.", dev_name = "dev_tmax.", monthly_mean = tmax_monthly_mean)
wbparams_foc = adjust_clim(wbparams_foc, base_name = "tmin.", dev_name = "dev_tmin.", monthly_mean = tmin_monthly_mean)
wbparams_foc = adjust_clim(wbparams_foc, base_name = "ppt.", dev_name = "dev_ppt.", monthly_mean = ppt_monthly_mean, multiply=TRUE)
wbparams_z07 = wbparams_foc %>% compute_tmean() %>% select(-starts_with("dev_"),-tot_ppt,-mean_temp) %>%
  mutate(scenario = "Zone_07")


#### Zone 11

tmin_monthly_mean = c(-0.329374099999999, 1.13836744679, 5.148311925, 10.62599795, 
                      16.103683975, 20.11362845321, 21.58137, 20.11362845321, 16.103683975, 
                      10.62599795, 5.14831192500001, 1.13836744679)
tmax_monthly_mean = c(8.4688886, 9.89140349160682, 13.77778645, 19.0866843, 24.39558215, 
                      28.2819651083932, 29.70448, 28.2819651083932, 24.39558215, 19.0866843, 
                      13.77778645, 9.89140349160683)
ppt_monthly_mean = c(25.15121, 34.207457688658, 58.9495865, 92.747963, 126.5463395, 
                     151.288468311342, 160.344716, 151.288468311342, 126.5463395, 
                     92.747963, 58.9495865, 34.2074576886581)

wbparams_foc = wbparams_rel
wbparams_foc = adjust_clim(wbparams_foc, base_name = "tmax.", dev_name = "dev_tmax.", monthly_mean = tmax_monthly_mean)
wbparams_foc = adjust_clim(wbparams_foc, base_name = "tmin.", dev_name = "dev_tmin.", monthly_mean = tmin_monthly_mean)
wbparams_foc = adjust_clim(wbparams_foc, base_name = "ppt.", dev_name = "dev_ppt.", monthly_mean = ppt_monthly_mean, multiply=TRUE)
wbparams_z11 = wbparams_foc %>% compute_tmean() %>% select(-starts_with("dev_"),-tot_ppt,-mean_temp) %>%
  mutate(scenario = "Zone_11")


#### Zone 14

tmin_monthly_mean = c(-0.463983899999999, 0.920737528368342, 4.703866825, 9.87171755, 
                      15.039568275, 18.8226975716317, 20.207419, 18.8226975716317, 
                      15.039568275, 9.87171755, 4.70386682500001, 0.920737528368345
)
tmax_monthly_mean = c(10.0245175, 11.4444681915213, 15.323845625, 20.62317375, 25.922501875, 
                      29.8018793084787, 31.22183, 29.8018793084787, 25.922501875, 20.62317375, 
                      15.323845625, 11.4444681915213)
ppt_monthly_mean = c(81.72432, 84.6661159631562, 92.703252, 103.682184, 114.661116, 
                     122.698252036844, 125.640048, 122.698252036844, 114.661116, 103.682184, 
                     92.703252, 84.6661159631562)

wbparams_foc = wbparams_rel
wbparams_foc = adjust_clim(wbparams_foc, base_name = "tmax.", dev_name = "dev_tmax.", monthly_mean = tmax_monthly_mean)
wbparams_foc = adjust_clim(wbparams_foc, base_name = "tmin.", dev_name = "dev_tmin.", monthly_mean = tmin_monthly_mean)
wbparams_foc = adjust_clim(wbparams_foc, base_name = "ppt.", dev_name = "dev_ppt.", monthly_mean = ppt_monthly_mean, multiply=TRUE)
wbparams_z14 = wbparams_foc %>% compute_tmean() %>% select(-starts_with("dev_"),-tot_ppt,-mean_temp) %>%
  mutate(scenario = "Zone_14")



#### Zone 29

tmin_monthly_mean = c(-22.0210038, -20.406454451481, -15.9954236, -9.9698434, -3.9442632, 
                      0.466767651481035, 2.081317, 0.466767651481039, -3.94426319999999, 
                      -9.96984339999999, -15.9954236, -20.406454451481)
tmax_monthly_mean = c(-7.638623, -6.20932285778816, -2.30440225, 3.0298185, 8.36403925, 
                      12.2689598577882, 13.69826, 12.2689598577882, 8.36403925, 3.0298185, 
                      -2.30440224999999, -6.20932285778816)
ppt_monthly_mean = c(11.68857, 16.5887950218327, 29.97645875, 48.2643475, 66.55223625, 
                     79.9398999781673, 84.840125, 79.9398999781673, 66.55223625, 48.2643475, 
                     29.97645875, 16.5887950218327)

wbparams_foc = wbparams_rel
wbparams_foc = adjust_clim(wbparams_foc, base_name = "tmax.", dev_name = "dev_tmax.", monthly_mean = tmax_monthly_mean)
wbparams_foc = adjust_clim(wbparams_foc, base_name = "tmin.", dev_name = "dev_tmin.", monthly_mean = tmin_monthly_mean)
wbparams_foc = adjust_clim(wbparams_foc, base_name = "ppt.", dev_name = "dev_ppt.", monthly_mean = ppt_monthly_mean, multiply=TRUE)
wbparams_z29 = wbparams_foc %>% compute_tmean() %>% select(-starts_with("dev_"),-tot_ppt,-mean_temp) %>%
  mutate(scenario = "Zone_29")


# 
# 
# 
# 
# 
# ## Aseasonal temp
# wbparams_consttemp = wbparams_base %>%
#   mutate_at(vars(starts_with("tmin.")), ~tmin_mean) %>%
#   mutate_at(vars(starts_with("tmean.")), ~tmean_mean) %>%
#   mutate_at(vars(starts_with("tmax.")), ~tmax_mean) %>%
#   mutate_at(vars(starts_with("td.")), ~td_mean) %>%
#   mutate(scenario = "consttemp")
# 
# ## Aseasonal precip
# ppt_mean = wbparams_base %>% select(starts_with("ppt.")) %>% rowMeans
# wbparams_constppt = wbparams_base %>%
#   mutate_at(vars(starts_with("ppt.")), ~ppt_mean) %>%
#   mutate(scenario = "constppt")
# 
# ## Aseasonal temp and precip
# wbparams_consttempppt = wbparams_base %>%
#   mutate_at(vars(starts_with("ppt.")), ~ppt_mean) %>%
#   mutate_at(vars(starts_with("tmin.")), ~tmin_mean) %>%
#   mutate_at(vars(starts_with("tmean.")), ~tmean_mean) %>%
#   mutate_at(vars(starts_with("tmax.")), ~tmax_mean) %>%
#   mutate_at(vars(starts_with("td.")), ~td_mean) %>%
#   mutate(scenario = "consttempppt")
# 
# ## Aseasonal temp and precip, hot dry
# wbparams_consttemppptHotdry = wbparams_base %>%
#   mutate_at(vars(starts_with("ppt.")), ~ppt_mean/2) %>%
#   mutate_at(vars(starts_with("tmin.")), ~tmin_mean+10) %>%
#   mutate_at(vars(starts_with("tmean.")), ~tmean_mean+10) %>%
#   mutate_at(vars(starts_with("tmax.")), ~tmax_mean+10) %>%
#   mutate_at(vars(starts_with("td.")), ~td_mean) %>%
#   mutate(scenario = "consttemppptHotdry")
# 
# ## Double precip
# wbparams_doubleppt = wbparams_base %>%
#   mutate_at(vars(starts_with("ppt.")), ~.*2) %>%
#   mutate(scenario = "doubleppt")
# 
# ## Double aseasonal precip, aseasonal temp
# wbparams_consttemppptdoubleppt = wbparams_consttempppt %>%
#   mutate_at(vars(starts_with("ppt.")), ~.*2) %>%
#   mutate(scenario = "consttemppptdoubleppt")


## LIST OF PARAMETERS DFs
inputs = list(base = wbparams_base,
              z04 = wbparams_z04,
              z07 = wbparams_z07,
              z11 = wbparams_z11,
              z14 = wbparams_z14,
              z29 = wbparams_z29)

# inputs = inputs[6:7]

wb = future_map(inputs,set_wb,PET.methods = PET.methods, PET.mods = PET.mods, AET.methods = AET.methods, monthly=FALSE, dobr.wb = TRUE,
               .progress = TRUE,
               .options = future_options(globals = c("PET.PT","Point_WB","monthlengths","SatVapPresSlope","ElevToPress","PET.mod.STD","PET.mod.cc025","AET.Wil150mm","AET.Wil","run_dob_wb","dobr_wb","snowmod","monthlyET0","aetmod"),
                                         packages = "dplyr"))



# Compile and save Data frames of inputs and outputs , plus rasters: one for each scenario

for(i in 1:length(inputs)) {
  scenario = names(inputs)[i]
  inputs_foc = inputs[[scenario]]
  outputs_foc = wb[[scenario]]
  in_out_foc = left_join(inputs_foc, outputs_foc,by="ID") %>%
    select(-scenario)
  write.csv(in_out_foc,paste0("data/wb_output/",scenario,".csv"),row.names=FALSE)
  
  
  
  ## Add all the modeled WB values/layers to the WB inputs DF so that outputs are stacked with inputs
  ## arrange the wb.nonwb DF with IDs in the same order as the raster cells
  order <- match(raster::values(id.layer),in_out_foc$ID)
  d.rast.compare <- in_out_foc[order,]
  
  wb.inputs.foc = wb.inputs
  
  for (j in 1:ncol(d.rast.compare)) {
    col <- d.rast.compare[,j]
    name <- names(d.rast.compare)[j]

    if(!(name %in% names(wb.inputs.foc))) { # if the layer we are adding is not a name in the raster, it was just added at the end with a genericc name, so rename it
      newlayer = wb.inputs.foc[[1]] * 0
      values(newlayer) = col
      names(newlayer) = name
      wb.inputs.foc = addLayer(wb.inputs.foc,newlayer)
    } else { # it already exists, so replace it
      wb.inputs.foc[[name]] <- col
    }
    
  }
  
  
  writeRaster(wb.inputs.foc,paste0("data/wb_output/rasters/",scenario,".grd"), overwrite=TRUE)
  
  
}

