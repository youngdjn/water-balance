# This script defines a function to model water balance for many points (e.g., a landscape). It is provided a single climate scenario worth of points at a time. It is called from scripts/01_run_water_balance.R
# It is more complicated than it needs to be for the GEB analysis because it is adapted from code that can apply multiple water balance models to a single dataset
# Columns returned ending in cc100 and cc025 indicate water balance values modeled assuming a PET coefficient (or crop coefficient, "cc") of 1.00 and 0.25 respectively.


# Load the functions for calling the Dobrowski water balance model
source("scripts/water_balance_functions/dobrowski_wb.R")

# Constant defining number of days per month
monthlengths = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

set_wb <- function(in.data, PET.methods, PET.mods, AET.methods) {

  nplots <- nrow(in.data)
  if(length(unique(in.data$ID)) != nplots)
    stop("Identifiers in column 'ID' must be unique")
  
  # Load all input variables
  
  T.m <- in.data[, grep("tmean.01", names(in.data))[1]:grep("tmean.12", names(in.data))]
  Tmin.m <- in.data[, grep("tmin.01", names(in.data))[1]:grep("tmin.12", names(in.data))]
  Tmax.m <- in.data[, grep("tmax.01", names(in.data))[1]:grep("tmax.12", names(in.data))]
  Tdew.m <- in.data[, grep("td.01", names(in.data))[1]:grep("td.12", names(in.data))]
  P.m <- in.data[, grep("ppt.01", names(in.data))[1]:grep("ppt.12", names(in.data))]
  R.m <- in.data[, grep("rad.01", names(in.data))[1]:grep("rad.12", names(in.data))]
  R.nldas.m <- in.data[, grep("rad.nldas.01",names(in.data))[1]:grep("rad.nldas.12",names(in.data))]
  wind.m <- in.data[, grep("wind.01", names(in.data))[1]:grep("wind.12", names(in.data))]
  L     <- in.data$Lat
  E     <- in.data$elev
  
  # Object type conversions
  T.m <- apply(T.m, 2, as.numeric)
  Tmin.m <- apply(Tmin.m, 2, as.numeric)
  Tmax.m <- apply(Tmax.m, 2, as.numeric)
  Tdew.m <- apply(Tdew.m, 2, as.numeric)
  P.m <- apply(P.m, 2, as.numeric)
  R.m <- apply(R.m, 2, as.numeric)
  R.nldas.m <- apply(R.nldas.m, 2, as.numeric)
  wind.m <- apply(wind.m, 2, as.numeric)
  colnames(T.m) <- colnames(P.m) <- colnames(R.m) <- colnames(R.nldas.m) <- colnames(Tmin.m) <- colnames(Tmax.m) <- colnames(Tdew.m) <- colnames(wind.m) <- NULL
  
  month.param <- 1:12
  
  # Manually specify soil water holding capacity = 150mm
  S.max <- rep(150,nrow(T.m))
  
  # Compute average wind speed across the landscape for each month, so that wind is not spatially variable, to facilitate interpretation of the drivers of water balance variation
  wind.avg.m <- apply(wind.m,2,mean,na.rm=TRUE)
  wind.avg.m <- matrix(wind.avg.m,nrow=nrow(T.m),ncol=ncol(T.m),byrow=TRUE)
  


  ##### Compute Dobrowski WB ###
  
  ## Convert rad to units expected by Dobrowski method (MJ/m^2/day)
  
  R.m <- R.m *.0036
  R.m <- t(apply(R.m, 1, function(x) {x / monthlengths}))
  
  R.nldas.m <- R.nldas.m *.0036
  R.nldas.m <- t(apply(R.nldas.m, 1, function(x) {x / monthlengths}))
  
  ## Run Dobrowski water balance for the supplied points, assuming PET coefficient (pet.mult) is 1
  # This function is defined in scripts/water_balance_functions/dobrowski_wb.R
  dob_wb_ret <- run_dob_wb(T.m,P.m,R.m,R.nldas.m,Tmin.m,Tmax.m,wind.avg.m,Tdew.m,L,E,S.max,pet.mult=1)
  
  PET.Dobr.cc100 <- dob_wb_ret$PET.Dobr
  AET.Dobr.cc100 <- dob_wb_ret$AET.Dobr
  Deficit.Dobr.cc100 <- dob_wb_ret$Deficit.Dobr
  
  ## Re-run Dobrowski water balance for the supplied points, assuming PET coefficient (pet.mult) is 0.25
  dob_wb_ret <- run_dob_wb(T.m,P.m,R.m,R.nldas.m,Tmin.m,Tmax.m,wind.avg.m,Tdew.m,L,E,S.max,pet.mult=0.25)
  
  PET.Dobr.cc025 <- dob_wb_ret$PET.Dobr
  AET.Dobr.cc025 <- dob_wb_ret$AET.Dobr
  Deficit.Dobr.cc025 <- dob_wb_ret$Deficit.Dobr
  
  # Compile the separate PET, AET, Deficit outputs into a single DF and return it
  
  PET.dob <- cbind(PET.Dobr.cc100,PET.Dobr.cc025)
  AET.dob <- cbind(AET.Dobr.cc100,AET.Dobr.cc025)
  Deficit.dob <- cbind(Deficit.Dobr.cc100,Deficit.Dobr.cc025)
  
  PET.dob <- as.data.frame(PET.dob)
  AET.dob <- as.data.frame(AET.dob)
  Deficit.dob <- as.data.frame(Deficit.dob)
  
  wb_output = bind_cols(PET.dob,AET.dob,Deficit.dob)
  
  wb_output$ID = in.data$ID
  
  return(wb_output)
}