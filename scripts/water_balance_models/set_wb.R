source("scripts/water_balance_models/point_wb.R")
source("scripts/water_balance_models/wb_utils.R")
source("scripts/water_balance_models/dobrowski_wb.R")

set_wb <- function(in.data, PET.methods, PET.mods, AET.methods, monthly=FALSE, dobr.wb = FALSE) {

  nplots <- nrow(in.data)
  if(length(unique(in.data$ID)) != nplots)
    stop("Identifiers in column 'ID' must be unique")
  
  #in.data <- in.data[!is.na(in.data$rad.01),] ## remove a row with no rad data

  T.m <- in.data[, grep("tmean.01", names(in.data))[1]:grep("tmean.12", names(in.data))]
  Tmin.m <- in.data[, grep("tmin.01", names(in.data))[1]:grep("tmin.12", names(in.data))]
  Tmax.m <- in.data[, grep("tmax.01", names(in.data))[1]:grep("tmax.12", names(in.data))]
  Tdew.m <- in.data[, grep("td.01", names(in.data))[1]:grep("td.12", names(in.data))]
  P.m <- in.data[, grep("ppt.01", names(in.data))[1]:grep("ppt.12", names(in.data))]
  R.m <- in.data[, grep("rad.01", names(in.data))[1]:grep("rad.12", names(in.data))]
  R.nldas.m <- in.data[, grep("rad.nldas.01",names(in.data))[1]:grep("rad.nldas.12",names(in.data))]
  wind.m <- in.data[, grep("wind.01", names(in.data))[1]:grep("wind.12", names(in.data))]    
  PET.BCM.m <- in.data[, grep("bcm.pet.01", names(in.data))[1]:grep("bcm.pet.12", names(in.data))]
  T.m <- apply(T.m, 2, as.numeric)
  Tmin.m <- apply(Tmin.m, 2, as.numeric)
  Tmax.m <- apply(Tmax.m, 2, as.numeric)
  Tdew.m <- apply(Tdew.m, 2, as.numeric)
  P.m <- apply(P.m, 2, as.numeric)
  R.m <- apply(R.m, 2, as.numeric)
  R.nldas.m <- apply(R.nldas.m, 2, as.numeric)
  wind.m <- apply(wind.m, 2, as.numeric)
  PET.BCM.m <- apply(PET.BCM.m, 2, as.numeric)
  colnames(T.m) <- colnames(P.m) <- colnames(R.m) <- colnames(R.nldas.m) <- colnames(Tmin.m) <- colnames(Tmax.m) <- colnames(Tdew.m) <- colnames(wind.m) <- colnames(PET.BCM.m) <- NULL
  
  month.param <- 1:12
  
  L     <- in.data$Lat
  E     <- in.data$elev
  S.max <- 150
  
  wind.avg.m <- apply(wind.m,2,mean,na.rm=TRUE)
  wind.avg.m <- matrix(wind.avg.m,nrow=nrow(T.m),ncol=ncol(T.m),byrow=TRUE)
  
  #this is the predicted soil whc from the regression
  S.max.reg <- -3.618513 - 0.017472*E + 3.309366*L
  
  AET = data.frame()
  PET = data.frame()
  Deficit = data.frame()
  
  browser()
  
  for(i in 1:nplots) {
    
    #run for cc100
    params <- list(T.m=T.m[i,], P.m=P.m[i,], R.m=R.m[i,], R.nldas.m=R.nldas.m[i,], L=L[i], E=E[i], S.max=150, PET.BCM.m=PET.BCM.m[i,], Tmin.m=Tmin.m[i,], Tmax.m=Tmax.m[i,], wind.m=wind.m[i,],month=month.param,Tdew.m=Tdew.m[i,])
    
    output <- Point_WB(params, PET.methods, PET.mods, AET.methods, monthly=monthly)
    
    PET = bind_rows(PET, output$PET)
    AET = bind_rows(AET, output$AET)
    Deficit = bind_rows(Deficit, output$Deficit)
    

    
    cat(paste("\rFinished plot number:",i," of ",nplots))
    
  }
  
  if(monthly) {
    methods.mult <- 12
  } else {
    methods.mult <- 1
  }

  
  rownames(PET) <- rownames(AET) <- rownames(Deficit) <- in.data$ID
  
  names(PET) <- paste0("PET.",names(PET))
  names(AET) <- paste0("AET.",names(AET))
  names(Deficit) <- paste0("Deficit.",names(Deficit))

  ##### Also compute Dobrowski WB if specified ###
  if(dobr.wb) {
    
    R.m <- R.m *.0036
    R.m <- t(apply(R.m, 1, function(x) {x / monthlengths()}))
    
    R.nldas.m <- R.nldas.m *.0036
    R.nldas.m <- t(apply(R.nldas.m, 1, function(x) {x / monthlengths()}))
    
    ### Fixed 150mm soil; fixed monthly wind

    
    # dob.inputs <- data.frame(T.m,P.m,R.m,R.nldas.m,Tmin.m,Tmax.m,wind.avg.m,Tdew.m,L,E)
    # compl <- complete.cases(dob.inputs)
    # dob.inputs$incomplete <- 1
    # dob.inputs[compl,"incomplete"] <- 0
    
    dob_wb_ret <- run_dob_wb(T.m,P.m,R.m,R.nldas.m,Tmin.m,Tmax.m,wind.avg.m,Tdew.m,L,E,rep(150,nrow(T.m)),pet.mult=1)
    
    PET.Dobr.cc100 <- dob_wb_ret$PET.Dobr
    AET.Dobr.cc100 <- dob_wb_ret$AET.Dobr
    Deficit.Dobr.cc100 <- dob_wb_ret$Deficit.Dobr
    
    dob_wb_ret <- run_dob_wb(T.m,P.m,R.m,R.nldas.m,Tmin.m,Tmax.m,wind.avg.m,Tdew.m,L,E,rep(150,nrow(T.m)),pet.mult=0.25)
    
    PET.Dobr.cc025 <- dob_wb_ret$PET.Dobr
    AET.Dobr.cc025 <- dob_wb_ret$AET.Dobr
    Deficit.Dobr.cc025 <- dob_wb_ret$Deficit.Dobr
    
    PET.dob <- cbind(PET.Dobr.cc100,PET.Dobr.cc025)
    AET.dob <- cbind(AET.Dobr.cc100,AET.Dobr.cc025)
    Deficit.dob <- cbind(Deficit.Dobr.cc100,Deficit.Dobr.cc025)
    
    PET.dob <- as.data.frame(PET.dob)
    AET.dob <- as.data.frame(AET.dob)
    Deficit.dob <- as.data.frame(Deficit.dob)
    
    rownames(PET.dob) <- rownames(AET.dob) <- rownames(Deficit.dob) <- in.data$ID
    
  }

  
  # merge dobrowski and non-dobrowski
  
  PET.all <- cbind(PET,PET.dob)
  AET.all <- cbind(AET,AET.dob)
  Deficit.all <- cbind(Deficit,Deficit.dob)
  
  wb_output = cbind(PET.all, AET.all, Deficit.all)
  wb_output$ID = in.data$ID
  
  return(wb_output)
}