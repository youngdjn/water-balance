# PET_methods:
#     Various ways of calculating PET from WB.params (see Point_WB.R for description).
#     Not all methods require all elements of WB.params to be populated.
#
#     Each method returns a 12-element vector of monthly PET values.
#

PET.PT <- function(WB.params) {
  
  T.m <- WB.params$T.m
  Rs.m <- WB.params$R.m * 0.0036 / monthlengths()
  Tmin.m <- WB.params$Tmin.m
  Tmax.m <- WB.params$Tmax.m
  cloudiness.m <- WB.params$cloudiness.m
  Tdew.m <- WB.params$Tdew.m
  month <- WB.params$month
  lat <- WB.params$L
  elev <- WB.params$E
  rad.nldas.m <- WB.params$R.nldas.m * 0.0036 / monthlengths()
  
  Tmin.m.k <- Tmin.m + 273.16
  Tmax.m.k <- Tmax.m + 273.16
  
  ea.m <- 0.6108 * exp((17.27*Tdew.m)/(Tdew.m + 237.3))
  
  
  ## Begin code for computing Rn.m (net rad, which is incoming minus outgoing, from Dobrowski)
  ### Compute cloudiness using, for clear-sky, empirical latitude- and date-based (following Dobrowski) and for true-sky, NLDAS
  
  # Calculate potential max solar radiation or clear sky radiation	
  daysinmonth=c(31,28,31,30,31,30,31,31,30,31,30,31)
  d2=c(31,59,90,120,151,181,212,243,273,304,334,365)
  d1=c(1,32,60,91,121,152,182,213,244,274,305,335)
  DoY <- (d1[month]+d2[month])/2 # use middle day of month to represent monthly average. 
  n_days <- daysinmonth[month]
  
  
  GSC=0.082      # MJ m -2 min-1 (solar constant)
  phi <- pi*lat/180 
  dr <- 1+0.033*cos(2*pi/365*DoY)      
  delt <- 0.409*sin(2*pi/365*DoY-1.39)     
  omegas <- acos(-tan(phi)*tan(delt)) 
  Ra <- 24*60/pi*GSC*dr*(omegas*sin(phi)*sin(delt) +cos(phi)*cos(delt)*sin(omegas))    # Daily extraterrestrial radiation
  Rso <- Ra*(0.75+2e-5*elev)     #For a cloudless day, Rs is roughly 75% of extraterrestrial radiation (Ra)
  
  #Rso from Dobrowski is per day. Convert back to P-T units
  Rso <- Rso
  
  
  # radfraction is a measure of relative shortwave radiation, or of
  # possible radiation (cloudy vs. clear-sky)
  radfraction <- rad.nldas.m/Rso
  radfraction[radfraction>1] <- 1
  
  cloudiness.m <- radfraction
  
  
  
  
  fcd.m <- 1.35*cloudiness.m-0.35
  fcd.m[fcd.m < 0.3] <- 0.3
  fcd.m[fcd.m > 1] <- 1
  
  
  Rns.m <- Rs.m - 0.23 * Rs.m
  
  Rnl.m <- 4.901e-9 * fcd.m * (0.34-0.14*sqrt(ea.m)) * ((Tmax.m.k^4 + Tmin.m.k^4)/2)
  
  Rn.m <- Rns.m - Rnl.m 
  
  ## End code to get Rn.m
  
  lambda <- 2.45
  s <- SatVapPresSlope(T.m)
  press <- ElevToPress(WB.params$E)
  gamma  <- 0.00163*press/lambda  
  alpha <- 1.26
  PET.m <- (1/lambda) * ((s * Rn.m) / (s+gamma)) * alpha * monthlengths()
  PET.m[PET.m < 0] <- 0
  
  return(PET.m)
}



