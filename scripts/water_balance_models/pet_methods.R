#
# PET_methods:
#     Various ways of calculating PET from WB.params (see Point_WB.R for description).
#     Not all methods require all elements of WB.params to be populated. See documentation for details.
#
#     Each method returns a 12-element vector of monthly PET values.
#

#source("utils.R")

PET.Thorn <- function(WB.params)
{
  
  T.m     <- WB.params$T.m
  T.m.frz <- T.m < 0
  T.m.hot <- T.m >= 26.5
  
  L <- as.radians(WB.params$L)
  I <- sum((T.m / 5)[!T.m.frz]^1.514)
  a <- 6.75e-7 * I^3 - 7.71e-5 * I^2 + 1.79e-2 * I + 0.49
  
  PET.m          <- 16 * (10 * T.m / I)^a                            #    0-26.5 degrees
  PET.m[T.m.frz] <- 0                                                #  < 0 degrees
  PET.m[T.m.hot] <- (-415.85 + 32.24 * T.m - 0.43 * T.m^2)[T.m.hot]  # >= 26.5 degrees
  
  #cat("\nPET.Thorn: ",PET.m[7])
  
  return(PET.m * (monthlengths() / 30) * (daylength(L) / 12))
}

PET.Ding <- function(WB.params)
{
  T.m <- WB.params$T.m
  L   <- as.radians(WB.params$L)
  
  # saturation vapor pressure
  Vsp <- as.numeric(0.611 * exp((17.3 * T.m) / (T.m + 237.3)))
  
  PET.m <- as.numeric(29.8 * daylength(L) * monthlengths() * (Vsp / (T.m + 273.2)))
  #cat("\nPET.Ding: ",PET.m[7])
  
  return(PET.m)
}

PET.Turc <- function(WB.params)
{
  T.m <- WB.params$T.m
  R.m <- WB.params$R.m
  
  R.m <- (R.m * 0.08598)/monthlengths() # Convert from Wh/sq. m/month to cal/sq. cm/day
  
  RH.term <- ifelse((WB.params$RH.m < 50),(1 + ((50 - WB.params$RH.m) / 70)),1)
  
  PET.m <- 0.013 * (T.m / (T.m + 15)) * (R.m + 50) * RH.term * monthlengths()
  PET.m[T.m < 0] <- 0
  
  #cat("\nPET.Turc2nlhvap: ",PET.m[7])
  
  return(PET.m)
}

# PET.Turc2lhvap <- function(WB.params)
# {
#   T.m <- WB.params$T.m
#   R.m <- (WB.params$R.m * 0.0036) / monthlengths() # convert Wh/sq. m/month to MJ/sq. m/day
#   lambda <- 2.501 - 0.002361*T.m          # latent heat of vaporization; MJ / kg
#   
#   RH.term <- ifelse((WB.params$RH.m < 50),(1 + ((50 - WB.params$RH.m) / 70)),1)
#   
#   PET.m <- 0.013 * (T.m / (T.m + 15)) * ((23.8856*R.m + 50)/lambda) * RH.term * monthlengths()
#   PET.m[PET.m < 0] <- 0
#   
#   #cat("\nPET.Turc2lhvap: ",PET.m[7])
#   
#   return(PET.m)
# }

#PET.PriestleyTaylor <- function(WB.params) {
#  T.m <- WB.params$T.m
#  R.m <- (WB.params$R.m * 0.0036) / monthlengths() # convert Wh/sq. m/month to MJ/sq. m/day
#	
#	lambda <- 2.45                 # latent heat of vaporization; MJ / kg
#	s      <- SatVapPresSlope(T.m) # slope of saturation vapor pressure; kPa / C
#	gamma  <- 0.066                # psychrometric constant; kPa / C
#	alpha  <- 1.26                 # Priestley-Taylor coefficient; unitless
#	PET.m  <- (1/lambda) * ((s * R.m) / (s + gamma)) * alpha * monthlengths()
#	PET.m[PET.m < 0] <- 0
#	
#	# Note that when all of the units cancel, you're left with kg water/m^2.
#	# Per the density of water, this is equivalent to mm of water.
#	return(PET.m)
#}

#Makkinik-Hansen, as presented in Cristea et al. 2013
PET.MH <- function(WB.params) {
  T.m <- WB.params$T.m
  R.m <- (WB.params$R.m * 0.0036) / monthlengths() # convert Wh/sq. m/month to MJ/sq. m/day
  press <- ElevToPress(WB.params$E)
  
  lambda <- 2.501 - 0.002361*T.m          # latent heat of vaporization; MJ / kg
  gamma  <- 0.00163*press/lambda          
  delta <- SatVapPresSlope(T.m)
  PET.m  <- (0.7/lambda) * (delta/(delta+gamma)) * R.m * monthlengths()
  PET.m[PET.m < 0] <- 0
  
  #cat("\nPET.MH: ",PET.m[7])
  return(PET.m)
}

#Makkinik 1957
PET.Mak <- function(WB.params) {
  T.m <- WB.params$T.m
  R.m <- (WB.params$R.m * 0.0036) / monthlengths() # convert Wh/sq. m/month to MJ/sq. m/day
  press <- ElevToPress(WB.params$E)
  
  lambda <- 2.501 - 0.002361*T.m          # latent heat of vaporization; MJ / kg
  gamma  <- 0.00163*press/lambda          
  delta <- SatVapPresSlope(T.m)
  PET.m  <- ((0.61/lambda) * (delta/(delta+gamma)) * R.m - 0.12) * monthlengths() 
  PET.m[PET.m < 0] <- 0
  
  #cat("\nPET.Mak: ",PET.m[7])
  return(PET.m)  
}

#Makkinik-Hansen, as adjusted by Cristea at al. 2013 assuming 2.25m/s wind
PET.MHadj <- function(WB.params) {
  
  T.m <- WB.params$T.m
  R.m <- (WB.params$R.m * 0.0036) / monthlengths() # convert Wh/sq. m/month to MJ/sq. m/day
  press <- ElevToPress(WB.params$E)
  RH <- WB.params$RH.m
  wind <- 2.25 # 2m wind speed; calc from NREL map sierra average, converted to 2m with Cristea formula
  Cadj <- 1.036 - 0.527*(RH/100) + 0.041 * wind
  lambda <- 2.501 - 0.002361*T.m          # latent heat of vaporization; MJ / kg
  gamma  <- 0.00163*press/lambda          
  delta <- SatVapPresSlope(T.m)
  PET.m  <- (Cadj/lambda) * (delta/(delta+gamma)) * R.m * monthlengths()
  PET.m[PET.m < 0] <- 0
  
  #cat("\nPET.MHadj: ",PET.m[7])
  return(PET.m)  
}  

#Valiantzas 2013
PET.Vali <- function(WB.params) {
  R.m <- (WB.params$R.m * 0.0036) / monthlengths() # convert Wh/sq. m/month to MJ/sq. m/day
  T.m <- WB.params$T.m
  RH.m <- WB.params$RH.m
  Lat.rad <- as.radians(WB.params$L)
  a <- ifelse(T.m <= -9.5,0.01,T.m+9.5) # do not let the term being rooted be negative
  PET.m <- (0.0393*R.m*sqrt(a)-0.19*(R.m^0.6)*(Lat.rad^0.15)+0.078*(T.m+20)*(1-(RH.m/100)))*monthlengths()
  PET.m[PET.m < 0] <- 0
  #cat("\nPET.Vali: ",PET.m[7])
  return(PET.m)
}

#Copais
PET.Cop <- function(WB.params) {
  
  R.m <- (WB.params$R.m * 0.0036) / monthlengths() # convert Wh/sq. m/month to MJ/sq. m/day
  T.m <- WB.params$T.m
  RH.m <- WB.params$RH.m
  
  m1 <- 0.057
  m2 <- 0.277
  m3 <- 0.643
  m4 <- 0.0124
  
  C1 <- 0.6416 - 0.00784*RH.m + 0.372*R.m - 0.00264*R.m*RH.m
  C2 <- -0.0033 + 0.00812*T.m + 0.0101*R.m + 0.00584*R.m*T.m
  
  PET.m <- (m1 + m2*C2 + m3*C1 + m4*C1*C2) * monthlengths()
  
  PET.m[PET.m < 0] <- 0
  
  #cat("\nPET.cop: ",PET.m[7])
  return(PET.m)
}

PET.Har <- function(WB.params) {
  R.m <- (WB.params$R.m * 0.0036) / monthlengths() # convert Wh/sq. m/month to MJ/sq. m/day
  T.m <- WB.params$T.m
  
  PET.m <- 0.0135 * 0.408*R.m*(T.m+17.8) * monthlengths()
  
  PET.m[PET.m < 0] <- 0
  
  #cat("\nPET.har: ",PET.m[7])
  return(PET.m)
}

PET.BCM <- function(WB.params) {
  return(WB.params$PET.BCM.m)
}

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
  
  lambda <- 2.45
  s <- SatVapPresSlope(T.m)
  press <- ElevToPress(WB.params$E)
  gamma  <- 0.00163*press/lambda  
  alpha <- 1.26
  PET.m <- (1/lambda) * ((s * Rn.m) / (s+gamma)) * alpha * monthlengths()
  PET.m[PET.m < 0] <- 0
  
  return(PET.m)
}



