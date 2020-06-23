#
# AET_methods:
#     Various ways of calculating AET given monthly PET and WB.params (see Point_WB.R for details).
#     Not all methods require all elements of WB.params to be populated. See documentation for details.
#
#     Each method returns a 12-element vector of monthly AET values.
#

# Test values
#T.m   <- c(-0.8,-0.6,0.7,3.1,7.1,12.8,17.4,18.1,13.5,9.4,4.7,-2.3)
#P.m   <- c(394,201,33,118,51,2,10,2,40,29,44,14)
#P.m   <- c(94,33,20,18,51,2,8,2,10,15,34,34)
#S.max <- 180
#PET.m <- c(13.79687, 14.52889, 19.75247, 25.19059, 36.73518, 51.10025, 65.13021, 59.57761, 42.96944, 28.98110, 17.42535, 13.51610)
#PET.m <- c(0.000000, 0.000000, 4.355127, 19.12804, 46.87502, 82.57885, 112.4007, 109.3416, 72.94102,  48.0792,  21.6015, 0.000000)

#source("utils.R")

AET.Thorn <- function(PET.m, WB.params, ret.soil.water=FALSE)
{
  
  
  
  P.m   <- WB.params$P.m
  S.max <- WB.params$S.max
  
  P.m     <- c(NA, P.m, P.m[1])
  PET.m   <- c(NA, PET.m, PET.m[1])
  Storage <- c(S.max, rep(NA, 13))
  AET.m   <- rep(NA, 14)
  
  cnt <- 0
  while(any(is.na(Storage)) | Storage[2] != Storage[14])
  {
    
    for(i in 2:14)
    {
      Avail.H2O <- Storage[i - 1] + P.m[i]
      
      
      Def <- PET.m[i] - Avail.H2O
      if(Def < 0)
      {
        Storage[i] <- Storage[i - 1] + P.m[i] - PET.m[i]
        if(Storage[i] > S.max)
          Storage[i] <- S.max
        Def <- 0
      }
      else
        Storage[i] <- 0
      AET.m[i] <- PET.m[i] - Def
    }
    Storage[1] <- Storage[13]
    cnt <- cnt + 1
    if(cnt > 10) {
      stop("Error in Thornthwaite AET")
    }
  }
  
  if(ret.soil.water == TRUE) {
    return(Storage[2:13])
  } else {  
    return(AET.m[2:13])
  }
}

AET.ThornReg <- function(PET.m, WB.params, ret.soil.water=FALSE) {
  WB.params$S.max <- WB.params$S.max.reg
  return(AET.Thorn(PET.m, WB.params, ret.soil.water))
}

AET.Thorn150mm <- function(PET.m, WB.params, ret.soil.water=FALSE) {
  WB.params$S.max <- 150
  return(AET.Thorn(PET.m, WB.params, ret.soil.water))
}

AET.Ding <- function(PET.m, WB.params, ret.soil.water=FALSE)
{
  T.m   <- WB.params$T.m
  P.m   <- WB.params$P.m
  S.max <- WB.params$S.max
  
  ### Water input, snow, rain, etc.
  ## Monthly melt factor
  MF.m <- rep(0, 12)
  partial.months <- T.m > 0 & T.m < 6
  full.months    <- T.m >= 6
  MF.m[partial.months] <- .167 * T.m[partial.months]
  MF.m[full.months]    <- 1
  
  ## Rain & Snow
  rain.m <- as.numeric(MF.m * P.m)
  snow.m <- as.numeric((1 - MF.m) * P.m)
  
  ## Snowpack & Melt & Water input:  uses August as 0 snowpack month to spin up for 29 months
  # Spin up precip, melt, and snowpack
  p.spin    <- as.numeric(c(P.m[c(8:12)], P.m, P.m))  # create Aug - Dec + 24 more months
  melt.spin <- as.numeric(c(MF.m[c(8:12)], MF.m, MF.m))
  pack.spin <- rep(0, 29)
  
  # First month, August, is left at 0. Run the rest.  
  pack.spin[2:29] <- (1 - melt.spin[2:29])^2 * p.spin[2:29]
  summand <- (1 - melt.spin[2:29])
  for(i in 2:29)
    pack.spin[i] <- pack.spin[i] + summand[i - 1] * pack.spin[i -1]
  
  melt.m  <- MF.m * pack.spin[c(17:28)]
  water.m <- rain.m + melt.m
  pack.m  <- pack.spin[c(18:29)]
  
  ### Soil water balance:  May is saturated month to start spin up for 30 months
  soil.water <- c(S.max, rep(0,31))
  water.spin <- c(water.m[c(5:12)], water.m, water.m)
  PET.spin   <- c(PET.m[c(5:12)], PET.m, PET.m)
  
  for(i in 2:32)
  {
    if(water.spin[i] >= PET.spin[i])
    {
      soil.water[i] <- min(S.max, water.spin[i] - PET.spin[i] + soil.water[i - 1])
    } else # PET.spin[i] > water.spin[i]
    {
      soil.water[i] <- soil.water[i - 1] - (soil.water[i - 1] * (1 - exp(-(PET.spin[i] - water.spin[i]) / S.max)))
    }
  }
  soil.m     <- soil.water[c(21:32)]
  soil.delta <- soil.m[c(12, 1:11)] - soil.m
  soil.delta[soil.delta < 0] <- 0
  
  AET.m <- ifelse((PET.m > water.m) & (T.m > 5), water.m[1:12] + soil.delta[1:12], PET.m[1:12])
  AET.m[AET.m > PET.m] <- PET.m[AET.m > PET.m]
  
  
  if(ret.soil.water == TRUE) {
    return(soil.m)
  } else {
    return(AET.m)    
  }
}

AET.DingReg <- function(PET.m, WB.params, ret.soil.water=FALSE) {
  WB.params$S.max <- WB.params$S.max.reg
  return(AET.Ding(PET.m, WB.params,ret.soil.water))
}

AET.Ding150mm <- function(PET.m, WB.params, ret.soil.water=FALSE) {
  WB.params$S.max <- 150
  return(AET.Ding(PET.m, WB.params,ret.soil.water))
}

AET.Wil <- function(PET.m, WB.params, ret.soil.water=FALSE)
{
  T.m   <- WB.params$T.m
  P.m   <- WB.params$P.m
  S.max <- WB.params$S.max
  
  if(S.max <= 0)
  {
    warning(paste("S.max must be positive to calculate Wilmott AET, but it is",S.max))
    return(rep(-1, 12))
  }
  
  # Set initial guesses for W.s and W
  W.s.0 <- 0
  W.0   <- S.max
  W.s   <- 1
  W     <- 0
  # Tolerance for convergence in mm
  tol   <- 1
  # Counter to keep track of how many times we've run through
  # If it gets too large, we'll abort and assume that something is wrong
  ctr   <- 0
  while(abs(W.s[1] - W.s.0) > tol | abs(W[1] - W.0) > tol)
  {
    # In this method, everything is calculated on a quasi-daily basis
    # So temp, precip, PET, etc. are all expanded to reflect this
    
    T <- unlist(mapply(function(x, y) rep(x, y), T.m, monthlengths()))
    P.s <- P.r <- unlist(mapply(function(x, y) rep(x, y) / y, P.m, monthlengths()))
    P.s[T >= -1] <- 0
    P.r[T <  -1] <- 0
    E.0 <- unlist(mapply(function(x, y) rep(x, y) / y, PET.m, monthlengths()))
    
    # Amount that melts each day
    M <- 2.63 + (2.55 * T) + (0.0912 * T * P.r)
    M[M < 0] <- 0
    
    # Amount of water stored as snow each day
    W.s <- c(W.s.0, rep(NA, 364))
    for(i in 2:365)
    {
      snow.cover <- W.s[i - 1] + P.s[i]
      if(M[i] > snow.cover)
        M[i] <- snow.cover
      
      W.s[i] <- W.s[i - 1] + P.s[i] - M[i]
    }
    
    # Evaporative demand: (+) is recharging, (-) is a demand
    D <- M + P.r - E.0
    
    # Evapotranspiration proportion: evapotranspiration slows down as the soil water becomes depleted
    beta.d <- rep(1, 365)
    if(D[1] < 0)
      beta.d[1] <- 1 - exp(-6.68)
    
    # Amount of water in the soil and surplus each day
    W <- c(W.0, rep(NA, 364))
    S <- rep(0, 365)
    for(i in 2:365)
    {
      if(D[i] < 0)
        beta.d[i] <- 1 - exp(-6.68 * W[i - 1] / S.max)
      W[i] <- W[i - 1] + beta.d[i]*D[i]
      if(W[i] > S.max)
      {
        S[i] <- W[i] - S.max
        W[i] <- S.max
      }
    }
    
    # Amount of water in the soil on the last day of each month
    W.30s <- W[cumsum(monthlengths())]
    # Difference between months
    delta.W <- c(W.30s[1] - W.30s[12], W.30s[2:12] - W.30s[1:11])
    
    # Dec 31 is the new initial condition for Jan 1
    W.s.0 <- W.s[365]
    W.0   <- W[365]
    ctr   <- ctr + 1
    if(ctr > 10)
      stop("Error in Wilmott AET")
  }
  
  monthends <- cumsum(monthlengths())
  monthstarts <- c(1, 1 + monthends)[1:12]
  monthdays <- mapply(function(x, y) x:y, monthstarts, monthends)
  
  # Reduce daily variables to monthly variables
  P.r <- sapply(monthdays, function(x) sum(P.r[x]))
  M   <- sapply(monthdays, function(x) sum(M[x]))
  S   <- sapply(monthdays, function(x) sum(S[x]))
  
  AET.m <- P.r + M - delta.W - S
  AET.m[AET.m > PET.m] <- PET.m[AET.m > PET.m]
  
  if(ret.soil.water==TRUE) {
    return(W.30s)
  } else {
    return(AET.m)
  }
}

AET.WilReg <- function(PET.m, WB.params,ret.soil.water=FALSE) {
  WB.params$S.max <- WB.params$S.max.reg
  return(AET.Wil(PET.m, WB.params, ret.soil.water))
}

AET.Wil150mm <- function(PET.m, WB.params, ret.soil.water=FALSE) {
  WB.params$S.max <- 150
  return(AET.Wil(PET.m, WB.params, ret.soil.water))
}