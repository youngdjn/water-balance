#
# PET_mods:
#     Various modifications to PET. Pass in monthly PET values and
#     WB.params (as described in Point_WB.R) and get back modified
#     monthly PET values.
#
#     For descriptions of particular modifications, see documentation.
#



PET.mod.STD <- function(PET.m, WB.params)
{
  return(PET.m)
}

PET.mod.cc050 <- function(PET.m, WB.params)
{
  return(PET.m  * .5)
}

PET.mod.cc025 <- function(PET.m, WB.params)
{
  return(PET.m  * 0.25)
}

PET.mod.SRI <- function(PET.m, WB.params)
{
  L   <- WB.params$L
  R.m <- WB.params$R.m
  R.mean.m <- WB.params$R.mean.m
  
  return(PET.m * (R.m / R.mean.m))
}

PET.mod.HLI <- function(PET.m, WB.params)
{
  
  Sl <- as.radians(WB.params$Sl)
  A  <- as.radians(abs(180 - abs(WB.params$A - 225)))
  L  <- as.radians(WB.params$L)
  
  # McCune & Keon 2002
  HLI <- 0.339 + 0.808 * (cos(L) * cos(Sl)) - 0.196 * (sin(L) * sin(Sl)) - 0.482 * (cos(A) * sin(Sl))
  return(PET.m * HLI)
}

#This one standardizes within each month
PET.mod.DL1 <- function(PET.m, WB.params)
{
  R.hr.m <- WB.params$R.hr.m
  L <- as.radians(WB.params$L)
  
  return(PET.m * (R.hr.m / (daylength(L) * monthlengths())))
}

#This one standardizes to the year
PET.mod.DL2 <- function(PET.m, WB.params)
{
  R.hr.m <- WB.params$R.hr.m
  L <- as.radians(WB.params$L)
  
  return(PET.m * (R.hr.m / (12 * monthlengths())))
}

#Removed this modifier because as originally conceptualized it was identical to SRI
#PET.mod.DL2 <- function(PET.m, WB.params)
#{
#
#	L <- as.radians(WB.params$L)
#	PET.m.SRI <- PET.mod.SRI(PET.m, WB.params)
#	
#	return(PET.m.SRI * (daylength(L) / 12))
#}
