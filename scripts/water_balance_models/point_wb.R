#
# Point_WB
#	Runs every permutation of selected water balance methods for a point.
#
# WB.params:
#	list of water balance parameters with the elements necessary for the
#	chosen PET/AET methods. Possible items in the list are:
#		T.m		monthly mean temperature (degrees C)
#		P.m		monthly precipitation (mm)
#		R.m		monthly total solar radiation (W hrs per sq. m)
#		R.hr.m	monthly total hours of direct radiation (hr)
#		Sl		slope (degrees)
#		A		aspect (degrees)
#		L		latitude (decimal degrees)
#		S.max		maximum soil water holding at 150 cm (mm)
#
# PET.methods:
#	character vector of PET methods to use. Choices are:
#		"Thornthwaite"	Thornthwaite 1948
#		"Dingman"		Dingman 2002
#		"Turc"		Turc 1961
#
# PET.mods:
#	character vector of PET modifiers to use. Choices are:
#		"STD"		Standard day length modifier
#		"SRI"		Solar radiation index
#		"HLI"		Heat load index
#		"DL1"		Topographically modified day length I: based on R.h instead of day length
#		"DL2"		Topographically modified day length II: Scale standard day length by SRI
#
# AET.methods:
#	character vector of AET methods to use. Choices are:
#		"Thornthwaite"	Thornthwaite 1948
#		"Dingman"		Dingman 2002
#		"Wilmott"		Wilmott 1985
#
# monthly:
#	boolean: if TRUE, report monthly values, if FALSE, report annual summary. Default FALSE.
#
# returns:
#     Big ugly list. Each list item represents one method permutation, and in turn contains three sub-lists:
#     PET, AET, and Deficit. Each sub-list contains values for each point.
#

source("scripts/water_balance_models/pet_methods.R")
source("scripts/water_balance_models/pet_mods.R")
source("scripts/water_balance_models/aet_methods.R")

Point_WB <- function(WB.params, PET.methods, PET.mods, AET.methods, monthly=TRUE)
{
  
  PET <- lapply(PET.methods, function(x)
  {
    #print(paste("Running ","PET", x, sep=" "))
    do.call(paste("PET", x, sep="."),list(WB.params=WB.params))
  })
  names(PET) <- PET.methods
  
  PET <- unlist(lapply(PET, function(pet)
  {
    lapply(PET.mods, function(mod)
    {
      #print(paste("Running ","PET.mods", mod, sep=" "))
      do.call(paste("PET.mod", mod, sep="."), list(PET.m=pet, WB.params=WB.params))
    })
  }), recursive=F)
  names(PET) <- apply(matrix(c(sapply(PET.methods, function(x) rep(x, length(PET.mods))), rep(PET.mods, length(PET.methods))), ncol=2), 1, function(x) paste(x[1], x[2], sep="."))
  
  #print(PET)
  
  AET <- unlist(lapply(PET, function(pet)
  {
    lapply(AET.methods, function(aet)
    {
      #print(paste("Running ","AET", aet, sep=" "))
      do.call(paste("AET", aet, sep="."), list(PET.m=pet, WB.params=WB.params))
    })
  }), recursive=F)
  if(length(AET.methods) == 1)
  {
    names(AET) <- sapply(names(AET), function(x) paste(x, AET.methods, sep="."))
  } else
  {
    names(AET) <- sapply(names(AET), function(x) paste(substr(x, 0, nchar(x) - 1), AET.methods[as.integer(substr(x, nchar(x), nchar(x)))], sep="."))
  }
  
  PET <- unlist(lapply(PET, function(x) rep(list(x), length(AET.methods))), recursive=F)
  names(PET) <- names(AET)
  
  Deficit <- mapply(function(x, y) x - y, PET, AET, SIMPLIFY=F)
  
  if(!monthly)
  {
    PET <- lapply(PET, sum)
    AET <- lapply(AET, sum)
    Deficit <- lapply(Deficit, sum)
  }
  
  return(list(PET=PET, AET=AET, Deficit=Deficit))
}
