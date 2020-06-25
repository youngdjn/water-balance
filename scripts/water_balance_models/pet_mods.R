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
