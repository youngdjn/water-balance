monthlengths <- function()
{
  return(c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
}

SatVapPresSlope <- function(temp) {
  sl <- (4098*0.6108*exp(17.27*temp/(temp+237.3)))/((temp+237.3)^2)
  return(sl)
}

#Bernoulli formula
ElevToPress <- function(elev) {
  
  101325*exp(-9.80665*0.0289644*elev/(8.31447*288.15))/1000
  
}