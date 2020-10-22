getEquallyLogSpacedIntervals <- function(start, end, n){
  distance_log <- (log(end) - log(start))/(n-1)
  x <- numeric(n)
  x[1] <- start
  for(i in 2:n){
    x[i] <- exp(log(x[i-1]) + distance_log) 
  }
  x <- round(x, digits=7)
  return(x)
}

getNneg <- function(conc, ntot=20000){
  nneg <- exp(log(ntot) - conc)
  nneg <- round(nneg)
  return(nneg)
}
