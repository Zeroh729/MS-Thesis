cal_concentration <- function(nneg, total, volDrp){
  res <- list()
  volSamp <- 20
  
  # print(paste("Negatives",nneg))
  # print(paste("Total",total))
  npos <- total - nneg
  
  lambda <- -log(nneg/total)
  lower <- lambda - 1.96 * sqrt((total - length(nneg))/(total * length(nneg)))
  upper <- lambda + 1.96 * sqrt((total - length(nneg))/(total * length(nneg)))
  res$lambda <- c(lambda, upper, lower)
  names(res$lambda) <- c("lambda", "upper", "lower")
  
  conc <- lambda * volSamp/volDrp * 1000
  conc_lo <- lower * volSamp/volDrp * 1000
  conc_up <- upper * volSamp/volDrp * 1000
  res$conc <- c(conc, conc_up, conc_lo)
  names(res$conc) <- c("conc", "upper", "lower")
  return(res)
}