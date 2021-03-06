setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/simulated2")
source("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Jacobs/umbrella-master/1D/Umbrella_1d_V1.R")
source("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/Cloudy-V2-04_classification.R")
source("../classifier_EM.R")
source("../classifier_EM_tskew.R")
library(doParallel)

getConcFromLabel <- function(label){
  if(label == "A") return(0.1)
  if(label == "B") return(0.4)
  if(label == "C") return(1)
  if(label == "D") return(2.5)
}

getMethodFilename <- function(method){
  if(method=="EM"){
    return("EM (ICL)")
  }else if(method=="EM_t"){
    return("EM_t (BIC)")
  }else if(method=="EM_tskew"){
    return("EM_tskew (BIC)")
  }
  return(method)
}

getResultList <- function(n_neg, n_pos, n_rain, filename, rep){
  n_tot <- n_pos + n_rain + n_neg
  estLambda <- -log(n_neg/n_tot)
  metadata <- strsplit(filename, "_")[[1]]
  concLabel <- gsub("conc", "", metadata[2])
  conc <- getConcFromLabel(concLabel)
  rain_setting <- gsub("R", "", metadata[3])
  exp <- gsub("(rep|.csv)", "", metadata[4])
  return(list(
    conc = conc, 
    rain_setting = rain_setting, 
    exp = exp, 
    rep = rep, 
    n_tot = n_tot, 
    n_neg = n_neg, 
    n_rain = n_rain, 
    n_pos = n_pos, 
    estLambda = estLambda
  ))
}

main <- function(method){
  setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/simulated2")
  if(method == "umbrella"){
    x <- mainUmbrella()
  }else{
    x <- mainEM_orCloudy(method)
  }
  y <- data.frame(t(sapply(x, function(z) return(unlist(z)))))
  filename <- paste0("../Estimates_simulated_",getMethodFilename(method),".csv")
  write.csv(y, filename, row.names = FALSE)
  writeLines(paste0("Saved! ", filename))
}

mainEM_orCloudy <- function(method){
  x <- foreach(f = list.files(pattern = ".csv"), .combine = c) %:% 
    foreach(t = 1:3) %dopar% {
      writeLines(paste0("f ", f))
      writeLines(paste0("t ", t))
      target <- read.csv(f)[,t]  # only read columns Target1-3, no need for NTC
      if(method == "EM"){
        emres <- emclassifier(target, volDrp = 0.85, crit = "ICL")
        classification <- emres$classification
      }else if(method == "EM_t"){
        emres <- emclassifier_t(target, volDrp = 0.85, crit = "BIC")
        classification <- emres$classification
      }else if(method == "EM_tskew"){
        emres <- emclassifier_tskew(target, volDrp = 0.85, crit = "BIC")
        classification <- emres$classification
      }else if(method == "cloudy"){
        res <- cloudyClassifier(target, showRain = TRUE)
        classification <- res
      }
      
      n_neg <- sum(classification == "neg")
      n_rain <- sum(classification == "rain")
      n_pos <- sum(classification == "pos")

      return(getResultList(n_neg = n_neg, n_pos = n_pos, n_rain = n_rain, filename = f, rep = t))
    }
  return(x)
}

mainUmbrella <- function(){
  x <- foreach(f = list.files(pattern = ".csv"), .combine = c) %dopar% {
    writeLines(paste0("f ", f))
    experiment <- read.csv(f)
    umbrellaRes <- Umbrella1d(experiment,NTC=c("NTC1","NTC2","NTC3"))
    resList <- list()
    for(t in c("Target1", "Target2", "Target3")){
      n_pos <- umbrellaRes[[t]]$tresh["pos","count"]
      n_rain <- umbrellaRes[[t]]$tresh["rain","count"]
      n_neg <- umbrellaRes[[t]]$tresh["neg","count"]
      
      resList[[(length(resList)+1)]] <- getResultList(n_neg = n_neg, n_pos = n_pos, n_rain = n_rain, filename = f, rep = (length(resList)+1))
    }
    return(resList)
  }
  return(x)
}
# 
# main("cloudy")
# main("umbrella")
# main("EM_tskew")
main("EM_t")
