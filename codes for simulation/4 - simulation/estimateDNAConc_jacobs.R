setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/4 - simulation/simulated/")
source("../getEquallyLogSpacedIntervals.R")
source("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Jacobs/umbrella-master/1D/Umbrella_1d_V1.R")
source("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/Cloudy-V2-04_classification.R")
source("../../../classifier_EM.R")
source("../../../classifier_EM_tskew.R")
source("../../../DropEmPCR.R")
library(doParallel)
registerDoParallel(cores = 3)

getConcFromLabel <- function(label){
  concs <- getEquallyLogSpacedIntervals(0.1, 2.5, 5)
  return(concs[which(label == LETTERS)])
}

getMethodFilename <- function(method){
  if(method=="EM"){
    return("EM")
  }else if(method=="EM_t"){
    return("EM_t")
  }else if(method=="EM_tskew"){
    return("EM_tskew")
  }
  return(method)
}

getResultList <- function(drp_classification, n_neg, n_pos, n_rain, filename, rep){
  n_tot <- n_pos + n_rain + n_neg
  estLambda <- -log(n_neg/n_tot)
  metadata <- strsplit(filename, "_")[[1]]
  concLabel <- gsub("conc", "", metadata[2])
  conc <- getConcFromLabel(concLabel)
  
  trueThres <- getNneg(conc)
  # Misclassification in the negative population
  TN <- sum(drp_classification[1:trueThres] == "neg")
  FP <- sum(drp_classification[1:trueThres] != "neg")
  
  # Misclassification in the positive population
  TP <- sum(drp_classification[(trueThres+1):length(drp_classification)] != "neg")
  FN <- sum(drp_classification[(trueThres+1):length(drp_classification)] == "neg")
  
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
    TN = TN, FP = FP, TP = TP, FN = FN,
    estLambda = estLambda
  ))
}

main <- function(method){
  setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/4 - simulation/simulated/")
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
  x <- foreach(f = list.files(pattern = ".csv"), .combine = c,  .export = c("getResultList", "getConcFromLabel")) %:% 
    foreach(t = 1:3) %dopar% {
      setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/4 - simulation/simulated/")
      source("../getEquallyLogSpacedIntervals.R")
      source("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/Cloudy-V2-04_classification.R")
      source("../../../classifier_EM.R")
      source("../../../classifier_EM_tskew.R")
      source("../../../DropEmPCR.R")
      writeLines(paste0("f ", f))
      writeLines(paste0("t ", t))
      target <- read.csv(f)[,t]  # only read columns Target1-3, no need for NTC
      if(method == "EM"){
        emres <- emclassifier(target, volDrp = 0.85, crit = "ICL")
        classification <- emres$classification[order(target)]
      }else if(method == "EM_t"){
        emres <- emclassifier_t(target, volDrp = 0.85, crit = "ICL")
        classification <- emres$classification[order(target)]
      }else if(method == "EM_tskew"){
        emres <- emclassifier_tskew(target, volDrp = 0.85, crit = "ICL")
        classification <- emres$classification[order(target)]
      }else if(method == "DropEmPCR_t"){
        emres <- DropEmPCR(target, volDrp = 0.85, distr = "mvt", maxGroups = 2)
        classification <- emres$classification[order(target)]
      }else if(method == "DropEmPCR_tskew"){
        emres <- DropEmPCR(target, volDrp = 0.85, distr = "mst", maxGroups = 2)
        classification <- emres$classification[order(target)]
      }else if(method == "cloudy"){
        res <- cloudyClassifier(target, showRain = TRUE)
        classification <- res
      }
      
      n_neg <- sum(classification == "neg")
      n_rain <- sum(classification == "rain")
      n_pos <- sum(classification == "pos")
      resultList <- getResultList(drp_classification = classification, n_neg = n_neg, n_pos = n_pos, n_rain = n_rain, filename = f, rep = t)
      return(resultList)
    }
  return(x)
}

mainUmbrella <- function(){
  x <- foreach(f = list.files(pattern = ".csv"), .combine = c,  .export = c("getResultList", "getConcFromLabel")) %dopar% {
    setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/4 - simulation/simulated/")
    source("../getEquallyLogSpacedIntervals.R")
    source("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Jacobs/umbrella-master/1D/Umbrella_1d_V1.R")
    writeLines(paste0("f ", f))
    experiment <- read.csv(f)
    umbrellaRes <- Umbrella1d(experiment,NTC=c("NTC1","NTC2","NTC3"))
    resList <- list()
    for(t in c("Target1", "Target2", "Target3")){
      n_pos <- umbrellaRes[[t]]$tresh["pos","count"]
      n_rain <- umbrellaRes[[t]]$tresh["rain","count"]
      n_neg <- umbrellaRes[[t]]$tresh["neg","count"]
      classification <- ifelse(umbrellaRes[[t]]$droppi >= .8, "neg", "pos")
      
      resList[[(length(resList)+1)]] <- getResultList(drp_classification = classification, n_neg = n_neg, n_pos = n_pos, n_rain = n_rain, filename = f, rep = (length(resList)+1))
    }
    return(resList)
  }
  return(x)
}
# 
# main("cloudy")
main("umbrella")
# main("EM_tskew")
# main("EM_t")
main("DropEmPCR_tskew")
main("DropEmPCR_t")
beepr::beep()