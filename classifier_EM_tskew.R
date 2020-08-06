setwd("D:/~Masters/~ MS-STAT/~THESIS/Code")
source("utils.R")
source("common.R")
library(EMCluster)
library(EMMIXskew)
library(dplyr)

emclassifier_tskew <- function(vecFluo, volDrp, crit="BIC"){
  res <- em_tskew(data.frame(Fluorescence=sort(vecFluo)), distr="mst",volDrp = volDrp, crit)
  res_info <- get_resInfo(res)
  return(res_info)
}

emclassifier_t <- function(vecFluo, volDrp, crit="BIC"){
  res <- em_tskew(data.frame(Fluorescence=sort(vecFluo)), distr = "mvt", volDrp = volDrp, crit)
  res_info <- get_resInfo(res)
  return(res_info)
}

get_resInfo <- function(res){
  res_info <- list()
  estParam <- data.frame(Mu = c(res$em$mu), Sigma = sqrt(c(res$em$sigma)), Df = res$em$dof, Skew = c(res$em$delta), MixProp = res$em$pro) %>%
    mutate(NegThres = rep(res$negThres, nrow(.))) %>%
    mutate(PosThres = rep(res$posThres, nrow(.))) %>% 
    arrange(Mu)
  res_info$est_parameter <- estParam
  res_info$classification <- res$classification
  res_info$G <- G <- res$G
  res_info$emres <- res
  
  title <- bquote("Neg~"~t~"("~v==.(round(estParam[1,"Df"],2))~","~mu==.(round(estParam[1,"Mu"],2))~","~sigma==.(round(estParam[1,"Sigma"],2))~","~delta==.(round(estParam[1,"Skew"],2))~")")
  subtitle <- bquote("Pos~"~t~"("~v==.(round(estParam[G,"Df"],2))~","~mu==.(round(estParam[G,"Mu"],2))~","~sigma==.(round(estParam[G,"Sigma"],2))~","~delta==.(round(estParam[G,"Skew"],2))~")")
  if(G == 3){
    res_info$desc <- bquote("Rain~"~t~"("~v==.(round(estParam[2,"Df"],2))~","~mu==.(round(estParam[2,"Mu"],2))~","~sigma==.(round(estParam[2,"Sigma"],2))~","~delta==.(round(estParam[2,"Skew"],2))~")")
  }
  
  res_info$title <- title
  res_info$subtitle <- subtitle
  return(res_info)
}

em_tskew <- function(drp, volDrp, distr, crit="BIC"){
  getClusMemberProb <- function(emres){
    clusterMem <- emres$tau
    means <- emres$mu
    posClust <- which.max(means)
    negClust <- which.min(means)
    rainClust <- if(length(means)==3) which(means == median(means)) else NA
    
    clusts <- c(posClust, negClust, rainClust)
    names(clusts) <- c("pos", "neg", "rain")
    clusts <- sort(clusts)
    classification <- factor(names(clusts[emres$clust]), levels = c("pos", "neg",  "rain"))
    
    print(paste("Negative Mu:", means[negClust], "Positive Mu:", means[posClust]))
    if(!is.na(rainClust))
      print(paste("Rain Mu:", means[rainClust]))
    posMemberProb <- clusterMem[,posClust]
    negMemberProb <- clusterMem[,negClust]
    rainMemberProb <- clusterMem[,rainClust]
    return(list(negMemberProb,posMemberProb, classification))
  }
  
  getInitParams <- function(drp, G){
    set.seed(1234)
    initEM <- init.EM(drp, nclass = G)
    initEMSkew <- list(mu = t(initEM$Mu),
                       sigma = array(c(initEM$LTSigma), c(1,1,G)),
                       pro = initEM$pi,
                       dof = rep(30, G),
                       delta = t(matrix(rep(0,G))))
    writeLines(paste("Initial parameters for G=", G))
    writeLines(paste("Mu = ", initEMSkew$mu))
    writeLines(paste("Sigma = ", initEMSkew$sigma))
    writeLines(paste("Pi = ", initEMSkew$pro))
    return(initEMSkew)
  }
  
  set.seed(1234)
  # Manual - https://cran.r-project.org/web/packages/EMMIXskew/EMMIXskew.pdf; ncov = 3 means general variance (1 & 2 gives me equal variances)
  emres_G2 <- EmSkew(drp, init = getInitParams(drp, G = 2), g = 2, nkmeans = 2, nrandom = 2, distr = distr, ncov = 3, initloop = 20, debug = FALSE) # MAIN FUNCTION
  emres_G3 <- EmSkew(drp, init = getInitParams(drp, G = 3), g = 3, nkmeans = 3, nrandom = 3, distr = distr, ncov = 3, initloop = 20, debug = FALSE) # MAIN FUNCTION
  
  if(crit == "BIC"){
    scoreG2 <- emres_G2$bic
    scoreG3 <- emres_G3$bic
  }else if(crit == "ICL"){
    scoreG2 <- emres_G2$ICL
    scoreG3 <- emres_G3$ICL
  }else if(crit == "AIC"){
    scoreG2 <- emres_G2$aic
    scoreG3 <- emres_G3$aic
  }
  
  if(scoreG3 < scoreG2){
    emres <- emres_G3
    emres_lower <- emres_G2
    G <- 3
  }else{
    emres <- emres_G2
    emres_lower <- emres_G3
    G <- 2
  }
  
  g(negMemberProb, posMemberProb, classification) %=% getClusMemberProb(emres)
  nneg <- sum(classification == "neg")
  
  res <- list()
  res$est <- cal_concentration(nneg, nrow(drp),volDrp)
  res$member <- list(negProb = negMemberProb, posProb = posMemberProb)
  res$G <- G
  res$em <- emres
  res$em_lower <- emres_lower
  res$crit <- crit
  res$critScore <- min(scoreG2, scoreG3)
  res$bestmodel <- paste(paste(crit, "of G=2 is",scoreG2), paste(crit, "of G=3 is", scoreG3), sep = "\n")
  res$classification <- classification
  res$negThres <- drp[max(which(classification == "neg")),]
  res$posThres <- drp[min(which(classification == "pos")),]
  return(res)
}

# setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")
# if(!exists("df_orig")) df_orig <- readRDS("Dataset_t_sampled.RDS")

# flou <- df_orig %>% filter(react.ID==13) %>% select(-c(1:11))
# flou <- flou[!is.na(flou)]
# flou <- sort(as.numeric(as.character((flou))))
# 
# emres <- emclassifier_tskew(flou, 0.85, crit="BIC")
# emres$emres$bestmodel
# emres$G
# emres$emres$em$mu
# 
# 
# emres$emres$em$aic
# emres$emres$em_lower$aic

# Comparison with 3 Targets in Umbrella tutorial data
#           Umbrella  (by thres) |  (by robust)  | Mine
# Target 1 :  6021  |  6008  |  6021
# Target 2 :  5893  |  5914  |  5111
# Target 3 :  5717  |  5697  |  4584


