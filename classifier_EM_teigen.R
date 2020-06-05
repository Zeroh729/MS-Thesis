# install.packages("teigen")
# install.packages("metRology")
# install.packages("extraDistr")
library(teigen)
library(metRology)
library(extraDistr)
library(dplyr)
library(magrittr)
source("D:/~Masters/~ MS-STAT/~THESIS/Code/utils.R")
source("D:/~Masters/~ MS-STAT/~THESIS/Code/common.R")


emclassifier_teigen <- function(vecFluo, volDrp, crit="BIC"){
  res <- em_teigen(vecFluo, volDrp, crit)
  return(res)
}


em_teigen <- function(drp, volDrp, crit="BIC"){ # "ICL" or "BIC"
  getClusMemberProb <- function(emres){
    means <- emres$parameters$mean
    posClust <- which.max(means)
    negClust <- which.min(means)
    rainClust <- if(length(means)==3) which(means == median(means)) else NA
    
    clusts <- c(posClust, negClust, rainClust)
    names(clusts) <- c("pos", "neg", "rain")
    clusts <- sort(clusts)
    classification <- names(clusts[emres$classification])
    
    print(paste("Negative Mu:", means[negClust], "Positive Mu:", means[posClust]))
    if(!is.na(rainClust))
      print(paste("Rain Mu:", means[rainClust]))
    posMemberProb <- emres$fuzzy[,posClust]
    negMemberProb <- emres$fuzzy[,negClust]
    rainMemberProb <- emres$fuzzy[,rainClust]
    return(list(negMemberProb,posMemberProb, classification))
  }
  
  res <- list()
  emres_orig <- teigen(x = drp, Gs=2:3, scale = FALSE, convstyle = "lop")   # MAIN FUNCTION
  
  if(crit=="ICL"){
    emres <- emres_orig$iclresults
  }else{
    emres <- emres_orig
  }
  
  g(negMemberProb, posMemberProb, classification) %=% getClusMemberProb(emres)
  nneg <- sum(classification == "neg")
  
  res$est <- cal_concentration(nneg, length(drp), volDrp)
  res$member <- list(negProb = negMemberProb, posProb = posMemberProb)
  res$G <- emres$G
  res$em <- emres
  res$crit <- crit
  res$critScore <- emres[[tolower(crit)]]
  res$bestmodel <- paste0(emres_orig$bestmodel,"\n", emres_orig$iclresults$bestmodel)
  res$classification <- classification
  res$negThres <- drp[max(which(classification == "neg"))]
  res$posThres <- drp[min(which(classification == "pos"))]
  return(res)
}


teigen_dist <- function(x, df, mu, sigma){
  p <- 1
  prob <- c()
  
  for(i in x){
    numerator <- gamma((df + p) / 2)  * det(as.matrix(sigma))^(-1/2)
    denominator <- (pi * df)^(p/2) * gamma(df/2) * (1 + (mahalanobis(i, center = mu, cov=as.matrix(sigma))/df))^((df+p)/2)
    prob <- c(prob,(numerator/denominator))
  }
  return(prob)
}

# x <- emclassifier_teigen(flou, 0.85, crit="BIC")

