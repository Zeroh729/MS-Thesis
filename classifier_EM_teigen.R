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


emclassifier_teigen <- function(vecFluo, volDrp){
  res <- em_teigen(vecFluo, volDrp)
  return(res)
}


em_teigen <- function(drp, volDrp){
  getClusMemberProb <- function(emres){
    means <- emres$parameters$mean
    posClust <- which.max(means)
    negClust <- which.min(means)
    rainClust <- if(length(means)==3) which(means == median(means)) else NA
    
    clusts <- c(posClust, negClust, rainClust)
    names(clusts) <- c("pos", "neg", "rain")
    #clusts <- sort(clusts)
    classification <- names(clusts[emres$classification])
    
    print(paste("Negative Mu:", means[negClust], "Positive Mu:", means[posClust]))
    if(!is.na(rainClust))
      print(paste("Rain Mu:", means[rainClust]))
    posMemberProb <- emres$fuzzy[,posClust]
    negMemberProb <- emres$fuzzy[,negClust]
    rainMemberProb <- emres$fuzzy[,rainClust]
    thres <- NA
    return(list(negMemberProb,posMemberProb, thres, classification))
  }
  
  res <- list()
  emres <- teigen(x = sort(drp), Gs=2:3, scale = FALSE, convstyle = "lop")   # MAIN FUNCTION
  
  g(negMemberProb, posMemberProb, thres, classification) %=% getClusMemberProb(emres)
  nneg <- sum(negMemberProb > posMemberProb)
  
  res$est <- cal_concentration(nneg, length(drp), volDrp)
  res$member <- list(negProb = negMemberProb, posProb = posMemberProb)
  res$em <- emres
  res$thres <- thres
  res$classification <- classification
  res$negThres <- max(which(x$classification == "neg"))
  res$posThres <- min(which(x$classification == "pos"))
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

