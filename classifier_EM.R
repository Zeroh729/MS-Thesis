source("D:/~Masters/~ MS-STAT/~THESIS/Code/utils.R")
source("D:/~Masters/~ MS-STAT/~THESIS/Code/common.R")

emclassifier <- function(vecFluo, volDrp, crit="BIC"){
  res <- em(data.frame(Fluorescence=sort(vecFluo)), volDrp, crit=crit)
  return(res)
}


## EM STARTS HERE
library(EMCluster)
em <- function(drp, volDrp, crit="BIC"){
  getClusMemberProb <- function(emres){
    clusterMem <- e.step(drp, emres, norm = T)$Gamma
    means <- emres$Mu
    posClust <- which.max(means)
    negClust <- which.min(means)
    rainClust <- if(length(means)==3) which(means == median(means)) else NA
    
    clusts <- c(posClust, negClust, rainClust)
    names(clusts) <- c("pos", "neg", "rain")
    clusts <- sort(clusts)
    classification <- names(clusts[emres$class])
    
    print(paste("Negative Mu:", means[negClust], "Positive Mu:", means[posClust]))
    if(!is.na(rainClust))
      print(paste("Rain Mu:", means[rainClust]))
    posMemberProb <- clusterMem[,posClust]
    negMemberProb <- clusterMem[,negClust]
    rainMemberProb <- clusterMem[,rainClust]
    return(list(negMemberProb,posMemberProb, classification))
  }
  
  set.seed(1234)
  emres_G2 <- emcluster(drp, init.EM(drp, nclass = 2), assign.class = TRUE) # MAIN FUNCTION
  emres_G3 <- emcluster(drp, init.EM(drp, nclass = 3), assign.class = TRUE) # MAIN FUNCTION

  if(crit == "BIC"){
    scoreG2 <- em.bic(drp, emres_G2)
    scoreG3 <- em.bic(drp, emres_G3)
  }else if(crit == "ICL"){
    scoreG2 <- em.icl(drp, emres_G2)
    scoreG3 <- em.icl(drp, emres_G3) 
  }
  
  if(scoreG3 < scoreG2){
    emres <- emres_G3
    G <- 3
  }else{
    emres <- emres_G2
    G <- 2
  }
  
  g(negMemberProb, posMemberProb, classification) %=% getClusMemberProb(emres)
  nneg <- sum(classification == "neg")
  
  res <- list()
  res$est <- cal_concentration(nneg, nrow(drp),volDrp)
  res$member <- list(negProb = negMemberProb, posProb = posMemberProb)
  res$G <- G
  res$em <- emres
  res$crit <- crit
  res$critScore <- min(scoreG2, scoreG3)
  res$bestmodel <- paste(paste(crit, "of G=2 is",scoreG2), paste(crit, "of G=3 is", scoreG3), sep = "\n")
  res$classification <- classification
  res$negThres <- drp[max(which(classification == "neg")),]
  res$posThres <- drp[min(which(classification == "pos")),]
  return(res)
}

# setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/Dataset - one csv per reactionID")
# drp <- read.csv("7_cru_338.csv")
# res <- em(drp, 0.85)
# thres <- res$thres

# funcShaded <- function(x, clus) {
#   y <- dnorm(x, mean = res$em$Mu[clus,], sd = sqrt(res$em$LTSigma[clus,]))
#   return(y)
# }
# ggplot(drp, aes(x=Fluorescence)) +
#   geom_histogram(aes(y =..density..),color="gray64", fill="gray88") +
#   geom_vline(xintercept = thres)+
#   stat_function(fun = dnorm, colour = "indianred1", args = list(res$em$Mu[1,], sqrt(res$em$LTSigma[1,]))) +
#   stat_function(fun = dnorm, colour = "green3", args = list(res$em$Mu[2,], sqrt(res$em$LTSigma[2,]))) +
#   stat_function(fun = funcShaded, geom="area", fill="indianred1", alpha=0.3, args=list(clus=1)) +
#   stat_function(fun = funcShaded, geom="area", fill="green3", alpha=0.3, args=list(clus=2)) +
#   labs(title="", x="", y = "") +
#   theme_minimal()

## TEST function
# drp <- read.csv("6_M88017_308.csv")
# 
library(magrittr)
# setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")
# if(!exists("df_orig")) df_orig <- read.csv("Dataset_t.csv")
# drp <- as.numeric(df_orig[1,-c(1:11)]) %>% 
#        .[!is.na(.)]
# write.csv(drp, "drp1.csv")
