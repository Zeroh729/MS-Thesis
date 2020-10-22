source("D:/~Masters/~ MS-STAT/~THESIS/Code/utils.R")
source("D:/~Masters/~ MS-STAT/~THESIS/Code/common.R")
library(EMMIXskew)
library(dplyr)

DropEmPCR <- function(vecFluo, volDrp, distr, maxGroups=3){
  #' distr = c("mvn", "msn", "mvt", "mst")
  #' 
  res <- em(vecFluo, distr=distr,volDrp = volDrp, maxGroups=maxGroups)
  res_info <- get_resInfo(res, distr)
  return(res_info)
}

get_resInfo <- function(res, distr){
  res_info <- list()
  G <- res$G
  
  if(distr == "mst"){
    estParam <- data.frame(Mu = c(res$em$modpts), Sigma = sqrt(c(res$em$sigma)), Df = res$em$dof, Skew = c(res$em$delta), MixProp = res$em$pro) %>%
      mutate(NegThres = rep(res$negThres, nrow(.))) %>%
      mutate(PosThres = rep(res$posThres, nrow(.))) %>% 
      arrange(Mu)
    
    title <- bquote("Neg~"~t~"("~v==.(round(estParam[1,"Df"],2))~","~mu==.(round(estParam[1,"Mu"],2))~","~sigma==.(round(estParam[1,"Sigma"],2))~","~delta==.(round(estParam[1,"Skew"],2))~")")
    if(G == 2){
      subtitle <- bquote("Pos~"~t~"("~v==.(round(estParam[G,"Df"],2))~","~mu==.(round(estParam[G,"Mu"],2))~","~sigma==.(round(estParam[G,"Sigma"],2))~","~delta==.(round(estParam[G,"Skew"],2))~")")
    }else{
      subtitle <- ""
    }
    if(G == 3){
      res_info$desc <- bquote("Rain~"~t~"("~v==.(round(estParam[2,"Df"],2))~","~mu==.(round(estParam[2,"Mu"],2))~","~sigma==.(round(estParam[2,"Sigma"],2))~","~delta==.(round(estParam[2,"Skew"],2))~")")
    }
  }else if(distr == "mvt"){
    estParam <- data.frame(Mu = c(res$em$modpts), Sigma = sqrt(c(res$em$sigma)), Df = res$em$dof, MixProp = res$em$pro) %>%
      mutate(NegThres = rep(res$negThres, nrow(.))) %>%
      mutate(PosThres = rep(res$posThres, nrow(.))) %>% 
      arrange(Mu)
    
    title <- bquote("Neg~"~t~"("~v==.(round(estParam[1,"Df"],2))~","~mu==.(round(estParam[1,"Mu"],2))~","~sigma==.(round(estParam[1,"Sigma"],2))~")")
    if(G == 2){
      subtitle <- bquote("Pos~"~t~"("~v==.(round(estParam[G,"Df"],2))~","~mu==.(round(estParam[G,"Mu"],2))~","~sigma==.(round(estParam[G,"Sigma"],2))~")")
    }else{
      subtitle <- ""
    }
    if(G == 3){
      res_info$desc <- bquote("Rain~"~t~"("~v==.(round(estParam[2,"Df"],2))~","~mu==.(round(estParam[2,"Mu"],2))~","~sigma==.(round(estParam[2,"Sigma"],2))~")")
    }
  }
  
  res_info$est_parameter <- estParam
  res_info$classification <- res$classification
  res_info$G <- G 
  res_info$emres <- res
  res_info$title <- title
  res_info$subtitle <- subtitle
  return(res_info)
}

em <- function(drp, volDrp, distr, maxGroups=3){
  getClusMemberProb <- function(drp, emres){
    if(length(emres$mu) > 1){
      clusterMem <- emres$tau
      means <- emres$modpts
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
    }else{
      # WARNING : This code assumes that the population detected is the negative distribution
      # If you want it to be dynamic, check if the min or max of x is nearer mu
      mean <- emres$mu[1,1]
      sigma <- emres$sigma[1,1,1]
      df <- emres$dof
      x_mayBePos <- drp[drp > mean]
      if(emres$distr == "mvt"){
        d <- EMMIXskew::ddmvt(matrix(x_mayBePos, ncol=1), n=length(x_mayBePos), p=1, mean=mean, cov=sigma, nu = df)
      }else if(emres$distr == "mst"){
        d <- EMMIXskew::ddmst(matrix(x_mayBePos, ncol=1), n=length(x_mayBePos), p=1, mean=mean, cov=sigma, nu = df, del =  emres$delta[1,1])
      }
      surenegs <- length(drp)-length(x_mayBePos)
      i_pos <- which(d < 0.0000001) + surenegs
      i_neg <- setdiff(1:length(drp), i_pos)
      
      negMemberProb <- c(rep(1, surenegs), d)
      posMemberProb <- 1 - negMemberProb
      
      classification <- character(length(drp))
      classification[i_neg] <- "neg"
      classification[i_pos] <- "pos"
      classification <- factor(classification, levels = c("pos", "neg"))
      return(list(negMemberProb,posMemberProb, classification))
    }
  }
  
  # modes <- getInitMus(drp, maxGroups = maxGroups)
  modes <- getInitMus2(drp)
  G <- length(modes)
  initEMSkew <- list(mu = modes,
                     sigma = array(rep(1000, G), c(1,1,G)),
                     pro = rep(1/G, G),
                     dof = rep(30, G),
                     delta = t(matrix(rep(0,G))))
  
  # Manual - https://cran.r-project.org/web/packages/EMMIXskew/EMMIXskew.pdf; ncov = 3 means general variance (1 & 2 gives me equal variances)
  df_drp <- data.frame(Fluorescence=sort(drp))
  emres <- EmSkew(df_drp, init = initEMSkew, g = G, distr = distr, ncov = 3, debug = FALSE) # MAIN FUNCTION
  
  .g(negMemberProb, posMemberProb, classification) %=% getClusMemberProb(drp, emres)
  nneg <- sum(classification == "neg")
  
  res <- list()
  res$est <- cal_concentration(nneg, nrow(df_drp),volDrp)
  res$member <- list(negProb = negMemberProb, posProb = posMemberProb)
  res$G <- G
  res$em <- emres
  res$classification <- classification
  if(any(classification == "pos")){
    res$negThres <- df_drp[max(which(classification == "neg")),]
    res$posThres <- df_drp[min(which(classification == "pos")),]
  }else{
    res$negThres <- 0
    res$posThres <- 0
  }
  return(res)
}

getInitMus <- function(x, maxGroups = 3){
  z <- hist(x, breaks = 50, plot = FALSE)
  mode_i <- which(diff(sign(diff(c(0, z$counts))))==-2)
  mode_i <- mode_i[z$counts[mode_i] >= 50]
  mode_i <- groupNearModes(z, mode_i)
  
  # If only one population is detected, try less smoothed histogram w/ Scott
  if(length(mode_i) < 2){
    z <- hist(x, breaks = "Scott", plot = FALSE)
    mode_i <- which(diff(sign(diff(c(0, z$counts))))==-2)
    mode_i <- mode_i[z$counts[mode_i] >= 50]
    mode_i <- groupNearModes(z, mode_i)
  }
  
  # If still there is one population, increase number of breaks until 3/4s of the number of droplets
  nbreaks <- length(z$breaks)
  repeat{
    if(length(mode_i) >= 2){
      break
    }
    if(nbreaks >= 200){
      warning("Only 1 population found")
      return(mode_i)
    }
    nbreaks <- nbreaks + 1
    z <- hist(x, breaks = nbreaks, plot = FALSE)
    mode_i <- which(diff(sign(diff(c(0, z$counts))))==-2)
    mode_i <- mode_i[z$counts[mode_i] >= 50]
    mode_i <- groupNearModes(z, mode_i)
  }
  
  # Get top distant maxGroups modes with highest counts
  modes <- z$mids[mode_i]
  if(length(mode_i) > maxGroups){
    cluster_membership <- cutree(hclust(dist(modes)), k = maxGroups)
    modes <- sapply(unique(cluster_membership), function(x){
      clusters_i <- mode_i[cluster_membership == x]
      z$mids[clusters_i[which.max(z$counts[clusters_i])]]
    })
  }
  plot(z)
  abline(v=modes)
  return(modes)
}

findpeaks <- function(vec, bw = 1, x.coo = c(1:length(vec))){  #where bw = is box width, setting the sensitivity of the search
  ###set all vectors to null
  pos.x.max <- NULL ;	pos.y.max <- NULL ;	pos.x.min <- NULL ;	pos.y.min <- NULL
  ###Start of for loop:    we walk down the vector with a window of size "bw"
  for(i in 1:(length(vec)-1)){
    #check if we have reached the end of the vector
    if((i+1+bw)>length(vec)){sup.stop <- length(vec)}else{sup.stop <- i+1+bw}
    #check if we are at beginning of the vector
    if((i-bw) < 1){inf.stop <- 1}else{inf.stop <- i-bw}
    #select window in two parts: values beyond i (superior), and values before i (inferior)
    subset.sup <- vec[(i+1):sup.stop]
    subset.inf <- vec[inf.stop:(i-1)]
    ##############################################################
    #are ALL trailing data smaller than i?
    is.max   <- sum(subset.inf > vec[i]) == 0
    #are ALL leading data smaller than i?
    is.nomin <- sum(subset.sup > vec[i]) == 0
    #are ALL trailing data larger than i?
    no.max   <- sum(subset.inf > vec[i]) == length(subset.inf)
    #are ALL leading data larger than i?
    no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup)
    ##############################################################
    #a maximum is found if  all data before and after i are smaller than i
    if(is.max & is.nomin){
      pos.x.max <- c(pos.x.max, x.coo[i])
      pos.y.max <- c(pos.y.max, vec[i])
    }
    #a maximum is found if  all data before and after i are larger than i
    if(no.max & no.nomin){
      pos.x.min <- c(pos.x.min, x.coo[i])
      pos.y.min <- c(pos.y.min, vec[i])
    }
  }#end of for loop
  ###Output
  return(list("max.X" = pos.x.max, "max.Y" = pos.y.max, "min.X" = pos.x.min, "min.Y" = pos.y.min))
}

groupNearModes <- function(z, mode_i){
  # Make sure modes are not 2 bins apart
  if(length(mode_i) > 2){
    for(i in 1:(length(mode_i))){
      if(i == 1){
        dist_modes <- diff(mode_i)
        clusters <- 1
        cluster_membership <- c()
      }
      cluster_membership <- c(cluster_membership, clusters)
      if(i < length(mode_i) && dist_modes[i] > 2){
        clusters <- clusters + 1
      }
    }
    mode_i <- sapply(unique(cluster_membership), function(x){
      clusters_i <- mode_i[cluster_membership == x]
      clusters_i[which.max(z$counts[clusters_i])]
    })
  }
  return(mode_i)
}

getInitMus2 <- function(x){
  bw <- bw.nrd0(x)
  if(bw < 50){
    bw <- 50
  } 
  krn <- density(x, bw = bw)
  krn <- rbind(krn$y, krn$x)
  
  piek <- findpeaks(krn[1, ], bw=20)$max.X
  piek <- rbind(piek, krn[1, piek])    #add peak heights
  piek <- rbind(piek, krn[2, piek[1,]])#add peak x-locations
  piek <- piek[, -which(piek[2,] < max(piek[2,])/100)] 
  if(is.null(dim(piek))){
    piek <- as.matrix(piek)
  }   
  piek <- piek[3,]
  # png(file.path("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/plot_peaks", paste0(rID, ".png")))
  # plot(density(x, bw = bw), main = paste("React ID",rID))
  # for(p in piek){
  # abline(v=p)
  # }
  # dev.off()
  return(piek)
}

# setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/4 - simulation/simulated")
# vecFluo <- read.csv("amplitude/Simulated_091_Amplitude.csv")
# getInitMus(vecFluo[,1], maxGroups = 3)

# df_orig <- readRDS("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/Dataset_t_sampled.RDS")
# for(r in 161:168){
#   print(r)
#   x <- prepareFlou(df_orig[df_orig$react.ID==r,-c(1:11)])
#   emres <- DropEmPCR(x, volDrp = 0.85, distr = "mst")
#   print(emres$emres$est$lambda[['lambda']])
# }

# df_jones <- readRDS("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/df_jones.RDS")
# for(r in 21:28){
#   x <- prepareFlou(df_jones[r,-c(1:3)])
#   emres <- DropEmPCR(x, volDrp = 0.91, distr = "mst")
#   print(emres$emres$est$lambda[['lambda']])
# }
