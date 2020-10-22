setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")

if(!exists("df_orig")) df_orig <- readRDS("Dataset_t_sampled.RDS")
library(ggplot2)
library(gridExtra)
library(dplyr)
library(magrittr)
library(EMMIXskew)
library(EMCluster)
prepareFlou <- function(flou){
  flou <- flou[!is.na(flou)]
  flou <- as.numeric(as.character((flou)))
  flou <- sort(flou)
  return(flou)
}

sim <- read.csv("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/3 - simulation/simulated3/sim_concD_R1_rep3.csv")
sim <- sim[,"Target2"]
sim <- sim[!is.na(sim)]
sim <- prepareFlou(sim)
findPopns(sim)

findPopns <- function(x, maxGroups = 3){
  z <- hist(x, plot = FALSE)
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
    if(nbreaks >= length(x)*0.75){
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

temp <- function(drp_i){
  x <- df_orig[drp_i,-(1:11)] %>% t()
  x <- x[!is.na(x)]
  x <- prepareFlou(x)
  findPopns(x)
}
for(i in 1:460){
  print(i)
  temp(i)
}
temp(463)
temp(150)

hist(x)

x <- df_orig[463,-(1:11)] %>% t()
x <- x[!is.na(x)]
x <- prepareFlou(x)
hist(x)
z2 <- DropEmPCR(sim, 0.85, "mst")
z2$emres$est

modes <- findPopns(x)
G <- length(modes)
initEMSkew <- list(mu = modes,
                   sigma = array(rep(1000, G), c(1,1,G)),
                   pro = rep(1/G, G),
                   dof = rep(30, G),
                   delta = t(matrix(rep(0,G))))

emres <- EmSkew(x, init = initEMSkew, g = G, distr = "mvt", ncov = 3, debug = FALSE)
-log(table(emres$clust)[[1]]/length(emres$clust))

