source("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/4 - simulation/getEquallyLogSpacedIntervals.R")
library(ghyp)
library(MixGHD)

filepath <- "D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/Dataset_t_sampled.RDS"

x <- readRDS(filepath)
x <- x[x$react.ID %in% c(393,394,395, 400),]

prepareFlou <- function(flou){
  flou <- flou[!is.na(flou)]
  flou <- as.numeric(as.character((flou)))
  flou <- sort(flou)
  return(flou)
}

gh <- list()
scaleConsts <- list()
for(i in 1:4){
  print(i)
  vecFluo <- prepareFlou(t(x[i,-c(1:11)]))
  plot(density(vecFluo))

  gh[[i]] <- MGHD(vecFluo, G=2, scale = TRUE)
  scaleConsts[[i]] <- list(m = mean(vecFluo), s = sd(vecFluo))
  coef(gh[[i]])
}


vecFlous <- sapply(1:4, function(i) prepareFlou(t(x[i,-c(1:11)])))

scaleConsts <- list(
  list(m=917.2365,s=258.3078),
  list(m=1142.409,s=465.2219),
  list(m=1660.784,s=914.1312),
  list(m=3557.501,s=2753.308)
)
n <- list(
  list(negs=0.4984655*15890, pos=0.5015345*15890),
  list(negs=0.4647437*14394, pos=0.5352563*14394),
  list(negs=0.4510775*11948, pos=0.5489225*11948),
  list(negs=0.5272043*14810, pos=0.4727957*14810)
)

# first y_ntc1 uses default max.iter and eps 
y_ntc1 <- rGHD(n[[1]][['negs']],1, mu=-1.031838,alpha=0.3236473,sigma=0.06334825,omega=1.065034,lambda=-1.614129)
y_ntc2 <- rGHD(n[[2]][['negs']],1, mu=-1.050779,alpha=0.1374008,sigma=0.0188096,omega=1.70541,lambda=-2.777454)
y_ntc3 <- rGHD(n[[3]][['negs']],1, mu=-1.094076,alpha=0.06105236,sigma=0.004694025,omega=2.092005,lambda=-2.988504)
y_ntc4 <- rGHD(n[[4]][['negs']],1, mu=-1.057021,alpha=0.02366937,sigma=0.00073866,omega=1.619113,lambda=-3.120553)

y_tar1 <- rGHD(n[[1]][['pos']],1, mu=1.213399-0.3,alpha=-0.4315404,sigma=0.4942149,omega=3.719954,lambda=-1.740423)
y_tar2 <- rGHD(n[[2]][['pos']],1, mu=1.178894-0.3,alpha=-0.3605672,sigma=0.2134476,omega=1.4123658,lambda=-0.7825815)
y_tar3 <- rGHD(n[[3]][['pos']],1, mu=1.019818-0.3,alpha=-0.132391,sigma=0.07584152,omega=0.4895525,lambda=-0.4393224)
y_tar4 <- rGHD(n[[4]][['pos']],1, mu=0.977734-0.3,alpha=-0.09873485,sigma=0.02468989,omega=0.3748432,lambda=-1.3308905)

hist(rGHD(7869, 1, mu=1.019818-0.3,alpha=-0.132391,sigma=0.07584152,omega=0.4895525,lambda=-0.4393224))

y1 <- (c(y_ntc1, y_tar1)*scaleConsts[[1]][['s']])+scaleConsts[[1]][['m']]
y2 <- (c(y_ntc2, y_tar2)*scaleConsts[[2]][['s']])+scaleConsts[[2]][['m']]
y3 <- (c(y_ntc3, y_tar3)*scaleConsts[[3]][['s']])+scaleConsts[[3]][['m']]
y4 <- (c(y_ntc4, y_tar4)*scaleConsts[[4]][['s']])+scaleConsts[[4]][['m']]

# heavy rain
plot(density(y1), ylim=c(0,0.0055), xlim=c(500, 1500))  
plot(density(vecFlous[[1]]), ylim=c(0,0.0055), xlim=c(500, 1500))
# mod rain
plot(density(y2), ylim=c(0,0.0035), xlim=c(300, 2000))
plot(density(vecFlous[[2]]), ylim=c(0,0.0035), , xlim=c(300, 2000))
# low rain
plot(density(y3), ylim=c(0,0.002), xlim=c(0, 3500))
plot(density(vecFlous[[3]]),ylim=c(0,0.002), xlim=c(0, 3500))
# no rain
plot(density(y4), ylim=c(0,0.0007), xlim=c(0, 7500))
plot(density(vecFlous[[4]]), ylim=c(0,0.0007), xlim=c(0, 7500))

g <- gh[[1]]
for(g in gh){
  print(g@gpar[[1]])
}
gh[[4]]@gpar[[2]]
