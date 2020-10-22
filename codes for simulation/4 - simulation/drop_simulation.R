source("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/4 - simulation/getEquallyLogSpacedIntervals.R")
library(ghyp)

simulateNTC <- function(R, n=20000){
  param <- list(
    list(mu=-1.057021, sigma=0.00073866,  alpha=0.02366937, omega=1.619113, lambda=-3.120553),  # 1 - No rain
    list(mu=-1.094076, sigma=0.004694025, alpha=0.06105236, omega=2.092005, lambda=-2.988504),  # 2 - Lo rain
    list(mu=-1.050779, sigma=0.0188096,   alpha=0.1374008,  omega=1.70541,  lambda=-2.777454),  # 3 - Mo rain
    list(mu=-1.031838, sigma=0.06334825,  alpha=0.3236473,  omega=1.065034, lambda=-1.614129)   # 4 - Hi rain
  )
  ntc <- rGHD(n, 1, 
              mu     = param[[R]][['mu']],
              alpha  = param[[R]][['alpha']],
              sigma  = param[[R]][['sigma']],
              omega  = param[[R]][['omega']],
              lambda = param[[R]][['lambda']])
  ntc <- (ntc*scaleConsts[[R]][['s']])+scaleConsts[[R]][['m']]
  return(ntc)
}

# MAIN
simulateTarget <- function(conc, R, ntot=20000){   
  # Arguments              | values
  # conc : Concentration   | 0.1, 0.4, 1, 2.5
  # R    : Rain settings   | 1, 2, 3, 4
  nneg <- getNneg(conc, ntot)
  npos <- ntot - nneg
  print(c(ntot, nneg, npos))
  
  # Step 1 : Get f0
  x_f0 <- simulateNTC(R, nneg)
  
  param <- list(
    list(mu=0.977734-1,      sigma=0.02468989,  alpha=-0.09873485, omega=0.3748432,  lambda=-1.3308905),  # 1 - Target
    list(mu=1.019818-1.5,  sigma=0.07584152,  alpha=-0.132391,   omega=0.4895525,  lambda=-0.4393224),  # 2 - Target
    list(mu=1.178894-1, sigma=0.2134476,   alpha=-0.3605672,  omega=1.4123658,  lambda=-0.7825815),  # 3 - Target
    list(mu=1.213399-1,  sigma=0.4942149,   alpha=-0.4315404,  omega=3.719954,   lambda=-1.740423)    # 4 - Target
  )
  
  # Step 2 : Get f2
  x_f1 <- c()
  repeat{
    n_remaining <- npos - length(x_f1)
    z_new <- rGHD(n_remaining, 1, mu=param[[R]][['mu']], alpha=param[[R]][['alpha']], sigma=param[[R]][['sigma']], omega=param[[R]][['omega']], lambda=param[[R]][['lambda']])
    x_new <- (z_new*scaleConsts[[R]][['s']])+scaleConsts[[R]][['m']]
    x_new <- x_new[which(x_new >= median(x_f0))]
    x_f1 <- c(x_f1, x_new)
    if(length(x_f1) == npos){
      break
    }
  }
  
  # Step 3 : append f0 and f1
  x_flou <- c(x_f0, x_f1)
  # hist(x_flou, n=50)
  
  invisible(x_flou)
}

concs <- getEquallyLogSpacedIntervals(0.1, 2.5, 5)

scaleConsts <- list(
  list(m=917.2365,s=258.3078),
  list(m=1142.409,s=465.2219),
  list(m=1660.784,s=914.1312),
  list(m=3557.501,s=2753.308)
)

setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/4 - simulation/simulated")
set.seed(12345)
for(conc in concs){
  for(r in 1:4){        # Rain settings
    for(i in 1:5){      # Replicates
      Targets <- lapply(1:3, function(x) simulateTarget(conc = conc, R = r))
      NTCs <- lapply(1:3, function(x) simulateNTC(R = r))

      simulated <- data.frame(c(Targets, NTCs))
      colnames(simulated) <- c(paste0("Target",1:3), paste0("NTC",1:3))
      
      label_conc <- LETTERS[which(conc == concs)]
            
      filename <- paste0("sim_conc",label_conc,"_R",r,"_rep",i,".csv")
      write.csv(simulated, file = filename, row.names = FALSE)
      # break
    }
    # break
  }
  # break
}

XXX <- 4
simulateTarget(conc = concs[1], R = XXX)
plot(sample(simulateTarget(conc = concs[1], R = XXX)))
simulateTarget(conc = concs[2], R = XXX)
plot(sample(simulateTarget(conc = concs[2], R = XXX)))
simulateTarget(conc = concs[3], R = XXX)
plot(sample(simulateTarget(conc = concs[3], R = XXX)))
simulateTarget(conc = concs[4], R = XXX)
plot(sample(simulateTarget(conc = concs[4], R = XXX)))
simulateTarget(conc = concs[5], R = XXX)
plot(sample(simulateTarget(conc = concs[5], R = XXX)))

