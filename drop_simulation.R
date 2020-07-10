library(ghyp)

getNneg <- function(conc, ntot=20000){
  nneg <- exp(log(ntot) - conc)
  nneg <- round(nneg)
  return(nneg)
}

simulateNTC <- function(n=20000){
  model_ntc <- ghyp(mu = 2000,         # 0
                    sigma = 100,      # 1
                    lambda = 0.5,   # 0.5 - kurtosis?
                    chi = 55,        # 0.5 - kurtosis?
                    gamma = 0       # 0 - skew?
  )
  # plot(model_ntc)
  ntc <- rghyp(n, model_ntc)
    
  return(ntc)
}

# MAIN
simulateTarget <- function(conc, R, ntot=20000){   
  # Arguments              | values
  # conc : Concentration   | 0.1, 0.4, 1, 2.5
  # R    : Rain settings   | 1, 2, 3, 4
  
  nneg <- getNneg(conc, ntot)
  npos <- ntot - nneg
  
  # Step 1 : Get f0
  x_f0 <- simulateNTC(nneg)
  
  param <- list(
    list(mu=8000, sigma=250, lambda=0.5, chi=55, gamma=0),  # Rain 1
    list(mu=8000, sigma=350, lambda=0.5, chi=55, gamma=-350),  # Rain 2
    list(mu=3500, sigma=450, lambda=0.5, chi=55, gamma=200),  # Rain 3
    list(mu=2500, sigma=550, lambda=0.5, chi=55, gamma=350)  # Rain 
  )
  
  model_f1 <- ghyp(mu = param[[R]]$mu,   # 0
                   sigma = param[[R]]$sigma,      # 1
                   lambda = param[[R]]$lambda,    # 0.5 - kurtosis?
                   chi = param[[R]]$chi,          # 0.5 - kurtosis?
                   gamma = param[[R]]$gamma       # 0 - skew? (upto -800)
  )
  # plot(model_f1)
  
  # Step 2 : Get f2
  repeat{
    x_f1 <- rghyp(npos, model_f1)
    if(min(x_f1) >= 0){
      # hist(x_f1)
      break
    }
  }
  
  # Step 3 : append f0 and f1
  x_flou <- c(x_f0, x_f1)
  # hist(x_flou)
  
  return(x_flou)
}

# setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/simulated")
# for(conc in c(0.1, 0.4, 1, 2.5)){
#   for(r in 1:4){        # Rain settings
#     for(i in 1:5){      # Replicates
#       Targets <- lapply(1:3, function(x) simulateTarget(conc = conc, R = r))
#       NTCs <- lapply(1:3, function(x) simulateNTC())
#       
#       simulated <- data.frame(c(Targets, NTCs))
#       colnames(simulated) <- c(paste0("Target",1:3), paste0("NTC",1:3))
#       
#       if(conc==0.1){ label_conc <- "A" }else 
#       if(conc==0.4){ label_conc <- "B" }else 
#       if(conc==1  ){ label_conc <- "C" }else 
#       if(conc==2.5){ label_conc <- "D" }
#       
#       filename <- paste0("sim_conc",label_conc,"_R",r,"_rep",i,".csv")
#       write.csv(simulated, file = filename, row.names = FALSE)
#       # break
#     }
#     # break
#   }
#   # break
# }
