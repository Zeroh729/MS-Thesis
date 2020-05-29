setwd("D:/~Masters/~ MS-STAT/~THESIS/Code")
source("classifier_EM.R")

setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")
source("mine_util.R")
source("Cloudy-V2-04.R")
library(zoo)

getEstConc <- function(dataName, lambda){
  Vsamp <- 20
  
  if(dataName == "lievens"){
    Vdrp <- 0.85
  }else if(dataName == "jones"){
    Vdrp <- 0.91
  }
  
  return(lambda * Vsamp/Vdrp * 1000)
}

mainLooped <- function(df, classifier="cloudy"){  # ("cloudy", "EM")
  res <- list(est=c(), est_lower=c(), est_upper=c())
  for(i in 1:length(df[,1])){
    if(i %in%  c(358, 383)){
      res[['est']] <- c(res[['est']], NA)
      res[['est_lower']] <- c(res[['est_lower']], Inf)
      res[['est_upper']] <- c(res[['est_upper']], Inf)
      next
    }
    
    print(paste0("In experiment ", i))
    drp <- as.numeric(df[i,])
    drp <- drp[!is.na(drp)]
    if(classifier == "EM"){
      lambda_0 <- em(as.matrix(drp))
    
    }else{ 
      #cloudy is default
      lambda_0 <- cloudy(drp, dVol = 0.85, sVol = 20, plots = FALSE, silent = TRUE, vec = FALSE)[['lambda']]
    }
    
    res[['est']] <- c(res[['est']], lambda_0[['lambda']])
    res[['est_lower']] <- c(res[['est_lower']], lambda_0[['lower']])
    res[['est_upper']] <- c(res[['est_upper']], lambda_0[['upper']])
  }
  return(res)
}

dataName <- "jones"     # "lievens"
classifier <- "cloudy"      # "EM"  ,  "cloudy"


factorCols <- NA
df_orig <- NA
if(dataName == "lievens"){
  setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")
  df_orig <- read.csv("Dataset_t.csv")
  factorCols <- c(1:11)
}else if(dataName == "jones"){
  setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Jones/definetherain-master/data/Albumin/sampled")
  df_orig <- read.csv("df_jones.csv")
  factorCols <- c(1:3)
}
res <- mainLooped(df_orig[,-factorCols], classifier=classifier)

df_summarized <- df_orig[,factorCols] %>% 
  mutate(estLambda = res$est) %>% 
  mutate(estLambda_l = res$est_lower) %>% 
  mutate(estLambda_u = res$est_upper) %>% 
  mutate(estConc = getEstConc(dataName, estLambda)) %>% 
  mutate(estConc_l = getEstConc(dataName, estLambda_l)) %>% 
  mutate(estConc_u = getEstConc(dataName, estLambda_u))
dir.create("summarized", showWarnings=FALSE)
write.csv(df_summarized, paste0("summarized/Estimates_",dataName,"_",classifier,".csv"), row.names = F)
print("Done")
