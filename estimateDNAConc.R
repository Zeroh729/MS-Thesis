setwd("D:/~Masters/~ MS-STAT/~THESIS/Code")
source("classifier_EM.R")
source("classifier_EM_teigen.R")
source("utils.R")

setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")
source("mine_util.R")
source("Cloudy-V2-04_classification.R")
library(zoo)
library(dplyr)

listOfList_toColumn <- function(res_list, colName){
  return(sapply(res_list, function(x) return(x[[colName]])))
}
getVolDrp <- function(dataName){
  if(dataName == "lievens"){
    Vdrp <- 0.85
  }else if(dataName == "jones"){
    Vdrp <- 0.91
  }
  return(Vdrp)
}
calculate_estConc <- function(dataName, lambda){
  Vsamp <- 20
  Vdrp <- getVolDrp(dataName)
  return(lambda * Vsamp/Vdrp * 1000)
}

calculate_lambda <- function(classification){ # string vector or "pos", "neg", or "rain"
  r <- list()
  r$n_neg <- sum(classification == "neg") 
  r$n_rain <- sum(classification == "rain")
  r$n_pos <- sum(classification == "pos") 
  r$n_total <- length(classification)
  r$estLambda <- -log(r$n_neg/r$n_tot)
  r$estLambda_l <- r$estLambda - 1.96 * sqrt((r$n_tot - r$n_neg)/(r$n_tot * r$n_neg))   #lower 95% confidence bound for lambda 1 (cpd)
  r$estLambda_u <- r$estLambda + 1.96 * sqrt((r$n_tot - r$n_neg)/(r$n_tot * r$n_neg))
  return(r)
}
mainLooped <- function(df, dataName, classifier="cloudy"){  # ("cloudy", "EM (BIC | ICF)", "EM_t (BIC | ICF)")
  res <- list()
  
  mainClassifier <- strsplit(classifier, " ")[[1]][1]
  criteria <- gsub("(\\(|\\))", "", strsplit(classifier, " ")[[1]][2])
  
  for(i in 1:length(df[,1])){
    if(i %in% c(358, 383)){
       r <- list(n_neg = NA,
                 n_rain = NA,
                 n_pos = NA,
                 n_total = NA,
                 estLambda = NA,
                 estLambda_l = NA,
                 estLambda_u = NA)
      if(mainClassifier == "EM"){
        r[['n_components']] <- NA
        r[['critScore']] <- NA
      }else if(mainClassifier == "EM_t"){
        r[['modelname']] <- NA
        r[['n_components']] <- NA
        r[['critScore']] <- NA
      }
      res[[length(res)+1]] <- r
      next
    }
    
    print(paste0("In experiment ", i))
    drp <- as.numeric(df[i,])
    drp <- drp[!is.na(drp)]
    res_extra <- list()
    if(mainClassifier == "cloudy"){
      # cloudy is default
      classification <- cloudyClassifier(drp, silent = TRUE, showRain = TRUE)
    }else {
      # only for EM or EM_t
      if(mainClassifier == "EM"){
        emres <- emclassifier(drp, volDrp = getVolDrp(dataName), crit = criteria)
      }else if(mainClassifier == "EM_t"){
        emres <- emclassifier_teigen(drp, volDrp = getVolDrp(dataName), crit = criteria)
        res_extra[['modelname']] <- emres$em$modelname
      }
      res_extra[['n_components']] <- emres$G
      res_extra[['critScore']] <- emres$critScore
      classification <- emres[['classification']]
    }
    
    res[[length(res)+1]] <- c(calculate_lambda(classification), res_extra)
  }
  return(res)
}


getDataset <- function(dataName){
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
  return(list(df_orig, factorCols))
}

run <- function(df_orig, factorCols, dataName, classifier, save=FALSE){
  res <- mainLooped(df_orig[,-factorCols], dataName, classifier=classifier)
  
  # 1 - Create table of results
  df_summarized <- df_orig[,factorCols]
  for(i in attributes(res[[1]])$names){
    df_summarized[[i]] <- listOfList_toColumn(res, i)
  }
  df_summarized <- df_summarized %>% 
    mutate(estConc = calculate_estConc(dataName, estLambda)) %>% 
    mutate(estConc_l = calculate_estConc(dataName, estLambda_l)) %>% 
    mutate(estConc_u = calculate_estConc(dataName, estLambda_u))
  
  # 2 - Save table to CSV
  if(save){
    dir.create("summarized", showWarnings=FALSE)
    write.csv(df_summarized, paste0("summarized/Estimates_",dataName,"_",classifier,".csv"), row.names = F)
    print("Done")    
  }
  return(df_summarized)
}

dataName <- "lievens"     # "lievens" , "jones
g(df_orig, factorCols) %=% getDataset(dataName)

classifier <- "cloudy"      # "EM (BIC | ICL)", "EM_t (BIC | ICL)"  ,  "cloudy"
run(df_orig, factorCols, dataName, classifier=classifier, save=TRUE)





