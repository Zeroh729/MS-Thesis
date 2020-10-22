setwd("D:/~Masters/~ MS-STAT/~THESIS/Code")
source("classifier_EM.R")
source("classifier_EM_tskew.R")
source("DropEmPCR.R")
source("utils.R")

setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")
source("mine_util.R")
source("Cloudy-V2-04_classification.R")
library(zoo)
library(dplyr)
library(beepr)
library(doParallel)
registerDoParallel(cores = 2)

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

# WARNING : will not work for EM or EM_t, replace their return functions to res_info first
calculate_lambda_em <- function(res_info){ 
  r <- list()
  classification <- res_info[['classification']]
  
  r$n_neg <- sum(classification == "neg") 
  r$n_rain <- sum(classification == "rain")
  r$n_pos <- sum(classification == "pos") 
  r$n_total <- length(classification)
  
  # Formula 2
  # if(res_info[['G']] == 3){
  #   # option 1 - All rain is counted as negative
  #   # r$n_neg <- r$n_neg + r$n_rain
  #   
  #   # option 2 - Check first if rain is closer to negative
  #   mus <- sort(c(res_info[['est_parameter']]$Mu))
  #   if(abs(mus[2] - mus[1]) < abs(mus[2] - mus[3])){
  #     writeLines(paste0("old neg : ", r$n_neg))
  #     r$n_neg <- r$n_neg + r$n_rain
  #     writeLines(paste0("new neg : ", r$n_neg))
  #   }
  # }
  
  r$estLambda <- -log(r$n_neg/r$n_tot)
  r$estLambda_l <- r$estLambda - 1.96 * sqrt((r$n_tot - r$n_neg)/(r$n_tot * r$n_neg))   #lower 95% confidence bound for lambda 1 (cpd)
  r$estLambda_u <- r$estLambda + 1.96 * sqrt((r$n_tot - r$n_neg)/(r$n_tot * r$n_neg))
  return(r)
}

mainLooped <- function(df, dataName, classifier="cloudy"){  # ("cloudy", "EM (BIC | ICF)", "EM_t (BIC | ICF)")
  res <- list()
  
  mainClassifier <- strsplit(classifier, " ")[[1]][1]
  criteria <- gsub("(\\(|\\))", "", strsplit(classifier, " ")[[1]][2])

  res <- foreach(i = 1:length(df[,1]), .export = c("calculate_lambda", "calculate_lambda_em", "calculate_estConc","getVolDrp")) %do% {
    source("D:/~Masters/~ MS-STAT/~THESIS/Code/classifier_EM.R")
    source("D:/~Masters/~ MS-STAT/~THESIS/Code/classifier_EM_tskew.R")
    source("D:/~Masters/~ MS-STAT/~THESIS/Code/DropEmPCR.R")
    source("D:/~Masters/~ MS-STAT/~THESIS/Code/utils.R")
    source("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/mine_util.R")
    source("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/Cloudy-V2-04_classification.R")
    
    print(paste0("In experiment ", i))
    drp <- as.numeric(df[i,])
    drp <- drp[!is.na(drp)]
    print(paste0("sum is ", sum(drp)))
    
    if(i %in% c(358, 383) || sum(drp) == 0){
       r <- list(n_neg = NA,
                 n_rain = NA,
                 n_pos = NA,
                 n_total = NA,
                 estLambda = NA,
                 estLambda_l = NA,
                 estLambda_u = NA)
      if(mainClassifier == "EM" || mainClassifier == "EM_t" || mainClassifier == "EM_tskew"){
        r[['n_components']] <- NA
        r[['critScore']] <- NA
      }
      
       if(mainClassifier == "EM_tskew"){
         r[['bic_G2']] <- NA
         r[['bic_G3']] <- NA
         r[['aic_G2']] <- NA
         r[['aic_G3']] <- NA
         r[['icl_G2']] <- NA
         r[['icl_G3']] <- NA
       }
       if(mainClassifier == "DropEmPCR_tskew" || mainClassifier == "DropEmPCR_t"){
         r[['G']] <- NA
       }
      i_res <- r
      return(i_res)
    }
    
    res_extra <- list()
    if(mainClassifier == "cloudy"){
      # cloudy is default
      classification <- cloudyClassifier(drp, silent = TRUE, showRain = TRUE)
      estLambda <- calculate_lambda(classification)
    }else if(grepl("DropEmPCR",mainClassifier)){
      if(mainClassifier == "DropEmPCR_t"){
        res_info <- DropEmPCR(drp, volDrp = getVolDrp(dataName), distr = "mvt")
        emres <- res_info$emres
        res_extra[['G']] <- emres$G
      }else if(mainClassifier == "DropEmPCR_tskew"){
        res_info <- DropEmPCR(drp, volDrp = getVolDrp(dataName), distr = "mst")
        emres <- res_info$emres
        res_extra[['G']] <- emres$G
      }
      classification <- emres[['classification']]
      estLambda <- calculate_lambda(classification)
    }else {
      # only for EM or EM_t
      if(mainClassifier == "EM"){
        emres <- emclassifier(drp, volDrp = getVolDrp(dataName), crit = criteria)
      }else if(mainClassifier == "EM_t"){
        res_info <- emclassifier_t(drp, volDrp = getVolDrp(dataName), crit = criteria)
      }else if(mainClassifier == "EM_tskew"){
        res_info <- emclassifier_tskew(drp, volDrp = getVolDrp(dataName), crit = criteria)
      }
      if(mainClassifier == "EM_t" || mainClassifier == "EM_tskew"){
        emres <- res_info$emres
        if(res_info$G == 2){
          temp_emres_G2 <- res_info$emres$em
          temp_emres_G3 <- res_info$emres$em_lower
        }else{
          temp_emres_G2 <- res_info$emres$em_lower
          temp_emres_G3 <- res_info$emres$em
        }
        res_extra[['bic_G2']] <- temp_emres_G2$bic
        res_extra[['bic_G3']] <- temp_emres_G3$bic
        res_extra[['icl_G2']] <- temp_emres_G2$ICL
        res_extra[['icl_G3']] <- temp_emres_G3$ICL
        res_extra[['aic_G2']] <- temp_emres_G2$aic
        res_extra[['aic_G3']] <- temp_emres_G3$aic
      }
      res_extra[['n_components']] <- emres$G
      res_extra[['critScore']] <- emres$critScore
      classification <- emres[['classification']]
      # estLambda <- calculate_lambda_em(res_info)
      estLambda <- calculate_lambda(classification)
    }
    
    i_res <- c(estLambda, res_extra)
    i_res
  }
  return(res)
}


getDataset <- function(dataName){
  factorCols <- NA
  df_orig <- NA
  if(dataName == "lievens"){
    setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")
    df_orig <- readRDS("Dataset_t_sampled.RDS")
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
    #  (formula2-opt1).csv
    write.csv(df_summarized, paste0("summarized/Estimates_",dataName,"_",classifier,".csv"), row.names = F)
    print("Done")    
  }
  beep()
  return(df_summarized)
}

dataName <- "lievens"     # "lievens" , "jones
.g(df_orig, factorCols) %=% getDataset(dataName)

# classifier <- "EM_tskew (ICL)"      # "EM (BIC | ICL)", "EM_t (BIC | ICL)", "EM_tskew (BIC | ICL | AIC)"  ,  "cloudy"
# classifier <- "DropEmPCR_t"
# run(df_orig, factorCols, dataName, classifier=classifier, save=TRUE)

classifier <- "DropEmPCR_tskew"
run(df_orig, factorCols, dataName, classifier=classifier, save=TRUE)

