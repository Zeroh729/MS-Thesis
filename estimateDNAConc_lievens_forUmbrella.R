setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/summarized")
source("D:/~Masters/~ MS-STAT/~THESIS/Code/utils.R")
library(dplyr)
library(doParallel)
library(magrittr)
df_orig <- readRDS("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/Dataset_t_sampled.RDS")
df <- mutate(df_orig, Enhancer = ifelse(is.na(Enhancer), "NA", Enhancer)) %>% 
  filter(react.ID != 358 & react.ID != 383)


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

estimateDNAConc <- function(df){
  source("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Jacobs/umbrella-master/1D/Umbrella_1d_V1.R")
  drp <- t(filter(df, react.ID != 358 & react.ID != 383)[,-(1:11)])
  drp <- data.frame(apply(drp, 2, as.numeric))
  
  # MAIN
  umbrellaRes <- Umbrella1d(drp)

  res_list <- list()
  for(t in colnames(drp)){
    classification <- ifelse(umbrellaRes[[t]]$droppi >= .8, "neg", "pos")
    classification <- na.omit(classification)
    res_list[[length(res_list)+1]] <- calculate_lambda(classification)
  }

  df_res <- data.frame(i = 1:length(res_list))
  for(i in attributes(res_list[[1]])$names){
    df_res[[i]] <- listOfList_toColumn(res_list, i)
  }
  print(df_res)
  
  return(df_res[,-1])
}

list_dfRes <- list()
for(i in unique(df[['plate.ID']])){
  df_plate <- df %>% filter(plate.ID == i)
  
  subFactor <- c("Target")
  if(i == 2){
    subFactor <- c(subFactor, "Primers")  
  }else if(i == 3){
    subFactor <- c(subFactor, "Conc")
  }else if(i == 4){
    subFactor <- c(subFactor,"Enhancer")
  }else if(i == 5){
    subFactor <- c(subFactor, "Conc", "Cycles")
  }else if(i == 6){
    subFactor <- c(subFactor, "Sonication")
  }
  
  list_groupComb <- list()
  for(f in subFactor){
     list_groupComb[[f]] <- unique(df_plate[[f]])
  }
  df_groupComb <- do.call(expand.grid, list_groupComb)
  
  list_resSample <- foreach(r = 1:nrow(df_groupComb)) %dopar% {
    df_sample <- df_plate
    for(f in colnames(df_groupComb)){
      df_sample <- df_sample %>% filter_at(vars(f), all_vars(. == df_groupComb[r, f]))
    }
    print(paste0("In Plate ", i, " Group ", paste0(df_groupComb[r,], collapse = " ")))
    res <- estimateDNAConc(df_sample)  # MAIN
    df_res <- cbind(df_sample[1:11], res)
  }
  for(resSample in list_resSample){
    list_dfRes[[length(list_dfRes)+1]] <- resSample
  }
}

df_res <- list_dfRes[[1]]
for(l in 2:length(list_dfRes)){
  df_res <- rbind(df_res, list_dfRes[[l]])
}
df_removed <- filter(df_orig, react.ID == 358 | react.ID == 383)[,1:11] %>% 
              mutate(n_neg=NA,
                     n_rain=NA,
                     n_pos=NA,
                     n_total=NA,
                     estLambda=NA,
                     estLambda_l=NA,
                     estLambda_u=NA)
df_res <- rbind(df_res, df_removed)

df_resFinal <- df_res %>% arrange(react.ID)

write.csv(df_resFinal, "Estimates_lievens_Umbrella.csv")