setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Trypsteen/amplitude_simulated4 - (.8) output_ddpcrquant/Final report CSV files")
source("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/4 - simulation/getEquallyLogSpacedIntervals.R")
library(magrittr)
library(dplyr)

concs <- getEquallyLogSpacedIntervals(0.1, 2.5, 5)
resFiles <- list.files(pattern = "*.csv")
list_conc <- list_rainSetting <- list_exp <- list_nPos <- list_nNeg <- list_nTot <- list_TP <- list_TN <- list_FP <- list_FN <- c()
for(i in resFiles){
  df <- read.csv(i, row.names = NULL)[-1, ]
  fnameParts <- strsplit(i, "_")[[1]] %>% gsub(".csv", "",., fixed = TRUE)
  conc <- gsub("conc", "", fnameParts[[3]]) 
  conc <- concs[which(conc == LETTERS)]
  rainSetting <- gsub("R", "", fnameParts[[4]]) 
  experiment <- gsub("rep", "", fnameParts[[5]])
  list_conc <- c(list_conc, rep(conc,3))
  list_rainSetting <- c(list_rainSetting, rep(rainSetting,3))
  list_exp <- c(list_exp, rep(experiment, 3))
  list_nPos <- c(list_nPos, pull(df, "positive.droplets"))
  list_nNeg <- c(list_nNeg, pull(df, "negative.droplets"))
  list_nTot <- c(list_nTot, pull(df, "total.droplets"))
  list_TP <- c(list_TP, pull(df, "TP"))
  list_TN <- c(list_TN, pull(df, "TN"))
  list_FP <- c(list_FP, pull(df, "FP"))
  list_FN <- c(list_FN, pull(df, "FN"))
  # break
}


df <- data.frame(conc = list_conc,
                 rain_setting = list_rainSetting,
                 exp = list_exp,
                 rep = rep(1:3, length(resFiles)),
                 n_tot = list_nTot,
                 n_neg = list_nNeg,
                 n_rain = 0,
                 n_pos = list_nPos,
                 TN = list_TN,
                 FP = list_FP,
                 TP = list_TP,
                 FN = list_FN,
                 estLambda = -log(list_nNeg/list_nTot))

write.csv(df, "Estimates_simulated_ddpcrquant.csv", row.names = FALSE)
