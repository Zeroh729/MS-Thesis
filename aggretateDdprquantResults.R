setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Trypsteen/amplitude_simulated2 - Final report CSV files")
library(magrittr)
library(dplyr)

list_conc <- list_rainSetting <- list_exp <- list_nPos <- list_nNeg <- list_nTot <- c()
resFiles <- list.files(pattern = "*.csv")
for(i in resFiles){
  df <- read.csv(i, row.names = NULL)[-1, ]
  fnameParts <- strsplit(i, "_")[[1]] %>% gsub(".csv", "",., fixed = TRUE)
  conc <- gsub("conc", "", fnameParts[[3]]) 
  conc <- case_when(conc == "A" ~ 0.1,
                    conc == "B" ~ 0.4, 
                    conc == "C" ~ 1,
                    conc == "D" ~ 2.5)
  rainSetting <- gsub("R", "", fnameParts[[4]]) 
  experiment <- gsub("rep", "", fnameParts[[5]])
  list_conc <- c(list_conc, rep(conc,3))
  list_rainSetting <- c(list_rainSetting, rep(rainSetting,3))
  list_exp <- c(list_exp, rep(experiment, 3))
  list_nPos <- c(list_nPos, pull(df, "positive.droplets"))
  list_nNeg <- c(list_nNeg, pull(df, "negative.droplets"))
  list_nTot <- c(list_nTot, pull(df, "total.droplets"))
  # break
}


df <- data.frame(conc = list_conc,
         rain_setting = list_rainSetting,
         exp = list_exp,
         rep = rep(1:3, length(resFiles)),
         n_tot = list_nTot,
         n_neg = list_nNeg,
         n_rain = 0,
         n_pos = list_nPos)

write.csv(df, "Estimates_simulated_ddpcrquant.csv", row.names = FALSE)
