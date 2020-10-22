setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/3 - simulation/simulated3")
source("../getEquallyLogSpacedIntervals.R")
library(tidyr)
library(magrittr)
library(dplyr)
library(ggplot2)

FILE_DIR <- "histograms/"
dir.create(FILE_DIR, recursive = TRUE, showWarnings = FALSE)

list_setting_ConcRain <- unique(substr(list.files(pattern = "*.csv"), start=5, stop = 12))
concs <- getEquallyLogSpacedIntervals(0.1, 2.5, 5)

for(setting in list_setting_ConcRain){
  files <- list.files(pattern = setting)
  conc <- substr(setting, 5, 5)
  conc <- concs[which(conc == LETTERS)]
  rain <- substr(setting, nchar(setting), nchar(setting))
    
  df <- NA
  for(file in files){
    rep <- substr(file, 17, 17)
    df_rep <- read.csv(file) %>% 
      pivot_longer(., colnames(.), names_to="sample") %>% 
      mutate(sample = factor(sample, levels = c("NTC1", "NTC2", "NTC3", "Target1", "Target2", "Target3"))) %>%
      mutate(rep = rep) %>% 
      arrange(sample)
    
    if(all(is.na(df))){
      df <- df_rep
    }else{
      df <- rbind(df_rep, df)
    }
  }
  
  # Graph NTCs and Targets in one image
  df <- mutate(df, type = as.factor(substr(sample, 1, ncwhguuuhar(as.character(df[['sample']]))-1)))
  p <- ggplot(df, aes(x=value, colour=type, fill=type)) +
    geom_histogram(bins=60, lwd=0.1, alpha=0.6) +
    labs(title=paste("Concentration:", conc), subtitle = paste("Rain setting:", rain), x = 'Samples', y = 'Rep') +
    facet_grid(rows = vars(rep), cols = vars(sample),scales = "free_y") +
    theme_light() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ggsave(p, filename=paste0(FILE_DIR, setting,".png"), device="png", width=14, height=7)
  
  
  # Graph Targets only in one image
  df2 <- filter(df, sample %in% c("Target1", "Target2", "Target3"))
  p2 <- ggplot(df2, aes(x=value)) +
    geom_histogram(bins=60, lwd=0.1, fill="#66D9DC", col="#00BFC4") +
    labs(title=paste("Concentration:", conc), subtitle = paste("Rain setting:", rain), x = 'Samples', y = 'Rep') +
    facet_grid(rows = vars(rep), cols = vars(sample),scales = "free_y") +
    theme_light() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ggsave(p2, filename=paste0(FILE_DIR, "Targets only - ", setting,".png"), device="png", width=14, height=7)
}
