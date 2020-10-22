setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/summarized")
source("D:/~Masters/~ MS-STAT/~THESIS/Code/utils.R")
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(magrittr)

out_dir <- "New Lievens Perf Evaluation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

list_methods <- list(
  list(method = "cloudy", filename = "Estimates_lievens_cloudy.csv"),
  list(method = "umbrella", filename = "Estimates_lievens_Umbrella.csv"),
  list(method = "em-t", filename = "Estimates_lievens_EM_t (BIC).csv"),
  list(method = "em-tskew", filename = "Estimates_lievens_EM_tskew (BIC).csv"))

for(m in list_methods){
  df <- read.csv(m[['filename']]) %>% 
    mutate(Enhancer = ifelse(is.na(Enhancer), "NA", Enhancer),
           Primers = as.character(Primers),
           Conc = as.character(Conc),
           Enhancer = as.character(Enhancer),
           Cycles = as.character(Cycles),
           Sonication = as.character(Sonication)) %>% 
    mutate(method = m[['method']]) %>% 
    filter(react.ID != 358 & react.ID != 383) %>% 
    select(method, react.ID,plate.ID,Target,Conc,Primers,Probe,Anneal,Touchdown,Enhancer,Cycles,Sonication,n_neg,n_rain,n_pos,n_total,estLambda,estLambda_l,estLambda_u) %>% 
    mutate(plate.ID2 = case_when(plate.ID != 5 ~ paste0(plate.ID),
                                 plate.ID == 5 ~ paste0(plate.ID, " (conc ", Conc, ")")))
    
  df2 <- df %>% mutate(subfactor_label = case_when(plate.ID %in% c(1, 7, 8, 9) ~ paste0(Target),
                                              plate.ID == 2 ~ paste(Target, "Primers", Primers),
                                              plate.ID == 3 ~ paste(Target, "Conc", Conc),
                                              plate.ID == 4 ~ paste(Target, "Enhancer", Enhancer),
                                              plate.ID == 5 ~ paste(Target, "Conc", Conc, "Cycles", Cycles),
                                              plate.ID == 6 ~ paste(Target, "Sonication", Sonication)
                                             ),
                       subfactors = case_when(plate.ID %in% c(1, 7, 8, 9) ~ "",
                                              plate.ID == 2 ~ "Primers",
                                              plate.ID == 3 ~ "Conc",
                                              plate.ID == 4 ~ "Enhancer",
                                              plate.ID == 5 ~ "Cycles",
                                              plate.ID == 6 ~ "Sonication"))
  if(list_methods[[1]][['filename']] == m[['filename']]){
    df_combined <- df2
  }else{
    df_combined <- rbind(df_combined, df2)
  }
  # break
}

plotData <- df_combined %>% group_by(method, plate.ID, subfactor_label) %>% 
  summarise(cv = sd(estLambda)/mean(estLambda)) %>% 
  mutate(cv_inv = 1-cv, 
         method = factor(method, levels = listOfList_toColumn(list_methods, "method"))) %>% 
  filter(plate.ID != 8)

p <- ggplot(plotData, aes(x=method, y=cv, color=method, fill = method)) +
  geom_violin(alpha=0.4, lwd=0.5) + 
  labs(title = "CV Distribution",
       subtitle = "For estimated concentrations of all experiments") +
  theme_bw()
ggsave(p, filename = "CV Distribution (overall).png", path = out_dir, width = 7, height = 4.5)

p <- p + facet_wrap(~ plate.ID, ncol = 3, scales = "free") +
  labs(subtitle = "For estimated concentrations per plate experiment") +
  theme_bw()
ggsave(p, filename = "CV Distribution (per plate).png", path = out_dir, width = 10, height = 7)

for(i in unique(df_combined$plate.ID2)){
  plotData <- filter(df_combined, plate.ID2 == i) %>% mutate(method = factor(method, levels = listOfList_toColumn(list_methods, "method")))
  max_n_reps <- plotData %>% group_by(method, subfactor_label) %>% summarise(n = n()) %>% pull() %>% max()
  if(max_n_reps <= 8) {
    pd <- position_dodge(0.75)
  }else if(max_n_reps >= 16){
    pd <- position_dodge(1)
  }else {
    pd <- position_dodge(0.5)
  }
  n_targets <- length(unique(plotData[['Target']]))
  if(n_targets == 2){
    facet_ncols <- 1
  }else if(max_n_reps >= 8){
    facet_ncols <- 2
  }else{
    facet_ncols <- 3
  }

  subfactor <- plotData[['subfactors']][[1]]
  if(subfactor == ""){
    p <- ggplot(plotData, aes_string(x = "method", y = "estLambda", color = "method", group="react.ID", label="react.ID"))
  }else{
    p <- ggplot(plotData, aes_string(x = "method", y = "estLambda", color = "method", group="react.ID", label="react.ID", shape = subfactor))
  }

  p <- p +geom_pointrange(aes( ymin = estLambda_l, ymax = estLambda_u), position=pd) +
    labs(title = paste("Plate", i), y = "Estimated DNA concentration (per partition)") +
    facet_wrap(~ Target, ncol=ifelse(n_targets == 2, 1, 3), scales = "free")+
    theme_light() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  # print(p)
  ggsave(p, filename = paste0("Plate ", i,".png"), path = out_dir, width = 11, height = 8)

  p_labeled <- p + geom_text_repel(direction = "x", segment.alpha = 0.3, size = 2.5, position = pd)
  ggsave(p_labeled, filename = paste0("Plate ", i," (w label).png"), path = out_dir, width = 11, height = 8)

  p <- ggplot(plotData, aes_string(x = "method", y = "estLambda", color = "method", label="react.ID")) + 
    geom_boxplot() +
    labs(title = paste("Plate", i), y = "Estimated DNA concentration (per partition)") +
    facet_wrap(~ Target, ncol=ifelse(n_targets == 2, 1, 3), scales = "free")+
    theme_light() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  ggsave(p, filename = paste0("Boxplot - Plate ", i,".png"), path = out_dir, width = 11, height = 8)
}

plotData <- df_combined %>% filter(plate.ID == 3) %>% mutate(logConc = log(as.numeric(Conc), base = 2) - (log(250, base = 2)-1)) %>% 
  filter(method != "em-tskew")

ggplot(plotData, aes(x = estLambda, y = logConc, color=method)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(rows = vars(method)) +
  ylim(c(min(plotData$logConc), max(plotData$logConc))) +
  theme_bw()
