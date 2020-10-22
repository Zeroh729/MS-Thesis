setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/3 - simulation")
source("getEquallyLogSpacedIntervals.R")
filesPerMethod <- list(
  "cloudy"     = "Estimates_simulated_cloudy.csv",
  "umbrella"   = "Estimates_simulated_umbrella.csv",
  "ddpcRquant" = "Estimates_simulated_ddpcrquant.csv",
  # "EM_t"       = "Estimates_simulated_EM_t (ICL).csv",
  # "EM_tskew"   = "Estimates_simulated_EM_tskew (ICL).csv",
  "DropEmPCR_t" = "Estimates_simulated_DropEmPCR_t(1).csv",
  "DropEmPCR_tskew" = "Estimates_simulated_DropEmPCR_tskew(1).csv"
)
library(dplyr)
library(magrittr)
library(ggplot2)

for(i in 1:length(filesPerMethod)){
  df_temp <- read.csv(filesPerMethod[[i]]) %>% mutate(method = names(filesPerMethod)[[i]])
  if(i == 1){
    df <- df_temp
  }else{
    if(all(colnames(df) == colnames(df_temp))){
      df <- rbind(df, df_temp)      
    }else{
      print(paste0("Not the same column names! ", names(filesPerMethod)[[i]]))
    }
  }
}
df[['method']] <- factor(df[['method']], levels = unique(df[['method']]))

getLabelData <- function(df, option=1){
  # Option 1 : Summarize method across all concs
  # Option 2 : Summarize method for each conc
  if(option == 1){
    labelData <- df %>% dplyr::select(method) %>% distinct()
  }else if(option == 2){
    labelData <- df %>% dplyr::select(rain_setting, method) %>% distinct()
  }
  
  list_res <- list()
  for(m in unique(df[['method']])){
    print(paste0("method ", m))
    x <- subset(df, method==m)
    
    if(option == 1){
      res <- getResultStatistics(x)
      res$method <- m
      list_res[[length(list_res)+1]] <- res
    }else if(option == 2){
      for(r in unique(x[['rain_setting']])){
        print(paste0("rain_setting ", r))
        res <- getResultStatistics(subset(x, rain_setting==r))    
        res$method <- m
        res$rain_setting <- r
        list_res[[length(list_res)+1]] <- res    
      }
    }
  }
  labelData <- data.frame(t(matrix(unlist(list_res), ncol=length(list_res))))
  colnames(labelData) <- names(list_res[[1]])
  labelData[['method']] <- factor(labelData[['method']], levels=unique(df[['method']]))
  labelData[['fitRsquared']] <- as.numeric(labelData[['fitRsquared']])
  labelData[['fitPval']] <- as.numeric(labelData[['fitPval']])
  labelData[['fitResSE']] <- as.numeric(labelData[['fitResSE']])
  labelData[['label']] <- apply(labelData, 1, function(x) {
    paste0(
      "y ~ ",round(as.numeric(x[['beta0']]), 4), " + ",  round(as.numeric(x[['beta1']]), 4), "X\n",
      # "beta0 = ", round(as.numeric(x[['beta0']]), 4),"\n",
      # "beta1 = ", round(as.numeric(x[['beta1']]), 4),"\n",
      "R^2 = ", round(as.numeric(x[['fitRsquared']]), 4),"\n",
      # "Pval  = ", round(as.numeric(x[['fitPval']]), 4),"\n",
      "Res. SE = ", round(as.numeric(x[['fitResSE']]), 4))
  })
  return(labelData)
}


getResultStatistics <- function(x){
  trueConc <- pull(x, conc)
  estConc <- pull(x, estLambda)
  fit <- lm(log(trueConc) ~ log(estConc))
  fstat <- summary(fit)$fstatistic
  
  fitFormula <-   paste0("y ~ ", round(coefficients(fit)[1],2), " + ", 
           paste(sprintf("%.2f * %s", 
                         coefficients(fit)[-1],  
                         names(coefficients(fit)[-1])), 
                 collapse=" + ")
    )
  beta0 <- coefficients(fit)[[1]]
  beta1 <- coefficients(fit)[[2]]
  fitRsquared <- summary(fit)$r.squared
  fitPval <- pf(fstat[['value']], df1 = fstat[['numdf']], df2 = fstat[['dendf']], lower.tail = FALSE)
  fitResSE <- summary(fit)$sigma
  return(list(beta0 = beta0, beta1 = beta1, fitFormula=fitFormula, fitRsquared=fitRsquared, fitPval=fitPval, fitResSE=fitResSE))
}


addTheme_Common <- function(p, labelData){
  log_space <- (log(2.5) - log(0.1))/4
  y_min <- exp(log(0.1) - log_space)
  y_max <- exp(log(2.5) + log_space)
  y_breaks <- c(y_min, getEquallyLogSpacedIntervals(0.1, 2.5, 5), y_max) %>% round(2)
  
  p <- p + geom_abline(slope = 1, color="gray60") +
    geom_point() +
    geom_smooth(method="lm") + 
    geom_text(data = labelData, mapping = aes(x = y_min, y = y_max, vjust=1, hjust=0, label=label), color="black", size=2.7) +
    scale_y_continuous(breaks=y_breaks, trans = "log", labels = c("", y_breaks[2:6], "")) +
    scale_x_continuous(breaks=y_breaks, trans = "log", labels = c("", y_breaks[2:6], "")) +
    labs(x="Estimated Concentration (per partition)", y = "True Concentration (per partition)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size = 6))
}


plotData <- df %>% dplyr::select(conc, rain_setting, estLambda, method) %>% mutate(concFactor = factor(conc)) %>% mutate(Total = "All rain settings")

labelData1 <- getLabelData(df, option=1)
p1 <- ggplot(plotData, aes(x = estLambda, y = conc, color=method)) %>% 
  addTheme_Common(labelData1) +
  facet_grid(rows = vars(method), cols = vars(Total)) +
  theme(legend.position = "none")
p1
ggsave(p1, width = 3, height = 7, filename = "Results summary - log linear (across all rain).png")

labelData2 <- getLabelData(df, option=2)
rainLabels <- c("No rain", "Low rain", "Moderate rain", "High rain")
names(rainLabels) <- unique(df[['rain_setting']])
p2 <- ggplot(plotData, aes(x = estLambda, y = conc, color=method)) %>% 
  addTheme_Common(labelData2) +
  facet_grid(rows = vars(method), cols = vars(rain_setting), labeller = labeller(rain_setting = rainLabels))

p2
ggsave(p2, width = 9, height = 7, filename = "Results summary - log linear (per rain).png")


p <- cowplot::plot_grid(p1, NULL, p2, rel_widths = c(0.22, 0.03, 0.75), nrow = 1)

ggsave(p, width = 12, height = 7, filename = "Results summary - log linear (combined and per rain).png")





