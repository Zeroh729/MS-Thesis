setwd("D:/~Masters/~ MS-STAT/~THESIS/Code")
source("utils.R")
source("common.R")
source("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/graph_utils.R")
setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/summarized")
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(magrittr)


list_methods <- list(
  list(method = "cloudy", filename = "Estimates_lievens_cloudy.csv"),
  list(method = "umbrella", filename = "Estimates_lievens_Umbrella.csv"),
  list(method = "em-t", filename = "Estimates_lievens_EM_t (BIC).csv"),
  list(method = "em-tskew", filename = "Estimates_lievens_EM_tskew (BIC).csv"))

main <- function(df_est, save = F){
  for(i in unique(df_est$plate.ID)){
    print(paste("In Plate", i))
    x <- df_est[df_est$plate.ID == i,]
    if(i %in% c(1, 2, 3, 4, 6, 7, 8, 9)){
      for(j in unique(x$Target)){
        print(paste("   In Target", j))
        y <- x[x$Target == j,]
        if(i == 1 || i == 9){
          replicates <- 4
          factorname <- "plate.ID"
        }else if(i == 2){
          replicates <- 2
          factorname <- "Primers"
        }else if(i %in% c(3, 6)){
          replicates <- 4
          if(i == 3) factorname <- "Conc"
          else if(i == 6) factorname <- "Sonication"
        }else if(i == 4){
          replicates <- 4
          factorname <- "Enhancer"
        }else if(i %in% c(7, 8)){
          replicates <- 1
          factorname <- "Anneal"
        }
        y$factors <- as.factor(y[[factorname]])
        y$replicate <- as.factor(rep(c(1:replicates), length(unique(y$factors))))
        
        graph_ErrorBar(y, custom_title(i, j), factorname, saveImg=save)
        # break
      }
    }else if(i == 5){
      for(j in unique(x$Target)){
        y <- x[x$Target == j,]
        for(k in unique(y$Conc)){
          z <- y[y$Conc == k,]
          replicates <- 4
          factorname <- "Cycles"
          z$factors <- as.factor(z[[factorname]])
          z$replicate <- as.factor(rep(c(1:replicates), length(unique(z$factors))))
          graph_ErrorBar(z, custom_title(i, paste(j,"Conc",k)), factorname, saveImg=save)
        }
      }
    }
    # break
  }
  print("done!")
}

graph_ErrorBar <- function(d, title, factorname, saveImg){ #data must have the columns "factors" & "replicates" 
  classEsts <- colSums(is.na(d[,estcols])) != nrow(d)
  classEsts <-  names(classEsts[classEsts])
  
  keepcols <- c("plate.ID", "Target", "react.ID", "factors", classEsts)
  
  plotData <- pivot_longer(d[,keepcols], cols = classEsts, names_to = "classifier", values_to = "estimates") %>% 
    mutate(classifier = sapply(strsplit(classifier,split="_"), function(x) x[[1]]))
  
  getBounds <- function (colnames_bounds, to_colname = "lower"){
    bounds <- pivot_longer(d[,c("react.ID", colnames_bounds)], cols = colnames_bounds, names_to = "classifier", values_to = to_colname) %>% 
      mutate(classifier = sapply(strsplit(classifier,split="_"), function(x) x[[1]]))
    return(bounds)
  }
  
  plotData <- inner_join(plotData, getBounds(estcols_l, to_colname="lower"), by = c("react.ID", "classifier")) %>% 
    inner_join(getBounds(estcols_u, to_colname="upper"), by = c("react.ID", "classifier")) %>% 
    mutate(classifier = factor(classifier, levels = sapply(list_methods, function(x) x$method)))
  
  if(ESTIMATE_TO_PLOT == "conc"){
    plotData[6:8] <- getDNAConc(plotData[6:8])
  }
  
  plotTitle <- title

  n_repl <- length(unique(plotData$react.ID))
  if(n_repl <= 8) {
    pd <- position_dodge(0.25)
  }else if(n_repl >= 20){
    pd <- position_dodge(0.6)
  }else {
    pd <- position_dodge(0.4)
  }
  
  p <- ggplot(plotData, aes(x=classifier, y=estimates, colour=factors, group=react.ID, label=react.ID)) + 
    geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
    labs(title=plotTitle, x="Classifiers", y = "Estimates", colour = factorname) +
    geom_point(position=pd) +
    theme_minimal()
  print(p)
  
  # Add hline of true conc if within the range of highest and lowest est. conc
  trueVal <- d[1,"Conc"][[1]]
  min_conc <- min(plotData[['lower']], na.rm = TRUE)
  max_conc <- max(plotData[['upper']], na.rm = TRUE)
  if(min_conc <= trueVal && trueVal <= max_conc){
    p <- p + geom_hline(yintercept = trueVal)
  }else{
    # Else, add note if under or overestimate
    note <- if(min_conc > trueVal) "overestimate" else "underestimate"
    p <- p + annotate("text", x=1, y=max_conc, label= note, colour = "gray")
  }
  
  
  p2 <- p + geom_text_repel(direction = "x", segment.alpha = 0.3, size = 2.5, position = pd)
  
  if(saveImg){
    filename <- title
    # directory for original plot
    filedir <- "plot_errorbars_combinedPlates"
    dir.create(filedir, showWarnings = FALSE)
    
    # directory for labeled plot
    filedir2 <- "plot_errorbars_combinedPlates (labeled)"
    dir.create(filedir2, showWarnings = FALSE)
    ggsave(paste0(filedir,"/",filename,".png"), p,"png", width = 9, height = 5, units = c("in"))
    ggsave(paste0(filedir2,"/",filename,".png"), p2,"png", width = 9, height = 5, units = c("in"))
  }
}

df <- getMergedEstimates(list_methods)
df <- addLambdaColumns(df, list_methods)


ESTIMATE_TO_PLOT <- "conc"  # "lambda" | "conc"

options(scipen = 7)
estcols_l <- colnames(df)[grepl("_lambda_l", colnames(df))]
estcols_u <- colnames(df)[grepl("_lambda_u", colnames(df))]
estcols <-  colnames(df)["_lambda" == substr(colnames(df), start = nchar(colnames(df))-6, stop = nchar(colnames(df)))]
  

main(df, save = TRUE)
