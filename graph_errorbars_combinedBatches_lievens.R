setwd("D:/~Masters/~ MS-STAT/~THESIS/Code")
source("utils.R")
source("common.R")
setwd("D:/~Masters/~ MS-STAT/~THESIS")
df_batch <- readxl::read_xlsx("Dataset Batches - Tabular.xlsx", sheet = "Precision_filtered (form2_op1)") %>% 
            rename(react.ID = `Reaction ID`) %>% 
            select(c("Batch", "Target", "react.ID", "Sample Type", "Conc"))

setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/summarized")

list_methods <- list(
  list(method = "cloudy", filename = "Estimates_lievens_cloudy.csv"),
  # list(method = "dtr", filename = "Estimates_definetherain.csv"), # Extract n_neg and n_total from RAW
  list(method = "em", filename = "Estimates_lievens_EM (ICL).csv"),
  list(method = "em-t", filename = "Estimates_lievens_EM_t (ICL).csv"),
  list(method = "em-tskew", filename = "Estimates_lievens_EM_tskew (BIC).csv")
)

df_ests <- getMergedEstimates(list_methods) %>% 
  addLambdaColumns(., list_methods) %>% 
  .[,-c(2:11)]

df <- inner_join(df_batch, df_ests, by="react.ID")

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(doParallel)
registerDoParallel(cores = 4)

save <- TRUE
ESTIMATE_TO_PLOT <- "conc"  # "lambda" | "conc"


getDNAConc <- function(lambda){
  return(lambda * 20 / 0.85 * 1000)
}
options(scipen = 7)
estcols_l <- colnames(df)[grepl("_lambda_l", colnames(df))]
estcols_u <- colnames(df)[grepl("_lambda_u", colnames(df))]
estcols <-  colnames(df)["_lambda" == substr(colnames(df), start = nchar(colnames(df))-6, stop = nchar(colnames(df)))]

for(batch in unique(df$Batch)){
  writeLines(paste0("In batch ", batch))
  foreach(target = unique(df[df$Batch == batch,]$Target), 
          .packages = c("tidyr", "dplyr", "ggplot2", "ggrepel")) %dopar% {
    source("D:/~Masters/~ MS-STAT/~THESIS/Code/utils.R")
    source("D:/~Masters/~ MS-STAT/~THESIS/Code/common.R")
    
    writeLines(paste0("     In Target ", target))
    d <- df[df$Batch == batch & df$Target == target,]
    classEsts <- colSums(is.na(d[,estcols])) != nrow(d)
    classEsts <-  names(classEsts[classEsts])
    
    keepcols <- c("Batch", "Target", "react.ID", "Sample Type", classEsts)
    
    plotData <- pivot_longer(d[,keepcols], cols = classEsts, names_to = "classifier", values_to = "estimates") %>% 
      mutate(classifier = sapply(strsplit(classifier,split="_"), function(x) x[[1]]))
    
    getBounds <- function (colnames_bounds, to_colname = "lower"){
      bounds <- pivot_longer(d[,c("react.ID", colnames_bounds)], cols = colnames_bounds, names_to = "classifier", values_to = to_colname) %>% 
        mutate(classifier = sapply(strsplit(classifier,split="_"), function(x) x[[1]]))
      return(bounds)
    }
    
    plotData <- inner_join(plotData, getBounds(estcols_l, to_colname="lower"), by = c("react.ID", "classifier")) %>% 
      inner_join(getBounds(estcols_u, to_colname="upper"), by = c("react.ID", "classifier"))
    
    if(ESTIMATE_TO_PLOT == "conc"){
      plotData[6:8] <- getDNAConc(plotData[6:8])
    }
    
    plotTitle <- paste("Target", target)
    plotSubtit <- paste("Test", batch)
    
    n_repl <- length(unique(plotData$react.ID))
    if(n_repl <= 8) {
      pd <- position_dodge(0.25)
    }else if(n_repl >= 20){
      pd <- position_dodge(0.6)
    }else {
      pd <- position_dodge(0.4)
    }

    p <- ggplot(plotData, aes(x=classifier, y=estimates, colour=`Sample Type`, group=react.ID, label=react.ID)) + 
      geom_errorbar(aes(ymin=lower, ymax=upper), width=.1, position=pd) +
      labs(title=plotTitle, subtitle = plotSubtit, x="Classifiers", y = "Estimates") +
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
    
    if(save){
      filename <- paste0("Test ",batch," - Target ", target)
      # directory for original plot
      filedir <- "plot_errorbars_combinedBatches"
      dir.create(filedir, showWarnings = FALSE)
      
      # directory for labeled plot
      filedir2 <- "plot_errorbars_combinedBatches (labeled)"
      dir.create(filedir2, showWarnings = FALSE)
      
      ggsave(paste0(filedir,"/",filename,".png"), p,"png", width = 9, height = 5, units = c("in"))
      ggsave(paste0(filedir2,"/",filename,".png"), p2,"png", width = 9, height = 5, units = c("in"))
    }
    # break
  }
  # break
}


# x[x$classifier == "cloudy",]
# x %>% merge(., aggregate(estimates ~ classifier, x, sd), by = x$classifier, sort = F)

