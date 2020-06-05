setwd("D:/~Masters/~ MS-STAT/~THESIS")
df <- readxl::read_xlsx("Dataset Batches - Tabular.xlsx", sheet = "Precision_filtered")
setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/summarized")

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

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
  for(target in unique(df[df$Batch == batch,]$Target)){
    writeLines(paste0("     In Target ", target))
    d <- df[df$Batch == batch & df$Target == target,]

    classEsts <- colSums(is.na(d[,estcols])) != nrow(d)
    classEsts <-  names(classEsts[classEsts])
    
    keepcols <- c("Batch", "Target", "Reaction ID", "Sample Type", classEsts)
    
    plotData <- pivot_longer(d[,keepcols], cols = classEsts, names_to = "classifier", values_to = "estimates") %>% 
      mutate(classifier = sapply(strsplit(classifier,split="_"), function(x) x[[1]]))
    
    getBounds <- function (colnames_bounds, to_colname = "lower"){
      bounds <- pivot_longer(d[,c("Reaction ID", colnames_bounds)], cols = colnames_bounds, names_to = "classifier", values_to = to_colname) %>% 
        mutate(classifier = sapply(strsplit(classifier,split="_"), function(x) x[[1]]))
      return(bounds)
    }
    
    plotData <- inner_join(plotData, getBounds(estcols_l, to_colname="lower"), by = c("Reaction ID", "classifier")) %>% 
      inner_join(getBounds(estcols_u, to_colname="upper"), by = c("Reaction ID", "classifier"))
    
    if(ESTIMATE_TO_PLOT == "conc"){
      plotData[6:8] <- getDNAConc(plotData[6:8])
    }
    
    plotTitle <- paste("Target", target)
    plotSubtit <- paste("Test", batch)
    
    n_repl <- length(unique(plotData$`Reaction ID`))
    if(n_repl <= 8) {
      pd <- position_dodge(0.25)
    }else if(n_repl >= 20){
      pd <- position_dodge(0.6)
    }else {
      pd <- position_dodge(0.4)
    }

    p <- ggplot(plotData, aes(x=classifier, y=estimates, colour=`Sample Type`, group=`Reaction ID`, label=`Reaction ID`)) + 
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
      filedir <- "plot_errorbars_combined"
      dir.create(filedir, showWarnings = FALSE)
      
      # directory for labeled plot
      filedir2 <- "plot_errorbars_combined (labeled)"
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

