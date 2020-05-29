library(ggplot2)
library(gridExtra)
library(dplyr)
library(magrittr)

calculate_SE <- function(data, title, classifier, factorname, saveImg){ #data must have the columns "factors" & "replicates" 
  
  est_summary <- unique(data[,c("Target","factors")])  %>% 
    mutate(factor = factorname) %>% 
    merge(., aggregate(estConc ~ factors, data, mean), by = "factors", sort = F) %>% 
    merge(., aggregate(estConc ~ factors, data, sd), by = "factors", sort = F) %>% 
    merge(., aggregate(estLambda ~ factors, data, mean), by = "factors", sort = F) %>% 
    merge(., aggregate(estLambda ~ factors, data, sd), by = "factors", sort = F) %>% 
    subset(., select=c(2,3,1,4,5,6,7))
  
  colnames(est_summary)[3:7] <- c("subfactor", "conc_mean", "conc_se", "lambda_mean", "lambda_se")

  if(saveImg){
    write.table(est_summary, outfilename, sep=",", append=T, row.names = F, col.names = !file.exists(outfilename))
  }else{
    print(est_summary)
  }
}

graph_ErrorBar <- function(data, title, classifier, factorname, saveImg){ #data must have the columns "factors" & "replicates" 
    p <- ggplot(data, aes(x=factors, y=estLambda, colour=Replicate, group=Replicate)) + 
          geom_errorbar(aes(ymin=estLambda_l, ymax=estLambda_u), width=.1, position=position_dodge(0.1)) +
          labs(title=title, x=factorname, y = "Estimate") +
          geom_point(position=position_dodge(0.1))
  
  print(p)
  if(saveImg){
    filedir <- paste0("plot_errorbars -",classifier)
    dir.create(filedir, showWarnings = FALSE)
    filename <- paste0(filedir,"/",title)
    ggsave(paste0(filename,".png"), p,"png", width = 12, height = 4.5, units = c("in"))
  }
}


graph_LogLinearReg <- function(data, xName, yName, classifier, saveImg=F){
  plotData <- data %>% 
              mutate(plotData_X = log(as.numeric(data[,xName])+1, base = 10)) %>% 
              mutate(plotData_Y = log(as.numeric(data[,yName])+1, base = 10))
                     
  linearMod <- lm(plotData_Y ~ plotData_X, data=plotData)
  print(summary(linearMod))
  
  pLimit <- max(max(plotData$plotData_X),max(plotData$plotData_Y))
  
  p <- ggplot(plotData, aes(x = plotData_X, y = plotData_Y)) +
        geom_point() +
        geom_smooth(method=lm, se=TRUE) +
        geom_abline(aes(intercept=0, slope=1)) +
        ggtitle(label=paste0("Log Linear Regression (", classifier, ")"), subtitle = paste0("R-squared : ", round(summary(linearMod)$r.squared, 4))) +
        scale_x_continuous(name="log10(estConc + 1)", limits=c(0, pLimit)) +
        scale_y_continuous(name="log10(Conc + 1)", limits=c(0, pLimit))
  print(p)
  if(saveImg){
    filedir <- paste0("plot_loglinear","-",classifier)
    title <- custom_title("", "Albumin", classifier)
    dir.create(filedir, showWarnings = FALSE)
    filename <- paste0(filedir,"/",title)
    ggsave(paste0(filename,".png"), p,"png", width = 4, height = 4, units = c("in"))
  }
}

graphPlateTarget <- function(plate_id, target, title, saveImg=F){
  data <- df_est[df_est$plate.ID == plate_id & df_est$Target == target,]
  replicates <- 2
  factorname <- "Primers"
  data$factors <- as.factor(data$Primers)
  data$replicate <- as.factor(rep(c(1:replicates), length(unique(data$factors))))
  
  graph_ErrorBar(data, title, factorname, saveImg=saveImg)
}

custom_title <- function(plate="", target="", classifier=""){
  title <- ""
  if(plate != ""){
    title <- paste("Plate",plate)
  }
  
  if(target != ""){
    title <- paste(title,"Target",target)
  }
  
  if(classifier != "") {
    title <- paste0(title," (",classifier,")")
  }
  
  return(title)
}

main <- function(classifier, func, save = F){
  y <- df_est
  
  factorname <- "Conc"
  y$factors <- as.factor(y[[factorname]])
  
  func(y, custom_title("", "Albumin", classifier), classifier, factorname, saveImg=save)
  print("done!")
}

classifier <- "cloudy"           #  [ "cloudy", "definetherain" ]
dataFilename <- paste0("Estimates_jones_",classifier,".csv")
setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Jones/definetherain-master/data/Albumin/sampled")
df_est <- read.csv(dataFilename)

# Classifiers : "cloudy", "cloudy_norain", "EM", "definetherain"
#
# OPTION 1 : Graph Error bars
main(classifier=classifier, graph_ErrorBar, save = F)
#
# OPTION 2 : Calculate Estimate Statistics
outfilename <- paste0("EstimateStats_",classifier,".csv")
main(classifier=classifier, calculate_SE, save = T)
#
# OPTION 3 : Graph Log Linear Regression (for concentration series data)
graph_LogLinearReg(df_est, "estConc", "Conc", classifier, saveImg=T)

# graphPlateTarget(2, "acp", custom_title(2, "acp", "cloudy"))
