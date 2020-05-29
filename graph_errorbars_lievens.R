library(ggplot2)
library(gridExtra)
library(dplyr)
library(magrittr)


calculate_SE <- function(data, title, classifier, factorname, saveImg){ #data must have the columns "factors" & "replicates" 
  est_summary <- unique(data[,c("plate.ID","Target","factors")])  %>% 
    mutate(factor = factorname) %>% 
    merge(., aggregate(est ~ factors, data, sd), by = "factors", sort = F) %>% 
    merge(., aggregate(est ~ factors, data, mean), by = "factors", sort = F) %>% 
    subset(., select=c(2,3,4,1,5,6))
  
  colnames(est_summary)[4:6] <- c("subfactor","se", "mean")
  
  if(saveImg){
    write.table(est_summary, outfilename, sep=",", append=T, row.names = F, col.names = !file.exists(outfilename))
  }else{
    print(est_summary)
  }
}


graph_ErrorBar <- function(data, title, classifier, factorname, saveImg){ #data must have the columns "factors" & "replicates" 
  pd <- position_dodge(0.1)
  p <- ggplot(data, aes(x=factors, y=est, colour=replicate, group=replicate)) + 
    geom_errorbar(aes(ymin=est_lower, ymax=est_upper), width=.1, position=pd) +
    labs(title=title, x=factorname, y = "Estimate") +
    geom_point(position=pd)
  if(saveImg){
    filedir <- paste0("plot_errorbars","-",classifier)
    dir.create(filedir, showWarnings = FALSE)
    filename <- paste0(filedir,"/",title)
    ggsave(paste0(filename,".png"), p,"png", width = 12, height = 4.5, units = c("in"))
  }else{
    print(p)
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

custom_title <- function(plate, target="", classifier=""){
  title <- paste("Plate",plate)
  
  title <- paste(title,"Target",target)
  
  if(classifier != "") {
    title <- paste0(title," (",classifier,")")
  }
  
  return(title)
}

main <- function(classifier, func, save = F){
  for(i in unique(df_est$plate.ID)){
    print(paste("In Plate", i))
    x <- df_est[df_est$plate.ID == i,]
    if(i %in% c(2, 3, 4, 6, 7, 8)){
      for(j in unique(x$Target)){
        print(paste("   In Target", j))
        y <- x[x$Target == j,]
        if(i == 2){
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
        func(y, custom_title(i, j, classifier), classifier, factorname, saveImg=save)
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
          func(z, custom_title(i, paste(j,"Conc",k), classifier), classifier, factorname, saveImg=save)
        }
      }
    }
  }
  print("done!")
}

setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/summarized")
df_est <- read.csv("Estimates_definetherain.csv")

# Classifiers : "cloudy", "cloudy_norain", "EM", "definetherain"
#
# OPTION 1 : Graph Error bars
# main(classifier="definetherain", graph_ErrorBar, save = T)
#
# OPTION 2 : Calculate Estimate Statistics
outfilename <- "EstimateStats_definetherain.csv"
main(classifier="definetherain", calculate_SE, save = T)


# graphPlateTarget(2, "acp", custom_title(2, "acp", "cloudy"))
