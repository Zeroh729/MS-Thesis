setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")

if(!exists("df_orig")) df_orig <- read.csv("Dataset_t_sampled.csv")

library(ggplot2)
library(gridExtra)
library(dplyr)
library(magrittr)
source("Cloudy-V2-04_classification.R")
source("graph_utils.R")
source("D:/~Masters/~ MS-STAT/~THESIS/Code/classifier_EM.R")
source("D:/~Masters/~ MS-STAT/~THESIS/Code/classifier_EM_teigen.R")
source("D:/~Masters/~ MS-STAT/~THESIS/Code/utils.R")

prepareFlou <- function(flou){
  flou <- flou[!is.na(flou)]
  flou <- as.numeric(as.character((flou)))
  flou <- sort(flou)
  return(flou)
}

vecFlou_toDf  <- function(flou, classifier=""){
  flou <- prepareFlou(flou)  
  obs <- c(1:length(flou))
  cl <- rep("drp", length(flou))
  res_info <- list(classifier=classifier)
  mainClassifier <- strsplit(classifier, " ")[[1]][1]
  criteria <- gsub("(\\(|\\))", "", strsplit(classifier, " ")[[1]][2])
  if(mainClassifier == "cloudy"){
    cl <- cloudyClassifier(flou)
    cl <- factor(cl, levels = c("pos", "neg",  "rain"))
  }else if(mainClassifier == "cloudy_rain"){
    cl <- cloudyClassifier(flou, showRain = TRUE)
    cl <- factor(cl, levels = c("pos", "neg",  "rain"))
  }else if(mainClassifier == "EM"){
    criteria <- if(is.na(criteria)) "BIC" else criteria
    res <- emclassifier(flou, volDrp=0.85, crit = criteria)
    cl <- res$classification
    cl <- factor(cl, levels = c("pos", "neg",  "rain"))
    
    # TODO - Move all this in emclassifier function
    estParam <- data.frame(Mu = res$em$Mu, Sigma = sqrt(res$em$LTSigma), MixProp = res$em$pi) %>% 
      mutate(NegThres = rep(res$negThres, nrow(.))) %>% 
      mutate(PosThres = rep(res$posThres, nrow(.))) %>% 
      arrange(Mu)
    res_info$est_parameter <- estParam
    G <- res$G
    
    title <- bquote("Neg~N("~mu==.(round(estParam[1,"Mu"],2))~","~sigma==.(round(estParam[1,"Sigma"],2))~")")
    subtitle <- bquote("Pos~N("~mu==.(round(estParam[G,"Mu"],2))~","~sigma==.(round(estParam[G,"Sigma"],2))~")")
    if(G == 3){
      res_info$desc <- bquote("Rain~N("~mu==.(round(estParam[2,"Mu"],2))~","~sigma==.(round(estParam[2,"Sigma"],2))~")")
    }
    
    res_info$G <- G
    res_info$title <- title
    res_info$subtitle <- subtitle
  }else if(mainClassifier == "EM_t"){
    criteria <- if(is.na(criteria)) "BIC" else criteria
    res <- emclassifier_teigen(flou, volDrp=0.85, crit = criteria)
    cl <- res$classification
    cl <- factor(cl, levels = c("pos", "neg",  "rain"))
    
    # TODO - Move all this in emclassifier_teigen function
    estParam <- data.frame(Mu = res$em$parameters$mean, Sigma = sqrt(c(res$em$parameters$sigma)), Df = res$em$parameters$df, MixProp = res$em$parameters$pig) %>%
      mutate(NegThres = rep(res$negThres, nrow(.))) %>%
      mutate(PosThres = rep(res$posThres, nrow(.))) %>% 
      arrange(Mu)
    res_info$est_parameter <- estParam
    G <- res$G
    
    title <- bquote("Neg~"~t~"("~v==.(round(estParam[1,"Df"],2))~","~mu==.(round(estParam[1,"Mu"],2))~","~sigma==.(round(estParam[1,"Sigma"],2))~")")
    subtitle <- bquote("Pos~"~t~"("~v==.(round(estParam[G,"Df"],2))~","~mu==.(round(estParam[G,"Mu"],2))~","~sigma==.(round(estParam[G,"Sigma"],2))~")")
    if(G == 3){
      res_info$desc <- bquote("Rain~"~t~"("~v==.(round(estParam[2,"Df"],2))~","~mu==.(round(estParam[2,"Mu"],2))~","~sigma==.(round(estParam[2,"Sigma"],2))~")")
    }
    res_info$G <- G
    res_info$title <- title
    res_info$subtitle <- subtitle
  }
  res_info$classifier <- mainClassifier
  res_info$criteria <- if(is.na(criteria)) "" else paste0("(",criteria,")")
  data <- data.frame(flou, obs, cl)
  return(list(data, res_info))
} 

getList_Scatter <- function(listData){
  listPlots <- list()
  n_repl <- length(listData)
  for(i in 1:n_repl){
    data <- listData[[i]]
    plot <- ggplot(data, aes(x=id, y=flou, color=cl)) +
      geom_point() +
      scale_color_manual(values=c("#00AFBB", "#F8786F", "#999999")) +
      labs(title="", x="", y = "")
    if(i != 1) { plot <- plot + theme(legend.position = "none") }
    listPlots[[length(listPlots)+1]] <- plot
  }
  
  p1 <- listPlots[[1]]
  listPlots[[n_repl+1]] <- get_legend(p1)
  listPlots[[1]] <- p1 + theme(legend.position = "none")
  widths <- c(rep(95/n_repl, n_repl),5)
  nrow <- 1
  ncol <- n_repl+1
  return(list(listPlots=listPlots, widths=widths, ncol=ncol, nrow=nrow, xlab="flourescence", ylab="reaction ID"))
}

getList_Hist <- function(listData, listResInfo, wDens = FALSE){
  listPlots <- list()
  n_repl <- length(listData)
  for(i in 1:n_repl){
    if(!wDens){
      listPlots <- getHist(i, listPlots, listData)
    }else{
      listPlots <- getHistWDensity(i, listPlots, listData, listResInfo)
    }
  }
  
  nrow <- 2
  # ncol <- as.integer(n_repl/2)
  # widths <- c(rep(100/ncol, ncol))
  
  ncol <- as.integer(n_repl/2)
  p1 <- listPlots[[1]]
  listPlots[(ncol+2):(n_repl+1)] <- listPlots[(ncol+1):(n_repl)]
  listPlots[[(ncol+1)]] <- get_legend(p1)
  listPlots[[1]] <- p1 + theme(legend.position = "none")
  widths <- c(rep(95/ncol, ncol),5)
  ncol <- ncol + 1
  
  return(list(listPlots=listPlots, widths=widths, ncol=ncol, nrow=nrow, xlab="flourescence", ylab="frequency"))
}

getHist <- function(i, listPlots, listData){
  data <- listData[[i]] %>% mutate(cl = rep("drp", nrow(.)))
  plot <- ggplot(data, aes(x=flou, color=cl)) +
    geom_histogram(color="darkblue", fill="lightblue") +
    geom_density(alpha=.4, fill="lightblue" ) +
    scale_color_manual(values=c("#00AFBB", "#F8786F", "#999999")) +
    labs(title="", x="", y = "")
  if(i != 1) { plot <- plot + theme(legend.position = "none") }
  listPlots[[length(listPlots)+1]] <- plot
  return(listPlots)
}

getHistWDensity <- function(i, listPlots, listData, listResInfo){
  data <- listData[[i]]
  resInfo <- listResInfo[[i]]
  estParam <- resInfo$est_parameter
  fills <- c("pos" = "green", "neg"="red", "rain"="gray30")
  funcShaded <- function(x, clus) {
    if(resInfo$classifier == "EM"){
      y <- dnorm(x, mean = estParam[clus,"Mu"], sd = estParam[clus,"Sigma"])     
    }else if(resInfo$classifier == "EM_t"){
      y <- teigen_dist(x, df=estParam[clus,"Df"], mu=estParam[clus,"Mu"], sigma=estParam[clus,"Sigma"]^2)
    }
    y <- y * estParam[clus, "MixProp"]
    return(y)
  }
  
  G <- resInfo$G
  xlabel <- if("desc" %in% names(resInfo)) resInfo$desc else ""
  title <- if("title" %in% names(resInfo)) resInfo$title else ""
  subtitle <- if("subtitle" %in% names(resInfo)) resInfo$subtitle else ""
  negThres <- estParam[1, "NegThres"]
  posThres <- if(G == 3) estParam[1, "PosThres"] else NA
  
  plot <- ggplot(data, aes(x=flou, fill=cl)) +
    geom_histogram(aes(y =..density..),color="gray64", fill="gray88") +
    stat_function(data = . %>% filter(cl=="neg"), fun = funcShaded, geom="area", alpha=0.3, args=list(clus=1)) +
    stat_function(data = . %>% filter(cl=="pos"), fun = funcShaded, geom="area", alpha=0.3, args=list(clus=G)) +
    scale_fill_manual(values=fills) +
    ggtitle(label = title, subtitle = subtitle) +
    xlab(xlabel) +
    ylab("") +
    theme_minimal() +
    theme(axis.title.x = element_text(color = "gray60", size = 8, hjust = 0),
          plot.title = element_text(size = 8.5), plot.subtitle = element_text(size = 8.5))
  
  if(G == 3){
    plot <- plot +
      stat_function(data = . %>% filter(cl=="rain"), fun = funcShaded, geom="area", alpha=0.3, args=list(clus=2))
  }
  
  for(thres in c(negThres, posThres)){
    if(!is.na(thres)){
      # xdens <- summary(density(data$flou)$x)
      thresLabel_offset <- (max(data$flou) - min(data$flou))/25#(xdens[['Max.']] - xdens[['Min.']])/25
      plot <- plot +
        geom_vline(xintercept = thres) +
        geom_text(aes(y=summary(density(flou)$y)[['Max.']]),x=thres+thresLabel_offset, label=round(thres,2), colour="gray45", text=element_text(size=6))
    }
  }
  
  if(i != 1) { plot <- plot + theme(legend.position = "none") }
  listPlots[[length(listPlots)+1]] <- plot
  return(listPlots)
}

classifyFlou <- function(list_vecFlou, classifier){
  listData <- list()
  listResInfo <- list()
  for(i in 1:nrow(list_vecFlou)){
    # returns id | flou | cl=["", "cloudy", "cloudy_rain"]
    g(data, res_info) %=% vecFlou_toDf(list_vecFlou[i,], classifier)
    listData[[length(listData)+1]] <- data
    listResInfo[[length(listResInfo)+1]] <- res_info
  }
  return(list(listData, listResInfo))
}

graphFacetFlou <- function(list_vecFlou, title, classifier = "", plot="scatter", saveImg = F){ #scatter|hist
  g(listData, listResInfo) %=% classifyFlou(list_vecFlou, classifier)
  if(plot == "scatter"){
    p <- getList_Scatter(listData)
  }else if(plot == "hist"){
    p <- getList_Hist(listData, listResInfo)
  }else if(plot == "histWDens"){
    p <- getList_Hist(listData, listResInfo, wDens = TRUE)
  }
  listPlots <- p[['listPlots']]
  widths <- p[['widths']]
  ncol <- p[['ncol']]
  nrow <- p[['nrow']]
  xlab <- p[['ncol']]
  ylab <- p[['nrow']]
  
  facet <- do.call(grid.arrange, c(listPlots, nrow = nrow, ncol=ncol, widths=list(widths), top=title, name="name"))
  if(saveImg){
    filedir <- paste0("plot_",plot,"-",classifier)
    dir.create(filedir, showWarnings = FALSE)
    filename <- paste0(filedir,"/",title)
    height <- if(plot == "scatter") 4.5 else 8
    ggsave(paste0(filename,".png"), facet,"png", width = 12, height = height, units = c("in"))
  }
}
graphFacetFlou_Reps <- function(list_vecFlou, plate, target, classifier, plot, factorID="", saveImg){
  for(repID in 1:4){
    title <- custom_title(plate, target, factorID, repID, classifier) #plate, target, factorID, replicateID, classifier
    sampledFlou <- filterReplicates(nrow(list_vecFlou), plate, factorID, repID)
    graphFacetFlou(list_vecFlou[sampledFlou,], title, classifier=classifier, plot=plot, saveImg=saveImg) # [TC1507, M88017]
  }
}

graphFacetFlou_Facs <- function(list_vecFlou, plate, target, classifier, plot, saveImg){
  if(plate == 5){
    for(factID in 1:2){
      graphFacetFlou_Reps(list_vecFlou, plate, target, classifier=classifier, plot=plot, factorID=factID, saveImg=saveImg)
    }
  }else if(plate == 4){
    for(factID in 1:3){
        title <- custom_title(plate, target, factID, replicateID="", classifier) #plate, target, factorID, replicateID, classifier
        sampledFlou <- filterReplicates(nrow(list_vecFlou), plate, factID)
        graphFacetFlou(list_vecFlou[sampledFlou,], title, classifier=classifier, plot=plot, saveImg=saveImg) # [TC1507, M88017]    
      }
  }
}

filterReplicates <- function(n_repl, plate, factorID=1, replicateID=1){
  n_repl <- 1:n_repl
  if(plate == 5){
    startoffset <- c(1, 17)[factorID] # 1 - Conc 8000; 17 - Conc 4000
    n_repl <- n_repl[startoffset : (startoffset + 15)]
  }
  if(plate == 4){
    startoffset <- c(1, 5, 9)[factorID] # 1 - NA; 5 - DMSO2%; 9 - Trehalose0.2M
    n_repl <- n_repl[startoffset : (startoffset + 3)]
  }
  if(plate == 3 || plate == 5 || plate == 6){
    skips <- 4 
    n_repl <- n_repl[which(n_repl %% skips == 0)] - skips + replicateID #change to 1,2,3,4 to offset replicates
  }
  return(n_repl)
}

main <- function(classifier, plot="scatter", saveImg = F){
  for(i in c(7)){ #unique(df_orig$plate.ID)){
    print(paste("In Plate", i))
    
    x <- df_orig[df_orig$plate.ID == i,]
    for(j in unique(x$Target)){
      if(i == 7 && (j == "acp" || j == "M88017")){
        next
      }
      print(paste("   In Target", j))
      y <- x[x$Target == j,-c(1:11)]
      if(i %in% c(3, 6)){
        graphFacetFlou_Reps(y, i, j, classifier = classifier, plot=plot, saveImg = saveImg)
      }else if(i %in% c(4, 5)){
        graphFacetFlou_Facs(y, i, j, classifier = classifier, plot=plot, saveImg = saveImg)
      }else{
        title <- custom_title(i, j, "", "", classifier) #plate, target, factorID, replicateID, classifier
        graphFacetFlou(y, title, classifier = classifier, plot=plot, saveImg = saveImg) 
      }
    }
  }
  print("done!")
}

graphPlateTarget <- function(plate_id, target, classifier, factorID="", replicateID="", plot="scatter", saveImg=F){
  title <- custom_title(plate_id, target, "", "", classifier)
  graphFacetFlou(df_orig[df_orig$plate.ID == plate_id & df_orig$Target == target,-c(1:11)], title, plot = plot, classifier = classifier, saveImg=saveImg)
}

# classifier : "", "cloudy", "cloudy_rain", "EM (BIC | ICL)", "EM_t (BIC | ICL)"
# plot : "hist", "histWDens", "scatter"
# main(classifier="EM_t", plot="histWDens", saveImg = TRUE)

graphPlateTarget(1,"acp", "EM (ICL)", plot="histWDens", saveImg = FALSE)
# graphPlateTarget(4,"TC1507", "EM", plot="histWDens", saveImg = FALSE)
# graphPlateTarget(7,"M88017", "EM", plot="histWDens", saveImg = FALSE)

# Previous Codes
# ["le1", "cru", hmg", "acp", "M810", "TC1507", "M88017", "M1445", "GTS4032", "M88701","M89788","GT73"]
# graphPlateTarget(3, "M88017", classifier = "", saveImg = T) # [TC1507, M88017]
# graphPlateTarget(4, "TC1507", classifier = "", saveImg = T) # [TC1507, M88017]