setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")
library(ggplot2)
library(gridExtra)
source("Cloudy-V2-04_classification.R")
source("graph_utils.R")
setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Jones/definetherain-master/data/Albumin/sampled")


vecFlou_toDf  <- function(flou, classifier=""){
  flou <- flou[!is.na(flou)]
  flou <- as.numeric(as.character((flou)))
  id <- c(1:length(flou))
  cl <- rep("drp", length(flou))
  if(classifier == "cloudy"){
    cl <- cloudyClassifier(flou)
    cl <- factor(cl, levels = c("pos", "neg",  "rain"))
  }else if(classifier == "cloudy_rain"){
    cl <- cloudyClassifier(flou, showRain = TRUE)
    cl <- factor(cl, levels = c("pos", "neg",  "rain"))
  }
  return(data.frame(flou, id, cl))
}

# graphFacetFlou <- function(list_vecFlou, title, classifier = "", saveImg = F){
#   listPlots <- list()
#   n_repl <- nrow(list_vecFlou)
#   p1 <- ""
#   for(i in 1:n_repl){
#     data <- vecFlou_toDf(list_vecFlou[i,], classifier)
#     plot <- ggplot(data, aes(x=id, y=flou, color=cl)) +
#       geom_point() +
#       scale_color_manual(values=c("#00AFBB", "#F8786F", "#999999")) +
#       labs(title="", x="", y = "")
#     if(i == 1) { p1 <- plot }
#     else { plot <- plot + theme(legend.position = "none") }
#     listPlots[[length(listPlots)+1]] <- plot
#   }
#   listPlots[[n_repl+1]] <- get_legend(p1)
#   listPlots[[1]] <- p1 + theme(legend.position = "none")
#   widths <- c(rep(95/n_repl, n_repl),5)
#   facet <- do.call(grid.arrange, c(listPlots, ncol=(n_repl+1), widths=list(widths), top=title, left="flourescence", bottom="reaction ID"))
#   if(saveImg){
#     filename <- paste0("plots -",classifier,"/",title)
#     ggsave(paste0(filename,".png"), facet,"png", width = 12, height = 4.5, units = c("in"))
#   }
# }


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

getList_Hist <- function(listData){
  listPlots <- list()
  n_repl <- length(listData)
  for(i in 1:n_repl){
    data <- listData[[i]]
    plot <- ggplot(data, aes(x=flou, color=cl)) +
      geom_histogram(aes(y=..density..), color="darkblue", fill="lightblue") +
      geom_density(alpha=.4, fill="lightblue" ) +
      scale_color_manual(values=c("#00AFBB", "#F8786F", "#999999")) +
      labs(title="", x="", y = "")
    if(i != 1) { plot <- plot + theme(legend.position = "none") }
    listPlots[[length(listPlots)+1]] <- plot
  }
  
  nrow <- 2
  ncol <- as.integer(n_repl/2)
  
  p1 <- listPlots[[1]]
  listPlots[(ncol+2):(n_repl+1)] <- listPlots[(ncol+1):(n_repl)]
  listPlots[[(ncol+1)]] <- get_legend(p1)
  listPlots[[1]] <- p1 + theme(legend.position = "none")
  widths <- c(rep(95/ncol, ncol),5)
  return(list(listPlots=listPlots, widths=widths, ncol=(ncol+1), nrow=nrow, xlab="flourescence", ylab="frequency"))
}

classifyFlou <- function(list_vecFlou, classifier){
  listData <- list()
  for(i in 1:nrow(list_vecFlou)){
    # returns id | flou | cl=["", "cloudy", "cloudy_rain"]
    data <- vecFlou_toDf(list_vecFlou[i,], classifier) 
    listData[[length(listData)+1]] <- data
  }
  return(listData)
}

graphFacetFlou <- function(list_vecFlou, title, classifier = "", plot="scatter", saveImg = F){ #scatter|hist
  listData <- classifyFlou(list_vecFlou, classifier)
  if(plot == "scatter"){
    p <- getList_Scatter(listData)
  }else{ #if(plot == "hist")
    p <- getList_Hist(listData)
  }
  listPlots <- p[['listPlots']]
  widths <- p[['widths']]
  ncol <- p[['ncol']]
  nrow <- p[['nrow']]
  xlab <- p[['ncol']]
  ylab <- p[['nrow']]
  
  facet <- do.call(grid.arrange, c(listPlots, nrow = nrow, ncol=ncol, widths=list(widths), top=title, left=xlab, bottom=ylab))
  if(saveImg){
    filedir <- paste0("plot_",plot,"-",classifier)
    dir.create(filedir, showWarnings = FALSE)
    filename <- paste0(filedir,"/",title)
    ggsave(paste0(filename,".png"), facet,"png", width = 12, height = 4.5, units = c("in"))
  }
}

targets <- list.files()
for (t in targets){
  df <- read.csv(t)
  title <- strsplit(t, ".csv")[[1]]
  graphFacetFlou(df, title, classifier="", plot ="hist", saveImg = T) # "", "cloudy", "cloudy_rain"
}
