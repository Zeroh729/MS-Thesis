# install.packages("MixGHD")
# install.packages("gmp")
library(MixGHD)
library(dplyr)
library(magrittr)

setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")
if(!exists("df_orig")) df_orig <- read.csv("Dataset_t_sampled.csv")
flou <- df_orig %>% filter(react.ID==5) %>% select(-c(1:11))
flou <- flou[!is.na(flou)]
flou <- as.numeric(as.character((flou)))

r <- MGHD(flou, G=3)

hist(flou)
flou[1:10]
r@data$V1[1:10]
ztrans <- (flou - mean(flou))/sd(flou)
ztrans[1:10]

(r@gpar[[3]]$mu*sd(flou)+mean(flou))

MixGHD::plot(r)
