setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master")

df_orig <- readRDS("Dataset_t_sampled.RDS")

x <- df_orig[df_orig$react.ID == 357, -(1:11)]

hist(as.numeric(x), breaks = 60,
     main = "Plate 7 acp",
     sub = "React ID 357 (replicate 5)",
     xlab = "Fluorescence")

plot(as.numeric(x[!is.na(x)]),
     main = "React ID 357",
     sub = "Plate 7 acp",
     ylab = "Fluorescence",
     xlab = "Droplet")

df_orig[df_orig$plate.ID == 7 & df_orig$Target == "acp", 1:5]

