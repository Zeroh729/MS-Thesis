setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/4 - simulation")
source("getEquallyLogSpacedIntervals.R")
library(ggplot2)
in_dir <- "simulated/"
files_firstRep <- list.files(in_dir, pattern = "_rep1.csv")
concs <- getEquallyLogSpacedIntervals(0.1, 2.5, 5)

for(f in files_firstRep){
  fname_parts <- strsplit(f, "_")[[1]]
  rain_setting <- substr(fname_parts[[3]], 2,2)
  conc <- substr(fname_parts[[2]], 5, 5)
  conc <- concs[which(conc == LETTERS)]
  fluo_rep1 <- read.csv(paste0(in_dir, f))[,1]
  d <- data.frame(fluo = fluo_rep1,
                     conc = factor(conc),
                     rain_setting = factor(rain_setting))
  if(f == files_firstRep[1]){
    df <- d
  }else{
    df <- rbind(df, d)
  }
}

label_rainSetting <- c("Low rain", "Moderate rain", "High rain", "V. High rain")
names(label_rainSetting) <- unique(df[['rain_setting']])
p <- ggplot(df, aes(x = fluo, fill=conc, col=conc)) +
  geom_histogram(bins=40, lwd=0.001, alpha=0.9) + 
  scale_x_continuous(n.breaks = 7) + 
  facet_wrap(conc ~ rain_setting, ncol = length(label_rainSetting), scales="free", labeller = labeller(rain_setting =label_rainSetting)) +
  labs(title = "All first sample replicates from Simulated data", y = "Frequency", x = "Droplet Fluoresence", color = "True Concentration", fill = "True Concentration") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 7))
p
ggsave(p, width=12.5, height = 7, filename = "Histogram of all first replicate.png")
