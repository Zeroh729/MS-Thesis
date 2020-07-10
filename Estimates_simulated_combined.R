df_est <- read.csv("D:/~Masters/~ MS-STAT/~THESIS/Code/Estimates_simulated_combined.csv")

library(ggplot2)
df_est$Attribute <- as.factor(df_est$Attribute) 
df_est$Attribute <- relevel(x = df_est$Attribute, ref = "umbrella")
df_est$conc2 <- as.factor(df_est$conc)


p <- ggplot(df_est, aes(x=Attribute, y=Value, colour=conc2)) +
  geom_boxplot() +
  labs(x = 'Method', y = 'DNA Concentration') +
  geom_hline(aes(yintercept=conc), color="black", linetype="dashed") +
  facet_grid(rows = vars(conc), cols = vars(rain_setting),scales = "free_y", ) +
  theme_light() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggsave(p, filename="Estimates_simulated_combined.png", device="png", width=10, height=7)



df_hists <- read.csv("D:/~Masters/~ MS-STAT/~THESIS/Code/test.csv")
df_hists$conc2 <- as.factor(df_hists$conc)
p2 <- ggplot(df_hists, aes(x=Target1, colour=conc2, fill=conc2)) +
  geom_histogram(alpha=0.6) +
  labs(x = 'Rain Setting', y = 'DNA Concentration') +
  facet_grid(rows = vars(conc), cols = vars(rain_setting),scales = "free_y", ) +
  theme_light() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggsave(p2, filename="Estimates_simulated_histograms.png", device="png", width=10, height=7)
