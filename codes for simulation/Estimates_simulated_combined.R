df_est <- read.csv("D:/~Masters/~ MS-STAT/~THESIS/Code/Estimates_simulated_combined.csv")

library(ggplot2)
df_est$Attribute <- as.factor(df_est$Attribute) 
df_est$Attribute <- relevel(x = df_est$Attribute, ref = "umbrella")
df_est$conc2 <- as.factor(df_est$conc)
df_est$rain_setting <- as.factor(df_est$rain_setting)
df_est$Value <- as.numeric(df_est$Value)


p <- ggplot(df_est, aes(x=Attribute, y=Value, colour=conc2)) +
  geom_boxplot() +
  labs(x = 'Method', y = 'DNA Concentration') +
  geom_hline(aes(yintercept=conc), color="black", linetype="dashed") +
  facet_grid(rows = vars(conc), cols = vars(rain_setting),scales = "free_y", ) +
  theme_light() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggsave(p, filename="Estimates_simulated_combined.png", device="png", width=10, height=12)