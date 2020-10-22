library(dplyr)
library(magrittr)

setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/codes for simulation/4 - simulation")
filesPerMethod <- list(
  "cloudy"     = "Estimates_simulated_cloudy.csv",
  "umbrella"   = "Estimates_simulated_umbrella.csv",
  "ddpcRquant" = "Estimates_simulated_ddpcrquant (.8).csv",
  "EM-T" = "Estimates_simulated_DropEmPCR_t.csv",
  "EM-skewT" = "Estimates_simulated_DropEmPCR_tskew.csv"
)

out_dir <- "Results w DropEmPCR"
dir.create(out_dir, recursive = TRUE)
df_res <- NA
for(m in names(filesPerMethod)){
  file <- filesPerMethod[[m]]
  d <- read.csv(file) %>% dplyr::select(conc, rain_setting, exp, rep, n_tot, TP, TN, FP, FN) %>% 
    mutate(method = factor(m, levels = names(filesPerMethod)),
           concFactor = factor(conc),
           rainSettingFactor = factor(rain_setting))
  if(all(is.na(df_res))) {
    df_res <- d
  }else{
    df_res <- rbind(df_res, d)
  }
}
df_res <-  df_res %>% 
  mutate(accuracy = (TP + TN)/n_tot,
         sensitivity = TP/(TP + FN),
         specificity = TN/(TN + FP),
         FPR = FP/(TN + FP))

df_res %>% group_by(method) %>% 
  summarise(mean_acc  = mean(accuracy),    sd_acc  = sd(accuracy),
            mean_sens = mean(sensitivity), sd_sens = sd(sensitivity),
            mean_spec = mean(specificity), sd_spec = sd(specificity),
            mean_fpr  = mean(FPR),         sd_fpr  = sd(FPR)) %>% 
  write.csv(paste0(out_dir, "/Results summary - Accuracy.csv"))

# Sensitivity - FPR graph
p <- ggplot(df_res, aes(x = FPR, y=sensitivity, col=rainSettingFactor)) + 
  geom_point() +
  geom_line(aes(group=rainSettingFactor)) +
  facet_grid(cols=vars(method), rows=vars(rainSettingFactor)) +
  scale_color_manual(values = RColorBrewer::brewer.pal(5, "YlOrRd")[2:5], 
                     labels = c("No rain", "Low rain", "Moderate rain", "High rain")) +
  theme_bw() +
  labs(color = "Rain setting")
p
ggsave(p, width = 12, height = 3, filename = "Results summary - sensitivity FPR curve.png", path = out_dir)

# p <- ggplot(subset(df_res, method %in% c("cloudy", "EM_t", "EM_tskew")), 
#             aes(x = FPR, y=sensitivity, col=rainSettingFactor)) + 
#   geom_point() +
#   geom_line(aes(group=rainSettingFactor)) +
#   facet_grid(cols=vars(method)) +
#   scale_color_manual(values = RColorBrewer::brewer.pal(5, "YlOrRd")[2:5], 
#                      labels = c("No rain", "Low rain", "Moderate rain", "High rain")) +
#   theme_bw() +
#   labs(color = "Rain setting") + 
#   coord_cartesian(ylim = c(min(sensitivity), 1))
# p
ggsave(p, width = 12, height = 3, filename = "Results summary - sensitivity FPR curve (trimmed).png", path = out_dir)
