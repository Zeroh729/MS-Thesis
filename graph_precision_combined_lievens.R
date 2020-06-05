setwd("D:/~Masters/~ MS-STAT/~THESIS")
df_cv <- read.csv("Dataset Batches - Tabular_CV.csv", stringsAsFactors = TRUE)

setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Lievens/ddPCR-master/summarized")

colnames(df_cv) <- c("Batch", "Target", "classifier", "CV")
df_cv$classifier <- as.factor(df_cv$classifier)
df_cv$Batch <- as.factor(df_cv$Batch)

ggplot(df_cv, aes(x = classifier, y = CV, colour = classifier)) +
  geom_violin(adjust = 1) +
  geom_jitter(size = 1, alpha = 0.5, width = 0.25, colour = 'black') +
  facet_grid( Batch ~ . ) +
  theme_bw()


p <- ggplot(df_cv, aes(x = classifier, y = CV, colour = classifier)) +
  geom_boxplot(notch = FALSE) +
  facet_grid( . ~ Batch ) +
  theme_bw()

ggsave("CV boxplots.png", p, width = 13, height = 7, device = "png")
