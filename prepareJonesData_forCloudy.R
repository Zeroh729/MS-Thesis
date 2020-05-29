setwd("D:/~Masters/~ MS-STAT/~THESIS/Papers/(Supplementary Files) Jones/definetherain-master/data/Albumin/sampled")
filenames <- paste0(c("Alb 10e0","Alb 10e1","Alb 10e2","Alb 10e3","Alb 10e4", "Alb 10e5","Alb Neg"), ".csv")

for(i in 1:length(filenames)){
  print(paste0("Reading ", filenames[i]))
  Conc <- "0"
  if(i != length(filenames))
    Conc <- format(10^(i-1), scientific=FALSE)
  
  df_sample <- read.csv(filenames[i]) %>% mutate(Conc = Conc) %>%  
    mutate(Replicate = seq(nrow(.))) %>% 
    mutate(Target = "Albumin")
  
  writeLines(paste0("Sample number of cols : ", ncol(df_sample)))
  
  if(i == 1){
    df_jones <<- df_sample
  }else{
    if(ncol(df_sample) > ncol(df_jones)){
      writeLines("df_sample has more columns.")
      df_jones[,setdiff(colnames(df_sample), colnames(df_jones))] <- NA
    }else{
      writeLines("df_jones has more columns.")
      df_sample[,setdiff(colnames(df_jones), colnames(df_sample))] <- NA
    }
    df_jones <<- rbind(df_sample, df_jones)  
  }
  writeLines(paste0("Merged number of cols : ", ncol(df_sample)))
}

df_jones_final <- df_jones %>% arrange(desc(Conc)) %>%  select(Target, Conc, Replicate, everything())
write.csv(df_jones_final, "df_jones.csv", row.names = FALSE)

