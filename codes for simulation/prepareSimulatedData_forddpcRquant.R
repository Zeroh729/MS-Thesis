setwd("D:/~Masters/~ MS-STAT/~THESIS/Code/simulated2")
library(magrittr)

files <- list.files(pattern="*.csv")
list_well <- stringr::str_pad(1:(length(files)*6), 3, side = "left", pad = "0")
list_assay <- character(length(list_well))
list_sample <- character(length(list_well))

dir.create("amplitude", showWarnings = FALSE)

i <- 1
for(f in files){
  df <- read.csv(f)
  for(col in colnames(df)){
    # print(col)
    
    list_assay[[i]] <- gsub(".csv", "", f) %>% substr(5, nchar(.))
    list_sample[[i]] <- paste0(col,"_", list_assay[[i]])
    
    df_amp <- data.frame("Assay1 Amplitude" = df[[col]], check.names = FALSE)
    write.csv(df_amp, paste0("amplitude/Simulated_",list_well[[i]],"_Amplitude.csv"), row.names = FALSE)
    
    i <- i + 1
    # break
  }
  # break
}


df_header <- data.frame(Well = list_well,
                        Sample = list_sample,
                        TypeAssay="Ch1Unknown",
                        Assay=list_assay)
write.csv(df_header, "amplitude/Simulated.csv", row.names = FALSE)
