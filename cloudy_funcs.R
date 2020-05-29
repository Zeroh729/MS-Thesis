extractFlourescence <- function(df, row){
  data <- as.numeric(df[1,-c(1:11)])
  data <- data[!is.na(data)]
  return(data)
}

findpeaks <- function(vec, bw = 1, x.coo = c(1:length(vec))){  
  # bw : box width/window size, sets the sensitivity of the search
  
  list_max_x <- NULL;	
  list_max_y <- NULL;	
  list_min_x <- NULL;	
  list_min_y <- NULL;
  
  for(i in 1:(length(vec)-1)){
    
    # Step 1 : Get superior value
    if((i+1+bw)>length(vec)){
      leading_i <- length(vec)
    }else{
      leading_i <- i+1+bw
    }
    
    # Step 2 : Get inferior value
    if((i-bw) < 1){
      trailing_i <- 1
    }else{
      trailing_i <- i-bw
    }
    
    # Step 3 : Select window in two parts: values beyond superior, and values before inferior
    leading <- vec[(i+1):leading_i]
    trailing <- vec[trailing_i:(i-1)]
    
    # Step 4 : Get Max and Min values using two conditions
    cnt_trailingGtVal <- length(which(trailing > vec[i]))
    cnt_leadingGtVal <- length(which(leading > vec[i]))
    
    isValMax_trailing <- cnt_trailingGtVal == 0
    isValMax_leading  <- cnt_leadingGtVal == 0
    isValMin_trailing <- cnt_trailingGtVal == length(trailing)
    isValMin_leading  <- cnt_leadingGtVal == length(leading)
    
    if(isValMax_trailing & isValMax_leading){
      list_max_x <- c(list_max_x, x.coo[i])
      list_max_y <- c(list_max_y, vec[i])
    }
    
    if(isValMin_trailing & isValMin_leading){
      list_min_x <- c(list_min_x, x.coo[i])
      list_min_y <- c(list_min_y, vec[i])
    }
  }
  
  return(list("max_X" = list_max_x, "max_Y" = list_max_y, "min_X" = list_min_x, "min_Y" = list_min_y))
}