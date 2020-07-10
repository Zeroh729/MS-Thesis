cal_concentration <- function(nneg, total, volDrp){
  res <- list()
  volSamp <- 20
  
  # print(paste("Negatives",nneg))
  # print(paste("Total",total))
  npos <- total - nneg
  
  lambda <- -log(nneg/total)
  lower <- lambda - 1.96 * sqrt((total - length(nneg))/(total * length(nneg)))
  upper <- lambda + 1.96 * sqrt((total - length(nneg))/(total * length(nneg)))
  res$lambda <- c(lambda, upper, lower)
  names(res$lambda) <- c("lambda", "upper", "lower")
  
  conc <- lambda * volSamp/volDrp * 1000
  conc_lo <- lower * volSamp/volDrp * 1000
  conc_up <- upper * volSamp/volDrp * 1000
  res$conc <- c(conc, conc_up, conc_lo)
  names(res$conc) <- c("conc", "upper", "lower")
  return(res)
}


addLambdaColumns <- function(df, list_methods){
  for(method in list_methods){
    colname_n_neg <- paste0(method$method, "_n_neg")
    colname_n_total <- paste0(method$method, "_n_total")
    colname_lambda <- paste0(method$method, "_lambda")
    colname_lambda_l <- paste0(method$method, "_lambda_l")
    colname_lambda_u <- paste0(method$method, "_lambda_u")
    g(lambda,lambda_l,lambda_u) %=% getLambda(df[[colname_n_neg]], df[[colname_n_total]])
    df <- df %>% mutate(!!colname_lambda := lambda) %>% mutate(!!colname_lambda_l := lambda_l) %>% mutate(!!colname_lambda_u := lambda_u)
  }
  return(df)
}

getMergedEstimates <- function(list_methods){
  for(i in 1:length(list_methods)){
    method <- list_methods[[i]]
    print(method)
    df_method <- read.csv(method[['filename']])
    
    if(i == 1){
      df <- df_method[,c(1:11)]
    }
    df[[paste0(method$method,"_n_neg")]] <- df_method$n_neg
    df[[paste0(method$method,"_n_total")]] <- df_method$n_total
  }
  return(df)
}



getDNAConc <- function(lambda){
  return(lambda * 20 / 0.85 * 1000)
}

getLambda <- function(n_neg, n_total){
  lambda <- -log(n_neg / n_total)
  lambda_l <-  lambda - 1.96 * sqrt((n_total - n_neg)/(n_total * n_neg))
  lambda_u <-  lambda + 1.96 * sqrt((n_total - n_neg)/(n_total * n_neg))
  return(list(lambda, lambda_l, lambda_u))
}