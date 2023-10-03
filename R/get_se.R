get_se <- function(inf_matrix, boot, biters, cluster) {
  
  if(boot){
    
    top_quant <- 0.75
    bot_quant <- 0.25
    if(!is.null(cluster)){
      cluster_n <- aggregate(cluster, by=list(cluster), length)[,2]
      inf_matrix <- fsum(inf_matrix, cluster) / cluster_n #the mean without 0 for each cluster of each setting
      inf_matrix[is.na(inf_matrix)|abs(inf_matrix) < sqrt(.Machine$double.eps)*10] <- 0 #some cancels out, but not exactly zero because floating point zz
    }
    
    n_cluster <- nrow(inf_matrix)
    
    boot_results <- sqrt(n_cluster)*BMisc::multiplier_bootstrap(inf_matrix, biters = biters) %>% as.data.table()
    
    boot_top <- boot_results[, lapply(.SD, function(x) quantile(x, top_quant, type=1, na.rm = TRUE)),]
    boot_bot <- boot_results[, lapply(.SD, function(x) quantile(x, bot_quant, type=1, na.rm = TRUE)),]
    
    dt_se <- rbind(boot_bot, boot_top) %>% transpose()
    names(dt_se) <- c("boot_bot", "boot_top")
    dt_se[, n_adjust := (n_cluster-1)/colSums(inf_matrix!=0)]
    se <- dt_se[,(boot_top-boot_bot)/(qnorm(top_quant) - qnorm(bot_quant))*n_adjust/sqrt(n_cluster)]
    
  } else {
    
    if(!is.null(cluster)){stop("clustering only available with bootstrap")}
    
    inf_matrix <- inf_matrix %>% as.data.table()
    se <- inf_matrix[, lapply(.SD, function(x) sd(x, na.rm = TRUE)*sqrt(length(x)-1)/length(x[x!=0]))] %>% as.vector() #should maybe use n-1 but did use n
    
  }
  return(unlist(se))
}

