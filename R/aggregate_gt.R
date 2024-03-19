aggregate_gt <- function(gt_result, aux, p){
  

  #release the stuff
  id_cohorts <- aux$dt_inv[, G]
  
  id_dt <- data.table(weight = aux$weights/sum(aux$weights), G = id_cohorts)
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group_time <- gt_result$gt |> merge(pg_dt, by = "G")
  
  setorder(group_time, time, G) #change the order to match the order in gtatt
  
  gt_result$inf_func <- as.matrix(gt_result$inf_func)
  
  if(p$result_type == "group_time"){
    
    #don't need to do anything
    targets <- group_time[, unique(G*max(time)+time)]
    inf_matrix <- gt_result$inf_func
    agg_att <- as.vector(gt_result$att)
    
  } else {
    
    #get which gt(s) is a part of the aggregated param
    agg_sch <- get_aggregate_scheme(group_time, p$result_type, aux$weights, id_cohorts, p$balanced_event_time)
    targets <- agg_sch$targets
    
    #aggregated att
    agg_att <- agg_sch$agg_weights %*% gt_result$att
    
    #get the influence from weight estimation
    #this needs to be optimized!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    inf_weights <- sapply(asplit(agg_sch$agg_weights, 1), function (x){
      get_weight_influence(x, gt_result$att, aux$weights, id_cohorts, group_time[, .(G, time)])
    })
    
    #aggregated influence function
    inf_matrix <- (gt_result$inf_func %*% t(agg_sch$agg_weights)) + inf_weights 
    
  }
  
  #get se from influence function
  agg_se <- get_se(inf_matrix, p$boot, p$biters, aux$cluster, p$clustervar)
  
  # post process
  result <- data.table(targets, agg_att, agg_se)
  names(result) <- c("target", "att", "se")
  result[,outcome := gt_result$outname]
  
  return(result)
}

get_aggregate_scheme <- function(group_time, result_type, weights, id_cohorts, balanced_event_time){
  
  #browser()
  
  agg_weights <- data.table()
  gt_count <- group_time[, .N]
  
  bool_to_pn <- function(x){ifelse(x, 1, -1)}
  
  #choose the target based on aggregation type
  if (result_type == "dynamic") {
    group_time[, target := time-G]
  } else if (result_type == "group") {
    group_time[, target := G*(bool_to_pn(time>=G))] # group * treated
  } else if (result_type == "time") {
    group_time[, target := time*(bool_to_pn(time>=G))] #calendar time * treated
  } else if (result_type == "simple") {
    group_time[, target := bool_to_pn(time>=G)] #treated / not treated
  }  
  
  targets <- group_time[, unique(target)]
  targets <- sort(as.integer(targets))
  
  #for balanced cohort composition in dynamic setting
  #a cohort us only used if it is seen for all dynamic time
  if(result_type == "dynamic" & !is.na(balanced_event_time)){

    cohorts <- group_time[, .(max_et = max(time-G),
                              min_et = min(time-G)), by = "G"]
    cohorts[, used := max_et >= balanced_event_time] #the max
    if(!cohorts[, any(used)]){stop("balanced_comp_range outside avalible range")}
    group_time[, used := G %in% cohorts[used == TRUE, G]]
    
    targets <- targets[targets <= balanced_event_time & targets >= cohorts[used == TRUE, min(min_et)]]
    
  } else{group_time[, used := TRUE]}
  
  for(tar in targets){ #the order matters
    
    group_time[, targeted := target == tar & used]
    
    total_pg <- group_time[targeted == TRUE, sum(pg)] #all gt that fits in the target
    group_time[, weight := ifelse(targeted, pg/total_pg, 0)] #weight is 0 if not a target
    target_weights <- group_time[, .(weight)] |> transpose()
    
    group_time[, targeted := NULL]
    
    agg_weights <- rbind(agg_weights, target_weights)
  }
  
  return(list(agg_weights = as.matrix(agg_weights), #a matrix of each target and gt's weight in it 
              targets = targets)) 
}

get_weight_influence <- function(agg_weights, gt_att, weights, id_cohorts, group) {

  keepers <- which(agg_weights > 0)
  
  id_dt <- data.table(weight = weights/sum(weights), G = id_cohorts)
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group <- group |> merge(pg_dt, by = "G")
  
  group[, time := as.integer(time)]
  group[, G := as.integer(G)]
  setorder(group, time, G)

  # effect of estimating weights in the numerator
  if1 <- sapply(keepers, function(k) {
    (weights*BMisc::TorF(id_cohorts == group[k,G]) - group[k,pg]) /
      sum(group[keepers,pg])
  })
  # effect of estimating weights in the denominator
  if2 <- base::rowSums(sapply(keepers, function(k) {
    weights*BMisc::TorF(id_cohorts == group[k,G]) - group[k,pg]
  })) %*%
    t(group[keepers,pg]/(sum(group[keepers,pg])^2))
  # return the influence function for the weights
  inf_weight <- (if1 - if2) %*% as.vector(gt_att[keepers])
  inf_weight[abs(inf_weight) < sqrt(.Machine$double.eps)*10] <- 0 #fill zero 
  return(inf_weight)
}

get_se <- function(inf_matrix, boot, biters, cluster, clustervar) {
  
  if(boot){
    
    top_quant <- 0.75
    bot_quant <- 0.25
    if(!allNA(clustervar)){
      #take average within the cluster
      cluster_n <- stats::aggregate(cluster, by=list(cluster), length)[,2]
      inf_matrix <- fsum(inf_matrix, cluster) / cluster_n #the mean without 0 for each cluster of each setting
    }
    
    boot_results <- BMisc::multiplier_bootstrap(inf_matrix, biters = biters) %>% as.data.table()
    
    boot_top <- boot_results[, lapply(.SD, function(x) stats::quantile(x, top_quant, type=1, na.rm = TRUE))]
    boot_bot <- boot_results[, lapply(.SD, function(x) stats::quantile(x, bot_quant, type=1, na.rm = TRUE))]
    
    dt_se <- rbind(boot_bot, boot_top) %>% transpose()
    names(dt_se) <- c("boot_bot", "boot_top")
    
    se <- dt_se[,(boot_top-boot_bot)/(qnorm(top_quant) - qnorm(bot_quant))]
    se[se < sqrt(.Machine$double.eps)*10] <- NA
    
  } else {

    inf_matrix <- inf_matrix  |> as.data.table()
    se <- inf_matrix[, lapply(.SD, function(x) sqrt(sum(x^2, na.rm = TRUE)/length(x)^2))] %>% as.vector() #should maybe use n-1 but did use n
    
  }
  return(unlist(se))
}

