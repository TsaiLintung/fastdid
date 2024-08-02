aggregate_gt <- function(all_gt_result, aux, p){
  rbindlist(lapply(all_gt_result, aggregate_gt_outcome, aux, p))
}

aggregate_gt_outcome <- function(gt_result, aux, p){

  agg_sch <- get_agg_sch(gt_result, aux, p)
  
  #get att and se
  agg_att <- get_agg_att(gt_result, agg_sch, p)
  inf_matrix <- get_agg_inf(gt_result, agg_sch, aux, p)
  agg_se_result <- get_se(inf_matrix, aux, p)
  
  # post process
  result <- data.table(agg_sch$targets, agg_att, agg_se_result$se)
  names(result) <- c("target", "att", "se")
  result[,`:=`(outcome = gt_result$outname,
              att_ciub = att+se*agg_se_result$crit_val,
              att_cilb = att-se*agg_se_result$crit_val)]
  
  return(result)
}

#scheme for aggregation
get_agg_sch <- function(gt_result, aux, p){
  
  #setup stuff
  weights <- aux$weights
  id_cohorts <- aux$dt_inv[, G]
  result_type <- p$result_type
  agg_weights <- data.table()
  
  #create group_time
  id_dt <- data.table(weight = weights/sum(weights), G = id_cohorts)
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group_time <- gt_result$gt |> merge(pg_dt, by = "G")
  setorder(group_time, time, G) #change the order to match the order in gtatt
  gt_count <- group_time[, .N]
  
  #nothing to do
  if(result_type == "group_time"){
    return(list(targets = group_time[, paste0(G, ".", time)]))
  }
  
  #choose the target based on aggregation type
  group_time[, post := as.numeric(ifelse(time >= G, 1, -1))]
  if (result_type == "dynamic") {
    group_time[, target := time-G]
  } else if (result_type == "group") {
    group_time[, target := G*post] # group * treated
  } else if (result_type == "time") {
    group_time[, target := time*post] #calendar time * treated
  } else if (result_type == "simple") {
    group_time[, target := post] #treated / not treated
  }  
  
  targets <- sort(group_time[, unique(target)])
  
  #for balanced cohort composition in dynamic setting
  #a cohort us only used if it is seen for all dynamic time
  if(result_type == "dynamic" & !is.na(p$balanced_event_time)){
    
    cohorts <- group_time[, .(max_et = max(time-G),
                              min_et = min(time-G)), by = "G"]
    cohorts[, used := max_et >= p$balanced_event_time] #the max
    if(!cohorts[, any(used)]){stop("balanced_comp_range outside avalible range")}
    group_time[, used := G %in% cohorts[used == TRUE, G]]
    
    targets <- targets[targets <= p$balanced_event_time & targets >= cohorts[used == TRUE, min(min_et)]]
    
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
              targets = targets,
              group_time = group_time))
}

#aggregated influence function
get_agg_inf <- function(gt_result, agg_sch, aux, p){
  
  if(p$result_type == "group_time"){return(gt_result$inf_func)}
  
  inf_weights <- sapply(asplit(agg_sch$agg_weights, 1), function (x){
    get_weight_influence(x, gt_result$att, aux$weights, aux$dt_inv[, G], agg_sch$group_time[, .(G, time)])
  })
  
  #aggregated influence function
  inf_matrix <- (gt_result$inf_func %*% t(agg_sch$agg_weights)) + inf_weights 
  return(inf_matrix)
}

#influence from weight calculation
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

#aggregated att
get_agg_att <- function(gt_result, agg_sch, p){
  if(p$result_type == "group_time"){
    return(as.vector(gt_result$att))
  } else {
    return(agg_sch$agg_weights %*% gt_result$att)
  }
}

#aggregated standard error
get_se <- function(inf_matrix, aux, p) {

  if(p$boot){
    
    cluster <- aux$cluster
    
    top_quant <- 0.75
    bot_quant <- 0.25
    if(!allNA(p$clustervar)){
      #take average within the cluster
      cluster_n <- stats::aggregate(cluster, by=list(cluster), length)[,2]
      inf_matrix <- fsum(inf_matrix, cluster) / cluster_n #the mean without 0 for each cluster of each setting
    }
    
    boot_results <- BMisc::multiplier_bootstrap(inf_matrix, biters = p$biters) %>% as.data.table()
    
    boot_top <- boot_results[, lapply(.SD, function(x) stats::quantile(x, top_quant, type=1, na.rm = TRUE))]
    boot_bot <- boot_results[, lapply(.SD, function(x) stats::quantile(x, bot_quant, type=1, na.rm = TRUE))]
    
    dt_se <- rbind(boot_bot, boot_top) %>% transpose()
    names(dt_se) <- c("boot_bot", "boot_top")
    
    #get sigma
    se <- dt_se[,(boot_top-boot_bot)/(qnorm(top_quant) - qnorm(bot_quant))]
    se[se < sqrt(.Machine$double.eps)*10] <- NA


    
    
  } else {

    inf_matrix <- inf_matrix  |> as.data.table()
    se <- inf_matrix[, lapply(.SD, function(x) sqrt(sum(x^2, na.rm = TRUE)/length(x)^2))] %>% as.vector() #should maybe use n-1 but did use n
    
  }
  
  #get critical value
  crit_val <- NA
  point_crit_val <- qnorm(1-p$alpha/2)
  if(p$cband){
    boot_tv <- apply(boot_results, 1, function(b){max(abs(b/se), na.rm = TRUE)})
    boot_tv <= boot_tv[is.finite(boot_tv)]
    crit_val <- quantile(boot_tv, 1-p$alpha, type = 1, na.rm = TRUE) #alp set at 0.95 for now
  } 
  if(is.na(crit_val)|is.infinite(crit_val)|crit_val < point_crit_val){
    crit_val <- point_crit_val
  }
  
  return(list(se = unlist(se), crit_val = crit_val))
}

