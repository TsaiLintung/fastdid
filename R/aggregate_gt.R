aggregate_gt <- function(gt_result, cohort_sizes, 
                         id_weights, id_cohorts,
                         result_type){
  
  gt_att <- gt_result$att
  gt_inf_func <- gt_result$inf_func
  gt <- gt_result$gt
  
  id_dt <- data.table(weight = id_weights/sum(id_weights), G = id_cohorts)
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group_time <- gt |> merge(pg_dt, by = "G")
  
  setorder(group_time, time, G) #change the order to match the order in gtatt
  
  gt_inf_func <- as.matrix(gt_inf_func)
  
  if(result_type == "group_time"){
    
    targets <- group_time[, unique(G*max(time)+time)]
    inf_matrix <- gt_inf_func
    agg_att <- as.vector(gt_att)
    
  } else {
    
    agg_sch <- get_aggregate_scheme(group_time, result_type, id_weights, id_cohorts)
    targets <- agg_sch$targets
    weights <- as.matrix(agg_sch$weights)
    
    #influence from estimating the weights
    inf_matrix <- gt_inf_func %*% t(weights) 
    
    #get the influence from weight estimation
    inf_weights <- sapply(asplit(weights, 1), function (x){
      get_weight_influence(x, gt_att, id_weights, id_cohorts, group_time[, .(G, time)])
    })
    inf_matrix <- inf_matrix + inf_weights 
    agg_att <- weights %*% gt_att
    
  }
  return(list(inf_matrix = inf_matrix, agg_att = agg_att, targets = targets))
}

get_aggregate_scheme <- function(group_time, result_type, id_weights, id_cohorts){
  
  weights <- data.table()
  gt_count <- group_time[, .N]
  
  bool_to_pn <- function(x){ifelse(x, 1, -1)}
  
  if (result_type == "dynamic") {
    group_time[, target := time-G]
  } else if (result_type == "group") {
    group_time[, target := G*(bool_to_pn(time>=G))]
  } else if (result_type == "time") {
    group_time[, target := time*(bool_to_pn(time>=G))]
  } else if (result_type == "simple") {
    group_time[, target := bool_to_pn(time>=G)]
  } 
  
  min_target<- group_time[, min(target)]
  max_target<- group_time[, max(target)]
  
  targets <- group_time[, unique(target)]
  targets <- sort(as.integer(targets))
  for(tar in targets){ #the order matters
    
    group_time[, agg_weight := 0]
    total_pg <- group_time[target == tar, sum(pg)]
    group_time[, weight := ifelse(target == tar, pg/total_pg, 0)]
    target_weights <- group_time[, .(weight)] |> transpose()
    
    weights <- rbind(weights, target_weights)
  }
  
  return(list(weights = weights, targets = targets))
}

get_weight_influence <- function(agg_weights, gt_att, id_weights, id_cohorts, group) {

  keepers <- which(agg_weights > 0)
  
  id_dt <- data.table(weight = id_weights/sum(id_weights), G = id_cohorts)
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group <- group |> merge(pg_dt, by = "G")
  
  group[, time := as.integer(time)]
  group[, G := as.integer(G)]
  setorder(group, time, G)

  # effect of estimating weights in the numerator
  if1 <- sapply(keepers, function(k) {
    (id_weights*BMisc::TorF(id_cohorts == group[k,G]) - group[k,pg]) /
      sum(group[keepers,pg])
  })
  # effect of estimating weights in the denominator
  if2 <- base::rowSums(sapply(keepers, function(k) {
    id_weights*BMisc::TorF(id_cohorts == group[k,G]) - group[k,pg]
  })) %*%
    t(group[keepers,pg]/(sum(group[keepers,pg])^2))
  # return the influence function for the weights
  inf_weight <- (if1 - if2) %*% as.vector(gt_att[keepers])
  inf_weight[abs(inf_weight) < sqrt(.Machine$double.eps)*10] <- 0 #fill zero 
  return(inf_weight)
}