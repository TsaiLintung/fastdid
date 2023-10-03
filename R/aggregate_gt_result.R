aggregate_gt_result <- function(gt_result, cohort_sizes,
                                result_type){
  
  gt_att <- gt_result$att
  gt_inf_func <- gt_result$inf_func
  gt <- gt_result$gt
  
  group_time <- gt |> merge(cohort_sizes, by = "G")
  
  setorder(group_time, time, G) #change the order to match the order in gtatt
  
  if(result_type == "group_time"){
    
    targets <- group_time[, unique(G*max(time)+time)]
    inf_matrix <- as.matrix(gt_inf_func)
    agg_att <- as.vector(gt_att)
    
  } else {
    
    agg_sch <- get_aggregate_scheme(group_time, result_type)
    targets <- agg_sch$targets
    weights <- as.matrix(agg_sch$weights)
    colnames(weights) <- NULL
    
    inf_matrix <- as.matrix(gt_inf_func) %*% t(weights)
    agg_att <- weights %*% gt_att
    
  }
  return(list(inf_matrix = inf_matrix, agg_att = agg_att, targets = targets))
}

get_aggregate_scheme <- function(group_time, result_type){
  
  weights <- data.table()
  gt_count <- group_time[, .N]
  
  bool_to_pn <- function(x){ifelse(x, 1, -1)}
  
  if (result_type == "dynamic") {
    group_time[, target := time-G]
  } else if (result_type == "group") {
    group_time[, target := G*(bool_to_pn(time>=G))]
  } else if (result_type == "time") {
    group_time[, target := time*(bool_to_pn(time>=G))]
  } 
  
  min_target<- group_time[, min(target)]
  max_target<- group_time[, max(target)]
  
  targets <- group_time[, unique(target)]
  for(tar in targets){ #the order matters
    
    group_time[, agg_weight := 0]
    total_size <- group_time[target == tar, sum(cohort_size)]
    group_time[, weight := ifelse(target == tar, cohort_size/total_size, 0)]
    target_weights <- group_time[, .(weight)] |> transpose()
    
    weights <- rbind(weights, target_weights)
  }
  
  return(list(weights = weights, targets = targets))
}