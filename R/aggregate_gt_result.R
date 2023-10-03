aggregate_gt_result <- function(gt_result, cohort_sizes,
                                result_type){
  
  gt_att <- gt_result$att
  gt_inf_func <- gt_result$inf_func
  
  group_time <- gt_att[,.(G, time)] |> merge(cohort_sizes, by = "G")
  
  if(result_type == "group_time"){
    
    targets <- group_time[, unique(G*max(time)+time)]
    inf_matrix <- as.matrix(gt_inf_func)
    agg_att <- as.vector(gt_att[, att])
    
  } else {
    
    agg_sch <- get_aggregate_scheme(group_time, result_type, balanced_composition)
    targets <- agg_sch$targets
    weights <- agg_sch$weights
    
    inf_matrix <- as.matrix(gt_inf_func) %*% t(as.matrix(weights))
    agg_att <- as.matrix(weights) %*% as.vector(gt_att[, att]) |> as.vector()
    
  }
  return(list(inf_matrix = inf_matrix, agg_att = agg_att, targets = targets))
}

get_aggregate_scheme <- function(group_time, result_type, balanced_composition){
  
  weights <- data.table()
  gt_count <- group_time[, .N]
  
  if(balanced_composition){stop("balanced_composition not implemented yet")}
  
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
  for(tar in targets){
    group_time[, agg_weight := 0]
    total_size <- group_time[target == tar, sum(cohort_size)]
    group_time[, weight := ifelse(target == tar, cohort_size/total_size, 0)]
    target_weights <- group_time[, .(weight)] |> transpose()
    names(target_weights) <- names(gt_inf_func)
    
    weights <- rbind(weights, target_weights)
  }
  return(list(weights = weights, targets = targets))
}