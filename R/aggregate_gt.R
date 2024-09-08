# high level -------------------------------------------------------------------

aggregate_gt <- function(all_gt_result, aux, p){
  results <- lapply(all_gt_result, aggregate_gt_outcome, aux, p)
  return(list(
    est = rbindlist(lapply(results, function(x) {x$result})),
    inf_func = lapply(results, function(x) {x$inf_func}),
    agg_weight_matrix = lapply(results, function(x) {x$weight_matrix})
  ))
}

aggregate_gt_outcome <- function(gt_result, aux, p){
  
  #get aggregation scheme from g-t to target parameters
  agg_sch <- get_agg_sch(gt_result, aux, p)

  if(!p$event_specific | is.na(p$cohortvar2)){
    att <- gt_result$att
    inf_func <- gt_result$inf_func %*% t(agg_sch$agg_weights)
  } else {
    att <- agg_sch$es_weight %*% gt_result$att
    inf_func <- gt_result$inf_func %*% t(agg_sch$agg_weights %*% agg_sch$es_weight)
  }
  
  #get att
  agg_att <- agg_sch$agg_weights %*% att
  
  #get influence function matrix
  
  inf_weights <- get_weight_influence(att, agg_sch, aux, p)
  inf_matrix <- inf_func + inf_weights 

  #get se
  agg_se <- get_se(inf_matrix, aux, p)
  
  # post process
  result <- data.table(agg_sch$targets, agg_att, agg_se$se)
  names(result) <- c("target", "att", "se")
  result[,`:=`(outcome = gt_result$outname,
              att_ciub = att+se*agg_se$crit_val,
              att_cilb = att-se*agg_se$crit_val)]
  
  return(list(result = result,
              inf_func = inf_matrix,
              weight_matrix = agg_sch$agg_weights))
}

# scheme ------------------------------------------------------------------------

#scheme for aggregation
get_agg_sch <- function(gt_result, aux, p){

  #create group_time
  id_dt <- data.table(weight = aux$weights/sum(aux$weights), G = aux$dt_inv[, G])
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group_time <- gt_result$gt |> merge(pg_dt, by = "G", sort = FALSE)
  group_time[, mg := ming(G)]
  group_time[, G1 := g1(G)]
  group_time[, G2 := g2(G)]
  setorder(group_time, time, mg, G1, G2) #change the order to match the order in gtatt
  if(!all(names(gt_result$att) == group_time[, paste0(G, ".", time)])){stop("some bug makes gt misaligned, please report this to the maintainer. Thanks.")}
  
  #get the event-specific matrix, and available ggts
  if(p$event_specific & !is.na(p$cohortvar2)){
    es <- get_es_scheme(group_time, aux, p)
    group_time <- es$group_time #some gt may not have availble effect (ex: g1 == g2)
    es_weight <- as.matrix(es$es_weight)
  } else {
    es_weight <- NULL
  }

  #choose the target based on aggregation type
  tg <- get_agg_targets(group_time, p)
  group_time <- tg$group_time
  targets <- tg$targets
  
  
  #get aggregation weights
  agg_weights <- data.table()
  for(tar in targets){ #the order matters

    group_time[, weight := 0] #weight is 0 if not a target
    group_time[target == tar & used, weight := pg/sum(pg)] 
    target_weights <- group_time[, .(weight)] |> transpose()

    agg_weights <- rbind(agg_weights, target_weights)
  }
  group_time[, pg := NULL]
  
  agg_weights <- as.matrix(agg_weights)
  
  return(list(agg_weights = agg_weights, #a matrix of each target and gt's weight in it 
              targets = targets,
              group_time = group_time,
              es_weight = es_weight))
}

#get the target parameters
get_agg_targets <- function(group_time, p){
  group_time[, post := as.numeric(ifelse(time >= g1(G), 1, -1))]
  switch(p$result_type,
    dynamic = group_time[, target := time-g1(G)],
    group = group_time[, target := g1(G)*post], # group * treated
    time =  group_time[, target := time*post],
    simple = group_time[, target := post],
    group_time = group_time[, target := paste0(g1(G), ".", time)],
    group_group_time = group_time[, target := paste0(G, ".", time)] ,
    dynamic_stagger = group_time[, target := paste0(time-g1(G), ".", g1(G)-g2(G))]
  )
  
  #allow custom aggregation scheme, this overides other stuff
  if(!is.na(p$exper$aggregate_scheme)){
    group_time[, target := eval(str2lang(p$exper$aggregate_scheme))]
  }

  targets <- group_time[, unique(target)]
  
  #for balanced cohort composition in dynamic setting
  #a cohort us only used if it is seen for all dynamic time
  if(p$result_type == "dynamic" & !is.na(p$balanced_event_time)){
    
    cohorts <- group_time[, .(max_et = max(target), #event time is target if in dynamic
                              min_et = min(target)), by = "G"]
    cohorts[, used := max_et >= p$balanced_event_time] #the max
    if(!cohorts[, any(used)]){stop("balanced_comp_range outside avalible range")}
    group_time[, used := G %in% cohorts[used == TRUE, G]]
    
    targets <- targets[targets <= p$balanced_event_time & targets >= cohorts[used == TRUE, min(min_et)]]
    
  } else{group_time[, used := TRUE]}
  
  return(list(group_time = group_time, targets = targets))
  
}

# influence function ------------------------------------------------------------

get_weight_influence <- function(att, agg_sch, aux, p){

  group <- agg_sch$group_time
  
  id_dt <- data.table(weight = aux$weights/sum(aux$weights), G = aux$dt_inv[, G])
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group <- group |> merge(pg_dt, by = "G", sort = FALSE)
  
  group[, time := as.integer(time)]
  
  if(is.na(p$cohortvar2)){
    group[, G := as.integer(G)]
    setorder(group, time, G)
  } else {
    group[, mg := ming(G)]
    group[, G1 := g1(G)]
    group[, G2 := g2(G)]
    setorder(group, time, mg, G1, G2) #sort
  }

  if(!p$parallel){
    inf_weights <- sapply(asplit(agg_sch$agg_weights, 1), function (x){
      get_weight_influence_param(x, group, att, aux, p)
    })
  } else {
    inf_weights <- matrix(unlist(mclapply(asplit(agg_sch$agg_weights, 1), function (x){
      get_weight_influence_param(x, group, att, aux, p)
    })), ncol = length(agg_sch$targets))
  }
  
  
  return(inf_weights)
  
}

#influence from weight calculation
get_weight_influence_param <- function(agg_weights, group, gt_att, aux, p) {

  keepers <- which(agg_weights != 0)
  group <- group[keepers,]
  #moving this outside will create a g*t*id matrix, not really worth the memory
  keepers_matrix <- as.matrix(aux$weights*sapply(1:nrow(group), function(g){as.integer(aux$dt_inv[, G] == group[g,G]) - group[g,pg]}))
  
  # gt weight = pgi / sum(pgi)
  if1 <- keepers_matrix/sum(group[,pg]) #numerator
  if2 <- rowSums(keepers_matrix) %*% t(group[,pg])/(sum(group[,pg])^2) #denominator
  
  # return the influence function for the weights
  inf_weight <- (if1 - if2) %*% as.vector(gt_att[keepers])
  inf_weight[abs(inf_weight) < sqrt(.Machine$double.eps)*10] <- 0 #fill zero 

  return(inf_weight)
}

# se -------------------------------------------------------------------

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

