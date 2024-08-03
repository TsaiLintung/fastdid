# high level -------------------------------------------------------------------

aggregate_gt <- function(all_gt_result, aux, p){
  rbindlist(lapply(all_gt_result, aggregate_gt_outcome, aux, p))
}

aggregate_gt_outcome <- function(gt_result, aux, p){
  
  agg_sch <- get_agg_sch(gt_result, aux, p)
  
  #get att
  agg_att <- agg_sch$agg_weights %*% gt_result$att
  #get influence function
  inf_weights <- sapply(asplit(agg_sch$agg_weights, 1), function (x){
    get_weight_influence(x, gt_result, aux, agg_sch, p)
  })
  inf_matrix <- (gt_result$inf_func %*% t(agg_sch$agg_weights)) + inf_weights 
  
  #get se
  agg_se <- get_se(inf_matrix, aux, p)
  
  # post process
  result <- data.table(agg_sch$targets, agg_att, agg_se$se)
  names(result) <- c("target", "att", "se")
  result[,`:=`(outcome = gt_result$outname,
              att_ciub = att+se*agg_se$crit_val,
              att_cilb = att-se*agg_se$crit_val)]

  
  return(result)
}

# scheme ------------------------------------------------------------------------

#scheme for aggregation
get_agg_sch <- function(gt_result, aux, p){

  #create group_time
  id_dt <- data.table(weight = aux$weights/sum(aux$weights), G = aux$dt_inv[, G])
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group_time <- gt_result$gt |> merge(pg_dt, by = "G")
  group_time[, mg := ming(G)]
  setorder(group_time, time, mg) #change the order to match the order in gtatt
  if(!all(names(gt_result$att) == group_time[, paste0(G, ".", time)])){stop("some bug makes gt misaligned, please report this to the maintainer. Thanks.")}
  
  #get the event-specific matrix, and available ggts
  if(p$exper$event_specific){
    es <- get_es_scheme(group_time, aux, p)
    group_time <- es$group_time #some gt may not have availble effect (ex: g1 == g2)
    es_weight <- es$es_weight
  }

  #choose the target based on aggregation type
  tg <- get_agg_targets(group_time, p)
  group_time <- tg$group_time
  targets <- tg$targets
  
  
  
  agg_weights <- data.table()
  for(tar in targets){ #the order matters

    group_time[, weight := 0] #weight is 0 if not a target
    group_time[target == tar & used, weight := pg/sum(pg)] 
    target_weights <- group_time[, .(weight)] |> transpose()

    agg_weights <- rbind(agg_weights, target_weights)
  }
  
  agg_weights <- as.matrix(agg_weights)
  
  
  #apply the transformation on the aggregation matrix
  if(p$exper$event_specific){
    agg_weights <- agg_weights %*% as.matrix(es_weight)
  }
  
  return(list(agg_weights = agg_weights, #a matrix of each target and gt's weight in it 
              targets = targets,
              group_time = group_time))
}

#get the target parameters
get_agg_targets <- function(group_time, p){
  group_time[, post := as.numeric(ifelse(time >= g1(G), 1, -1))]
  if (p$result_type == "dynamic") {
    group_time[, target := time-g1(G)]
  } else if (p$result_type == "group") {
    group_time[, target := g1(G)*post] # group * treated
  } else if (p$result_type == "time") {
    group_time[, target := time*post] #calendar time * treated
  } else if (p$result_type == "simple") {
    group_time[, target := post] #treated / not treated
  } else if(p$result_type == "group_time"){
    group_time[, target := paste0(g1(G), ".", time)]
  } else if(p$result_type == "group_group_time"){
    group_time[, target := paste0(G, ".", time)]
  } else if(p$result_type == "dynamic_sq"){
    group_time[, target := paste0(time-g1(G), ".", time-g2(G))]
  }
  
  targets <- group_time[, unique(target)]
  
  #for balanced cohort composition in dynamic setting
  #a cohort us only used if it is seen for all dynamic time
  if(p$result_type == "dynamic" & !is.na(p$balanced_event_time)){
    
    cohorts <- group_time[, .(max_et = max(time-G),
                              min_et = min(time-G)), by = "G"]
    cohorts[, used := max_et >= p$balanced_event_time] #the max
    if(!cohorts[, any(used)]){stop("balanced_comp_range outside avalible range")}
    group_time[, used := G %in% cohorts[used == TRUE, G]]
    
    targets <- targets[targets <= p$balanced_event_time & targets >= cohorts[used == TRUE, min(min_et)]]
    
  } else{group_time[, used := TRUE]}
  
  return(list(group_time = group_time, targets = targets))
  
}

#the scheme for getting event-specific effect
get_es_scheme <- function(group_time, aux, p){

  es_group_time <- copy(group_time) #group_time with available es effect
  es_weight <- data.table()
  for(ggt in seq_len(nrow(group_time))){

    group_time <- get_es_ggt_weight(group_time, ggt, aux, p)
    
    if(group_time[, all(weight == 0)]){ #no available stuff
      t <- group_time[ggt, time]
      gg <- group_time[ggt, G]
      es_group_time <- es_group_time[!(time == t & G == gg)] #remove the ggt from new group time
    } else {
      es_ggt_weights <- group_time[, .(weight)] |> transpose()
      es_weight <- es_weight |> rbind(es_ggt_weights)
    }

  }
  
  return(list(group_time = es_group_time, es_weight = es_weight))

}

get_es_ggt_weight <- function(group_time, ggt, aux, p){

  group_time[, weight := 0] #reset 
  t <- group_time[ggt, time]
  g1 <- group_time[ggt, g1(G)]
  g2 <- group_time[ggt, g2(G)]
  gg <- group_time[ggt, G]
  
  if(t < g2){ #direct pure effect

    group_time[ggt, weight := 1] #just use the observed effect
    
  } else if(g1 < g2) { #imputation = treat-pre + (control-post - control-pre)
    
    base_period <- g2 - 1
    if(base_period == t){return(group_time)} #need to be treated already in base period
    
    #get the cohorts
    tb <- group_time[,G == gg & time == base_period]
    c <-  group_time[,g1(G) == g1 & g2(G) > max(t,base_period) & g2(G) != g2]
    cp <- group_time[, c & time == t]
    cb <- group_time[, c & time == base_period]
    
    #if any group have no available cohort, skip
    if(sum(tb) == 0 | sum(cp) == 0 | sum(cb) == 0){return(group_time)}
    
    #assign the weights
    group_time[tb, weight := pg/sum(pg)]
    group_time[cp, weight := pg/sum(pg)]
    group_time[cb, weight := -pg/sum(pg)]

    
  } else if (g1 > g2) { #double did = (treat-post - treat-base) - (control-post - control-pre)

    base_period <- g1 - 1
    if(base_period == t){return(group_time)}

    #get the cohorts
    tp <- group_time[,.I == ggt]
    tb <- group_time[,G == gg & time == base_period]
    c <-  group_time[,g2(G) == g2 & g1(G) > max(t,base_period) & g1(G) != g1]
    cp <- group_time[, c & time == t]
    cb <- group_time[, c & time == base_period]
    
    #if any group have no available cohort, skip
    if(sum(tp) == 0 | sum(tb) == 0 | sum(cp) == 0 | sum(cb) == 0){return(group_time)}
    
    #assign the weights
    group_time[tp, weight := pg/sum(pg)]
    group_time[tb, weight := -pg/sum(pg)]
    group_time[cp, weight := -pg/sum(pg)]
    group_time[cb, weight := pg/sum(pg)]

  } 

  return(group_time)
  
}

# influence function ------------------------------------------------------------

#influence from weight calculation
get_weight_influence <- function(agg_weights, gt_result, aux, agg_sch, p) {
  
  group <- gt_result$gt[,.(G,time)]
  gt_att <- gt_result$att
  weights <- aux$weights
  id_cohorts <- aux$dt_inv[, G]
  keepers <- which(agg_weights != 0)
  
  id_dt <- data.table(weight = weights/sum(weights), G = id_cohorts)
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group <- group |> merge(pg_dt, by = "G")
  
  group[, time := as.integer(time)]
  
  if(is.na(p$exper$cohortvar2)){
    group[, G := as.integer(G)]
    setorder(group, time, G)
  } else {
    group[, mg := ming(G)]
    group[, G1 := g1(G)]
    group[, G2 := g2(G)]
    setorder(group, time, mg, G1, G2) #sort
  }

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

