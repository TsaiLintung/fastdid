#2025-05-09
message('loading fastdid source ver. ver: 1.0.5 date: 2025-05-09')
require(data.table);
 require(stringr);
 require(BMisc);
 require(collapse);
 require(dreamerr);
 require(parglm);
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
  
  att <- gt_result$att
  inf_func <- gt_result$inf_func 
  
  if(p$event_specific & !is.na(p$cohortvar2)){
    es_weight <- agg_sch$es_sto_weight+agg_sch$es_det_weight
    es_inf_weights <- get_weight_influence(att, agg_sch$pre_es_group_time, agg_sch$es_sto_weight, aux, p)
    att <- (es_weight) %*% att
    inf_func <- (inf_func %*% t(es_weight)) + es_inf_weights
  }
  
  #get att
  agg_att <- agg_sch$agg_weights %*% att
  
  #get influence function matrix
  
  inf_weights <- get_weight_influence(att, agg_sch$group_time, agg_sch$agg_weights, aux, p)
  inf_matrix <- (inf_func %*% t(agg_sch$agg_weights)) + inf_weights 

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
    pre_es_group_time <- group_time
    pre_es_group_time[, pg := NULL]
    group_time <- es$group_time #some gt may not have availble effect (ex: g1 == g2)
    es_det_weight <- as.matrix(es$es_det_weight)
    es_sto_weight <- as.matrix(es$es_sto_weight)
  } else {
    es_det_weight <- NULL
    es_sto_weight <- NULL
    pre_es_group_time <- NULL
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
              pre_es_group_time = pre_es_group_time,
              es_det_weight = es_det_weight,
              es_sto_weight = es_sto_weight))
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

get_weight_influence <- function(att, group, agg_weights, aux, p){
  
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
    inf_weights <- sapply(asplit(agg_weights, 1), function (x){
      get_weight_influence_param(x, group, att, aux, p)
    })
  } else {
    inf_weights <- matrix(unlist(mclapply(asplit(agg_weights, 1), function (x){
      get_weight_influence_param(x, group, att, aux, p)
    })), ncol = dim(agg_weights)[1])
  }
  
  return(inf_weights)
  
}

#influence from weight calculation
get_weight_influence_param <- function(agg_weights, group, gt_att, aux, p) {
  keepers <- which(agg_weights != 0)
  group <- group[keepers,]
  if(nrow(group) == 0){return(rep(0, length(aux$weights)))} #for direct double did 

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


# auxilary steps in the main fastdid function

get_exper_default <- function(exper, exper_args){
  for(arg in exper_args){
    if(is.null(exper[[arg]])){
      exper[[arg]] <- NA
    }
  }
  
  if(!is.na(exper$only_balance_2by2) & exper$only_balance_2by2){ #will create this col in the get_aux part
    exper$filtervar <- "no_na"
    exper$filtervar_post <- "no_na"
  }
  
  return(exper)
}

coerce_dt <- function(dt, p){
  
  if(!is.na(p$cohortvar2)){return(coerce_dt_doub(dt, p))} #in doubledid.R
  
  #chcek if there is availble never-treated group
  if(!is.infinite(dt[, max(G)])){
    if(p$control_option == "both"){warning("no never-treated availble, effectively using not-yet-but-eventually-treated as control")}
    if(p$control_option == "never"){stop("no never-treated availble.")}
  }
  
  if(p$allow_unbalance_panel){
    dt_inv_raw <- dt[dt[, .I[1], by = unit]$V1]
    setorder(dt_inv_raw, G)
    dt_inv_raw[, new_unit := seq_len(.N)] #let unit start from 1 .... N, useful for knowing which unit is missing
    dt <- dt |> merge(dt_inv_raw[,.(unit, new_unit)], by = "unit", sort = FALSE)
    dt[, unit := new_unit]
  }
  
  setorder(dt, time, G, unit) #sort the dataset essential for the sort-once-quick-access 
  
  #deal with time, coerice time to 1,2,3,4,5.......
  time_periods <- dt[, unique(time)]
  time_size <- length(time_periods)
  
  #TODO: this part is kinda ugly
  time_offset <- min(time_periods) - 1 #assume time starts at 1, first is min after sort :)
  gcol <- str_subset(names(dt), ifelse(is.na(p$cohortvar2), "G", "G1|G2")) 
  if(time_offset != 0){
    dt[, G := G-time_offset]
    
    dt[, time := time-time_offset]
    time_periods <- time_periods - time_offset
  }
  
  time_step <- 1 #time may not jump at 1
  if(any(time_periods[seq(2,length(time_periods),1)] - time_periods[seq_len(length(time_periods))-1] != 1)){
    time_step <- time_periods[2]-time_periods[1]
    time_periods <- (time_periods-1)/time_step+1
    if(any(time_periods[seq(2,length(time_periods),1)] - time_periods[seq_len(length(time_periods))-1] != 1)){stop("time step is not uniform")}
    dt[G != 1, G := (G-1)/time_step+1]
    
    dt[time != 1, time := (time-1)/time_step+1]
  }
  
  #add the information to t
  t <- list()
  t$time_step <- time_step
  t$time_offset <- time_offset
  
  if(nrow(dt) == 0){
    stop("no data after coercing the dataset")
  }
  
  return(list(dt = dt, p = p, t = t))
  
}

get_auxdata <- function(dt, p){
  
  time_periods <- dt[, unique(time)]
  id_size <- dt[, uniqueN(unit)]
  
  #construct the outcomes list for fast access later
  #loop for multiple outcome
  outcomes_list <- list()
  for(outcol in p$outcomevar){
    outcomes <- list()
    
    if(!p$allow_unbalance_panel){
      for(i in time_periods){
        start <- (i-1)*id_size+1
        end <- i*id_size
        outcomes[[i]] <- dt[seq(start,end), get(outcol)]
      }
    } else {
      
      for(i in time_periods){
        #populate a outcome vector of length N with outcome data in the right place
        #NA is data gone or missing, will be addressed in estimate_did_rc
        outcome_period <- rep(NA, id_size)
        data_pos <- dt[time == i, unit] #units observed in i
        outcome_period[data_pos] <- dt[time == i, get(outcol)]
        outcomes[[i]] <- outcome_period
      }
      
    }
    
    outcomes_list[[outcol]] <- outcomes
  }
  
  #the time-invariant parts 
  if(!p$allow_unbalance_panel){
    dt_inv <- dt[seq_len(id_size)]
  } else {
    dt_inv <- dt[dt[, .I[1], by = unit]$V1] #the first observation
    setorder(dt_inv, unit) #can't move this outside
  }
  
  cohorts <- dt_inv[, unique(G)]
  cohort_sizes <- dt_inv[, .(cohort_size = .N) , by = G]
  
  # the optional columns
  varycovariates <- list()
  if(!allNA(p$varycovariatesvar)){
    for(i in time_periods){
      start <- (i-1)*id_size+1
      end <- i*id_size
      varycovariates[[i]] <- dt[seq(start,end, by = 1), .SD, .SDcols = p$varycovariatesvar]
    }
  } else {
    varycovariates <- NA
  }
  
  #create na indicator for filtering
  if(!is.na(p$exper$only_balance_2by2) & p$exper$only_balance_2by2){
    if("no_na" %in% names(dt)){stop("no_na is already in dt, consider using another column name")}
    varnames <- unlist(p[str_ends(names(p), "var")], recursive = TRUE) #get all the argument that ends with "var"
    varnames <- varnames[!varnames %in% c(p$timevar, p$unitvar, p$cohortvar) & !is.na(varnames) & !is.null(varnames)]
    dt[, no_na := TRUE]
    for(col in varnames){ #nothing can be 
      if(is.na(col)){next}
      dt[is.na(get(col)), no_na := TRUE]
    }
  }
  
  # filters
  filters <- list()
  if(!is.na(p$exper$filtervar)){
    for(t in time_periods){
      #filters[[t]] <- rep(NA, id_size)
      #data_pos <- dt[time == t, unit] #units observed in i
      filters[[t]] <- unlist(dt[time == t,  .SD, .SDcols = p$exper$filtervar])
      if(p$allow_unbalance_panel){stop("unbalance panel not supported with filtervar")}
    }
  } else {
    filters <- NA
  }
  
  if(!allNA(p$covariatesvar)){
    covariates <- dt_inv[,.SD, .SDcols = p$covariatesvar]
  } else {
    covariates <- NA
  }
  
  if(!is.na(p$clustervar)){
    cluster <- dt_inv[, .SD, .SDcols = p$clustervar] |> unlist()
  } else {
    cluster <- NA
  }
  
  if(!is.na(p$weightvar)){
    weights <- dt_inv[, .SD, .SDcols = p$weightvar] |> unlist()
    weights <- weights/mean(weights) #normalize
  } else {
    weights <- rep(1, id_size)
  }
  
  aux <- as.list(environment())
  aux$dt <- NULL
  aux$p <- NULL
  class(aux) <- "locked"
  
  return(aux)
  
}

convert_targets <- function(results, p, t){
  
  if(!is.na(p$exper$aggregate_scheme)){return(results)}  #no conversion back if use custom
  
  switch(p$result_type,
         dynamic = {
           results[, event_time := target]
           setcolorder(results, "event_time", before = 1)
         },
         group = {
           results[, type := ifelse(target >= 0, "post", "pre")]
           results[, cohort := recover_time(abs(target), t)]
           setcolorder(results, "cohort", before = 1)
         },
         time = {
           results[, type := ifelse(target >= 0, "post", "pre")]
           results[, time := recover_time(abs(target), t)]
           setcolorder(results, "time", before = 1)
         },
         simple = {
           results[, type := ifelse(target >= 0, "post", "pre")]
         },
         group_time = {
           results[, cohort := as.numeric(str_split_i(target, "\\.", 1))]
           results[, time :=  as.numeric(str_split_i(target, "\\.", 2))]
           
           #recover the time
           results[, cohort := recover_time(cohort, t)]
           results[, time := recover_time(time, t)]
           
         },
         group_group_time = {
           results[, cohort := str_split_i(target, "\\.", 1)]
           results[, time :=  as.numeric(str_split_i(target, "\\.", 2))]
           
           results[, cohort1 := g1(cohort)]
           results[, cohort2 := g2(cohort)]
           
           results[, cohort1 := recover_time(cohort1, t)]
           results[, cohort2 := recover_time(cohort2, t)]
           results[, time := recover_time(time, t)]
           results[, `:=`(cohort = NULL)]
         },
         dynamic_stagger = {
           results[, event_time_1 :=  as.numeric(str_split_i(target, "\\.", 1))]
           results[, event_stagger :=  as.numeric(str_split_i(target, "\\.", 2))]
         }
  )
  
  results[, target := NULL]
  
  return(results)
}

# small stuff ---------

recover_time <- function(time, t){
  return(((time-1)*t$time_step)+1+t$time_offset)
}

#locked list
#from: https://stackoverflow.com/questions/58776481/make-r-function-return-a-locked-immutable-list
.S3method("[[<-", "locked", function(value) {stop("Can't assign into locked object")})
.S3method("[<-", "locked", function(value) {stop("Can't assign into locked object")})
.S3method("$<-", "locked", function(value) {stop("Can't assign into locked object")})
# a <- list(b = 1, c = 2)
# class(a) <- c("locked")



# utils -----------------------------------

#g11 <- c(1,1,2,3,4)
#g22 <- c(2,1,3,2,4)
#GG <- as.factor(paste0(g11, ".", g22))

g1 <- function(GG){
  if(is.numeric(GG)){return(GG)}
  return(as.numeric(str_split_i(GG, "-", 1)))
}

g2 <- function(GG){
  return(as.numeric(str_split_i(GG, "-", 2)))
}

ming <- function(GG){
  if(is.numeric(GG)){return(GG)}
  else {pmin(g1(GG), g2(GG))}
}

# overiden function -------------------------------------------------

coerce_dt_doub <- function(dt, p){
  
  setnames(dt, "G", "G1")
  dt[, mg := pmin(G1, G2)]
  setorder(dt, time, mg, G1, G2, unit)  #for sort one quick access
 
  #check if there is availble never-treated group
  if(!is.infinite(dt[, max(mg)])){
    if(p$control_option == "both"){warning("no never-treated availble, effectively using not-yet-but-eventually-treated as control")}
    if(p$control_option == "never"){stop("no never-treated availble.")}
  }
  
  if(p$allow_unbalance_panel){ #let unit start from 1 .... N, useful for knowing which unit is missing
    dt_inv_raw <- dt[dt[, .I[1], by = unit]$V1]
    setorder(dt_inv_raw, mg, G1, G2)
    dt_inv_raw[, new_unit := 1:.N] 
    dt <- dt |> merge(dt_inv_raw[,.(unit, new_unit)], by = "unit", sort = FALSE)
    dt[, unit := new_unit]
  }
  
  #deal with time, coerice time to 1,2,3,4,5.......
  time_periods <- dt[, unique(time)]
  time_size <- length(time_periods)

  time_offset <- min(time_periods) - 1 #assume time starts at 1, first is min after sort :)
  gcol <- c("G1", "G2")
  if(time_offset != 0){
    dt[, c(gcol) := .SD-time_offset, .SDcols = gcol]
    dt[, time := time-time_offset]
    time_periods <- time_periods - time_offset
  }

  time_step <- 1 #time may not jump at 1
  if(any(time_periods[2:length(time_periods)] - time_periods[1:length(time_periods)-1] != 1)){
    time_step <- time_periods[2]-time_periods[1]
    time_periods <- (time_periods-1)/time_step+1
    if(any(time_periods[2:length(time_periods)] - time_periods[1:length(time_periods)-1] != 1)){stop("time step is not uniform")}
    
    for(g in gcol){
      dt[get(g) != 1, c(g) := (get(g)-1)/time_step+1]
    }
    
    dt[time != 1, time := (time-1)/time_step+1]
  }
  dt[, mg := pmin(G1, G2)]
  dt[, G := paste0(G1, "-", G2)] #create G once its finalized
  
  #add the information to t
  t <- list()
  t$time_step <- time_step
  t$time_offset <- time_offset
  
  if(nrow(dt) == 0){
    stop("no data after coercing the dataset")
  }
  
  return(list(dt = dt, p = p, t = t))
  
}

# aggregation scheme -----------------------------------------------------------

#the scheme for getting event-specific effect
get_es_scheme <- function(group_time, aux, p){
  
  es_group_time <- copy(group_time) #group_time with available es effect
  #create lookup
  es_group_time[, mg := ming(G)]
  es_group_time[, G1 := g1(G)]
  es_group_time[, G2 := g2(G)]
  es_weight_list <- list()
  
  ggt <- as.list(seq_len(nrow(group_time)))
  if(!p$parallel){
    es_weight_list <- lapply(ggt, get_es_ggt_weight, group_time, aux, p)
  } else {
    es_weight_list <- mclapply(ggt, get_es_ggt_weight, group_time, aux, p, mc.cores = getDTthreads())
  }
  
  valid_ggt <- which(!sapply(es_weight_list, is.null))
  es_group_time <- es_group_time[valid_ggt] #remove the ones without
  es_weight_list <- es_weight_list[valid_ggt]
  es_det_weight <- do.call(rbind, lapply(es_weight_list, \(x){x$det})) 
  es_sto_weight <- do.call(rbind, lapply(es_weight_list, \(x){x$sto}))
  
  return(list(group_time = es_group_time, es_det_weight = es_det_weight, es_sto_weight = es_sto_weight))
  
}

#get the scheme for retriving group-group-time estimates
get_es_ggt_weight <- function(ggt, group_time, aux, p){
  
  group_time <- copy(group_time) #avoid accidental modification
  
  group_time[, det_weight := 0] #reset 
  group_time[, sto_weight := 0] #reset 
  t <- group_time[ggt, time]
  g1 <- group_time[ggt, G1]
  g2 <- group_time[ggt, G2]
  gg <- group_time[ggt, G]
  
  if(is.infinite(g1)){return(NULL)}
  
  if(t < g2){ #direct pure effect
    
    group_time[ggt, det_weight := 1] #just use the observed effect
    
  } else if(g1 < g2) { #imputation = treat-pre + (control-post - control-pre)
    
    base_period <- g2 - 1
    if(base_period == t){return(NULL)}
    min_control_cohort <- ifelse(p$double_control_option == "never", Inf, max(t,base_period)+1)
    
    #get the cohorts
    tb <- group_time[,G == gg & time == base_period]
    c <-  group_time[, G1 == g1 & G2 >= min_control_cohort & G2 != g2]
    if(p$control_option == "notyet"){
      c[group_time[, is.infinite(G2)]] <- FALSE
    }
    cp <- group_time[, c & time == t]
    cb <- group_time[, c & time == base_period]
    
    #if any group have no available cohort, skip
    if(sum(tb) == 0 | sum(cp) == 0 | sum(cb) == 0){return(NULL)}
    
    #assign the weights
    group_time[tb, det_weight := 1]
    group_time[cp, sto_weight := pg/sum(pg)]
    group_time[cb, sto_weight := -pg/sum(pg)]
    
    
  } else if (g1 > g2) { #double did = (treat-post - treat-base) - (control-post - control-pre)
    
    base_period <- g1 - 1
    if(base_period == t){return(NULL)}
    min_control_cohort <- ifelse(p$double_control_option == "never", Inf, max(t,base_period)+1)
    
    #get the cohorts
    tp <- group_time[,.I == ggt]
    tb <- group_time[,G == gg & time == base_period]
    c <-  group_time[,G2 == g2 & G1 >= min_control_cohort & G1 != g1]
    if(p$control_option == "notyet"){
      c[group_time[, is.infinite(G1)]] <- FALSE
    }
    cp <- group_time[, c & time == t]
    cb <- group_time[, c & time == base_period]
    
    #if any group have no available cohort, skip
    if(sum(tp) == 0 | sum(tb) == 0 | sum(cp) == 0 | sum(cb) == 0){return(NULL)}
    
    #assign the weights
    group_time[tp, det_weight := 1]
    group_time[tb, det_weight := -1]
    group_time[cp, sto_weight := -pg/sum(pg)]
    group_time[cb, sto_weight := pg/sum(pg)]
    
  } 
  
  if(all(group_time[, det_weight+sto_weight] == 0)){return(NULL)} #not redundant!
  return(list(det = group_time[, det_weight], sto = group_time[, sto_weight]))
  
}

estimate_did <- function(dt_did, covvars, p, cache){
  
  #estimate did
  param <- as.list(environment())
  if(!p$allow_unbalance_panel){
    result <- do.call(estimate_did_bp, param)
  } else {
    result <- do.call(estimate_did_rc, param)
  }
  return(result)
}

estimate_did_bp <- function(dt_did, covvars, p, cache){
  
  # preprocess --------
  oldn <- dt_did[, .N]
  data_pos <-  which(dt_did[, !is.na(D)])
  dt_did <- dt_did[data_pos]
  n <- dt_did[, .N]
  
  if(is.matrix(covvars)){
    ipw <- p$control_type %in% c("ipw", "dr") 
    or <- p$control_type %in% c("reg", "dr")
    covvars <- covvars[data_pos,] 
  } else {
    ipw <- FALSE
    or <- FALSE
  }
  
  # ipw --------

  if(ipw){
    if(is.null(cache)){ #if no cache, calcuate ipw
      #estimate the logit
      prop_score_est <- suppressWarnings(parglm::parglm.fit(covvars, dt_did[, D],
                                                    family = stats::binomial(), 
                                                    weights = dt_did[, weights],
                                                    control = parglm.control(nthreads = ifelse(p$parallel, 1, getDTthreads())), #no parallel if already parallel
                                                    intercept = FALSE))
      class(prop_score_est) <- "glm" #trick the vcov function to think that this is a glm object to dispatch the write method
      #const is implicitly put into the ipw formula, need to incorporate it manually

      logit_coef <-  prop_score_est$coefficients
      
      if(anyNA(logit_coef)){
        warning("some propensity score estimation resulted in NA coefficients, likely cause by perfect colinearity")
      }
      
      logit_coef[is.na(logit_coef)|abs(logit_coef) > 1e10] <- 0 #put extreme value and na to 0
      prop_score_fit <- fitted(prop_score_est)
      if(max(prop_score_fit) >= 1-1e-10){warning(paste0("extreme propensity score: ", max(prop_score_fit), ", support overlap is likely to be violated"))} #<=0 (only in control) is fine for ATT since it is just not used 
      prop_score_fit <- pmin(1-1e-10, prop_score_fit) #for the ipw

      hess <- stats::vcov(prop_score_est) * n #for the influence function
      hess[is.na(hess)|abs(hess) > 1e10] <- 0
      
      cache <- list(hess = hess, ps = prop_score_fit)

    } else { #when using multiple outcome, ipw cache can be reused
      hess <- cache$hess
      prop_score_fit <- cache$ps
      logit_coef <- NA #won't be needing the approximate cache
    }

    #get the results into the main did dt
    dt_did[, ps := prop_score_fit]
    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*ps*(1-D)/(1-ps)]

  } else {

    prop_score_fit <- rep(1,n)
    logit_coef <- NA
    hess <- NA

    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*(1-D)]

  }

  #delta y is needed for PR
  dt_did[, delta_y := post.y-pre.y]

  # or --------

  if(or){

    #TODO: this should be optimized with better backend and some caching

    #should change to speedlm or something

    control_bool <- dt_did[, D==0]
    reg_coef <- stats::coef(stats::lm.wfit(x = covvars[control_bool,], y = dt_did[control_bool,delta_y],
                                           w = dt_did[control_bool,weights]))

    if(anyNA(reg_coef)){
      stop("some outcome regression resulted in NA coefficients, likely cause by perfect colinearity")
    }

    #the control function from outcome regression
    dt_did[, or_delta := as.vector(tcrossprod(reg_coef, covvars))]

  } else {
    dt_did[, or_delta := 0]
  }

  #did --------

  dt_did[, att_treat := treat_ipw_weight*(delta_y-or_delta)] #minus the OR adjust
  dt_did[, att_cont := cont_ipw_weight*(delta_y-or_delta)]

  weighted_treat_delta <- dt_did[,sum(att_treat)/sum(treat_ipw_weight)]
  weighted_cont_delta <- dt_did[,sum(att_cont)/sum(cont_ipw_weight)]

  att <- weighted_treat_delta - weighted_cont_delta
  
  # influence --------

  # influence from ipw
  if(ipw){

    M2 <- colMeans(dt_did[, cont_ipw_weight*(delta_y-weighted_cont_delta-or_delta)] * covvars)

    score_ps <- dt_did[, weights*(D-ps)] * covvars
    
    asym_linear_ps <- score_ps %*% hess

    #ipw for control
    inf_cont_ipw <- asym_linear_ps %*% as.matrix(M2)

  } else {inf_cont_ipw <- 0}

  if(or){

    M1 <- colMeans(dt_did[, treat_ipw_weight] * covvars)
    M3 <- colMeans(dt_did[, cont_ipw_weight] * covvars)

    or_x <- dt_did[, weights*(1-D)] * covvars
    or_ex <- dt_did[, weights*(1-D)*(delta_y - or_delta)] * covvars
    XpX <- crossprod(or_x, covvars)/n

    #calculate alrw = eX (XpX)^-1 by solve XpX*alrw = ex, much faster since avoided inv
    asym_linear_or <- t(solve(XpX, t(or_ex)))# |> reverse_col()

    #or for treat
    inf_treat_or <- -asym_linear_or %*% M1 #a negative sign here, since or_delta is subtracted from the att

    #or for control
    inf_cont_or <- -asym_linear_or %*% M3

  } else {
    inf_treat_or <- 0
    inf_cont_or <- 0
  }

  # influence from did
  inf_cont_did <- dt_did[, att_cont - cont_ipw_weight*weighted_cont_delta]
  inf_treat_did <-  dt_did[, (att_treat - treat_ipw_weight*weighted_treat_delta)]


  #get overall influence function
  #if(dt_did[, mean(cont_ipw_weight)] < 1e-10){warning("little/no overlap in covariates between control and treat group, estimates are unstable.")}
  inf_cont <- (inf_cont_did+inf_cont_ipw+inf_cont_or)/dt_did[, mean(cont_ipw_weight)]
  inf_treat <- (inf_treat_did+inf_treat_or)/dt_did[,mean(treat_ipw_weight)]
  inf_func_no_na <- inf_treat - inf_cont

  #post process (fill zeros for irrelevant ones)
  inf_func <- rep(0, oldn) #the default needs to be 0 for the matrix multiplication
  inf_func_no_na <- inf_func_no_na * oldn / n #adjust the value such that mean over the whole id size give the right result
  inf_func[data_pos] <- inf_func_no_na
  
  return(list(att = att, inf_func = inf_func, cache = list(ps = prop_score_fit, hess = hess))) #for next outcome
}

estimate_did_rc <- function(dt_did, covvars, p, cache){
  
  #TODO: skip if not enough valid data
  
  # preprocess --------
  
  
  oldn <- dt_did[, .N]
  data_pos <-  which(dt_did[, !is.na(D)])
  dt_did <- dt_did[data_pos]
  n <- dt_did[, .N]
  
  #separate the dataset into pre and post
  dt_did[, inpre := as.numeric(!is.na(pre.y))]
  dt_did[, inpost := as.numeric(!is.na(post.y))]
  n_pre <- dt_did[, sum(!is.na(pre.y))]
  n_post <- dt_did[, sum(!is.na(post.y))]
  
  sum_weight_pre <- dt_did[, sum(inpre*weights)]
  sum_weight_post <- dt_did[, sum(inpost*weights)]
  
  if(is.matrix(covvars)){
    ipw <- p$control_type %in% c("ipw", "dr") 
    or <- p$control_type %in% c("reg", "dr")
    covvars <- covvars[data_pos,] 
  } else {
    ipw <- FALSE
    or <- FALSE
  }
  
  # ipw --------
  
  if(ipw){
    
    #no caching since composition changes by period
    
    #estimate the logit
    prop_score_est <- suppressWarnings(parglm.fit(covvars, dt_did[, D],
                                                  family = stats::binomial(),
                                                  weights = dt_did[, weights*(inpre+inpost)*n/(n_pre+n_post)], #when seen in both pre and post have double weight
                                                  control = parglm.control(nthreads = ifelse(p$parallel, 1, getDTthreads())),
                                                  intercept = FALSE)) #*(inpre+inpost)
    class(prop_score_est) <- "glm" #trick the vcov function to think that this is a glm object to dispatch the write method
    #const is implicitly put into the ipw formula, need to incorporate it manually
    
    #for the influence, will be cached
    hess <- stats::vcov(prop_score_est) * n #for the influence function
    
    logit_coef <-  prop_score_est$coefficients 
    logit_coef[is.na(logit_coef)|abs(logit_coef) > 1e10] <- 0 #put extreme value and na to 0
    prop_score_fit <- fitted(prop_score_est)
    if(max(prop_score_fit) >= 1){warning(paste0("support overlap condition violated for some group_time"))}
    prop_score_fit <- pmin(1-1e-16, prop_score_fit) #for the ipw
    
    #get the results into the main did dt
    dt_did[, ps := prop_score_fit]
    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*ps*(1-D)/(1-ps)]
    
  } else {
    
    prop_score_fit <- rep(1,n)
    logit_coef <- NA
    hess <- NA
    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*(1-D)]
    
  }
  
  # or --------
  
  if(or){
    
    control_bool_post <- dt_did[, D==0 & inpost] #control group and have obs in post period
    control_bool_pre <- dt_did[, D==0 & inpre]
    reg_coef_post <- stats::coef(stats::lm.wfit(x = covvars[control_bool_post,], y = dt_did[control_bool_post,post.y],
                                                w = dt_did[control_bool_post,weights]))

    reg_coef_pre <- stats::coef(stats::lm.wfit(x = covvars[control_bool_pre,], y = dt_did[control_bool_pre,pre.y],
                                               w = dt_did[control_bool_pre,weights]))

    if(anyNA(reg_coef_post)|anyNA(reg_coef_pre)){
      stop("some outcome regression resulted in NA coefficients, likely cause by perfect colinearity")
    }

    #the control function from outcome regression
    dt_did[, or_delta_post := as.vector(tcrossprod(reg_coef_post, covvars))]
    dt_did[, or_delta_pre := as.vector(tcrossprod(reg_coef_pre, covvars))]
    
  } else {
    dt_did[, or_delta_post := 0]
    dt_did[, or_delta_pre := 0]
  }
  
  #did --------
  
  #mean weight
  mean_wcpo <-  dt_did[,sum(cont_ipw_weight*inpost)/n_post]
  mean_wtpo <-  dt_did[,sum(treat_ipw_weight*inpost)/n_post]
  mean_wcpr <-  dt_did[,sum(cont_ipw_weight*inpre)/n_pre]
  mean_wtpr <-  dt_did[,sum(treat_ipw_weight*inpre)/n_pre]
  
  #delta y is needed for PR
  dt_did[, att_treat_post := treat_ipw_weight*(post.y-or_delta_post)/mean_wtpo] #minus the OR adjust
  dt_did[, att_cont_post :=  cont_ipw_weight*(post.y-or_delta_post)/mean_wcpo]
  dt_did[, att_treat_pre := treat_ipw_weight*(pre.y-or_delta_pre)/mean_wtpr] #minus the OR adjust
  dt_did[, att_cont_pre := cont_ipw_weight*(pre.y-or_delta_pre)/mean_wcpr]
  
  weighted_treat_post <- dt_did[,mean(att_treat_post, na.rm = TRUE)]
  weighted_cont_post <- dt_did[,mean(att_cont_post, na.rm = TRUE)]
  weighted_treat_pre <- dt_did[,mean(att_treat_pre, na.rm = TRUE)]
  weighted_cont_pre <- dt_did[,mean(att_cont_pre, na.rm = TRUE)]
  
  att <- (weighted_treat_post - weighted_treat_pre) - (weighted_cont_post - weighted_cont_pre)
  
  # influence --------
  
  # influence from ipw
  if(ipw){
    
    # a bit unsure about this part
    M2_post <- colSums(dt_did[, inpost*cont_ipw_weight*(post.y-weighted_cont_post-or_delta_post)/n] * covvars, na.rm = TRUE) / mean_wcpo
    M2_pre <- colSums(dt_did[, inpre*cont_ipw_weight*(pre.y-weighted_cont_pre-or_delta_pre)/n] * covvars, na.rm = TRUE) / mean_wcpr
    
    #not sure about /2
    score_ps <- dt_did[, weights*(inpre+inpost)*n/(n_pre+n_post)*(D-ps)] * covvars#weight is doubled for observed in both post and pre
    asym_linear_ps <- score_ps %*% hess 
    
    #ipw for control
    inf_cont_ipw_post <- asym_linear_ps %*% M2_post 
    inf_cont_ipw_pre  <- asym_linear_ps %*% M2_pre
    
    
  } else {
    inf_cont_ipw_post <- 0
    inf_cont_ipw_pre <- 0
  }
  
  
  if(or){

    M1_post <- colSums(dt_did[, inpost*treat_ipw_weight/n] * covvars, na.rm = TRUE) / mean_wtpo
    M1_pre <- colSums(dt_did[, inpre*treat_ipw_weight/n] * covvars, na.rm = TRUE) / mean_wtpr
    M3_post <- colSums(dt_did[, inpost*cont_ipw_weight/n] * covvars, na.rm = TRUE) / mean_wcpo
    M3_pre <- colSums(dt_did[, inpre*cont_ipw_weight/n] * covvars, na.rm = TRUE) / mean_wcpr

    or_x_post <- dt_did[, inpost*weights*(1-D)] * covvars
    or_x_pre <- dt_did[, inpre*weights*(1-D)] * covvars
    or_ex_post <- dt_did[, inpost*weights*(1-D)*(post.y - or_delta_post)] * covvars
    or_ex_pre <- dt_did[, inpre*weights*(1-D)*(pre.y - or_delta_pre)] * covvars
    XpX_post <- crossprod(or_x_post, covvars)/n_post
    XpX_pre <- crossprod(or_x_pre, covvars)/n_pre

    #calculate alrw = eX (XpX)^-1 by solve XpX*alrw = ex, much faster since avoided inv
    asym_linear_or_post <- t(solve(XpX_post, t(or_ex_post)))
    asym_linear_or_pre <- t(solve(XpX_pre, t(or_ex_pre)))

    #or for treat
    inf_treat_or_post <- -asym_linear_or_post %*% M1_post #a negative sign here, since or_delta is subtracted from the att, THE PROBLEM
    inf_treat_or_pre <- -asym_linear_or_pre %*% M1_pre

    #or for control
    inf_cont_or_post <- -asym_linear_or_post %*% M3_post
    inf_cont_or_pre <- -asym_linear_or_pre %*% M3_pre
    
  } else {
    inf_treat_or_post <- 0
    inf_treat_or_pre <- 0
    inf_cont_or_post <- 0
    inf_cont_or_pre <- 0
  }
  
  # influence from did
  inf_cont_did_post <- dt_did[, att_cont_post - cont_ipw_weight*inpost*weighted_cont_post/mean_wcpo]
  inf_treat_did_post <-  dt_did[, att_treat_post - treat_ipw_weight*inpost*weighted_treat_post/mean_wtpo]
  inf_cont_did_pre <- dt_did[, att_cont_pre - cont_ipw_weight*inpre*weighted_cont_pre/mean_wcpr]
  inf_treat_did_pre <-  dt_did[, att_treat_pre -  treat_ipw_weight*inpre*weighted_treat_pre/mean_wtpr]
  
  #fill zero to avoid NA from addition
  inf_cont_did_post[is.na(inf_cont_did_post)] <- 0
  inf_treat_did_post[is.na(inf_treat_did_post)] <- 0
  inf_cont_did_pre[is.na(inf_cont_did_pre)] <- 0
  inf_treat_did_pre[is.na(inf_treat_did_pre)] <- 0
  
  #get overall influence function
  inf_cont_post <- inf_cont_did_post+inf_cont_ipw_post+inf_cont_or_post
  inf_treat_post <- inf_treat_did_post+inf_treat_or_post
  inf_cont_pre <- inf_cont_did_pre+inf_cont_ipw_pre+inf_cont_or_pre
  inf_treat_pre <- inf_treat_did_pre+inf_treat_or_pre
  
  #post process
  inf_func_no_na_post <- (inf_treat_post - inf_cont_post) * oldn / n_post #adjust the value such that mean over the whole id size give the right result
  inf_func_no_na_post[is.na(inf_func_no_na_post)] <- 0 #fill 0 for NA part (no influce if not in this gt)
  
  inf_func_no_na_pre <- (inf_treat_pre - inf_cont_pre) * oldn / n_pre #adjust the value such that mean over the whole id size give the right result
  inf_func_no_na_pre[is.na(inf_func_no_na_pre)] <- 0
  
  inf_func <- rep(0, oldn) #the default needs to be 0 for the matrix multiplication
  inf_func[data_pos] <- inf_func_no_na_post - inf_func_no_na_pre

  return(list(att = att, inf_func = inf_func, cache = list(ps = prop_score_fit, hess = hess))) #for next outcome
}

# utilities ------

reverse_col <- function(x){
  return(x[,ncol(x):1])
}


estimate_gtatt <- function(aux, p){
  caches <- list()
  outcome_results <- list()
  for(outcol in p$outcomevar){
    y <- aux$outcomes_list[[outcol]]
    out_result <- estimate_gtatt_outcome(y, aux, p, caches)
    out_result$est$outname <- outcol
    outcome_results[[outcol]] <- out_result$est
    caches <- out_result$caches
  }
  return(outcome_results)
}

#gtatt for each outcome
estimate_gtatt_outcome <- function(y, aux, p, caches) {
    
    treated_cohort <- aux$cohorts[!is.infinite(ming(aux$cohorts))] #otherwise would try to calculate the pre-period of nevertreated in varying base period lol
    gt_all <- expand.grid(g = treated_cohort, t = aux$time_periods, stringsAsFactors = FALSE) |> transpose() |> as.list() #first loop t then g
    
    #main estimation 
    if(!p$parallel){
      gt_results <- lapply(gt_all, estimate_gtatt_outcome_gt, y, aux, p, caches)
    } else {
      gt_results <- mclapply(gt_all, estimate_gtatt_outcome_gt, y, aux, p, caches)
    }
    
    #post process
    gt_results <- gt_results[which(!sapply(gt_results, is.null))] #remove the ones with no valid didsetup
    if(length(gt_results) == 0){stop("no valid group-times att to compute")}
    
    gt <- lapply(gt_results, function(x) {x$gt}) |> as.data.table() |> transpose()
    names(gt) <- c("G", "time")
    gt[, time := as.integer(time)]
    
    gt_att <- lapply(gt_results, function(x) {x$result$att})
    gt_inf_func <- lapply(gt_results, function(x) {x$result$inf_func})
    caches <- lapply(gt_results, function(x) {x$result$cache})
    
    gt_names <- gt[,paste0(G,".",time)]
    names(gt_att) <- gt_names
    names(gt_inf_func) <- gt_names
    names(caches) <- gt_names
    
    gt_inf_func <- do.call(cbind, gt_inf_func)
    gt_att <- do.call(cbind, gt_att) |> t()

  return(list(est = list(gt = gt, att = gt_att, inf_func = gt_inf_func), caches = caches))    
}

#gtatt for each outcome, each gt
estimate_gtatt_outcome_gt <- function(gt, y, aux, p, caches){
  
  g <- gt[1]
  t <- as.numeric(gt[2])
  
  #find base time
  gt_name <- paste0(g,".",t)
  base_period <- get_base_period(g,t,p)
  if(t == base_period | #no treatment effect for the base period
     !base_period %in% aux$time_period){ #base period out of bounds
    return(NULL)
  } 
  #find treatment and control group
  did_setup <- get_did_setup(g,t, base_period, aux, p)
  valid_tc_groups <- any(did_setup == 1) & any(did_setup == 0) #if takes up too much time, consider use collapse anyv, but right now quite ok
  if(!isTRUE(valid_tc_groups)){return(NULL)} #no treatment group or control group #isTRUE for na as false
  
  #covariates matrix
  covvars <- get_covvars(base_period, t, aux, p)
  
  #the 2x2 dataset
  cohort_did <- data.table(did_setup, y[[t]], y[[base_period]], aux$weights)
  names(cohort_did) <- c("D", "post.y", "pre.y", "weights")
  
  # estimate --------------------
  result <- tryCatch(estimate_did(dt_did = cohort_did, covvars, p, caches[[gt_name]]),
                     error = function(e){stop("2-by-2 DiD failed for internal group-time ",g , 
                                              "-", t, ": ", e)})
  return(list(gt = gt, result = result))
  
}


get_base_period <- function(g,t,p){
  g <- ming(g) #for two period
  if(p$base_period == "universal"){
    base_period <- g-1-p$anticipation
  } else {
    base_period <- ifelse(t>=g, g-1-p$anticipation, t-1)
  }
  return(base_period)
}

get_did_setup <- function(g, t, base_period, aux, p){

  treated_cohorts <- aux$cohorts[!is.infinite(ming(aux$cohorts))]

  min_control_cohort <- ifelse(p$control_option == "never", Inf, max(t, base_period)+p$anticipation+1)
  max_control_cohort <- ifelse(p$control_option == "notyet", max(ming(treated_cohorts)), Inf) 

  if(!is.na(p$exper$max_control_cohort_diff)){
    max_control_cohort <- min(g+p$exper$max_control_cohort_diff, max_control_cohort)
  } 
  
  #select the control and treated cohorts
  did_setup <- rep(NA, aux$id_size)
  did_setup[get_control_pos(aux$cohort_sizes, min_control_cohort, max_control_cohort)] <- 0
  did_setup[get_treat_pos(aux$cohort_sizes, g)] <- 1 #treated cannot be controls, assign treated after control to overwrite
  
  if(!is.na(p$exper$filtervar)){
    did_setup[!aux$filters[[base_period]]] <- NA #only use units with filter == TRUE at base period
   
  }
  if(!is.na(p$exper$filtervar_post)){
    did_setup[!aux$filters[[t]]] <- NA #only use units with filter == TRUE at target period
  }
  
  if(length(did_setup) != aux$id_size){stop("internal bug: something wrong with did setup (again?)")}
  
  return(did_setup)
}

get_control_pos <- function(cohort_sizes, start_cohort, end_cohort = start_cohort){
  start <- cohort_sizes[ming(G) < start_cohort, sum(cohort_size)]+1 
  end <- cohort_sizes[ming(G) <= end_cohort, sum(cohort_size)]
  if(start > end){return(c())}
  return(seq(start, end, by = 1))
}

get_treat_pos <- function(cohort_sizes, treat_cohort){ #need to separate for double did to match exact g-g-t
  index <- which(cohort_sizes[,G] == treat_cohort)
  start <- ifelse(index == 1, 1, cohort_sizes[1:(index-1), sum(cohort_size)]+1)
  end <- cohort_sizes[1:index, sum(cohort_size)]
  if(start > end){return(c())}
  return(seq(start, end, by = 1))
}

get_covvars <- function(base_period, t, aux, p){
  
  if(all(is.na(p$covariatesvar)) & all(is.na(p$varycovariatesvar))){return(NA)}
  covvars <- data.table()
  
  #add time-varying covariates
  if(!all(is.na(p$varycovariatesvar))){
    
    precov <- aux$varycovariates[[base_period]]
    names(precov) <- paste0("pre_", names(precov))
    postcov <- aux$varycovariates[[t]]-aux$varycovariates[[base_period]]
    names(postcov) <- paste0("post_", names(aux$varycovariates[[t]]))
    covvars <- cbind(covvars, cbind(precov, postcov))
  }
  
  #add time-invariant covariates
  if(!all(is.na(p$covariatesvar))){
    covvars <- cbind(aux$covariates, covvars)
  }
  
  #add constant
  covvars <- as.matrix(cbind(const = -1, covvars))
  return(covvars)
}

#' Fast Staggered DID Estimation
#'
#' Performs Difference-in-Differences (DID) estimation.
#'
#' @param data data.table, the dataset.
#' @param timevar character, name of the time variable.
#' @param cohortvar character, name of the cohort (group) variable.
#' @param unitvar character, name of the unit (id) variable.
#' @param outcomevar character vector, name(s) of the outcome variable(s).
#' @param control_option character, control units used for the DiD estimates, options are "both", "never", or "notyet". 
#' @param result_type character, type of result to return, options are "group_time", "time", "group", "simple", "dynamic" (time since event), "group_group_time", or "dynamic_stagger". 
#' @param balanced_event_time number, max event time to balance the cohort composition.
#' @param control_type character, estimator for controlling for covariates, options are "ipw" (inverse probability weighting), "reg" (outcome regression), or "dr" (doubly-robust).
#' @param allow_unbalance_panel logical, allow unbalance panel as input or coerce dataset into one. 
#' @param boot logical, whether to use bootstrap standard error. 
#' @param biters number, bootstrap iterations. Default is 1000.
#' @param cband logical, whether to use uniform confidence band or point-wise.
#' @param alpha number, the significance level. Default is 0.05.
#' @param weightvar character, name of the weight variable.
#' @param clustervar character, name of the cluster variable. 
#' @param covariatesvar character vector, names of time-invariant covariate variables. 
#' @param varycovariatesvar character vector, names of time-varying covariate variables.
#' @param copy logical, whether to copy the dataset. 
#' @param validate logical, whether to validate the dataset. 
#' @param anticipation number, periods with anticipation. 
#' @param exper list, arguments for experimental features. 
#' @param base_period character, type of base period in pre-preiods, options are "universal", or "varying".
#' @param full logical, whether to return the full result (influence function, call, weighting scheme, etc,.). 
#' @param parallel logical, whether to use parallization on unix system. 
#' @param cohortvar2 character, name of the second cohort (group) variable.
#' @param event_specific logical, whether to recover target treatment effect or use combined effect.
#' @param double_control_option character, control units used for the double DiD, options are "both", "never", or "notyet". 
#' 
#' @import data.table stringr dreamerr ggplot2
#' @importFrom stats quantile vcov sd binomial fitted qnorm rnorm as.formula weighted.mean
#' @importFrom collapse allNA fnrow whichNA fnunique fsum na_insert
#' @importFrom parallel mclapply
#' @importFrom BMisc multiplier_bootstrap
#' @importFrom parglm parglm.fit parglm.control
#' @return A data.table containing the estimated treatment effects and standard errors or a list of all results when `full == TRUE`.
#' @export
#'
#' @details
#' `balanced_event_time` is only meaningful when `result_type == "dynamic`.
#' 
#' `result_type` as `group-group-time` and `dynamic staggered` is only meaningful when using double did.
#' 
#' `biter` and `clustervar` is only used when `boot == TRUE`.
#'
#' @examples
#' # simulated data
#' simdt <- sim_did(1e+02, 10, cov = "cont", second_cov = TRUE, second_outcome = TRUE, seed = 1)
#' dt <- simdt$dt
#' 
#' #basic call
#' result <- fastdid(data = dt, timevar = "time", cohortvar = "G", 
#'                   unitvar = "unit", outcomevar = "y",  
#'                   result_type = "group_time")
#'
#' @keywords difference-in-differences fast computation panel data estimation did
fastdid <- function(data,
                    timevar, cohortvar, unitvar, outcomevar, 
                    control_option="both",result_type="group_time", balanced_event_time = NA,
                    control_type = "ipw", allow_unbalance_panel = FALSE, boot=FALSE, biters = 1000, cband = FALSE, alpha = 0.05,
                    weightvar=NA,clustervar=NA, covariatesvar = NA, varycovariatesvar = NA, 
                    copy = TRUE, validate = TRUE,
                    anticipation = 0,  base_period = "universal",
                    exper = NULL, full = FALSE, parallel = FALSE, 
                    cohortvar2 = NA, event_specific = TRUE, double_control_option="both"){
  
  # preprocess --------------------------------------------------------
  
  if(!is.data.table(data)){
    warning("coercing input into a data.table.")
    data <- as.data.table(data)
  } 
  if(copy){dt <- copy(data)} else {dt <- data}
  
  # validate arguments
  p <- as.list(environment()) #collect everything besides data
  p$data <- NULL
  p$dt <- NULL
  
  exper_args <- c("filtervar", "filtervar_post", "only_balance_2by2",
                  "aggregate_scheme", "max_control_cohort_diff")
  p$exper <- get_exper_default(p$exper, exper_args)
  class(p) <- "locked" #no more changes!
  validate_argument(dt, p)
  
  #change name for main columns
  setnames(dt, c(timevar, cohortvar, unitvar), c("time", "G", "unit"))
  if(!is.na(p$cohortvar2)){setnames(dt, p$cohortvar2, "G2")}
  
  # validate and throw away not legal data 
  dt <- validate_dt(dt, p)
  
  #make dt conform to the WLOG assumptions of fastdid
  coerce_result <- coerce_dt(dt, p) #also changed dt
  
  # get auxiliary data
  aux <- get_auxdata(coerce_result$dt, p)

  # main estimation  -------------------------------------------------

  gt_result_list <- estimate_gtatt(aux, p)
  agg_result <- aggregate_gt(gt_result_list, aux, p)
  
  # post process -------------------------------------------
  
  #convert "targets" back to meaningful parameters
  est_results <- convert_targets(agg_result$est, p, coerce_result$t) 
  
  if(!p$full){
    return(est_results)
  } else {
    full_result <- list(
      call = p,
      estimate = est_results,
      gt_estimate = gt_result_list,
      agg_inf_func = agg_result$inf_func,
      agg_weight_matrix = agg_result$agg_weight_matrix
    )
    class(full_result) <- c("fastdid_result", class(full_result))
    return(full_result)
  }
}

#' Plot event study
#'
#' Plot event study results.
#'
#' @param x A data table generated with [fastdid] with one-dimensional index.
#' @param margin character, the x-axis of the plot
#'
#' @return A ggplot2 object
#' @examples
#' 
#' # simulated data
#' simdt <- sim_did(1e+02, 10, seed = 1)
#' dt <- simdt$dt
#' 
#' #estimation
#' result <- fastdid(data = dt, timevar = "time", cohortvar = "G", 
#'                   unitvar = "unit", outcomevar = "y",  
#'                   result_type = "dynamic")
#' 
#' #plot
#' plot_did_dynamics(result)
#' 
#' @export
plot_did_dynamics <-function(x, margin = "event_time"){
  
  #find the base_period
  if(margin == "event_time"){
    et_range <- min(x[, event_time]):max(x[, event_time])
    base_time <- et_range[!et_range %in% x[, unique(event_time)]]
    if(length(base_time)!=1){stop("missing more than one period")}
    
    #add the base period
    if("outcome" %in% names(x)){
      base_row <- data.table(att = 0, se = 0, event_time = base_time, outcome = x[, unique(outcome)], att_ciub = 0, att_cilb = 0)
    } else {
      base_row <- data.table(att = 0, se = 0, event_time = base_time, att_ciub = 0, att_cilb = 0)
    }
    x <- x |> rbind(base_row, fill = TRUE)
  } else {
    x <- x[type == "post"]
  }
  
  plot <- x |> 
     ggplot() +
     geom_hline(yintercept = 0, linetype = 'dashed', col = 'red') + 
     geom_point( aes(x = eval(str2lang(margin)), y = att), color = "black") + #point est
     geom_errorbar( aes(x = eval(str2lang(margin)), ymin = att_cilb, ymax = att_ciub), 
                width = 0.1, linetype = "dashed") + #CI
     labs(x = margin)
  
  if(margin == "event_time"){
    plot <- plot + geom_line( aes(x = eval(str2lang(margin)), y = att), color = "black") #point est
  }

  return(plot)
  
}

NULL

# quiets concerns of R CMD check re: the .'s that appear in pipelines and data.table variables
utils::globalVariables(c('.','agg_weight','att','att_cont','att_treat','attgt','cohort','cohort_size','conf_lwb','conf_upb',
                         'const','cont_ipw_weight','count','delta_y','element_rect','element_text','event_time','pg','placeholder',
                         'post.y','pre.y','ps','s','se','target','tau','time_fe',
                         'treat_ipw_weight','treat_latent','type','unit','unit_fe','weight','x','x2',
                         'x_trend','y','y0','y1','y2', 'time', 'weights', 'outcome', "G", "D", 'xvar',
                         'V1','att_cont_post','att_cont_pre','att_treat_post','att_treat_pre','inpost','inpre','max_et','min_et','new_unit','or_delta','or_delta_post','or_delta_pre','targeted','used',
                         "timevar", "cohortvar", "unitvar", "outcomevar", "control_option", "result_type", "balanced_event_time", "control_type",
                         "allow_unbalance_panel", "boot", "biters", "weightvar", "clustervar", "covariatesvar", "varycovariatesvar", "filtervar",
                         "copy", "validate", "max_control_cohort_diff", "anticipation", "min_control_cohort_diff", "base_period", "post", "att_ciub", "att_cilb", "cband", "alpha",
                         "G2", "G1", "mg", "cohort1", "cohort2", "event_time_1", "event_time_2",
                         "D2", "attgt2", "event", "atu2", "y01", "y10", "y11", "tau2", "parallel",
                         "tp", "cp", "tb", "cb", "no_na", "event_stagger", "double_control_option",
                         "det_weight", "sto_weight"))


#' Simulate a Difference-in-Differences (DiD) dataset
#'
#' Simulates a dataset for a Difference-in-Differences analysis with various customizable options.
#'
#' @param sample_size The number of units in the dataset.
#' @param time_period The number of time periods in the dataset.
#' @param untreated_prop The proportion of untreated units.
#' @param epsilon_size The standard deviation for the error term in potential outcomes.
#' @param cov The type of covariate to include ("no", "int", or "cont").
#' @param hetero The type of heterogeneity in treatment effects ("all" or "dynamic").
#' @param second_outcome Whether to include a second outcome variable.
#' @param second_cov Whether to include a second covariate.
#' @param na Whether to generate missing data ("none", "y", "x", or "both").
#' @param balanced Whether to balance the dataset by random sampling.
#' @param seed Seed for random number generation.
#' @param stratify Whether to stratify the dataset based on a binary covariate.
#' @param treatment_assign The method for treatment assignment ("latent" or "uniform").
#' @param vary_cov include time-varying covariates
#' @param second_cohort include confounding events
#' @param confound_ratio extent of event confoundedness
#' @param second_het heterogeneity of the second event
#'
#' @return A list containing the simulated dataset (dt) and the treatment effect values (att).
#'
#' @examples
#' # Simulate a DiD dataset with default settings
#' data <- sim_did(sample_size = 100, time_period = 5)
#'
#' @export
sim_did <- function(sample_size, time_period, untreated_prop = 0.3, epsilon_size = 0.001,
                    cov = "no", hetero = "all", second_outcome = FALSE, second_cov = FALSE, vary_cov = FALSE, na = "none", 
                    balanced = TRUE, seed = NA, stratify = FALSE, treatment_assign = "latent", second_cohort = FALSE, confound_ratio = 1, second_het = "all"){
  
  if(!is.na(seed)){set.seed(seed)}
  
  #unit  -------------
  dt_i <- data.table(unit = 1:sample_size)
  if(cov == "int"){
    dt_i[, x := sample.int(5, sample_size, replace = TRUE)] 
  } else if (cov == "no"){
    dt_i[, x := 1] 
  } else if (cov == "cont"){
    dt_i[, x := rnorm(sample_size)]
  }
  
  if(second_cov){
    dt_i[, x2 := rnorm(sample_size)]
  } else {dt_i[, x2 := 1]}
  
  if(stratify){
    dt_i[, s := fifelse(rnorm(sample_size) > 0, 1, 2)]
  } else {
    dt_i[, s := 1]
  }

  #treatment assignment ---------------------
  
  #assign treated group based on a latent related to X
  
  if(treatment_assign == "latent"){
    ep1 <- rnorm(sample_size)
    dt_i[, treat_latent := x*0.2 + x2*0.2 + ep1] #unit with larger X tend to be treated and treated earlier
    untreated_thres <- quantile(dt_i$treat_latent, untreated_prop)
    dt_i[treat_latent <= untreated_thres, G := Inf] #unit with low latent is never treated
    
    cohort_prop <- (1-untreated_prop)/(time_period-1)
    last_treat_thres <- untreated_thres
    for(t in time_period:2){ #unit start getting treated in t = 2
      treat_thres <- quantile(dt_i$treat_latent, untreated_prop + cohort_prop*(time_period - t + 1))
      dt_i[treat_latent <= treat_thres & treat_latent > last_treat_thres, G := t]
      last_treat_thres <- treat_thres
    }
    rm(t)
  } else if (treatment_assign == "uniform"){
    #when treatment is set to 'uniform', untreated propensity is fixed
    dt_i[,G := floor((unit-1)/(sample_size/time_period))]
    dt_i[G < 2, G := Inf]
  }
  
  if(second_cohort){
    setnames(dt_i, "G", "G2")
    if(treatment_assign == "latent"){
      dt_i[, treat_latent := x*0.2 + x2*0.2 + ep1*confound_ratio + rnorm(sample_size)] #unit with larger X tend to be treated and treated earlier
      untreated_thres <- quantile(dt_i$treat_latent, untreated_prop)
      dt_i[treat_latent <= untreated_thres, G := Inf] #unit with low latent is never treated
      
      cohort_prop <- (1-untreated_prop)/(time_period-1)
      last_treat_thres <- untreated_thres
      for(t in time_period:2){ #unit start getting treated in t = 2
        treat_thres <- quantile(dt_i$treat_latent, untreated_prop + cohort_prop*(time_period - t + 1))
        dt_i[treat_latent <= treat_thres & treat_latent > last_treat_thres, G := t]
        last_treat_thres <- treat_thres
      }
      rm(t)
    } else if (treatment_assign == "uniform"){
      #when treatment is set to 'uniform', untreated propensity is fixed
      dt_i[,G := floor((unit-1)/(sample_size/time_period))]
      dt_i[G < 2, G := Inf]
    }
  }

  #assign unit FE
  dt_i[, unit_fe := rnorm(sample_size)]
  
  #time ------------------
  
  dt_t <- data.table(time = 1:time_period)
  dt_t[, time_fe := rnorm(time_period)]
  dt_t[, x_trend := rnorm(time_period)]
  
  #panel --------------------------
  dt <- CJ(unit = 1:sample_size, time = 1:time_period)
  dt <- dt |> merge(dt_i, by = "unit")
  dt <- dt |> merge(dt_t, by = "time")
  
  #add time_varying covariates
  if(vary_cov){
    dt[, xvar := pmin(G, time_period+4)*time^(1/3)*0.1+rnorm(sample_size*time_period, 0,10)] #should be confounding....?
  } else {
    dt[, xvar := 1]
  }
  
  #untreated potential outcomes
  dt[, y0 := unit_fe + time_fe + (x+x2)*x_trend + xvar + rnorm(sample_size*time_period, sd = epsilon_size)]
  
  dt[, D := as.integer(time >= G)]
  #generate gtatt
  att <- CJ(G = 1:time_period, time = 1:time_period)
  if(hetero == "all"){
    att[, attgt := rnorm(time_period*time_period, mean = 2, sd = 1)]
  } else if (hetero == "dynamic"){
    for(event_t in 0:max(att[,time-G])){
      att[time - G == event_t, attgt := rnorm(1, mean = 2, sd = 1)]
    }
  }
  att[time < G, attgt := 0] #no anticipation
  dt <- dt |> merge(att, by = c("time", "G"), all.x = TRUE, all.y = FALSE)
  dt[is.na(attgt), attgt := 0]
  dt[, tau := attgt*s]
  
  if(second_cohort){
    dt[, D2 := as.integer(time >= G2)]
    #generate gtatt
    att2 <- CJ(G2 = 1:time_period, time = 1:time_period)
    if(second_het == "no"){
      att2[, attgt2 := 10]
    } else {
      if(hetero == "all"){
        att2[, attgt2 := rnorm(time_period*time_period, mean = 2, sd = 1)]
      } else if (hetero == "dynamic"){
        for(event_t in 0:max(att2[,time-G2])){
          att2[time - G2 == event_t, attgt2 := rnorm(1, mean = 2, sd = 1)]
        }
      }
    }
    att2[time < G2, attgt2 := 0] #no anticipation
    
    #add att2 to att
    att[, event := 1]
    att2[, event := 2]
    att <- rbind(att, att2, fill = TRUE)
    
    dt <- dt |> merge(att2, by = c("time", "G2"), all.x = TRUE, all.y = FALSE)
    dt[is.na(attgt2), attgt2 := 0]
    dt[, tau2 := attgt2*s]
    
    #potential outcome
    dt[, y10 := y0 + tau]
    dt[, y01 := y0 + tau2]
    dt[, y11 := y0 + tau + tau2]
    dt[, y := y0*(1-D)*(1-D2)+ y10*D*(1-D2) + y01*(1-D)*D2 + y11*D*D2]
    cols <- c("time", "G", "G2", "unit", "x", "x2", "y", "s", "xvar")
    
  } else {
    #potential outcome
    dt[, y1 := y0 + tau]
    dt[, y := y1*D + y0*(1-D)]
    cols <- c("time", "G", "unit", "x", "x2", "y", "s", "xvar")
  }
  
  dt <- dt[, .SD, .SDcols = cols]
  
  #additional -----------------
  
  if(na == "y"){
    dt[, y := na_insert(y)]
  } else if (na == "x") {
    dt[, x := na_insert(x)]
  } else if( na == "both"){
    dt[, y := na_insert(y)]
    dt[, x := na_insert(x)]
  }
  
  if(balanced == FALSE){
    size <- fnrow(dt)
    dt <- dt[sample(1:size, size*0.99)]
  }
  if(second_outcome == TRUE){  dt[, y2 := y + 1 + rnorm(fnrow(dt), 0, 0.1)]}
  
  return(list(dt = dt, att = att))
  
}





validate_argument <- function(dt, p){
  
  if(!p$validate){
    return(NULL)
  }
  
  dt_names <- names(dt)
  
  #release p
  for(name in names(p)){
    assign(name, p[[name]])
  }

  name_message <- "__ARG__ must be a character scalar and a name of a column from the dataset."
  check_set_arg(timevar, unitvar, cohortvar, "match", .choices = dt_names, .message = name_message, .up = 1)
  
  covariate_message <- "__ARG__ must be NA or a character vector which are all names of columns from the dataset."
  check_set_arg(varycovariatesvar, covariatesvar, outcomevar, 
                "NA | multi match", .choices = dt_names, .message = covariate_message, .up = 1)
  
  checkvar_message <- "__ARG__ must be NA or a character scalar if a name of columns from the dataset."
  check_set_arg(weightvar, clustervar, "NA | match", .choices = dt_names, .message = checkvar_message, .up = 1)
  
  check_set_arg(control_option, double_control_option, "match", .choices = c("both", "never", "notyet"), .up = 1) #kinda bad names since did's notyet include both notyet and never
  check_set_arg(control_type, "match", .choices = c("ipw", "reg", "dr"), .up = 1) 
  check_set_arg(base_period, "match", .choices = c("varying", "universal"), .up = 1)
  check_arg(copy, validate, boot, allow_unbalance_panel, cband, parallel, "scalar logical", .up = 1)
  check_arg(anticipation, alpha, "scalar numeric", .up = 1)
  
  if(!is.na(balanced_event_time)){
    if(result_type != "dynamic"){stop("balanced_event_time is only meaningful with result_type == 'dynamic'")}
    check_arg(balanced_event_time, "numeric scalar", .up = 1)
  }
  if(allow_unbalance_panel == TRUE & control_type == "dr"){
    stop("fastdid does not support DR when allowing for unbalanced panels.")
  }
  # if(allow_unbalance_panel == TRUE & !allNA(varycovariatesvar)){
  #   stop("fastdid currently only supprts time varying covariates when not allowing for unbalanced panels.")
  # }
  if(any(covariatesvar %in% varycovariatesvar) & !allNA(varycovariatesvar) & !allNA(covariatesvar)){
    stop("time-varying var and invariant var have overlaps.")
  }
  if(!boot & (!allNA(clustervar)|cband == TRUE)){
    stop("clustering and uniform confidence interval only available with bootstrap")
  }
  
  if(parallel){
    if(.Platform$OS.type != "unix"){stop("parallel option only available on unix sysytems")}
    if(!requireNamespace("parallel")){stop("parallel requires the parallel package")}
  }
  
  # varname collision
  varnames <- unlist(p[str_subset(names(p), "var")])
  varnames <- varnames[!is.na(varnames)]
  if(any(duplicated(varnames))){stop("-var arguments can not have duplicated names. (no need to specicify cluster on unit-level, it is automatically done.)")}
}

validate_dt <- function(dt, p){

  varnames <- unlist(p[str_ends(names(p), "var")], recursive = TRUE) #get all the argument that ends with "var"
  varnames <- varnames[!varnames %in% c(p$timevar, p$unitvar, p$cohortvar) & !is.na(varnames) & !is.null(varnames)]
  
  #change to int 
  uniquecols <- c("G", "time", "unit")
  for(col in uniquecols){
      if(!dt[, is.numeric(get(col))]){stop(col, " needs to be numeric.")}
      dt[!is.infinite(get(col)), c(col) := as.integer(get(col))] #yeah sometimes floating point can be annoying
  }
  
  raw_unit_size <- dt[, uniqueN(unit)]
  raw_time_size <- dt[, uniqueN(time)]
  
  if(!is.na(p$balanced_event_time)){
    if(p$balanced_event_time > dt[, max(time-G)]){stop("balanced_event_time is larger than the max event time in the data")}
  }

  #doesn't allow missing value
  if(is.na(p$exper$only_balance_2by2) | !p$exper$only_balance_2by2){
    for(col in varnames){
      na_obs <- whichNA(dt[,get(col)])
      if(length(na_obs) != 0){
        warning("missing values detected in ", col, ", removing ", length(na_obs), " observation.")
        dt <- dt[!na_obs]
      }
    }
  }
  
  if(!allNA(p$covariatesvar) && uniqueN(dt, by = c("unit", p$covariatesvar)) > raw_unit_size){
    warning("some covariates is time-varying, fastdid only use the first observation for covariates.")
  }
  
  
  if(!allNA(p$covariatesvar)|!allNA(p$varycovariatesvar)){
    for(cov in c(p$covariatesvar, p$varycovariatesvar)){
      if(is.na(cov)){next}
      #check covaraites is not constant  
      if(fnunique(dt[, get(cov)[1], by = "unit"][, V1]) == 1)stop(cov, " have no variation")
      if(!(is.numeric(dt[, get(cov)])|is.integer(dt[, get(cov)]))){stop(cov, " is not numeric or integer, do not support fixed effects.")}
    }
  }
  
  #check balanced panel
  #check if any is dup
  if(anyDuplicated(dt[, .(unit, time)])){
    dup_id <- dt[duplicated(dt[,.(unit, time)]), unique(unit)]
    stop(length(dup_id), " units is observed more than once in a period.")
  }
  
  #check if any is missing
  if(!p$allow_unbalance_panel){
    unit_count <- dt[, .(count = .N), by = unit]
    if(any(unit_count[, count < raw_time_size])){
      mis_unit <- unit_count[count < raw_time_size]
      warning(nrow(mis_unit), " units is missing in some periods, enforcing balanced panel by dropping them")
      dt <- dt[!unit %in% mis_unit[, unit]]
    }
  }
  
  # drop always_treated units
  always_treated <- dt[G <= min(time), unique(unit)]
  if(length(always_treated) > 0){
    warning(length(always_treated), " units is treated in the first period, dropping them")
    dt <- dt[!unit %in% always_treated]
  }
  
  #for double did part
  if(!is.na(p$cohortvar2)){
    always_treated <- dt[G2 <= min(time), unique(unit)]
    if(length(always_treated) > 0){
      warning(length(always_treated), " units is treated in the first period, dropping them")
      dt <- dt[!unit %in% always_treated]
    }
  }

  
  
  return(dt)

}

