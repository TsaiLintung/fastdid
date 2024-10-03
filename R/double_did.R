

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
  es_weight <- do.call(rbind, es_weight_list) #not sure if the dim is right
  
  return(list(group_time = es_group_time, es_weight = es_weight))
  
}

#get the scheme for retriving group-group-time estimates
get_es_ggt_weight <- function(ggt, group_time, aux, p){
  
  group_time <- copy(group_time) #avoid accidental modification
  
  group_time[, weight := 0] #reset 
  t <- group_time[ggt, time]
  g1 <- group_time[ggt, G1]
  g2 <- group_time[ggt, G2]
  gg <- group_time[ggt, G]
  
  if(is.infinite(g1)){return(NULL)}
  
  if(t < g2){ #direct pure effect
    
    group_time[ggt, weight := 1] #just use the observed effect
    
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
    group_time[tb, weight := pg/sum(pg)]
    group_time[cp, weight := pg/sum(pg)]
    group_time[cb, weight := -pg/sum(pg)]
    
    
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
    group_time[tp, weight := pg/sum(pg)]
    group_time[tb, weight := -pg/sum(pg)]
    group_time[cp, weight := -pg/sum(pg)]
    group_time[cb, weight := pg/sum(pg)]
    
  } 
  
  if(all(group_time[, weight] == 0)){return(NULL)}
  return(group_time[, weight])
  
}
