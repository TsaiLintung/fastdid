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
    for(col in varnames){
      if(is.na(col)){next}
      dt[is.na(get(col)), no_na := FALSE]
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