#2026-09-15
message('loading fastdid source ver. ver: 1.1.0 date: 2026-09-15')
require(data.table);
 require(stringr);
 require(BMisc);
 require(collapse);
 require(dreamerr);
# high level -------------------------------------------------------------------

aggregate_gt <- function(all_gt_result, aux, p) {
  results <- lapply(all_gt_result, aggregate_gt_outcome, aux, p)
  return(list(
    est = rbindlist(lapply(results, function(x) {
      x$result
    })),
    inf_func = lapply(results, function(x) {
      x$inf_func
    }),
    agg_weight_matrix = lapply(results, function(x) {
      x$weight_matrix
    }),
    es_weight_matrix = lapply(results, function(x) {
      x$es_weight
    }),
    effect_diag = lapply(results, function(x) {
      x$effect_diag
    })
  ))
}

aggregate_gt_outcome <- function(gt_result, aux, p) {
  cells <- get_cell_table(gt_result, aux, p)
  att <- gt_result$att
  inf_func <- gt_result$inf_func

  # the second stage runs before the aggregation, on the first-stage cells
  es_weight <- NULL
  effect_diag <- NULL
  if (p$event_specific && !allNA(p$cohortvar2)) {
    ss <- second_stage(cells, att, inf_func, aux, p)
    cells <- ss$cells # some gt may not have an identified effect (ex: g1 == g2)
    att <- ss$att
    inf_func <- ss$inf_func
    es_weight <- ss$weight
    effect_diag <- ss$diag
  }

  # get aggregation scheme from cells to target parameters
  agg_sch <- get_agg_sch(cells, p)

  # get att
  agg_att <- agg_sch$agg_weights %*% att

  # get influence function matrix

  inf_weights <- get_weight_influence(att, agg_sch$group_time, agg_sch$agg_weights, aux, p)
  inf_matrix <- (inf_func %*% t(agg_sch$agg_weights)) + inf_weights

  # get se
  agg_se <- get_se(inf_matrix, aux, p)

  # post process
  result <- data.table(agg_sch$targets, agg_att, agg_se$se)
  names(result) <- c("target", "att", "se")
  result[, `:=`(
    outcome = gt_result$outname,
    att_ciub = att + se * agg_se$crit_val,
    att_cilb = att - se * agg_se$crit_val
  )]

  return(list(
    result = result,
    inf_func = inf_matrix,
    weight_matrix = agg_sch$agg_weights,
    es_weight = es_weight,
    effect_diag = effect_diag
  ))
}

# scheme ------------------------------------------------------------------------

#' Scheme for aggregation.
#'
#' @param group_time the cell table, after the second stage when there is one.
#' @return the targets, the weight of each cell in each target, and the table.
#' @noRd
get_agg_sch <- function(group_time, p) {
  # choose the target based on aggregation type
  tg <- get_agg_targets(group_time, p)
  group_time <- tg$group_time
  targets <- tg$targets

  # get aggregation weights
  agg_weights <- data.table()
  for (tar in targets) { # the order matters

    group_time[, weight := 0] # weight is 0 if not a target
    group_time[target == tar & used, weight := pg / sum(pg)]
    target_weights <- group_time[, .(weight)] |> transpose()
    agg_weights <- rbind(agg_weights, target_weights)
  }
  group_time[, pg := NULL]

  agg_weights <- as.matrix(agg_weights)
  return(list(
    agg_weights = agg_weights, # a matrix of each target and gt's weight in it
    targets = targets,
    group_time = group_time
  ))
}

#' Get the target parameters.
#' @noRd
get_agg_targets <- function(group_time, p) {
  group_time[, post := as.numeric(ifelse(time >= g1(G), 1, -1))]
  switch(p$result_type,
    dynamic = group_time[, target := time - g1(G)],
    group = group_time[, target := g1(G) * post], # group * treated
    time = group_time[, target := time * post],
    simple = group_time[, target := post],
    group_time = group_time[, target := paste0(g1(G), ".", time)],
    group_group_time = group_time[, target := paste0(G, ".", time)],
    dynamic_stagger = group_time[, target := paste0(time - g1(G), ".", g1(G) - gprime(G))],
    dynamic_event = group_time[, target := e] # the event time of the modeled event
  )

  # allow custom aggregation scheme, this overides other stuff
  if (!is.na(p$exper$aggregate_scheme)) {
    group_time[, target := eval(str2lang(p$exper$aggregate_scheme))]
  }

  targets <- group_time[, unique(target)]

  # for balanced cohort composition in dynamic setting
  # a cohort us only used if it is seen for all dynamic time
  if (p$result_type == "dynamic" && !is.na(p$balanced_event_time)) {
    cohorts <- group_time[, .(
      max_et = max(target), # event time is target if in dynamic
      min_et = min(target)
    ), by = "G"]
    cohorts[, used := max_et >= p$balanced_event_time] # the max
    if (!cohorts[, any(used)]) {
      stop("balanced_comp_range outside available range")
    }
    group_time[, used := G %in% cohorts[used == TRUE, G]]

    min_event_time <- cohorts[used == TRUE, min(min_et)]
    max_event_time <- p$balanced_event_time
    
    if (min_event_time > max_event_time) {
      stop("Invalid balanced_event_time: The minimum available event time (", min_event_time, 
           ") is greater than balanced_event_time (", max_event_time, "). ",
           "Please specify a balanced_event_time >= ", min_event_time, 
           " or use a smaller value that matches your data structure.")
    }
    
    targets <- targets[targets <= p$balanced_event_time & targets >= min_event_time]
  } else {
    group_time[, used := TRUE]
  }

  return(list(group_time = group_time, targets = targets))
}

# influence function ------------------------------------------------------------

get_weight_influence <- function(att, group, agg_weights, aux, p, by_period = FALSE) {
  id_dt <- data.table(weight = aux$weights / sum(aux$weights), G = aux$dt_inv[, G])
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  group <- group |> merge(pg_dt, by = "G", sort = FALSE)

  group[, time := as.integer(time)]

  if (allNA(p$cohortvar2)) {
    group[, G := as.integer(G)]
    setorder(group, time, G)
  } else {
    M <- 1L + length(p$cohortvar2)
    gcol_w <- paste0("G", seq_len(M))
    group[, mg := ming(G)]
    for(d in seq_len(M)){
      group[, (paste0("G", d)) := gd(G, d)]
    }
    sortcols <- c("time", "mg", gcol_w)
    if ("event" %in% names(group)) sortcols <- c(sortcols, "event") # one row per event in a stacked fit
    do.call(setorderv, c(list(group), list(sortcols))) # sort
  }

  if (!p$parallel) {
    inf_weights <- sapply(asplit(agg_weights, 1), function(x) {
      get_weight_influence_param(x, group, att, aux, p, by_period)
    })
  } else {
    inf_weights <- matrix(unlist(mclapply(asplit(agg_weights, 1), function(x) {
      get_weight_influence_param(x, group, att, aux, p, by_period)
    })), ncol = dim(agg_weights)[1])
  }

  return(inf_weights)
}

#' Influence from the weight calculation.
#'
#' The weight of a cell is the cohort share pgi / sum(pgi), so the estimated
#' weight adds a term to the influence function. The plain aggregation weights
#' are positive and share one denominator. The double DiD weights are a signed
#' pair, and each period has its own denominator, so `by_period` splits the sum
#' into one block for each period and keeps the sign.
#'
#' @param agg_weights numeric vector, the weight of each row of `group`.
#' @param group the group-time table, with the cohort share `pg`.
#' @param gt_att the g-t estimates.
#' @param by_period logical, normalize each period on its own.
#' @return a column of the influence function for the weights.
#' @noRd
get_weight_influence_param <- function(agg_weights, group, gt_att, aux, p, by_period = FALSE) {
  keepers <- which(agg_weights != 0)
  if (length(keepers) == 0) {
    return(rep(0, length(aux$weights)))
  } # for direct double did
  group <- group[keepers, ]
  signs <- sign(agg_weights[keepers])
  att_keep <- as.vector(gt_att)[keepers]

  # moving this outside will create a g*t*id matrix, not really worth the memory
  keepers_matrix <- as.matrix(aux$weights * sapply(seq_len(nrow(group)), function(g) {
    as.integer(aux$dt_inv[, G] == group[g, G]) - group[g, pg]
  }))

  blocks <- if (by_period) group[, time] else rep(1L, nrow(group))

  # one normalized share for each block: d(pgi / sum(pgi))
  inf_weight <- rep(0, length(aux$weights))
  for (b in unique(blocks)) {
    idx <- which(blocks == b)
    pg_b <- group[idx, pg]
    block_matrix <- keepers_matrix[, idx, drop = FALSE]
    if1 <- block_matrix / sum(pg_b) # numerator
    if2 <- rowSums(block_matrix) %*% t(pg_b) / (sum(pg_b)^2) # denominator
    inf_weight <- inf_weight + (if1 - if2) %*% (signs[idx] * att_keep[idx])
  }

  inf_weight[abs(inf_weight) < sqrt(.Machine$double.eps) * 10] <- 0 # fill zero
  return(inf_weight)
}

# se -------------------------------------------------------------------

#' Aggregated standard error.
#' @noRd
get_se <- function(inf_matrix, aux, p) {
  if (p$boot) {
    cluster <- aux$cluster

    top_quant <- 0.75
    bot_quant <- 0.25
    if (!allNA(p$clustervar)) {
      # take average within the cluster
      cluster_n <- stats::aggregate(cluster, by = list(cluster), length)[, 2]
      inf_matrix <- fsum(inf_matrix, cluster) / cluster_n # the mean without 0 for each cluster of each setting
    }

    boot_results <- BMisc::multiplier_bootstrap(inf_matrix, biters = p$biters) |> as.data.table()

    boot_top <- boot_results[, lapply(.SD, function(x) stats::quantile(x, top_quant, type = 1, na.rm = TRUE))]
    boot_bot <- boot_results[, lapply(.SD, function(x) stats::quantile(x, bot_quant, type = 1, na.rm = TRUE))]

    dt_se <- rbind(boot_bot, boot_top) |> transpose()
    names(dt_se) <- c("boot_bot", "boot_top")

    # get sigma
    se <- dt_se[, (boot_top - boot_bot) / (qnorm(top_quant) - qnorm(bot_quant))]
    se[se < sqrt(.Machine$double.eps) * 10] <- NA
  } else {
    inf_matrix <- inf_matrix |> as.data.table()
    se <- inf_matrix[, lapply(.SD, function(x) sqrt(sum(x^2, na.rm = TRUE) / length(x)^2))] |> as.vector() # divides by n, matching the did package (see TODO.md)
  }

  # get critical value
  crit_val <- NA
  point_crit_val <- qnorm(1 - p$alpha / 2)
  if (p$cband) {
    boot_tv <- apply(boot_results, 1, function(b) {
      max(abs(b / se), na.rm = TRUE)
    })
    boot_tv <- boot_tv[is.finite(boot_tv)]
    crit_val <- quantile(boot_tv, 1 - p$alpha, type = 1, na.rm = TRUE)
  }
  if (is.na(crit_val) || is.infinite(crit_val) || crit_val < point_crit_val) {
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
  
  if(!is.na(exper$only_balance_2by2) && exper$only_balance_2by2){ #will create this col in the get_aux part
    exper$filtervar <- "no_na"
    exper$filtervar_post <- "no_na"
  }
  if(is.na(exper$effect_tol)){
    exper$effect_tol <- 1e-7
  }
  
  return(exper)
}

coerce_dt <- function(dt, p){
  
  if(!allNA(p$cohortvar2)){return(coerce_dt_doub(dt, p))} #in doubledid.R

  if(nrow(dt) == 0){
    stop("no data after coercing the dataset")
  }

  #check if there is available never-treated group
  if(!is.infinite(dt[, max(G)])){
    if(p$control_option == "both"){warning("no never-treated available, effectively using not-yet-but-eventually-treated as control")}
    if(p$control_option == "never"){stop("no never-treated available.")}
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

  time_offset <- min(time_periods) - 1 #assume time starts at 1, first is min after sort :)
  if(time_offset != 0){
    dt[, G := G-time_offset]
    
    dt[, time := time-time_offset]
    time_periods <- time_periods - time_offset
  }
  
  time_step <- 1 #time may not jump at 1
  if(length(time_periods) > 1 && any(time_periods[2:length(time_periods)] - time_periods[seq_len(length(time_periods)-1)] != 1)){
    # Calculate all intervals between consecutive time periods
    all_intervals <- time_periods[2:length(time_periods)] - time_periods[seq_len(length(time_periods)-1)]

    
    
    # Check if all intervals are identical (uniform step)
    if(length(unique(all_intervals)) > 1){
      stop("Time step is not uniform. Time periods: ", paste(head(time_periods, 10), collapse = ", "),
           if(length(time_periods) > 10) "..." else "",
           ". Intervals between periods: ", paste(unique(all_intervals), collapse = ", "),
           ". fastdid requires uniformly-spaced time periods.")
    }
    
    time_step <- all_intervals[1]
    time_periods <- (time_periods-1)/time_step+1
    
    # Verify that normalization worked (should always be consecutive integers now)
    if(any(time_periods[2:length(time_periods)] - time_periods[seq_len(length(time_periods)-1)] != 1)){
      stop("Internal error: time normalization failed. Please report this issue.")
    }
    
    dt[G != 1, G := (G-1)/time_step+1]
    dt[time != 1, time := (time-1)/time_step+1]
  }
  
  #add the information to t
  t <- list()
  t$time_step <- time_step
  t$time_offset <- time_offset

  return(list(dt = dt, p = p, t = t))
  
}

get_auxdata <- function(dt, p){
  
  time_periods <- dt[, unique(time)]
  id_size <- dt[, uniqueN(unit)]
  
  # Validate basic data structure
  if(id_size <= 0){
    stop("Invalid data structure: id_size is ", id_size, ". Dataset must contain at least one unit.")
  }
  if(length(time_periods) == 0){
    stop("Invalid data structure: no time periods found in the dataset.")
  }
  
  #construct the outcomes list for fast access later
  #loop for multiple outcome
  outcomes_list <- list()
  for(outcol in p$outcomevar){
    outcomes <- list()
    
    if(!p$allow_unbalance_panel){
      for(i in time_periods){
        start <- (i-1)*id_size+1
        end <- i*id_size
        if(start > end){
          stop("Invalid sequence when extracting outcomes: start (", start, ") > end (", end, ") for time period ", i, 
               ". This suggests an issue with the data structure (id_size=", id_size, ").")
        }
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
      if(start > end){
        stop("Invalid sequence when extracting varying covariates: start (", start, ") > end (", end, ") for time period ", i, 
             ". This suggests an issue with the data structure (id_size=", id_size, ").")
      }
      varycovariates[[i]] <- dt[seq(start,end, by = 1), .SD, .SDcols = p$varycovariatesvar]
    }
  } else {
    varycovariates <- NA
  }
  
  #create na indicator for filtering
  if(!is.na(p$exper$only_balance_2by2) && p$exper$only_balance_2by2){
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

           M <- 1L + length(p$cohortvar2)
           for(d in seq_len(M)){
             col_name <- paste0("cohort", d)
             results[, (col_name) := recover_time(gd(cohort, d), t)]
           }

           results[, time := recover_time(time, t)]
           results[, `:=`(cohort = NULL)]
         },
         dynamic_stagger = {
           results[, event_time_1 :=  as.numeric(str_split_i(target, "\\.", 1))]
           results[, event_stagger :=  as.numeric(str_split_i(target, "\\.", 2))]
         },
         dynamic_event = {
           results[, event_time := target]
           setcolorder(results, "event_time", before = 1)
         }
  )
  
  results[, target := NULL]
  
  if(p$add_base_period){
    results <- results |> rbind(data.table(event_time = -1 - p$anticipation, att = 0, se = 0, outcome = p$outcomevar, att_ciub = 0, att_cilb = 0))
  }
  
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



# utils -----------------------------------

g1 <- function(GG){
  if(is.numeric(GG)){return(GG)}
  return(as.numeric(str_split_i(GG, "-", 1)))
}

g2 <- function(GG){
  return(as.numeric(str_split_i(GG, "-", 2)))
}

#' Extract the d-th event timing from the G string.
#' @noRd
gd <- function(GG, d){
  if(is.numeric(GG)){return(GG)}
  return(as.numeric(str_split_i(GG, "-", d)))
}

#' Minimum over all confounding events d != 1 (g' in the paper).
#' @noRd
gprime <- function(GG){
  if(is.numeric(GG)){return(Inf)}
  M <- n_events(GG)
  if(M < 2){return(Inf)}
  result <- gd(GG, 2)
  if(M >= 3){
    for(d in 3:M){
      result <- pmin(result, gd(GG, d))
    }
  }
  result
}

#' Minimum over ALL events (ming = min(g1, gprime)).
#' @noRd
ming <- function(GG){
  if(is.numeric(GG)){return(GG)}
  pmin(g1(GG), gprime(GG))
}

#' Number of events M from a G string vector.
#' @noRd
n_events <- function(GG){
  if(is.numeric(GG)){return(1L)}
  length(str_split(as.character(GG[1]), "-")[[1]])
}

# overridden function -------------------------------------------------

coerce_dt_doub <- function(dt, p){

  M <- 1L + length(p$cohortvar2)      # total number of events
  gcol <- paste0("G", seq_len(M))     # c("G1", "G2", ..., "GM")

  setnames(dt, "G", "G1")
  dt[, mg := Reduce(pmin, .SD), .SDcols = gcol]
  do.call(setorderv, c(list(dt), list(c("time", "mg", gcol, "unit"))))

  #check if there is available never-treated group
  if(!is.infinite(dt[, max(mg)])){
    if(p$control_option == "both"){warning("no never-treated available, effectively using not-yet-but-eventually-treated as control")}
    if(p$control_option == "never"){stop("no never-treated available.")}
  }

  if(p$allow_unbalance_panel){ #let unit start from 1 .... N, useful for knowing which unit is missing
    dt_inv_raw <- dt[dt[, .I[1], by = unit]$V1]
    do.call(setorderv, c(list(dt_inv_raw), list(c("mg", gcol))))
    dt_inv_raw[, new_unit := seq_len(.N)]
    dt <- dt |> merge(dt_inv_raw[,.(unit, new_unit)], by = "unit", sort = FALSE)
    dt[, unit := new_unit]
  }

  #deal with time, coerce time to 1,2,3,4,5.......
  time_periods <- dt[, unique(time)]

  if(length(time_periods) == 0){
    stop("no data after coercing the dataset")
  }

  time_offset <- min(time_periods) - 1 #assume time starts at 1, first is min after sort :)
  if(time_offset != 0){
    dt[, c(gcol) := lapply(.SD, function(x) x - time_offset), .SDcols = gcol]
    dt[, time := time - time_offset]
    time_periods <- time_periods - time_offset
  }

  time_step <- 1 #time may not jump at 1
  if(any(time_periods[2:length(time_periods)] - time_periods[seq_len(length(time_periods)-1)] != 1)){
    # Calculate all intervals between consecutive time periods
    all_intervals <- time_periods[2:length(time_periods)] - time_periods[seq_len(length(time_periods)-1)]

    # Check if all intervals are identical (uniform step)
    if(length(unique(all_intervals)) > 1){
      stop("Time step is not uniform. Time periods: ", paste(head(time_periods, 10), collapse = ", "),
           if(length(time_periods) > 10) "..." else "",
           ". Intervals between periods: ", paste(unique(all_intervals), collapse = ", "),
           ". fastdid requires uniformly-spaced time periods.")
    }

    time_step <- all_intervals[1]
    time_periods <- (time_periods-1)/time_step+1

    # Verify that normalization worked (should always be consecutive integers now)
    if(any(time_periods[2:length(time_periods)] - time_periods[seq_len(length(time_periods)-1)] != 1)){
      stop("Internal error: time normalization failed. Please report this issue.")
    }

    for(g in gcol){
      dt[get(g) != 1, c(g) := (get(g)-1)/time_step+1]
    }

    dt[time != 1, time := (time-1)/time_step+1]
  }
  dt[, mg := Reduce(pmin, .SD), .SDcols = gcol]       # recompute after normalization
  dt[, G := do.call(paste, c(.SD, list(sep="-"))), .SDcols = gcol]   # create G string

  #add the information to t
  t <- list()
  t$time_step <- time_step
  t$time_offset <- time_offset

  if(nrow(dt) == 0){
    stop("no data after coercing the dataset")
  }

  return(list(dt = dt, p = p, t = t))

}

# effect model scheme ----------------------------------------------------------

#' The scheme for a formula-based second stage.
#'
#' A cell equals the pure target effect plus the confounding effect. The user
#' models one component as a linear function of cell features. The model is fit
#' by weighted least squares on the cells where the component is observed alone
#' (the clean cells), and imputed into the confounded cells. A cell is estimable
#' when its regressor row lies in the row space of the fit design. The imputed
#' value is then a linear combination of first-stage cells, so the weights slot
#' into the same contract as `get_es_scheme()`.
#'
#' Three fit modes: `"separate"` fits one component on its clean cells,
#' `"joint"` stacks the events additively and fits on every cell, `"ordered"`
#' does the same for the k-th occurrence of one event kind, where the stacked
#' sum holds by telescoping.
#'
#' @param cells the cell table from `get_cell_table()`.
#' @param att,inf_func the first-stage estimates and influence functions.
#' @return a list: `group_time` (the reported rows), `es_det_weight`,
#'   `es_sto_weight` (zero), and `diag` (weights, cells, fit, overid).
#' @noRd
get_effect_scheme <- function(cells, att, inf_func, aux, p) {
  tab <- build_effect_table(cells, p)
  K <- nrow(tab)

  # the direct rows: a clean cell is the pure target effect, pre-periods included
  direct <- tab[clean1 == TRUE & is.finite(G1)]
  direct[, `:=`(event = 1L, e = e1, component = "direct", estimable = TRUE, leverage = 0, e_gap = 0)]
  rep_rows <- direct
  det <- unit_rows(direct[, cell], K)
  fits <- list()

  if (p$effect_kind != "unrestricted") {
    long <- if (p$effect_fit == "separate") build_component_rows(tab, p) else build_event_rows(tab, p)
    long[, lid := .I] # a cell has one row per modeled event, so the cell index is not unique
    X <- build_design(p$effect_model, long)

    # fit each component, and collect the weight of each modeled row
    for (comp in unique(long[, component])) {
      in_comp <- long[, component == comp]
      lr <- long[in_comp]
      Xc <- X[in_comp, , drop = FALSE]
      fit_long <- lr[, fit]
      tr <- lr[target == TRUE]
      if (nrow(tr) == 0) next
      # a component with no clean cell identifies nothing
      if (!any(fit_long)) {
        tr[, `:=`(estimable = FALSE, leverage = NA_real_, e_gap = NA_real_)]
        if (comp == "confound") tr[, `:=`(event = 1L, e = e1)]
        rep_rows <- rbind(rep_rows, tr, fill = TRUE)
        det <- rbind(det, matrix(0, nrow(tr), K))
        next
      }
      # a stacked design sums the event rows of a cell
      Xf <- rowsum(Xc[fit_long, , drop = FALSE], lr[fit_long, cell])
      fit_cells <- as.integer(rownames(Xf))
      fit <- wls_fit(Xf, as.vector(att)[fit_cells], tab[fit_cells, pg], p$exper$effect_tol)
      fit$cells <- fit_cells
      fit$component <- comp
      fits[[comp]] <- fit

      # the event-time gap to the last fit cell of the same cohort and event
      own_max <- lr[fit == TRUE, .(emax = max(e)), by = .(G, event)]
      tr <- merge(tr, own_max, by = c("G", "event"), all.x = TRUE, sort = FALSE)
      tr[, e_gap := pmax(e - emax, 0)]
      Xt <- Xc[match(tr[, lid], lr[, lid]), , drop = FALSE]

      det_comp <- matrix(0, nrow(tr), K)
      est <- logical(nrow(tr))
      lev <- rep(NA_real_, nrow(tr))
      for (i in seq_len(nrow(tr))) {
        x <- Xt[i, ]
        est[i] <- is_estimable(x, fit, p$exper$effect_tol)
        if (!est[i]) next
        lev[i] <- as.numeric(t(x) %*% fit$pinv %*% x)
        w <- impute_weights(x, fit)
        det_comp[i, fit$cells] <- w
        if (comp == "confound") { # the interacted target effect is the cell minus the confounding effect
          det_comp[i, ] <- -det_comp[i, ]
          det_comp[i, tr[i, cell]] <- det_comp[i, tr[i, cell]] + 1
        }
      }
      tr[, `:=`(estimable = est, leverage = lev)]
      # the reported row is the effect of event 1, whichever component was modeled
      if (comp == "confound") tr[, `:=`(event = 1L, e = e1)]
      rep_rows <- rbind(rep_rows, tr, fill = TRUE)
      det <- rbind(det, det_comp)
    }
  }

  return(finish_effect_scheme(rep_rows, det, tab, fits, p, att, inf_func))
}

#' Order the reported rows, drop the cells that are not estimable, and build
#' the diagnostics.
#' @noRd
finish_effect_scheme <- function(rep_rows, det, tab, fits, p, att, inf_func) {
  M <- 1L + length(p$cohortvar2)
  gcol <- paste0("G", seq_len(M))
  rep_rows[, ri := .I]
  # a cohort treated and confounded in the same period is tried under both
  # components; the target component wins when both are estimable
  rep_rows[, prio := match(component, c("direct", "target", "confound", "joint"))]
  setorder(rep_rows, cell, event, prio)
  rep_rows[, keep := estimable & !duplicated(paste(cell, event, estimable))]
  do.call(setorderv, c(list(rep_rows), list(c("time", "mg", gcol, "event", "prio"))))
  det <- det[rep_rows[, ri], , drop = FALSE]

  # a cohort treated and confounded in the same period has no clean cell, so
  # only the other cells are worth a warning
  surprise <- rep_rows[estimable == FALSE & G1 != gprime(G)]
  if (nrow(surprise) > 0) {
    warning(nrow(surprise), " cell(s) are not estimable under the effect model and are dropped. ",
            "see `effect_diag$cells` with `full = TRUE`.")
  }
  keep <- rep_rows[, keep]
  if (!any(keep)) {
    stop("no cell is estimable under the effect model.")
  }
  if (!rep_rows[keep, any(time >= G1 - p$anticipation)]) {
    warning("no event-specific post-period effect is identified for any cohort. ",
            "check that some cohort has the first event at a different time than the confounding events.")
  }

  det_keep <- det[keep, , drop = FALSE]
  colnames(det_keep) <- tab[, paste0(G, ".", time)]
  rownames(det_keep) <- rep_rows[keep, paste0(G, ".", time, ".", event)]

  diag <- list(
    weights = det_keep,
    cells = rep_rows[, .(G, time, event, component, estimable, leverage, e_gap)],
    fit = lapply(fits, function(f) fit_summary(f, att, tab)),
    overid = lapply(fits, get_effect_overid, att, inf_func)
  )

  out <- rep_rows[keep, c("G", "time", "pg", "mg", gcol, "event", "e", "component"), with = FALSE]
  return(list(group_time = out, es_det_weight = det_keep, es_sto_weight = det_keep * 0, diag = diag))
}

# cell features ----------------------------------------------------------------

#' The cell table with the features that a formula can use.
#'
#' Adds the cell index, the event times `e1..eM`, the earliest confounding
#' date `gp`, the confounding profile `gconf`, and the windows: `active1` and
#' `clean1` for the target effect, `activec` and `cleanc` for the confounding
#' effect. The windows use the anticipation of each event.
#' @noRd
build_effect_table <- function(cells, p) {
  tab <- copy(cells)
  M <- 1L + length(p$cohortvar2)
  tab[, cell := .I]
  tab[, ord := .I]
  tab[, gvec := G]
  tab[, t := as.numeric(time)]
  for (d in seq_len(M)) {
    tab[, (paste0("e", d)) := t - get(paste0("G", d))]
  }
  tab[, gp := gprime(G)]
  tab[, gconf := str_split_fixed(G, "-", 2)[, 2]]
  a1 <- p$anticipation
  a2 <- p$anticipation2
  tab[, active1 := t >= G1 - a1]
  tab[, clean1 := t < gp - a2]
  tab[, activec := t >= gp - a2]
  tab[, cleanc := t < G1 - a1]
  return(tab)
}

#' One row per cell and modeled component, for the separate fit.
#'
#' The target component has one row per cell where event 1 is active. Its fit
#' set is the clean cells. The confounding component has one row per cell where
#' a confounding event is active. Its fit set is the cells before event 1.
#' `gactive` is the confounding profile active at `t`, so a control cohort with
#' an extra active event gets its own column.
#' @noRd
build_component_rows <- function(tab, p) {
  M <- 1L + length(p$cohortvar2)
  a2 <- p$anticipation2

  target <- tab[active1 == TRUE & is.finite(G1)]
  target[, `:=`(component = "target", event = 1L, gown = as.character(G1), gactive = as.character(G1),
                ghist = as.character(G1), e = e1)]
  target[, fit := clean1]
  # a cohort treated and confounded in the same period is tried under both components
  target[, target := !clean1 & G1 <= gp]

  confound <- tab[activec == TRUE]
  confound[, `:=`(component = "confound", event = 2L, gown = gconf, ghist = gconf, e = t - gp)]
  active_dates <- lapply(2:M, function(d) {
    confound[, ifelse(t >= get(paste0("G", d)) - a2, as.character(get(paste0("G", d))), "Inf")]
  })
  confound[, gactive := do.call(paste, c(active_dates, list(sep = "-")))]
  confound[, fit := cleanc]
  confound[, target := active1 & G1 >= gp]

  rows <- rbind(target, confound)
  return(rows)
}

#' One row per cell and active event, for the stacked fits.
#'
#' The design row of a cell is the sum of its event rows. `ghist` is the
#' cohort vector truncated at the event, so a formula can stratify on the
#' history of earlier events.
#' @noRd
build_event_rows <- function(tab, p) {
  M <- 1L + length(p$cohortvar2)
  a <- c(p$anticipation, rep(p$anticipation2, M - 1))
  rows <- rbindlist(lapply(seq_len(M), function(d) {
    r <- tab[t >= get(paste0("G", d)) - a[d]]
    if (nrow(r) == 0) return(NULL)
    r[, `:=`(component = "joint", event = as.integer(d), gown = as.character(get(paste0("G", d))),
             e = get(paste0("e", d)))]
    r[, gactive := gown]
    r[, ghist := do.call(paste, c(lapply(seq_len(d), function(k) as.character(get(paste0("G", k)))), list(sep = "-")))]
    r
  }))
  rows[, fit := TRUE]
  # a clean cell of event 1 is reported as it is
  rows[, target := !(event == 1L & clean1)]
  return(rows)
}

#' The design matrix of the long rows.
#'
#' One `model.matrix()` call over every row, so the factor levels and the
#' polynomial bases are shared between the fit rows and the target rows.
#' @noRd
build_design <- function(formula, long) {
  vars <- all.vars(formula)
  allowed <- names(long)
  bad <- setdiff(vars, allowed)
  if (length(bad) > 0) {
    stop("the effect model uses ", paste(bad, collapse = ", "), ". the available variables are: ",
         "gvec, t, event, e, gown, gactive, ghist, g1..gM, e1..eM.")
  }
  X <- tryCatch({
    mf <- stats::model.frame(formula, data = long, na.action = stats::na.pass)
    stats::model.matrix(formula, mf)
  }, error = function(err) {
    stop("the effect model cannot be built: ", conditionMessage(err),
         ". an event time `e2..eM` is infinite for a cohort without that event; use `e` and `gown`.")
  })
  if (any(!is.finite(X))) {
    stop("the effect model has a non-finite regressor. an event time `e2..eM` is infinite for a cohort ",
         "without that event; use `e` and `gown`, which refer to the modeled event.")
  }
  return(X)
}

# linear algebra ---------------------------------------------------------------

#' Weighted least squares through the SVD of the normal matrix.
#'
#' The singular values below `tol` times the largest are treated as zero. The
#' right singular vectors that remain span the row space of the design, which
#' the estimability check reads.
#' @return a list: `pinv`, `rank`, `V`, `X`, `w`, `y`, `coef`, `fitted`, `resid`.
#' @noRd
wls_fit <- function(X, y, w, tol) {
  A <- crossprod(X, X * w)
  s <- svd(A)
  keep <- s$d > tol * max(s$d, .Machine$double.eps)
  V <- s$v[, keep, drop = FALSE]
  pinv <- V %*% (t(V) / s$d[keep])
  coef <- as.vector(pinv %*% crossprod(X, w * y))
  fitted <- as.vector(X %*% coef)
  list(pinv = pinv, rank = sum(keep), V = V, X = X, w = w, y = y,
       coef = stats::setNames(coef, colnames(X)), fitted = fitted, resid = y - fitted)
}

#' Is the regressor row in the row space of the fit design?
#' @noRd
is_estimable <- function(x, fit, tol) {
  r <- x - fit$V %*% crossprod(fit$V, x)
  sqrt(sum(r^2)) <= sqrt(tol) * max(1, sqrt(sum(x^2)))
}

#' The weight of each fit cell in the imputed value `x' beta`.
#' @noRd
impute_weights <- function(x, fit) {
  as.vector(crossprod(fit$pinv %*% x, t(fit$X) * rep(fit$w, each = ncol(fit$X))))
}

#' Rows of the identity, one per direct cell.
#' @noRd
unit_rows <- function(cells, K) {
  m <- matrix(0, length(cells), K)
  m[cbind(seq_along(cells), cells)] <- 1
  m
}

# diagnostics ------------------------------------------------------------------

#' A summary of one fit: the formula, the size, the rank, the coefficients,
#' and the residual of every fit cell.
#' @noRd
fit_summary <- function(fit, att, tab) {
  list(
    component = fit$component,
    n_fit = length(fit$cells),
    rank = fit$rank,
    coef = fit$coef,
    residuals = data.table(G = tab[fit$cells, G], time = tab[fit$cells, time],
                           att = fit$y, fitted = fit$fitted, resid = fit$resid)
  )
}

#' The over-identification test of one fit.
#'
#' The residual of the fit cells is a linear combination of first-stage cells,
#' so its variance comes from the influence functions. The statistic is the
#' residual quadratic form in the pseudo-inverse of that variance, with the
#' degrees of freedom of the null space. The bootstrap does not cover it.
#' @noRd
get_effect_overid <- function(fit, att, inf_func) {
  n <- length(fit$cells)
  H <- fit$X %*% fit$pinv %*% (t(fit$X) * rep(fit$w, each = ncol(fit$X)))
  R <- diag(n) - H
  IF <- inf_func[, fit$cells, drop = FALSE] %*% t(R)
  V <- crossprod(IF) / nrow(IF)^2
  r <- fit$resid
  df <- n - fit$rank
  if (df <= 0) return(list(stat = 0, df = 0, pvalue = NA_real_))
  s <- svd(V)
  keep <- s$d > 1e-10 * max(s$d)
  Vinv <- s$v[, keep, drop = FALSE] %*% (t(s$u[, keep, drop = FALSE]) / s$d[keep])
  stat <- as.numeric(t(r) %*% Vinv %*% r)
  list(stat = stat, df = df, pvalue = stats::pchisq(stat, df, lower.tail = FALSE))
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
                                                    control = parglm::parglm.control(nthreads = ifelse(p$parallel, 1, getDTthreads())),
                                                    intercept = FALSE))
      class(prop_score_est) <- "glm" #trick the vcov function to think that this is a glm object to dispatch the write method
      #const is implicitly put into the ipw formula, need to incorporate it manually

      logit_coef <-  prop_score_est$coefficients
      
      if(anyNA(logit_coef)){
        warning("some propensity score estimation resulted in NA coefficients, likely cause by perfect colinearity")
      }
      
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
    if(max(dt_did[, cont_ipw_weight]) > 100){warning("extreme IPW weights detected (max: ", round(max(dt_did[, cont_ipw_weight]), 2), "), estimates may be unstable")}

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
    asym_linear_or <- t(solve(XpX, t(or_ex)))

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

  # the residuals are centered on the estimated group means, which deflates the
  # plug-in variance by (m-1)/m for a group of m effective units. inflate by the
  # Kish effective size of each group, else small cells under-cover.
  ess_treat <- dt_did[, sum(treat_ipw_weight)^2/sum(treat_ipw_weight^2)]
  ess_cont <- dt_did[, sum(cont_ipw_weight)^2/sum(cont_ipw_weight^2)]

  # a group of one effective unit has a zero residual, so its variance is not
  # estimable and the cell must be skipped
  if(!is.finite(ess_treat) || !is.finite(ess_cont) || ess_treat < 2 || ess_cont < 2){
    stop("a group has fewer than 2 effective units, the variance is not estimable")
  }
  inf_treat_did <- inf_treat_did * sqrt(ess_treat/(ess_treat-1))
  inf_cont_did <- inf_cont_did * sqrt(ess_cont/(ess_cont-1))

  #get overall influence function
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

  if(n_pre == 0 || n_post == 0){
    warning("No observations in pre or post period; skipping this 2x2 DiD")
    return(list(att = NA, inf_func = rep(0, oldn), cache = NULL))
  }

  # check for enough treated and control observations in each period
  any_zero <- dt_did[, sum(D == 0 & inpre)==0] | dt_did[, sum(D == 0 & inpost)==0] | dt_did[, sum(D == 1 & inpre)==0] | dt_did[, sum(D == 1 & inpost)==0]
  if(any_zero){
    stop("Not enough treated or control observations in pre or post period")
  }

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
    prop_score_est <- suppressWarnings(parglm::parglm.fit(covvars, dt_did[, D],
                                                  family = stats::binomial(),
                                                  weights = dt_did[, weights*(inpre+inpost)*n/(n_pre+n_post)], #when seen in both pre and post have double weight
                                                  control = parglm::parglm.control(nthreads = ifelse(p$parallel, 1, getDTthreads())),
                                                  intercept = FALSE)) #*(inpre+inpost)
    class(prop_score_est) <- "glm" #trick the vcov function to think that this is a glm object to dispatch the write method
    #const is implicitly put into the ipw formula, need to incorporate it manually
    
    #for the influence, will be cached
    hess <- stats::vcov(prop_score_est) * n #for the influence function
    
    logit_coef <-  prop_score_est$coefficients 
    prop_score_fit <- fitted(prop_score_est)
    if(max(prop_score_fit) >= 1-1e-10){warning(paste0("extreme propensity score: ", max(prop_score_fit), ", support overlap is likely to be violated"))}
    prop_score_fit <- pmin(1-1e-10, prop_score_fit) #for the ipw
    
    #get the results into the main did dt
    dt_did[, ps := prop_score_fit]
    dt_did[, treat_ipw_weight := weights*D]
    dt_did[, cont_ipw_weight := weights*ps*(1-D)/(1-ps)]
    if(max(dt_did[, cont_ipw_weight]) > 100){warning("extreme IPW weights detected (max: ", round(max(dt_did[, cont_ipw_weight]), 2), "), estimates may be unstable")}
    
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

    # ensure sufficient observations for outcome regression
    if(sum(control_bool_post) <= ncol(covvars) || sum(control_bool_pre) <= ncol(covvars)){
      stop("Not enough control observations to estimate outcome regression")
    }
    reg_coef_post <- stats::coef(stats::lm.wfit(x = covvars[control_bool_post,], y = dt_did[control_bool_post,post.y],
                                                w = dt_did[control_bool_post,weights]))

    reg_coef_pre <- stats::coef(stats::lm.wfit(x = covvars[control_bool_pre,], y = dt_did[control_bool_pre,pre.y],
                                               w = dt_did[control_bool_pre,weights]))

    if(anyNA(reg_coef_post) || anyNA(reg_coef_pre)){
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

  # small-group inflation, see estimate_did_bp: each group-period mean deflates
  # the plug-in variance by (m-1)/m for m effective units. a group-period of one
  # effective unit has a zero residual, so the cell must be skipped
  ess_all <- sapply(list(dt_did[, treat_ipw_weight*inpost], dt_did[, cont_ipw_weight*inpost],
                         dt_did[, treat_ipw_weight*inpre], dt_did[, cont_ipw_weight*inpre]),
                    function(w) sum(w)^2/sum(w^2))
  if(any(!is.finite(ess_all)) || any(ess_all < 2)){
    stop("a group-period has fewer than 2 effective units, the variance is not estimable")
  }
  inf_treat_did_post <- inf_treat_did_post * sqrt(ess_all[1]/(ess_all[1]-1))
  inf_cont_did_post <- inf_cont_did_post * sqrt(ess_all[2]/(ess_all[2]-1))
  inf_treat_did_pre <- inf_treat_did_pre * sqrt(ess_all[3]/(ess_all[3]-1))
  inf_cont_did_pre <- inf_cont_did_pre * sqrt(ess_all[4]/(ess_all[4]-1))
  
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

#' gtatt for each outcome.
#' @noRd
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

    # cells skipped for a too-small group are reported once, not one warning per cell
    ess_skip <- vapply(gt_results, function(x) isTRUE(x$skip_ess), logical(1))
    if(any(ess_skip)){
      warning(sum(ess_skip), " group-time(s) skipped: a group has fewer than 2 effective units, so the variance is not estimable")
      gt_results <- gt_results[!ess_skip]
    }
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

#' gtatt for each outcome, each gt.
#' @noRd
estimate_gtatt_outcome_gt <- function(gt, y, aux, p, caches){
  
  g <- gt[1]
  t <- as.numeric(gt[2])
  
  # skip g-t outside requested event-time range
  event_time <- t - ming(g)
  if (!is.na(p$exper$only_est_min) && event_time < p$exper$only_est_min) return(NULL)
  if (!is.na(p$exper$only_est_max) && event_time > p$exper$only_est_max) return(NULL)

  #find base time
  gt_name <- paste0(g,".",t)
  base_period <- get_base_period(g,t,p)
  if(t == base_period || #no treatment effect for the base period
     !base_period %in% aux$time_periods){ #base period out of bounds
    return(NULL)
  } 
  #find treatment and control group
  did_setup <- get_did_setup(g,t, base_period, aux, p)
  valid_tc_groups <- any(did_setup == 1) && any(did_setup == 0) #if takes up too much time, consider use collapse anyv, but right now quite ok
  if(!isTRUE(valid_tc_groups)){return(NULL)} #no treatment group or control group #isTRUE for na as false
  
  #covariates matrix
  covvars <- get_covvars(base_period, t, aux, p)
  
  #the 2x2 dataset
  cohort_did <- data.table(did_setup, y[[t]], y[[base_period]], aux$weights)
  names(cohort_did) <- c("D", "post.y", "pre.y", "weights")

  # skip if no valid observation in either period when allowing for
  # unbalanced panels
  if(sum(!is.na(cohort_did$pre.y)) == 0 || sum(!is.na(cohort_did$post.y)) == 0){
    return(NULL)
  }

  # estimate --------------------
  result <- tryCatch(estimate_did(dt_did = cohort_did, covvars, p, caches[[gt_name]]),
                     error = function(e){
                       # small-group skips are counted and reported once by the caller
                       if(grepl("fewer than 2 effective units", e$message, fixed = TRUE)){
                         return("skip_ess")
                       }
                       warning("Skipping group-time ", g, "-", t,
                               ": ", e$message)
                       return(NULL)
                     })
  if(is.null(result)){return(NULL)}
  if(identical(result, "skip_ess")){return(list(gt = gt, skip_ess = TRUE))}
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
  if(allNA(p$cohortvar2) || p$anticipation == p$anticipation2){
    # one horizon for every event, so the cohorts in range are contiguous after the sort
    control_pos <- get_control_pos(aux$cohort_sizes, min_control_cohort, max_control_cohort)
  } else {
    # each event has its own horizon, so screen the cohorts one by one
    control_pos <- get_control_pos_event(aux$cohort_sizes, t, base_period, max_control_cohort, p)
  }
  did_setup[control_pos] <- 0
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
  
  # Validate sequence parameters
  if(start > end || start <= 0 || end <= 0){
    return(c())  # Return empty vector when no valid control cohorts
  }
  return(seq(start, end, by = 1))
}

#' Control positions when the events have different anticipation horizons.
#'
#' A cohort is a valid not-yet-treated control only if every event is further
#' away than the horizon of that event. The first event uses `anticipation`, the
#' confounding events use `anticipation2`. The cohorts in range are not
#' contiguous after the sort, so screen them one by one.
#'
#' @param cohort_sizes the cohort table, in the order of the unit array.
#' @param t,base_period the two periods of the 2x2.
#' @param max_control_cohort the not-yet-treated upper bound.
#' @return the positions of the control units in the unit array.
#' @noRd
get_control_pos_event <- function(cohort_sizes, t, base_period, max_control_cohort, p){
  GG <- cohort_sizes[, G]
  last <- max(t, base_period)

  if(p$control_option == "never"){
    keep <- is.infinite(ming(GG))
  } else {
    keep <- g1(GG) > last + p$anticipation
    M <- 1L + length(p$cohortvar2)
    for(d in 2:M){
      keep <- keep & (gd(GG, d) > last + p$anticipation2)
    }
  }
  keep <- keep & (ming(GG) <= max_control_cohort)
  if(!any(keep)){return(c())}

  # each cohort is one contiguous block of the unit array
  end <- cumsum(cohort_sizes[, cohort_size])
  start <- end - cohort_sizes[, cohort_size] + 1
  return(unlist(lapply(which(keep), function(i) seq(start[i], end[i]))))
}

get_treat_pos <- function(cohort_sizes, treat_cohort){ #need to separate for double did to match exact g-g-t
  index <- which(cohort_sizes[,G] == treat_cohort)
  if(length(index) == 0){
    stop("Cohort ", treat_cohort, " not found in cohort_sizes")
  }
  start <- ifelse(index == 1, 1, cohort_sizes[1:(index-1), sum(cohort_size)]+1)
  end <- cohort_sizes[1:index, sum(cohort_size)]
  
  # Validate sequence parameters
  if(start > end || start <= 0 || end <= 0){
    return(c())  # Return empty vector when no valid treat positions
  }
  return(seq(start, end, by = 1))
}

get_covvars <- function(base_period, t, aux, p){
  
  if(allNA(p$covariatesvar) && allNA(p$varycovariatesvar)){return(NA)}
  covvars <- data.table()
  
  #add time-varying covariates
  if(!allNA(p$varycovariatesvar)){
    
    precov <- aux$varycovariates[[base_period]]
    names(precov) <- paste0("pre_", names(precov))
    postcov <- aux$varycovariates[[t]]-aux$varycovariates[[base_period]]
    names(postcov) <- paste0("post_", names(aux$varycovariates[[t]]))
    covvars <- cbind(covvars, cbind(precov, postcov))
  }
  
  #add time-invariant covariates
  if(!allNA(p$covariatesvar)){
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
#' @param result_type character, type of result to return, options are "group_time", "time", "group", "simple", "dynamic" (time since event), "group_group_time", "dynamic_stagger", or "dynamic_event" (the average effect of every event at each event time, for `effect_fit = "ordered"`).
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
#' @param anticipation2 number, periods with anticipation for the second event.
#' @param exper list, arguments for experimental features. Supported options:
#'   \describe{
#'     \item{`only_est_min`}{numeric scalar, minimum event time to estimate (`result_type == "dynamic"` only, not compatible with double DiD).}
#'     \item{`only_est_max`}{numeric scalar, maximum event time to estimate (`result_type == "dynamic"` only, not compatible with double DiD).}
#'     \item{`filtervar`}{character, name of a logical column; only units with TRUE at the base period are used.}
#'     \item{`filtervar_post`}{character, name of a logical column; only units with TRUE at the post period are used.}
#'     \item{`only_balance_2by2`}{logical, keep only units observed in both periods of each 2x2 DiD.}
#'     \item{`aggregate_scheme`}{character, a custom aggregation expression evaluated as `group_time[, target := <expr>]`.}
#'     \item{`max_control_cohort_diff`}{numeric, maximum cohort difference between treated and control groups.}
#'     \item{`effect_tol`}{numeric, the relative tolerance of the rank and estimability checks of an effect model. Default is 1e-7.}
#'   }
#' @param base_period character, type of base period in pre-preiods, options are "universal", or "varying".
#' @param full logical, whether to return the full result (influence function, call, weighting scheme, etc,.).
#' @param parallel logical, whether to use parallization on unix system.
#' @param cohortvar2 character or character vector, name(s) of the confounding event cohort variable(s). For M>2 events, provide a vector of length M-1 (e.g., `c("G2", "G3")` for M=3 events).
#' @param event_specific logical, whether to recover target treatment effect or use combined effect.
#' @param double_control_option character, control units used for the double DiD, options are "both", "never", or "notyet". "notyet" keeps a control cohort only if every confounding event is finite and later than the period, so the control must be confounded eventually by all of them. With M >= 3 events this can leave very few control cohorts, and "both" is the recommended option.
#' @param add_base_period logical, whether to add a placeholder base period in dynamic results.
#' @param effect_model the second-stage model of the treatment effects, for multiple events. `"parallel"` (the default) is the parallel treatment effects estimator. `"unrestricted"` reports the clean cells only. A one-sided formula models one component of the cell as a linear function of the cell features, see Details.
#' @param effect_fit character, how a formula is fit. `"separate"` fits the modeled component on its clean cells. `"joint"` stacks the events additively and fits on every cell. `"ordered"` does the same for the k-th occurrence of one event kind, where the cohort dates must increase.
#'
#' @import data.table stringr dreamerr ggplot2
#' @importFrom stats quantile vcov sd binomial fitted qnorm rnorm as.formula weighted.mean model.frame model.matrix na.pass pchisq setNames
#' @importFrom parglm parglm.fit parglm.control
#' @importFrom collapse allNA fnrow whichNA fnunique fsum na_insert
#' @importFrom parallel mclapply
#' @importFrom BMisc multiplier_bootstrap
#' @importFrom utils head
#' @return A data.table containing the estimated treatment effects and standard errors or a list of all results when `full == TRUE`.
#' @export
#'
#' @details
#' `balanced_event_time`, `add_base_period`, and the `exper` options `only_est_min`/`only_est_max` are only meaningful when `result_type == "dynamic"`.
#'
#' `result_type` as `"group_group_time"` and `"dynamic_stagger"` are only meaningful when using double DiD (`cohortvar2` is set).
#'
#' `cohortvar2` accepts a character vector of length M-1 to support M>2 treatment events.
#'
#' `biters` and `clustervar` are only used when `boot == TRUE`.
#'
#' **Effect models.** With multiple events, a first-stage cell (cohort vector, period) identifies the combined effect of every active event. The second stage recovers the effect of event 1. An `effect_model` formula states a linear model for one component of the cell, fit by weighted least squares on the cells where that component is observed alone, and imputed into the other cells. A cell is reported only when it is estimable: its regressor row lies in the row space of the fit design. The formula can use these cell features:
#'   \describe{
#'     \item{`gvec`}{the cohort vector, as the string `"g1-g2-...-gM"`}
#'     \item{`t`}{the period}
#'     \item{`g1`, ..., `gM`, `e1`, ..., `eM`}{the date and the event time of each event}
#'     \item{`event`}{the modeled event: 1 for the target, 2 for the confounding bundle, k in the stacked fits}
#'     \item{`gown`}{the date of the modeled event (the confounding profile for the bundle)}
#'     \item{`gactive`}{the dates of the modeled events that are active at `t`}
#'     \item{`ghist`}{the cohort vector truncated at the modeled event}
#'     \item{`e`}{the event time of the modeled event}
#'   }
#' The model `~ factor(gvec) + factor(gactive):factor(t)` is the parallel treatment effects assumption of `"parallel"`, fit on every clean cell at once. `~ factor(gvec) + poly(e, 2)` is a quadratic profile in event time with a cohort level. With `effect_fit = "separate"`, the target component is modeled for cohorts treated before they are confounded, and the confounding component for cohorts confounded before they are treated; the reported effect is the pure effect in the first case and the interacted effect in the second. With `effect_fit = "joint"` or `"ordered"`, the design row of a cell is the sum over its active events, so a formula can share parameters across events, for example `~ 0 + factor(gvec):factor(event) + factor(e)`. Use `~ 0 + ...` in a stacked fit, because an intercept counts once per active event. `result_type = "dynamic_event"` averages the effect of every event at each event time, with cohort-size weights. The weights of the fit are treated as fixed in the influence function. `double_control_option` and `control_option = "notyet"` do not restrict the fit set of a formula. With `full = TRUE`, `effect_diag` returns the weight of every first-stage cell in every reported cell, the estimability of every candidate cell, the fit residuals, and an over-identification test.
#'
#' @examples
#' # simulated data
#' simdt <- sim_did(1e+02, 10, cov = "cont", second_cov = TRUE, second_outcome = TRUE, seed = 1)
#' dt <- simdt$dt
#'
#' # basic call
#' result <- fastdid(
#'   data = dt, timevar = "time", cohortvar = "G",
#'   unitvar = "unit", outcomevar = "y",
#'   result_type = "group_time"
#' )
#'
#' @keywords difference-in-differences fast computation panel data estimation did
fastdid <- function(data,
                    timevar, cohortvar, unitvar, outcomevar,
                    control_option = "both", result_type = "group_time", balanced_event_time = NA,
                    control_type = "ipw", allow_unbalance_panel = FALSE, boot = FALSE, biters = 1000, cband = FALSE, alpha = 0.05,
                    weightvar = NA, clustervar = NA, covariatesvar = NA, varycovariatesvar = NA,
                    copy = TRUE, validate = TRUE,
                    anticipation = 0, anticipation2 = 0, base_period = "universal",
                    exper = NULL, full = FALSE, parallel = FALSE,
                    cohortvar2 = NA, event_specific = TRUE, double_control_option = "both", add_base_period = FALSE,
                    effect_model = "parallel", effect_fit = "separate") {
  # preprocess --------------------------------------------------------

  if (!is.data.table(data)) {
    warning("coercing input into a data.table.")
    data <- as.data.table(data)
  }
  if (copy) {
    dt <- copy(data)
  } else {
    dt <- data
  }

  # the presets are matched here, so that the rest of the code reads one value
  if (is.character(effect_model)) {
    effect_model <- match.arg(effect_model, c("parallel", "unrestricted"))
  }
  effect_kind <- if (is.character(effect_model)) effect_model else "formula"

  # validate arguments
  p <- as.list(environment()) # collect everything besides data
  p$data <- NULL
  p$dt <- NULL

  exper_args <- c(
    "filtervar", "filtervar_post", "only_balance_2by2",
    "aggregate_scheme", "max_control_cohort_diff",
    "only_est_min", "only_est_max", "effect_tol"
  )
  p$exper <- get_exper_default(p$exper, exper_args)
  class(p) <- "locked" # no more changes!
  validate_argument(dt, p)

  # change name for main columns
  setnames(dt, c(timevar, cohortvar, unitvar), c("time", "G", "unit"))
  if (!allNA(p$cohortvar2)) {
    new_names <- paste0("G", seq(2L, 1L + length(p$cohortvar2)))
    setnames(dt, p$cohortvar2, new_names)
  }

  # validate and throw away not legal data
  dt <- validate_dt(dt, p)

  # make dt conform to the WLOG assumptions of fastdid
  coerce_result <- coerce_dt(dt, p) # also changed dt

  # get auxiliary data
  aux <- get_auxdata(coerce_result$dt, p)

  # main estimation  -------------------------------------------------

  gt_result_list <- estimate_gtatt(aux, p)
  agg_result <- aggregate_gt(gt_result_list, aux, p)

  # post process -------------------------------------------

  # convert "targets" back to meaningful parameters
  est_results <- convert_targets(agg_result$est, p, coerce_result$t)

  if (!p$full) {
    return(est_results)
  } else {
    full_result <- list(
      call = p,
      estimate = est_results,
      gt_estimate = gt_result_list,
      agg_inf_func = agg_result$inf_func,
      agg_weight_matrix = agg_result$agg_weight_matrix,
      es_weight_matrix = agg_result$es_weight_matrix,
      effect_diag = agg_result$effect_diag
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
#' # estimation
#' result <- fastdid(
#'   data = dt, timevar = "time", cohortvar = "G",
#'   unitvar = "unit", outcomevar = "y",
#'   result_type = "dynamic"
#' )
#'
#' # plot
#' plot_did_dynamics(result)
#'
#' @export
plot_did_dynamics <- function(x, margin = "event_time") {
  # find the base_period
  if (margin == "event_time") {
    et_range <- min(x[, event_time]):max(x[, event_time])
    base_time <- et_range[!et_range %in% x[, unique(event_time)]]
    if (length(base_time) > 1) {
      stop("missing more than one period")
    }

    # add the base period; a result with post periods only has none to add
    if (length(base_time) == 1) {
      if ("outcome" %in% names(x)) {
        base_row <- data.table(att = 0, se = 0, event_time = base_time, outcome = x[, unique(outcome)], att_ciub = 0, att_cilb = 0)
      } else {
        base_row <- data.table(att = 0, se = 0, event_time = base_time, att_ciub = 0, att_cilb = 0)
      }
      x <- x |> rbind(base_row, fill = TRUE)
    }
  } else {
    x <- x[type == "post"]
  }

  plot <- x |>
    ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", col = "red") +
    geom_point(aes(x = eval(str2lang(margin)), y = att), color = "black") + # point est
    geom_errorbar(aes(x = eval(str2lang(margin)), ymin = att_cilb, ymax = att_ciub),
      width = 0.1, linetype = "dashed"
    ) + # CI
    labs(x = margin)

  if (margin == "event_time") {
    plot <- plot + geom_line(aes(x = eval(str2lang(margin)), y = att), color = "black") # point est
  }

  return(plot)
}

NULL

# quiets concerns of R CMD check re: the .'s that appear in pipelines and data.table variables
utils::globalVariables(c(
    ".", "agg_weight", "att", "att_cont", "att_treat", "attgt", "cohort", "cohort_size", "conf_lwb", "conf_upb",
    "const", "cont_ipw_weight", "count", "delta_y", "element_rect", "element_text", "event_time", "pg", "placeholder",
    "post.y", "pre.y", "ps", "s", "se", "target", "tau", "time_fe",
    "treat_ipw_weight", "treat_latent", "type", "unit", "unit_fe", "weight", "x", "x2",
    "x_trend", "y", "y0", "y1", "y2", "time", "weights", "outcome", "G", "D", "xvar",
    "V1", "att_cont_post", "att_cont_pre", "att_treat_post", "att_treat_pre", "inpost", "inpre", "max_et", "min_et", "new_unit", "or_delta", "or_delta_post", "or_delta_pre", "targeted", "used",
    "timevar", "cohortvar", "unitvar", "outcomevar", "control_option", "result_type", "balanced_event_time", "control_type",
    "allow_unbalance_panel", "boot", "biters", "weightvar", "clustervar", "covariatesvar", "varycovariatesvar", "filtervar",
    "copy", "validate", "max_control_cohort_diff", "anticipation", "anticipation2", "min_control_cohort_diff", "base_period", "post", "att_ciub", "att_cilb", "cband", "alpha",
    "G2", "G1", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10",
    "mg", "cohort1", "cohort2", "cohort3", "cohort4", "cohort5", "event_time_1", "event_time_2",
    "D2", "attgt2", "event", "atu2", "y01", "y10", "y11", "tau2", "parallel",
    "tp", "cp", "tb", "cb", "no_na", "event_stagger", "double_control_option",
    "det_weight", "sto_weight", "add_base_period", "cohortvar2", "exper",
    "effect_model", "effect_fit", "effect_kind", "cell", "ord", "gvec", "t", "gp", "gconf",
    "active1", "clean1", "activec", "cleanc", "component", "estimable", "leverage", "e_gap",
    "gown", "gactive", "ghist", "e", "e1", "fit", "emax", "resid", "ri", "prio", "keep", "lid"
))

# cell table -------------------------------------------------------------------

#' Build the cell table from the first-stage result.
#'
#' One row per first-stage cell, in the order of the first-stage estimates. The
#' table carries the cohort share `pg`, the first event `mg`, and one column per
#' event date `G1..GM`. The second stage and the aggregation read the cells from
#' this table only.
#'
#' @param gt_result the first-stage result of one outcome: `gt`, `att`, `inf_func`.
#' @return the cell table.
#' @noRd
get_cell_table <- function(gt_result, aux, p) {
  id_dt <- data.table(weight = aux$weights / sum(aux$weights), G = aux$dt_inv[, G])
  pg_dt <- id_dt[, .(pg = sum(weight)), by = "G"]
  cells <- gt_result$gt |> merge(pg_dt, by = "G", sort = FALSE)
  cells[, mg := ming(G)]
  M <- if (allNA(p$cohortvar2)) 1L else 1L + length(p$cohortvar2)
  gcol <- paste0("G", seq_len(M))
  for (d in seq_len(M)) {
    cells[, (paste0("G", d)) := gd(G, d)]
  }
  do.call(setorderv, c(list(cells), list(c("time", "mg", gcol)))) # match order in gtatt
  if (!all(names(gt_result$att) == cells[, paste0(G, ".", time)])) {
    stop("some bug makes gt misaligned, please report this to the maintainer. Thanks.")
  }
  return(cells)
}

# second stage -----------------------------------------------------------------

#' Apply the second stage to the first-stage cells.
#'
#' The second stage is a linear combination of first-stage cells. The scheme
#' gives one weight row per identified cell, split into a deterministic part and
#' a part that depends on the estimated cohort shares. The estimate is the
#' weighted sum, and its influence function is the same combination of the
#' first-stage influence functions plus the influence of the shares.
#'
#' @param cells the cell table from `get_cell_table()`.
#' @param att,inf_func the first-stage estimates and influence functions.
#' @return a list: the identified cells, their `att` and `inf_func`, and the
#'   weight matrix over the first-stage cells.
#' @noRd
second_stage <- function(cells, att, inf_func, aux, p) {
  if (p$effect_kind == "parallel") {
    sch <- get_es_scheme(cells, aux, p)
  } else {
    sch <- get_effect_scheme(cells, att, inf_func, aux, p)
  }
  det_weight <- as.matrix(sch$es_det_weight)
  sto_weight <- as.matrix(sch$es_sto_weight)
  es_weight <- det_weight + sto_weight

  # the share weights are signed, and each period is normalized on its own
  if (any(sto_weight != 0)) {
    pre_cells <- copy(cells)
    pre_cells[, pg := NULL] # get_weight_influence merges the shares itself
    es_inf_weights <- get_weight_influence(att, pre_cells, sto_weight, aux, p, by_period = TRUE)
  } else {
    es_inf_weights <- 0
  }

  out_cells <- sch$group_time
  att <- es_weight %*% att
  inf_func <- (inf_func %*% t(es_weight)) + es_inf_weights

  # a stacked fit reports every event; the event-1 rows serve the other result types
  if ("event" %in% names(out_cells) && p$result_type != "dynamic_event" && is.na(p$exper$aggregate_scheme)) {
    keep <- out_cells[, event == 1L]
    out_cells <- out_cells[keep]
    att <- att[keep, , drop = FALSE]
    inf_func <- inf_func[, keep, drop = FALSE]
    es_weight <- es_weight[keep, , drop = FALSE]
  }

  return(list(
    cells = out_cells,
    att = att,
    inf_func = inf_func,
    weight = es_weight,
    diag = sch$diag
  ))
}

# parallel treatment effects scheme -------------------------------------------

#' The scheme for the event-specific effect.
#' @noRd
get_es_scheme <- function(group_time, aux, p){

  es_group_time <- copy(group_time) #group_time with available es effect
  #create lookup (columns already populated by get_cell_table)
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

  # a cohort with g1 == g' carries no separable effect, so the post-periods can all be gone
  if(nrow(es_group_time) == 0 || !es_group_time[, any(time >= G1 - p$anticipation)]){
    warning("no event-specific post-period effect is identified for any cohort. ",
            "check that some cohort has the first event at a different time than the confounding events.")
  }

  es_det_weight <- do.call(rbind, lapply(es_weight_list, \(x){x$det}))
  es_sto_weight <- do.call(rbind, lapply(es_weight_list, \(x){x$sto}))

  return(list(group_time = es_group_time, es_det_weight = es_det_weight, es_sto_weight = es_sto_weight))

}

#' Keep the control cohorts that are available at both periods.
#'
#' A first-stage cell can be missing, for example after an estimation failure or
#' with an unbalanced panel. The two periods then normalize over different
#' populations. This function restricts both to the common cohorts.
#'
#' @param group_time the group-time table.
#' @param cp,cb logical vectors, the control rows at t and at the base period.
#' @param gg,t,base_period the target cohort and the two periods, for the message.
#' @return a list with the two restricted logical vectors, or NULL if no cohort is common.
#' @noRd
intersect_control <- function(group_time, cp, cb, gg, t, base_period){
  common <- intersect(group_time[cp, G], group_time[cb, G])
  if(length(common) == 0){
    warning("the control cohorts at ", t, " and at ", base_period,
            " do not overlap for cohort ", gg, ". fastdid skips the cell.")
    return(NULL)
  }
  in_common <- group_time[, G %in% common]
  return(list(cp = cp & in_common, cb = cb & in_common))
}

#' The scheme for the group-group-time estimates.
#' Implements Theorem 3 of Tsai (2026) for M >= 2 events.
#' @noRd
get_es_ggt_weight <- function(ggt, group_time, aux, p){

  group_time <- copy(group_time) #avoid accidental modification

  group_time[, det_weight := 0] #reset
  group_time[, sto_weight := 0] #reset
  t     <- group_time[ggt, time]
  g1_val <- group_time[ggt, G1]
  gg    <- group_time[ggt, G]

  if(is.infinite(g1_val)){return(NULL)}

  M_val <- 1L + length(p$cohortvar2)   # total number of events
  gp    <- gprime(gg)                  # g' = min_{d!=1}(g^d), earliest confounding event

  # the anticipation of the confounding event contaminates the periods before g', so
  # the direct case stops one anticipation horizon earlier
  if(t < gp - p$anticipation2){ # Case 1: direct pure effect (before any confounding event)

    group_time[ggt, det_weight := 1]

  } else if(g1_val < gp) { # Case 2: imputation (treated before confounded)
    # C^imp = {h : h^1 = g^1, for all d!=1, h^d > t}

    base_period <- gp - 1 - p$anticipation2
    if(base_period == t){return(NULL)}
    min_control_cohort <- ifelse(p$double_control_option == "never", Inf, max(t,base_period)+p$anticipation2+1)

    tb <- group_time[, G == gg & time == base_period]

    # control: same G1, and ALL confounding events not yet occurred (each h^d >= min_control_cohort)
    c <- group_time[, G1 == g1_val]
    for(d in 2:M_val){
      Gd_vals <- group_time[[paste0("G", d)]]
      c <- c & (Gd_vals >= min_control_cohort)
    }
    if(p$control_option == "notyet"){
      # exclude never-confounded units: require each h^d < Inf
      for(d in 2:M_val){
        Gd_vals <- group_time[[paste0("G", d)]]
        c <- c & !is.infinite(Gd_vals)
      }
    }
    cp <- group_time[, c & time == t]
    cb <- group_time[, c & time == base_period]

    #if any group have no available cohort, skip
    if(sum(tb) == 0 | sum(cp) == 0 | sum(cb) == 0){return(NULL)}

    common <- intersect_control(group_time, cp, cb, gg, t, base_period)
    if(is.null(common)){return(NULL)}
    cp <- common$cp
    cb <- common$cb

    #assign the weights
    group_time[tb, det_weight := 1]
    group_time[cp, sto_weight := pg/sum(pg)]
    group_time[cb, sto_weight := -pg/sum(pg)]

  } else if (g1_val > gp) { # Case 3: double DiD (confounded before treated)
    # C^did = {h : h^1 > t, h^d = g^d if g^d <= t, h^d > max(t, base) if g^d > t}

    # the theorem covers the post-period of the first event only
    if(t < g1_val - p$anticipation){return(NULL)}

    base_period <- g1_val - 1 - p$anticipation
    if(base_period == t){return(NULL)}
    min_control_cohort <- ifelse(p$double_control_option == "never", Inf, max(t,base_period)+p$anticipation+1)
    min_conf_cohort <- max(t, base_period) + p$anticipation2 + 1

    tp <- group_time[,.I == ggt]
    tb <- group_time[,G == gg & time == base_period]

    # control: h^1 not yet treated, and for every confounding event either the same
    # timing as the target (g^d <= t), or no event at all yet (g^d > t). without the
    # second rule a control can carry an event the target does not have.
    c <- group_time[, G1 >= min_control_cohort & G1 != g1_val]
    for(d in 2:M_val){
      g_d <- gd(gg, d)
      Gd_vals <- group_time[[paste0("G", d)]]
      if(g_d <= t){  # this confounding event has already occurred for the target cohort
        c <- c & (Gd_vals == g_d)
      } else {       # the target is not confounded by it, so the control must not be either
        c <- c & (Gd_vals >= min_conf_cohort)
      }
    }
    if(p$control_option == "notyet"){
      c[group_time[, is.infinite(G1)]] <- FALSE
    }
    cp <- group_time[, c & time == t]
    cb <- group_time[, c & time == base_period]

    #if any group have no available cohort, skip
    if(sum(tp) == 0 || sum(tb) == 0 || sum(cp) == 0 || sum(cb) == 0){return(NULL)}

    common <- intersect_control(group_time, cp, cb, gg, t, base_period)
    if(is.null(common)){return(NULL)}
    cp <- common$cp
    cb <- common$cb

    #assign the weights
    group_time[tp, det_weight := 1]
    group_time[tb, det_weight := -1]
    group_time[cp, sto_weight := -pg/sum(pg)]
    group_time[cb, sto_weight := pg/sum(pg)]

  }

  if(all(group_time[, det_weight+sto_weight] == 0)){return(NULL)} #not redundant!
  return(list(det = group_time[, det_weight], sto = group_time[, sto_weight]))

}

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
                    balanced = TRUE, seed = NA, stratify = FALSE, treatment_assign = "latent", second_cohort = FALSE, confound_ratio = 1, second_het = "all") {
  if (!is.na(seed)) {
    set.seed(seed)
  }

  # unit  -------------
  dt_i <- data.table(unit = 1:sample_size)
  if (cov == "int") {
    dt_i[, x := sample.int(5, sample_size, replace = TRUE)]
  } else if (cov == "no") {
    dt_i[, x := 1]
  } else if (cov == "cont") {
    dt_i[, x := rnorm(sample_size)]
  }

  if (second_cov) {
    dt_i[, x2 := rnorm(sample_size)]
  } else {
    dt_i[, x2 := 1]
  }

  if (stratify) {
    dt_i[, s := fifelse(rnorm(sample_size) > 0, 1, 2)]
  } else {
    dt_i[, s := 1]
  }

  # treatment assignment ---------------------

  # assign treated group based on a latent related to X

  if (treatment_assign == "latent") {
    ep1 <- rnorm(sample_size)
    dt_i[, treat_latent := x * 0.2 + x2 * 0.2 + ep1] # unit with larger X tend to be treated and treated earlier
    untreated_thres <- quantile(dt_i$treat_latent, untreated_prop)
    dt_i[treat_latent <= untreated_thres, G := Inf] # unit with low latent is never treated

    cohort_prop <- (1 - untreated_prop) / (time_period - 1)
    last_treat_thres <- untreated_thres
    for (t in time_period:2) { # unit start getting treated in t = 2
      treat_thres <- quantile(dt_i$treat_latent, untreated_prop + cohort_prop * (time_period - t + 1))
      dt_i[treat_latent <= treat_thres & treat_latent > last_treat_thres, G := t]
      last_treat_thres <- treat_thres
    }
    rm(t)
  } else if (treatment_assign == "uniform") {
    # when treatment is set to 'uniform', untreated propensity is fixed
    dt_i[, G := floor((unit - 1) / (sample_size / time_period))]
    dt_i[G < 2, G := Inf]
  }

  if (second_cohort) {
    setnames(dt_i, "G", "G2")
    if (treatment_assign == "latent") {
      dt_i[, treat_latent := x * 0.2 + x2 * 0.2 + ep1 * confound_ratio + rnorm(sample_size)] # unit with larger X tend to be treated and treated earlier
      untreated_thres <- quantile(dt_i$treat_latent, untreated_prop)
      dt_i[treat_latent <= untreated_thres, G := Inf] # unit with low latent is never treated

      cohort_prop <- (1 - untreated_prop) / (time_period - 1)
      last_treat_thres <- untreated_thres
      for (t in time_period:2) { # unit start getting treated in t = 2
        treat_thres <- quantile(dt_i$treat_latent, untreated_prop + cohort_prop * (time_period - t + 1))
        dt_i[treat_latent <= treat_thres & treat_latent > last_treat_thres, G := t]
        last_treat_thres <- treat_thres
      }
      rm(t)
    } else if (treatment_assign == "uniform") {
      # when treatment is set to 'uniform', untreated propensity is fixed
      dt_i[, G := floor((unit - 1) / (sample_size / time_period))]
      dt_i[G < 2, G := Inf]
    }
  }

  # assign unit FE
  dt_i[, unit_fe := rnorm(sample_size)]

  # time ------------------

  dt_t <- data.table(time = 1:time_period)
  dt_t[, time_fe := rnorm(time_period)]
  dt_t[, x_trend := rnorm(time_period)]

  # panel --------------------------
  dt <- CJ(unit = 1:sample_size, time = 1:time_period)
  dt <- dt |> merge(dt_i, by = "unit")
  dt <- dt |> merge(dt_t, by = "time")

  # add time_varying covariates
  if (vary_cov) {
    dt[, xvar := pmin(G, time_period + 4) * time^(1 / 3) * 0.1 + rnorm(sample_size * time_period, 0, 10)]
  } else {
    dt[, xvar := 1]
  }

  # untreated potential outcomes
  dt[, y0 := unit_fe + time_fe + (x + x2) * x_trend + xvar + rnorm(sample_size * time_period, sd = epsilon_size)]

  dt[, D := as.integer(time >= G)]
  # generate gtatt
  att <- CJ(G = 1:time_period, time = 1:time_period)
  if (hetero == "all") {
    att[, attgt := rnorm(time_period * time_period, mean = 2, sd = 1)]
  } else if (hetero == "dynamic") {
    for (event_t in 0:max(att[, time - G])) {
      att[time - G == event_t, attgt := rnorm(1, mean = 2, sd = 1)]
    }
  }
  att[time < G, attgt := 0] # no anticipation
  dt <- dt |> merge(att, by = c("time", "G"), all.x = TRUE, all.y = FALSE)
  dt[is.na(attgt), attgt := 0]
  dt[, tau := attgt * s]

  if (second_cohort) {
    dt[, D2 := as.integer(time >= G2)]
    # generate gtatt
    att2 <- CJ(G2 = 1:time_period, time = 1:time_period)
    if (second_het == "no") {
      att2[, attgt2 := 10]
    } else {
      if (hetero == "all") {
        att2[, attgt2 := rnorm(time_period * time_period, mean = 2, sd = 1)]
      } else if (hetero == "dynamic") {
        for (event_t in 0:max(att2[, time - G2])) {
          att2[time - G2 == event_t, attgt2 := rnorm(1, mean = 2, sd = 1)]
        }
      }
    }
    att2[time < G2, attgt2 := 0] # no anticipation

    # add att2 to att
    att[, event := 1]
    att2[, event := 2]
    att <- rbind(att, att2, fill = TRUE)

    dt <- dt |> merge(att2, by = c("time", "G2"), all.x = TRUE, all.y = FALSE)
    dt[is.na(attgt2), attgt2 := 0]
    dt[, tau2 := attgt2 * s]

    # potential outcome
    dt[, y10 := y0 + tau]
    dt[, y01 := y0 + tau2]
    dt[, y11 := y0 + tau + tau2]
    dt[, y := y0 * (1 - D) * (1 - D2) + y10 * D * (1 - D2) + y01 * (1 - D) * D2 + y11 * D * D2]
    cols <- c("time", "G", "G2", "unit", "x", "x2", "y", "s", "xvar")
  } else {
    # potential outcome
    dt[, y1 := y0 + tau]
    dt[, y := y1 * D + y0 * (1 - D)]
    cols <- c("time", "G", "unit", "x", "x2", "y", "s", "xvar")
  }

  dt <- dt[, .SD, .SDcols = cols]

  # additional -----------------

  if (na == "y") {
    dt[, y := na_insert(y)]
  } else if (na == "x") {
    dt[, x := na_insert(x)]
  } else if (na == "both") {
    dt[, y := na_insert(y)]
    dt[, x := na_insert(x)]
  }

  if (balanced == FALSE) {
    size <- fnrow(dt)
    dt <- dt[sample(1:size, size * 0.99)]
  }
  if (second_outcome == TRUE) {
    dt[, y2 := y + 1 + rnorm(fnrow(dt), 0, 0.1)]
  }

  return(list(dt = dt, att = att))
}

validate_argument <- function(dt, p) {
  if (!p$validate) {
    return(NULL)
  }

  dt_names <- names(dt)

  # release p
  for (name in names(p)) {
    assign(name, p[[name]])
  }

  name_message <- "__ARG__ must be a character scalar and a name of a column from the dataset."
  check_set_arg(timevar, unitvar, cohortvar, "match", .choices = dt_names, .message = name_message, .up = 1)

  covariate_message <- "__ARG__ must be NA or a character vector which are all names of columns from the dataset."
  check_set_arg(varycovariatesvar, covariatesvar, outcomevar,
    "NA | multi match",
    .choices = dt_names, .message = covariate_message, .up = 1
  )

  checkvar_message <- "__ARG__ must be NA or a character scalar if a name of columns from the dataset."
  check_set_arg(weightvar, clustervar, "NA | match", .choices = dt_names, .message = checkvar_message, .up = 1)

  check_set_arg(control_option, double_control_option, "match", .choices = c("both", "never", "notyet"), .up = 1) # kinda bad names since did's notyet include both notyet and never
  check_set_arg(control_type, "match", .choices = c("ipw", "reg", "dr"), .up = 1)
  check_set_arg(base_period, "match", .choices = c("varying", "universal"), .up = 1)
  check_arg(copy, validate, boot, allow_unbalance_panel, cband, parallel, add_base_period, "scalar logical", .up = 1)
  check_arg(anticipation, anticipation2, alpha, "scalar numeric", .up = 1)
  
  if (anticipation < 0) {
    stop("anticipation must be non-negative (>= 0), got: ", anticipation)
  }
  if (anticipation2 < 0) {
    stop("anticipation2 must be non-negative (>= 0), got: ", anticipation2)
  }

  if (!is.na(balanced_event_time)) {
    if (result_type != "dynamic") {
      stop("balanced_event_time is only meaningful with result_type == 'dynamic'")
    }
    check_arg(balanced_event_time, "numeric scalar", .up = 1)
    if (balanced_event_time < 0) {
      stop("balanced_event_time must be non-negative (>= 0), got: ", balanced_event_time)
    }
  }

  if (add_base_period == TRUE) {
    if (!result_type %in% c("dynamic", "dynamic_event")) {
      stop("add_base_period is only possible with result_type == 'dynamic' or 'dynamic_event'")
    }
  }

  # the effect model of the second stage
  check_set_arg(effect_model, "os formula | match", .choices = c("parallel", "unrestricted"), .up = 1)
  check_set_arg(effect_fit, "match", .choices = c("separate", "joint", "ordered"), .up = 1)
  if (effect_kind != "parallel") {
    if (allNA(cohortvar2)) {
      stop("effect_model needs multiple events: set cohortvar2.")
    }
    if (!event_specific) {
      stop("effect_model needs event_specific = TRUE.")
    }
    if (base_period != "universal") {
      stop("effect_model needs base_period = 'universal'. a varying base period makes the pre-period cells short differences.")
    }
    if (double_control_option != "both" || control_option == "notyet") {
      warning("effect_model does not restrict the fit set: double_control_option and control_option = 'notyet' do not apply to the second stage.")
    }
  }
  if (effect_fit != "separate") {
    if (effect_kind != "formula") {
      stop("effect_fit = '", effect_fit, "' needs a formula in effect_model.")
    }
  }
  if (effect_fit == "ordered" && anticipation != anticipation2) {
    stop("effect_fit = 'ordered' models one event kind, so anticipation and anticipation2 must be equal.")
  }
  if (result_type == "dynamic_event" && effect_fit != "ordered") {
    stop("result_type 'dynamic_event' needs effect_fit = 'ordered'.")
  }
  
  # Validate only_est_min / only_est_max
  has_est_min <- !is.na(exper$only_est_min)
  has_est_max <- !is.na(exper$only_est_max)
  if (has_est_min || has_est_max) {
    if (result_type != "dynamic") {
      stop("only_est_min/only_est_max can only be used with result_type == 'dynamic'")
    }
    if (!allNA(cohortvar2)) {
      stop("only_est_min/only_est_max cannot be used with double DiD (cohortvar2)")
    }
    if (has_est_min && (!is.numeric(exper$only_est_min) || length(exper$only_est_min) != 1))
      stop("only_est_min must be a numeric scalar")
    if (has_est_max && (!is.numeric(exper$only_est_max) || length(exper$only_est_max) != 1))
      stop("only_est_max must be a numeric scalar")
    if (has_est_min && has_est_max && exper$only_est_min > exper$only_est_max) {
      stop("only_est_min (", exper$only_est_min, ") must be <= only_est_max (", exper$only_est_max, ")")
    }
  }

  # Validate result_type for double DiD
  check_set_arg(result_type, "match", .choices = c("group_time", "time", "group", "simple", "dynamic", "group_group_time", "dynamic_stagger", "dynamic_event"), .up = 1)
  if (result_type %in% c("group_group_time", "dynamic_stagger", "dynamic_event")) {
    if (allNA(cohortvar2)) {
      stop("result_type '", result_type, "' can only be used with double DiD (cohortvar2 must be specified)")
    }
  }
  
  if (allow_unbalance_panel == TRUE && control_type == "dr") {
    stop("fastdid does not support DR when allowing for unbalanced panels.")
  }

  if(allow_unbalance_panel == TRUE & !allNA(varycovariatesvar)){
     stop("fastdid currently only supports time varying covariates when not allowing for unbalanced panels.")
  }
  
  if (any(covariatesvar %in% varycovariatesvar) && !allNA(varycovariatesvar) && !allNA(covariatesvar)) {
    stop("time-varying var and invariant var have overlaps.")
  }
  if (!boot && (!allNA(clustervar) || cband == TRUE)) {
    stop("clustering and uniform confidence interval only available with bootstrap")
  }

  if (parallel) {
    if (.Platform$OS.type != "unix") {
      stop("parallel option only available on unix systems")
    }
    if (!requireNamespace("parallel")) {
      stop("parallel requires the parallel package")
    }
  }

  # varname collision
  varnames <- unlist(p[str_subset(names(p), "var")])
  varnames <- varnames[!is.na(varnames)]
  if (any(duplicated(varnames))) {
    stop("-var arguments can not have duplicated names. (no need to specify cluster on unit-level, it is automatically done.)")
  }
}

validate_dt <- function(dt, p) {
  varnames <- unlist(p[str_ends(names(p), "var")], recursive = TRUE) # get all the argument that ends with "var"
  varnames <- varnames[!varnames %in% c(p$timevar, p$unitvar, p$cohortvar) & !is.na(varnames) & !is.null(varnames)]

  # the confounding events of double did, already renamed to G2 ... GM
  gcol2 <- character(0)
  if (!allNA(p$cohortvar2)) {
    gcol2 <- paste0("G", seq(2L, 1L + length(p$cohortvar2)))
  }

  # screen the confounding cohorts: a missing value silently drops a unit from the
  # control sets, and a fractional value corrupts the cohort labels
  for (col in gcol2) {
    if (!dt[, is.numeric(get(col))]) {
      stop(col, " needs to be numeric.")
    }
    na_units <- dt[is.na(get(col)), unique(unit)]
    if (length(na_units) > 0) {
      warning(length(na_units), " units have a missing value in ", col, ". fastdid drops them.")
      dt <- dt[!unit %in% na_units]
    }
    frac <- dt[!is.infinite(get(col)) & get(col) %% 1 != 0, .N]
    if (frac > 0) {
      stop(col, " must be a whole number or Inf. ", frac, " observations are not.")
    }
  }
  if (length(gcol2) > 0 && nrow(dt) == 0) {
    stop("no observations remain after the confounding cohorts are screened.")
  }

  # the k-th occurrence of one event: the dates must increase, and Inf is a suffix
  if (p$effect_fit == "ordered") {
    gcols <- c("G", gcol2)
    for (k in seq_len(length(gcols) - 1)) {
      a <- dt[[gcols[k]]]
      b <- dt[[gcols[k + 1]]]
      ok <- is.infinite(b) | (is.finite(a) & b > a)
      if (any(!ok)) {
        stop("effect_fit = 'ordered' needs g1 < g2 < ... for every unit, with Inf only after every finite date. ",
             dt[!ok, uniqueN(unit)], " units break the order between event ", k, " and event ", k + 1, ".")
      }
    }
  }

  # change to int
  uniquecols <- c("G", "time", "unit", gcol2)
  for (col in uniquecols) {
    if (!dt[, is.numeric(get(col))]) {
      stop(col, " needs to be numeric.")
    }
    dt[!is.infinite(get(col)), c(col) := as.integer(get(col))] # yeah sometimes floating point can be annoying
  }

  raw_unit_size <- dt[, uniqueN(unit)]
  raw_time_size <- dt[, uniqueN(time)]

  # Validate balanced_event_time against actual data
  if (!is.na(p$balanced_event_time)) {
    # Early check: ensure balanced_event_time doesn't exceed max possible event time in data
    max_event_time <- dt[, max(time - G)]
    if (p$balanced_event_time > max_event_time) {
      stop("balanced_event_time (", p$balanced_event_time, 
           ") is larger than the maximum event time in the data (", max_event_time, "). ",
           "Please specify a value between 0 and ", max_event_time, ".")
    }
  }

  # doesn't allow missing value
  if (is.na(p$exper$only_balance_2by2) || !p$exper$only_balance_2by2) {
    for (col in varnames) {
      na_obs <- whichNA(dt[, get(col)])
      if (length(na_obs) != 0) {
        warning("missing values detected in ", col, ", removing ", length(na_obs), " observation.")
        # whichNA returns integer row indices; remove those rows with negative indexing
        dt <- dt[-na_obs]
      }
    }
  }

  if (!allNA(p$covariatesvar) && uniqueN(dt, by = c("unit", p$covariatesvar)) > raw_unit_size) {
    warning("some covariates is time-varying, fastdid only use the first observation for covariates.")
  }

  if (!is.na(p$weightvar) && dt[,.(uniqueN(get(p$weightvar))), by = "unit"] [V1>1,] |> nrow() > 0) {
    stop("weightvar is time-varying, fastdid does not support time-varying weights.")
  }

  if (!is.na(p$weightvar) && nrow(dt[get(p$weightvar) == 0, ] > 0)) {
    stop("some weights are zero.")
  }

  if (!allNA(p$covariatesvar) || !allNA(p$varycovariatesvar)) {
    for (cov in c(p$covariatesvar, p$varycovariatesvar)) {
      if (is.na(cov)) {
        next
      }
      # check covaraites is not constant
      # for time-invariant covariates only the first obs per unit is used, so check variation there;
      # for time-varying covariates variation can come from any observation (within or across units)
      if (cov %in% p$varycovariatesvar) {
        novar <- fnunique(dt[, get(cov)]) == 1
      } else {
        novar <- fnunique(dt[, get(cov)[1], by = "unit"][, V1]) == 1
      }
      if (novar) stop(cov, " have no variation")
      if (!(is.numeric(dt[, get(cov)]) || is.integer(dt[, get(cov)]))) {
        stop(cov, " is not numeric or integer, do not support fixed effects.")
      }
    }
  }

  # check balanced panel
  # check if any is dup
  if (anyDuplicated(dt[, .(unit, time)])) {
    dup_id <- dt[duplicated(dt[, .(unit, time)]), unique(unit)]
    stop(length(dup_id), " units is observed more than once in a period.")
  }

  # check if any is missing
  if (!p$allow_unbalance_panel) {
    unit_count <- dt[, .(count = .N), by = unit]
    
    if (any(unit_count[, count < raw_time_size])) {
      mis_unit <- unit_count[count < raw_time_size]
      warning(nrow(mis_unit), " units is missing in some periods, enforcing balanced panel by dropping them")
      dt <- dt[!unit %in% mis_unit[, unit]]
      
      # Validate we still have data after dropping
      if (nrow(dt) == 0) {
        stop("No observations remain after enforcing balanced panel. Consider setting allow_unbalance_panel = TRUE")
      }
    }
  }

  # drop always_treated units
  if (nrow(dt) > 0) {
    min_time <- dt[, min(time)]
    always_treated <- dt[G <= min_time, unique(unit)]
    if (length(always_treated) > 0) {
      warning(length(always_treated), " units is treated in the first period, dropping them")
      dt <- dt[!unit %in% always_treated]
    }
  }

  # for double did part: check all confounding event columns
  if (!allNA(p$cohortvar2) && nrow(dt) > 0) {
    min_time <- dt[, min(time)]
    M <- 1L + length(p$cohortvar2)
    for(d in 2L:M){
      Gd_col <- paste0("G", d)
      always_treated <- dt[get(Gd_col) <= min_time, unique(unit)]
      if (length(always_treated) > 0) {
        warning(length(always_treated), " units is treated in the first period by event ", d, ", dropping them")
        dt <- dt[!unit %in% always_treated]
      }
    }
  }

  # test is weights 


  return(dt)
}

