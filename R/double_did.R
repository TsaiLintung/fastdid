

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

# aggregation scheme -----------------------------------------------------------

#' The scheme for the event-specific effect.
#' @noRd
get_es_scheme <- function(group_time, aux, p){

  es_group_time <- copy(group_time) #group_time with available es effect
  #create lookup (columns already populated by get_agg_sch)
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
