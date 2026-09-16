

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
