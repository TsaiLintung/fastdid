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
