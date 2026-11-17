# high level -------------------------------------------------------------------

#' The second stage and the aggregation of every outcome, one first stage each.
#' @noRd
aggregate_gt <- function(fs_list, p) {
  results <- lapply(fs_list, run_double_did, p)
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

  # get aggregation weights, one row per target in the order of the targets
  agg_weights <- matrix(0, length(targets), nrow(group_time))
  for (i in seq_along(targets)) {
    sel <- group_time[, which(target == targets[i] & used)]
    agg_weights[i, sel] <- group_time[sel, pg / sum(pg)]
  }
  group_time[, pg := NULL]
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

#' The influence of the estimated cohort shares in the weights.
#'
#' @param group the cells of the weights, in the order of `att`.
#' @param agg_weights the weights, one row per target.
#' @param fs the first stage: the shares `pg` by cohort, and the cohort labels.
#' @return a matrix of the cohorts times the targets: the influence of a unit is
#'   its weight times the entry of its cohort (see `lin_inf()`).
#' @noRd
get_weight_influence <- function(att, group, agg_weights, fs, p, by_period = FALSE) {
  group <- group |> merge(fs$shares, by = "G", sort = FALSE)

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

  n_cohort <- length(fs$cohorts)
  rows <- asplit(agg_weights, 1)
  if (!p$parallel) {
    inf_weights <- vapply(rows, function(x) {
      get_weight_influence_param(x, group, att, fs$cohorts, by_period)
    }, numeric(n_cohort))
  } else {
    inf_weights <- unlist(mclapply(rows, function(x) {
      get_weight_influence_param(x, group, att, fs$cohorts, by_period)
    }))
  }

  return(matrix(inf_weights, n_cohort, nrow(agg_weights)))
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
#' @param cohorts the cohort labels of the first stage.
#' @param by_period logical, normalize each period on its own.
#' @return the influence of a unit of weight 1 in each cohort.
#' @noRd
get_weight_influence_param <- function(agg_weights, group, gt_att, cohorts, by_period = FALSE) {
  n_cohort <- length(cohorts)
  keepers <- which(agg_weights != 0)
  if (length(keepers) == 0) {
    return(rep(0, n_cohort))
  } # for direct double did
  group <- group[keepers, ]
  signs <- sign(agg_weights[keepers])
  att_keep <- as.vector(gt_att)[keepers]

  # the share of a cell moves with the units of its cohort, so the influence of
  # a unit depends only on its cohort
  cell_cohort <- match(group[, G], cohorts)

  blocks <- if (by_period) group[, time] else rep(1L, nrow(group))

  # one normalized share for each block: d(pgi / sum(pgi))
  inf_weight <- rep(0, n_cohort)
  for (b in unique(blocks)) {
    idx <- which(blocks == b)
    pg_b <- group[idx, pg]
    v <- signs[idx] * att_keep[idx]
    # the sums over the cells of the block that belong to each cohort
    own_v <- tabulate_sum(cell_cohort[idx], v, n_cohort)
    own_n <- tabulate_sum(cell_cohort[idx], rep(1, length(idx)), n_cohort)
    pv <- sum(pg_b * v)
    inf_weight <- inf_weight + (own_v - pv) / sum(pg_b) - (own_n - sum(pg_b)) * pv / sum(pg_b)^2
  }

  inf_weight[abs(inf_weight) < sqrt(.Machine$double.eps) * 10] <- 0 # fill zero
  return(inf_weight)
}

#' The sum of `x` by an integer index in `1..n`, zero for an absent index.
#' @noRd
tabulate_sum <- function(index, x, n) {
  out <- numeric(n)
  agg <- rowsum(x, index)
  out[as.integer(rownames(agg))] <- agg[, 1]
  out
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
