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
