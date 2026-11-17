# the first stage ----------------------------------------------------------------

#' Build the first stage of one outcome from the cells of `estimate_gtatt()`.
#'
#' The first stage is the only input of the second stage and the aggregation:
#' the cell table, the cell estimates, the cohort shares, and the uncertainty
#' of the estimates. The uncertainty is one of four classes:
#' - `unit_inf`: the influence functions of the units, one column per cell
#' - `strata_inf`: the strata path, which forms the influence function of a
#'   combination from the means (`strata.R`)
#' - `strata_vcov`: the strata path with a factor of the group moments in
#'   place of the units (`first_stage()`)
#' - a matrix: the covariance of the estimates (`as_first_stage()`)
#'
#' @param gt_result the first-stage result of one outcome: `gt`, `att`, `inf_func`.
#' @param t the time information of `coerce_dt()`.
#' @return a `fastdid_first_stage` object.
#' @noRd
new_first_stage <- function(gt_result, aux, p, t) {
  cells <- get_cell_table(gt_result, aux, p)
  cohorts <- aux$cohort_sizes[, G]
  shares <- data.table(w = aux$weights / sum(aux$weights), G = aux$dt_inv[, G])[, .(pg = sum(w)), by = "G"]
  unc <- gt_result$inf_func
  if (!inherits(unc, "strata_inf")) {
    unc <- structure(list(inf = unc, cohort = match(aux$dt_inv[, G], cohorts), w = aux$weights,
                          cluster = aux$cluster), class = "unit_inf")
  }
  evs <- c(p$cohortvar, if (!allNA(p$cohortvar2)) p$cohortvar2)
  xcols <- if (allNA(p$covariatesvar)) character(0) else p$covariatesvar
  # lab is the cohort label of fastdid, apart from the cohort columns of the data
  est <- data.table(kind = "cell", id = seq_len(nrow(cells)), lab = cells$G, tt = as.integer(cells$time),
                    estimate = as.vector(gt_result$att))
  if (inherits(unc, "strata_inf")) {
    # the means of each group: a cohort row r and a stratum x
    sx <- unc$sx
    n_k <- as.vector(sx$n_rs)
    means <- CJ(k = seq_len(sx$R * sx$S), tt = seq_len(ncol(unc$sum_y)))[n_k[k] > 0]
    means[, `:=`(kind = "mean", id = k, lab = cohorts[(k - 1L) %% sx$R + 1L],
                 estimate = unc$sum_y[cbind(k, tt)] / n_k[k])]
    if (length(xcols)) means <- cbind(means, sx$strata[(means$k - 1L) %/% sx$R + 1L])
    est <- rbind(est, means[, k := NULL], use.names = TRUE, fill = TRUE)
  }
  est <- cbind(est, data_keys(est$lab, est$tt, evs, p$timevar, t))
  fs <- list(cells = cells, att = gt_result$att, unc = unc, shares = shares, cohorts = cohorts,
             estimates = est, call = p, time_info = t, outname = gt_result$outname,
             cohortvar = evs, timevar = p$timevar, covariatesvar = xcols)
  return(structure(fs, class = "fastdid_first_stage"))
}

#' The cohort columns and the period on the scale of the data.
#' @noRd
data_keys <- function(G, tt, evs, timevar, t) {
  out <- data.table(tt = recover_time(rep_len(tt, length(G)), t))
  setnames(out, "tt", timevar)
  for (d in seq_along(evs)) set(out, j = evs[d], value = recover_time(gd(G, d), t))
  return(out)
}

#' The first stage of a fastdid result
#'
#' Returns the first stage of `fastdid(full = TRUE)`: the cells, their
#' uncertainty, and the cohort shares. `double_did()` runs the second stage
#' and the aggregation on it. On the strata path the first stage also holds
#' the means of the outcome by cohort, stratum, and period, and its
#' uncertainty is a factor of the moments within each group, so it holds no
#' unit data. On the other paths it holds the influence functions of the
#' units.
#'
#' @param result the output of `fastdid(full = TRUE)`.
#' @param outcome the outcome, a name or a position.
#' @return a `fastdid_first_stage` object. Its `estimates` table has `kind`
#'   ("cell" or "mean"), the cohort columns, the period column, the covariate
#'   columns (for a mean), and `estimate`.
#' @examples
#' dt <- sim_did(1e+03, 5, seed = 1)$dt
#' res <- fastdid(dt, "time", "G", "unit", "y", control_type = "strata", full = TRUE)
#' fs <- first_stage(res)
#' @export
first_stage <- function(result, outcome = 1) {
  if (!inherits(result, "fastdid_result")) stop("first_stage needs the output of fastdid(full = TRUE).")
  fs <- result$first_stage[[outcome]]
  if (inherits(fs$unc, "strata_inf")) fs$unc <- strata_factor(fs$unc)
  return(fs)
}

#' A first stage from the cells of another estimator
#'
#' Builds a first stage from cell estimates and their uncertainty, so that
#' `double_did()` runs the second stage of fastdid on the cells of another
#' first-stage estimator. Give the uncertainty either as a covariance matrix
#' or as the influence functions of the units. Without the cohort of each
#' unit, the cohort shares are fixed in the standard errors.
#'
#' @param cells a data.table with the cohort columns, the period column, and
#'   `estimate`: one row per cell.
#' @param shares a data.table with the cohort columns and `share`: the share
#'   of each cohort in the population, for every cohort of `cells`.
#' @param cohortvar the names of the cohort columns: the event of interest
#'   first, then the confounding events.
#' @param timevar the name of the period column.
#' @param vcov the covariance matrix of the estimates, in the order of the
#'   rows of `cells`.
#' @param inf_func the influence functions of the units, one column per row
#'   of `cells`, with the convention `estimate - truth = mean(inf_func)`.
#' @param unit_cohort optional, with `inf_func`: a data.table with the cohort
#'   columns of each unit, in the order of the rows of `inf_func`. It gives
#'   the influence of the estimated shares.
#' @return a `fastdid_first_stage` object.
#' @export
as_first_stage <- function(cells, shares, cohortvar, timevar, vcov = NULL, inf_func = NULL,
                           unit_cohort = NULL) {
  if (is.null(vcov) == is.null(inf_func)) stop("give one of vcov and inf_func.")
  cells <- as.data.table(cells)
  shares <- as.data.table(shares)
  miss <- setdiff(c(cohortvar, timevar, "estimate"), names(cells))
  if (length(miss)) stop("cells needs the columns ", paste(miss, collapse = ", "), ".")
  miss <- setdiff(c(cohortvar, "share"), names(shares))
  if (length(miss)) stop("shares needs the columns ", paste(miss, collapse = ", "), ".")
  C <- nrow(cells)
  if (!is.null(vcov) && !all(dim(vcov) == C)) stop("vcov must be a square matrix with one row per cell.")
  if (!is.null(inf_func) && ncol(inf_func) != C) stop("inf_func must have one column per cell.")

  # the periods on the scale 1, 2, ... of coerce_dt()
  periods <- sort(unique(cells[[timevar]]))
  step <- if (length(periods) > 1) min(diff(periods)) else 1
  t <- list(time_offset = min(periods) - 1, time_step = step)
  norm <- function(x) (x - 1 - t$time_offset) / t$time_step + 1
  label <- function(x) {
    g <- lapply(cohortvar, \(v) norm(x[[v]]))
    if (length(g) == 1) return(g[[1]])
    return(do.call(paste, c(g, list(sep = "-"))))
  }

  M <- length(cohortvar)
  tab <- data.table(G = label(cells), time = as.integer(round(norm(cells[[timevar]]))), pos = seq_len(C))
  for (d in seq_len(M)) set(tab, j = paste0("G", d), value = norm(cells[[cohortvar[d]]]))
  tab[, mg := Reduce(pmin, .SD), .SDcols = paste0("G", seq_len(M))]
  sh <- data.table(G = label(shares), pg = shares$share)
  if (!all(tab$G %in% sh$G)) stop("shares has no share for some cohorts of cells.")
  tab <- merge(tab, sh, by = "G", sort = FALSE)
  setorderv(tab, c("time", "mg", paste0("G", seq_len(M))))
  ord <- tab$pos
  att <- matrix(cells$estimate[ord], ncol = 1, dimnames = list(paste0(tab$G, ".", tab$time), NULL))

  cohorts <- sh$G
  if (!is.null(vcov)) {
    unc <- as.matrix(vcov)[ord, ord, drop = FALSE]
  } else {
    unit <- if (is.null(unit_cohort)) NULL else match(label(as.data.table(unit_cohort)), cohorts)
    if (anyNA(unit)) stop("unit_cohort names a cohort that is not in shares.")
    unc <- structure(list(inf = as.matrix(inf_func)[, ord, drop = FALSE], cohort = unit,
                          w = rep(1, nrow(inf_func)), cluster = NULL), class = "unit_inf")
  }
  p <- default_params(cohortvar, timevar)
  est <- data.table(kind = "cell", id = seq_len(C), lab = tab$G, tt = tab$time, estimate = as.vector(att))
  est <- cbind(est, data_keys(est$lab, est$tt, cohortvar, timevar, t))
  fs <- list(cells = tab[, c("G", "time", "pg", "mg", paste0("G", seq_len(M))), with = FALSE],
             att = att, unc = unc, shares = sh, cohorts = cohorts, estimates = est, call = p,
             time_info = t, outname = "estimate", cohortvar = cohortvar, timevar = timevar,
             covariatesvar = character(0))
  return(structure(fs, class = "fastdid_first_stage"))
}

#' The dense covariance matrix of the cells of a first stage
#'
#' @param object a `fastdid_first_stage` object.
#' @param ... unused.
#' @return a matrix with one row and one column per cell.
#' @export
vcov.fastdid_first_stage <- function(object, ...) {
  nc <- nrow(object$cells)
  return(lin_vcov(object, list(W = diag(nc))))
}

# the parameters of the second stage ---------------------------------------------

#' The parameters of a first stage that no fastdid() call made: the defaults of fastdid().
#' @noRd
default_params <- function(cohortvar, timevar) {
  f <- formals(fastdid)
  p <- lapply(f, \(x) if (is.symbol(x) && !nzchar(as.character(x))) NULL else eval(x))
  p$cohortvar <- cohortvar[1]
  p$cohortvar2 <- if (length(cohortvar) > 1) cohortvar[-1] else NA
  p$timevar <- timevar
  return(second_stage_params(p, list()))
}

#' The parameters of the first stage, with the arguments of the second stage in place.
#' @noRd
second_stage_params <- function(base, args) {
  allowed <- c("result_type", "balanced_event_time", "event_specific", "effect_model", "effect_fit",
               "double_control_option", "anticipation2", "control_option", "exper", "boot",
               "biters", "cband", "alpha", "parallel")
  bad <- setdiff(names(args), allowed)
  if (length(bad)) stop("double_did() does not take the arguments ", paste(bad, collapse = ", "), ".")
  p <- unclass(base)
  p[names(args)] <- args
  if (is.character(p$effect_model)) p$effect_model <- match.arg(p$effect_model, c("parallel", "unrestricted"))
  p$effect_kind <- if (is.character(p$effect_model)) p$effect_model else "formula"
  p$exper <- get_exper_default(p$exper, exper_arg_names())
  class(p) <- "locked"
  return(p)
}

# the second stage on the first stage ----------------------------------------------

#' The second stage and the targets of a first stage
#'
#' Runs the second stage of fastdid on a first stage, then the targets: the
#' aggregation of `result_type`, or the linear combinations of `w_cell` and
#' `w_mean`. `fastdid()` calls the same code on its own first stage.
#'
#' @param fs a `fastdid_first_stage` object, from `first_stage()` or
#'   `as_first_stage()`.
#' @param ... the arguments of `fastdid()` for the second stage and the
#'   targets: `result_type`, `balanced_event_time`, `event_specific`,
#'   `effect_model`, `effect_fit`, `double_control_option`, `anticipation2`,
#'   `control_option`, `exper`, `boot`, `biters`, `cband`, `alpha`, and
#'   `parallel`. An argument that is not given keeps the value of the first
#'   stage.
#' @param w_cell an optional data.table of targets with their own weights,
#'   with the columns `target`, the cohort columns, the period column, and
#'   `w`. The weights apply to the cells after the second stage. A weight on
#'   a cell before the first event of its cohort has no effect: that cell is
#'   zero by no anticipation.
#' @param w_mean an optional data.table with the same columns, and the
#'   covariate columns of the call: the weight of each target on a mean of
#'   the first stage.
#' @return with `w_cell` or `w_mean`, a data.table with `target`, `estimate`,
#'   and `se`, with the covariance matrix of the targets as the attribute
#'   `vcov`. Else the estimates of `result_type`, as `fastdid()` returns them.
#' @export
double_did <- function(fs, ..., w_cell = NULL, w_mean = NULL) {
  if (!inherits(fs, "fastdid_first_stage")) stop("double_did needs a first stage: see first_stage().")
  p <- second_stage_params(fs$call, list(...))
  out <- run_double_did(fs, p, w_cell, w_mean)
  if (!is.null(w_cell) || !is.null(w_mean)) return(out$result)
  return(convert_targets(out$result, p, fs$time_info))
}

#' The second stage, then the targets, on one first stage.
#'
#' @return a list: the result table, the influence functions of the targets
#'   (NULL without unit data), the aggregation weights, the second-stage
#'   weights, and the diagnostics of the effect model.
#' @noRd
run_double_did <- function(fs, p, w_cell = NULL, w_mean = NULL) {
  cells <- fs$cells
  att <- fs$att
  lin <- list()
  es_weight <- NULL
  effect_diag <- NULL
  if (p$event_specific && !allNA(p$cohortvar2)) {
    ss <- second_stage(fs, p)
    cells <- ss$cells # some gt may not have an identified effect (ex: g1 == g2)
    att <- ss$att
    es_weight <- ss$weight
    effect_diag <- ss$diag
    lin <- lin_compose(lin, ss$weight, ss$cc)
  }
  if (!is.null(w_cell) || !is.null(w_mean)) return(custom_targets(fs, cells, att, lin, p, w_cell, w_mean))

  agg_sch <- get_agg_sch(cells, p)
  agg_att <- agg_sch$agg_weights %*% att
  cc <- get_weight_influence(att, agg_sch$group_time, agg_sch$agg_weights, fs, p)
  lin <- lin_compose(lin, agg_sch$agg_weights, cc)
  s <- lin_se(fs, lin, p)

  result <- data.table(agg_sch$targets, agg_att, s$se)
  names(result) <- c("target", "att", "se")
  result[, `:=`(
    outcome = fs$outname,
    att_ciub = att + se * s$crit_val,
    att_cilb = att - se * s$crit_val
  )]
  return(list(result = result, inf_func = s$inf, weight_matrix = agg_sch$agg_weights,
              es_weight = es_weight, effect_diag = effect_diag))
}

#' The targets of user weights on the cells and on the means.
#' @noRd
custom_targets <- function(fs, cells, att, lin, p, w_cell, w_mean) {
  evs <- fs$cohortvar
  tv <- fs$timevar
  map <- function(w, keys, cols) {
    miss <- setdiff(c("target", cols, "w"), names(w))
    if (length(miss)) stop("the weights need the columns ", paste(miss, collapse = ", "), ".")
    w <- merge(as.data.table(w)[, c("target", cols, "w"), with = FALSE], keys, by = cols,
               all.x = TRUE, sort = FALSE)
    if (anyNA(w$pos)) stop(sum(is.na(w$pos)), " weights name an estimate that is not in the data.")
    return(w)
  }
  targets <- unique(c(w_cell$target, w_mean$target))
  J <- length(targets)
  dense <- function(w, n) {
    w <- w[, .(w = sum(w)), by = .(j = match(target, targets), pos)]
    out <- matrix(0, J, n)
    out[cbind(w$j, w$pos)] <- w$w
    return(out)
  }

  keys <- data_keys(cells$G, cells$time, evs, tv, fs$time_info)[, pos := .I]
  if (is.null(w_cell)) {
    wc <- keys[0, c(evs, tv), with = FALSE][, `:=`(target = character(0), w = numeric(0))]
  } else {
    wc <- as.data.table(w_cell)
  }
  # a cohort that is not in the data stops before the cells before the first event drop
  coh <- unique(data_keys(fs$cohorts, 1L, evs, tv, fs$time_info)[, evs, with = FALSE])[, ok := TRUE]
  bad <- merge(wc[, evs, with = FALSE], coh, by = evs, all.x = TRUE)[is.na(ok)]
  if (nrow(bad)) stop(nrow(bad), " weights name a cohort that is not in the data.")
  wc <- wc[wc[[tv]] >= do.call(pmin, unname(as.list(wc[, evs, with = FALSE])))]
  Wt <- dense(map(wc, keys, c(evs, tv)), nrow(cells))
  estimate <- as.vector(Wt %*% att)
  lin <- lin_compose(lin, Wt)

  if (!is.null(w_mean)) {
    means <- fs$estimates[fs$estimates$kind == "mean"]
    if (!nrow(means)) stop("the first stage has no means: only the strata path has them.")
    mkeys <- means[, c(evs, tv, fs$covariatesvar), with = FALSE][, pos := .I]
    Wm <- dense(map(w_mean, mkeys, c(evs, tv, fs$covariatesvar)), nrow(means))
    estimate <- estimate + as.vector(Wm %*% means$estimate)
    lin$Wm <- Wm
  }
  s <- lin_se(fs, lin, p)
  V <- s$vcov
  dimnames(V) <- list(targets, targets)
  out <- data.table(target = targets, estimate = estimate, se = s$se)
  setattr(out, "vcov", V)
  inf <- s$inf
  if (!is.null(inf)) colnames(inf) <- targets
  return(list(result = out, inf_func = inf))
}

# linear maps of the first stage ---------------------------------------------------

#' Compose a linear map of the first stage with the weights `W2`.
#'
#' A map gives each target as `W` times the first-stage cells, plus `Wm`
#' times the first-stage means, plus `cc`, the influence of the estimated
#' cohort shares, which is constant within a cohort (cohorts times targets).
#' An empty map is the identity on the cells. The influence function of the
#' targets is formed once, at the end, so no step builds a matrix of the
#' units times the cells.
#'
#' @param lin a map: a list with `W`, `Wm`, and `cc`, each NULL if absent.
#' @param W2 the weights of the new targets on the targets of `lin`.
#' @param cc2 the share influence of the new targets, or NULL.
#' @return the composed map.
#' @noRd
lin_compose <- function(lin, W2, cc2 = NULL) {
  out <- list(W = if (is.null(lin$W)) W2 else W2 %*% lin$W,
              Wm = if (is.null(lin$Wm)) NULL else W2 %*% lin$Wm,
              cc = if (is.null(lin$cc)) cc2 else lin$cc %*% t(W2))
  if (!is.null(lin$cc) && !is.null(cc2)) out$cc <- out$cc + cc2
  return(out)
}

#' The influence functions of the targets of a map, one column per target.
#' @noRd
lin_inf <- function(fs, lin) {
  unc <- fs$unc
  if (inherits(unc, "unit_inf")) {
    out <- if (is.null(lin$W)) unc$inf else unc$inf %*% t(lin$W)
    if (!is.null(lin$cc) && !is.null(unc$cohort)) out <- out + unc$w * lin$cc[unc$cohort, , drop = FALSE]
    return(out)
  }
  if (inherits(unc, "strata_inf")) {
    means <- fs$estimates[fs$estimates$kind == "mean"]
    J <- nrow(lin$W)
    out <- vapply(seq_len(J), \(j) strata_unit_eval(unc, strata_lin_coef(unc, lin, j, means)), numeric(unc$N))
    return(matrix(out, unc$N, J))
  }
  stop("the first stage holds no influence function of the units.")
}

#' The covariance matrix of the targets of a map.
#' @noRd
lin_vcov <- function(fs, lin) {
  unc <- fs$unc
  if (inherits(unc, c("unit_inf", "strata_inf"))) {
    inf <- lin_inf(fs, lin)
    return(crossprod(inf) / nrow(inf)^2)
  }
  if (inherits(unc, "strata_vcov")) {
    means <- fs$estimates[fs$estimates$kind == "mean"]
    thetas <- lapply(seq_len(nrow(lin$W)), function(j) {
      cf <- strata_lin_coef(unc, lin, j, means)
      return(cbind(cf$B, cf$A))
    })
    return(moment_quad(unc$L, thetas, unc$N))
  }
  # a covariance matrix holds no units, so the shares are fixed
  if (!is.null(lin$Wm)) stop("the first stage has no means.")
  return(lin$W %*% unc %*% t(lin$W))
}

#' The standard errors of the targets of a map.
#'
#' With the influence functions of the units, `get_se()` gives the standard
#' errors, with the bootstrap and the clusters. Without them, the standard
#' errors are the square roots of the diagonal of the covariance.
#'
#' @return a list: `se`, `crit_val`, `inf` (NULL without units), and `vcov`.
#' @noRd
lin_se <- function(fs, lin, p) {
  if (inherits(fs$unc, c("unit_inf", "strata_inf"))) {
    inf <- lin_inf(fs, lin)
    s <- get_se(inf, list(cluster = fs$unc$cluster), p)
    return(list(se = s$se, crit_val = s$crit_val, inf = inf, vcov = crossprod(inf) / nrow(inf)^2))
  }
  if (p$boot) stop("the bootstrap needs the influence functions of the units, which this first stage does not hold.")
  V <- lin_vcov(fs, lin)
  return(list(se = sqrt(diag(V)), crit_val = stats::qnorm(1 - p$alpha / 2), inf = NULL, vcov = V))
}

#' The coefficients of target j of a map on the strata path.
#'
#' @param obj a `strata_inf` or `strata_vcov` object.
#' @param means the mean rows of the estimates, with `id` (the group) and `tt`.
#' @return the coefficients of `strata_coef()`, with the share influence in `B`.
#' @noRd
strata_lin_coef <- function(obj, lin, j, means) {
  wc <- lin$W[j, ]
  idx <- which(wc != 0)
  wm <- NULL
  if (!is.null(lin$Wm)) {
    m <- which(lin$Wm[j, ] != 0)
    wm <- data.table(k = means$id[m], tt = means$tt[m], w = lin$Wm[j, m])
  }
  cf <- strata_coef(obj, idx, wc[idx], wm)
  # a group is a cohort row and a stratum, and the share influence is by cohort row
  if (!is.null(lin$cc)) cf$B <- cf$B + rep(lin$cc[, j], obj$sx$S)
  return(cf)
}

# the factored covariance of the strata path ------------------------------------------

#' The strata uncertainty without the units: a factor of the group moments.
#' @noRd
strata_factor <- function(obj) {
  out <- list(cells = obj$cells, sx = obj$sx[c("R", "S", "first", "n_rs", "pn")],
              sum_y = obj$sum_y, py = obj$py, N = obj$N, L = strata_moments(obj))
  return(structure(out, class = "strata_vcov"))
}

#' A factor of the moments of the outcome within each group of the strata path.
#'
#' For a group k, the moments are `Q_k = sum_i z_i z_i'` over its units, with
#' `z = (1, Y_1 - M_k1, ..., Y_T - M_kT)` and `M_k` the means of the group.
#' The factor stacks one block `F_k` per group with `F_k' F_k = Q_k`, from the
#' eigenvectors of the nonzero eigenvalues, so a group has at most
#' `min(n_k, T + 1)` rows.
#'
#' @param obj a `strata_inf` object.
#' @return a list: `F`, the stacked factor (`T + 1` columns), `k`, the group
#'   of each row, and `M`, the means by group and period.
#' @noRd
strata_moments <- function(obj) {
  K <- obj$sx$R * obj$sx$S
  uk <- obj$sx$unit_k
  # the uncentered moments lose the digits of the variance to the level of the outcome
  mean_y <- obj$sum_y / pmax(as.vector(obj$sx$n_rs), 1)
  z <- c(list(rep(1, obj$N)), lapply(seq_along(obj$y), \(t) obj$y[[t]] - mean_y[uk, t]))
  M <- length(z)
  Q <- array(0, c(K, M, M))
  for (a in seq_len(M)) {
    for (b in a:M) {
      q <- tabulate_sum(uk, z[[a]] * z[[b]], K)
      Q[, a, b] <- q
      Q[, b, a] <- q
    }
  }
  blocks <- lapply(which(as.vector(obj$sx$n_rs) > 0), function(k) {
    e <- eigen(Q[k, , ], symmetric = TRUE)
    keep <- e$values > max(e$values) * 1e-12
    f <- t(e$vectors[, keep, drop = FALSE] %*% diag(sqrt(e$values[keep]), sum(keep)))
    return(list(F = f, k = rep(k, nrow(f))))
  })
  return(list(F = do.call(rbind, lapply(blocks, `[[`, "F")),
              k = unlist(lapply(blocks, `[[`, "k")), M = mean_y))
}

#' The covariance of influence functions that are linear in z within each group.
#'
#' The influence of a unit of group k on target j is `theta_kj' z_i`, so the
#' sum over the units of the products of two targets is
#' `sum_k theta_kj' F_k' F_k theta_kl`. Each row of the factor acts as one
#' pseudo unit, and the covariance divides by `N^2`, as `get_se()` does.
#'
#' @param L the factor of `strata_moments()`.
#' @param thetas a list, one `K x (T + 1)` coefficient matrix for each target.
#' @return the covariance matrix of the targets.
#' @noRd
moment_quad <- function(L, thetas, N) {
  psi <- vapply(thetas, function(th) {
    # the constant absorbs the means, because z holds the centered outcome
    th[, 1] <- th[, 1] + rowSums(th[, -1, drop = FALSE] * L$M)
    return(rowSums(L$F * th[L$k, , drop = FALSE]))
  }, numeric(nrow(L$F)))
  psi <- matrix(psi, nrow(L$F))
  return(crossprod(psi) / N^2)
}
