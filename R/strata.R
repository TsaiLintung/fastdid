# the strata first stage --------------------------------------------------------

#' The units, the cohorts, and the strata of the strata path.
#'
#' `control_type = "strata"` computes every cell from the means of the outcome
#' by cohort, stratum, and period. A stratum is one value of `covariatesvar`
#' (one combination, for several columns), and without covariates every unit is
#' in one stratum. A group is a cohort row of `aux$cohort_sizes` and a stratum.
#' The cohort rows are sorted by the first event, so the control cohorts of a
#' cell are a contiguous range of rows, as in `get_control_pos()`.
#'
#' @return a list: the number of rows `R` and strata `S`, the first event of
#'   each row, the group of each unit (`r + (x - 1) R`), the group sizes as an
#'   `R x S` matrix, and their prefix sums over the rows.
#' @noRd
get_strata_setup <- function(aux, p) {
  cs <- aux$cohort_sizes
  R <- nrow(cs)
  unit_row <- match(aux$dt_inv[, G], cs[, G])
  if (allNA(p$covariatesvar)) {
    unit_x <- rep(1L, aux$id_size)
  } else {
    unit_x <- GRP(aux$covariates, sort = TRUE)$group.id
  }
  S <- max(unit_x)
  unit_k <- unit_row + (unit_x - 1L) * R
  n_rs <- matrix(tabulate(unit_k, R * S), R, S)
  strata <- if (allNA(p$covariatesvar)) NULL else GRP(aux$covariates, sort = TRUE)$groups
  list(R = R, S = S, first = ming(cs[, G]), unit_k = unit_k, n_rs = n_rs,
       pn = rbind(0, col_cumsum(n_rs)), strata = as.data.table(strata))
}

#' The cells of one outcome from the stratified means.
#'
#' A cell is the DiD within each stratum, averaged with the stratum shares of
#' the treated cohort,
#' \deqn{ATT(g, t) = \sum_x \pi_{gx} [\Delta\bar Y_{gx} - \Delta\bar Y_{Cx}],}
#' where `C` is the control set of `get_control_range()`. This equals
#' `control_type = "reg"` with one dummy for each stratum. The cells follow
#' the rules of `estimate_gtatt_outcome()`: the same base periods, the same
#' cells skipped, and the same order.
#'
#' @param y the outcome, a list of one vector per period.
#' @param sx the output of `get_strata_setup()`.
#' @return a list: `gt`, `att`, and `inf_func`, a `strata_inf` object that
#'   stores no influence function per cell (see `apply_inf()`).
#' @noRd
estimate_gtatt_outcome_strata <- function(y, sx, aux, p) {
  R <- sx$R
  S <- sx$S
  TT <- length(aux$time_periods)

  # the sums of the outcome by group and period, and their prefix sums over the rows
  sum_y <- vapply(seq_len(TT), function(tt) tabulate_sum(sx$unit_k, y[[tt]], R * S), numeric(R * S))
  sum_y <- matrix(sum_y, R * S, TT)
  py <- rbind(0, col_cumsum(matrix(sum_y, R, S * TT)))

  # every cohort and period, first loop t then g, as in estimate_gtatt_outcome()
  treated <- which(is.finite(sx$first))
  cells <- CJ(t = aux$time_periods, r = treated, sorted = FALSE)
  setorder(cells, t, r)
  cells[, first := sx$first[r]]
  cells[, b := get_base_period(first, t, p)]
  cells <- cells[t != b & b %in% aux$time_periods]
  event_time <- cells[, t - first]
  if (!is.na(p$exper$only_est_min)) cells <- cells[event_time >= p$exper$only_est_min]
  event_time <- cells[, t - first]
  if (!is.na(p$exper$only_est_max)) cells <- cells[event_time <= p$exper$only_est_max]

  # the control cohorts: rows ra to rb, less the treated cohort if it is in range
  range <- get_control_range(cells$first, cells$t, cells$b, aux, p)
  cells[, ra := findInterval(range$min, sx$first, left.open = TRUE) + 1L]
  cells[, rb := findInterval(range$max, sx$first)]
  cells[, excl := r >= ra & r <= rb]

  parts <- strata_cell_parts(cells, sx, sum_y, py)
  cells[, `:=`(att = parts$att, n_g = parts$n_g, n_c = parts$n_c, overlap = parts$overlap)]

  # the rules of the 2x2 path: no treated or no control unit, then a treated
  # stratum with no control unit (a coefficient of the regression is NA), then
  # a group of fewer than 2 units
  cells <- cells[n_g > 0 & n_c > 0]
  if (any(!cells$overlap)) {
    warning(sum(!cells$overlap), " group-time(s) skipped: a stratum of the treated cohort has no control unit")
    cells <- cells[overlap == TRUE]
  }
  small <- cells[, n_g < 2 | n_c < 2]
  if (any(small)) {
    warning(sum(small), " group-time(s) skipped: a group has fewer than 2 effective units, so the variance is not estimable")
    cells <- cells[!small]
  }
  if (nrow(cells) == 0) stop("no valid group-times att to compute")

  G <- aux$cohort_sizes[, G][cells$r]
  gt <- data.table(G = G, time = as.integer(cells$t))
  att <- matrix(cells$att, ncol = 1, dimnames = list(paste0(G, ".", cells$t), NULL))
  inf_func <- structure(list(cells = cells[, .(r, t, b, ra, rb, excl, att, n_g, n_c)],
                             y = y, sx = sx, sum_y = sum_y, py = py, N = aux$id_size,
                             cohorts = aux$cohort_sizes[, G], cluster = aux$cluster),
                        class = "strata_inf")
  list(gt = gt, att = att, inf_func = inf_func)
}

#' The stratum parts of a set of cells.
#'
#' @param cells a table with the treated row `r`, the periods `t` and `b`, and
#'   the control rows `ra`, `rb`, `excl`.
#' @param keep_strata keep the `m x S` matrices, else only the vectors.
#' @return a list: the vectors `att`, `n_g`, `n_c`, and `overlap` (every
#'   treated stratum has a control unit), and with `keep_strata` the `m x S`
#'   matrices `n_gx` (treated units), `n_cx` (control units), `d` (the change of
#'   the treated), and `c` (the change of the controls).
#' @noRd
strata_cell_parts <- function(cells, sx, sum_y, py, keep_strata = FALSE) {
  # chunks of about a million stratum cells keep the index vectors small
  size <- max(1L, floor(1e6 / sx$S))
  chunks <- split(seq_len(nrow(cells)), ceiling(seq_len(nrow(cells)) / size))
  out <- lapply(chunks, function(i) strata_cell_parts_chunk(cells[i], sx, sum_y, py, keep_strata))
  vec <- c("att", "n_g", "n_c", "overlap")
  res <- lapply(setNames(vec, vec), function(v) unlist(lapply(out, `[[`, v), use.names = FALSE))
  if (keep_strata) {
    mat <- c("n_gx", "n_cx", "d", "c")
    res <- c(res, lapply(setNames(mat, mat), function(v) do.call(rbind, lapply(out, `[[`, v))))
  }
  res
}

#' @noRd
strata_cell_parts_chunk <- function(cells, sx, sum_y, py, keep_strata) {
  R <- sx$R
  S <- sx$S
  m <- nrow(cells)
  xs <- rep(seq_len(S), each = m)
  r <- rep(cells$r, S)
  ra <- rep(cells$ra, S)
  rb <- rep(cells$rb, S)
  excl <- rep(cells$excl, S)
  k <- r + (xs - 1L) * R

  n_gx <- sx$n_rs[cbind(r, xs)]
  n_cx <- sx$pn[cbind(rb + 1L, xs)] - sx$pn[cbind(ra, xs)] - excl * n_gx
  # the sum of the controls at a period: the prefix sums of the range, less the treated
  csum <- function(tt) {
    tt <- rep(tt, S)
    col <- xs + (tt - 1L) * S
    py[cbind(rb + 1L, col)] - py[cbind(ra, col)] - excl * sum_y[cbind(k, tt)]
  }
  d <- (sum_y[cbind(k, rep(cells$t, S))] - sum_y[cbind(k, rep(cells$b, S))]) / n_gx
  cc <- (csum(cells$t) - csum(cells$b)) / n_cx

  n_gx <- matrix(n_gx, m, S)
  n_cx <- matrix(n_cx, m, S)
  d <- matrix(d, m, S)
  cc <- matrix(cc, m, S)
  n_g <- rowSums(n_gx)
  # a stratum with no treated unit has no weight in the cell
  gap <- ifelse(n_gx > 0, n_gx * (d - cc), 0)
  out <- list(att = rowSums(gap) / n_g, n_g = n_g, n_c = rowSums(n_cx),
              overlap = rowSums(n_gx > 0 & n_cx == 0) == 0)
  if (keep_strata) out <- c(out, list(n_gx = n_gx, n_cx = n_cx, d = d, c = cc))
  out
}

#' The cumulative sums of the columns of a matrix.
#' @noRd
col_cumsum <- function(x) {
  if (nrow(x) == 1L) return(x)
  apply(x, 2, cumsum)
}

# the influence function without a column per cell ------------------------------

#' The influence functions of linear combinations of the cells.
#'
#' The first stage returns `inf_func`, a matrix of the units times the cells,
#' or a `strata_inf` object on the strata path. The second stage and the
#' aggregation only multiply it on the right by a weight matrix, so the strata
#' path forms the influence function of each combination directly.
#'
#' @param inf_func the influence functions of the cells.
#' @param weights a matrix of the cells times the combinations.
#' @return a matrix of the units times the combinations.
#' @noRd
apply_inf <- function(inf_func, weights) {
  if (!inherits(inf_func, "strata_inf")) return(inf_func %*% weights)
  weights <- as.matrix(weights)
  out <- matrix(0, inf_func$N, ncol(weights))
  for (j in seq_len(ncol(weights))) {
    nz <- which(weights[, j] != 0)
    out[, j] <- strata_inf_target(inf_func, nz, weights[nz, j])
  }
  out
}

#' The number of cells of an influence function.
#' @noRd
n_cells <- function(inf_func) {
  if (inherits(inf_func, "strata_inf")) nrow(inf_func$cells) else ncol(inf_func)
}

#' The influence function of one linear combination on the strata path.
#'
#' For the cell of the treated cohort g, a unit i of the cohort in stratum x
#' has the influence `(N / n_g) infl_g (dY_i - dY_Cx - ATT)`, and a control unit
#' j in stratum x has `-N [pi_gx / n_Cx + (infl_C - 1) / n_C] (dY_j - dY_Cx)`.
#' These are the terms of the regression path of `estimate_did_bp()` with one
#' dummy for each stratum, with its small-group inflation `infl`. The influence
#' of a unit is then linear in its outcomes, with coefficients that depend only
#' on its group, so the function builds the coefficients by group and period,
#' and reads the outcomes once. The control coefficients cover a range of rows,
#' so they go through a difference array over the rows.
#'
#' @param obj a `strata_inf` object.
#' @param idx,w the cells of the combination and their weights.
#' @param wm an optional weight on the group means: a table with the group `k`,
#'   the period `tt`, and `w`.
#' @return the influence function, one value per unit.
#' @noRd
strata_inf_target <- function(obj, idx, w, wm = NULL) {
  sx <- obj$sx
  R <- sx$R
  S <- sx$S
  K <- R * S
  TT <- ncol(obj$sum_y)
  N <- obj$N
  A <- numeric(K * TT)       # the coefficient on the outcome, by group and period
  B <- numeric(K)            # the constant, by group
  DA <- numeric((R + 1) * S * TT)
  DB <- numeric((R + 1) * S)

  if (length(idx) > 0) {
    cs <- obj$cells[idx]
    pr <- strata_cell_parts(cs, sx, obj$sum_y, obj$py, keep_strata = TRUE)
    m <- nrow(cs)
    xs <- rep(seq_len(S), each = m)
    k <- rep(cs$r, S) + (xs - 1L) * R
    t_ <- rep(cs$t, S)
    b_ <- rep(cs$b, S)
    cx <- as.vector(pr$c)

    # the treated cohort, in each stratum with a treated unit
    a <- rep(w * N * sqrt(cs$n_g / (cs$n_g - 1)) / cs$n_g, S)
    has_g <- as.vector(pr$n_gx) > 0
    A <- add_at(A, c(k + (t_ - 1L) * K, k + (b_ - 1L) * K)[c(has_g, has_g)], c(a, -a)[c(has_g, has_g)])
    B <- add_at(B, k[has_g], (-a * (cx + rep(cs$att, S)))[has_g])

    # the control cohorts, rows ra to rb in each stratum with a control unit
    infl_c <- sqrt(cs$n_c / (cs$n_c - 1))
    q <- -rep(w, S) * N * (as.vector(pr$n_gx / cs$n_g / pr$n_cx) + rep((infl_c - 1) / cs$n_c, S))
    has_c <- as.vector(pr$n_cx) > 0
    q[!has_c] <- 0
    cx[!has_c] <- 0
    ra <- rep(cs$ra, S)
    rb1 <- rep(cs$rb, S) + 1L
    da <- function(row, tt) row + (xs - 1L) * (R + 1L) + (tt - 1L) * (R + 1L) * S
    DA <- add_at(DA, c(da(ra, t_), da(rb1, t_), da(ra, b_), da(rb1, b_)), c(q, -q, -q, q))
    db <- function(row) row + (xs - 1L) * (R + 1L)
    DB <- add_at(DB, c(db(ra), db(rb1)), c(-q * cx, q * cx))
    # the treated cohort is not its own control
    ex <- rep(cs$excl, S) & has_c
    A <- add_at(A, c(k + (t_ - 1L) * K, k + (b_ - 1L) * K)[c(ex, ex)], c(-q, q)[c(ex, ex)])
    B <- add_at(B, k[ex], (q * cx)[ex])
  }

  # a weight on the group means: (N / n_k) (Y_it - M_kt)
  if (!is.null(wm) && nrow(wm) > 0) {
    n_k <- as.vector(sx$n_rs)[wm$k]
    v <- wm$w * N / n_k
    A <- add_at(A, wm$k + (wm$tt - 1L) * K, v)
    B <- add_at(B, wm$k, -v * obj$sum_y[cbind(wm$k, wm$tt)] / n_k)
  }

  A <- A + as.vector(col_cumsum(matrix(DA, R + 1L, S * TT))[seq_len(R), , drop = FALSE])
  B <- B + as.vector(col_cumsum(matrix(DB, R + 1L, S))[seq_len(R), , drop = FALSE])

  # read the outcomes once, in the periods with a coefficient
  A <- matrix(A, K, TT)
  psi <- B[sx$unit_k]
  for (tt in which(colSums(A != 0) > 0)) {
    psi <- psi + A[sx$unit_k, tt] * obj$y[[tt]]
  }
  psi
}

#' Add the values `val` at the positions `idx` of `vec`, duplicates summed.
#' @noRd
add_at <- function(vec, idx, val) {
  if (length(idx) == 0) return(vec)
  vec + tabulate_sum(idx, val, length(vec))
}

# linear combinations of the cells ----------------------------------------------

#' Linear combinations of the cells of the strata path
#'
#' Estimates targets that are linear in the cells of
#' `fastdid(control_type = "strata", full = TRUE)`, and in the means of the
#' outcome by cohort, stratum, and period, with their influence functions.
#' The strata path stores no influence function per cell, so a target with
#' its own weights goes through this function instead of `agg_inf_func`. The
#' influence function of a target is the same combination of the influence
#' functions of the cells that `fastdid()` aggregates. The weights are fixed.
#'
#' @param result the output of `fastdid(control_type = "strata", full = TRUE)`.
#' @param w_cell a data.table with the columns `target`, the cohort columns of
#'   the call (`cohortvar` and `cohortvar2`), the period column (`timevar`),
#'   and `w`, all on the scale of the data. A weight on a cell before the first
#'   event of its cohort has no effect: that cell is zero by no anticipation.
#' @param w_mean an optional data.table with the same columns, and the columns
#'   of `covariatesvar` if the call has strata: the weight of each target on the
#'   mean of a cohort and a stratum at a period.
#' @param outcome the outcome, a name or a position.
#' @param return_if if `TRUE`, also return the influence functions, one column
#'   per target, with the convention `estimate - truth = mean(psi)`.
#' @return a data.table with `target`, `estimate`, and `se`; with `return_if`,
#'   a list of that table and the matrix `inf_func`.
#' @examples
#' dt <- sim_did(1e+03, 5, seed = 1)$dt
#' res <- fastdid(dt, "time", "G", "unit", "y", control_type = "strata", full = TRUE)
#' w <- data.table::data.table(target = "sum", G = 3, time = 3:4, w = 1)
#' combine_cells(res, w)
#' @importFrom collapse GRP
#' @export
combine_cells <- function(result, w_cell, w_mean = NULL, outcome = 1, return_if = FALSE) {
  obj <- if (inherits(result, "fastdid_result")) result$gt_estimate[[outcome]]$inf_func else NULL
  if (!inherits(obj, "strata_inf")) {
    stop("combine_cells needs the output of fastdid(control_type = \"strata\", full = TRUE).")
  }
  p <- result$call
  ti <- result$time_info
  sx <- obj$sx
  evs <- c(p$cohortvar, if (!allNA(p$cohortvar2)) p$cohortvar2)
  cols <- c("target", evs, p$timevar, "w")
  norm <- function(x) (x - 1 - ti$time_offset) / ti$time_step + 1

  # the cohort rows on the scale of the data
  keys <- data.table(r = seq_len(sx$R))
  for (d in seq_along(evs)) keys[, (evs[d]) := recover_time(gd(obj$cohorts, d), ti)]
  map <- function(w, extra = character(0)) {
    miss <- setdiff(c(cols, extra), names(w))
    if (length(miss)) stop("the weights need the columns ", paste(miss, collapse = ", "), ".")
    w <- merge(as.data.table(w)[, c(cols, extra), with = FALSE], keys, by = evs, all.x = TRUE)
    if (anyNA(w$r)) stop(sum(is.na(w$r)), " weights name a cohort that is not in the data.")
    w[, tt := as.integer(round(norm(get(p$timevar))))]
    w
  }

  wc <- map(w_cell)
  targets <- unique(c(wc$target, if (!is.null(w_mean)) w_mean$target))
  wc <- wc[tt >= sx$first[r]] # the cells before the first event are zero
  wc[obj$cells[, .(r, t, cell = .I)], cell := i.cell, on = c(r = "r", tt = "t")]
  if (anyNA(wc$cell)) stop(sum(is.na(wc$cell)), " weights name a cell that fastdid() did not estimate.")

  if (!is.null(w_mean)) {
    xcols <- if (allNA(p$covariatesvar)) character(0) else p$covariatesvar
    wm <- map(w_mean, xcols)
    if (length(xcols)) {
      xkey <- copy(sx$strata)[, x := .I]
      wm <- merge(wm, xkey, by = xcols, all.x = TRUE)
      if (anyNA(wm$x)) stop(sum(is.na(wm$x)), " mean weights name a stratum that is not in the data.")
    } else {
      wm[, x := 1L]
    }
    wm[, k := r + (x - 1L) * sx$R]
  }

  psi <- matrix(0, obj$N, length(targets), dimnames = list(NULL, targets))
  est <- numeric(length(targets))
  for (j in seq_along(targets)) {
    x <- wc[target == targets[j], .(w = sum(w)), by = cell]
    y <- if (!is.null(w_mean)) wm[target == targets[j], .(w = sum(w)), by = .(k, tt)] else NULL
    est[j] <- sum(x$w * obj$cells$att[x$cell]) +
      if (is.null(y)) 0 else sum(y$w * obj$sum_y[cbind(y$k, y$tt)] / as.vector(sx$n_rs)[y$k])
    psi[, j] <- strata_inf_target(obj, x$cell, x$w, y)
  }
  se <- get_se(psi, list(cluster = obj$cluster), p)$se
  res <- data.table(target = targets, estimate = est, se = se)
  if (return_if) list(estimates = res, inf_func = psi) else res
}
