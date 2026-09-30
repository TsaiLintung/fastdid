# cells without covariates ----------------------------------------------------

#' First-stage cells from cohort means, without covariates
#'
#' Computes the not-yet-treated first-stage cells on the cohort vector, for any
#' number of event columns, from the cohort-by-period means of the outcome. The
#' cohort of a unit is the vector of its event dates, and a cell is
#' \deqn{ATT(h, s) = [\bar Y_{h,s} - \bar Y_{h,m}] - [\bar Y_{c,s} - \bar Y_{c,m}],}
#' with the universal base period \eqn{m = \min h - 1} and the comparison
#' units \eqn{c} with no event by \eqn{s}. A cell before \eqn{\min h} is zero
#' by no anticipation. The function stores no influence function of a cell:
#' [combine_cells()] computes the influence function of a linear combination of
#' cells in closed form, so the memory does not grow with the number of cells.
#' This is the mode for a large panel with many cohorts and no covariates.
#'
#' @param data a balanced panel, a data.table or an object that coerces to one.
#' @param timevar,unitvar,outcomevar the column names of the period (consecutive
#'   integers), the unit, and the outcome.
#' @param cohortvars the column names of the event dates, one per event. A date
#'   is `Inf` if the unit never has the event. The dates must be constant within
#'   a unit.
#' @return an object of class `fastdid_cells`: `cohorts` (one row per cohort,
#'   with its dates, the first date `first`, and the size `n`), `cells` (the
#'   cohort row, the period `s`, and `att`), and the data that [combine_cells()]
#'   needs.
#' @examples
#' dt <- sim_did(1e+03, 5, seed = 1)$dt
#' fc <- fastdid_cells(dt, "time", "unit", "y", "G")
#' head(fc$cells)
#' @export
fastdid_cells <- function(data, timevar, unitvar, outcomevar, cohortvars) {
  check_arg(timevar, unitvar, outcomevar, "scalar charin", .choices = names(data))
  check_arg(cohortvars, "vector charin", .choices = names(data))
  dt <- as.data.table(data)[, c(unitvar, timevar, outcomevar, cohortvars), with = FALSE]
  setnames(dt, c(unitvar, timevar, outcomevar), c("unit", "time", "y"))
  if (anyNA(dt)) stop("fastdid_cells needs no missing value in the unit, the period, the outcome, or the event dates.")
  per <- sort(unique(dt$time))
  if (!identical(per, seq(min(per), max(per)))) stop("fastdid_cells needs consecutive integer periods.")
  if (dt[, .N, by = unit][, any(N != length(per))]) stop("fastdid_cells needs a balanced panel.")
  if (nrow(unique(dt[, c("unit", cohortvars), with = FALSE])) != uniqueN(dt$unit)) {
    stop("fastdid_cells needs event dates that are constant within a unit.")
  }

  # units, their cohort, and their first event
  setkey(dt, unit, time)
  units <- dt[time == min(per), c("unit", cohortvars), with = FALSE]
  units[, first := do.call(pmin, .SD), .SDcols = cohortvars]
  drop <- units[, first <= min(per)]
  if (any(drop)) {
    warning(sum(drop), " units have an event in the first period and have no base period. fastdid_cells drops them.")
    dt <- dt[unit %in% units[!drop, unit]]
    units <- units[!drop]
  }
  Y <- matrix(dt$y, ncol = length(per), byrow = TRUE) # units in key order, periods in order
  cohorts <- units[, .(n = .N), by = c(cohortvars, "first")]
  setorderv(cohorts, c("first", cohortvars))
  cohorts[, row := .I]
  units <- merge(units, cohorts[, c(cohortvars, "row"), with = FALSE], by = cohortvars, sort = FALSE)
  setkey(units, unit)

  # cohort-by-period means, and the shares of the comparison units at each period
  M <- unname(rowsum(Y, units$row, reorder = TRUE) / cohorts$n)
  S <- vapply(seq_along(per), function(k) {
    w <- cohorts$n * (cohorts$first > per[k])
    if (sum(w) == 0) rep(NA_real_, nrow(cohorts)) else w / sum(w)
  }, numeric(nrow(cohorts)))
  S <- matrix(S, nrow = nrow(cohorts))

  # the cells: the change of the cohort from its base minus the change of the comparison units
  col <- function(s) s - min(per) + 1
  cells <- cohorts[is.finite(first), .(row, first)][, .(s = per), by = .(row, first)]
  cells[, m := first - 1]
  ms <- unique(cells[s >= first, .(m, s)])
  ms[, ctrl := colSums(S[, col(s), drop = FALSE] * (M[, col(s), drop = FALSE] - M[, col(m), drop = FALSE]))]
  cells <- merge(cells, ms, by = c("m", "s"), all.x = TRUE)
  setorder(cells, row, s)
  cells[, att := 0]
  post <- cells[, s >= first]
  cells[post, att := M[cbind(row, col(s))] - M[cbind(row, col(m))] - ctrl]
  if (cells[post, anyNA(att)]) warning("some cells have no comparison unit: the last periods need never-treated units.")

  structure(list(cohorts = cohorts[], cells = cells[, .(row, s, att)], M = M, S = S, Y = Y,
                 unit_row = units$row, periods = per, cohortvars = cohortvars),
            class = "fastdid_cells")
}

#' Linear combinations of cells, with analytic standard errors
#'
#' Estimates targets that are linear in the cells of [fastdid_cells()], and in
#' the cohort-by-period means, with their influence functions. Each cell is a
#' linear function of the cohort means, so a target is too, and its influence
#' function is
#' \deqn{\psi_i = \frac{N}{n_{g_i}}\sum_s W_{g_i, s}\,(Y_{is} - \bar Y_{g_i, s}),}
#' where \eqn{W} is the weight of the target on the cohort means. The cohort
#' sizes are held fixed: the standard errors cover the sampling error of the
#' outcome means, not of the cohort shares.
#'
#' @param fc the output of [fastdid_cells()].
#' @param w_cell a data.table with the columns `target`, `row`, `s`, and `w`:
#'   the weight of each target on the cell of cohort `row` at period `s`. A
#'   weight on a cell before the first event of its cohort has no effect, since
#'   that cell is zero.
#' @param w_mean an optional data.table with the same columns: the weight of
#'   each target on the cohort mean of cohort `row` at period `s`.
#' @param return_if if `TRUE`, also return the influence functions, one column
#'   per target, with the convention `estimate - truth = mean(psi)`.
#' @return a data.table with `target`, `estimate`, and `se`; with `return_if`,
#'   a list of that table and the matrix `inf_func`.
#' @export
combine_cells <- function(fc, w_cell, w_mean = NULL, return_if = FALSE) {
  if (!inherits(fc, "fastdid_cells")) stop("fc must be the output of fastdid_cells().")
  per <- fc$periods
  col <- function(s) s - min(per) + 1
  R <- nrow(fc$cohorts)
  first <- fc$cohorts$first
  w_cell <- as.data.table(w_cell)[, .(target, row, s, w)]
  w_cell <- w_cell[s >= first[row]]                   # the cells before the first event are zero
  if (!is.null(w_mean)) w_mean <- as.data.table(w_mean)[, .(target, row, s, w)]
  targets <- unique(c(w_cell$target, if (!is.null(w_mean)) w_mean$target))
  N <- length(fc$unit_row)
  n_unit <- fc$cohorts$n[fc$unit_row]
  resid <- fc$Y - fc$M[fc$unit_row, , drop = FALSE]

  out <- vector("list", length(targets))
  psi_all <- if (return_if) matrix(NA_real_, N, length(targets), dimnames = list(NULL, targets)) else NULL
  for (k in seq_along(targets)) {
    W <- matrix(0, R, length(per))
    x <- w_cell[target == targets[k], .(w = sum(w)), by = .(row, s)]
    if (nrow(x) > 0) {
      x[, m := first[row] - 1]
      # the cohort change: + at s, - at the base m
      W <- W + as.matrix(Matrix_add(R, length(per), x$row, col(x$s), x$w))
      W <- W - as.matrix(Matrix_add(R, length(per), x$row, col(x$m), x$w))
      # the comparison change: - shares at s, + shares at m
      a <- x[, .(w = sum(w)), by = s]
      W[, col(a$s)] <- W[, col(a$s), drop = FALSE] - fc$S[, col(a$s), drop = FALSE] %*% diag(a$w, nrow(a))
      B <- x[, .(w = sum(w)), by = .(m, s)]
      for (j in seq_len(nrow(B))) {
        W[, col(B$m[j])] <- W[, col(B$m[j])] + B$w[j] * fc$S[, col(B$s[j])]
      }
    }
    if (!is.null(w_mean)) {
      y <- w_mean[target == targets[k], .(w = sum(w)), by = .(row, s)]
      if (nrow(y) > 0) W <- W + as.matrix(Matrix_add(R, length(per), y$row, col(y$s), y$w))
    }
    if (anyNA(W)) stop("target ", targets[k], " uses a cell with no comparison unit.")
    est <- sum(W * fc$M)
    psi <- rowSums(resid * W[fc$unit_row, , drop = FALSE]) * N / n_unit
    out[[k]] <- data.table(target = targets[k], estimate = est, se = sqrt(sum(psi^2)) / N)
    if (return_if) psi_all[, k] <- psi
  }
  res <- rbindlist(out)
  if (return_if) list(estimates = res, inf_func = psi_all) else res
}

#' A dense R x T matrix with the values `w` added at (`i`, `j`), duplicates summed.
#' @noRd
Matrix_add <- function(nr, nc, i, j, w) {
  z <- matrix(0, nr, nc)
  agg <- rowsum(w, (j - 1) * nr + i)
  z[as.integer(rownames(agg))] <- agg[, 1]
  z
}
