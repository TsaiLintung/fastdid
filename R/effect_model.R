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
    long <- switch(p$effect_fit,
                   separate = build_component_rows(tab, p),
                   state = build_state_rows(tab, p),
                   build_event_rows(tab, p))
    long[, lid := .I] # a cell has one row per modeled event, so the cell index is not unique
    if (p$effect_fit == "state") { # the state before event k is the row of event k - 1
      long[, lid_prev := lid[match(paste(cell, event - 1L), paste(cell, event))]]
    }
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
      if (comp == "state") { # the marginal effect of event k is the state after it minus the state before it
        prev <- match(tr[, lid_prev], lr[, lid])
        Xp <- Xc[ifelse(is.na(prev), 1L, prev), , drop = FALSE]
        Xp[is.na(prev), ] <- 0
        Xt <- Xt - Xp
      }

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
  rep_rows[, prio := match(component, c("direct", "target", "confound", "joint", "state"))]
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

#' One row per cell and number of active events, for the state fit.
#'
#' The cell is a function of its state: the number of active events `nact`,
#' the parity `status` (1 when an odd number of events is active, the "on"
#' state of an on-and-off treatment), and the time since the last active
#' event `e`. The row of event `k` describes the state after the first `k`
#' events. The actual state of the cell is its fit row, and the marginal
#' effect of event `k` is the difference between the states after `k` and
#' after `k - 1` events.
#' @noRd
build_state_rows <- function(tab, p) {
  M <- 1L + length(p$cohortvar2)
  a <- p$anticipation
  nact <- Reduce(`+`, lapply(seq_len(M), function(d) as.integer(tab[, t >= get(paste0("G", d)) - a])))
  tab[, nact_cell := nact]
  rows <- rbindlist(lapply(seq_len(M), function(k) {
    r <- tab[nact_cell >= k]
    if (nrow(r) == 0) return(NULL)
    r[, `:=`(component = "state", event = as.integer(k), nact = as.integer(k), status = as.integer(k %% 2L),
             gown = as.character(get(paste0("G", k))), e = get(paste0("e", k)))]
    r[, gactive := gown]
    r[, ghist := do.call(paste, c(lapply(seq_len(k), function(j) as.character(get(paste0("G", j)))), list(sep = "-")))]
    r
  }))
  rows[, fit := event == nact_cell]
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
         "gvec, t, event, e, gown, gactive, ghist, g1..gM, e1..eM, and status, nact in the state fit.")
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
