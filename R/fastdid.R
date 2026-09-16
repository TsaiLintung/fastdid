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
#' @param effect_fit character, how a formula is fit. `"separate"` models one block (the target effect, or the confounding bundle) and leaves the other unmodeled, so the fit set is that block's clean cells; both orders are tried. `"joint"` models the sequential marginal effect of every event, in the order of the cohort columns, and fits on every cell; the cell is the sum of the marginal effects of its active events, so no additivity is assumed. `"ordered"` is the same fit for the k-th occurrence of one event kind, and checks that the cohort dates increase. `"state"` models the cell as a function of its state (the number of active events `nact`, its parity `status`, and the time since the last event `e`), and reports the marginal effect of event k as the state after k events minus the state after k-1 events; the dates must increase.
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
#'     \item{`nact`, `status`}{in the state fit: the number of active events, and its parity (1 for the "on" state of an on-and-off treatment)}
#'   }
#' The model `~ factor(gvec) + factor(gactive):factor(t)` is the parallel treatment effects assumption of `"parallel"`, fit on every clean cell at once. `~ factor(gvec) + poly(e, 2)` is a quadratic profile in event time with a cohort level. With `effect_fit = "separate"`, the target component is modeled for cohorts treated before they are confounded, and the confounding component for cohorts confounded before they are treated; the reported effect is the pure effect in the first case and the interacted effect in the second. With `effect_fit = "joint"` or `"ordered"`, each event block models the marginal effect of that event given the events before it, the design row of a cell is the sum over its active events, and a formula can share parameters across events, for example `~ 0 + factor(gvec):factor(event) + factor(e)`. Use `~ 0 + ...` in a stacked fit, because an intercept counts once per active event. An on-and-off treatment is the ordered case with alternating switches: `effect_fit = "state"` with `~ 0 + factor(status):factor(e)` gives one profile in time since the switch for the on state and one for the off state, and the reported effect of switch k is the change of state that it causes. `result_type = "dynamic_event"` averages the effect of every event at each event time, with cohort-size weights. The weights of the fit are treated as fixed in the influence function. `double_control_option` and `control_option = "notyet"` do not restrict the fit set of a formula. With `full = TRUE`, `effect_diag` returns the weight of every first-stage cell in every reported cell, the estimability of every candidate cell, the fit residuals, and an over-identification test.
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
