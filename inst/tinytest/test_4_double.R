tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+03, 5, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1,
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE, second_cohort = TRUE)
dt <- simdt$dt

# drop singleton cross-cohorts: their variance is not estimable, so every call would warn
inv <- unique(dt[, .(unit, G, G2)])
keep_units <- inv[, if(.N >= 2) .SD, by = .(G, G2)][, unit]
dt <- dt[unit %in% keep_units]

# basic tests ------------------------------------------------------------------

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2"),
              info = "double call")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", event_specific = FALSE),
              info = "double, combined effect")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", allow_unbalance_panel = TRUE),
              info = "double, unbalanced panel")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", control_option = "never"),
              info = "double, only never")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", double_control_option = "never"),
              info = "double call, double never")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", double_control_option = "notyet"),
              info = "double cal, double notyet")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      covariatesvar = "x",
                      cohortvar2 = "G2"),
              info = "double, covariates")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", anticipation = 0, anticipation2 = 1),
              info = "double, one horizon for each event")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_group_time",
                      cohortvar2 = "G2"),
              info = "double, ggt")

expect_silent(fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "dynamic_stagger",
                      cohortvar2 = "G2"),
              info = "double, dynamic_sq")

# non-standard data ---------------------------------------------------------

dt2 <- data.table::copy(dt)
dt2 <- dt2[G < 4]
expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", allow_unbalance_panel = TRUE, control_option = "notyet"),
              info = "double, limited treated group")


#all G2 > G
dt2 <- data.table::copy(dt)
dt2 <- dt2[G < G2]
expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", allow_unbalance_panel = TRUE, control_option = "notyet"),
              info = "double, G2 > G")

#all G > G2
dt2 <- data.table::copy(dt)
dt2 <- dt2[G2 < G]
expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", allow_unbalance_panel = TRUE, control_option = "notyet"),
              info = "double, G > G2")


dt2 <- data.table::copy(dt)
keep <- sample(c(rep(TRUE, 15),FALSE), dt2[,.N], TRUE)
dt2 <- dt2[keep]
# dropping rows leaves some group-period with one effective unit, which warns
expect_warning(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                       cohortvar2 = "G2", allow_unbalance_panel = TRUE),
               pattern = "fewer than 2 effective units",
               info = "double, unbalanced panel, unbalanced data")

dt2 <- data.table::copy(dt)
dt2[, time := time*2 + 3]
dt2[, G := G*2 + 3]
dt2[, G2 := G2*2 + 3]

expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2"),
             info = "time offset")

if(.Platform$OS.type == "unix" & at_home() & requireNamespace("parallel")){
  expect_silent(fastdid(dt, timevar = "time",cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                        cohortvar2 = "G2", parallel = TRUE),
                info = "parallel double")
}

# M=3 tests ----------------------------------------------------------------

set.seed(42)
N3 <- 500; TT3 <- 6
g3_main  <- sample(c(2L, 3L, 4L, Inf), N3, replace = TRUE, prob = c(0.2, 0.2, 0.2, 0.4))
g3_conf1 <- sample(c(2L, 3L, 5L, Inf), N3, replace = TRUE, prob = c(0.2, 0.2, 0.2, 0.4))
g3_conf2 <- sample(c(3L, 4L, 6L, Inf), N3, replace = TRUE, prob = c(0.2, 0.2, 0.2, 0.4))
dt_m3 <- data.table::data.table(
  unit = rep(1:N3, each = TT3),
  time = rep(1:TT3, N3),
  G    = rep(g3_main,  each = TT3),
  G2   = rep(g3_conf1, each = TT3),
  G3   = rep(g3_conf2, each = TT3),
  y    = rnorm(N3 * TT3)
)

# drop singleton triple-cohorts: their variance is not estimable, so every call would warn
inv3 <- unique(dt_m3[, .(unit, G, G2, G3)])
keep3 <- inv3[, if(.N >= 2) .SD, by = .(G, G2, G3)][, unit]
dt_m3 <- dt_m3[unit %in% keep3]

expect_silent(
  fastdid(dt_m3, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
          result_type = "group_time", cohortvar2 = c("G2", "G3")),
  info = "M=3 basic run"
)

expect_silent(
  fastdid(dt_m3, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
          result_type = "group_time", cohortvar2 = c("G2", "G3"), event_specific = FALSE),
  info = "M=3 combined effect"
)

expect_silent(
  fastdid(dt_m3, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
          result_type = "dynamic", cohortvar2 = c("G2", "G3")),
  info = "M=3 dynamic"
)

expect_silent(
  fastdid(dt_m3, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
          result_type = "group_time", cohortvar2 = c("G2", "G3"), double_control_option = "never"),
  info = "M=3 double never"
)

expect_silent(
  fastdid(dt_m3, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
          result_type = "group_group_time", cohortvar2 = c("G2", "G3")),
  info = "M=3 group_group_time"
)
res_m3_ggt <- fastdid(dt_m3, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                      result_type = "group_group_time", cohortvar2 = c("G2", "G3"))
expect_true(all(c("cohort1", "cohort2", "cohort3", "time") %in% names(res_m3_ggt)),
            info = "M=3 group_group_time has cohort1/2/3 columns")

# correctness --------------------------------------------------------------

#' A panel with unit and time fixed effects, known event effects, and small noise.
#' Event 1 is flat at 1. Event 2 grows by 0.5 a period. Event 3 grows by 3 a period.
#' Each element of `cohorts` is the triple (g1, g2, g3) of one cohort.
make_double_dt <- function(cohorts, TT = 6, n = 500, seed = 7) {
  set.seed(seed)
  dt <- data.table::rbindlist(lapply(seq_along(cohorts), function(i) {
    g <- cohorts[[i]]
    data.table::CJ(unit = (i - 1) * n + seq_len(n), time = seq_len(TT))[
      , `:=`(G = g[1], G2 = g[2], G3 = g[3])]
  }))
  dt[, y := unit / 100 + time * 0.3 + stats::rnorm(.N, 0, 0.01) +
       data.table::fifelse(time >= G, 1, 0) +
       data.table::fifelse(time >= G2, 0.5 * (time - G2 + 1), 0) +
       data.table::fifelse(time >= G3, 3 * (time - G3 + 1), 0)]
  dt[]
}

# imputation case: the first event comes before the confounding one (g1 < g')
dt_imp <- make_double_dt(list(c(2, 4, Inf), # target
                              c(2, Inf, Inf), # C^imp control, same g1, never confounded
                              c(Inf, Inf, Inf))) # never treated, for the first stage
res_imp <- fastdid(dt_imp, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                   result_type = "group_group_time", cohortvar2 = "G2")
att_imp <- res_imp[cohort1 == 2 & cohort2 == 4 & time >= 4, att]
expect_equal(length(att_imp), 3L, info = "double, imputation case has all post periods")
expect_equal(att_imp, rep(1, length(att_imp)), tolerance = tol,
             info = "double, imputation case recovers ATT^1")

# did case: the confounding event comes before the first one (g1 > g')
dt_did <- make_double_dt(list(c(4, 2, Inf), # target
                              c(Inf, 2, Inf), # C^did control, same g2, not yet treated
                              c(Inf, Inf, Inf)))
res_did <- fastdid(dt_did, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                   result_type = "group_group_time", cohortvar2 = "G2")
att_did <- res_did[cohort1 == 4 & cohort2 == 2 & time >= 4, att]
expect_equal(length(att_did), 3L, info = "double, did case has all post periods")
expect_equal(att_did, rep(1, length(att_did)), tolerance = tol,
             info = "double, did case recovers ATT^1")

# M=3 did case: a control with the same g2 but an extra third event must not be used
dt_con <- make_double_dt(list(c(4, 2, Inf), # target
                              c(Inf, 2, 3), # same g2, but event 3 is active: not a valid control
                              c(Inf, 2, Inf), # the valid control
                              c(Inf, Inf, Inf)))
res_con <- fastdid(dt_con, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                   result_type = "group_group_time", cohortvar2 = c("G2", "G3"))
att_con <- res_con[cohort1 == 4 & cohort2 == 2 & is.infinite(cohort3) & time >= 4, att]
expect_equal(length(att_con), 3L, info = "M=3 did case has all post periods")
expect_equal(att_con, rep(1, length(att_con)), tolerance = tol,
             info = "M=3 did case excludes the contaminated control")

