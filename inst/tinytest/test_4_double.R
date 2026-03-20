tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+03, 5, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE, second_cohort = TRUE)
dt <- simdt$dt

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
expect_silent(fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                      cohortvar2 = "G2", allow_unbalance_panel = TRUE),
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
dt_m3 <- data.table(
  unit = rep(1:N3, each = TT3),
  time = rep(1:TT3, N3),
  G    = rep(g3_main,  each = TT3),
  G2   = rep(g3_conf1, each = TT3),
  G3   = rep(g3_conf2, each = TT3),
  y    = rnorm(N3 * TT3)
)

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

