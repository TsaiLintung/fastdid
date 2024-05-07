# setup ----------------------

library(did)
#set.seed(1)
tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+03, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE,
                 seed = 1, stratify = FALSE,
                 second_cov = TRUE)
dt <- simdt$dt

est_diff_ratio <- function(result, did_result){
  did_result_dt <- data.table(cohort = did_result$group, time = did_result$t, did_att = did_result$att, did_se = did_result$se)
  did_result_dt <- did_result_dt[did_att != 0]
  compare <- did_result_dt |> merge(result, by = c("cohort", "time"), all = TRUE) 
  att_diff_per <- compare[, sum(abs(did_att-att), na.rm = TRUE)/sum(did_att, na.rm = TRUE)]
  se_diff_per <- compare[, sum(abs(did_se-se), na.rm = TRUE)/sum(did_se, na.rm = TRUE)]
  return(c(att_diff_per , se_diff_per))
}

# simple -------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time")
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                          #xformla = ~x,
                          control_group = "notyettreated",
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "simple")
rm(result, did_result)

# never -------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  control_option = "never")
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                          #xformla = ~x,
                          control_group = "nevertreated",
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "never treated")
rm(result, did_result)

# covariates -------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  control_type = "ipw",
                  covariatesvar = c("x", "x2"))
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                          xformla = ~x+x2,
                          control_group = "notyettreated",
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "covariates ipw")
rm(result, did_result)

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  control_type = "reg",
                  covariatesvar = c("x", "x2"))
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "reg",cband = FALSE,
                          xformla = ~x+x2,
                          control_group = "notyettreated",
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "covariates reg")
rm(result, did_result)


result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  control_type = "dr",
                  covariatesvar = c("x"))
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "dr",cband = FALSE,
                          xformla = ~x,
                          control_group = "notyettreated",
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "covariates dr")
rm(result, did_result)

# anticipation -----------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  control_type = "dr",
                  covariatesvar = c("x"), anticipation = 2)
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "dr",cband = FALSE,
                          xformla = ~x,
                          control_group = "notyettreated",
                          clustervars = NULL,
                          bstrap = FALSE, anticipation = 2)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "anticipation")
rm(result, did_result)

# alternative base period -----------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  control_type = "dr",
                  covariatesvar = c("x"), base_period = "varying")
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "varying",est_method = "dr",cband = FALSE,
                          xformla = ~x,
                          control_group = "notyettreated",
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "varying base period")
rm(result, did_result)

# bootstrap --------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  boot = TRUE, biters = 10000)
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                          control_group = "notyettreated",
                          clustervars = NULL,
                          bstrap = TRUE, biters = 10000)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "bootstrap")
rm(result, did_result)

# bootstrap clustered --------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  clustervar = "x",
                  boot = TRUE, biters = 10000)
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                          control_group = "notyettreated",
                          clustervars = "x",
                          bstrap = TRUE, biters = 10000)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "bootstrap clustered")
rm(result, did_result)

#multiple outcome ------------------------------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = c("y", "y2"),  result_type = "group_time")
did_result <- did::att_gt(yname = "y2",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                          #xformla = ~x,
                          control_group = "notyettreated",
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result[outcome == "y2"], did_result), c(0,0), tolerance = tol,
             info = "simple multiple outcome")
rm(result, did_result)


result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = c("y", "y2"),  result_type = "group_time",
                  covariatesvar = "x")
did_result <- did::att_gt(yname = "y2",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                          xformla = ~x,
                          control_group = "notyettreated",
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result[outcome == "y2"], did_result), c(0,0), tolerance = tol,
             info = "multiple outcome with covariates")
rm(result, did_result)
#clean up environment
detach("package:did", unload=TRUE)