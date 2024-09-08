# setup ----------------------

library(did)
#set.seed(1)
tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+03, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE,
                 seed = 1, stratify = FALSE,
                 second_cov = TRUE)
dt <- simdt$dt

est_diff_ratio <- function(result, did_result){
  did_result_dt <- data.table::data.table(cohort = did_result$group, time = did_result$t, did_att = did_result$att, did_se = did_result$se)
  did_result_dt <- did_result_dt[did_att != 0]
  compare <- did_result_dt |> merge(result, by = c("cohort", "time"), all = TRUE) 
  att_diff_per <- compare[, sum(abs(did_att-att), na.rm = TRUE)/sum(did_att, na.rm = TRUE)]
  se_diff_per <- compare[, sum(abs(did_se-se), na.rm = TRUE)/sum(did_se, na.rm = TRUE)]
  return(c(att_diff_per , se_diff_per))
}

est_diff_ratio_agg <- function(result, did_result){
  names(result) <- c("target", "att", "se", "outcome", "atub", "atlb")
  did_result_dt <- data.table::data.table(target = did_result$egt, did_att = did_result$att.egt, did_se = did_result$se.egt)
  compare <- did_result_dt |> merge(result, by = c("target"), all = TRUE) 
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

# simple notyet only ----------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time", control_option = "notyet")
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt[!is.infinite(G)],base_period = "universal",est_method = "ipw",cband = FALSE,
                          #xformla = ~x,
                          control_group = "notyettreated",
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "notyet")
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



# dynamic -------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "dynamic")
did_result_gt <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                             #xformla = ~x,
                             control_group = "notyettreated",
                             clustervars = NULL,
                             bstrap = FALSE)
did_result <- did::aggte(did_result_gt, type = "dynamic")
expect_equal(est_diff_ratio_agg(result, did_result), c(0,0), tolerance = tol,
             info = "dynamic")
rm(result, did_result)

# group -------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group")
result <- result[type == "post",]
result[, type := NULL]
did_result_gt <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                             #xformla = ~x,
                             control_group = "notyettreated",
                             clustervars = NULL,
                             bstrap = FALSE)
did_result <- did::aggte(did_result_gt, type = "group")
expect_equal(est_diff_ratio_agg(result, did_result), c(0,0), tolerance = tol,
             info = "group")
rm(result, did_result)

# time --------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "time")
result <- result[type == "post",]
result[, type := NULL]
did_result_gt <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                             #xformla = ~x,
                             control_group = "notyettreated",
                             clustervars = NULL,
                             bstrap = FALSE)
did_result <- did::aggte(did_result_gt, type = "calendar")
expect_equal(est_diff_ratio_agg(result, did_result), c(0,0), tolerance = tol,
             info = "time")
rm(result, did_result)

# bootstrap clustered --------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "simple")
did_result_gt <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                             #xformla = ~x,
                             control_group = "notyettreated",
                             clustervars = NULL,
                             bstrap = FALSE)
did_result <- did::aggte(did_result_gt, type = "simple")
diffs <- c((did_result$overall.att-result[type == "post", att])/did_result$overall.att,
           (did_result$overall.se-result[type == "post", se])/did_result$overall.se)
expect_equal(diffs, c(0,0), tolerance = tol,
             info = "simple group")
rm(result, did_result)

# dynamic with balanced cohort -------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "dynamic",
                  balanced_event_time = 4)
did_result_gt <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                             #xformla = ~x,
                             control_group = "notyettreated",
                             clustervars = NULL,
                             bstrap = FALSE)
did_result <- did::aggte(did_result_gt, type = "dynamic",
                         balance_e = 4)
expect_equal(est_diff_ratio_agg(result, did_result), c(0,0), tolerance = tol,
             info = "dynamic with balanced cohort")
rm(result, did_result)

# unbalanced ------------------------------

dt2 <- data.table::copy(dt)
keep <- sample(c(rep(TRUE, 15),FALSE), dt2[,.N], TRUE)
dt2 <- dt2[keep]

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  allow_unbalance_panel = TRUE)
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                          #xformla = ~x,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = TRUE,
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "unbalanced method, balance panel, simple")
rm(result, did_result)

result <- fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  allow_unbalance_panel = TRUE)
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt2,base_period = "universal",est_method = "ipw",cband = FALSE,
                          #xformla = ~x,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = TRUE,
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "unbalanced method, unbalance panel, simple")
rm(result, did_result)

# unbalanced, ipw  ------------------------------------------

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  control_type = "ipw",
                  covariatesvar = c("x", "x2"),
                  allow_unbalance_panel = TRUE)
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",cband = FALSE,
                          est_method = "ipw",xformla = ~x+x2,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = TRUE,
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "unbalanced method, balance panel, ipw")
rm(result, did_result)

result <- fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  
                  result_type = "group_time",
                  control_type = "ipw",
                  covariatesvar = c("x", "x2"),
                  allow_unbalance_panel = TRUE)
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt2,
                          base_period = "universal",cband = FALSE,
                          est_method = "ipw",xformla = ~x+x2,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = TRUE,
                          clustervars = NULL,
                          bstrap = FALSE)

expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "unbalanced method,unbalanced panel, data unbalanced, ipw")
rm(result, did_result)

# unbalanced, reg --------------------------------------------------------------------------

# reg and dr unbalanced will be used after the OR influence issues is resolved

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                  control_type = "reg",
                  covariatesvar = c("x", "x2"),
                  allow_unbalance_panel = TRUE)
did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",cband = FALSE,
                          est_method = "reg",xformla = ~x+x2,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = TRUE,
                          clustervars = NULL,
                          bstrap = FALSE)
expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
             info = "unbalanced panel reg, but data balanced")
rm(result, did_result)

# estimate not the same, see coverage test
# result <- fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
#                   control_type = "reg",
#                   covariatesvar = c("x", "x2"),
#                   allow_unbalance_panel = TRUE)
# did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt2 ,base_period = "universal",cband = FALSE,
#                           est_method = "reg",xformla = ~x+x2,
#                           control_group = "notyettreated",
#                           allow_unbalanced_panel = TRUE,
#                           clustervars = NULL,
#                           bstrap = FALSE)
# 
# expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
#              info = "unbalanced panel reg, but data unbalanced")
# rm(result, did_result)

#clean up environment -------------------------
rm(dt, simdt, tol, est_diff_ratio_agg)
detach("package:did", unload=TRUE)