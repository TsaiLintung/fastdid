# setup ----------------------

library(did)

tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e03, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, stratify = FALSE)
dt <- simdt$dt

est_diff_ratio_agg <- function(result, did_result){
  names(result) <- c("target", "att", "se", "outcome", "atub", "atlb")
  did_result_dt <- data.table::data.table(target = did_result$egt, did_att = did_result$att.egt, did_se = did_result$se.egt)
  compare <- did_result_dt |> merge(result, by = c("target"), all = TRUE) 
  att_diff_per <- compare[, sum(abs(did_att-att), na.rm = TRUE)/sum(did_att, na.rm = TRUE)]
  se_diff_per <- compare[, sum(abs(did_se-se), na.rm = TRUE)/sum(did_se, na.rm = TRUE)]
  return(c(att_diff_per , se_diff_per))
}

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

#clean up environment -------------------------
rm(dt, simdt, tol, est_diff_ratio_agg)
detach("package:did", unload=TRUE)

