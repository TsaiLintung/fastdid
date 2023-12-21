# setup ----------------------

library(did)

tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e03, 10, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, stratify = FALSE,
                 second_cov = TRUE)
dt <- simdt$dt
est_diff_ratio <- function(result, did_result){
  did_result_dt <- data.table(cohort = did_result$group, time = did_result$t, did_att = did_result$att, did_se = did_result$se)
  compare <- did_result_dt |> merge(result, by = c("cohort", "time"), all = TRUE) 
  att_diff_per <- compare[, sum(abs(did_att-att), na.rm = TRUE)/sum(did_att, na.rm = TRUE)]
  se_diff_per <- compare[, sum(abs(did_se-se), na.rm = TRUE)/sum(did_se, na.rm = TRUE)]
  return(c(att_diff_per , se_diff_per))
}
# unbalanced ------------------------------

dt2 <- copy(dt)
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

# unbalanced method with balanced panel ------------------------------------------

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

# reg and dr unbalanced will be used after the OR influence issues is resolved

# result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
#                   control_type = "reg",
#                   covariatesvar = c("x", "x2"),
#                   allow_unbalance_panel = TRUE)
# did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",cband = FALSE,
#                           est_method = "reg",xformla = ~x+x2,
#                           control_group = "notyettreated",
#                           allow_unbalanced_panel = TRUE,
#                           clustervars = NULL,
#                           bstrap = FALSE)
# 
# expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
#              info = "unbalanced panel reg, but data balanced")
# rm(result, did_result)
# 
# result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
#                   control_type = "dr",
#                   covariatesvar = c("x", "x2"),
#                   allow_unbalance_panel = TRUE)
# did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",cband = FALSE,
#                           est_method = "dr",xformla = ~x+x2,
#                           control_group = "notyettreated",
#                           allow_unbalanced_panel = TRUE,
#                           clustervars = NULL,
#                           bstrap = FALSE)
# 
# expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
#              info = "unbalanced panel dr, but data balanced")
# rm(result, did_result)
# 
# 
# # unbalanced panel unbalanced method ------------------------------------------------------------------
# 
# 
# 
# 
# result <- fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
#                   control_type = "reg",
#                   covariatesvar = c("x", "x2"),
#                   allow_unbalance_panel = TRUE)
# did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt2,base_period = "universal",cband = FALSE,
#                           est_method = "reg",xformla = ~x+x2,
#                           control_group = "notyettreated",
#                           allow_unbalanced_panel = TRUE,
#                           clustervars = NULL,
#                           bstrap = FALSE)
# 
# expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
#              info = "unbalanced panel reg, data unbalanced")
# rm(result, did_result)
# 
# result <- fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
#                   control_type = "dr",
#                   covariatesvar = c("x", "x2"),
#                   allow_unbalance_panel = TRUE)
# did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt2,base_period = "universal",cband = FALSE,
#                           est_method = "dr",xformla = ~x+x2,
#                           control_group = "notyettreated",
#                           allow_unbalanced_panel = TRUE,
#                           clustervars = NULL,
#                           bstrap = FALSE)
# 
# expect_equal(est_diff_ratio(result, did_result), c(0,0), tolerance = tol,
#              info = "unbalanced panel dr, data unbalanced")
# rm(result, did_result)

#clean up environment
detach("package:did", unload=TRUE)