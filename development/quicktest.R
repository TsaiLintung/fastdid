
library(did)
library(DRDID)

tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(5e02, 3, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE, seed = 1, stratify = FALSE,
                 second_cov = TRUE)
dt <- simdt$dt
dt2 <- copy(dt)
keep <- sample(c(rep(TRUE, 15),FALSE), dt2[,.N], TRUE)
dt2 <- dt2[keep]


est_diff_ratio <- function(result, did_result){
  did_result_dt <- data.table(cohort = did_result$group, time = did_result$t, did_att = did_result$att, did_se = did_result$se)
  compare <- did_result_dt |> merge(result, by = c("cohort", "time"), all = TRUE) 
  att_diff_per <- compare[, sum(abs(did_att-att), na.rm = TRUE)/sum(did_att, na.rm = TRUE)]
  se_diff_per <- compare[, sum(abs(did_se-se), na.rm = TRUE)/sum(did_se, na.rm = TRUE)]
  return(c(att_diff_per , se_diff_per))
}
# unbalanced ------------------------------



balanced_panel_rc <- function(method){
  
  result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                    control_type = method,
                    allow_unbalance_panel = TRUE,
                    covariatesvar = c("x", "x2"))
  did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",
                            est_method = method,#method,
                            cband = FALSE,
                            xformla = ~x+x2,
                            allow_unbalanced_panel = TRUE,
                            control_group = "notyettreated",
                            clustervars = NULL,
                            bstrap = FALSE)
  print(est_diff_ratio(result, did_result))
}
  


unbalanced_panel_rc <- function(method){
  
  result <- fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                    control_type = method,
                    allow_unbalance_panel = TRUE,
                    covariatesvar = c("x", "x2"))
  did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt2,base_period = "universal",
                            est_method = method,#my_reg_did_rc
                            cband = FALSE,
                            xformla = ~x+x2,
                            allow_unbalanced_panel = TRUE,
                            control_group = "notyettreated",
                            clustervars = NULL,
                            bstrap = FALSE)
  print(est_diff_ratio(result, did_result))
}

simple_rc <- function(){
  
  result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                    control_type = "ipw",
                    allow_unbalance_panel = TRUE)
  did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                            allow_unbalanced_panel = TRUE,
                            control_group = "notyettreated",
                            clustervars = NULL,
                            bstrap = FALSE)
  print(est_diff_ratio(result, did_result))
}

test_again <- function(method){
  source("interactive/custom_func.R")
  load_all()
  simple_rc()
  balanced_panel_rc(method)
  unbalanced_panel_rc(method)
}