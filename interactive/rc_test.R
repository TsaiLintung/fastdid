
ipw_unbalanced_data <- function(method){
  result <- fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                    control_type = method,
                    covariatesvar = c("x", "x2"),
                    allow_unbalance_panel = TRUE)
  did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt2,
                            base_period = "universal",cband = FALSE,
                            est_method = drdid_rc1,xformla = ~x+x2,
                            control_group = "notyettreated",
                            allow_unbalanced_panel = TRUE,
                            clustervars = NULL,
                            bstrap = FALSE)
  est_diff_ratio(result, did_result)
}

ipw_balanced_data <- function(method){
  result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                    control_type = method,
                    covariatesvar = c("x", "x2"),
                    allow_unbalance_panel = TRUE)
  did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,
                            base_period = "universal",cband = FALSE,
                            est_method = drdid_rc,xformla = ~x+x2,
                            control_group = "notyettreated",
                            allow_unbalanced_panel = TRUE,
                            clustervars = NULL,
                            bstrap = FALSE)
  est_diff_ratio(result, did_result)
}

simple_unbalanced <- function(){
  result <- fastdid(dt2, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
                    allow_unbalance_panel = TRUE)
  did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt2,base_period = "universal",est_method = "ipw",cband = FALSE,
                            #xformla = ~x,
                            control_group = "notyettreated",
                            allow_unbalanced_panel = TRUE,
                            clustervars = NULL,
                            bstrap = FALSE)
  
  est_diff_ratio(result, did_result)
}

run_new <- function(method){
  load_all()
  print(simple_unbalanced())
  print(ipw_balanced_data(method))
  print(ipw_unbalanced_data(method))
}
