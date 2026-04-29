# setup ----------------------

library(devtools)
library(fastdid)
library(did)

set.seed(1)
tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+03, 3, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE,
                 seed = 1, stratify = FALSE,
                 second_cov = TRUE)
dt <- simdt$dt

dt[, weight := pmax(rnorm(.N, 0.8, 0.2), 0.3)] # weight is time-varying

did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "reg",cband = FALSE,
                          #xformla = ~x,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = TRUE,
                          clustervars = NULL,
                          faster_mode = TRUE,
                          weightsname = "weight",
                          bstrap = FALSE)

did_result2 <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "reg",cband = FALSE,
                          #xformla = ~x,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = TRUE,
                          clustervars = NULL,
                          faster_mode = FALSE,
                          weightsname = "weight",
                          bstrap = FALSE)

did_result3 <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "reg",cband = FALSE,
                          #xformla = ~x,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = TRUE,
                          clustervars = NULL,
                          faster_mode = TRUE,
                          bstrap = FALSE)

did_result4 <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "reg",cband = FALSE,
                          #xformla = ~x,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = TRUE,
                          clustervars = NULL,
                          faster_mode = FALSE,
                          bstrap = FALSE)

did_result5 <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "reg",cband = FALSE,
                          #xformla = ~x,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = TRUE,
                          clustervars = "x",
                          faster_mode = TRUE,
                          biters = 20000,
                          bstrap = TRUE)

did_result6<- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "reg",cband = FALSE,
                          #xformla = ~x,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = TRUE,
                          clustervars = "x",
                          faster_mode = FALSE,
                          biters = 20000,
                          bstrap = TRUE)