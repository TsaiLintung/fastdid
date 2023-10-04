
rm(list = ls())
gc()

library(devtools)
library(profvis)

setwd("~/GitHub/EventStudyCode")


# setup --------------------------------

load_all()
#load_all("~/Github/did")
library(did)

tol <- 0.01 #allow 1% different between estimates
simdt <- sim_did(1000, 10, cov = "no", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, stratify = FALSE, 
                 epsilon_size = 1)
dt <- simdt$dt

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", result_type = "dynamic")
did_result_gt <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "ipw",cband = FALSE,
                             #xformla = ~x,
                             control_group = "notyettreated",
                             clustervars = NULL,
                             bstrap = FALSE)
did_result <- did::aggte(did_result_gt, type = "dynamic")

names(result) <- c("target", "att", "se")
did_result_dt <- data.table(target = did_result$egt, did_att = did_result$att.egt, did_se = did_result$se.egt)
compare <- did_result_dt |> merge(result, by = c("target"), all = TRUE) 