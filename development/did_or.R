library(did)
# install.packages("devtools")
#devtools::install_github("TsaiLintung/fastdid")
library(fastdid)

simdt <- sim_did(1e+03, 5, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = TRUE,
                  stratify = FALSE,
                 second_cov = TRUE)
dt <- simdt$dt

did_result <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "dr",cband = FALSE,
                          xformla = ~x+x2,
                          control_group = "notyettreated",
                          allow_unbalanced_panel = FALSE, #this is the only differece
                          clustervars = NULL,
                          bstrap = FALSE)
did_result2 <- did::att_gt(yname = "y",gname = "G",idname = "unit",tname = "time",data = dt,base_period = "universal",est_method = "dr",cband = FALSE,
                           xformla = ~x+x2,
                           control_group = "notyettreated",
                           allow_unbalanced_panel = TRUE, #this is the only differece
                           clustervars = NULL,
                           bstrap = FALSE)

sum(did_result$att, na.rm = TRUE)/sum(did_result2$att, na.rm = TRUE)
sum(did_result$se, na.rm = TRUE)/sum(did_result2$se, na.rm = TRUE) #should be 1 if the SE is the same, but not??????

#very weird
#DRDID:drdid_rc1, line 203: inf.cont <- inf.cont.post - inf.cont.pre + inf.cont.ps + inf.cont.or
#DRDID:drdid_panel, line 169: inf.control <- (inf.cont.1 + inf.cont.2 - inf.cont.3) / mean(w.cont)
#influence from outcome have different sign??