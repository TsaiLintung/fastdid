did_result <- att_gt(yname = "y",
                     gname = "G",
                     idname = "unit",
                     tname = "time",
                     data = dt,
                     #xformla = ~x,
                     base_period = "universal",
                     control_group = "notyettreated",
                     est_method = "ipw",
                     clustervars = NULL,
                     bstrap = boot)

did_result2 <- att_gt(yname = "y",
                      gname = "G",
                      idname = "unit",
                      tname = "time",
                      data = dt,
                      #xformla = ~x,
                      base_period = "universal",
                      control_group = "notyettreated",
                      est_method = "ipw",
                      clustervars = NULL,
                      bstrap = boot)


did_result_dt <- data.table(G = did_result$group, time = did_result$t, did_att = did_result$att, did_se = did_result$se, did_se2 = did_result2$se)
did_result_dt[,target := G*max(time)+time]
compare <- did_result_dt |> merge(result, by = c("target"), all = TRUE) 
compare <- compare[!is.na(did_se)]
compare[, se_diff := se-did_se]
compare[, se_rand := did_se2-did_se]
compare[, att_diff := att-did_att]
compare[, se_diff_ratio := abs(se_diff/did_se)]
compare[, att_diff_ratio := abs(att_diff/did_att)]
compare[, se_rand_ratio := abs(se_rand/did_se)]
View(compare)


timetaken(started.at)
result