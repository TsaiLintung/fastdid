boot <- TRUE
clust <- "G"

did_result <- att_gt(yname = "y",
                     gname = "G",
                     idname = "unit",
                     tname = "time",
                     data = dt,
                     #xformla = ~x,
                     base_period = "universal",
                     control_group = "notyettreated",
                     est_method = "ipw",
                     clustervars = clust,
                     bstrap = boot,
                     biters = 10000)

did_result2 <- att_gt(yname = "y",
                      gname = "G",
                      idname = "unit",
                      tname = "time",
                      data = dt,
                      #xformla = ~x,
                      base_period = "universal",
                      control_group = "notyettreated",
                      est_method = "ipw",
                      clustervars = clust,
                      bstrap = boot,
                      biters = 10000)


did_result_dt <- data.table(cohort = did_result$group, time = did_result$t, did_att = did_result$att, did_se = did_result$se, did_se2 = did_result2$se)

setorder(did_result_dt, time, cohort)


compare <- did_result_dt |> merge(result, by = c("cohort", "time"), all = TRUE) 
compare[, ratio := did_se/se]
compare[, ratio2 := did_se/did_se2]

compare[, mean(ratio, na.rm = TRUE)]
compare[, mean(ratio2, na.rm = TRUE)]
# 
# compare <- compare[!is.na(did_se)]
# compare[, se_diff := se-did_se]
# compare[, se_rand := did_se2-did_se]
# compare[, att_diff := att-did_att]
# compare[, se_diff_ratio := abs(se_diff/did_se)]
# compare[, att_diff_ratio := abs(att_diff/did_att)]
# compare[, se_rand_ratio := abs(se_rand/did_se)]
# View(compare)
# 
# 
# timetaken(started.at)
# result