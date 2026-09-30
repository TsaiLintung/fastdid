# fastdid_cells() and combine_cells(): the no-covariate cells and their linear combinations

simdt <- sim_did(2000, 6, seed = 3, second_cohort = TRUE)
dt <- simdt$dt

# the cells match fastdid() on one event, with the universal base and no covariates
fc <- fastdid_cells(dt, "time", "unit", "y", "G")
ref <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
               result_type = "group_time", base_period = "universal", control_option = "both")
ref <- ref[time >= cohort]
own <- merge(fc$cells, fc$cohorts[, .(row, G)], by = "row")[s >= G]
cmp <- merge(own, ref[, .(G = cohort, s = time, att_ref = att, se_ref = se)], by = c("G", "s"))
expect_equal(nrow(cmp), nrow(ref), info = "every post cell of fastdid() is a cell")
expect_equal(cmp$att, cmp$att_ref, tolerance = 1e-8, info = "the cells equal fastdid()")

# one cell as a target: the same estimate, and the same standard error up to fastdid's small-group inflation
w <- own[, .(target = paste(G, s), row, s, w = 1)]
est <- combine_cells(fc, w)
cmp2 <- merge(cmp[, .(target = paste(G, s), att, se_ref)], est, by = "target")
expect_equal(cmp2$estimate, cmp2$att, tolerance = 1e-8, info = "a one-cell target returns the cell")
expect_true(all(abs(cmp2$se / cmp2$se_ref - 1) < 0.02), info = "the analytic se matches fastdid()")

# a linear combination is the combination of the cells, and the influence function sums to zero
w2 <- rbind(w[1:3][, target := "sum"], w[4][, `:=`(target = "sum", w = -2)])
r <- combine_cells(fc, w2, return_if = TRUE)
expect_equal(r$estimates$estimate, sum(own$att[1:3]) - 2 * own$att[4], tolerance = 1e-8)
expect_equal(mean(r$inf_func[, 1]), 0, tolerance = 1e-8)

# a weight on the means: the mean change of cohort 1 from its first to its last period
wm <- data.table(target = "mean", row = 1L, s = c(1L, 6L), w = c(-1, 1))
r3 <- combine_cells(fc, w2[0], w_mean = wm)
expect_equal(r3$estimate, fc$M[1, 6] - fc$M[1, 1], tolerance = 1e-10)

# several event columns: the cohort is the vector of dates
fc2 <- fastdid_cells(dt, "time", "unit", "y", c("G", "G2"))
expect_equal(nrow(fc2$cohorts), uniqueN(dt[, .(G, G2)][, first := pmin(G, G2)][first > 1]),
             info = "one cohort per observed vector of dates")

# defenses
expect_error(fastdid_cells(dt[!(unit == 1 & time == 2)], "time", "unit", "y", "G"), "balanced")
