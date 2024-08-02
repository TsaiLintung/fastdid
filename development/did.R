library(did)
library(fastdid)

rm(list = ls());gc()

data(mpdta)
mpdta$x <- rnorm(nrow(mpdta))

#the unbalanced panel OR influence function calculation (DRDID::drdid_rc.R, 230, should be post - pre - or), since control function is subtracted

# did ---------

# repeated cross section with 
out <- att_gt(yname="lemp",
              tname="year",
              idname="countyreal",
              gname="first.treat",
              data=mpdta,
              bstrap = FALSE,
              allow_unbalanced_panel = TRUE)
agg <- aggte(out, type = "simple")
agg$overall.se

out <- att_gt(yname="lemp",
              tname="year",
              idname="countyreal",
              gname="first.treat",
              data=mpdta,
              bstrap = FALSE,
              allow_unbalanced_panel = FALSE)
agg <- aggte(out, type = "simple")
agg$overall.se

for(est in c("dr", "ipw", "reg")){
  for(p in c(FALSE, TRUE)){
    out <- att_gt(yname="lemp",
                     tname="year",
                     idname="countyreal",
                     gname="first.treat",
                     xformla=~x,
                     data=mpdta,
                     allow_unbalanced_panel = p,
                     bstrap = FALSE,
                     est_method = est)
    agg <- aggte(out, type = "simple")
    message("se for allow_unbalance_panel = ", p, " est = ", est, ": ", agg$overall.se)
    rm(out, agg)
  }
}

# fastdid ------------
mpdta <- as.data.table(mpdta)
setnames(mpdta, "first.treat", "firsttreat")
result <- fastdid(mpdta, timevar = "year", cohortvar = "firsttreat", unitvar = "countyreal",outcomevar = "lemp",  result_type = "simple",
                  allow_unbalance_panel = FALSE, covariatesvar = "x")
