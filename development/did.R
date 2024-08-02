library(did)

rm(list = ls());gc()

data(mpdta)
mpdta$x <- rnorm(nrow(mpdta))

#the unbalanced panel OR influence function calculation (DRDID::drdid_rc.R, 230, should be post - pre - or), since control function is subtracted
#simultaneous valid band, should subset t>g  before getting the new critical value? (compute.aggte 376-383)


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
