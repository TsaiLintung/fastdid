#setwd("~/Documents/GitHub/fastdid")

library(here)
library(devtools)
library(tinytest)
library(roxygen2)
library(profvis)

setwd(here())

load_all()



tol <- 1e-2 #allow 1% different between estimates
simdt <- sim_did(1e+03, 4, cov = "cont", hetero = "all", balanced = TRUE, second_outcome = FALSE, seed = 1, 
                 stratify = FALSE, second_cov = TRUE, vary_cov = TRUE, second_cohort = TRUE)
dt <- simdt$dt
res <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit",outcomevar = "y",  result_type = "group_time",
        cohortvar2 = "G2", control_option = "notyet")

# did ------

library(did)
library(data.table)

rm(list = ls());gc()

#suspected problem: the unbalanced panel OR influence function calculation (DRDID::drdid_rc.R, 230, should be post - pre - or), since control function is subtracted

#data
dt <- data.table(
  time = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3),
  G = c(2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, Inf, Inf, Inf, Inf, Inf, Inf, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, Inf, Inf, Inf, Inf, Inf, Inf, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, Inf, Inf, Inf, Inf, Inf, Inf),
  unit = c(6, 9, 11, 13, 14, 16, 17, 1, 3, 4, 8, 10, 15, 19, 2, 5, 7, 12, 18, 20, 6, 9, 11, 13, 14, 16, 17, 1, 3, 4, 8, 10, 15, 19, 2, 5, 7, 12, 18, 20, 6, 9, 11, 13, 14, 16, 17, 1, 3, 4, 8, 10, 15, 19, 2, 5, 7, 12, 18, 20),
  x = c(4, 5, 4, 4, 4, 4, 1, 1, 4, 2, 1, 2, 4, 5, 1, 4, 3, 1, 2, 4, 4, 5, 4, 4, 4, 4, 1, 1, 4, 2, 1, 2, 4, 5, 1, 4, 3, 1, 2, 4, 4, 5, 4, 4, 4, 4, 1, 1, 4, 2, 1, 2, 4, 5, 1, 4, 3, 1, 2, 4),
  y = c(-3.18649152, -3.25869527, -1.64841563, -2.43408372, -0.09610498, -1.65460665, 1.81631785, -0.53701069, -2.98435867, 0.50109073, 0.20436527, -0.24096772, -1.50765524, -3.44179477, 0.56591236, -3.48124692, -2.36647953, 2.25948062, -1.02782941, -0.82430706, 0.72018109, 1.34004336, 2.25809458, 1.47164831, 3.80850612, 2.25153806, 3.64015287, -0.02881929, -0.39501198, 1.70162601, 0.71358084, 0.96240644, 1.08049450, -0.16098213, 1.07408065, -0.89082004, -0.47072410, 2.76987894, 0.17469733, 1.76498883, -5.71884992, -6.37625509, -4.18192234, -4.96832686, -2.63206375, -4.18432324, 1.03219718, -2.16749889, -6.36580419, -1.71299118, -1.42506617, -2.45314625, -4.88839472, -7.40582661, -2.24368809, -8.03861697, -6.34331587, -0.55011068, -4.42010676, -5.38291640)
)

# did with 
for(est in c("dr", "ipw", "reg")){
  for(p in c(FALSE, TRUE)){
    out <- att_gt(yname="y",
                  tname="time",
                  idname="unit",
                  gname="G",
                  xformla=~x,
                  data=dt,
                  allow_unbalanced_panel = p,
                  bstrap = FALSE,
                  est_method = est)
    agg <- aggte(out, type = "simple")
    message("se for allow_unbalance_panel = ", p, " est = ", est, ": ", agg$overall.se)
    rm(out, agg)
  }
}

# test filter -------

rm(list = ls()); gc()
library(devtools)
load_all()
library(data.table)
library(ggplot2)
cohort_n <- 500

dt <- data.table(cohort = rep(2010:2015, each = 16, cohort_n),
                yr = rep(2005:2020, 6*cohort_n),
                id = rep(1:(cohort_n*6), each = 16)
                )
dt[, randmn := rep(ceiling(runif(cohort_n*6)*10)%%5, each = 16)]
dt[, event := ifelse(id%%7==1, "a1", "a2")]

dt[, relative_t := yr - cohort]
dt[event =="a2", diff_preg_born := cohort+randmn]
dt[, inc := 50000]
dt[event =="a1" & relative_t>=0 , inc:= 10000]
dt[event =="a2" & relative_t>=0, inc:= 12000]

# dt[event =="a2" & randmn>=1 & relative_t>=0, inc:= 12000]
# dt[event =="a2" & randmn>=2 & relative_t>=0, inc:= 12000]
# dt[event =="a2" & randmn>=3 & relative_t>=0, inc:= 12000]
# dt[event =="a2" & randmn>=4 & relative_t>=0, inc:= 12000]
#dt[event =="a2" & yr>=diff_preg_born , inc:= 10000]

dt[, before_birth := yr==cohort-1]
dt[, reg_event_year := as.numeric(cohort)]
dt[!event == "a1", reg_event_year := Inf]


dt[,.(inc = mean(inc)), by = c("cohort", "yr", "before_birth", "event")] |> 
    ggplot(aes(x = yr, color = as.factor(event), group = cohort, size = ifelse(before_birth, 5, 1), y = inc)) +
    geom_point() + facet_wrap(~event+cohort)


# fastdid

table(dt$event, dt$cohort)
fastdid_func <- function( dt){
 dt[, nID := rleid(id)]
 res <- fastdid(dt,#& marr_before_preg==m], #& first_marr<=2015 &first_marr>=2010
                timevar = "yr",
                cohortvar = "reg_event_year",
                unitvar = "nID",
                outcomevar = "inc",
                result_type = "dynamic",
                control_option = "never", allow_unbalance_panel = F,
                validate = T,
                exper = list(filtervar = "before_birth"),
                copy = T)
 print(plot_did_dynamics(res))

 return(res)
}

est_1_event <-data.table()
# for(i in c(0:4)){
res <- fastdid_func( dt)#[event == "a1" |(event=="a2"&randmn==i)]


res <- fastdid(dt[event == "a1"],#& marr_before_preg==m], #& first_marr<=2015 &first_marr>=2010
                timevar = "yr",
                cohortvar = "cohort",
                unitvar = "nID",
                outcomevar = "inc",
                result_type = "dynamic",
                control_option = "notyet", allow_unbalance_panel = F,
                validate = T,
                copy = T)
 print(plot_did_dynamics(res))

res <- fastdid(dt[event == "a2"],#& marr_before_preg==m], #& first_marr<=2015 &first_marr>=2010
                timevar = "yr",
                cohortvar = "cohort",
                unitvar = "nID",
                outcomevar = "inc",
                result_type = "dynamic",
                control_option = "notyet", allow_unbalance_panel = F,
                validate = T,
                copy = T)
print(plot_did_dynamics(res))