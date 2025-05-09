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