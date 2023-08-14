## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=T, eval=FALSE, message=FALSE----------------------------------------
#  devtools::install_github("setzler/DiDforBigData")

## ----echo=T, eval=T, message=FALSE--------------------------------------------
library(DiDforBigData)

## ----echo=T, eval=T, message=FALSE--------------------------------------------
sim = SimDiD(sample_size = 400, seed=123)

# true ATTs in the simulation
print(sim$true_ATT)

# simulated data
simdata = sim$simdata
print(simdata)

## ----echo=T, eval=T, message=FALSE--------------------------------------------
varnames = list()
varnames$time_name = "year" 
varnames$outcome_name = "Y"
varnames$cohort_name = "cohort"
varnames$id_name = "id"

## ----echo=T, eval=T, message=FALSE--------------------------------------------
did_2010 = DiDge(inputdata = simdata, varnames = varnames, 
             cohort_time = 2010, event_postperiod = 3)

print(did_2010)

## ----echo=T, eval=T, message=FALSE--------------------------------------------
did_all = DiD(inputdata = simdata, varnames = varnames, min_event = -3, max_event = 5)

## ----echo=T, eval=T, message=FALSE--------------------------------------------
print(did_all$results_average)

## ----echo=T, eval=T, message=FALSE--------------------------------------------
print(did_all$results_cohort[EventTime==1 | EventTime==2])

## ----echo=T, eval=T, message=FALSE--------------------------------------------
did_all = DiD(inputdata = simdata, varnames = varnames, min_event = -3, max_event = 5, 
              Esets = list(c(1,2), c(1,2,3)))

## ----echo=T, eval=T, message=FALSE--------------------------------------------
print(did_all$results_Esets)

