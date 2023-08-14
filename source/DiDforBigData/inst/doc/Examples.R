## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=T, eval=FALSE, message=FALSE----------------------------------------
#  install.packages("DiDforBigData")

## ----echo=T, eval=FALSE, message=FALSE----------------------------------------
#  devtools::install_github("setzler/DiDforBigData")

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
library(DiDforBigData)

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
sim = SimDiD(seed=123, sample_size = 1000)
simdata = sim$simdata
print(simdata)

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
print(sim$true_ATT[cohort==2007])

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
varnames = list()
varnames$time_name = "year" 
varnames$outcome_name = "Y"
varnames$cohort_name = "cohort"
varnames$id_name = "id"

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
did = DiD(inputdata = simdata, varnames = varnames, min_event = -3, max_event=5)

print(did$results_cohort[Cohort==2007])

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
did = DiD(inputdata = simdata, varnames = varnames, 
          control_group = "never-treated", min_event = -3, max_event=5)

print(did$results_cohort[Cohort==2007])

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
did = DiD(inputdata = simdata, varnames = varnames, 
          control_group = "future-treated", min_event = -3, max_event=5)

print(did$results_cohort[Cohort==2007])

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
sim = SimDiD(seed=123, sample_size = 200, anticipation = 2)
simdata = sim$simdata
print(simdata)

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
print(sim$true_ATT[cohort=="Average"])

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
varnames = list()
varnames$time_name = "year" 
varnames$outcome_name = "Y"
varnames$cohort_name = "cohort"
varnames$id_name = "id"

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
did = DiD(inputdata = simdata, varnames = varnames, min_event = -3, max_event=3)

print(did$results_average)

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
did = DiD(inputdata = simdata, varnames = varnames, 
          base_event = -3, min_event = -3, max_event=3)

print(did$results_average) 

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------

sim = SimDiD(sample_size=1000, time_covars=TRUE)
simdata = sim$simdata
print(simdata)

print(sim$true_ATT[cohort==2007])

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
varnames = list()
varnames$time_name = "year" 
varnames$outcome_name = "Y"
varnames$cohort_name = "cohort"
varnames$id_name = "id"

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
did = DiD(inputdata = simdata, varnames = varnames, min_event = -3, max_event=5)

print(did$results_cohort[Cohort==2007])

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
varnames$covariate_names = c("X1","X2")

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
did = DiD(inputdata = simdata, varnames = varnames, min_event = -3, max_event=5)

print(did$results_cohort[Cohort==2007])

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------

sim = SimDiD(sample_size=1000, clusters = TRUE)
simdata = sim$simdata
print(simdata)

print(sim$true_ATT[cohort=="Average"])

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
varnames = list()
varnames$time_name = "year" 
varnames$outcome_name = "Y"
varnames$cohort_name = "cohort"
varnames$id_name = "id"

## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
did = DiD(inputdata = simdata, varnames = varnames, min_event = -1, max_event=3)

print(did$results_average)


## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------

varnames$cluster_names = "cluster" 

did = DiD(inputdata = copy(simdata), varnames = varnames, min_event = -1, max_event=3)

print(did$results_average)


## ----echo=TRUE, eval = TRUE, message=FALSE------------------------------------
sim = SimDiD(seed=123, sample_size = 1000)
simdata = sim$simdata

## ----echo=TRUE, eval=FALSE, message=FALSE-------------------------------------
#  varnames = list()
#  varnames$time_name = "year"
#  varnames$outcome_name = "Y"
#  varnames$cohort_name = "cohort"
#  varnames$id_name = "id"

## ----echo=TRUE, eval=FALSE, message=FALSE-------------------------------------
#  did = DiD(inputdata = copy(simdata), varnames = varnames, min_event = -1, max_event=3)

## ----echo=TRUE, eval=FALSE, message=FALSE-------------------------------------
#  did = DiD(inputdata = copy(simdata), varnames = varnames, min_event = -1, max_event=3, parallel_cores = 2)

## ----echo=TRUE, eval=TRUE, message=FALSE--------------------------------------
sim = SimDiD(seed=123, sample_size = 1000)
simdata = sim$simdata

## ----echo=TRUE, eval=TRUE, message=FALSE--------------------------------------
varnames = list()
varnames$time_name = "year" 
varnames$outcome_name = "Y"
varnames$cohort_name = "cohort"
varnames$id_name = "id"

## ----echo=TRUE, eval=TRUE, message=FALSE--------------------------------------
did = DiD(inputdata = copy(simdata), varnames = varnames, min_event = -1, max_event=3, Esets = list(c(1,2,3), c(1,3)))

print(did)

