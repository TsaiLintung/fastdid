% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastdid.R
\name{fastdid}
\alias{fastdid}
\title{Fast Staggered DID Estimation}
\usage{
fastdid(
  data,
  timevar,
  cohortvar,
  unitvar,
  outcomevar,
  control_option = "both",
  result_type = "group_time",
  balanced_event_time = NULL,
  control_type = "ipw",
  allow_unbalance_panel = FALSE,
  boot = FALSE,
  biters = 1000,
  weightvar = NULL,
  clustervar = NULL,
  covariatesvar = NULL,
  copy = TRUE,
  validate = TRUE
)
}
\arguments{
\item{data}{A data table containing the panel data.}

\item{timevar}{The name of the time variable.}

\item{cohortvar}{The name of the cohort (group) variable.}

\item{unitvar}{The name of the unit (id) variable.}

\item{outcomevar}{The name of the outcome variable.}

\item{control_option}{The control units used for the DiD estimates. Default is "both".}

\item{result_type}{A character string indicating the type of result to be returned. Default is "group_time".}

\item{balanced_event_time}{A numeric scalar that indicates the max event time to balance the cohort composition, only meaningful when result_type == "dynamic". Default is NULL}

\item{control_type}{The method for controlling for covariates. "ipw" for inverse probability weighting, "reg" for outcome regression, or "dr" for doubly-robust}

\item{allow_unbalance_panel}{Whether allow unbalance panel as input (if false will coerce the dataset to a balanced panel). Default is FALSE}

\item{boot}{Logical, indicating whether bootstrapping should be performed. Default is FALSE.}

\item{biters}{The number of bootstrap iterations. Only relevant if boot = TRUE. Default is 1000.}

\item{weightvar}{The name of the weight variable (optional).}

\item{clustervar}{The name of the cluster variable, can only be used when boot == TRUE (optional).}

\item{covariatesvar}{A character vector containing the names of covariate variables (optional).}

\item{copy}{whether to copy the dataset before processing, set to false to speed up the process, but the input data will be altered.}

\item{validate}{whether to validate the dataset before processing.}
}
\value{
A data.table containing the estimated treatment effects and standard errors.
}
\description{
Performs Difference-in-Differences (DID) estimation fast.
}
\examples{
# simulated data
simdt <- sim_did(1e+03, 10, cov = "cont", second_cov = TRUE, second_outcome = TRUE)
dt <- simdt$dt

#basic call
result <- fastdid(dt, timevar = "time", cohortvar = "G", 
                  unitvar = "unit", outcomevar = "y",  
                  result_type = "group_time")

#control for covariates
result2 <- fastdid(dt, timevar = "time", cohortvar = "G", 
                   unitvar = "unit", outcomevar = "y",  
                   result_type = "group_time",
                   covariatesvar = c("x", "x2"))
                  
#bootstrap and clustering
result3 <- fastdid(dt, timevar = "time", cohortvar = "G", 
                   unitvar = "unit", outcomevar = "y",  
                   result_type = "group_time",
                   boot = TRUE, clustervar = "x")

#estimate for multiple outcomes
result4 <- fastdid(dt, #the dataset
                   timevar = "time", cohortvar = "G", unitvar = "unit", 
                   outcomevar = c("y", "y2"), #name of the outcome columns
                   result_type = "group_time") 

}
\keyword{computation}
\keyword{data}
\keyword{did}
\keyword{difference-in-differences}
\keyword{estimation}
\keyword{fast}
\keyword{panel}
