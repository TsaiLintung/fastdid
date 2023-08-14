#' DiD data simulator with staggered treatment.
#' @description
#' Simulate data from the model Y_it =  alpha_i + mu_t + ATT*(t >= G_i) + epsilon_it,
#' where i is individual, t is year, and G_i is the cohort.
#' The ATT formula is ATTat0 + EventTime*ATTgrowth + \*cohort_counter\*ATTcohortdiff,
#' where cohort_counter is the order of treated cohort (first, second, etc.).
#' @param seed Set the random seed. Default is seed=1.
#' @param sample_size Number of individuals. Default is sample_size=100.
#' @param cohorts Vector of years at which treatment onset occurs. Default is cohorts=c(2007,2010,2012).
#' @param ATTat0 Treatment effect at event time 0. Default is 1.
#' @param ATTgrowth Increment in the ATT for each event time after 0. Default is 1.
#' @param ATTcohortdiff Incrememnt in the ATT for each cohort. Default is 0.5.
#' @param anticipation Number of years prior to cohort to allow 50% treatment effects. Default is anticipation=0.
#' @param minyear Minimum calendar year to include in the data. Default is minyear=2003.
#' @param maxyear Maximum calendar year to include in the data. Default is maxyear=2013.
#' @param idvar Variance of individual fixed effects (alpha_i). Default is idvar=1.
#' @param yearvar Variance of year effects (mu_i). Default is yearvar=1.
#' @param shockvar Variance of idiosyncratic shocks (epsilon_it). Default is shockvar=1.
#' @param indivAR1 Each individual's shocks follow an AR(1) process. Default is FALSE.
#' @param time_covars Add 2 time-varying covariates, called "X1" and "X2". Default is FALSE.
#' @param clusters Add 10 randomly assigned clusters, with cluster-specific AR(1) shocks. Default is FALSE.
#' @param markets Add 10 randomly assigned markets, with market-specific shocks that are systematically greater for markets that are treated earlier. Default is FALSE.
#' @param randomNA If TRUE, randomly assign the outcome variable with missing values (NA) in some cases. Default is FALSE.
#' @param missingCohorts If set to a particular cohort (or vector of cohorts), all of the outcomes for that cohort at event time -1 will be set to missing. Default is NULL.
#' @return A list with two data.tables.
#' The first data.table is simulated data with variables (id, year, cohort, Y), where Y is the outcome variable.
#' The second data.table contains the true ATT values, both at the (event,cohort) level and by event averaging across cohorts.
#' @examples
#' # simulate data with default options
#' SimDiD()
#' @export
SimDiD <- function(seed=1,sample_size=100, cohorts=c(2007,2010,2012), ATTat0=1, ATTgrowth=1, ATTcohortdiff=0.5, anticipation=0, minyear=2003, maxyear=2013, idvar=1, yearvar=1, shockvar=1, indivAR1=FALSE, time_covars=FALSE, clusters=FALSE, markets=FALSE, randomNA=FALSE, missingCohorts=NULL){
  # seed=1; sample_size=1000; cohorts=c(2007,2010,2012); ATTat0=1; ATTgrowth=1; ATTcohortdiff=0.5; anticipation=0; minyear=2003; maxyear=2013; idvar=1; yearvar=1; shockvar=1; indivAR1=FALSE; time_covars=FALSE; bin_covars_nobias=FALSE; bin_covars_bias=FALSE
  set.seed(seed)

  # create id-by-year data
  simdata = setDT(expand.grid(id=1:sample_size,year=minyear:maxyear))[order(id,year)]

  # simulate unobservables
  simdata[, shock := rnorm(nrow(simdata),sd=shockvar)]
  if(indivAR1){
    numyears = simdata[,length(unique(year))]
    simdata[, shock := as.numeric(arima.sim(list(order=c(shockvar,0,0), ar=.5), n=numyears)), id]
  }
  simdata[, individualFE := rnorm(1,sd=idvar), id]
  simdata[, yearFE := rnorm(1,sd=yearvar), year]
  simdata[, individualFE_ecdf := ecdf(individualFE)(individualFE)]

  # simulate treatment cohort
  simcohorts = c(cohorts,Inf)
  cutoffs = rep(1,length(simcohorts))/length(simcohorts)
  cutoffs = c(0,cumsum(cutoffs))
  for(ii in 1:(length(cutoffs)-1)){
    simdata[individualFE_ecdf >= cutoffs[ii] & individualFE_ecdf <= cutoffs[ii+1], cohort := simcohorts[ii]]
  }
  simdata[, event := year - cohort]

  # simulate ATT
  simdata[, cohort_counter := 0.0]
  for(ii in 1:length(cohorts)){
    simdata[cohort==cohorts[ii], cohort_counter := (ii-1)]
  }
  simdata[, ATT := 0.0]
  simdata[year >= cohort, ATT := (ATTat0 + event*ATTgrowth + cohort_counter*ATTcohortdiff)]
  simdata[year < cohort & year >= (cohort-anticipation), ATT := 0.5*ATTat0]

  # extract true ATT
  true_ATTge = simdata[(year >= (cohort-anticipation)),list(ATTge = mean(ATT)), list(cohort, event)][order(cohort,event)][, cohort := as.character(cohort)]
  true_ATTe = simdata[(year >= (cohort-anticipation)),list(cohort="Average",ATTge = mean(ATT)), list(event)][order(event)]
  true_ATT = rbindlist(list(true_ATTge,true_ATTe),use.names=T)

  # simulate outcome
  keep_vars = c("id","year","cohort","Y")
  simdata[, Y := 10 + individualFE + yearFE + ATT + shock]
  if(time_covars){
    simdata[, X1 := year - mean(year) + rnorm(nrow(simdata))] # this causes no issues, since it differs out on average
    simdata[, X2 := 0.0]
    simdata[event > 0, X2 := event]
    simdata[, X2 := log(1 + cohort_counter + X2) + rnorm(nrow(simdata))] # this causes problems because there is an event-cohort interaction
    simdata[, covariate_term := -1.0*X1 + 1.0*X2]
    simdata[, Y := Y + covariate_term]
    keep_vars = c(keep_vars, "X1", "X2")
  }
  if(clusters | markets){
    simdata[, cohort_counter := 0.0]
    for(ii in 1:length(cohorts)){
      simdata[cohort==cohorts[ii], cohort_counter := (1+length(cohorts)-ii)]
    }
    id_cohort = unique(simdata[,list(id,cohort_counter)])
    bin_assignments = NULL
    for(bin in 1:10){
      bin_assignment = copy(id_cohort)
      bin_assignment[, bin_counter := bin]
      bin_assignments = rbindlist(list(bin_assignments, bin_assignment))
    }
    # we assume low-numbered markets (e.g. market=1, market=2) tend to receive treatment earlier than high-numbered markets
    bin_assignments[, bin_index := exp(0.1*cohort_counter*bin_counter)]
    bin_assignments[, bin_prob := sum(bin_index), id]
    bin_assignments[, bin_prob := bin_index/bin_prob]
    bin_draws = NULL
    for(ii in 1:sample_size){
      p = bin_assignments[id==ii, bin_prob]
      bin_draw = sample(1:10, 1, prob = p, replace = TRUE)
      bin_draws = rbindlist(list(bin_draws,data.table(id=ii,bin=bin_draw)))
    }
    simdata = merge(simdata, bin_draws, by=c("id"))
    if(clusters){
      # the shocks are randomly assigned to the bin-years, so they do not induce bias
      numyears = simdata[,length(unique(year))]
      binyears = unique(simdata[,list(bin,year)])[order(bin,year)]
      binyears[, bin_shock := 0.0]
      binyears[, bin_shock := arima.sim(list(order=c(1,0,0), ar=.5), n=numyears), list(bin)] # each bin gets an AR(1) process with persistence 0.5
      simdata = merge(simdata, binyears, by=c("bin","year"))
      keep_vars = c(keep_vars, "cluster")
      setnames(simdata,'bin','cluster')
    }
    if(markets){
      # the shocks are systematically more positive for bins that are treated earlier
      simdata[, bin_shock := - cohort_counter*bin/10 ]
      keep_vars = c(keep_vars, "market")
      setnames(simdata,'bin','market')
    }
    simdata[, Y := Y + bin_shock]
  }
  if(randomNA){
    simdata[, change_to_NA := runif(nrow(simdata)) < 0.1] # 10% selected randomly
    simdata[change_to_NA==TRUE, Y := NA]
  }
  if(!is.null(missingCohorts)){
    for(cc in missingCohorts){
      simdata[cohort==cc & year==(cc-1), Y := NA]
    }
  }
  simdata = simdata[order(id,year)]
  simdata = simdata[,.SD,.SDcols=keep_vars]
  return(list(simdata=simdata,true_ATT=true_ATT))
}
