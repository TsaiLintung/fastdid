
DiDe <- function(inputdata, varnames, control_group = "all", base_event=-1, min_event=NULL, max_event=NULL, return_data=FALSE, return_ATTs_only=TRUE, parallel_cores=1){

  # set up variable names
  time_name = varnames$time_name
  outcome_name = varnames$outcome_name
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name
  covariate_names = varnames$covariate_names
  cluster_names = varnames$cluster_names
  fixedeffect_names = varnames$fixedeffect_names

  # check packages
  check_parallel = requireNamespace("parallel", quietly=TRUE, warn.conflicts = FALSE)
  check_progress = requireNamespace("progress", quietly=TRUE, warn.conflicts = FALSE)

  # set up and checks
  cohorts = sort(inputdata[, unique(get(cohort_name))])
  nevertreated_exist = sum(is.infinite(cohorts))>0
  cohorts = cohorts[!is.infinite(cohorts)]

  upperbound_outcomeyear = Inf

  if(control_group=="never-treated"){
    if(!nevertreated_exist){
      stop("You specified control_group='never-treated', but there are no never-treated observations. Note: you must code never-treated observations as infinity (Inf).")
    }
  }

  if(control_group=="future-treated"){
    if(length(cohorts) <= 1){
      stop("You specified control_group='future-treated', but there is only one treatment cohort with finite cohort time, so no control units are available based on future-treated observations.")
    }
    upperbound_outcomeyear = max(cohorts)
    cohorts = cohorts[!(cohorts == max(cohorts))]
  }

  # drop cohorts that are useless since they were treated before the base_event
  min_time_actual = inputdata[,min(get(time_name))]
  too_early_cohorts = cohorts - base_event < min_time_actual
  if(sum(too_early_cohorts) > 0){
    warning(sprintf("We cannot provide ATT estimates for cohort %s due to the absence of a base_preperiod in the data.", paste0(cohorts[too_early_cohorts], collapse=",")))
    cohorts_todrop = cohorts[too_early_cohorts]
    cohorts = cohorts[!too_early_cohorts]
    inputdata = inputdata[!(get(cohort_name) %in% cohorts_todrop)]
  }

  # cohort-specific estimation
  get_cohort_results <- function(cc){
    results_this_cohort = data.table()
    data_this_cohort = data.table()
    # set up cohort-specific event times
    times_for_cohort = sort(inputdata[get(cohort_name) == cc, unique(get(time_name))])
    times_for_cohort = times_for_cohort[times_for_cohort < upperbound_outcomeyear]
    event_periods = times_for_cohort - cc
    # drop event times outside of the specified min/max range, for this specific cohort
    if(!is.null(min_event)){
      event_periods = event_periods[event_periods >= min_event]
    }
    if(!is.null(max_event)){
      event_periods = event_periods[event_periods <= max_event]
    }
    # loop DiD over the event times for this cohort
    for(event_postperiod in event_periods){
      res = DiDge(inputdata, cohort_time = cc, event_postperiod = event_postperiod, base_event = base_event,
                  varnames=varnames, control_group = control_group, return_data=TRUE, return_ATTs_only=FALSE)
      results_this_cohort = rbindlist(list(results_this_cohort,res$results))
      # this is key -- collect cohort-event-specific control-treatment-paired data to estimate DiDe later
      data_this_cohort = rbindlist(list(data_this_cohort, res$data))
    }
    # finish
    return(list(results_cohort=results_this_cohort,data_cohort=data_this_cohort))
  }

  # apply the cohort-specific estimation in parallel
  if(!check_parallel){
    the_results = lapply(cohorts,get_cohort_results)
  }
  if(check_parallel){
    if(!check_progress){
      the_results = parallel::mclapply(cohorts,get_cohort_results,mc.cores=parallel_cores)
    }
    if(check_progress){
      the_results = Errorhandle.mclapply(cohorts,get_cohort_results,mc.cores=parallel_cores)
    }
  }

  # extract the results from the output list
  extract_results <- function(index_count,subtype){
    the_results[[index_count]][[subtype]]
  }
  # gather the cohort-year DiDge estimates and SEs
  results_cohort = rbindlist(lapply(1:length(cohorts), extract_results, subtype="results_cohort"))
  results_cohort = results_cohort[order(Cohort,EventTime)]
  # gather the cohort-event-specific control-treatment-paired data to estimate DiDe later
  data_cohort = rbindlist(lapply(1:length(cohorts), extract_results, subtype="data_cohort"))
  data_cohort = data_cohort[order(Cohort,EventTime)]

  # take the average across cohorts
  results_cohort[, Ntreated_event := sum(Ntreated), by="EventTime"]
  results_cohort[, cohort_weights := Ntreated/Ntreated_event]
  results_average = results_cohort[Ntreated > 1 & Ncontrol > 1,
                                   list(
                                     ATTe_simple=sum(cohort_weights * ATTge_simple),
                                     ATTe_SE_simple=sqrt(sum(cohort_weights^2 * ATTge_SE_simple^2)), # these are the SEs if DiDge were uncorrelated across g
                                     Etreated_post=sum(cohort_weights*Etreated_post),
                                     Etreated_pre=sum(cohort_weights*Etreated_pre),
                                     Etreated_SE=sqrt(sum(cohort_weights^2*Etreated_SE^2)),
                                     Econtrol_post=sum(cohort_weights*Econtrol_post),
                                     Econtrol_pre=sum(cohort_weights*Econtrol_pre),
                                     Econtrol_SE=sqrt(sum(cohort_weights^2*Econtrol_SE^2)),
                                     Ntreated=sum(Ntreated),
                                     Ncontrol=sum(Ncontrol)
                                   ),
                                   list(EventTime,BaseEvent)][order(EventTime,BaseEvent)]

  # collect the SEs that account for correlation in DiDge across g, merge them to the results
  ATTe_SEs = DiD_getSEs_EventTime(data_cohort=data_cohort,varnames=varnames,base_event=base_event)

  results_average = merge(results_average, ATTe_SEs, by="EventTime")[order(EventTime,BaseEvent)]

  # clean up the results
  ordered_names = c("EventTime","BaseEvent",
                    "ATTe","ATTe_SE","ATTe_simple","ATTe_SE_simple",
                    "Etreated_post","Etreated_pre","Etreated_SE",
                    "Econtrol_post","Econtrol_pre","Econtrol_SE",
                    "Ntreated","Ncontrol")
  results_average = results_average[,.SD,.SDcols=ordered_names]

  if(return_data){
    return(list(results_cohort=results_cohort, results_average=results_average, data_cohort=data_cohort))
  }

  if(return_ATTs_only){
    results_cohort=results_cohort[,.SD,.SDcols=c("Cohort", "EventTime", "BaseEvent", "CalendarTime", "ATTge", "ATTge_SE", "Ncontrol", "Ntreated")]
    results_average=results_average[,.SD,.SDcols=c("EventTime", "BaseEvent", "ATTe", "ATTe_SE", "Ncontrol", "Ntreated")]
    return(list(results_cohort=results_cohort, results_average=results_average))
  }
  return(list(results_cohort=results_cohort, results_average=results_average))
}


#' Combine DiD estimates across cohorts and event times.
#'
#' @description
#' Estimate DiD for all possible cohorts and event time pairs (g,e), as well as the average across cohorts for each event time (e).
#'
#' @param inputdata A data.table.
#' @param varnames A list of the form varnames = list(id_name, time_name, outcome_name, cohort_name), where all four arguments of the list must be a character that corresponds to a variable name in inputdata.
#' @param control_group There are three possibilities: control_group="never-treated" uses the never-treated control group only; control_group="future-treated" uses those units that will receive treatment in the future as the control group; and control_group="all" uses both the never-treated and the future-treated in the control group. Default is control_group="all".
#' @param base_event This is the base pre-period that is normalized to zero in the DiD estimation. Default is base_event=-1.
#' @param min_event This is the minimum event time (e) to estimate. Default is NULL, in which case, no minimum is imposed.
#' @param max_event This is the maximum event time (e) to estimate. Default is NULL, in which case, no maximum is imposed.
#' @param Esets If a list of sets of event times is provided, it will loop over those sets, computing the average ATT_e across event times e. Default is NULL.
#' @param return_ATTs_only Return only the ATT estimates and sample sizes. Default is TRUE.
#' @param parallel_cores Number of cores to use in parallel processing. If greater than 1, it will try to run library(parallel), so the "parallel" package must be installed. Default is 1.
#' @return A list with two components: results_cohort is a data.table with the DiDge estimates (by event e and cohort g), and results_average is a data.table with the DiDe estimates (by event e, average across cohorts g). If the Esets argument is specified, a third component called results_Esets will be included in the list of output.
#' @examples
#' # simulate some data
#' simdata = SimDiD(sample_size=200, ATTcohortdiff = 2)$simdata
#'
#' # define the variable names as a list()
#' varnames = list()
#' varnames$time_name = "year"
#' varnames$outcome_name = "Y"
#' varnames$cohort_name = "cohort"
#' varnames$id_name = "id"
#'
#' # estimate the ATT for all cohorts at event time 1 only
#' DiD(simdata, varnames, min_event=1, max_event=1)
#'
#' @export
DiD <- function(inputdata, varnames, control_group = "all", base_event=-1, min_event=NULL, max_event=NULL, Esets=NULL, return_ATTs_only=TRUE, parallel_cores=1){
  # case without averaging event sets
  if(is.null(Esets)){
    results = DiDe(inputdata=inputdata, varnames=varnames, control_group=control_group, base_event=base_event, min_event=min_event, max_event=max_event, return_data=FALSE, return_ATTs_only=return_ATTs_only, parallel_cores=parallel_cores)
    return(results)
  }
  # case with averaging event sets
  if(!is.null(Esets)){
    results = DiDe(inputdata=inputdata, varnames=varnames, control_group=control_group, base_event=base_event, min_event=min_event, max_event=max_event, return_data=TRUE, parallel_cores=parallel_cores)
    data_cohort = results$data_cohort
    results_Esets = data.table()
    for(Eset in Esets){
      results_Esets = rbindlist(list(
        results_Esets,
        DiD_getSEs_multipleEventTimes(data_cohort,varnames,Eset=Eset,min_event=min_event,max_event=max_event)
      ))
    }
    results$results_Esets = results_Esets
    results$data_cohort = NULL

    if(return_ATTs_only){
      results$results_cohort=results$results_cohort[,.SD,.SDcols=c("Cohort", "EventTime", "BaseEvent", "CalendarTime", "ATTge", "ATTge_SE", "Ncontrol", "Ntreated")]
      results$results_average=results$results_average[,.SD,.SDcols=c("EventTime", "BaseEvent", "ATTe", "ATTe_SE", "Ncontrol", "Ntreated")]
    }

    return(results)
  }
}

