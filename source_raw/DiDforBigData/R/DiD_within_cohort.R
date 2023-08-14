
DiDge_main <- function(inputdata, varnames, cohort_time, event_postperiod, base_event = -1, control_group = "all", return_data=FALSE){

  # set up variable names
  time_name = varnames$time_name
  outcome_name = varnames$outcome_name
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name
  covariate_names = varnames$covariate_names
  cluster_names = varnames$cluster_names
  fixedeffect_names = varnames$fixedeffect_names
  keep_vars = c(outcome_name,covariate_names)
  all_keep_vars = c(outcome_name,covariate_names,cluster_names,fixedeffect_names)
  all_keep_vars = all_keep_vars[all_keep_vars != id_name]

  # check if fixest is available
  check_fixest = requireNamespace("fixest", quietly=TRUE, warn.conflicts = FALSE)
  if(!is.null(fixedeffect_names)){
    if(!check_fixest){
      stop("Since varnames$fixedeffect_names is non-missing, you must install the 'fixest' package, which estimates fixed-effects.")
    }
  }

  # prepare time periods
  pre_time = cohort_time + base_event
  post_time = cohort_time + event_postperiod
  time_set = sort(inputdata[get(cohort_name) == cohort_time, unique(get(time_name))])
  if(!(pre_time %in% time_set)){
    stop(print(sprintf("error: for cohort %s, preperiod %s is unavailable.",cohort_time,base_event)))
  }
  if(!(post_time %in% time_set)){
    stop(print(sprintf("error: for cohort %s, postperiod %s is unavailable.",cohort_time,event_postperiod)))
  }

  # restrict to pre and post time periods of interest, then restrict to observations present in both time periods
  inputdata = inputdata[get(time_name)==pre_time | get(time_name)==post_time]
  inputdata = inputdata[,.SD,.SDcols=c(id_name, time_name, cohort_name, all_keep_vars)]
  nn0 = nrow(inputdata)
  inputdata = na.omit(inputdata) # remove any rows with missing values
  nn1 = nrow(inputdata)
  if(nn1 < nn0){
    warning(sprintf("%s out of %s observations dropped due to missing values",(nn0-nn1),nn0))
  }
  inputdata[, nobs := .N, by=id_name]
  inputdata = inputdata[nobs==2 | (get(time_name)==pre_time & get(time_name)==post_time)] # the second condition keeps the base period
  gc()

  # define treated data and get means
  treated_data_prepost = merge(
    inputdata[(get(cohort_name) == cohort_time) & (get(time_name) == pre_time), .SD, .SDcols=c(id_name,all_keep_vars)],
    inputdata[(get(cohort_name) == cohort_time) & (get(time_name) == post_time), .SD, .SDcols=c(id_name,keep_vars)],
    by=c(id_name)
  )
  names(treated_data_prepost) = gsub("\\.x","_pre",names(treated_data_prepost))
  names(treated_data_prepost) = gsub("\\.y","_post",names(treated_data_prepost))
  setnames(treated_data_prepost, paste0(keep_vars,"_pre"), paste0("treated_",keep_vars,"_pre"))
  setnames(treated_data_prepost, paste0(keep_vars,"_post"), paste0("treated_",keep_vars,"_post"))
  treated_data_prepost[, treated := 1.0]

  # define control data and get means
  control_data_prepost = NULL
  if(control_group == "all"){
    control_data_prepost = merge(
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & (get(time_name) == pre_time), .SD, .SDcols=c(id_name,all_keep_vars)],
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & (get(time_name) == post_time), .SD, .SDcols=c(id_name,keep_vars)],
      by=c(id_name)
    )
  }
  if(control_group == "never-treated"){
    control_data_prepost = merge(
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.infinite(get(cohort_name)) & (get(time_name) == pre_time), .SD, .SDcols=c(id_name,all_keep_vars)],
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.infinite(get(cohort_name)) & (get(time_name) == post_time), .SD, .SDcols=c(id_name,keep_vars)],
      by=c(id_name)
    )
  }
  if(control_group == "future-treated"){
    control_data_prepost = merge(
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.finite(get(cohort_name)) & (get(time_name) == pre_time), .SD, .SDcols=c(id_name,all_keep_vars)],
      inputdata[(get(cohort_name) > max(post_time,cohort_time)) & is.finite(get(cohort_name)) & (get(time_name) == post_time), .SD, .SDcols=c(id_name,keep_vars)],
      by=c(id_name)
    )
  }
  names(control_data_prepost) = gsub("\\.x","_pre",names(control_data_prepost))
  names(control_data_prepost) = gsub("\\.y","_post",names(control_data_prepost))
  setnames(control_data_prepost, paste0(keep_vars,"_pre"), paste0("control_",keep_vars,"_pre"))
  setnames(control_data_prepost, paste0(keep_vars,"_post"), paste0("control_",keep_vars,"_post"))
  control_data_prepost[, treated := 0.0]

  # levels
  treated_pre_outcome = paste0("treated_",outcome_name,"_pre")
  treated_post_outcome = paste0("treated_",outcome_name,"_post")
  treated_pre = treated_data_prepost[, list(Ntreated_pre=sum(!is.na(get(treated_pre_outcome))),Etreated_pre=mean(get(treated_pre_outcome)))]
  treated_post = treated_data_prepost[, list(Ntreated_post=sum(!is.na(get(treated_post_outcome))),Etreated_post=mean(get(treated_post_outcome)),Etreated_var=var(get(treated_post_outcome)))]
  treated_results = cbind(treated_pre,treated_post)
  treated_results$treated_diff_var = treated_data_prepost[, var(get(treated_post_outcome) - get(treated_pre_outcome))]
  treated_results[, Ntreated_pre := NULL]
  setnames(treated_results,"Ntreated_post","Ntreated")

  control_pre_outcome = paste0("control_",outcome_name,"_pre")
  control_post_outcome = paste0("control_",outcome_name,"_post")
  control_pre = control_data_prepost[, list(Ncontrol_pre=sum(!is.na(get(control_pre_outcome))),Econtrol_pre=mean(get(control_pre_outcome)))]
  control_post = control_data_prepost[, list(Ncontrol_post=sum(!is.na(get(control_post_outcome))),Econtrol_post=mean(get(control_post_outcome)),Econtrol_var=var(get(control_post_outcome)))]
  control_results = cbind(control_pre,control_post)
  control_results$control_diff_var = control_data_prepost[, var(get(control_post_outcome) - get(control_pre_outcome))]
  control_results[, Ncontrol_pre := NULL]
  setnames(control_results,"Ncontrol_post","Ncontrol")

  # set up results
  results = cbind(control_results,treated_results)
  results[, Cohort := cohort_time]
  results[, Preperiod := pre_time]
  results[, CalendarTime := post_time]
  results[, BaseEvent := base_event]
  results[, EventTime := event_postperiod]
  results[, Econtrol_SE := sqrt(Econtrol_var/(Ncontrol))]
  results[, Etreated_SE := sqrt(Etreated_var/(Ntreated))]
  results[, pred_Etreated_post := Etreated_pre + (Econtrol_post - Econtrol_pre)]
  results[, ATTge_simple := (Etreated_post - pred_Etreated_post)]
  results[, ATTge_SE_simple := sqrt(treated_diff_var/(Ntreated) + control_diff_var/(Ncontrol))]

  results_variables_order = c("Cohort","EventTime","BaseEvent","CalendarTime",
                              "ATTge","ATTge_SE","ATTge_simple","ATTge_SE_simple",
                              "Econtrol_pre","Econtrol_post","Econtrol_SE",
                              "Etreated_pre","Etreated_post","Etreated_SE",
                              "Ncontrol","Ntreated")

  # set up differences
  data_prepost = NULL
  names(treated_data_prepost) = gsub("treated_","",names(treated_data_prepost))
  names(control_data_prepost) = gsub("control_","",names(control_data_prepost))
  data_prepost = rbindlist(list(treated_data_prepost,control_data_prepost))
  for(ii in keep_vars){
    data_prepost[, (paste0(ii,"_diff")) := get(paste0(ii,"_post")) - get(paste0(ii,"_pre"))]
  }
  data_prepost = data_prepost[,.SD,.SDcols=c(id_name,"treated",paste0(keep_vars,"_diff"),cluster_names,fixedeffect_names)]


  # the case in which ATT is mechanically 0
  if(pre_time==post_time){
    results[, ATTge := NA]
    results[, ATTge_SE := NA]
    if(results$Ntreated > 0 & results$Ncontrol > 0){
      results[, ATTge := 0.0]
      results[, ATTge_SE := 0.0]
    }
    results = results[,.SD,.SDcols=results_variables_order]
    results = results[order(Cohort,EventTime)]
    if(!return_data){
      return(results)
    }
    if(return_data){
      data_prepost[, Cohort := cohort_time]
      data_prepost[, EventTime := event_postperiod]
      export_vars = c(id_name,"Cohort","EventTime","treated",paste0(keep_vars,"_diff"),cluster_names)
      data_prepost = data_prepost[,.SD,.SDcols=export_vars]
      return(list(results=results,data_prepost=data_prepost))
    }
  }

  # execute the regression
  OLSlm = NULL
  if(is.null(fixedeffect_names)){
    OLSformula = paste0(paste0(outcome_name,"_diff"), paste0(" ~ treated"))
    if(!is.null(covariate_names)){ # reg formula, with covariates
      OLSformula = paste0(OLSformula, " + ", paste0(paste0(covariate_names,"_diff"),collapse=" + "))
    }
    if(check_fixest){ # prefer feols() if installed
      OLSlm = fixest::feols(as.formula(OLSformula),data=data_prepost)
    }
    if(!check_fixest){ # use lm() if feols() not installed
      OLSlm = lm(as.formula(OLSformula),data=data_prepost)
    }
  }
  if(!is.null(fixedeffect_names)){ # fixed effects
    OLSformula = paste0(paste0(outcome_name,"_diff"), paste0(" ~ -1 + treated"))
    if(!is.null(covariate_names)){ # reg formula, with covariates
      OLSformula = paste0(OLSformula, " + ", paste0(paste0(covariate_names,"_diff"),collapse=" + "))
    }
    OLSformula = paste0(OLSformula, " | ", paste0(fixedeffect_names, collapse=" + "))
    OLSlm = fixest::feols(as.formula(OLSformula),data=data_prepost)
  }
  newATT = as.numeric(OLSlm$coefficients["treated"])

  # check if the treated coefficient is missing
  if(is.na(newATT)){
    results[, ATTge := NA]
    results[, ATTge_SE := NA]
  }

  # if non-missing ATT, get standard error
  if(!is.na(newATT)){
    results[, ATTge := newATT]
    if(!is.null(cluster_names)){
      data_prepost <<- copy(data_prepost) # due to a well-known scoping bug in R's base lm.predict that no one will fix despite years of requests, this redundancy is necessary!
      CLformula = as.formula(paste0(" ~ ", paste0(cluster_names, collapse=" + ")))
      if(!check_fixest){
        OLSvcov = vcovCL(OLSlm, cluster = CLformula, type = "HC1")
      }
      if(check_fixest){
        OLSvcov = vcov(OLSlm, cluster = cluster_names, type = "HC1")
      }
    }
    if(is.null(cluster_names)){
      if(!check_fixest){
        OLSvcov = vcovHC(OLSlm, type = "HC1")
      }
      if(check_fixest){
        OLSvcov = vcov(OLSlm, vcov="hetero")
      }
    }
    OLSvcov = OLSvcov["treated", "treated"]
    newATTSE = sqrt(as.numeric(OLSvcov))
    results[, ATTge_SE := newATTSE]
  }

  # combine means into an output table
  results = results[,.SD,.SDcols=results_variables_order]
  results = results[order(Cohort,EventTime)]
  if(return_data){
    data_prepost[, Cohort := cohort_time]
    data_prepost[, EventTime := event_postperiod]
    export_vars = unique(c(id_name,"Cohort","EventTime","treated",paste0(keep_vars,"_diff"),cluster_names,fixedeffect_names))
    data_prepost = data_prepost[,.SD,.SDcols=export_vars]
    return(list(results=results,data_prepost=data_prepost))
  }
  return(results)
}

#' Estimate DiD for a single cohort (g) and a single event time (e).
#'
#' @param inputdata A data.table.
#' @param varnames A list of the form varnames = list(id_name, time_name, outcome_name, cohort_name), where all four arguments of the list must be a character that corresponds to a variable name in inputdata.
#' @param cohort_time The treatment cohort of reference.
#' @param event_postperiod Number of time periods after the cohort time at which to estimate the DiD.
#' @param control_group There are three possibilities: control_group="never-treated" uses the never-treated control group only; control_group="future-treated" uses those units that will receive treatment in the future as the control group; and control_group="all" uses both the never-treated and the future-treated in the control group. Default is control_group="all".
#' @param base_event This is the base pre-period that is normalized to zero in the DiD estimation. Default is base_event=-1.
#' @param return_data If true, this returns the treated and control differenced data. Default is FALSE.
#' @param return_ATTs_only Return only the ATT estimates and sample sizes. Default is TRUE.
#' @return A single-row data.table() containing the estimates and various statistics such as sample size. If `return_data=TRUE`, it instead returns a list in which the `data_prepost` entry is the previously-mentioned single-row data.table(), and the other argument `data_prepost`  contains the constructed data that should be provided to OLS.
#' @examples
#' # simulate some data
#' simdata = SimDiD(sample_size=200)$simdata
#'
#' # define the variable names as a list()
#' varnames = list()
#' varnames$time_name = "year"
#' varnames$outcome_name = "Y"
#' varnames$cohort_name = "cohort"
#' varnames$id_name = "id"
#'
#' # estimate the ATT for cohort 2007 at event time 1
#' DiDge(simdata, varnames, cohort_time=2007, event_postperiod=1)
#'
#' # change the base period to -3
#' DiDge(simdata, varnames, base_event=-3, cohort_time=2007, event_postperiod=1)
#'
#' # use only the never-treated control group
#' DiDge(simdata, varnames, control_group = "never-treated", cohort_time=2007, event_postperiod=1)
#'
#' @export
DiDge <- function(inputdata, varnames, cohort_time, event_postperiod, base_event = -1, control_group = "all", return_data=FALSE, return_ATTs_only=TRUE){

  results = DiDge_main(inputdata=inputdata, varnames=varnames, cohort_time=cohort_time, event_postperiod=event_postperiod, base_event=base_event, control_group = control_group, return_data=return_data)

  if(return_ATTs_only){
    return(results[,.SD,.SDcols=c("Cohort","EventTime","BaseEvent","CalendarTime","ATTge","ATTge_SE","Ncontrol","Ntreated")])
  }

  return(results)
}

