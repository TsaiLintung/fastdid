
DiD_getSEs_EventTime <- function(data_cohort,varnames,base_event){

  # one should always cluster on id_name in a stacked regression in which the same unit can appear multiple times
  varnames$cluster_names = unique(c(varnames$cluster_names, varnames$id_name))

  # check that there are treated and control units available
  data_cohort[,available_treated := sum(treated), list(Cohort,EventTime)]
  data_cohort[,available_control := sum(1-treated), list(Cohort,EventTime)]
  data_cohort = data_cohort[available_treated > 1 & available_control > 1]
  data_cohort[,available_treated := NULL]
  data_cohort[,available_control := NULL]

  # set up variable names
  time_name = varnames$time_name
  outcome_name = paste0(varnames$outcome_name,"_diff")
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name
  cluster_names = varnames$cluster_names
  covariate_names = NULL
  if(!is.null(varnames$covariate_names)){
    covariate_names = paste0(varnames$covariate_names,"_diff")
  }
  fixedeffect_names = varnames$fixedeffect_names

  # check if fixest is available
  check_fixest <- requireNamespace("fixest")
  if(length(check_fixest) == 0){
    # for some reason, this magically fixes the problems with devtools::check when it can't find requireNamespace
  }
  if(!is.null(fixedeffect_names)){
    if(!check_fixest){
      stop("Since varnames$fixedeffect_names is non-missing, you must install the 'fixest' package, which estimates fixed-effects.")
    }
  }

  # get the SE
  all_events = data_cohort[,sort(unique(EventTime))]
  ATTe_SEs = data.table()
  for(ee in all_events){
    if(ee==base_event){
      ATTe_SEs = rbindlist(list(ATTe_SEs, data.table(EventTime=ee, ATTe=0.0, ATTe_SE=0.0)))
      next
    }
    intercepts = NULL
    treateds = NULL
    treated_weights = NULL
    covariates = NULL
    data_event = data_cohort[EventTime==ee]
    all_cohorts = data_event[,unique(Cohort)]
    for(cc in all_cohorts){
      cc_name = copy(cc)
      if(cc<0){
        cc_name = paste0("neg",abs(cc))
      }
      if(check_fixest){
        intercepts = c(intercepts, sprintf("intercept_%s",cc_name))
        data_event[, (sprintf("intercept_%s",cc_name)) := as.numeric((Cohort==cc))]
      }
      treateds = c(treateds, sprintf("treated_%s",cc_name))
      data_event[, (sprintf("treated_%s",cc_name)) := as.numeric((treated==1)*(Cohort==cc))]
      treated_weights = c(treated_weights, data_event[, sum(get(sprintf("treated_%s",cc_name))) ] )
      if(!is.null(covariate_names)){
        for(vv in covariate_names){
          covariates = c(covariates, sprintf("%s_%s",vv,cc_name))
          data_event[, (sprintf("%s_%s",vv,cc_name)) := as.numeric((Cohort==cc)*get(vv))]
        }
      }
    }

    # regression
    OLSlm = NULL
    OLSformula = paste0(outcome_name," ~ -1 + ", paste0(treateds, collapse=" + "))
    if(!is.null(covariate_names)){
      OLSformula = paste0(OLSformula, " + ", paste0(covariates, collapse=" + "))
    }
    if(is.null(fixedeffect_names)){
      if(check_fixest){ # prefer feols() if installed
        OLSformula = paste0(OLSformula, " | Cohort")
        OLSlm = fixest::feols(as.formula(OLSformula),data=data_event)
      }
      if(!check_fixest){ # use lm() if feols() not installed
        print("warning: fixest not found")
        OLSformula = paste0(OLSformula," + ",paste0(intercepts, collapse=" + "))
        OLSlm = lm(as.formula(OLSformula),data=data_event)
      }
    }
    if(!is.null(fixedeffect_names)){ # the case with fixed-effects
      OLSformula = paste0(OLSformula, " | Cohort + ", paste0(fixedeffect_names, collapse=" + "))
      OLSlm = fixest::feols(as.formula(OLSformula),data=data_event)
    }
    OLSmeans = as.numeric(OLSlm$coefficients[treateds])

    # covariance matrix for errors
    OLSvcov = NULL
    if(!is.null(cluster_names)){
      data_event <<- copy(data_event) # due to a well-known scoping bug in R's base lm.predict that no one will fix despite years of requests, this redundancy is necessary!
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

    # finish
    OLSvcov = OLSvcov[treateds, treateds]
    OLSvcov = as.matrix(OLSvcov)
    treated_weights = treated_weights/sum(treated_weights)
    ATTe = as.numeric(t(treated_weights) %*% OLSmeans)
    ATTe_SE = sqrt(as.numeric((t(treated_weights) %*% OLSvcov) %*% treated_weights))
    ATTe_SEs = rbindlist(list(ATTe_SEs, data.table(EventTime=ee, ATTe=ATTe, ATTe_SE=ATTe_SE)))
  }

  return(ATTe_SEs)
}


DiD_getSEs_multipleEventTimes <- function(data_cohort,varnames,Eset,min_event,max_event){

  if(max(Eset) > max_event){
    stop(sprintf("The largest value in Eset is %s. The max_event is %s. Eset should be within the range of events [min_event, max_event].",max(Eset),max_event))
  }
  if(min(Eset) < min_event){
    stop(sprintf("The smallest value in Eset is %s. The max_event is %s. Eset should be within the range of events [min_event, max_event].",min(Eset),min_event))
  }


  # one should always cluster on id_name in a stacked regression in which the same unit can appear multiple times
  varnames$cluster_names = unique(c(varnames$cluster_names, varnames$id_name))

  # check that there are treated and control units available
  data_cohort[,available_treated := sum(treated), list(Cohort,EventTime)]
  data_cohort[,available_control := sum(1-treated), list(Cohort,EventTime)]
  data_cohort = data_cohort[available_treated > 1 & available_control > 1]
  data_cohort[,available_treated := NULL]
  data_cohort[,available_control := NULL]

  # set up variable names
  time_name = varnames$time_name
  outcome_name = paste0(varnames$outcome_name,"_diff")
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name
  cluster_names = varnames$cluster_names
  covariate_names = NULL
  if(!is.null(varnames$covariate_names)){
    covariate_names = paste0(varnames$covariate_names,"_diff")
  }
  fixedeffect_names = varnames$fixedeffect_names

  # check if fixest is available
  check_fixest <- requireNamespace("fixest")
  if(length(check_fixest) == 0){
    # for some reason, this magically fixes the problems with devtools::check when it can't find requireNamespace
  }
  if(!is.null(fixedeffect_names)){
    if(!check_fixest){
      stop("Since varnames$fixedeffect_names is non-missing, you must install the 'fixest' package, which estimates fixed-effects.")
    }
  }

  # get the SE
  data_event = data_cohort[EventTime %in% Eset]
  cohortevents = data.table()
  for(ee in Eset){
    cohortevents = rbindlist(list(cohortevents, data.table(cohort=data_event[EventTime==ee,sort(unique(Cohort))],event=ee)))
  }
  cohortevents = cohortevents[order(cohort,event)]
  numce = nrow(cohortevents)
  cohortevents_covmat = matrix(0, numce, numce)
  cohortevents_weights = rep(NA,numce/2)
  cohortevents_means = rep(NA,numce)

  intercepts = c()
  treateds = c()
  covariates = c()
  treated_weights = c()
  for(rowiter in 1:numce){
    cc = cohortevents[rowiter][,cohort]
    ee = cohortevents[rowiter][,event]
    cc_name = copy(cc)
    if(cc<0){
      cc_name = paste0("neg",abs(cc))
    }
    ee_name = copy(ee)
    if(ee<0){
      ee_name = paste0("neg",abs(ee))
    }
    if(check_fixest){
      intercepts = c(intercepts, sprintf("intercept_%s_%s",cc_name,ee_name))
      data_event[, (sprintf("intercept_%s_%s",cc_name,ee_name)) := (Cohort==cc)*(EventTime==ee)]
    }
    treateds = c(treateds, sprintf("treated_%s_%s",cc_name,ee_name))
    data_event[, (sprintf("treated_%s_%s",cc_name,ee_name)) := (treated==1)*(Cohort==cc)*(EventTime==ee)]
    treated_weights = c(treated_weights, data_event[, sum(get(sprintf("treated_%s_%s",cc_name,ee_name))) ] )
    if(!is.null(covariate_names)){
      for(vv in covariate_names){
        covariates = c(covariates, sprintf("%s_%s_%s",vv,cc_name,ee_name))
        data_event[, (sprintf("%s_%s_%s",vv,cc_name,ee_name)) := (Cohort==cc)*(EventTime==ee)*get(vv)]
      }
    }
  }

  # regression
  OLSlm = NULL
  OLSformula = paste0(outcome_name," ~ -1 + ", paste0(treateds, collapse=" + "))
  if(!is.null(covariate_names)){
    OLSformula = paste0(OLSformula, " + ", paste0(covariates, collapse=" + "))
  }

  if(is.null(fixedeffect_names)){
    if(check_fixest){ # prefer feols() if installed
      OLSformula = paste0(OLSformula, " | Cohort")
      OLSlm = fixest::feols(as.formula(OLSformula),data=data_event)
    }
    if(!check_fixest){ # use lm() if feols() not installed
      OLSformula = paste0(OLSformula," + ",paste0(intercepts, collapse=" + "))
      OLSlm = lm(as.formula(OLSformula),data=data_event)
    }
  }
  if(!is.null(fixedeffect_names)){
    OLSformula = paste0(OLSformula, " | Cohort + ", paste0(fixedeffect_names, collapse=" + "))
    OLSlm = fixest::feols(as.formula(OLSformula),data=data_event)
  }
  OLSmeans = as.numeric(OLSlm$coefficients[treateds])

  # covariance matrix for errors
  OLSvcov = NULL
  if(!is.null(cluster_names)){
    data_event <<- copy(data_event) # due to a well-known scoping bug in R's base lm.predict that no one will fix despite years of requests, this redundancy is necessary!
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

  # standard errors
  OLSvcov = OLSvcov[treateds, treateds]
  treated_weights = treated_weights/sum(treated_weights)
  ATT_Eset = as.numeric(t(treated_weights) %*% OLSmeans)
  ATT_Eset_SE = sqrt(as.numeric((t(treated_weights) %*% OLSvcov) %*% treated_weights))
  return(data.table(Eset = paste0(Eset,collapse=","), ATT_Eset=ATT_Eset, ATT_Eset_SE=ATT_Eset_SE))
}

