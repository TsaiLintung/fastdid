#' Create Event Data
#'
#' Generates event data for treated households based on the input parameters.
#'
#' @param maindata A data.table containing the main data.
#' @param timevar The name of the column representing time.
#' @param unitvar The name of the column representing unit (e.g., household).
#' @param cohortvar The name of the column representing the period when the event of interest arises.
#' @param control_group Options: "both", "never", and "later".
#' @param base_time The base time for analysis (default is -1).
#' @param anycohortvar The name of the column representing the period when ANY kind of event arises, used to find control cohorts (default is NULL).
#' @param covariate_base_stratify Vector of variable names; treatment effects are stratified by their joint unique values (default is NULL).
#' @param covariate_base_balance Vector of variable names to include in finding exact matches for treated units (default is NULL).
#' @param covariate_base_balance_linear Vector of variable names to include linearly in propensity score estimator (default is NULL).
#' @param covariate_base_balance_linear_subset Vector of variable names within which the linear propensity score regression will be estimated (default is NULL).
#' @param covariate_base_support Vector of variable names on which to impose common support across treated and control units (default is NULL).
#' @param balanced_panel Logical, whether to enforce a balanced panel (default is TRUE).
#' @param check_not_treated If TRUE, only allows controls to be used for a cohort if they are observed in the dataset to be untreated in that cohort (default is FALSE).
#' @param stratify_by_cohort Logical, whether to stratify by cohort (default is FALSE).
#' @param lower_event_time Earliest period (relative to treatment time) for which to estimate effects (default is -Inf).
#' @param upper_event_time Latest period (relative to treatment time) for which to estimate effects (default is Inf).
#' @param min_control_gap Controls must be treated at least this many years away from matched treated units (default is NULL).
#' @param max_control_gap Controls must be treated no more than this many years away from matched treated units (default is NULL).
#' @param treat_criteria Replacement of onsetagevar, provide a character column that must == TRUE for the treated (default is NULL).
#' @param base_restrict A single variable, restricting to households for which the variable == 1 in the base period (default is NULL).
#' @param base_restrict_treated A single variable, restricting TREATED GUYS to households for which the variable == 1 in the base period (default is NULL).
#' @param verbose Logical, whether to display messages (default is TRUE).
#' @param copy_dataset Logical, whether to create a copy of the dataset to avoid changing the original input (default is TRUE).
#' @param validate Logical, whether to validate the data (default is TRUE).
#' @param combine Logical, whether to combine the resulting event data or return a list (default is TRUE).
#'
#' @return Event data as a data.table.
#'
#' @details This function generates event data for treated households based on the provided input parameters.
#' 
#' @examples
#' 
#' create_event_data(maindata, "time", "id", "cohort", "both")
#' 
#'
#' @export
create_event_data<-function(maindata,
                            timevar,
                            unitvar,
                            cohortvar, #period when the particular event of interest arises
                            control_group, #options: "both", "never", and "later"
                            base_time = -1,
                            anycohortvar = NULL, #period when ANY kind of event arises, used to find control cohorts
                            
                            
                            #covaraites
                            covariate_base_stratify=NULL, #vector of variable names, treatment effects are stratified by their joint unique values.
                            covariate_base_balance=NULL, #vector of variable names to include in finding exact matches for treated units
                            covariate_base_balance_linear=NULL, #vector of variable names to include linearly in propensity score estimator
                            covariate_base_balance_linear_subset=NULL, #vector of variable names within which the linear propensity score regression will be estimated
                            covariate_base_support=NULL, #vector of variable names on which to impose common support across treated and control units
                            
                            #checks and control
                            balanced_panel = TRUE,
                            check_not_treated = FALSE, #If TRUE, only allows controls to be used for cohort o if they are observed in the dataset to be untreated in year o.                      
                            stratify_by_cohort= FALSE,
                            lower_event_time = -Inf, #Earliest period (relative to treatment time) for which to estimate effects
                            upper_event_time = Inf, #Latest period (relative to treatment time) for which to estimate effects
                            min_control_gap = NULL, #controls must be treated at least this many years away from matched treated units
                            max_control_gap = NULL, #controls must be treated no more than this many years away from matched treated units
                            treat_criteria = NULL, #replacement of onsetagevar, provide a character col that must == TRUE for the treated
                            base_restrict = NULL, #a single var, restricting to households for which var == 1 in base period,
                            base_restrict_treated = NULL, #a single var, restricting TREATED GUYS to households for which var == 1 in base period
                            
                            #function behavior
                            verbose = TRUE,
                            copy_dataset = TRUE,
                            validate = TRUE,
                            combine = TRUE

) 
{

  # argument validation -----------------------------------
  
  
  if(!is.data.table(maindata)) stop("rawdata must be a data.table")
  
  dt_names <- names(maindata)
  name_message <- "__ARG__ must be a character scalar and a name of a column from the dataset."
  check_arg(timevar, unitvar, cohortvar, "scalar charin", .choices = dt_names, .message = name_message)

  covariate_message <- "__ARG__ must be NULL or a character vector which are all names of columns from the dataset."
  check_arg(covariate_base_stratify, covariate_base_balance, covariate_base_balance, covariate_base_balance_linear, 
            covariate_base_balance_linear_subset, covariate_base_support,
            "NULL | multi charin", .choices = dt_names, .message = covariate_message)
  
  checkvar_message <- "__ARG__ must be NULL or a character scalar which are all names of columns from the dataset."
  check_arg(anycohortvar, base_restrict, base_restrict_treated, treat_criteria,
            "NULL | scalar charin", .choices = dt_names, .message = checkvar_message)

  check_arg(control_group, "scalar charin", .choices = c("both", "never", "later"))
  check_arg(balanced_panel, check_not_treated, stratify_by_cohort, "scalar logical")

  if(lower_event_time > base_time) stop("lower_event_time must lie below base_time")
  if("anycohort" %in% dt_names) warning("column name 'anycohort' detected in input data. The input may already be changed by past function call, avoid this by calling with copy == TRUE.")

  # name change ----------------------------------

  if(copy_dataset){
    maindata <- copy(maindata)
    if(verbose) message("copying the dataset to avoid changing the original input.")
  }
  
  if(is.null(anycohortvar)){
    maindata[,anycohortvar:=get(cohortvar)]
    anycohortvar<-"anycohortvar"
  }
  setnames(maindata,c(timevar,unitvar,cohortvar,anycohortvar),c("time","id","cohort","anycohort"))
  setkey(maindata, id)
  if(!is.null(base_restrict)) setnames(maindata,base_restrict,"base_restrict")
  if(!is.null(base_restrict_treated)) setnames(maindata,base_restrict_treated,"base_restrict_treated")
  if(is.null(covariate_base_balance_linear_subset)) covariate_base_balance_linear_subset <- covariate_base_balance
  
  # data validation
  if(validate){
    covariates <-  c(covariate_base_stratify, covariate_base_balance, covariate_base_support, covariate_base_balance_linear)
    maindata <- maindata |> validate_eventdata(covariates, balanced_panel)
  }
  
  # main part --------------------------------------

  maindata <- maindata |> covariate_to_factor(base_time, 
                                              covariate_base_stratify, covariate_base_balance, covariate_base_support, covariate_base_balance_linear_subset,
                                              stratify_by_cohort)
  
  
  stacked_list <- maindata |> stack_for_cohort(control_group, base_time,
                                               treat_criteria, lower_event_time, upper_event_time, 
                                               min_control_gap, max_control_gap, check_not_treated) 

  factor_cols <- c("id", "cohort", "time_pair", "time")

  event_list <- stacked_list |> lapply(function (x) check_stacked_data(x, base_time,
                                                                       balanced_panel, base_restrict, base_restrict_treated, covariate_base_balance_linear) |> 
                                                    stack_for_event_time(base_time) |>  
                                                    var_to_factor(factor_cols) |> 
                                                    estimate_ipw(covariate_base_stratify, covariate_base_balance, covariate_base_balance_linear))
  
  
  if(combine){
    return(rbindlist(event_list))
  } else {
    return(event_list)
  }
}

# Functions for important steps --------------------------------------------------------------

validate_eventdata <- function(maindata, covariates, balanced_panel){
  
  if(any(is.na(maindata$cohort))) stop("cohort variable should not be missing (it can be infinite instead)")
  if(any(is.na(maindata$anycohort))) stop("anycohort variable should not be missing (it can be infinite instead)")
  
  if(balanced_panel){
    
    #check if any is dup
    if(anyDuplicated(maindata[, .(id, time)])){
      dup <- duplicated(maindata[,.(id, time)])
      warning(nrow(dup_id), " units is observed more than once in some periods, enforcing balanced panel by dropping them")
      maindata <- maindata[!id %in% dup]
    }
    
    #check if any is missing
    id_count <- maindata[, .(count = .N), by = id]
    time_period <- maindata[, uniqueN(time)]
    if(any(id_count[, count < time_period])){
      mis_id <- id_count[count < time_period]
      warning(nrow(mis_id), " units is missing in some periods, enforcing balanced panel by dropping them")
      maindata <- maindata[!id %in% mis_id$id]
    }
    
  }

  #check covariate is time-invariant
  for(out in covariates){
    if(is.null(out))next
    
    cov_count <- maindata[, .(count = uniqueN(c(out))), by = "id"]
    if(any(cov_count[, count] > 1)){stop("covariates need to be time invariant")}

  }
  
  return(maindata)
  
}

covariate_to_factor <- function(dt, base_time,
                                covariate_base_stratify, covariate_base_balance, covariate_base_support, covariate_base_balance_linear_subset,
                                stratify_by_cohort){
  
  #turn covariates interaction to factors
  if(!is.null(covariate_base_stratify)){
    stratifyvars <- ifelse(stratify_by_cohort, c(covariate_base_stratify, "cohort"), covariate_base_stratify)
  } else {stratifyvars <- covariate_base_stratify}
  
  for(covariate_type in c("stratify", "balancevars", "balancevars_linear_subset", "supportvars")){
    
    cov_vars <- switch(covariate_type, 
                       stratify = stratifyvars,
                       balancevars = covariate_base_balance,
                       balancevars_linear_subset = covariate_base_balance_linear_subset,
                       supportvars = covariate_base_support)
    if(!is.null(cov_vars)){
      dt[,(covariate_type) :=  do.call(finteraction, dt[, cov_vars, with = FALSE])]
    } else {
      dt[,(covariate_type) := factor(1,levels=c(1,"OMIT"))]
    }
  }
  return(dt)
  
}

stack_for_cohort <- function(maindata, 
                             control_group, base_time,
                             treat_criteria, lower_event_time, upper_event_time, 
                             min_control_gap, max_control_gap, check_not_treated){
  
  setkey(maindata, id)
  
  #stacking control cohorts.
  dt_list <- list()
  cohorts <- maindata[, unique(cohort)]
  cohorts <- cohorts[cohorts != maindata[, min(time)] & !is.infinite(cohorts)] #the min cohort won't have valid base_time control
  for(o in cohorts){
    
    treatcohort<-maindata[cohort == o,] #cohort treated on the first period is not useful
    
    if(!is.null(treat_criteria)){
      treatcohort <- treatcohort[get(treat_criteria == TRUE)]
    }
    
    treatcohort[,treated:=1]
    treatcohort[,event_time:=time-cohort]
    treatcohort<-treatcohort[event_time >= lower_event_time & event_time <= upper_event_time,]
    
    
    #find the relevant control group
    if(control_group=="both") {
      controlcohort <- maindata[(anycohort > o | is.infinite(anycohort)) ,]
    }
    if(control_group=="later"){
      controlcohort <- maindata[(anycohort > o & !is.infinite(anycohort)) ,]
    }
    if(control_group=="never") {
      controlcohort <- maindata[is.infinite(anycohort) ,]
    }
    
    #enfornce min/max control gap
    if (!is.null(max_control_gap)) controlcohort <- controlcohort[anycohort - o <=  max_control_gap,]
    if (!is.null(min_control_gap))  controlcohort <- controlcohort[anycohort - o >=  min_control_gap,]
    
    controlcohort[,cohort := o]
    controlcohort[, event_time := time-cohort]
    
    #a control would only be useful if it is observed in the base period
    controlcohort[, useful_control := any(event_time == base_time), by = "id"]
    controlcohort <- controlcohort[useful_control == TRUE]
    controlcohort[, useful_control := NULL]
    
    controlcohort[, treated := 0]
    
    controlcohort <- controlcohort[anycohort - cohort > event_time]
    
    if(!is.infinite(-lower_event_time)|!is.infinite(upper_event_time)){
      controlcohort <- controlcohort[event_time >= lower_event_time & event_time <= upper_event_time]
    }
    
    if(check_not_treated){
      controlcohort[ ,obscohort := any(time == cohort),by=id]
      controlcohort<-controlcohort[obscohort==1,]
      controlcohort[,obscohort:=NULL]
    }
    
    controlcohort[, anycohort := NULL]
    treatcohort[, anycohort := NULL]
    
    dt_list <- c(dt_list, list(list(treat = treatcohort, control = controlcohort)))
    
  }
  
  return(dt_list)
  
} 

check_stacked_data<- function(dt_list, base_time,
                              balanced_panel, base_restrict, base_restrict_treated, covariate_base_balance_linear){
  
  treatdata <- dt_list$treat
  controldata <- dt_list$control
  
  if(!balanced_panel){
    
    #only needed when panel is not already balanced
    treatdata[,obsbase:=sum(event_time==base_time),by=.(id,cohort)]
    controldata[,obsbase:=sum(event_time==base_time),by=.(id,cohort)]
    if(max(treatdata$obsbase)>1) stop("Error: some treated units are observed more than once in the reference period")
    if(max(controldata$obsbase)>1) stop("Error: some control units are observed more than once in the reference period")
    treatdata <- treatdata[obsbase==1,]
    controldata <- controldata[obsbase==1,]
    treatdata[,obsbase := NULL]
    controldata[,obsbase := NULL]
    
  }
  
  #check base-restrict
  if(!is.null(base_restrict)){
    controldata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
    controldata <- controldata[base_restrict == 1,]
    treatdata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
    treatdata <- treatdata[base_restrict == 1,]
  }
  if(!is.null(base_restrict_treated)){
    treatdata[,base_restrict_treated := max(base_restrict_treated * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
    treatdata <- treatdata[base_restrict_treated == 1,]
  }
  
  #check common support
  if(!is.null(covariate_base_balance_linear)){
    
    #common support checking is only needed when there is linear regression - interpolation and extrapolation
    #o.w. a propensity score within 0,1 means it has common support
    
    treatdata[,temp:=finteraction(balancevars,supportvars,stratify,event_time,cohort)]
    controldata[,temp:=finteraction(balancevars,supportvars,stratify,event_time,cohort)]
    commonvals<-intersect(unique(treatdata$temp),unique(controldata$temp))
    treatdata<-treatdata[temp%in%commonvals,]
    controldata<-controldata[temp%in%commonvals,]
    
    treatdata[,temp:=NULL]
    controldata[,temp:=NULL]
    rm(commonvals)
    gc()
    
  }
  
  dt_list$treat <- treatdata
  dt_list$control <- controldata  
  
  return(dt_list)
  
}

stack_for_event_time <- function(eventdata, base_time){
  
  controldata <- eventdata$control
  treatdata <- eventdata$treat
 
  
  controldata[,`:=`(max_event_time = max(event_time),
                  min_event_time = min(event_time)) , by = "id"]
  treatdata[,`:=`(max_event_time = max(event_time),
                    min_event_time = min(event_time)) , by = "id"]
   
  event_times<-treatdata[treated == 1,funique(event_time)] #equals to 0, 
  
  setkey(treatdata, event_time)
  setkey(controldata, event_time)
  
  double_stack_list <- list()
  for(t in event_times){
    
    if(t == -1){next}
    
    event_time_treat <- treatdata[t <= max_event_time & t >= min_event_time & (event_time == -1 | event_time == t)] #already made sure every unit have -1 in it
    event_time_treat[,time_pair := t]
    
    event_time_control <- controldata[t <= max_event_time & t >= min_event_time & (event_time == -1 | event_time == t)] #already made sure every unit have -1 in it
    event_time_control[,time_pair := t]
    
    #remove uncessary cols ASAP
    event_time_treat[, max_event_time := NULL]
    event_time_treat[, min_event_time := NULL]
    event_time_control[, max_event_time := NULL]
    event_time_control[, min_event_time := NULL]
    
    double_stack_list<-c(double_stack_list, list(event_time_treat), list(event_time_control))
    
  }
  
  double_stack_dt <- rbindlist(double_stack_list, use.names=TRUE)
  
  
  
  return(double_stack_dt)
  
}

var_to_factor <- function(x, factor_cols){
    #turn the cols into factors
    for(col in factor_cols){
      x[, (col) := qF(get(col))]
    }
    #keep a numeric version for later comparisons
    x[, event_time_fact := qF(event_time)]
}


estimate_ipw<- function(eventdata, covariate_base_stratify, covariate_base_balance, covariate_base_balance_linear){
  
  if(is.null(covariate_base_stratify)){
    eventdata[, stratify := 1]
  }
  
  if(is.null(covariate_base_balance)){
    eventdata[, balancevars := 1]
  }
  
  #construct the call  
  if(is.null(covariate_base_balance_linear)){
    
    call <- "treated ~ 1 | finteraction(cohort,time_pair,event_time_fact,stratify,balancevars)" 
    
  } else {
    call <- paste0("treated ~",paste0(covariate_base_balance_linear,collapse="+",
                                      sep="finteraction(cohort,time_pair,event_time_fact,stratify,balancevars_linear_subset,)"),
                   "| finteraction(cohort,event_time_fact,stratify,balancevars)")
  }
  
  #estimate propensity score
  #in a specific cohort-time_pair estimation, propensity to get treated at event_time, given stratify and balance
  eventdata[,pval:= feols(as.formula(call), data = eventdata, lean = FALSE)$fitted.values]
  
  #only keep propensity score between 0,1 is equivalent to checking common support
  eventdata <- eventdata[pval < 1 & pval > 0]
  eventdata[,pweight := ifelse(treated == 1, 1, pval/(1-pval))]
  eventdata[,pval:=NULL]
  
  eventdata[, treated := qF(treated)]#can't do it before feols call
  
  return(eventdata)
  
}





