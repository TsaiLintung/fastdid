
# create event data -----------------------------------

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
                            combine = FALSE
                            
) 
{
  
  # argument validation -----------------------------------
  
  
  if(!is.data.table(maindata)) stop("rawdata must be a data.table")
  
  dt_names <- names(maindata)
  name_message <- "__ARG__ must be a character scalar and a name of a column from the dataset."
  check_arg(timevar, unitvar, cohortvar, "scalar charin", .choices = dt_names, .message = name_message)
  
  covariate_message <- "__ARG__ must be NULL or a character vector which are all names of columns from the dataset."
  check_arg(covariate_base_stratify, covariate_base_balance, covariate_base_balance, covariate_base_balance_linear, 
            covariate_base_balance_linear_subset, 
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
    # if(verbose) message("copying the dataset to avoid changing the original input.")
  }
  
  if(is.null(anycohortvar)){
    maindata[,anycohortvar:=get(cohortvar)]
    anycohortvar<-"anycohortvar"
  }
  setnames(maindata,c(timevar,unitvar,cohortvar,anycohortvar),c("time","id","cohort","anycohort"))
  setkey(maindata, id)
  if(!is.null(base_restrict)) setnames(maindata,base_restrict,"base_restrict")
  if(!is.null(base_restrict_treated)) setnames(maindata,base_restrict_treated,"base_restrict_treated")
  if(!is.null(treat_criteria)) setnames(maindata,treat_criteria,"treat_criteria")
  if(is.null(covariate_base_balance_linear_subset)) covariate_base_balance_linear_subset <- covariate_base_balance
  
  # data validation
  if(validate){
    covariates <-  c(covariate_base_stratify, covariate_base_balance, covariate_base_balance_linear)
    maindata <- maindata |> validate_eventdata(covariates, balanced_panel)
  }
  
  # main part --------------------------------------
  
  maindata <- maindata |> covariate_to_factor(base_time, 
                                              covariate_base_stratify, covariate_base_balance, covariate_base_balance_linear_subset,
                                              stratify_by_cohort)
  
  
  stacked_list <- maindata |> stack_for_cohort(control_group, base_time,
                                               treat_criteria, lower_event_time, upper_event_time, 
                                               min_control_gap, max_control_gap, check_not_treated) 
  
  rm(maindata)
  
  event_list <- stacked_list |> lapply(function (x) check_stacked_data(x, base_time,
                                                                       balanced_panel, base_restrict, base_restrict_treated, covariate_base_balance_linear) |> 
                                         stack_for_event_time(base_time))
  
  rm(stacked_list)
  
  event_list <- unlist(event_list, recursive = FALSE)
  
  factor_cols <- c("cohort", "time_pair")
  event_list <- event_list |> lapply(function (x) var_to_factor(x, factor_cols) |> 
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
                                covariate_base_stratify, covariate_base_balance, covariate_base_balance_linear_subset,
                                stratify_by_cohort){
  
  #turn covariates interaction to factors
  if(!is.null(covariate_base_stratify)){
    stratifyvars <- ifelse(stratify_by_cohort, c(covariate_base_stratify, "cohort"), covariate_base_stratify)
  } else {stratifyvars <- covariate_base_stratify}
  
  for(covariate_type in c("stratify", "balancevars", "balancevars_linear_subset")){
    
    cov_vars <- switch(covariate_type, 
                       stratify = stratifyvars,
                       balancevars = covariate_base_balance,
                       balancevars_linear_subset = covariate_base_balance_linear_subset)
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
  cohorts <- maindata[, unique(anycohort)]
  cohorts <- cohorts[cohorts != maindata[, min(time)] & !is.infinite(cohorts)] #the min cohort won't have valid base_time control
  
  if(control_group == "later" | maindata[is.infinite(anycohort), .N] == 0){
    cohorts <- cohorts[cohorts != max(cohorts)] #if no never treated, the last cohort won't have never treated
  }
  
  potential_controldata <- maindata[anycohort > time]
  
  for(o in cohorts){
    
    treatcohort<-maindata[anycohort == o,] #cohort treated on the first period is not useful
    
    if(!is.null(treat_criteria)){
      treatcohort <- treatcohort[treat_criteria == TRUE]
    }
    
    treatcohort[,treated:=1]
    treatcohort[,event_time:=time-cohort]
    treatcohort<-treatcohort[event_time >= lower_event_time & event_time <= upper_event_time,]
    
    
    #find the relevant control group
    if(control_group=="both") {
      controlcohort <- potential_controldata[(anycohort > o | is.infinite(anycohort)) ,]
    }
    if(control_group=="later"){
      controlcohort <- potential_controldata[(anycohort > o & !is.infinite(anycohort)) ,]
    }
    if(control_group=="never") {
      controlcohort <- potential_controldata[is.infinite(anycohort) ,]
    }
    
    #enfornce min/max control gap
    if (!is.null(max_control_gap)) controlcohort <- controlcohort[anycohort - o <=  max_control_gap,]
    if (!is.null(min_control_gap))  controlcohort <- controlcohort[anycohort - o >=  min_control_gap,]
    
    controlcohort[,cohort := o]
    controlcohort[,event_time := time-cohort]
    
    #a control would only be useful if it is observed in the base period
    controlcohort[, useful_control := any(event_time == base_time), by = "id"]
    controlcohort <- controlcohort[useful_control == TRUE]
    controlcohort[, useful_control := NULL]
    
    controlcohort[, treated := 0]
    
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
    
    if(nrow(controlcohort)==0|nrow(treatcohort)==0){next}
    
    cohortdata <- rbind(treatcohort, controlcohort)
    
    dt_list <- c(dt_list, list(cohortdata))
    
  }
  
  return(dt_list)
  
} 

check_stacked_data<- function(eventdata, base_time,
                              balanced_panel, base_restrict, base_restrict_treated, covariate_base_balance_linear){
  
  if(!balanced_panel){
    
    #only needed when panel is not already balanced
    eventdata[,obsbase:=sum(event_time==base_time),by=.(id,cohort)]
    if(max(eventdata$obsbase)>1) warning("Error: some units are observed more than once in the reference period")
    eventdata <- eventdata[obsbase==1,]
    eventdata[,obsbase := NULL]
    
  }
  
  #check base-restrict
  if(!is.null(base_restrict)){
    eventdata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
    eventdata <- eventdata[base_restrict == 1,]
  }
  if(!is.null(base_restrict_treated)){
    eventdata[,base_restrict_treated := max(base_restrict_treated * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
    eventdata <- eventdata[base_restrict_treated == 1 | treated == 0,] #only limit it for treated units
  }
  
  #check common support
  if(!is.null(covariate_base_balance_linear)){
    
    #common support checking is only needed when there is linear regression - interpolation and extrapolation
    #o.w. a propensity score within 0,1 means it has common support
    
    eventdata[,temp:=finteraction(balancevars,stratify,event_time,cohort)]
    commonvals<-intersect(eventdata[treated == 1, unique(temp)],eventdata[treated == 0, unique(temp)])
    eventdata<-eventdata[temp%in%commonvals,]
    eventdata[,temp:=NULL]
    rm(commonvals)
    gc()
    
  }
  
  return(eventdata)
  
}

stack_for_event_time <- function(eventdata, base_time){
  
  eventdata[,`:=`(max_event_time = max(event_time),
                  min_event_time = min(event_time)) , by = "id"]
  
  event_times<-eventdata[treated == 1,funique(event_time)] #equals to 0, 
  
  double_stack_list <- list()
  for(t in event_times){
    
    if(t == -1){next}
    
    event_time_data <- eventdata[t <= max_event_time & t >= min_event_time & (event_time == -1 | event_time == t)] #already made sure every unit have -1 in it
    
    event_time_data[, time_pair := t]
    
    #remove uncessary cols ASAP
    event_time_data[, max_event_time := NULL]
    event_time_data[, min_event_time := NULL]
    
    if(event_time_data[treated == 1, .N]==0|event_time_data[treated == 0, .N]==0){next}
    
    new_list <- list(event_time_data)
    names(new_list) <- as.character(t)
    double_stack_list<-c(double_stack_list, new_list)
    
  }
  
  cohort <- eventdata[, as.character(first(cohort))]
  names(double_stack_list) <- str_c(names(double_stack_list), paste0(".", cohort))
  
  return(double_stack_list)
  
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
  
  
  fe_var <- "event_time_fact"
  
  if(!is.null(covariate_base_stratify)){
    fe_var <- paste0(fe_var, ",stratify")
  }
  
  if(!is.null(covariate_base_balance)){
    fe_var_wb <- paste0(fe_var, ",balancevars")
  } else {fe_var_wb <- fe_var}
  
  #construct the call  
  if(is.null(covariate_base_balance_linear)){
    
    call <- paste0("treated ~ 1 | finteraction(", fe_var_wb,")") 
    
  } else {
    call <- paste0("treated ~ 1 | ",paste0(paste0(paste0("finteraction(", fe_var, ",balancevars_linear_subset)["),
                                                  covariate_base_balance_linear, "]"),collapse="+"),
                   paste0("+ finteraction(",fe_var_wb,")"))
  }
  
  #estimate propensity score
  #in a specific cohort-time_pair estimation, propensity to get treated at event_time, given stratify and balance
  eventdata[,pval:= feols(as.formula(call), data = eventdata, lean = FALSE, note = FALSE)$fitted.values]
  
  #only keep propensity score between 0,1 is equivalent to checking common support
  eventdata <- eventdata[pval < 1 & pval > 0]
  eventdata[,pweight := ifelse(treated == 1, 1, pval/(1-pval))]
  eventdata[,pval:=NULL]
  
  eventdata[, treated := qF(treated)]#can't do it before feols call
  
  return(eventdata)
  
}

# get event result ------------------------------------------------


get_event_result <- function(eventdata,
                             variable,
                             result_type,
                             covariate = NULL,
                             clustervar = "id",
                             weights = "pweight",
                             base_time = -1,
                             trends = FALSE,
                             mem.clean = FALSE,
                             cohort_obs = NULL) {
  
  if(!is.data.table(eventdata[[1]])){stop("please provide a data.table")}
  
  if(result_type %in% c("cohort_event_time", "dynamic")){
    all_results <- data.table()
    for(i in seq(1, length(eventdata))){
      
      results <- compute_event_result(eventdata = eventdata[[i]],
                                      variable = variable,
                                      result_type = "cohort_event_time",
                                      covariate = covariate,
                                      clustervar=clustervar, 
                                      weights=weights,
                                      base_time = base_time,
                                      trends=trends,
                                      mem.clean = mem.clean)
      all_results <- rbind(all_results, results)
    }
    
    if(result_type == "cohort_event_time") return(all_results)
    else { #dynamic
      
      all_results <- all_results |> merge(cohort_obs, by = "cohort")
      all_results[, weight := count/sum(count), by = c("event_time", "stratify", "outcome")]
      
      all_results <- all_results[, .(Estimate = sum(Estimate*weight),
                                     `Std. Error` = sqrt(sum(`Std. Error`^2*weight^2)),
                                     obs = sum(obs)) , by = c("event_time", "stratify", "outcome")]
      
      return(all_results)
      
    }
    
  } else {
    
    eventdata <- rbindlist(eventdata)
    
    results <- compute_event_result(eventdata = eventdata,
                                    variable = variable,
                                    result_type = result_type,
                                    covariate = covariate,
                                    clustervar=clustervar, 
                                    weights=weights,
                                    base_time = base_time,
                                    trends=trends,
                                    mem.clean = mem.clean)
    
    return(eventdata)
    
  }
  
}

compute_event_result <- function(eventdata,
                                 variable,
                                 result_type,
                                 covariate = NULL,
                                 clustervar = "id",
                                 weights = "pweight",
                                 base_time = -1,
                                 trends = FALSE,
                                 mem.clean = FALSE) {
  
  #allow list input (useful when not combined)
  if(!is.data.table(eventdata)){stop("input panel must be data.table")} 
  
  #input validation
  dt_names <- names(eventdata)
  check_arg(variable,"multi charin", .choices = dt_names)
  check_arg(result_type,"charin", .choices = c("cohort_event_time", "dynamic", "pooled", "means", "covariate"))
  check_arg(clustervar, weights,"charin", .choices = dt_names)
  if(is.null(covariate) & result_type == "covariates"){stop("please provide covariate when getting covariates result")}
  if(nrow(eventdata)==0){ 
    stop("event panel is empty")
  }
  if(any(eventdata[,lapply(.SD, stats::var), .SDcols = variable] == 0)){
    stop(variable[eventdata[,lapply(.SD, stats::var), .SDcols = variable] == 0], " have no variation")
  }
  
  eventdata <- eventdata[!is.na(get(weights))]
  
  eventdata <- eventdata |> construct_event_variable(result_type, variable, covariate, base_time)
  
  call <- get_estimate_call(result_type, variable, trends, covariate)
  
  results <- eventdata |> estimate_att(call, weights, clustervar, result_type,mem.clean )
  
  return(results)
  
}

# helper functions --------------------------------------------

construct_event_variable <- function(eventdata, result_type, variable, covariate, base_time){
  
  base_event_stratify <- paste0(c(base_time,1),collapse=".")
  
  if(result_type == "cohort_event_time"){
    
    eventdata[,unitfe := finteraction(cohort, time_pair,treated,stratify)] #the level difference between treated and untreated units within each time_pair and stratify
    
    #cohort_event_time_stratify
    eventdata[,cohort_event_time_stratify:= finteraction(event_time_fact,cohort,stratify)]
    eventdata[event_time==base_time, cohort_event_time_stratify := base_event_stratify]
    eventdata[,cohort_event_time_stratify:=relevel(cohort_event_time_stratify,ref = base_event_stratify)]
    
    #treated_cohort_event_time_stratify
    eventdata[,treated_cohort_event_time_stratify := cohort_event_time_stratify]
    eventdata[event_time==base_time | treated == 0 ,treated_cohort_event_time_stratify := base_event_stratify]
    eventdata[,treated_cohort_event_time_stratify:=relevel(treated_cohort_event_time_stratify,ref = base_event_stratify)]
    
  }else if(result_type == "dynamic"){
    
    eventdata[,unitfe := finteraction(time_pair,treated,stratify)] #the level difference between treated and untreated units within each time_pair and stratify
    
    #event_time_stratify
    eventdata[,event_time_stratify:= finteraction(event_time_fact,stratify)]
    eventdata[event_time==base_time, event_time_stratify := base_event_stratify]
    eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = base_event_stratify)]
    
    #treated_event_time_stratify
    eventdata[,treated_event_time_stratify := event_time_stratify]
    eventdata[event_time==base_time | treated == 0 ,treated_event_time_stratify := base_event_stratify]
    eventdata[,treated_event_time_stratify:=relevel(treated_event_time_stratify,ref = base_event_stratify)]
    
  }else if(result_type == "pooled"){
    
    eventdata[,post:=event_time >= 0]
    
    eventdata[,unitfe := finteraction(time_pair,treated,stratify)] #the level difference between treated and untreated units within each time_pair and stratify
    
    #event_time_stratify
    eventdata[,event_time_stratify:= finteraction(event_time_fact,stratify)]
    eventdata[event_time==base_time, event_time_stratify := base_event_stratify]
    eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = base_event_stratify)]
    
    #treated_post_stratify
    eventdata[,treated_post := qF((treated == 1) * (post == 1))]
    eventdata[,treated_post_stratify := finteraction(treated_post,stratify)]
    eventdata[treated_post == 0 ,treated_post_stratify := paste0(c(0,1),collapse=".")]
    eventdata[,treated_post_stratify:=relevel(treated_post_stratify,ref = paste0(c(0,1),collapse="."))]
    
    #treated_pre_stratify
    eventdata[,treated_pre := qF((treated == 1) * (post == 0))]
    eventdata[,treated_pre_stratify := finteraction(treated_pre,stratify)]
    eventdata[treated_pre==0 | event_time==base_time,treated_pre_stratify := 1]
    
  }else if(result_type == "means"){
    
    #stratify
    eventdata[, stratify:= factor(eventdata$stratify, levels = c(levels(eventdata$stratify), "empty"))]
    
    #subset
    eventdata <- eventdata[treated == 1 & event_time == base_time,]
    
  }else if(result_type == "covariate"){
    
    eventdata[,int_event_time_stratify:=interaction(event_time,stratify)]
    estimate_data <- eventdata[as.numeric(as.character(event_time)) == as.numeric(as.character(time_pair)) | as.numeric(as.character(time_pair)) == 0,]
    
  }
}

get_estimate_call <- function(result_type, variable, trends, covariate = NULL){
  
  if(result_type == "covariate"){
    
    call <- paste0("c(", paste0(covariate,collapse=",") , ") ~ int_event_time_stratify + int_event_time_stratify:treated- 1" )
    
  } else if (result_type == "means"){
    
    call <- paste0("c(",paste0(variable,collapse=",") ,") ~ interaction(stratify) - 1")
    
  } else {
    
    outcomes_call <- paste0("c(", paste0(variable,collapse=","), ")")
    
    if(result_type == "cohort_event_time"){
      event_stratify_call <- ifelse(trends, "~ cohort_event_time_stratify + treated_cohort_event_time_stratify |", 
                                    "~ treated_cohort_event_time_stratify | cohort_event_time_stratify +")
    }else if(result_type == "dynamic"){
      event_stratify_call <- ifelse(trends, "~ event_time_stratify + treated_event_time_stratify |", 
                                    "~ treated_event_time_stratify | event_time_stratify +")
    }else if(result_type == "pooled"){
      event_stratify_call <- ifelse(trends, "~ event_time_stratify + treated_pre_stratify + treated_post_stratify |", 
                                    "~ treated_pre_stratify + treated_post_stratify | event_time_stratify +")
    }
    
    call <- paste0(outcomes_call, event_stratify_call, " unitfe" )
    
  }
  
  return(call)
  
}


estimate_att <- function(eventdata, call, weights, clustervar, result_type, mem.clean){
  all_results_dt <- data.table()
  
  if(!eventdata[, uniqueN(stratify)] > 1){
    
    results<-feols(as.formula(call),
                   data = eventdata,
                   weights= eventdata[,get(weights)],
                   cluster=clustervar, lean = TRUE, mem.clean = mem.clean, notes = FALSE)
    
    all_results_dt <- parse_event_result(results, result_type)
    
  } else {
    
    for(stratify_type in eventdata[, unique(stratify)]){
      
      strat_eventdata <- eventdata[stratify == stratify_type]
      
      strat_results<-feols(as.formula(call),
                           data = strat_eventdata,
                           weights= strat_eventdata[,get(weights)],
                           cluster=clustervar, lean = TRUE, mem.clean = mem.clean, notes = FALSE)
      strat_results_dt <- parse_event_result(strat_results, result_type)
      strat_results_dt[, stratify := stratify_type]
      all_results_dt <- rbind(all_results_dt, strat_results_dt)
      
    }
    
  }
  
  return(all_results_dt)
}

get_event_time <- function(x){
  start <- str_locate(x, "stratify")[2]
  end <- str_locate(x, "\\.")[1]
  return(as.numeric(str_sub(x, start + 1, end - 1)))
}

get_cohort <- function(x){
  
  dot_pos <- str_locate_all(x, "\\.")
  start <- dot_pos[[1]][1, 1]
  end <- dot_pos[[1]][2, 1]
  return(as.numeric(str_sub(x, start + 1, end - 1)))
  
}

get_stratify <- function(x){
  
  start <- str_locate(x, "\\.(?=[^.]*$)")[1]
  end <-str_length(x)
  return(str_sub(x, start + 1, end))
  
}

parse_event_result <- function(results, result_type, depth = 1){
  
  if(is.null(results$coeftable)){
    dt <- data.table()
    if(depth >= 3){stop("empty list")}
    for(result in results){
      new_dt <- parse_event_result(result, result_type, depth + 1)
      dt <- rbind(dt, new_dt)
    }
    return(dt)
  }
  
  est <- results$coeftable
  dt <- data.table(variable = rownames(est), est, obs = results$nobs, outcome = deparse(results$fml[[2]]))
  
  
  dt[,result:=result_type]
  
  if(result_type %in% c("dynamic", "cohort_event_time")){
    dt[,event_time:= lapply(variable, get_event_time)]
    dt[,event_time := unlist(event_time)]
  }
  
  if(result_type == "cohort_event_time"){
    dt[,cohort:= lapply(variable, get_cohort)]
    dt[,cohort := unlist(cohort)]
  }
  
  dt[,stratify:= lapply(variable, get_stratify)]
  dt[,stratify := unlist(stratify)]
  
  return(dt)
  
}

# plot event dynamics -----------------------------------------------

plot_event_dynamics <-function(raw_dt, 
                               graphname = "event study plot", note = "", base_time = -1, significance_level = 0.05, stratify_offset =0.1){
  
  dt <- copy(raw_dt)
  
  suppressWarnings(dt[, c("variable", "t value", "Pr(>|t|)", "obs", "result") := NULL])
  setnames(dt, c("Estimate", "Std. Error"), c("att", "se"))
  
  #add the base period
  base_row <- data.table(att = 0, se = 0, event_time = base_time)
  dt <- dt |> add_base_row(base_row)
  
  
  
  dt[, conf_upb := att + qnorm(1-significance_level/2) * se]
  dt[, conf_lwb := att - qnorm(1-significance_level/2) * se]
  
  #add some offset
  if("stratify" %in% names(dt)){
    dt[, event_time := event_time + stratify_offset*(as.integer(stratify)-1)]
  }
  
  figure <- dt %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'red')
  
  if("stratify" %in% names(dt)){
    figure <- figure + geom_line(aes(x = event_time, y = att, color = stratify)) + 
      geom_point(aes(x = event_time, y = att, color = stratify)) +
      geom_errorbar(aes(x = event_time, ymin = conf_lwb, ymax = conf_upb), 
                    width = 0.1, linetype = "dashed")
  } else {
    figure <- figure + geom_line(aes(x = event_time, y = att), color = "black") + 
      geom_point(aes(x = event_time, y = att), color = "black") +
      geom_errorbar(aes(x = event_time, ymin = conf_lwb, ymax = conf_upb), 
                    width = 0.1, linetype = "dashed")
  }
  
  if("outcome" %in% names(dt)){
    figure <- figure + facet_wrap(~ outcome, scales = "free")
  }
  if("cohort" %in% names(dt)){
    figure <- figure + facet_wrap(~ cohort, scales = "free")
  }
  
  figure <- figure +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.background = element_rect(linetype = "dashed", color = "black"),
          legend.box = "horizontal",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0)) +
    labs(title = graphname, subtitle = note)
  
  return(figure)
  
}

add_base_row <- function(dt, base_row){
  
  col_names <- names(dt)
  opt_cols <- c("outcome", "stratify", "cohort")
  
  for(opt_col in opt_cols){
    if(opt_col %in% col_names){
      
      row_list <- list()
      
      for (col_name in dt[, unique(get(opt_col))]) {
        base_row_sub <- copy(base_row)
        base_row_sub[, c(opt_col) := col_name]
        row_list <- c(row_list, list(base_row_sub))
      }
      
    } 
    base_row <- rbindlist(row_list)
  }
  
  dt <- dt |> rbind(base_row)
  
  return(dt)
}

# utils ---------------------------------

get_cohort_size <- function(dt, cohortvar, unitvar) {
  cohort_obs <- dt[!is.infinite(get(cohortvar)), .(count = uniqueN(get(unitvar))) , by = cohortvar]
  setnames(cohort_obs, "G", "cohort")
  return(cohort_obs)
}


set_max_thread <- function(){
  setDTthreads(0)
  options(kit.nThread = getDTthreads())
  setFixest_nthreads(getDTthreads())
}

load_event_code_package <- function(){
  require(collapse)
  require(data.table)
  require(dreamerr)
  require(ggplot2)
  require(magrittr)
  require(rlist)
  require(stringr)
  require(fixest)
}

sim_did <- function(sample_size, time_period, untreated_prop = 0.3, 
                    cov = "no", hetero = "all", second_outcome = FALSE, na = "none", 
                    balanced = TRUE, seed = NA, stratify = FALSE, treatment_assign = "latent"){
  
  if(!is.na(seed)){set.seed(seed)}
  
  #unit  -------------
  dt_i <- data.table(unit = 1:sample_size)
  if(cov == "int"){
    dt_i[, x := sample.int(5, sample_size, replace = TRUE)] 
  } else if (cov == "no"){
    dt_i[, x := 1] 
  } else if (cov == "cont"){
    dt_i[, x := rnorm(sample_size)]
  }
  
  if(stratify){
    dt_i[, s := fifelse(rnorm(sample_size) > 0, 1, 2)]
  } else {
    dt_i[, s := 1]
  }
  
  #treatment assignment ---------------------
  
  #assign treated group based on a latent related to X
  
  if(treatment_assign == "latent"){
    dt_i[, treat_latent := x*0.2 + rnorm(sample_size)] #unit with larger X tend to be treated and treated earlier
    untreated_thres <- quantile(dt_i$treat_latent, untreated_prop)
    dt_i[treat_latent <= untreated_thres, G := Inf] #unit with low latent is never treated
    
    cohort_prop <- (1-untreated_prop)/(time_period-1)
    last_treat_thres <- untreated_thres
    for(t in time_period:2){ #unit start getting treated in t = 2
      treat_thres <- quantile(dt_i$treat_latent, untreated_prop + cohort_prop*(time_period - t + 1))
      dt_i[treat_latent <= treat_thres & treat_latent > last_treat_thres, G := t]
      last_treat_thres <- treat_thres
    }
    rm(t)
  } else if (treatment_assign == "uniform"){
    #when treatment is set to 'uniform', untreated propensity is fixed
    dt_i[,G := floor((unit-1)/(sample_size/time_period))]
    dt_i[G < 2, G := Inf]
  }
  
  #assign unit FE
  dt_i[, unit_fe := rnorm(sample_size)]
  
  #time ------------------
  
  dt_t <- data.table(time = 1:time_period)
  dt_t[, time_fe := rnorm(time_period)]
  dt_t[, x_trend := rnorm(time_period)]
  
  #panel --------------------------
  dt <- CJ(unit = 1:sample_size, time = 1:time_period)
  dt <- dt |> merge(dt_i, by = "unit")
  dt <- dt |> merge(dt_t, by = "time")
  
  dt[, D := as.integer(time >= G)]
  
  #untreated potential outcomes
  dt[, y0 := unit_fe + time_fe + x*x_trend + rnorm(sample_size*time_period, sd = 0.001)]
  
  #generate gtatt
  att <- CJ(G = 1:time_period, time = 1:time_period)
  if(hetero == "all"){
    att[, attgt := rnorm(time_period*time_period, mean = 2, sd = 1)]
  } else if (hetero == "dynamic"){
    for(event_t in 0:max(att[,time-G])){
      att[time - G == event_t, attgt := rnorm(1, mean = 2, sd = 1)]
    }
  }
  att[time < G, attgt := 0] #no anticipation
  
  dt <- dt |> merge(att, by = c("time", "G"), all.x = TRUE, all.y = FALSE)
  dt[is.na(attgt), attgt := 0]
  dt[, tau := attgt*s]
  
  #potential outcome
  dt[, y1 := y0 + tau]
  dt[, y := y1*D + y0*(1-D)]
  dt <- dt[, .SD, .SDcols = c("time", "G", "unit", "x", "y", "s")]
  
  #additional -----------------
  
  if(na == "y"){
    dt[, y := na_insert(y)]
  } else if (na == "x") {
    dt[, x := na_insert(x)]
  } else if( na == "both"){
    dt[, y := na_insert(y)]
    dt[, x := na_insert(x)]
  }
  
  if(balanced == FALSE){
    size <- fnrow(dt)
    dt <- dt[sample(1:size, size*0.99)]
  }
  if(second_outcome == TRUE){  dt[, y2 := y + 1 + rnorm(fnrow(dt), 0, 0.1)]}
  
  return(list(dt = dt, att = att))
  
}