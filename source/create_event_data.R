# eventcode

#Treated households:
create_event_data<-function(maindata,
                            #if no covariates for balancing or stratification, just specify a constant (as follows):
                            covariate_base_stratify=1, #vector of variable names, treatment effects are stratified by their joint unique values.
                            covariate_base_balance=1, #vector of variable names to include in finding exact matches for treated units
                            covariate_base_balance_linear=1, #vector of variable names to include linearly in propensity score estimator
                            covariate_base_balance_linear_subset=1, #vector of variable names within which the linear propensity score regression will be estimated
                            covariate_base_support=1, #vector of variable names on which to impose common support across treated and control units
                            min_control_gap = 1, #controls must be treated at least this many years away from matched treated units
                            max_control_gap = Inf, #controls must be treated no more than this many years away from matched treated units
                            instrument=NA,
                            covs_instrument_base_balance=NA, #Balance these characteristics jointly across BOTH treatment and instrument status.
                            instrument_exposure="full",
                            base_restrict = 1, #a single var, restricting to households for which var == 1 in base period,
                            base_restrict_treated = 1, #a single var, restricting TREATED GUYS to households for which var == 1 in base period,
                            timevar,
                            unitvar,
                            cohortvar, #period when the particular event of interest arises
                            anycohortvar = NULL, #period when ANY kind of event arises, used to find control cohorts
                            stratify_by_cohort=FALSE,
                            #(for whom no event of any kind has yet to happen)
                            #Only specify if it differs from cohortvar
                            onset_agevar=NULL, #variable indicating age at which LTC needs arose, for person affected by them
                            base_time = -1, #reference value for event time; the period from which differences over time are taken
                            onset_minimum=-Inf, #Drop treated guys whose onset happens BEFORE this time period
                            onset_maximum=Inf, #Drop treated guys whose onset happens AFTER this time period
                            never_treat_action = "both", #options: "both", "only", and "exclude"
                            #Both: include both never-treated and later-treated units as controls
                            #Only: use only never-treated units as controls
                            #Exclude: use only later-treated units as controls
                            lower_event_time = -Inf, #Earliest period (relative to treatment time) for which to estimate effects
                            upper_event_time = Inf, #Latest period (relative to treatment time) for which to estimate effects
                            balanced_panel = FALSE, #If TRUE, keep only units observed over full interval [lower_event_time, upper_event_time]
                            stratify_balance_val = NA, #Set to "mean" to balance each strata on covariates to look like the overall treated sample
                            #Set to a value of "stratify" to balance all other strata to look like the given one.
                            check_not_treated = FALSE #If TRUE, only allows controls to be used for cohort o if they are observed in the dataset to be untreated in year o.
                            ) 
{

  # name change ----------------------------------
  
  if(is.null(anycohortvar)) {
    maindata[,anycohortvar:=get(cohortvar)]
    anycohortvar<-"anycohortvar"
  }
  if(is.numeric(base_time)){
    maindata[,base_time:=base_time]
  } else {
    setnames(maindata,base_time,"base_time")
  }
  if(!is.null(onset_agevar)) {
    maindata[,onset_age_v:=get(onset_agevar)]
  } else {
    maindata[,onset_age_v:=Inf]
  }
  setnames(maindata,c(timevar,unitvar,cohortvar,anycohortvar),c("time","id","cohort","anycohort"))
  if(base_restrict != 1) setnames(maindata,base_restrict,"base_restrict")
  if(base_restrict_treated != 1) setnames(maindata,base_restrict_treated,"base_restrict_treated")
  if(any(is.na(maindata$cohort))) stop("cohort variable should not be missing (it can be infinite instead)")
  if(any(is.na(maindata$anycohort))) stop("anycohort variable should not be missing (it can be infinite instead)")
  if(!is.character(covariate_base_balance_linear_subset)) covariate_base_balance_linear_subset <- covariate_base_balance
  
  # validation -----------------------------------
  
  #if(balanced_panel==TRUE & !is.na(instrument) & !instrument_exposure%in%c("full","base")) stop("If imposing a balanced panel in an IV design, set instrument_exposure to full or base.")
  if(lower_event_time > base_time) stop("lower_event_time must lie below base_time")
  if(!is.data.table(maindata)) stop("rawdata must be a data.table")
  if(!all(is.na(stratify_balance_val)) & all(covariate_base_stratify == 1)) stop("It makes no sense to specify stratify_balance_val without specifying covariate_base_stratify")
  if(!instrument_exposure%in%c("full","partial","all","base")) stop("instrument_exposure, must be set to either full, partial, all, or base")

  if(balanced_panel){
    
    #check if any is dup
      
    
    if(anyDuplicated(maindata[, .(id, time)])){
      dup <- duplicated(maindata[,.(id, time)])
      warning(nrow(dup_id), " units is observed more than once in some periods, enforcing balanced panel by dropping them")
      maindata <- maindata[!id %fin% dup]
    }
    
    #check if any is missing
    id_count <- maindata[, .(count = .N), by = id]
    time_period <- maindata[, uniqueN(time)]
    if(any(id_count[, count < time_period])){
      mis_id <- id_count[count < time_period]
      warning(nrow(mis_id), " units is missing in some periods, enforcing balanced panel by dropping them")
      maindata <- maindata[!id %fin% mis_id$id]
    }
    
   #can't just check count = period because one may be missing in one period and observed multiple time in another
    
  }

  # stacking for cohort -----------------------------------------------------------------
  
  treatdata<-copy(maindata[onset_age_v>=onset_minimum  & ! is.infinite(cohort) & !is.infinite(anycohort) ,])
  treatdata[,treated:=1]
  treatdata[,event_time:=time-cohort]
  treatdata<-treatdata[event_time >= lower_event_time & event_time <= upper_event_time,]
  treatdata[,treatgroup:="treated"]
  
  #stacking control cohorts.
  #I assume people who never suffer the event have a value cohort = Inf
  control_list <- list()
  for(o in unique(treatdata$cohort)){
    
    #find the relevant control group
    if(never_treat_action=="both") {
      controlcohort <- maindata[(anycohort > o | is.infinite(anycohort)) ,]
      controlcohort[is.infinite(anycohort),treatgroup:="never-treated"]
      controlcohort[!is.infinite(anycohort),treatgroup:="later-treated"]
    }
    if(never_treat_action=="exclude"){
      controlcohort <- maindata[(anycohort > o & !is.infinite(anycohort)) ,]
      controlcohort[,treatgroup:="later-treated"]
    }
    if(never_treat_action=="only") {
      controlcohort <- maindata[is.infinite(anycohort) ,]
      controlcohort[,treatgroup:="never-treated"]
    }
    
    if (!is.infinite(max_control_gap)) controlcohort <- controlcohort[anycohort - o <=  max_control_gap,]
    if (min_control_gap != 1)  controlcohort <- controlcohort[anycohort - o >=  min_control_gap,]
    
    controlcohort[,cohort := o]
    controlcohort[,event_time := time - cohort]
    
    #Make sure people in the control cohort are actually observed in that period
    #(to verify they don't belong to the cohort)
    if(check_not_treated == TRUE){
    controlcohort[ ,obscohort := max(time == o),by=id]
    controlcohort<-controlcohort[obscohort==1,]
    controlcohort[,obscohort:=NULL]
    }
    #drop someone from the control cohort when they get treated:
    controlcohort<-controlcohort[anycohort - cohort > event_time ,]
    
    control_list <- c(control_list, list(controlcohort))
    
  }
  rm(controlcohort)
  
  controldata<-rbindlist(control_list)
  
  gc()
  
  controldata[,treated:=0]
  controldata<-controldata[event_time >= lower_event_time & event_time <= upper_event_time,]
  
  # covariate to factor -----------------------------------------------------------
  
  #turn covariates to factor
  covariates <- c(covariate_base_stratify, covariate_base_balance, covariate_base_support)
  covariates <- covariates[covariates %!=% 1]
  for(out in covariates){
    treatdata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
    treatdata[get(out) >= 9e9,eval(out) := NA, ]
    treatdata[,eval(out) := qF(get(out))]
    controldata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
    controldata[get(out) >= 9e9,eval(out) := NA, ]
    controldata[,eval(out) := qF(get(out))]
  }
  
  #turn covariates interaction to factors
  for(covariate_type in c("stratify", "balancevars", "balancevars_linear_subset", "supportvars")){
    
    cov_vars <- switch(covariate_type, 
                       stratify = ifelse(stratify_by_cohort, c(covariate_base_stratify, "cohort"), covariate_base_stratify),
                       balancevars = covariate_base_balance,
                       balancevars_linear_subset = covariate_base_balance_linear_subset,
                       supportvars = covariate_base_support)
    if(is.character(cov_vars)){
      treatdata[,(covariate_type) :=  do.call(finteraction, treatdata[, cov_vars, with = FALSE])]
      controldata[,(covariate_type) :=  do.call(finteraction, controldata[, cov_vars, with = FALSE])]
    } else {
      treatdata[,(covariate_type) := factor(1,levels=c(1,"OMIT"))]
      controldata[,(covariate_type) := factor(1,levels=c(1,"OMIT"))]
    }
  }
  
  # checking observation -----------------------------------------------------------------

  #If base_time varies across units, reassigning it to a common reference value:
  #This is relevant, for instance, with a dataset that moves from annual to bi-annual
  treatdata[event_time == base_time,event_time :=max(base_time)]
  controldata[event_time == base_time,event_time :=max(base_time)]
  
  base_time <- max(treatdata[, max(base_time)], controldata[, max(base_time)])
  
  treatdata[,base_time := NULL]
  controldata[,base_time := NULL]
  
  #check if any unit is observed more then once in the base period
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
  if(base_restrict != 1){
    controldata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
    controldata <- controldata[base_restrict == 1,]
    treatdata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
    treatdata <- treatdata[base_restrict == 1,]
    treatdata[,base_restrict_treated := max(base_restrict_treated * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
    treatdata <- treatdata[base_restrict_treated == 1,]
    controldata <- controldata[base_restrict_treated == 1]
  }
  
  #check common support
  if(is.character(covariate_base_balance_linear)){
    
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
  
  # stacking for event_time ----------------------------------------------------------------

  #if is balanced panel, after knowing its max and min, can be sure it is observed when in the middle
  # TODO: make sure the estimates is valid if there are missing value within the min max (FEOLS)
  
  controldata[, `:=`(min_event_time = min(event_time),
                     max_event_time = max(event_time)), by = .(id, cohort)]
  treatdata[, `:=`(min_event_time = min(event_time),
                   max_event_time = max(event_time)), by = .(id, cohort)]
  
  event_times<-treatdata[,funique(event_time)]
  data_list <- list()
  for(t in event_times){
    
    pair_treat_data <- treatdata[t >= min_event_time & t <= max_event_time & (event_time == t | event_time == base_time)]
    pair_control_data <- controldata[t >= min_event_time & t <= max_event_time & (event_time == t | event_time == base_time)]
    
    pair_treat_data[,time_pair := t]
    pair_control_data[,time_pair := t]
    
    data_list<-c(data_list, list(pair_treat_data), list(pair_control_data))
    
  }
  eventdata <- rbindlist(data_list,use.names=TRUE)
  
  # estimating ipw ----------------------------------------------------------------------------------
  
  if(is.null(eventdata)) {stop("eventdata is empty!")}

  eventdata[,anycohort:=NULL]
  
  rm(treatdata)
  rm(controldata)
  gc()
  
  eventdata[,post:=event_time >= 0]
  
  #turn the cols into factors
  factor_cols <- c("id", "treatgroup", "cohort", "time_pair", "time")
  for(col in factor_cols){
    eventdata[, (col) := qF(get(col))]
  }
  #keep a numeric version for later comparisons
  eventdata[, event_time_fact := qF(event_time)]

  #construct the call  
  if(!is.character(covariate_base_balance_linear)){
    call <- "treated ~ 1 | finteraction(cohort,event_time_fact,time_pair,stratify,balancevars)"
  } else {
    call <- paste0("treated ~",paste0(covariate_base_balance_linear,collapse="+",
                                      sep="finteraction(cohort,event_time_fact,time_pair,stratify,balancevars_linear_subset,)"),
                   "| finteraction(cohort,event_time_fact,time_pair,stratify,balancevars)")
  }
  
  #estimate ipw
  eventdata[,pval:= feols(as.formula(call), data = eventdata, lean = FALSE)$fitted.values]

  #only keep propensity score between 0,1 is equivalent to checking common support
  eventdata <- eventdata[pval < 1 & pval > 0]
  eventdata[,pweight := ifelse(treated == 1, 1, pval/(1-pval))]
  eventdata[,pval:=NULL]

  #  deal with extensions ----------------------------------------------
  
  if(!is.na(instrument)){
    eventdata <- eventdata |> process_iv(instrument, instrument_exposure, covs_instrument_base_balance, saturate)
  }
  
  if(!is.na(stratify_balance_val)){
    eventdata <- eventdata |> process_stratify(stratify_balance_val)
  }
  
  eventdata[,treated:=qF(treated)]
  
  # construct_event_variables --------------------------------------------
  
  if(!is.na(instrument)){
    eventdata <- construct_event_variables_iv(eventdata)
  }
  
  return(eventdata)

}





