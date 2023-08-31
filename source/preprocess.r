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
  
  # validation -----------------------------------
  
  if(balanced_panel==TRUE & !is.na(instrument) & !instrument_exposure%in%c("full","base")) stop("If imposing a balanced panel in an IV design, set instrument_exposure to full or base.")
  if(lower_event_time > base_time) stop("lower_event_time must lie below base_time")
  if(!is.data.table(maindata)) stop("rawdata must be a data.table")
  if(!all(is.na(stratify_balance_val)) & all(covariate_base_stratify == 1)) stop("It makes no sense to specify stratify_balance_val without specifying covariate_base_stratify")
  if(!instrument_exposure%in%c("full","partial","all","base")) stop("instrument_exposure, must be set to either full, partial, all, or base")

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
  
  # stacking for cohort -----------------------------------------------------------------
  
  treatdata<-copy(maindata[onset_age_v>=onset_minimum  & ! is.infinite(cohort) & !is.infinite(anycohort) ,])
  treatdata[,treated:=1]
  treatdata[,event_time:=time-cohort]
  treatdata<-treatdata[event_time >= lower_event_time & event_time <= upper_event_time,]
  treatdata[,treatgroup:="treated"]
  #stacking control cohorts.
  #I assume people who never suffer the event have a value cohort = Inf
  controldata<-NULL
  for(o in unique(treatdata$cohort)){
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
    
    controldata<-rbind(controldata,controlcohort)
  }
  rm(controlcohort)
  gc()
  
  controldata[,treated:=0]
  controldata<-controldata[event_time >= lower_event_time & event_time <= upper_event_time,]
  
  
  # checking observation -----------------------------------------------------------------
  
  
  #check if any unit is observed more then once in the base period
  treatdata[,obsbase:=sum(event_time==base_time),by=.(id,cohort)]
  controldata[,obsbase:=sum(event_time==base_time),by=.(id,cohort)]
  if(max(treatdata$obsbase)>1) stop("Error: some treated units are observed more than once in the reference period")
  if(max(controldata$obsbase)>1) stop("Error: some control units are observed more than once in the reference period")
  treatdata <- treatdata[obsbase==1,]
  controldata <- controldata[obsbase==1,]
  gc()
  
  #If base_time varies across units, reassigning it to a common reference value:
  #This is relevant, for instance, with a dataset that moves from annual to bi-annual
  treatdata[event_time == base_time,event_time :=max(base_time)]
  treatdata[,base_time :=max(base_time)]
  controldata[event_time == base_time,event_time :=max(base_time)]
  controldata[,base_time :=max(base_time)]
  
  treatdata[,obscount:=1]
  controldata[,obscount:=1]
  treatdata[,obscount:=sum(obscount),by=.(id,cohort)]
  controldata[,obscount:=sum(obscount),by=.(id,cohort)]
  
  if(balanced_panel==TRUE){
    numperiods<-length(unique(treatdata$event_time))
    treatdata<- treatdata[obscount== numperiods,]
    controldata <- controldata[obscount>=numperiods,]
  }
  
  #turn covariates to factor
  covariates <- c(covariate_base_stratify, covariate_base_balance, covariate_base_support)
  covariates <- covariates[covariates %!=% 1]
  for(out in covariates){
    controldata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
    controldata[get(out) >= 9e9,eval(out) := NA, ]
    controldata[,eval(out) := as.factor(get(out))]
    treatdata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
    treatdata[get(out) >= 9e9,eval(out) := NA, ]
    treatdata[,eval(out) := as.factor(get(out))]
  }

  controldata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
  controldata <- controldata[base_restrict == 1,]
  treatdata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
  treatdata <- treatdata[base_restrict == 1,]

  treatdata[,base_restrict_treated := max(base_restrict_treated * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
  treatdata <- treatdata[base_restrict_treated == 1,]
  controldata <- controldata[, base_restrict_treated := 1]
  
  #turn covariates interaction to factors
  for(covariate_type in c("stratify", "balancevars", "balancevars_linear_subset", "supportvars")){
    
    cov_vars <- switch(covariate_type, 
                       stratify = ifelse(stratify_by_cohort, c(covariate_base_stratify, "cohort"), covariate_base_stratify),
                       balancevars = covariate_base_balance,
                       balancevars_linear_subset = covariate_base_balance_linear_subset,
                       supportvars = covariate_base_support)
    
    if(is.character(cov_vars)){
      treatdata[,(covariate_type) :=interaction(treatdata[, cov_vars, with = FALSE], drop=TRUE)]
      controldata[,(covariate_type) :=interaction(controldata[, cov_vars, with = FALSE], drop=TRUE)]
    } else {
      treatdata[,(covariate_type) := factor(1,levels=c(1,"OMIT"))]
      controldata[,(covariate_type) := factor(1,levels=c(1,"OMIT"))]
    }
  
  }
  
  #check common support
  treatdata[,temp:=interaction(balancevars,supportvars,stratify,event_time,cohort,drop=TRUE)]
  controldata[,temp:=interaction(balancevars,supportvars,stratify,event_time,cohort,drop=TRUE)]
  commonvals<-intersect(unique(treatdata$temp),unique(controldata$temp))
  treatdata<-treatdata[temp%in%commonvals,]
  controldata<-controldata[temp%in%commonvals,]
  
  treatdata[,temp:=NULL]
  controldata[,temp:=NULL]
  rm(commonvals)
  gc()
  
  # stacking for event_time -----------------------------------------------------------------

  event_times<-treatdata[,unique(event_time)]
  eventdata<-NULL

  for(t in event_times){
    treatdata[,obst:=sum(event_time==t),by=.(id,cohort)]
    controldata[,obst:=sum(event_time==t),by=.(id,cohort)]
    
    pairdata<-rbind(treatdata[obsbase==1 & obst==1 & base_time != t & (event_time == t | event_time == base_time),],
                    controldata[obsbase==1 & obst==1 & base_time != t & (event_time == t | event_time == base_time),])
    pairdata[,time_pair:= t]
    eventdata<-rbind(eventdata, pairdata)
  }

  # estimating ipw ----------------------------------------------------------------------------------
  
  if(is.null(eventdata)==T){stop("eventdata is empty!!!")}

  eventdata[,obst:=NULL]
  eventdata[,obsbase:=NULL]
  eventdata[,anycohort:=NULL]
  
  rm(treatdata)
  rm(controldata)
  gc()
  
  eventdata[,post:=event_time >= 0]
  
  #turn the cols into factors
  factor_cols <- c("id", "treatgroup", "cohort", "event_time", "time_pair", "time")
  eventdata[, (factor_cols) := lapply(.SD, as.factor), .SDcols = factor_cols]

  basefactor<-unique(eventdata[event_time==base_time,event_time])
  basefactor<-basefactor[length(basefactor)]
  eventdata[,event_time:=relevel(event_time,ref=as.character(basefactor))]
  
  #calculating weights so that controls match treated households on characteristics
  #Note: this may have a lot of fixed effects, and may need to be broken down into 
  #multiple smaller regressions:

  eventdata[,pweight:=NULL]
  if(is.character(covariate_base_balance_linear))
    {eventdata[,pval:=feols(as.formula(paste0("treated ~",paste0(covariate_base_balance_linear,collapse="+",
                                                                 sep=":interaction(cohort,event_time,time_pair,stratify,balancevars_linear_subset, drop = TRUE)"),
                                              "| interaction(cohort,event_time,time_pair,stratify,balancevars, drop = TRUE)"
    )),data = eventdata, lean = FALSE)$fitted.values]
    }
  else{
    eventdata[,pval:=feols(treated ~ 1 | interaction(cohort,event_time,time_pair,stratify,balancevars, drop = TRUE),
                                  data = eventdata, lean = FALSE)$fitted.values]
  }
  
  eventdata[treated==1 & pval < 1 & pval > 0,pweight:=1]
  eventdata[treated==0 & pval < 1 & pval > 0,pweight:=pval/(1-pval)]
  eventdata[,pval:=NULL]
  eventdata<- eventdata[!is.na(pweight),]

  #  deal with extensions ----------------------------------------------
  
  if(!is.na(instrument)){
    eventdata <- process_iv(eventdata, instrument, instrument_exposure, covs_instrument_base_balance, saturate)
  }
  
  if(!is.na(stratify_balance_val)){
    eventdata <- process_stratify(eventdata)
  }
  
  eventdata[,treated:=as.factor(treated)]
  
  return(eventdata)

}


construct_event_variables<-function(eventdata,saturate=FALSE,IV=FALSE,response=NULL){
  #These regressions should work identically if the fixed effects (after the "|") were replaced with:
  # interaction(time_pair,id,cohort)
  eventdata[,treated_post := as.factor((treated == 1) * (post == 1))]
  eventdata[,treated_pre := as.factor((treated == 1) * (post == 0))]
  eventdata[,treated_event_time := event_time]
  eventdata[treated==0,treated_event_time := 1] #1 is the base level
  
  eventdata[,event_time_stratify:=interaction(event_time,stratify, drop = TRUE)]
  eventdata[,treated_post_stratify := interaction(treated_post,stratify, drop = TRUE)]
  
  eventdata[,treated_pre_stratify := interaction(treated_pre,stratify,drop=TRUE)]
  eventdata[treated_pre==0,treated_pre_stratify := 1]
  eventdata[event_time==base_time,treated_pre_stratify := 1]
  
  if(IV==FALSE) eventdata[,unitfe := interaction(time_pair,treated,stratify, drop = TRUE)]
  else eventdata[,unitfe := interaction(time_pair,id,treated,cohort,stratify, drop = TRUE)]
  
  eventdata[,treated_event_time_stratify := interaction(event_time,stratify, drop = TRUE)]

  #Omitting base year for all levels of --stratify--:
  eventdata[event_time==base_time,event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  
  #Omitting base year for all levels of --stratify--, for treated people
  eventdata[treated == 0 ,treated_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  eventdata[event_time==base_time,treated_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  eventdata[,treated_event_time_stratify:=relevel(treated_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  

  #Omitting effect for untreated people or observations in pre-period:
  eventdata[treated_post == 0 ,treated_post_stratify := paste0(c(0,1),collapse=".")]
  eventdata[,treated_post_stratify:=relevel(treated_post_stratify,ref = paste0(c(0,1),collapse="."))]
  
  # IV -------------------------------------------------------------------------------
  
  if(IV==TRUE){
  eventdata[,Z:=instrument_group[event_time != base_time],by=.(unitfe)]
  eventdata[,R:=get(response)[event_time != base_time],by=.(unitfe)]
  eventdata[,response_event_time_stratify := interaction(event_time,stratify, drop = TRUE)]
  eventdata[,instrument_event_time_stratify := interaction(event_time,stratify, drop = TRUE)]
  #Omitting base year for all levels of --stratify--, for response interaction
  eventdata[R == 0 ,response_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  eventdata[event_time==base_time,response_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  eventdata[,response_event_time_stratify:=relevel(response_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  #Omitting base year for all levels of --stratify--, for instrument interaction
  eventdata[Z == 0 ,instrument_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  eventdata[event_time==base_time,instrument_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  eventdata[,instrument_event_time_stratify:=relevel(instrument_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  #Defining response-event-time interactions, setting to omitted value if response == 0:
  eventdata[,treated_response_event_time_stratify := treated_event_time_stratify]
  eventdata[R==0,treated_response_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  
  #Defining pooled response interactions, setting to omitted value if response == 0:
  eventdata[,treated_response_post_stratify := treated_post_stratify]
  eventdata[R==0,treated_response_post_stratify := paste0(c(0,1),collapse=".")]
  eventdata[,treated_response_pre_stratify := treated_pre_stratify]
  eventdata[R==0,treated_response_pre_stratify := 1]
  
  #Defining instrument-event-time interactions, setting to omitted value if instrument == 0:
  eventdata[,treated_instrument_event_time_stratify := treated_event_time_stratify]
  eventdata[Z==0,treated_instrument_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  
  
  #Defining pooled instrument interactions, setting to omitted value if instrument == 0:
  eventdata[,treated_instrument_post_stratify := treated_post_stratify]
  eventdata[Z==0,treated_instrument_post_stratify := paste0(c(0,1),collapse=".")]
  eventdata[,treated_instrument_pre_stratify := treated_pre_stratify]
  eventdata[Z==0,treated_instrument_pre_stratify := 1]
  
  if(saturate==TRUE) {
    eventdata[,event_time_stratify:=interaction(event_time_stratify,balancevars, drop=TRUE)]
    eventdata[,instrument_event_time_stratify:=interaction(instrument_event_time_stratify,balancevars, drop=TRUE)]
    eventdata[,treated_event_time_stratify:=interaction(treated_event_time_stratify,balancevars, drop=TRUE)]
    eventdata[,treated_pre_stratify:=interaction(treated_pre_stratify,balancevars, drop=TRUE)]
    eventdata[,treated_post_stratify:=interaction(treated_post_stratify,balancevars, drop=TRUE)]
    # if(fullsaturate==TRUE){
    # eventdata[,treated_instrument_event_time_stratify:=interaction(treated_instrument_event_time_stratify,balancevars, drop=TRUE)]
    # eventdata[,treated_instrument_event_time_stratify:=relevel(treated_instrument_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1,0),collapse="."))]
    # eventdata[event_time==base_time | Z == 0 | treated==0,treated_instrument_event_time_stratify:=levels(treated_instrument_event_time_stratify)[1]]
    # eventdata[,treated_instrument_post_stratify:=interaction(treated_instrument_post_stratify,balancevars, drop=TRUE)]
    # eventdata[,treated_instrument_pre_stratify:=interaction(treated_instrument_pre_stratify,balancevars, drop=TRUE)]
    # }
  }
  }
  
  return(eventdata)
}

process_stratify <- function(eventdata){
  if(stratify_balance_val != "mean"){
    eventdata[,treated_base := treated == 1 & stratify == stratify_balance_val]
    
    stratvals<-unique(eventdata$stratify)
    for(strat in stratvals){
      if(strat != stratify_balance_val){
        stratbalmodel <- feols(treated_base ~ 1 | interaction(cohort,event_time,time_pair,balancevars, drop = TRUE),
                               data = eventdata[(treated == 1 & stratify == stratify_balance_val) | 
                                                  (treated == 1 & stratify == strat),], lean = TRUE, mem.clean=TRUE)
        eventdata[, eval(paste0("pval_1",strat)) := predict(stratbalmodel,
                                                            eventdata)]
        
        eventdata[treated == 1 & stratify == strat, pweight_stratbal := get(paste0("pval_1",strat))/(1-get(paste0("pval_1",strat)))]
      }
      else    eventdata[treated == 1 & stratify == strat, pweight_stratbal := 1]
      
      stratbalmodel <- feols(treated_base ~ 1 | interaction(cohort,event_time,time_pair,balancevars, drop = TRUE),
                             data = eventdata[(treated == 1 & stratify == stratify_balance_val) | 
                                                (treated == 0 & stratify == strat),], lean = TRUE, mem.clean=TRUE)
      
      eventdata[, eval(paste0("pval_0",strat)) := predict(stratbalmodel,
                                                          eventdata)]
      eventdata[treated == 0 & stratify == strat, pweight_stratbal := get(paste0("pval_0",strat))/(1-get(paste0("pval_0",strat)))]
      
      
    }
    
    checkvars <- names(eventdata)[grepl("pval_",names(eventdata))]
    for(var in checkvars){
      eventdata[get(var) == 1 | get(var) == 0 | is.na(get(var)),  pweight_stratbal := NA ]
    }
  }
  if(stratify_balance_val == "mean"){
    
    stratvals<-unique(eventdata$stratify)
    for(strat in stratvals){
      eventdata[,treated_base := treated == 1 & stratify == strat]
      stratbalmodel <- feols(treated_base ~ 1 | interaction(cohort,event_time,time_pair,balancevars, drop = TRUE),
                             data = eventdata, lean = TRUE, mem.clean=TRUE)
      eventdata[, eval(paste0("pval_1",strat)) := predict(stratbalmodel,
                                                          eventdata)]
      
      eventdata[treated == 1 & stratify == strat, pweight_stratbal :=1/(get(paste0("pval_1",strat)))]
      
      eventdata[,treated_base := treated == 0 & stratify == strat]
      stratbalmodel <- feols(treated_base ~ 1 | interaction(cohort,event_time,time_pair,balancevars, drop = TRUE),
                             data = eventdata, lean = TRUE, mem.clean=TRUE)
      
      eventdata[, eval(paste0("pval_0",strat)) := predict(stratbalmodel,
                                                          eventdata)]
      eventdata[treated == 0 & stratify == strat, pweight_stratbal := 1/(get(paste0("pval_0",strat)))]
      
      
      
    }
    #Need to rebalance everyone to look like distribution of TREATED guys:
    stratbalmodel <- feols(treated ~ 1 | interaction(cohort,event_time,time_pair,balancevars, drop = TRUE),
                           data = eventdata, lean = TRUE, mem.clean=TRUE)
    eventdata[, pval_treat := predict(stratbalmodel,
                                      eventdata)]
    eventdata[,pweight_stratbal:=pweight_stratbal*pval_treat]
    
    
    checkvars <- names(eventdata)[grepl("pval_",names(eventdata))]
    for(var in checkvars){
      eventdata[get(var) == 1 | get(var) == 0 | is.na(get(var)),  pweight_stratbal := NA ]
    }
  }
  return(eventdata)
}

process_iv <- function(eventdata, instrument, instrument_exposure, covs_instrument_base_balance, saturate){
  
  eventdata[,instrument_group_now:=max(get(instrument) == 1 & as.numeric(as.character(time_pair)) == as.numeric(as.character(event_time))),by=.(id,cohort,time_pair)]
  eventdata[,instrument_group_base:=max(get(instrument) == 1 & as.numeric(as.character(base_time)) == as.numeric(as.character(event_time))),by=.(id,cohort,time_pair)]
  if(instrument_exposure=="full") eventdata <- eventdata[instrument_group_now == instrument_group_base,]
  if(instrument_exposure=="partial") eventdata <- eventdata[(instrument_group_now == 0 & instrument_group_base == 0) | (instrument_group_base == 0 & instrument_group_now == 1),]
  if(instrument_exposure=="all") {
    eventdata <- eventdata[!is.na(instrument_group_now) & !is.na(instrument_group_base),]
    eventdata[,instrument_group_now := pmax(instrument_group_now,instrument_group_base)]
  }
  if(instrument_exposure=="base") {
    eventdata <- eventdata[!is.na(instrument_group_base),]
    eventdata[,instrument_group_now := instrument_group_base]
  }
  #After restricting on instrument status, need to re-impose balance:
  eventdata[,obscount:=1]
  eventdata[,obscount:=sum(ifelse(event_time==base_time,0,obscount)),by=.(id,cohort)]
  if(balanced_panel==TRUE) {
    eventdata <- eventdata[obscount==numperiods-1,]
  }
  
  eventdata[,instrument_group:=instrument_group_now]
  eventdata[,instrument_group_now:=NULL]
  eventdata[,instrument_group_base:=NULL]
  if(length(unique(eventdata$instrument_group)) == 1) 
    stop("No variation in instrument left in the panel. Maybe because you're trying to impose panel balance?")
  
  
  if(!is.na(covs_instrument_base_balance)) eventdata[,instrument_balancevars:=interaction(eventdata[,covs_instrument_base_balance,with=FALSE], drop=TRUE)]
  else eventdata[,instrument_balancevars:=factor(1,levels=c(1,"OMIT"))]
  
  eventdata[,ivgroup:=interaction(treated,instrument_group)]
  ivlevels<-unique(eventdata$ivgroup)
  for(strat in ivlevels){
    eventdata[,treated_base := ivgroup == strat]
    stratbalmodel <- feols(treated_base ~ 1 | interaction(cohort,event_time,time_pair,instrument_balancevars, drop = TRUE),
                           data = eventdata, lean = TRUE, mem.clean=TRUE,
                           weights = eventdata$pweight)
    eventdata[, eval(paste0("pval_",strat)) := predict(stratbalmodel,
                                                       eventdata)]
    
    eventdata[ivgroup == strat, pweight_new :=pweight/(get(paste0("pval_",strat)))]
    
  }
  #Need to rebalance everyone to look like distribution of TREATED guys:
  stratbalmodel <- feols(treated ~ 1 | interaction(cohort,event_time,time_pair,instrument_balancevars, drop = TRUE),
                         data = eventdata, lean = TRUE, mem.clean=TRUE,
                         weights = eventdata$pweight)
  eventdata[, pval_treat := predict(stratbalmodel,
                                    eventdata)]
  
  eventdata[,pweight_new:=pweight_new*pval_treat]
  eventdata[,pval_treat:=NULL]
  eventdata[,pweight:=pweight_new]
  eventdata[,pweight_new:=NULL]
  
  checkvars <- names(eventdata)[grepl("pval_",names(eventdata))]
  for(var in checkvars){
    eventdata[get(var) == 1 | get(var) == 0 | is.na(get(var)),  pweight := NA ]
    eventdata[,eval(var):=NULL]
  }
  eventdata<-eventdata[!is.na(pweight),]
  if(length(unique(eventdata$instrument_group)) == 1) stop("No variation in instrument left in the panel after imposing common support on observables
                                                             across instrument X treatment cells.")
  return(eventdata)
  
}



