# eventcode

library(data.table)
library(tidyverse)
library(fixest)
library(Hmisc)
library(lfe)
library(rlist)
library(car)

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
  if(balanced_panel==TRUE & !is.na(instrument) & !instrument_exposure%in%c("full","base")) stop("If imposing a balanced panel in an IV design, set instrument_exposure to full or base.")
  if(lower_event_time > base_time) stop("lower_event_time must lie below base_time")
  if(!is.data.table(maindata)) stop("rawdata must be a data.table")
  if(!all(is.na(stratify_balance_val)) & all(covariate_base_stratify == 1)) stop("It makes no sense to specify stratify_balance_val without specifying covariate_base_stratify")
  if(!instrument_exposure%in%c("full","partial","all","base")) stop("instrument_exposure, must be set to either full, partial, all, or base")
  #maindata<-copy(rawdata)
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
  
  #constructing the regression dataset, stacking each non-base year 
  #with a copy of the base year observation.
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
  
  if(is.character(covariate_base_stratify)){
    for(out in covariate_base_stratify){
      controldata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
      controldata[get(out) >= 9e9,eval(out) := NA, ]
      controldata[,eval(out) := as.factor(get(out))]
      treatdata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
      treatdata[get(out) >= 9e9,eval(out) := NA, ]
      treatdata[,eval(out) := as.factor(get(out))]
    }
  }
  if(is.character(covariate_base_balance)){
    for(out in covariate_base_balance){
      controldata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
      controldata[get(out) >= 9e9,eval(out) := Inf, ]
      controldata[,eval(out) := as.factor(get(out))]
      treatdata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
      treatdata[get(out) >= 9e9,eval(out) := Inf, ]
      treatdata[,eval(out) := as.factor(get(out))]
    }
  }
  
  if(is.character(covariate_base_support)){
    for(out in covariate_base_support){
      controldata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
      controldata[get(out) >= 9e9,eval(out) := Inf, ]
      controldata[,eval(out) := as.factor(get(out))]
      treatdata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
      treatdata[get(out) >= 9e9,eval(out) := Inf, ]
      treatdata[,eval(out) := as.factor(get(out))]
    }
  }
  

  controldata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
  controldata <- controldata[base_restrict == 1,]
  treatdata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
  treatdata <- treatdata[base_restrict == 1,]

  treatdata[,base_restrict_treated := max(base_restrict_treated * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
  treatdata <- treatdata[base_restrict_treated == 1,]
  controldata <- controldata[, base_restrict_treated := 1]
  
  if(is.character(covariate_base_stratify) & stratify_by_cohort == FALSE) {
    treatdata[,stratify:=interaction(treatdata[,covariate_base_stratify,with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
    controldata[,stratify:=interaction(controldata[,covariate_base_stratify,with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
  } else if(is.character(covariate_base_stratify) & stratify_by_cohort == TRUE) {
    treatdata[,stratify:=interaction(treatdata[,c(covariate_base_stratify,"cohort"),with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
    controldata[,stratify:=interaction(controldata[,c(covariate_base_stratify,"cohort"),with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
  } else {
    treatdata[,stratify:=factor(1,levels=c(1,"OMIT"))]
    controldata[,stratify:=factor(1,levels=c(1,"OMIT"))]
  }
  if(is.character(covariate_base_balance)) {
    treatdata[,balancevars:=interaction(treatdata[,covariate_base_balance,with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
    controldata[,balancevars:=interaction(controldata[,covariate_base_balance,with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
  }else {
    treatdata[,balancevars:=factor(1,levels=c(1,"OMIT"))]
    controldata[,balancevars:=factor(1,levels=c(1,"OMIT"))]
  }
  if(is.character(covariate_base_balance_linear_subset)) {
    treatdata[,balancevars_linear_subset:=interaction(treatdata[,covariate_base_balance_linear_subset,with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
    controldata[,balancevars_linear_subset:=interaction(controldata[,covariate_base_balance_linear_subset,with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
  }else {
    treatdata[,balancevars:=factor(1,levels=c(1,"OMIT"))]
    controldata[,balancevars:=factor(1,levels=c(1,"OMIT"))]
  }
  
  if(is.character(covariate_base_support)) {
    treatdata[,supportvars:=interaction(treatdata[,covariate_base_support,with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
    controldata[,supportvars:=interaction(controldata[,covariate_base_support,with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
  }else {
    treatdata[,supportvars:=factor(1,levels=c(1,"OMIT"))]
    controldata[,supportvars:=factor(1,levels=c(1,"OMIT"))]
  }
  
  treatdata[,temp:=interaction(balancevars,supportvars,stratify,event_time,cohort,drop=TRUE)]
  controldata[,temp:=interaction(balancevars,supportvars,stratify,event_time,cohort,drop=TRUE)]
  
  commonvals<-intersect(unique(treatdata$temp),unique(controldata$temp))
  
  treatdata<-treatdata[temp%in%commonvals,]
  controldata<-controldata[temp%in%commonvals,]
  
  treatdata[,temp:=NULL]
  controldata[,temp:=NULL]
  rm(commonvals)
  gc()
  
  event_times<-treatdata[,unique(event_time)]
  eventdata<-NULL
  #For the final dataset, we must stack each pairwise combo of (base year, other year)
  #for each household. The reason? We will assign control households weights that vary
  #based on the other year, as households enter/exit the sample.
#  if(balanced_panel==FALSE){
    for(t in event_times){
      treatdata[,obst:=sum(event_time==t),by=.(id,cohort)]
      controldata[,obst:=sum(event_time==t),by=.(id,cohort)]
      
      pairdata<-rbind(treatdata[obsbase==1 & obst==1 & base_time != t & (event_time == t | event_time == base_time),],
                      controldata[obsbase==1 & obst==1 & base_time != t & (event_time == t | event_time == base_time),])
      pairdata[,time_pair:= t]
      eventdata<-rbind(eventdata,
                       pairdata)
    }
#  }
  #This stacking is unnecessary if we restrict to a balanced panel. Then no one enters
  #or exits the panel, and weights are fixed at the household level.
  # if(balanced_panel==TRUE){
  #   numperiods<-length(unique(treatdata$event_time))
  #   eventdata<-rbind(eventdata,
  #                    treatdata[obscount== numperiods,],
  #                    controldata[obscount>=numperiods,]
  #   )
  #   eventdata[,time_pair:=1]
  #   
  # }
  if(is.null(eventdata)==T){
    return(data.table())
  }else{
    
    eventdata[,obst:=NULL]
    eventdata[,obsbase:=NULL]
    eventdata[,anycohort:=NULL]
    
    rm(treatdata)
    rm(controldata)
    gc()
    
    eventdata[,post:=event_time >= 0]
    eventdata[,id:=as.factor(id)]
    eventdata[,treatgroup:=as.factor(treatgroup)]
    eventdata[,cohort:=as.factor(cohort)]  
    eventdata[,event_time:=as.factor(event_time)]
    eventdata[,time_pair:=as.factor(time_pair)]
    eventdata[,time:=as.factor(time)]
    
    basefactor<-unique(eventdata[event_time==base_time,event_time])
    basefactor<-basefactor[length(basefactor)]
    eventdata[,event_time:=relevel(event_time,ref=as.character(basefactor))]
    
    #calculating weights so that controls match treated households on characteristics
    #Note: this may have a lot of fixed effects, and may need to be broken down into 
    #multiple smaller regressions:

    eventdata[,pweight:=NULL]
    if(is.character(covariate_base_balance_linear)) eventdata[,pval:=
      feols(as.formula(paste0("treated ~", 
      paste0(covariate_base_balance_linear,collapse="+",sep=":interaction(cohort,event_time,time_pair,stratify,balancevars_linear_subset, drop = TRUE)"),
      "| interaction(cohort,event_time,time_pair,stratify,balancevars, drop = TRUE)"
      )),
                           data = eventdata, lean = FALSE)$fitted.values]
    else     eventdata[,pval:=feols(treated ~ 1 | interaction(cohort,event_time,time_pair,stratify,balancevars, drop = TRUE),
                                    data = eventdata, lean = FALSE)$fitted.values]

    
    eventdata[treated==1 & pval < 1 & pval > 0,pweight:=1]
    eventdata[treated==0 & pval < 1 & pval > 0,pweight:=pval/(1-pval)]
    eventdata[,pval:=NULL]
    eventdata<- eventdata[!is.na(pweight),]
    
    if(!is.na(instrument)){
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
      
    }

    #calculating weights to match control and treated households on characteristics, across strata
    if(!is.na(stratify_balance_val) & stratify_balance_val != "mean"){
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
    if(!is.na(stratify_balance_val) & stratify_balance_val == "mean"){
      
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
    
    eventdata[,treated:=as.factor(treated)]
    
    return(eventdata)
  }
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



#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================


event_ATTs_dynamic<-function(eventdata,
                             outcomes,#vector of variable names
                             clustervar="id", 
                             weights="pweight",
                             keep_trends=TRUE,
                             se="cluster",
                             lean=TRUE,
                             ssc=NULL){
  
  eventdata[,event_time_stratify:=as.factor(as.character(event_time_stratify))]
  eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  eventdata[,treated_event_time_stratify:=as.factor(as.character(treated_event_time_stratify))]
  eventdata[,treated_event_time_stratify:=relevel(treated_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  
  
  if(keep_trends == TRUE){
    results<-feols(as.formula(paste0("c(",
                                     paste0(outcomes,collapse=",") , ") ~ event_time_stratify + treated_event_time_stratify | unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = TRUE, mem.clean = TRUE)
  }
  else{
    results<-feols(as.formula(paste0("c(",
                                     paste0(outcomes,collapse=",") , ") ~ treated_event_time_stratify | event_time_stratify + unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = TRUE, mem.clean = TRUE)
    
  }
  
  return(list(dynamic = results))
}


event_ATTs_means<-function(eventdata,
                           outcomes,#vector of variable names
                           clustervar="id", 
                           weights="pweight",
                           keep_trends=TRUE,
                           se="cluster",
                           lean=TRUE,
                           ssc=NULL){
  
  eventdata$stratify <- factor(eventdata$stratify, levels = c(levels(eventdata$stratify), "empty"))
  
  results_means<-feols(as.formula(paste0("c(",
                                         paste0(outcomes,collapse=",") , 
                                         ") ~ interaction(stratify) - 1"
  )),
  data = eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),],
  weights= eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),get(weights)],
  cluster=clustervar, lean = TRUE, mem.clean = TRUE)
  
  return(list(means = results_means))
}


event_ATTs_pooled<-function(eventdata,
                            outcomes,#vector of variable names
                            clustervar="id", 
                            weights="pweight",
                            keep_trends=TRUE,
                            se="cluster",
                            lean=TRUE,
                            ssc=NULL){
  
  eventdata[,event_time_stratify:=as.factor(as.character(event_time_stratify))]
  eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  eventdata[,treated_event_time_stratify:=as.factor(as.character(treated_event_time_stratify))]
  eventdata[,treated_event_time_stratify:=relevel(treated_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  eventdata[,treated_pre_stratify:=as.factor(as.character(treated_pre_stratify))]
  eventdata[,treated_post_stratify:=as.factor(as.character(treated_post_stratify))]
  eventdata[,treated_post_stratify:=relevel(treated_post_stratify,ref = paste0(c(0,1),collapse="."))]
  
  if(keep_trends == TRUE){
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ event_time_stratify + treated_pre_stratify + treated_post_stratify | unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = TRUE, mem.clean = TRUE)
  }
  else{
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ treated_pre_stratify + treated_post_stratify | event_time_stratify + unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = TRUE, mem.clean = TRUE
    )
  }
  return(list(pooled = results_pooled))
}


event_ATTs_covariates<-function(eventdata,
                     outcomes=NULL,#vector of variable names
                     covariates=NULL, #vector of covariates on which to check balance
                     clustervar="id", 
                     weights="pweight",
                     keep_trends=TRUE,
                     se="cluster",
                     lean=TRUE,
                     ssc=NULL){

  #eventdata$stratify <- factor(eventdata$stratify, levels = c(levels(eventdata$stratify), "empty"))
  eventdata[,int_event_time_stratify:=interaction(event_time,stratify)]
  
  if(!is.null(covariates)) results_covariates<-feols(as.formula(paste0("c(",
                                                                       paste0(covariates,collapse=",") , ") ~ int_event_time_stratify + int_event_time_stratify:treated- 1"
  )),
  data = eventdata[as.numeric(as.character(event_time)) == as.numeric(as.character(time_pair)) | as.numeric(as.character(time_pair)) == 0,],
  weights= eventdata[as.numeric(as.character(event_time)) == as.numeric(as.character(time_pair)) | as.numeric(as.character(time_pair)) == 0,get(weights)],
  se=se,
  cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  else results_covariates=NULL
  
  return(list(covariates = results_covariates))
}


event_ATTs<-function(eventdata,
                     outcomes,#vector of variable names
                     covariates=NULL, #vector of covariates on which to check balance
                     clustervar="id", 
                     weights="pweight",
                     keep_trends=TRUE,
                     se="cluster",
                     lean=TRUE,
                     ssc=NULL){
  #These regressions should work identically if the fixed effects (after the "|") were replaced with:
  # interaction(time_pair,id,cohort)
  if(se != "cluster") clustervar<-NULL
  if(is.null(ssc)) ssc<-ssc()

  construct_event_variables(eventdata)
  eventdata$stratify <- factor(eventdata$stratify, levels = c(levels(eventdata$stratify), "empty"))
  
  if(keep_trends == TRUE){
    results<-feols(as.formula(paste0("c(",
                                     paste0(outcomes,collapse=",") , ") ~ event_time_stratify + treated_event_time_stratify | unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    se=se,
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
    
    
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ event_time_stratify + treated_pre_stratify + treated_post_stratify | unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    se=se,
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  }else{
    results<-feols(as.formula(paste0("c(",
                                     paste0(outcomes,collapse=",") , ") ~ treated_event_time_stratify | event_time_stratify + unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    se=se,
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
    
    
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ treated_pre_stratify + treated_post_stratify | event_time_stratify + unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    se=se,
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE
    )
  }
  
  results_means<-feols(as.formula(paste0("c(",
                                         paste0(outcomes,collapse=",") , 
                                         ") ~ interaction(stratify) - 1"
  )),
  data = eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),],
  weights= eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),get(weights)],
  se=se,
  cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  
  if(!is.null(covariates)) results_covariates<-feols(as.formula(paste0("c(",
                                   paste0(covariates,collapse=",") , ") ~ interaction(event_time,stratify) + interaction(event_time,stratify):treated- 1"
  )),
  data = eventdata[as.numeric(as.character(event_time)) == as.numeric(as.character(time_pair)) | as.numeric(as.character(time_pair)) == 0,],
  weights= eventdata[as.numeric(as.character(event_time)) == as.numeric(as.character(time_pair)) | as.numeric(as.character(time_pair)) == 0,get(weights)],
  se=se,
  cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  else results_covariates=NULL
  
  return(list(pooled = results_pooled,
              dynamic = results,
              means = results_means,
              covariates = results_covariates))
}


event_IVs<-function(eventdata,
                     outcomes,#vector of variable names
                     response,#binary, endogenous response variable
#                     instrument,#instrument we use to estimate causal effect of the response variable
                    saturate=FALSE,
#                    fullsaturate=FALSE,
                     clustervar="id", 
                     weights="pweight",
                     keep_trends=TRUE,
                     lean = TRUE,
                      ssc = NULL){
  if(is.null(ssc)) ssc<-ssc()

  construct_event_variables(eventdata,saturate=saturate,IV=TRUE,response=response)
  
  eventdata$stratify <- factor(eventdata$stratify, levels = c(levels(eventdata$stratify), "empty"))
  
  if(keep_trends == TRUE){
    results<-feols(as.formula(paste0("c(",
                                     paste0(outcomes,collapse=",") , ") ~ event_time_stratify + instrument_event_time_stratify + treated_event_time_stratify | unitfe |  treated_response_event_time_stratify  ~ treated_instrument_event_time_stratify"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
    
    
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ event_time_stratify + instrument_event_time_stratify + treated_pre_stratify + treated_post_stratify  | unitfe | treated_response_pre_stratify  + treated_response_post_stratify ~ treated_instrument_pre_stratify + treated_instrument_post_stratify"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  }else{
    results<-feols(as.formula(paste0("c(",
                                     paste0(outcomes,collapse=",") , ") ~ treated_event_time_stratify | event_time_stratify + instrument_event_time_stratify + unitfe | treated_response_event_time_stratify ~  treated_instrument_event_time_stratify "
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
    
    
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ treated_pre_stratify + treated_post_stratify | event_time_stratify + instrument_event_time_stratify + unitfe |  treated_response_pre_stratify +  treated_response_post_stratify ~ treated_instrument_pre_stratify + treated_instrument_post_stratify"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE
    )
  }
  
  
  results_means<-feols(as.formula(paste0("c(",
                                         paste0(outcomes,collapse=",") , 
                                         ") ~ interaction(stratify) - 1"
  )),
  data = eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),],
  weights= eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),get(weights)],
  cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  
  return(list(pooled = results_pooled,
              dynamic = results,
              means = results_means))
}

se_adjusted<-function(semats,eventdata,pw="pweight"){
  construct_event_variables(eventdata)
  baset<-max(eventdata$base_time)
  
  V<-model.matrix( ~event_time_stratify -1 
                   ,data =eventdata[event_time!=base_time])
  V<-as.matrix(V[,!grepl(paste0("event_time_stratify",baset),colnames(V))])
  V2<-V
  V2[eventdata[event_time!=base_time,treated==0],]<-0
  V<-cbind(V,V2)
  V<-V*eventdata[event_time!=base_time,get(pw)]
  
  
  D<-diag(1/eventdata[event_time!=base_time,2*get(pw)])
  
  correctingmat<-t(V)%*%D%*%V
  ZZ_inv<-solve(semats$ZZ-correctingmat)
  
  cov<-list()
  for(outcome in names(semats$sandwiches)){
    cov<-c(cov,list(ZZ_inv%*%semats$sandwich[[outcome]]%*%ZZ_inv))
  }
  names(cov)<-names(semats$sandwiches)
  return(cov)
}

se_iv_adjusted<-function(semats,eventdata,response,pw="pweight"){
  construct_event_variables(eventdata,saturate=FALSE,IV=TRUE,response=response)
  eventdata[,Ztreated:=Z*as.numeric(as.character(treated))]
  eventdata[,Rtreated:=R*as.numeric(as.character(treated))]
  baset<-max(eventdata$base_time)
  
  V <- model.matrix( ~event_time_stratify -1 
                      ,data =eventdata[event_time!=base_time])
  V<-as.matrix(V[,!grepl(paste0("event_time_stratify",baset),colnames(V))])
  V<-V*eventdata[event_time!=base_time,get(pw)]
  
  VZ <-V
  VX <-V
  for(vb in c("Z","treated","Ztreated")){
  V2<-V
  V2[eventdata[event_time!=base_time,get(vb)==0],]<-0
  VZ<-cbind(VZ,V2)
  }
  for(vb in c("Z","treated","Rtreated")){
  V2<-V
  V2[eventdata[event_time!=base_time,get(vb)==0],]<-0
  VX<-cbind(VX,V2)
  }

  D<-diag(1/eventdata[event_time!=base_time,2*get(pw)])
  
  correctingmatZ<-t(VZ)%*%D%*%VZ
  correctingmatX<-t(VZ)%*%D%*%VX
  ZZ<-semats$ZZ-correctingmatZ
  ZX<-semats$ZX-correctingmatX
  ZZ_inv <- solve(semats$ZZ-correctingmatZ)
  ZX_inv <- solve(semats$ZX-correctingmatX)
  
#  cov<-list()
#  for(outcome in names(semats$sandwiches)){
#  cov<-c(cov,list(ZX_inv%*%semats$sandwiches[[outcome]]%*%t(ZX_inv)))
#  }
#  names(cov)<-names(semats$sandwiches)
  cov<-ZX_inv%*%semats$sandwich%*%t(ZX_inv)
  return(cov)
#  ses<- solve(t(ZZ_inv%*%ZX)%*%ZX)%*%t(ZZ_inv%*%ZX)%*%semats$sandwich%*%ZZ_inv%*%ZX%*%solve(t(ZZ_inv%*%ZX)%*%ZX)
}

event_ATT_ses<-function(outcomes,
                       results,
                       eventdata,
                       includeFEs=FALSE,
                       clustervar="id", 
                       weights="pweight"){
  #Define useful variables:
  construct_event_variables(eventdata)
  
  
if(includeFEs==TRUE){
  Z<-model.matrix( ~event_time_stratify + treated_event_time_stratify + unitfe
                   ,data =eventdata)
  ndof<-0
}
  else{
    Z<-model.matrix( ~event_time_stratify + treated_event_time_stratify
                     ,data =eventdata)
    ndof<-nrow(eventdata)/2
  }
  Z<-Z[,which(apply(Z,MARGIN =2,FUN=var)!=0)]
  ZZ <- t(Z)%*%(Z*eventdata[,get(weights)])
  
  sandwiches<-list()
  for(outcome in outcomes){
  if(length(outcomes==1)) result <- results$dynamic
  if(length(outcomes> 1)) {
    outpos<-which(names(results$dynamic)==outcome)
    result <- results$dynamic[[outpos]]
  }
  eventdata[,e:=get(outcome) - predict(result,eventdata)]
  Ze <- Z*eventdata[,e*get(weights)]
  Ze<-as.data.table(Ze)
  Ze[,id:=eventdata[,get(clustervar)]]
  Ze<-as.matrix(Ze[,lapply(.SD,sum),by=id][,id:=NULL])
  
  
  sandwich<-t(Ze)%*%Ze
  sandwiches<-c(sandwiches,list(sandwich))
  }
  dof<-c(obs=nrow(eventdata),
         cl=length(unique(eventdata[,get(clustervar)]))
         )  
  #  dof<- (nrow(eventdata)-1)/(nrow(eventdata)-ncol(Ze)+1-ndof)
  #  dof.cl<-length(unique(eventdata[,get(clustervar)]))/(length(unique(eventdata[,get(clustervar)]))-1)
  names(sandwiches)<-outcomes
  
  return(list(ZZ=ZZ,
              sandwiches=sandwiches,
              dof=dof)
  )
}

event_IV_ses<-function(outcome,
                       response,
                 results,
                 eventdata,
                 includeFEs=FALSE,
                 clustervar="id", 
                 weights="pweight"){
  
  #Define useful variables:
  construct_event_variables(eventdata,saturate=FALSE,IV=TRUE,response=response)
  
  
if(includeFEs==TRUE) {
  Z<-model.matrix( ~event_time_stratify + instrument_event_time_stratify + treated_event_time_stratify  + treated_instrument_event_time_stratify + unitfe
              ,data =eventdata)
  
  X<-model.matrix( ~event_time_stratify + instrument_event_time_stratify + treated_event_time_stratify  + treated_response_event_time_stratify + unitfe
                   ,data =eventdata)
  ndof<-0
}
  else{
    Z<-model.matrix( ~event_time_stratify + instrument_event_time_stratify + treated_event_time_stratify  + treated_instrument_event_time_stratify
                     ,data =eventdata)
    X<-model.matrix( ~event_time_stratify + instrument_event_time_stratify + treated_event_time_stratify  + treated_response_event_time_stratify 
                     ,data =eventdata)
    
    ndof<-nrow(eventdata)/2
  }
Z<-Z[,which(apply(Z,MARGIN =2,FUN=var)!=0)]
X<-X[,which(apply(X,MARGIN =2,FUN=var)!=0)]
ZZ <- t(Z)%*%(Z*eventdata[,get(weights)])
ZX <- t(Z)%*%(X*eventdata[,get(weights)])

varnames<-rownames(results)
#correcting name on the IV coefficient(s):
varnames<-gsub("fit_","",varnames)
if(!all(colnames(X)%in%varnames)) stop("The variables we need to predict a person's outcome (columns of matrix X) do not all appear in the given coeftable 'results'")

  #sandwiches<-list()
  #eventdata[,e:=get(outcome) - predict(result,eventdata)]
 e<-eventdata[,get(outcome)-mean(get(outcome)),by=unitfe]$V1
 Xres<-as.data.table(X)
 Xres[,unitfe:=eventdata[,unitfe]]
 Xres<-Xres[,.SD-lapply(.SD,mean),by=unitfe]
 Xres[,unitfe:=NULL]
 Xres<-as.matrix(Xres)
 e<-e - apply(Xres[,varnames], MARGIN=1,FUN=function(x) sum(x*results[,"Estimate"]))
 Xres<-NULL
 Ze <- Z*e*eventdata[,get(weights)]
 Ze<-as.data.table(Ze)
 Ze[,id:=eventdata[,get(clustervar)]]
 Ze<-as.matrix(Ze[,lapply(.SD,sum),by=id][,id:=NULL])
 
 
 sandwich<-t(Ze)%*%Ze
# sandwiches<-c(sandwiches,list(sandwich))

 dof<-c(obs=nrow(eventdata),
        cl=length(unique(eventdata[,get(clustervar)]))
 )  
 #  dof<- (nrow(eventdata)-1)/(nrow(eventdata)-ncol(Ze)+1-ndof)
 #  dof.cl<-length(unique(eventdata[,get(clustervar)]))/(length(unique(eventdata[,get(clustervar)]))-1)
 #names(sandwiches)<-outcome

  return(list(ZZ=ZZ,
              ZX=ZX,
              sandwich=sandwich,
              dof=dof)
  )
}

#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================

get_result_dynamic<-function(eventdata_panel,variable,trends=TRUE){
  
  if(nrow(eventdata_panel)==0){
    dt<-data.table()
    return(dt)
  }else if ((eventdata_panel[,var(get(variable))]==0)==T){
    dt<-data.table()
    return(dt)
  }else{
    
    results<-event_ATTs_dynamic(eventdata_panel,outcomes = c(variable),keep_trends = trends)
    dt<-data.table(variable = row.names(results$dynamic$coeftable),model=i,results$dynamic$coeftable,obs=results$dynamic$nobs)
    dt<-dt[,result:="dynamic"]
    rm(eventdata_panel)
    rm(results)
    gc()
    return(dt)
  }
}

get_result_pooled<-function(eventdata_panel,variable,trends=TRUE){
  
  if(nrow(eventdata_panel)==0){
    dt<-data.table()
    return(dt)
  }else if ((eventdata_panel[,var(get(variable))]==0)==T){
    dt<-data.table()
    return(dt)
  }else{
    
    results<-event_ATTs_pooled(eventdata_panel,outcomes = c(variable),keep_trends = trends)
    dt<-data.table(variable = row.names(results$pooled$coeftable),model=i,results$pooled$coeftable,obs=results$pooled$nobs)
    dt<-dt[,result:="pooled"]
    rm(eventdata_panel)
    rm(results)
    gc()
    return(dt)
    
  }
}

get_result_means<-function(eventdata_panel,variable,trends=TRUE){
  
  results<-event_ATTs_means(eventdata_panel,outcomes = c(variable),keep_trends = trends)
  dt<-data.table(variable = row.names(results$means$coeftable),model=i,results$means$coeftable,obs=results$means$nobs)
  dt<-dt[,result:="means"]
  
  return(dt)
  
}

get_result_covariates<-function(eventdata_panel,covariate,variable=NULL,trends=TRUE){
  
  results<-event_ATTs_covariates(eventdata_panel,covariates=covariate,outcomes = c(variable),keep_trends = trends)
  dt<-data.table(variable = row.names(results$covariates$coeftable),results$covariates$coeftable,obs=results$covariates$nobs)
  dt<-dt[,result:="covariates"][,model:=paste0("additional_covariate--",covariate)]
  
  return(dt)
  
}



#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================


print_ATT<-function(results,
                    outcomes=NULL,
                    outcome_names=NULL,
                    event_name="Event",
                    pooled_tables=TRUE,
                    dynamic_plots=TRUE,
                    pooled_tables_name=NULL,
                    base_time = NULL,
                    dynamic_pdfname=NULL,
                    stratify_values = NULL,
                    stratify_names = NULL,
                    decimals=0,
                    plot_pval=0.05){
  if(!is.null(stratify_values) & is.null(stratify_names)) stop("Must provide a name for each entry of stratify_values in stratify_names")
  if(pooled_tables==TRUE){
    
    if(is.null(outcome_names)) outcome_names <- names(results$pooled)
    if(is.null(outcomes)) outcomes <- names(results$pooled)
    if(!is.null(stratify_values)){
      tab<-TR(c("","Pre-period Mean","Treatment Effect"),cspan=c(1,length(stratify_values),length(stratify_values)))
      tab<-tab + TR(c("Outcome",stratify_names,stratify_names))
    }
    else{
      tab<-TR(c("Outcome","Pre-period Mean","Treatment Effect"))
      stratify_values<-1
    }
    tab<-tab + midrulep(list(c(2,length(stratify_values)+1),c(length(stratify_values)+2,2*length(stratify_values)+1))) 
    p<-1
    for(outcome in outcomes){
      pos<-which(names(results$pooled)==outcome)
      tab<-tab + TR(outcome_names[p]) %:%
        #Adding in the pre-period means:
        TR(results$means[[pos]]$coefficients[paste0("interaction(stratify)",stratify_values)],
           pvalues = ifelse(is.na(results$means[[pos]]$coeftable[paste0("interaction(stratify)",stratify_values),"Pr(>|t|))"]),
                            1,
                            results$means[[pos]]$coeftable[paste0("interaction(stratify)",stratify_values),"Pr(>|t|))"]
           ), dec = decimals) %:%
        #Adding in the treatment effects:
        TR(results$pooled[[pos]]$coefficients[paste0("treated_post_stratify1.",stratify_values)],
           pvalues = ifelse(is.na(results$pooled[[pos]]$coeftable[paste0("treated_post_stratify1.",stratify_values),"Pr(>|t|))"]),
                            1,
                            results$pooled[[pos]]$coeftable[paste0("treated_post_stratify1.",stratify_values),"Pr(>|t|))"]
           ), dec = decimals)
      #STANDARD ERRORS:
      tab<- tab + TR("") %:%
        #For pre-period means
        TR(results$means[[pos]]$se[paste0("interaction(stratify)",stratify_values)],se=TRUE,dec = decimals)%:%
        #For treatment effects:
        TR(results$pooled[[pos]]$se[paste0("treated_post_stratify1.",stratify_values)],se=TRUE,dec = decimals)
      
      p<-p+1
    }
    TS(tab, file=paste0(paste(c(pooled_tables_name,"pooled_table"),collapse="_")),
       output_path=".",
       pretty_rules=T,
       header=c('r',rep('c',2*length(stratify_values))))
    
  }
  
  if(dynamic_plots==TRUE){
    p<-1
    for(outcome in outcomes){
      pos<-which(names(results$pooled)==outcome)
      
      plot_table<-data.table(
        coef = results$dynamic[[pos]]$coefficients[
          grepl("treated_event_time_stratify",names(results$dynamic[[pos]]$coefficients))
        ],
        se = results$dynamic[[pos]]$se[
          grepl("treated_event_time_stratify",names(results$dynamic[[pos]]$se))
        ],
        event_time = as.numeric(as.character(gsub("\\..*","", 
                                                  gsub("treated_event_time_stratify", "",
                                                       names(results$dynamic[[pos]]$coefficients)[grepl("treated_event_time_stratify",names(results$dynamic[[pos]]$se))])
        ))),
        stratify_value = as.numeric(as.character(gsub(".*\\.","", 
                                                      gsub("treated_event_time_stratify", "",
                                                           names(results$dynamic[[pos]]$coefficients)[grepl("treated_event_time_stratify",names(results$dynamic[[pos]]$se))])
        )))
      )
      #adding in omitted year:
      #Looking for the first missing period in the interval of plotted periods
      #If the interval of plotted periods has no gap, I assume refernce period is first
      #period in the data.
      if(is.null(base_time)){ 
        missingperiods<-(min(plot_table$event_time):max(plot_table$event_time))[!(min(plot_table$event_time):max(plot_table$event_time))%in%unique(plot_table$event_time)]
        if(length(missingperiods)==0) refperiod <- min(plot_table$event_time)-1
        else refperiod <- min(missingperiods)
      }
      else refperiod <- base_time
      plot_table<-rbind(plot_table, data.table(coef=0,
                                               se = 0,
                                               event_time = refperiod,
                                               stratify_value = unique(plot_table$stratify_value)
      ))
      
      plot_table[,upper:= coef + abs(qt(plot_pval/2,
                                        df = results$dynamic[[pos]]$nobs - results$dynamic[[pos]]$nparams))*se]
      plot_table[,lower:= coef - abs(qt(plot_pval/2,
                                        df = results$dynamic[[pos]]$nobs - results$dynamic[[pos]]$nparams))*se]
      pdf(paste0("./",paste(c(dynamic_pdfname,outcome),collapse="_"),".pdf"))
      for(stratval in stratify_values){
        print(ggplot(data=plot_table[stratify_value==stratval,]) + 
                geom_ribbon(aes(ymin = lower, ymax=upper, x=event_time), fill="grey50", alpha=0.5) + 
                geom_hline(yintercept = 0, linetype="dashed") + 
                geom_vline(xintercept = 0, linetype="dashed") + 
                geom_line(aes(x = event_time,
                              y = coef)) +
                scale_x_continuous( breaks = pretty_breaks(12)) +
                ylim(c(min(plot_table$lower),max(plot_table$upper))) +
                labs(y=varnames[p], x = paste0("Years from ",event_name))
        )
      }
      dev.off()
      p<-p+1
    }
  }
  
}

#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================

counter_eventdata <- function(eventdata){
  eventdata_counter<-eventdata[treated==1 & pweight ==1 ,.N, by=.(stratify,event_time)]
  return(eventdata_counter)
}


summary_eventdata1 <- function(eventdata,
                               summarylist,
                               summarylevel){
  
  eventdata<-eventdata[,stratify2:=paste(get(summarylevel),post,treatgroup,sep=".")]
  weights<-"pweight"
  summarytable1<-data.table()
  for(h in summarylist){
    
    results_means<-feols(as.formula(paste0(h, "~ interaction(stratify2) - 1")),
                         data = eventdata, cluster = "id", weights = eventdata[,get(weights)],
                         lean = TRUE, mem.clean = TRUE)
    tab<-data.table(variable = row.names(results_means$coeftable),results_means$coeftable,outcome=h,obs=results_means$nobs)
    summarytable1<-rbind(summarytable1,tab)
    console.log(h)
  }
  tab_N<-eventdata[,.N,by=.(stratify2)][,variable:=paste0("interaction(stratify2)",stratify2)]
  summarytable1<-tab_N[summarytable1,on="variable"]
  
  return(summarytable1)
  rm(results_means)
  rm(summarytable1)
  gc()
  
}


summary_eventdata2 <- function(eventdata,
                               summarylist,
                               summarylevel){
  
  eventdata<-eventdata[,stratify2:=paste(get(summarylevel),post,treatgroup,sep=".")]
  weights<-"pweight"
  summarytable2<-data.table()
  for(h in summarylist){
    tab1<-eventdata[,.(weighted.mean(get(h),pweight),wtd.var(get(h),pweight),mean(get(h)),var(get(h))),by=.(stratify2)]%>%rename(weighted_mean=V1,
                                                                                                                                 weighted_var=V2,
                                                                                                                                 mean=V3,
                                                                                                                                 var=V4)
    tab1<-tab1[,outcome:=h]
    summarytable2<-rbind(summarytable2,tab1,fill=T)
    console.log(h)
    gc()
    
  }
  tab_N<-eventdata[,.N,by=.(stratify2)][,variable:=paste0("interaction(stratify2)",stratify2)]
  summarytable2<-tab_N[summarytable2,on="stratify2"]
  return(summarytable2) 
  rm(list=ls(pattern="^tab"))
}


#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================
#=======================================================================================================================================================

# How to deal with memory issue ?

# dynamic : 

#eventdata[,.N,by=.(time_pair)] -> find out lower bound and upper bound
#dt_dynamic<-data.table()
#for(i in -13:-2){
#  gc()
#  dt_dynamic<-get_result(eventdata,i,i,"outcome",dt_dynamic,trends=T)
#  gc()
#  print(i)
#}
#for(i in 0:12){
#  gc()
#  dt_dynamic<-get_result(eventdata,i,i,"outcome",dt_dynamic,trends=T)
#  gc()
#  print(i)
#}
#fwrite(dt_dynamic,"")

#get_result_dynamic_old<-function(eventdata_panel,start,end,variable,table,results=list(),pos=1,trends=TRUE){
  
#  b = end-start+1
#  for(eventtime in start:end){
#    results[[pos]] <- event_ATTs_dynamic(eventdata_panel[time_pair==eventtime,],outcomes=c(variable),keep_trends=trends)
#    pos<-pos+1
#  }
#  for(i in 1:b){
#    dt<-data.table(variable = row.names(results[[i]]$dynamic$coeftable),model=i,results[[i]]$dynamic$coeftable,obs=results[[i]]$dynamic$nobs)
#    table<-rbind(dt,table)
#  }
#  return(table)
#}


#event_ATTs_head<-function(eventdata,
#                          outcomes,#vector of variable names
#                          clustervar="id", 
#                          weights="pweight",
#                          keep_trends=TRUE){
  #These regressions should work identically if the fixed effects (after the "|") were replaced with:
  # interaction(time_pair,id,cohort)
#  eventdata[,treated_post := as.factor((treated == 1) * (post == 1))]
#  eventdata[,treated_pre := as.factor((treated == 1) * (post == 0))]
#  eventdata[,treated_event_time := event_time]
#  eventdata[treated==0,treated_event_time := 1] #1 is the base level
  
#  eventdata[,event_time_stratify:=interaction(event_time,stratify, drop = TRUE)]
#  eventdata[,treated_post_stratify := interaction(treated_post,stratify, drop = TRUE)]
  
#  eventdata[,treated_pre_stratify := interaction(treated_pre,stratify,drop=TRUE)]
#  eventdata[treated_pre==0,treated_pre_stratify := 1]
#  eventdata[event_time==base_time,treated_pre_stratify := 1]
  
#  eventdata[,unitfe := interaction(time_pair,treated,stratify, drop = TRUE)]
  
#  eventdata[,treated_event_time_stratify := interaction(event_time,stratify, drop = TRUE)]
  
  
  #Omitting base year for all levels of --stratify--:
#  eventdata[event_time==base_time,event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
#  eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  
  #Omitting base year for all levels of --stratify--, for treated people
#  eventdata[treated == 0 ,treated_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
#  eventdata[event_time==base_time,treated_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
#  eventdata[,treated_event_time_stratify:=relevel(treated_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  
  
  #Omitting effect for untreated people or observations in pre-period:
#  eventdata[treated_post == 0 ,treated_post_stratify := paste0(c(0,1),collapse=".")]
#  eventdata[,treated_post_stratify:=relevel(treated_post_stratify,ref = paste0(c(0,1),collapse="."))]
  
#  return(eventdata)
  
#}

#======================================================

console.log <- function ( ... ) { cat(format(Sys.time(),"(%Y/%b/%d) %X"),...,"\n") }