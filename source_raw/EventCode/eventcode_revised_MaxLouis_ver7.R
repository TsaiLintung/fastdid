# eventcode


#Treated households:
create_event_data<-function(maindata,
                            #if no covariates for balancing or stratification, just specify a constant (as follows):
                            covariate_base_stratify=1, #vector of variable names, treatment effects are stratified by their joint unique values.
                            covariate_base_balance=1, #vector of variable names to include in finding exact matches for treated units
                            base_restrict = 1, #a single var, restricting to households for which var == 1 in base period,
                            treated_restrict=1,
                            timevar,
                            unitvar,
                            cohortvar, #period when the particular event of interest arises
                            anycohortvar = NULL, #period when ANY kind of event arises, used to find control cohorts
                            #(for whom no event of any kind has yet to happen)
                            #Only specify if it differs from cohortvar
                            onset_agevar, #variable indicating age at which LTC needs arose, for person affected by them
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
                            stratify_balance_val = NA) 
{
  
  if(lower_event_time > base_time) stop("lower_event_time must lie below base_time")
  if(!is.data.table(maindata)) stop("rawdata must be a data.table")
  #maindata<-copy(rawdata)
  if(is.null(anycohortvar)) {
    maindata[,anycohortvar:=get(cohortvar)]
    anycohortvar<-"anycohortvar"
  }
  if(is.numeric(base_time)) maindata[,base_time:=base_time]
  else setnames(maindata,base_time,"base_time")
  setnames(maindata,c(timevar,unitvar,cohortvar,anycohortvar),c("time","id","cohort","anycohort"))
  if(base_restrict != 1) setnames(maindata,base_restrict,"base_restrict")
  if(treated_restrict != 1) setnames(maindata,treated_restrict,"treated_restrict")
  
  if(any(is.na(maindata$cohort))) stop("cohort variable should not be missing (it can be infinite instead)")
  if(any(is.na(maindata$anycohort))) stop("anycohort variable should not be missing (it can be infinite instead)")
  
  treatdata<-copy(maindata[onset_agevar>=onset_minimum  & ! is.infinite(cohort) & !is.infinite(anycohort) ,])
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
    
    controlcohort[,cohort := o]
    controlcohort[,event_time := time - cohort]
    
    #Make sure people in the control cohort are actually observed in that period
    #(to verify they don't belong to the cohort)
    controlcohort[ ,obscohort := max(time == o),by=id]
    controlcohort<-controlcohort[obscohort==1,]
    controlcohort[,obscohort:=NULL]
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
  controldata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
  controldata <- controldata[base_restrict == 1,]
  treatdata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
  treatdata <- treatdata[base_restrict == 1,]
  treatdata <- treatdata[treated_restrict == 1,]
  
  if(is.character(covariate_base_stratify)) {
    treatdata[,stratify:=interaction(treatdata[,covariate_base_stratify,with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
    controldata[,stratify:=interaction(controldata[,covariate_base_stratify,with=FALSE], drop=TRUE)]#, mem.clean=TRUE)]
  }else {
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
  treatdata[,temp:=interaction(balancevars,stratify,event_time,cohort,drop=TRUE)]
  controldata[,temp:=interaction(balancevars,stratify,event_time,cohort,drop=TRUE)]
  
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
  if(balanced_panel==FALSE){
    for(t in event_times){
      treatdata[,obst:=sum(event_time==t),by=.(id,cohort)]
      controldata[,obst:=sum(event_time==t),by=.(id,cohort)]
      
      pairdata<-rbind(treatdata[obsbase==1 & obst==1 & base_time != t & (event_time == t | event_time == base_time),],
                      controldata[obsbase==1 & obst==1 & base_time != t & (event_time == t | event_time == base_time),])
      pairdata[,time_pair:= t]
      eventdata<-rbind(eventdata, 
                       pairdata)
    }
  }

  
  eventdata[,obst:=NULL]
  eventdata[,obsbase:=NULL]
  eventdata[,anycohort:=NULL]

  #This stacking is unnecessary if we restrict to a balanced panel. Then no one enters
  #or exits the panel, and weights are fixed at the household level.
  if(balanced_panel==TRUE){
    eventdata<-rbind(eventdata,
                     treatdata[obscount== max(treatdata$obscount),],
                     controldata[obscount>=max(treatdata$obscount),]
    )
    
  }
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
  eventdata[,pweight:=NA]
  eventdata[,pval:=feols(treated ~ 1 | interaction(cohort,event_time,time_pair,stratify,balancevars, drop = TRUE),
                         data = eventdata, lean = FALSE)$fitted.values]
  
  eventdata[treated==1 & pval < 1 & pval > 0,pweight:=1.00]
  eventdata[, pweight := as.double(pweight)] #TSAI addition
  eventdata[treated==0 & pval < 1 & pval > 0,pweight:=pval/(1-pval)]
  eventdata[,pval:=NULL]
  eventdata<- eventdata[!is.na(pweight),]
  
  #calculating weights to match control and treated households on characteristics, across strata
  if(!is.na(stratify_balance_val)){
    eventdata[,treated_base := treated == 1 & stratify == stratify_balance_val]
    
    stratvals<-levels(eventdata$stratify)
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
  
  eventdata[,treated:=as.factor(treated)]
  
  return(eventdata)
}

event_ATTs_head<-function(eventdata,
                          outcomes,#vector of variable names
                          clustervar="id", 
                          weights="pweight",
                          keep_trends=TRUE){
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
  
  eventdata[,unitfe := interaction(time_pair,id,treated,cohort,stratify, drop = TRUE)]
  
  
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
  
  return(eventdata)
  
}


event_ATTs_dynamic<-function(eventdata,
                             outcomes,#vector of variable names
                             clustervar="id", 
                             weights="pweight",
                             keep_trends=TRUE){
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
                           keep_trends=TRUE){
  
  results_means<-feols(as.formula(paste0("c(",
                                         paste0(outcomes,collapse=",") , 
                                         ") ~ interaction(stratify) - 1"
  )),
  data = eventdata[treated == 1 & event_time == base_time & get(weights) == 1,],
  cluster=clustervar, lean = TRUE, mem.clean = TRUE)
  
  return(list(means = results_means))
}

event_ATTs_pooled<-function(eventdata,
                            outcomes,#vector of variable names
                            clustervar="id", 
                            weights="pweight",
                            keep_trends=TRUE){
  
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

event_ATTs<-function(eventdata,
                     outcomes,#vector of variable names
                     clustervar="id", 
                     weights="pweight",
                     keep_trends=TRUE){
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
  
  eventdata[,unitfe := interaction(time_pair,id,treated,cohort,stratify, drop = TRUE)]
  
  
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
  
  if(keep_trends == TRUE){
    results<-feols(as.formula(paste0("c(",
                                     paste0(outcomes,collapse=",") , ") ~ event_time_stratify + treated_event_time_stratify | unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = TRUE, mem.clean = TRUE)
    
    
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ event_time_stratify + treated_pre_stratify + treated_post_stratify | unitfe"
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
    
    
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ treated_pre_stratify + treated_post_stratify | event_time_stratify + unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = TRUE, mem.clean = TRUE
    )
  }
  
  results_means<-feols(as.formula(paste0("c(",
                                         paste0(outcomes,collapse=",") , 
                                         ") ~ interaction(stratify) - 1"
  )),
  data = eventdata[treated == 1 & event_time == base_time & get(weights) == 1,],
  cluster=clustervar, lean = TRUE, mem.clean = TRUE)
  
  return(list(pooled = results_pooled,
              dynamic = results,
              means = results_means))
}

event_ATTs_pooledmeans<-function(eventdata,
                                 outcomes,#vector of variable names
                                 clustervar="id", 
                                 weights="pweight",
                                 keep_trends=TRUE){
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
    cluster=clustervar, lean = TRUE, mem.clean = TRUE)
  }
  
  results_means<-feols(as.formula(paste0("c(",
                                         paste0(outcomes,collapse=",") , 
                                         ") ~ interaction(stratify) - 1"
  )),
  data = eventdata[treated == 1 & event_time == base_time & get(weights) == 1,],
  cluster=clustervar, lean = TRUE, mem.clean = TRUE)

  return(list(pooled = results_pooled,means = results_means))
}


get_result_dynamic<-function(eventdata_panel,start,end,variable,table,results=list(),pos=1,trends=TRUE){
  
  b = end-start+1
  for(eventtime in start:end){
    results[[pos]] <- event_ATTs_dynamic(eventdata_panel[time_pair==eventtime,],outcomes=c(variable),keep_trends=trends)
    pos<-pos+1
  }
  for(i in 1:b){
    dt<-data.table(variable = row.names(results[[i]]$dynamic$coeftable),model=i,results[[i]]$dynamic$coeftable,obs=results[[i]]$dynamic$nobs)
    table<-rbind(dt,table)
  }
  return(table)
}

get_result_pooledmeans<-function(eventdata_panel,variable,trends=TRUE){
  
  results<-event_ATTs_pooledmeans(eventdata_panel,outcomes = c(variable),keep_trends = trends)
  dt1<-data.table(variable = row.names(results$pooled$coeftable),model=i,results$pooled$coeftable,obs=results$pooled$nobs)
  dt2<-data.table(variable = row.names(results$means$coeftable),model=i,results$means$coeftable,obs=results$means$nobs)
  dt<-rbind(dt1[,result:="pooled"],dt2[,result:="means"],fill=T)
  
  return(dt)
  
}

get_result_pooled<-function(eventdata_panel,variable,trends=TRUE){
  
  results<-event_ATTs_pooled(eventdata_panel,outcomes = c(variable),keep_trends = trends)
  dt1<-data.table(variable = row.names(results$pooled$coeftable),model=i,results$pooled$coeftable,obs=results$pooled$nobs)
  #dt2<-data.table(variable = row.names(results$means$coeftable),model=i,results$means$coeftable,obs=results$means$nobs)
  dt<-rbind(dt1[,result:="pooled"],fill=T)
  
  return(dt)
  
}

get_result_means<-function(eventdata_panel,variable,trends=TRUE){
  
  results<-event_ATTs_means(eventdata_panel,outcomes = c(variable),keep_trends = trends)
  #dt1<-data.table(variable = row.names(results$pooled$coeftable),model=i,results$pooled$coeftable,obs=results$pooled$nobs)
  dt2<-data.table(variable = row.names(results$means$coeftable),model=i,results$means$coeftable,obs=results$means$nobs)
  dt<-rbind(dt2[,result:="means"],fill=T)
  
  return(dt)
  
}
