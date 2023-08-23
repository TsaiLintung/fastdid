#by LinTung Tsai to be used with Max's event code
#some wrappers
library(fixest)

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
                            stratify_balance_val = NA,
                            parallel = FALSE) 
{
  
  #deal with with name
  if(is.null(anycohortvar)) {
    maindata[,anycohortvar:=get(cohortvar)]
    anycohortvar<-"anycohortvar"
  }
  setnames(maindata,c(timevar,unitvar,cohortvar,anycohortvar),c("time","id","cohort","anycohort"))
  if(base_restrict != 1) setnames(maindata,base_restrict,"base_restrict")
  if(treated_restrict != 1) setnames(maindata,treated_restrict,"treated_restrict")
  
  #data validation
  if(lower_event_time > base_time) stop("lower_event_time must lie below base_time")
  if(!is.data.table(maindata)) stop("rawdata must be a data.table")
  if(any(is.na(maindata$cohort))) stop("cohort variable should not be missing (it can be infinite instead)")
  if(any(is.na(maindata$anycohort))) stop("anycohort variable should not be missing (it can be infinite instead)")
  if(any(fduplicated(maindata[,.(time, id)]))){stop ("some unit-time is observed more then once")}
  
  maindata[, id := charToFact(as.character(id))]
  
  #turn the covariates into a single interacted factor
  if(is.character(covariate_base_stratify)) {
    maindata[, stratify:= charToFact(as.character(.GRP)), by =covariate_base_stratify]
  }else {
    maindata[,stratify:="base_strat"]
  }
  if(is.character(covariate_base_balance)) {
    maindata[, balancevars:= charToFact(as.character(.GRP)), by =covariate_base_balance]
  }else {
    maindata[,balancevars:="base_balance"]
  }
  
  #create treated data
  treatdata<-copy(maindata[onset_agevar>=onset_minimum  & ! is.infinite(cohort) & !is.infinite(anycohort) ,])
  treatdata[,treated:=1]
  treatdata[,event_time:=time-cohort]
  treatdata<-treatdata[event_time >= lower_event_time & event_time <= upper_event_time,]
  treatdata[,treatgroup:="treated"]
  
  #stack control cohorts for each treated cohort
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
  
  

  #If base_time varies across units, reassigning it to a common reference value:
  #This is relevant, for instance, with a dataset that moves from annual to bi-annual
  #treatdata[event_time == base_time, event_time :=max(base_time)]
  #treatdata[,base_time :=max(base_time)]
  #controldata[event_time == base_time,event_time :=max(base_time)]
  #controldata[,base_time :=max(base_time)]
  
  #count how many times a id is viewed
  treatdata[,obscount:=1]
  controldata[,obscount:=1]
  treatdata[,obscount:=sum(obscount),by=.(id,cohort)]
  controldata[,obscount:=sum(obscount),by=.(id,cohort)]
  
  #if a covariate is > 9e9 and , it is set to NA
  if(is.character(covariate_base_stratify)){
    for(out in covariate_base_stratify){
      controldata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
      controldata[get(out) >= 9e9,eval(out) := NA, ]
      treatdata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
      treatdata[get(out) >= 9e9,eval(out) := NA, ]
    }
  }
  if(is.character(covariate_base_balance)){
    for(out in covariate_base_balance){
      controldata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
      controldata[get(out) >= 9e9,eval(out) := Inf, ]
      treatdata[,eval(out) := min(get(out) + 9e9 *(event_time != base_time)), by=.(id,cohort)]
      treatdata[get(out) >= 9e9,eval(out) := Inf, ]
    }
  }
  
  
  controldata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
  controldata <- controldata[base_restrict == 1,]
  treatdata[,base_restrict := max(base_restrict * (event_time == base_time), na.rm=TRUE),by=.(id, cohort)]
  treatdata <- treatdata[base_restrict == 1,]
  treatdata <- treatdata[treated_restrict == 1,]
  
  #only keep observations that have common support on covariates
  #treat_unique <- funique(treatdata[, .(balancevars,stratify,event_time, cohort)]) |> as.data.table()
  #control_unique <-funique(controldata[, .(balancevars,stratify,event_time, cohort)]) |> as.data.table()
  #common_unique <- fintersect(treat_unique, control_unique)
  
  #treatdata[,temp:=interaction(balancevars,stratify,event_time,cohort,drop=TRUE)]
  #controldata[,temp:=interaction(balancevars,stratify,event_time,cohort,drop=TRUE)]
  #commonvals<-intersect(unique(treatdata$temp),unique(controldata$temp))
  #treatdata<-treatdata[temp%in%commonvals,]
  #controldata<-controldata[temp%in%commonvals,]
  
  #treatdata[,temp:=NULL]
  #controldata[,temp:=NULL]
  #rm(commonvals)
  gc()
  
  #before stacking turn cohort and event time into fact
  treatdata[, event_time_fact := charToFact(as.character(event_time))]
  treatdata[, cohort_fact := charToFact(as.character(cohort))]
  controldata[, event_time_fact := charToFact(as.character(event_time))]
  controldata[, cohort_fact := charToFact(as.character(cohort))]
  
  #make sure observations are only observed once in the base period
  #treatdata[,obsbase:=sum(event_time==base_time),by=.(id,cohort)]
  #controldata[,obsbase:=sum(event_time==base_time),by=.(id,cohort)]
  #if(max(treatdata$obsbase)>1) stop("Error: some treated units are observed more than once in the reference period")
  #if(max(controldata$obsbase)>1) stop("Error: some control units are observed more than once in the reference period")
  #treatdata <- treatdata[obsbase==1,]
  #controldata <- controldata[obsbase==1,]
  
  #For the final dataset, we must stack each pairwise combo of (base year, other year)
  #for each household. The reason? We will assign control households weights that vary
  #based on the other year, as households enter/exit the sample.
  event_times<- as.integer(as.character(treatdata[,unique(event_time)]))
  event_times <- event_times[event_times %!=% -1]
  
  stack_eventtime <- function(t, treatdata, controldata){
    treatdata[,obst:=anyv(event_time,t),by=.(id,cohort)]
    controldata[,obst:=anyv(event_time,t),by=.(id,cohort)]
    
    pair_treat_data <- treatdata[obst == TRUE & event_time %in% c(t, base_time),]
    pair_control_data <- controldata[obst == TRUE & event_time %in% c(t, base_time),]
    
    fact_t <- charToFact(as.character(t))
    
    pair_treat_data[,time_pair := fact_t]
    pair_control_data[,time_pair := fact_t]
    
    #return this
    return(c(list(pair_treat_data), list(pair_control_data)))
  }
  
  data_list <- list()
  
  if(parallel == TRUE){
    options(future.globals.maxSize= 1000*1024*1024)
    doFuture::registerDoFuture()
    plan("multisession")
    data_list <- foreach(t = event_times) %dopar% stack_eventtime(t, treatdata, controldata)
  } else {
    data_list <- foreach(t = event_times) %do% stack_eventtime(t, treatdata, controldata)
  }


  data_list <- flatten(data_list)
  eventdata <- rbindlist(data_list, use.names=TRUE)

  eventdata[,obst:=NULL]
  eventdata[,anycohort:=NULL]
  
  rm(treatdata)
  rm(controldata)
  gc()
  
  eventdata[,post:= event_time >= 0]
  
  #calculate ipw
  #Note: this may have a lot of fixed effects, and may need to be broken down into multiple smaller regressions:
  eventdata[,pweight:=NA]
  
  eventdata[,pval:= feols(treated ~ 1 | cohort_fact^event_time_fact^time_pair^stratify^balancevars,
                          data = eventdata, lean = FALSE, combine.quick = TRUE)$fitted.values]
  
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
  
  return(eventdata)
}

event_ATTs_head<-function(eventdata,
                          outcomes,#vector of variable names
                          clustervar="id", 
                          weights="pweight",
                          keep_trends=TRUE,
                          base_time = -1){
  
  #These regressions should work identically if the fixed effects (after the "|") were replaced with:
  # interaction(time_pair,id,cohort)
  #eventdata[,treated_post := charToFact(as.character((treated == 1) * (post == 1)))]
  #eventdata[,treated_pre := charToFact(as.character((treated == 1) * (post == 0)))]
  #eventdata[,treated_event_time := event_time]
  #eventdata[treated==0,treated_event_time := 1] #1 is the base level

  
  eventdata[,event_time_stratify := interaction(event_time_fact,stratify, drop=TRUE)]
  #eventdata[,treated_post_stratify := interaction(treated_post,stratify, drop = TRUE)]
  
  #eventdata[,treated_pre_stratify := interaction(treated_pre,stratify,drop=TRUE)]
  #eventdata[treated_pre==0,treated_pre_stratify := 1]
  #eventdata[event_time==base_time,treated_pre_stratify := 1]
  
  eventdata[,unitfe := .GRP, by = .(time_pair,id,treated,cohort_fact,stratify)]
  eventdata[,treated_event_time_stratify := interaction(event_time_fact,stratify, drop = TRUE)]
  
  #Omitting base year for all levels of --stratify--:
  
  base_stratify <-  paste0(c(base_time,1),collapse=".")
  
  eventdata[event_time==base_time,event_time_stratify := base_stratify]
  eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = base_stratify)]
  
  #Omitting base year for all levels of --stratify--, for treated people
  eventdata[treated == 0 ,treated_event_time_stratify :=base_stratify]
  eventdata[event_time==base_time,treated_event_time_stratify :=base_stratify]
  eventdata[,treated_event_time_stratify:=relevel(treated_event_time_stratify,ref = base_stratify)]
  
  #Omitting effect for untreated people or observations in pre-period:
  #eventdata[treated_post == 0 ,treated_post_stratify := paste0(c(0,1),collapse=".")]
  #eventdata[,treated_post_stratify:=relevel(treated_post_stratify,ref = paste0(c(0,1),collapse="."))]
  
  return(eventdata)
  
}


get_result_dynamic<-function(eventdata_panel,start,end,variable,table = data.table(), results=list(),pos=1,trends=TRUE){
  
  b = end-start+ifelse(-1 %in% start:end, 0, 1)
  for(eventtime in start:end){
    if(eventtime == -1){next}
    
    results[[pos]] <- event_ATTs_dynamic(eventdata_panel[as.character(time_pair)==eventtime,],outcomes=c(variable),keep_trends=trends)
    
    
    pos<-pos+1
  }
  for(i in 1:b){
    dt<-data.table(variable = row.names(results[[i]]$dynamic$coeftable),model=i,results[[i]]$dynamic$coeftable,obs=results[[i]]$dynamic$nobs)
    table<-rbind(dt,table)
  }
  return(table)
}

event_ATTs_dynamic<-function(eventdata,
                             outcomes,#vector of variable names
                             clustervar="id", 
                             weights="pweight",
                             keep_trends=TRUE){
  trend_call <- ifelse(keep_trends, "+ event_time_stratify |", "| event_time_stratify + ")
  call <- paste0("c(", paste0(outcomes,collapse=","), ") ~ treated_event_time_stratify ", trend_call, " unitfe")
  results<-feols(as.formula(call),
                 data = eventdata,
                 weights= eventdata[,get(weights)],
                 cluster=clustervar, lean = TRUE, mem.clean = FALSE)
  
  return(list(dynamic = results))
}

# wrapers ---------------------------------------------------------------------------------------------------------

estimate_event_dynamics <- function(panel, start, end, outcomes, control = c(), stratify = c(), use_never_treat = TRUE,
                                    timevar = "time", unitvar = "id", cohortvar = "event_time"){
  
  event_panel <- copy(panel) #copying so that the original does not change
  
  setnames(event_panel, cohortvar, "event_time")
  onset <- event_panel[, min(event_time) - 1]
  event_panel[, onset_time := onset]
  event_panel <- event_panel %>% create_event_data(timevar = timevar, unitvar = unitvar, 
                                                   cohortvar = "event_time",
                                                   onset_agevar = "onset_time",
                                                   covariate_base_balance = control,
                                                   covariate_base_stratify = stratify,
                                                   never_treat_action = ifelse(use_never_treat, "both", "exclude"),
                                                   balanced_panel = FALSE)
  
  event_panel <- event_ATTs_head(event_panel, outcomes)
  
  message("preprocessing done.")
  
  dt_dynamic <- data.table()
  
  for(outcome in outcomes){
    dt_dynamic_outcome <-data.table()
    
    dt_dynamic_outcome <- suppressMessages(get_result_dynamic(event_panel,start,-2,outcome,dt_dynamic_outcome, trends = FALSE)) #add new result to dt_dynamic
    dt_dynamic_outcome <- suppressMessages(get_result_dynamic(event_panel,0,end,outcome,dt_dynamic_outcome, trends = FALSE)) #add new result to dt_dynamic
    
    dt_dynamic_outcome[, outcome_type := outcome]
    dt_dynamic <- rbind(dt_dynamic, dt_dynamic_outcome)
    message(outcome, " done.")
  }
  
  return(dt_dynamic)
}

plot_event_study <-function(dt, graphname, note = ""){
  
  # Author : 冠儒, Hsiao, Don
  # Motified from graph 2 and graph2 0720
  
  #add the base period
  for (outcome in dt[, unique(outcome_type)]) {
    dt <- add_row(dt, variable = "treated_event_time_stratify-1.0", 
                  model = 1, Estimate = 0, "Std. Error" = 0, "t value" = 0, "Pr(>|t|)" = 0, obs = 196890,
                  outcome_type = outcome)
    dt <- add_row(dt, variable = "treated_event_time_stratify-1.1", 
                  model = 1, Estimate = 0, "Std. Error" = 0, "t value" = 0, "Pr(>|t|)" = 0, obs = 196890,
                  outcome_type = outcome)
  }
  
  significance_level <- 0.05
  dt[, c("first", "stratify_dummy") := tstrsplit(variable, ".", fixed = TRUE)]
  dt[, coefficient := str_extract(first, '([a-z\\_]+)')]
  dt[, event_time := str_extract(first, '([-0-9]+)')]
  dt[, event_time := as.numeric(event_time)]
  dt[, conf_upb := Estimate + qnorm(1-significance_level/2) * `Std. Error`]
  dt[, conf_lwb := Estimate - qnorm(1-significance_level/2) * `Std. Error`]
  
  figure <- dt[str_sub(variable, 1, 7) == "treated"] %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'red') +
    geom_line(aes(x = event_time, y = Estimate, color = "black")) + 
    geom_point(aes(x = event_time, y = Estimate, color = "black")) +
    geom_errorbar(aes(x = event_time, ymin = conf_lwb, ymax = conf_upb), 
                  width = 0.1, linetype = "dashed") +
    facet_wrap(~ outcome_type, scales = "free") +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.background = element_rect(linetype = "dashed", color = "black"),
          legend.box = "horizontal",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0)) +
    #scale_y_continuous(name = ylabel, limits = yscale, breaks = ybks) +
    #scale_x_continuous(name = xlabel, limits = xscale, breaks = xbks) +
    labs(title = graphname, subtitle = note)
  
  return(figure)
  
}


# BELOW IS NOT USED -------------------------------------------------------------------------------------------------------