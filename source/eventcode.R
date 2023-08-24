#by LinTung Tsai to be used with Max's event code
#some wrappers
library(fixest)

#Treated households:
create_event_data<-function(maindata,
                            #if no covariates for balancing or stratification, just specify a constant (as follows):
                            covariate_base_stratify=c(), #vector of variable names, treatment effects are stratified by their joint unique values.
                            covariate_base_balance=c(), #vector of variable names to include in finding exact matches for treated units
                            timevar,
                            unitvar,
                            outcomevar,
                            clustervar = NULL,
                            cohortvar, #period when the particular event of interest arises
                            base_time = -1, #reference value for event time; the period from which differences over time are taken
                            never_treat_action = "both", #options: "both", "only", and "exclude"
                            lower_event_time = -Inf, #Earliest period (relative to treatment time) for which to estimate effects
                            upper_event_time = Inf, #Latest period (relative to treatment time) for which to estimate effects
                            data_validate = "check"
  ) 
{
  
  #handle names
  setnames(maindata,c(timevar,unitvar,cohortvar),c("time","id","cohort"))
  if(is.null(clustervar)){clustervar <- "id"}
  
  if(data_validate == "check"){
    
    if(!is.data.table(maindata)) stop("rawdata must be a data.table")
    if(any(is.na(maindata$cohort))) stop("cohort variable should not be missing (it can be infinite instead)")
    if(any(fduplicated(maindata[,.(time, id)]))){stop ("some unit-time is observed more than once")}
    if(any(missing_cases(maindata, cols = c("cohort", "time", "id")))){stop ("some units have missing cohort / time / id. Set cohort to Inf if it is never-treated")}
    if(any(missing_cases(maindata, cols = c(covariate_base_stratify, covariate_base_balance)))){stop ("some balance / stratify cov is missing")}
    if(any(missing_cases(maindata, cols = outcomevar))){stop ("some outcome is missing")}
    
  } else if (data_validate == "fix") { #fix the maindata as much as possible
    
    #not implemented yet
  
  } else if (data_validate != "trust") { #trust the data, perform no checks
    stop("pick either check / fix / trust for data validation strategy")
  }
  
  maindata[, id := charToFact(as.character(id))]
  
  #turn the covariates into a single interacted factor
  if(length(covariate_base_stratify) > 0) {
    maindata[, stratify:= charToFact(as.character(.GRP)), by = covariate_base_stratify]
  }else {
    maindata[,stratify:="base_strat"]
  }
  if(length(covariate_base_balance) > 0){
    maindata[, balancevars:= charToFact(as.character(.GRP)), by =covariate_base_balance]
  }else {
    maindata[,balancevars:="base_balance"]
  }

  #create treated data
  treatdata<-copy(maindata[!is.infinite(cohort),])
  treatdata[,treated:=1]
  treatdata[,event_time:=time-cohort]
  treatdata[,cohort_pair := cohort] #the pair
  treatdata<-treatdata[event_time >= lower_event_time & event_time <= upper_event_time,]
  treatdata[,treatgroup:="treated"]
  
  #stack control cohorts for each treated cohort
  #I assume people who never suffer the event have a value cohort = Inf
  control_list <- list()

  for(o in unique(treatdata$cohort_pair)){
    if(never_treat_action=="both") {
      controlcohort <- maindata[(cohort > o | is.infinite(cohort)) ,]
      controlcohort[is.infinite(cohort),treatgroup:="never-treated"]
      controlcohort[!is.infinite(cohort),treatgroup:="later-treated"]
    }
    if(never_treat_action=="exclude"){
      controlcohort <- maindata[(cohort > o & !is.infinite(cohort)) ,]
      controlcohort[,treatgroup:="later-treated"]
    }
    if(never_treat_action=="only") {
      controlcohort <- maindata[is.infinite(cohort) ,]
      controlcohort[,treatgroup:="never-treated"]
    }
    
    controlcohort[,cohort_pair := o]
    controlcohort[,event_time := time - cohort_pair]

    control_list<- c(control_list, list(controlcohort))
    
  }
  
  rm(controlcohort)
  gc()
  
  controldata <- rbindlist(control_list)
  
  controldata[,treated:=0]

  #drop someone from the control cohort when they get treated:
  controldata<-controldata[time < cohort ,]
  
  #drop control for event time outside needed
  controldata<-controldata[event_time >= lower_event_time & event_time <= upper_event_time,]
  
  #before stacking turn cohort and event time into fact
  treatdata[, event_time_fact := charToFact(as.character(event_time))]
  treatdata[, cohort_pair_fact := charToFact(as.character(cohort_pair))]
  controldata[, event_time_fact := charToFact(as.character(event_time))]
  controldata[, cohort_pair_fact := charToFact(as.character(cohort_pair))]

  #For the final dataset, we must stack each pairwise combo of (base year, other year)
  #for each household. The reason? We will assign control households weights that vary
  #based on the other year, as households enter/exit the sample.
  event_times<- as.integer(as.character(treatdata[,unique(event_time)]))
  event_times <- event_times[event_times != -1 & event_times >= lower_event_time & event_times <= upper_event_time]
  
  controldata[, `:=`(min_event_time = min(event_time),
                  max_event_time = max(event_time)), by = .(id, cohort_pair)]
  treatdata[, `:=`(min_event_time = min(event_time),
                  max_event_time = max(event_time)), by = .(id, cohort_pair)]
  
  stack_eventtime <- function(t, treatdata, controldata){
    
    #t is the event time
    
    pair_treat_data <- treatdata[t >= min_event_time & t <= max_event_time][event_time %in% c(t, base_time),]
    pair_control_data <- controldata[t >= min_event_time & t <= max_event_time][event_time %in% c(t, base_time),]
    
    pair_treat_data[,time_pair := t]
    pair_control_data[,time_pair := t]
    
    #return this
    return(c(list(pair_treat_data), list(pair_control_data)))
    
  }

  data_list <- foreach(t = event_times) %do% stack_eventtime(t, treatdata, controldata)

  data_list <- flatten(data_list)
  eventdata <- rbindlist(data_list, use.names=TRUE)
  
  eventdata[, time_pair := charToFact(as.character(time_pair))]
  
  rm(treatdata)
  rm(controldata)
  gc()
  
  #calculate ipw
  eventdata[,pweight:=NA]
  eventdata[,pval:= feols(treated ~ 1 | interaction(cohort_pair_fact,event_time_fact,time_pair,stratify,balancevars),
                          data = eventdata, lean = FALSE, combine.quick = TRUE)$fitted.values]
  
  eventdata[treated==1 & pval < 1 & pval > 0,pweight:=1.00] #att is estimated
  eventdata[, pweight := as.double(pweight)] #TSAI addition
  eventdata[treated==0 & pval < 1 & pval > 0,pweight:=pval/(1-pval)]
  
  eventdata <- eventdata[!is.na(pweight),] #ppl with pval outside 0 and 1 <-> no common support is dropped
  
  #calculating weights to match control and treated households on characteristics, across strata
  if(length(covariate_base_stratify) > 0 ){
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
  
  #original event att head
  eventdata[,event_time_stratify := interaction(event_time_fact,stratify, drop=TRUE)]
  eventdata[,unitfe := .GRP, by = .(time_pair,id,treated,cohort_pair_fact,stratify)]
  eventdata[,treated_event_time_stratify := interaction(event_time_fact,stratify, drop = TRUE)]
  
  base_stratify <- paste0(c(base_time,1),collapse=".")
  
  eventdata[event_time==base_time,event_time_stratify := base_stratify]
  eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = base_stratify)]
  
  #Omitting base year for all levels of --stratify--, for treated people
  eventdata[treated == 0 ,treated_event_time_stratify := base_stratify]
  eventdata[event_time==base_time,treated_event_time_stratify :=base_stratify]
  eventdata[,treated_event_time_stratify:=relevel(treated_event_time_stratify,ref = base_stratify)]
  
  keep_columns <- c(outcomevar, clustervar, "event_time_stratify", "treated_event_time_stratify", "unitfe", "time_pair", "pweight")
  eventdata <- eventdata[, .SD, .SDcols = keep_columns]
  gc()
  
  return(eventdata)
  
}

get_result_dynamic<-function(eventdata_panel, variable, clustervar = "id", weights = "pweight"){

  call <- paste0("c(", paste0(variable,collapse=","), ") ~ treated_event_time_stratify | event_time_stratify + unitfe")
  results<-feols(as.formula(call),
                 data = eventdata_panel,
                 weights= eventdata_panel[,get(weights)],
                 split = "time_pair",
                 cluster=clustervar, lean = TRUE, mem.clean = FALSE)
  table <- data.table()
  for(result in results){
    dt<-data.table(variable = row.names(result$coeftable),result$coeftable,obs=result$nobs)
    table<-rbind(dt,table)
  }
  table[, event_time := as.integer(str_remove_all(str_extract(variable, "y(.*?)\\."), "y|\\."))]
  setnames(table, c("Estimate", "Std. Error"), c("att", "se"))
  setorder(table, event_time)
  table <- table[,.(event_time, att, se, obs)]
  return(table)
  
}

# wrapers ---------------------------------------------------------------------------------------------------------

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