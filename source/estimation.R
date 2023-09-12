
# estimation ---------------------------------------

event_ATTs_cohort_event_time <-function(eventdata,
                                  outcomes, #vector of variable names
                                  clustervar="id", 
                                  weights="pweight",
                                  keep_trends=TRUE,
                                  se="cluster",
                                  lean=TRUE,
                                  ssc=NULL,
                                  mem.clean = TRUE){
  

  base_event_stratify <- paste0(c(eventdata[, max(base_time)],1),collapse=".")

  eventdata[,unitfe := finteraction(cohort,time_pair,treated,id)]
  
  eventdata[,cohort_event_time_event_time_stratify:= finteraction(event_time_fact,cohort,stratify)]
  eventdata[,treated_cohort_event_time_event_time_stratify := cohort_event_time_event_time_stratify]
  
  #Omitting base year for all levels of --stratify--:
  eventdata[event_time==base_time, event_time_stratify := base_event_stratify]
  eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = base_event_stratify)]
  
  #Omitting base year for all levels of --stratify--, for treated people
  eventdata[event_time==base_time | treated == 0 ,treated_cohort_event_time_event_time_stratify := base_event_stratify]
  eventdata[,treated_cohort_event_time_event_time_stratify:=relevel(treated_cohort_event_time_event_time_stratify,ref = base_event_stratify)]
  
  #construct the call
  outcomes_call <- paste0("c(", paste0(outcomes,collapse=","), ")")
  event_stratify_call <- ifelse(keep_trends, "~ cohort_event_time_event_time_stratify + treated_cohort_event_time_event_time_stratify |", 
                                "~ treated_cohort_event_time_event_time_stratify | cohort_event_time_event_time_stratify +")
  call <- paste0(outcomes_call, event_stratify_call, " unitfe" )
  
  results<-feols(as.formula(call),
                 data = eventdata,
                 weights= eventdata[,get(weights)],
                 #split = "cohort",
                 cluster=clustervar, lean = TRUE, mem.clean = mem.clean)
  
  return(list(cohort_event_time = results))
  
}



event_ATTs_dynamic<-function(eventdata,
                             outcomes,#vector of variable names
                             clustervar="id", 
                             weights="pweight",
                             keep_trends=TRUE,
                             se="cluster",
                             lean=TRUE,
                             ssc=NULL,
                             mem.clean = TRUE){

  #construct the call
  outcomes_call <- paste0("c(", paste0(outcomes,collapse=","), ")")
  event_stratify_call <- ifelse(keep_trends, "~ event_time_stratify + treated_event_time_stratify |", 
                                "~ treated_event_time_stratify | event_time_stratify +")
  call <- paste0(outcomes_call, event_stratify_call, " unitfe" )

  results<-feols(as.formula(call),
                 data = eventdata,
                 weights= eventdata[,get(weights)],
                 cluster=clustervar, lean = TRUE, mem.clean = mem.clean)
  
  return(list(dynamic = results))
  
}


event_ATTs_means<-function(eventdata,
                           outcomes,#vector of variable names
                           clustervar="id", 
                           weights="pweight",
                           keep_trends=TRUE,
                           se="cluster",
                           lean=TRUE,
                           ssc=NULL,
                           mem.clean = TRUE){
  
  eventdata[, stratify:= factor(eventdata$stratify, levels = c(levels(eventdata$stratify), "empty"))]
  
  estimate_data <- eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),]
  
  call <- paste0("c(",paste0(outcomes,collapse=",") ,") ~ interaction(stratify) - 1")

  results_means<-feols(as.formula(call),
  data = estimate_data,
  weights= estimate_data[,get(weights)],
  cluster=clustervar, lean = TRUE, mem.clean = mem.clean)
  
  return(list(means = results_means))
}


event_ATTs_pooled<-function(eventdata,
                            outcomes,#vector of variable names
                            clustervar="id", 
                            weights="pweight",
                            keep_trends=TRUE,
                            se="cluster",
                            lean=TRUE,
                            ssc=NULL,
                            mem.clean = TRUE){
  
  eventdata[,treated_post := qF((treated == 1) * (post == 1))]
  eventdata[,treated_pre := qF((treated == 1) * (post == 0))]
  
  eventdata[,treated_post_stratify := finteraction(treated_post,stratify)]
  eventdata[,treated_pre_stratify := finteraction(treated_pre,stratify)]
  
  eventdata[treated_pre==0 | event_time==base_time,treated_pre_stratify := 1]
  
  #Omitting effect for untreated people or observations in pre-period:
  eventdata[treated_post == 0 ,treated_post_stratify := paste0(c(0,1),collapse=".")]
  eventdata[,treated_post_stratify:=relevel(treated_post_stratify,ref = paste0(c(0,1),collapse="."))]
  
  outcomes_call <- paste0("c(", paste0(outcomes,collapse=","), ")")
  event_stratify_call <- ifelse(keep_trends, "~ event_time_stratify + treated_pre_stratify + treated_post_stratify |", 
                                "~ treated_pre_stratify + treated_post_stratify | event_time_stratify +")
  call <- paste0(outcomes_call, event_stratify_call, " unitfe" )
  
  results_pooled <-feols(as.formula(call),
                 data = eventdata,
                 weights= eventdata[,get(weights)],
                 cluster=clustervar, lean = TRUE, mem.clean = mem.clean)

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
                                ssc=NULL,
                                mem.clean = TRUE){
  
  if(is.null(covariates)){
    stop("covariates is null")
  }
  
  #eventdata$stratify <- factor(eventdata$stratify, levels = c(levels(eventdata$stratify), "empty"))
  eventdata[,int_event_time_stratify:=interaction(event_time,stratify)]
  
  call <- paste0("c(", paste0(covariates,collapse=",") , ") ~ int_event_time_stratify + int_event_time_stratify:treated- 1" )
  
  estimate_data <- eventdata[as.numeric(as.character(event_time)) == as.numeric(as.character(time_pair)) | as.numeric(as.character(time_pair)) == 0,]
  
  results_covariates <- feols(as.formula(call),
      data = estimate_data,
      weights= estimate_data[,get(weights)],
      se=se,
      cluster=clustervar, lean = lean, ssc=ssc, mem.clean = mem.clean)

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
  
  results_dynamic <- event_ATTs_dynamic(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  results_pooled <- event_ATTspooled(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  results_means <- event_ATTs_means(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  results_covariates <- event_ATTs_covariates(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  
  return(list(pooled = results_pooled,
              dynamic = results,
              means = results_means,
              covariates = results_covariates))
}

# get result -------------------

get_result_cohort_event_time<-function(eventdata_panel,variable,trends=TRUE,mem.clean = TRUE){
  
  validate_eventdata(eventdata_panel, variable)
  
  results<-event_ATTs_cohort_event_time(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  
  dt <- parse_event_ATTs(results, variable, "cohort_event_time")
  
  return(dt)
  
}

get_result_dynamic<-function(eventdata_panel,variable,trends=TRUE, mem.clean = TRUE){
  
  validate_eventdata(eventdata_panel, variable)

  results<-event_ATTs_dynamic(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  
  dt <- parse_event_ATTs(results, variable, "dynamic")

  return(dt)
  
}

get_result_pooled<-function(eventdata_panel,variable,trends=TRUE,
                            mem.clean = TRUE){
  
  validate_eventdata(eventdata_panel, variable)
  
  results<-event_ATTs_pooled(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  
  dt <- parse_event_ATTs(results, variable, "pooled")
 
}

get_result_means<-function(eventdata_panel,variable,trends=TRUE,
                           mem.clean = TRUE){
  
  validate_eventdata(eventdata_panel, variable)
  
  results<-event_ATTs_means(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  
  dt <- parse_event_ATTs(results, variable, "means")
  
  return(dt)
  
}

get_result_covariates<-function(eventdata_panel,covariate,variable=NULL,trends=TRUE,
                                mem.clean = TRUE){
  
  validate_eventdata(eventdata_panel, variable)
  
  results<-event_ATTs_covariates(eventdata_panel,covariates=covariate,outcomes = c(variable),keep_trends = trends)
  
  dt <- parse_event_ATTs(results, variable, "covariates")
  
  dt[,model:=paste0("additional_covariate--",covariate)]
  
  return(dt)
  
}

# helper functions --------------------------------------------

validate_eventdata <- function(eventdata_panel, variable){
  if(nrow(eventdata_panel)==0){ 
    stop("event panel is empty, returning empty data.table")
  }
  
  if(any(eventdata_panel[,lapply(.SD, stats::var), .SDcols = variable] == 0)){
    stop("some outcome have no variation, returning empty data.table")
  }
}

parse_event_ATTs <- function(results, variable, result_type){
  
  results <- results[[result_type]]
  
  if(length(variable) > 1){
    dt <- data.table()
    for(i in 1:length(variable)){
      
      result <- results[[i]]
      dt <- rbind(dt, data.table(variable = row.names(result$coeftable),result$coeftable,obs=result$nobs,outcome = variable[i]))
      
    }
  } else {dt<-data.table(variable = row.names(results$coeftable),results$coeftable,obs=results$nobs)}
  
  dt[,result:=result_type]
  
  get_event_time <- function(x){
    start <- str_locate(x, "stratify")[2]
    end <- str_locate(x, "\\.")[1]
    return(as.numeric(str_sub(x, start + 1, end - 1)))
  }
  
  dt[,event_time:= lapply(variable, get_event_time)]
  dt[,event_time := unlist(event_time)]
}






