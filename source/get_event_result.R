
# estimation ---------------------------------------

get_event_result<-function(eventdata,
                           variable,
                           result_type,
                           covariate = NULL,
                           clustervar="id", 
                           weights="pweight",
                           base_time = -1,
                           trends=FALSE,
                           mem.clean = TRUE,
                           separate_stratify = TRUE)
{
  
  if(is.null(covariate) & result_type == "covariates"){stop("please provide covariate when getting covariates result")}
  eventdata |> validate_eventdata(variable)
  
  eventdata <- eventdata[!is.na(get(weights))]
  
  eventdata <- eventdata |> construct_event_variable(result_type, variable, covariate, base_time)
  
  call <- get_estimate_call(result_type, variable, trends, covariate)
  
  if(eventdata[, uniqueN(stratify)] > 1 & separate_stratify){
    
    dt <- data.table()
    for(stratify_type in eventdata[, unique(stratify)]){
      
      strat_eventdata <- eventdata[stratify == stratify_type]
      
      strat_results<-feols(as.formula(call),
                     data = strat_eventdata,
                     weights= strat_eventdata[,get(weights)],
                     cluster=clustervar, lean = TRUE, mem.clean = mem.clean)
      
      strat_dt <- parse_event_result(strat_results, variable, result_type)
      strat_dt[, stratify := stratify_type]
      dt <- rbind(dt, strat_dt)
    }
    
  } else {
    
    results<-feols(as.formula(call),
                   data = eventdata,
                   weights= eventdata[,get(weights)],
                   cluster=clustervar, lean = TRUE, mem.clean = mem.clean)
    dt <- parse_event_result(results, variable, result_type)
    
  }
  
  return(dt)
  
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

validate_eventdata <- function(eventdata_panel, variable){
  if(nrow(eventdata_panel)==0){ 
    stop("event panel is empty, returning empty data.table")
  }
  
  if(any(eventdata_panel[,lapply(.SD, stats::var), .SDcols = variable] == 0)){
    stop("some outcome have no variation, returning empty data.table")
  }
}

get_event_time <- function(x){
  start <- str_locate(x, "stratify")[2]
  end <- str_locate(x, "\\.")[1]
  return(as.numeric(str_sub(x, start + 1, end - 1)))
}

parse_event_result <- function(results, variable, result_type){

  if(length(variable) == 1){
    
    dt<-data.table(variable = row.names(results$coeftable),results$coeftable,obs=results$nobs)
    
  } else {
    
    dt <- data.table()
    for(i in 1:length(variable)){
      
      result <- results[[i]]
      dt <- rbind(dt, data.table(variable = row.names(result$coeftable),result$coeftable,obs=result$nobs,outcome = variable[i]))
    }  
    
  }
  
  dt[,result:=result_type]
  dt[,event_time:= lapply(variable, get_event_time)]
  dt[,event_time := unlist(event_time)]
  
  return(dt)
  
}






