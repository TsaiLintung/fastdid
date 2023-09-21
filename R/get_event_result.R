#' Get Event Result
#'
#' Calculate various types of event results from an event dataset.
#'
#' @param eventdata The event dataset in data.table format.
#' @param variable The variable of interest for the analysis.
#' @param result_type The type of result to calculate. Options include
#'        "cohort_event_time", "dynamic", "pooled", "means", and "covariate".
#' @param covariate An optional covariate variable to include in the analysis.
#' @param clustervar The clustering variable for grouped data. Default is "id".
#' @param weights The weight variable for weighted analysis. Default is "pweight".
#' @param base_time The base time for dynamic analysis. Default is -1.
#' @param trends Logical, indicating whether to include trends in the analysis.
#' @param mem.clean Logical, indicating whether to perform memory cleanup during FEOLS.
#' @param separate_stratify Logical, indicating whether to separate stratified cohorts.
#' @param separate_cohort_time Logical, indicating whether to separate cohort time.
#'
#' @return A data.table containing the calculated event results.
#'
#' @examples
#' # Example usage:
#' result <- get_event_result(eventdata, "variable_name", "cohort_event_time")
#'
#' @import data.table
#'
#' @export
get_event_result <- function(eventdata,
                             variable,
                             result_type,
                             covariate = NULL,
                             clustervar = "id",
                             weights = "pweight",
                             base_time = -1,
                             trends = FALSE,
                             mem.clean = FALSE,
                             separate_stratify = TRUE,
                             separate_cohort_time = TRUE) {
  #allow list input (useful when not combined)
  if(!is.data.table(eventdata)){
    if(!is.data.table(eventdata[[1]])){stop("please provide a data.table")}
    if(result_type != "cohort_event_time"){stop("by-cohort list is only available when using 'cohort_event_time'")}
    
    all_results <- data.table()
    for(i in seq(1, length(eventdata))){
      results <- get_event_result(eventdata = eventdata[[i]],
                                  variable = variable,
                                  result_type = result_type,
                                  covariate = covariate,
                                  clustervar=clustervar, 
                                  weights=weights,
                                  base_time = base_time,
                                  trends=trends,
                                  mem.clean = mem.clean,
                                  separate_stratify = separate_stratify,
                                  separate_cohort_time = separate_cohort_time)
      all_results <- rbind(all_results, results)
    }
    
    return(all_results)
  } 
  
  
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
  
  results <- eventdata |> estimate_att(call, weights, clustervar, result_type, 
                                       mem.clean, separate_cohort_time, separate_stratify)
  
  return(results)
  
}

# helper functions --------------------------------------------

construct_event_variable <- function(eventdata, result_type, variable, covariate, base_time){
  
  base_event_stratify <- paste0(c(base_time,1),collapse=".")
  
  eventdata[, splitvar := 1]
  
  if(result_type == "cohort_event_time"){
    
    eventdata[,unitfe := finteraction(cohort, time_pair,treated,stratify)] #the level difference between treated and untreated units within each time_pair and stratify
    eventdata[, splitvar := finteraction(cohort, time_pair)]
    
    
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
    eventdata[,splitvar := time_pair]
    
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


estimate_att <- function(eventdata, call, weights, clustervar, result_type, 
                         mem.clean,separate_cohort_time, separate_stratify){
  all_results_dt <- data.table()
  
  if(!eventdata[, uniqueN(stratify)] > 1 | !separate_stratify){
    
    if(!separate_cohort_time){
      results<-feols(as.formula(call),
                     data = eventdata,
                     weights= eventdata[,get(weights)],
                     split = ~splitvar,
                     cluster=clustervar, lean = TRUE, mem.clean = mem.clean, notes = FALSE)
      
      all_results_dt <- parse_event_result(results, result_type)
      
    } else {
      
      #split by cohort / cohort+event_time
      for(split in eventdata[, unique(splitvar)]){
        
        eventdata_split <- eventdata[splitvar == split]
        results_split <- feols(as.formula(call),
                               data = eventdata_split,
                               weights= eventdata_split[,get(weights)],
                               split = ~splitvar,
                               cluster=clustervar, lean = TRUE, mem.clean = mem.clean, notes = FALSE)
        results_split_dt <- parse_event_result(results_split, result_type)
        all_results_dt <- rbind(all_results_dt, results_split_dt)
        
      }      
      
    }
    
  } else {
    
    for(stratify_type in eventdata[, unique(stratify)]){
      
      strat_eventdata <- eventdata[stratify == stratify_type]
      
      if(!separate_cohort_time){
        strat_results<-feols(as.formula(call),
                             data = strat_eventdata,
                             weights= strat_eventdata[,get(weights)],
                             split = ~splitvar,
                             cluster=clustervar, lean = TRUE, mem.clean = mem.clean, notes = FALSE)
        strat_results_dt <- parse_event_result(strat_results, result_type)
        strat_results_dt[, stratify := stratify_type]
        all_results_dt <- rbind(all_results_dt, strat_results_dt)
        
      }  else {
        
        #split by cohort / cohort+event_time
        strat_results <- list()
        for(split in strat_eventdata[, unique(splitvar)]){
          
          strat_eventdata_split <- strat_eventdata[splitvar == split]
          strat_results_split <- feols(as.formula(call),
                                 data = strat_eventdata_split,
                                 split = ~splitvar,
                                 weights= strat_eventdata_split[,get(weights)],
                                 cluster=clustervar, lean = TRUE, mem.clean = mem.clean, notes = FALSE)
          strat_results_split_dt <- parse_event_result(strat_results_split, result_type)
          strat_results_split_dt[, stratify := stratify_type]
          all_results_dt <- rbind(all_results_dt, strat_results_split_dt)
          
        }      
 
      }

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
  
  dot_pos <- str_locate_all(x, "\\.")
  start <- ifelse(length(dot_pos[[1]])==0,str_locate(x, "stratify")[2], dot_pos[[1]][length(dot_pos), 1])
  end <-str_length(x)
  return(str_sub(x, start + 1, end))
  
}

parse_event_result <- function(results, result_type){
  
  dt <- data.table()
  for(j in seq(1, length(results))){
    est <- results[[j]]$coeftable
    dt_row <- data.table(variable = rownames(est), est, obs = results[[j]]$nobs, outcome = deparse(results[[j]]$fml[[2]]))
    dt <- rbind(dt, dt_row)
  }

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
  if(all(dt[, stratify == 1], na.rm = TRUE)){dt[, stratify := NULL]}
  
  return(dt)
  
}






