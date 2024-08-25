#' @export
plot.fastdid_result <- function(x,...){
  x <- x$estimate
  #dispatch to functions by result_type
  if("event_time" %in% names(x)){
    plot <- plot_did_dynamics(x)
  } else if ("cohort" %in% names(x) & "time" %in% names(x)){
    plot <- plot_did_group_time(x)
  } else if ("cohort" %in% names(x)){
    plot <- plot_did_dynamics(x, "cohort")
  } else if ("time" %in% names(x)){
    plot <- plot_did_dynamics(x, "time")
  }
  
  if(x[,uniqueN(outcome)] > 1){
    plot <- plot + facet_wrap(~outcome, scales = "free")
  }
  
  return(plot)
}

#' Create plot for Difference-in-Differences (DiD) analysis.
#'
#' This function generates an plot based on the results of a DiD analysis.
#'
#' @param x A data table containing the results of the DiD analysis. It should include columns for 'att' (average treatment effect), 'se' (standard error), and 'event_time' (time points).
#' @param margin the x-axis of the plot
#'
#' @return A ggplot2 object representing the event study plot.
#' @export
plot_did_dynamics <-function(x, margin = "event_time"){
  
  #find the base_period
  if(margin == "event_time"){
    et_range <- min(x[, event_time]):max(x[, event_time])
    base_time <- et_range[!et_range %in% x[, unique(event_time)]]
    if(length(base_time)!=1){stop("missing more than one period")}
    
    #add the base period
    if("outcome" %in% names(x)){
      base_row <- data.table(att = 0, se = 0, event_time = base_time, outcome = x[, unique(outcome)], att_ciub = 0, att_cilb = 0)
    } else {
      base_row <- data.table(att = 0, se = 0, event_time = base_time, att_ciub = 0, att_cilb = 0)
    }
    x <- x |> rbind(base_row, fill = TRUE)
  }
  
  plot <- x |> 
     ggplot() +
     geom_hline(yintercept = 0, linetype = 'dashed', col = 'red') + 
     geom_point( aes(x = eval(str2lang(margin)), y = att), color = "black") + #point est
     geom_errorbar( aes(x = eval(str2lang(margin)), ymin = att_cilb, ymax = att_ciub), 
                width = 0.1, linetype = "dashed") + #CI
     labs(x = margin)
  
  if(margin == "event_time"){
    plot <- plot + geom_line( aes(x = eval(str2lang(margin)), y = att), color = "black") #point est
  }

  return(plot)
  
}


#' @export
plot_did_group_time <- function(x, range = c(-1, 1),...) {
  x$gtexpo[time >= cohort] |>  ggplot( aes(x = as.factor(time-cohort), y = as.factor(cohort), fill = gamma)) +  geom_tile() + 
     scale_fill_distiller(palette = "RdBu", limits = range) + 
     labs(y = "cohort", x = "event_time") 
}
