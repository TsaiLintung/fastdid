#' Plot event study
#'
#' Plot event study results.
#'
#' @param x A data table generated with [fastdid] with one-dimensional index.
#' @param margin character, the x-axis of the plot
#'
#' @return A ggplot2 object
#' @examples
#' 
#' # simulated data
#' simdt <- sim_did(1e+02, 10, seed = 1)
#' dt <- simdt$dt
#' 
#' #estimation
#' result <- fastdid(data = dt, timevar = "time", cohortvar = "G", 
#'                   unitvar = "unit", outcomevar = "y",  
#'                   result_type = "dynamic")
#' 
#' #plot
#' plot_did_dynamics(result)
#' 
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
  } else {
    x <- x[type == "post"]
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