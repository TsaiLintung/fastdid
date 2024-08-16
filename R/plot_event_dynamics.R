#' Create an event study plot for Difference-in-Differences (DiD) analysis.
#'
#' This function generates an event study plot based on the results of a DiD analysis.
#'
#' @param dt A data table containing the results of the DiD analysis. It should include columns for 'att' (average treatment effect), 'se' (standard error), and 'event_time' (time points).
#'
#' @return A ggplot2 object representing the event study plot.
#' @export
plot_did_dynamics <-function(dt ){
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("The ggplot2 package must be installed to use plotting functions")   #Either exit or do something without rgl
    return(NULL)
  }
  
  #find the base_period
  et_range <- min(dt[, event_time]):max(dt[, event_time])
  base_time <- et_range[!et_range %in% dt[, unique(event_time)]]
  if(length(base_time)!=1){stop("missing more then one period")}
  
  #add the base period
  if("outcome" %in% names(dt)){
    base_row <- data.table(att = 0, se = 0, event_time = base_time, outcome = dt[, unique(outcome)], att_ciub = 0, att_cilb = 0)
  } else {
    base_row <- data.table(att = 0, se = 0, event_time = base_time, att_ciub = 0, att_cilb = 0)
  }
  
  #dt <- dt |> add_base_row(base_row)
  dt <- dt |> rbind(base_row, fill = TRUE)
  

  figure <- dt |> 
    ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = 'dashed', col = 'red')
  
  figure <- figure + ggplot2::geom_line(ggplot2::aes(x = event_time, y = att), color = "black") + 
    ggplot2::geom_point(ggplot2::aes(x = event_time, y = att), color = "black") +
    ggplot2::geom_errorbar(ggplot2::aes(x = event_time, ymin = att_cilb, ymax = att_ciub), 
                width = 0.1, linetype = "dashed")

  if("outcome" %in% names(dt)){
    figure <- figure + ggplot2::facet_wrap(~outcome, scales = "free")
  }

  return(figure)
  
}
