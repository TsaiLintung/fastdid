#' Create an event study plot for Difference-in-Differences (DiD) analysis.
#'
#' This function generates an event study plot based on the results of a DiD analysis.
#'
#' @param dt A data table containing the results of the DiD analysis. It should include columns for 'att' (average treatment effect), 'se' (standard error), and 'event_time' (time points).
#' @param graphname A character string specifying the title of the plot (default is "event study plot").
#' @param note A character string for adding additional notes or comments to the plot (default is empty).
#' @param base_time The time point representing the base period (default is -1).
#' @param significance_level The significance level for confidence intervals (default is 0.05).
#'
#' @return A ggplot2 object representing the event study plot.
#' @export

plot_did_dynamics <-function(dt, 
                             graphname = "event study plot", note = "", base_time = -1, significance_level = 0.05#, 
                             #stratify_offset =0.1
){
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("The ggplot2 package must be installed to use plotting functions")
    #Either exit or do something without rgl
    return(NULL)
  }
  

  #add the base period
  if("outcome" %in% names(dt)){
    base_row <- data.table(att = 0, se = 0, event_time = base_time, outcome = dt[, unique(outcome)])
  } else {
    base_row <- data.table(att = 0, se = 0, event_time = base_time)
  }
  
  #dt <- dt |> add_base_row(base_row)
  dt <- dt |> rbind(base_row, fill = TRUE)
  
  dt[, conf_upb := att + qnorm(1-significance_level/2) * se]
  dt[, conf_lwb := att - qnorm(1-significance_level/2) * se]
  

  figure <- dt |> 
    ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = 'dashed', col = 'red')
  
  figure <- figure + ggplot2::geom_line(ggplot2::aes(x = event_time, y = att), color = "black") + 
    ggplot2::geom_point(ggplot2::aes(x = event_time, y = att), color = "black") +
    ggplot2::geom_errorbar(ggplot2::aes(x = event_time, ymin = conf_lwb, ymax = conf_upb), 
                width = 0.1, linetype = "dashed")

  if("outcome" %in% names(dt)){
    figure <- figure + ggplot2::facet_wrap(~outcome, scales = "free")
  }

  figure <- figure + ggplot2::theme_classic()

  return(figure)
  
}
