plot_event_dynamics <-function(raw_dt, 
                               graphname = "event study plot", note = "", base_time = -1, significance_level = 0.05, stratify_offset =0.1){
  
  dt <- copy(raw_dt)
  
  dt[, c("variable", "t value", "Pr(>|t|)", "obs", "result") := NULL]
  setnames(dt, c("Estimate", "Std. Error"), c("att", "se"))
  
  #add the base period
  base_row <- data.table(att = 0, se = 0, event_time = base_time)
  dt <- dt |> add_base_row(base_row)
  
 

  dt[, conf_upb := att + qnorm(1-significance_level/2) * se]
  dt[, conf_lwb := att - qnorm(1-significance_level/2) * se]
  
  #add some offset
  if("stratify" %in% names(dt)){
    dt[, event_time := event_time + stratify_offset*(as.integer(stratify)-1)]
  }
  
  figure <- dt %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'red')
  
  if("stratify" %in% names(dt)){
    figure <- figure + geom_line(aes(x = event_time, y = att, color = stratify)) + 
      geom_point(aes(x = event_time, y = att, color = stratify)) +
      geom_errorbar(aes(x = event_time, ymin = conf_lwb, ymax = conf_upb), 
                    width = 0.1, linetype = "dashed")
  } else {
    figure <- figure + geom_line(aes(x = event_time, y = att), color = "black") + 
    geom_point(aes(x = event_time, y = att), color = "black") +
    geom_errorbar(aes(x = event_time, ymin = conf_lwb, ymax = conf_upb), 
                  width = 0.1, linetype = "dashed")
  }
  
  if("outcome" %in% names(dt)){
    figure <- figure + facet_wrap(~ outcome, scales = "free")
  }
  if("cohort" %in% names(dt)){
    figure <- figure + facet_wrap(~ cohort, scales = "free")
  }

  figure <- figure +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.background = element_rect(linetype = "dashed", color = "black"),
          legend.box = "horizontal",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0)) +
    labs(title = graphname, subtitle = note)
  
  return(figure)
  
}

add_base_row <- function(dt, base_row){
  
  col_names <- names(dt)
  opt_cols <- c("outcome", "stratify", "cohort")
  
  for(opt_col in opt_cols){
    if(opt_col %in% col_names){
      
      row_list <- list()
      
      for (col_name in dt[, unique(get(opt_col))]) {
        base_row_sub <- copy(base_row)
        base_row_sub[, c(opt_col) := col_name]
        row_list <- c(row_list, list(base_row_sub))
      }
      
    } 
    base_row <- rbindlist(row_list)
  }

  dt <- dt |> rbind(base_row)
  
  return(dt)
}