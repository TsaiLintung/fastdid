plot_event_dynamics <-function(dt, graphname = "event study plot", note = ""){
  
  #add the base period
  for (out_name in dt[, unique(outcome)]) {
    dt <- rbind(dt, data.table(outcome = out_name, event_time = -1, att = 0, se = 0, obs = 0))
  }
  
  significance_level <- 0.05
  dt[, conf_upb := att + qnorm(1-significance_level/2) * se]
  dt[, conf_lwb := att - qnorm(1-significance_level/2) * se]
  
  figure <- dt %>%
    ggplot() +
    geom_hline(yintercept = 0, linetype = 'dashed', col = 'red') +
    geom_line(aes(x = event_time, y = att, color = "black")) + 
    geom_point(aes(x = event_time, y = att, color = "black")) +
    geom_errorbar(aes(x = event_time, ymin = conf_lwb, ymax = conf_upb), 
                  width = 0.1, linetype = "dashed") +
    facet_wrap(~ outcome, scales = "free") +
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