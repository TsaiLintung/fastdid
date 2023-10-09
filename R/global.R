NULL

if(FALSE){
  text <- ". agg_weight att att_cont att_treat attgt cohort cohort_size
    conf_lwb conf_upb const cont_ipw_weight count delta_y element_rect
    element_text event_time pg placeholder post.y pre.y ps s se target
    tau time_fe treat_ipw_weight treat_latent type unit unit_fe weight x
    x2 x_trend y y0 y1 y2 time"
  text <- text |> str_remove_all("\\\n") |>str_split(" ") |> unlist()
  text <- text[text!=""]
  text <- text |> str_flatten(collapse = "','")
  text <- paste0("c('", text, "')")
}

# quiets concerns of R CMD check re: the .'s that appear in pipelines and data.table variables
utils::globalVariables(c('.','agg_weight','att','att_cont','att_treat','attgt','cohort','cohort_size','conf_lwb','conf_upb',
                         'const','cont_ipw_weight','count','delta_y','element_rect','element_text','event_time','pg','placeholder',
                         'post.y','pre.y','ps','s','se','target','tau','time_fe',
                         'treat_ipw_weight','treat_latent','type','unit','unit_fe','weight','x','x2','x_trend','y','y0','y1','y2', 'time', 'weights'))