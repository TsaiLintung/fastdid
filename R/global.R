NULL

# quiets concerns of R CMD check re: the .'s that appear in pipelines and data.table variables
utils::globalVariables(c('.','agg_weight','att','att_cont','att_treat','attgt','cohort','cohort_size','conf_lwb','conf_upb',
                         'const','cont_ipw_weight','count','delta_y','element_rect','element_text','event_time','pg','placeholder',
                         'post.y','pre.y','ps','s','se','target','tau','time_fe',
                         'treat_ipw_weight','treat_latent','type','unit','unit_fe','weight','x','x2',
                         'x_trend','y','y0','y1','y2', 'time', 'weights', 'outcome', "G", "D", 'xvar',
                         'V1','att_cont_post','att_cont_pre','att_treat_post','att_treat_pre','inpost','inpre','max_et','min_et','new_unit','or_delta','or_delta_post','or_delta_pre','targeted','used',
                         "timevar", "cohortvar", "unitvar", "outcomevar", "control_option", "result_type", "balanced_event_time", "control_type",
                         "allow_unbalance_panel", "boot", "biters", "weightvar", "clustervar", "covariatesvar", "varycovariatesvar", "filtervar",
                         "copy", "validate", "max_control_cohort_diff", "anticipation", "min_control_cohort_diff", "base_period", "post", "att_ciub", "att_cilb", "cband", "alpha",
                         "G2", "G1", "mg", "cohort1", "cohort2", "event_time_1", "event_time_2",
                         "D2", "attgt2", "event", "atu2", "y01", "y10", "y11", "tau2", "parallel",
                         "tp", "cp", "tb", "cb", "no_na", "event_stagger", "double_control_option"))

