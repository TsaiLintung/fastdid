

library(data.table)
library(did)

sim_did <- function(sample_size, time_period, untreated_prop = 0.3, cov = "no", hetero = "dynamic",
                    second_outcome = TRUE, na = "none", balanced = TRUE){
  
  #unit  -------------
  dt_i <- data.table(unit = 1:sample_size)
  if(cov == "int"){
    dt_i[, x := sample.int(5, sample_size, replace = TRUE)] #event code can't take double 
  } else if (cov == "no"){
    dt_i[, x := 1] 
  } else if (cov == "cont"){
    dt_i[, x := rnorm(sample_size)]
  }

  dt_i[, treat_latent := x*0.2 + rnorm(sample_size)]
  
  #treatment assignment ---------------------
  
  #assign treated group based on a latent related to X
  #unit with larger X tend to be treated and treated earlier
  #unit start getting treated in t = 2
  untreated_thres <- quantile(dt_i$treat_latent, untreated_prop)
  dt_i[treat_latent <= untreated_thres, G := Inf]
  cohort_prop <- (1-untreated_prop)/(time_period-1)
  last_treat_thres <- untreated_thres
  for(t in time_period:2){
    treat_thres <- quantile(dt_i$treat_latent, untreated_prop + cohort_prop*(time_period - t + 1))
    dt_i[treat_latent <= treat_thres & treat_latent > last_treat_thres, G := t]
    last_treat_thres <- treat_thres
  }
  rm(t)
  
  #assign unit FE
  dt_i[, unit_fe := rnorm(sample_size)]
  
  #time ------------------
  dt_t <- data.table(time = 1:time_period)
  dt_t[, time_fe := rnorm(time_period)]
  dt_t[, x_trend := rnorm(time_period)]
  
  #panel --------------------------
  dt <- CJ(unit = 1:sample_size, time = 1:time_period)
  dt <- dt |> merge(dt_i, by = "unit")
  dt <- dt |> merge(dt_t, by = "time")
  
  dt[, D := as.integer(time >= G)]
  
  #untreated potential outcomes
  dt[, y0 := unit_fe + time_fe + x*x_trend + rnorm(sample_size*time_period, sd = 0.1)]
  
  #generate gtatt
  att <- CJ(G = 1:time_period, time = 1:time_period)
  if(hetero == "all"){
    att[, attgt := rnorm((time_period-1)*(time_period-1), mean = 2, sd = 0.5)]
  } else if (hetero == "dynamic"){
    for(event_t in 0:max(att[,time-G])){
      att[time - G == event_t, attgt := rnorm(1, mean = 2, sd = 0.5)]
    }
  }
  att[time < G, attgt := 0] #no anticipation
  
  dt <- dt |> merge(att, by = c("time", "G"), all.x = TRUE, all.y = FALSE)
  dt[is.na(attgt), attgt := 0]
  dt[, tau := attgt + rnorm(n = sample_size*time_period, sd = 0.1)]
  
  #potential outcome
  dt[, y1 := y0 + tau]
  dt[, y := y1*D + y0*(1-D)]
  dt <- dt[, .(time, G, unit, x, y)]
  
  #additional -----------------
  
  if(na == "y"){
    dt[, y := na_insert(y)]
  } else if (na == "x") {
    dt[, x := na_insert(x)]
  } else if( na == "both"){
    dt[, y := na_insert(y)]
    dt[, x := na_insert(x)]
  }
  
  if(balanced == FALSE){
    size <- fnrow(dt)
    dt <- dt[sample(1:size, size*0.95)]
  }
  if(second_outcome == TRUE){  dt[, y2 := y + 1 + rnorm(fnrow(dt), 0, 0.1)]}
  
  return(list(dt = dt, att = att))
  
}

validate_att_est <- function(true_att, att_hat, att_se_hat, type = "all"){
  
  
  if(type == "all"){
    
    #input vector like G,time = 2,2, 2,3,......, 3,2
    
    att <- CJ(G = min(true_att$G):max(true_att$G), time = min(true_att$time):max(true_att$time))
    att[, att_hat := att_hat]
    att[, att_se_hat := att_se_hat]
    att <- att |> merge(true_att, by = c("G", "time"))
    
  } else if (type == "dynamic"){
    
    true_att[, event_time := time-G]
    true_att <- true_att[, .(attgt = mean(attgt)), by = "event_time"]
    setorder(true_att, event_time)

    
    att <- true_att
    att <- att[event_time != -1 & event_time != max(event_time)] #base time is not estimated
    att[, att_hat := att_hat]
    att[, att_se_hat := att_se_hat]
    
  }
  

  
  att[, ci_ub := att_hat+att_se_hat*1.96]
  att[, ci_lb := att_hat-att_se_hat*1.96]
  att[, par_in_ci := (attgt <= ci_ub & attgt >= ci_lb)]

  ratio <- round(att[, mean(par_in_ci)], 2)
  
  message(ratio, " of true parameter is in the CI")
  
  return(att)
  
}




