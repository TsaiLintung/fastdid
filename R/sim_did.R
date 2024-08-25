#' Simulate a Difference-in-Differences (DiD) dataset
#'
#' Simulates a dataset for a Difference-in-Differences analysis with various customizable options.
#'
#' @param sample_size The number of units in the dataset.
#' @param time_period The number of time periods in the dataset.
#' @param untreated_prop The proportion of untreated units.
#' @param epsilon_size The standard deviation for the error term in potential outcomes.
#' @param cov The type of covariate to include ("no", "int", or "cont").
#' @param hetero The type of heterogeneity in treatment effects ("all" or "dynamic").
#' @param second_outcome Whether to include a second outcome variable.
#' @param second_cov Whether to include a second covariate.
#' @param na Whether to generate missing data ("none", "y", "x", or "both").
#' @param balanced Whether to balance the dataset by random sampling.
#' @param seed Seed for random number generation.
#' @param stratify Whether to stratify the dataset based on a binary covariate.
#' @param treatment_assign The method for treatment assignment ("latent" or "uniform").
#' @param vary_cov include time-varying covariates
#' @param second_cohort include confounding events
#' @param confound_ratio extent of event confoundedness
#' @param second_het heterogeneity of the second event
#'
#' @return A list containing the simulated dataset (dt) and the treatment effect values (att).
#'
#' @examples
#' # Simulate a DiD dataset with default settings
#' data <- sim_did(sample_size = 100, time_period = 5)
#'
#' # Simulate a DiD dataset with customized settings
#' data <- sim_did(sample_size = 200, time_period = 8, cov = "int", hetero = "dynamic")
#'
#' @export
sim_did <- function(sample_size, time_period, untreated_prop = 0.3, epsilon_size = 0.001,
                    cov = "no", hetero = "all", second_outcome = FALSE, second_cov = FALSE, vary_cov = FALSE, na = "none", 
                    balanced = TRUE, seed = NA, stratify = FALSE, treatment_assign = "latent", second_cohort = FALSE, confound_ratio = 1, second_het = "all"){
  
  if(!is.na(seed)){set.seed(seed)}
  
  #unit  -------------
  dt_i <- data.table(unit = 1:sample_size)
  if(cov == "int"){
    dt_i[, x := sample.int(5, sample_size, replace = TRUE)] 
  } else if (cov == "no"){
    dt_i[, x := 1] 
  } else if (cov == "cont"){
    dt_i[, x := rnorm(sample_size)]
  }
  
  if(second_cov){
    dt_i[, x2 := rnorm(sample_size)]
  } else {dt_i[, x2 := 1]}
  
  if(stratify){
    dt_i[, s := fifelse(rnorm(sample_size) > 0, 1, 2)]
  } else {
    dt_i[, s := 1]
  }

  #treatment assignment ---------------------
  
  #assign treated group based on a latent related to X
  
  if(treatment_assign == "latent"){
    ep1 <- rnorm(sample_size)
    dt_i[, treat_latent := x*0.2 + x2*0.2 + ep1] #unit with larger X tend to be treated and treated earlier
    untreated_thres <- quantile(dt_i$treat_latent, untreated_prop)
    dt_i[treat_latent <= untreated_thres, G := Inf] #unit with low latent is never treated
    
    cohort_prop <- (1-untreated_prop)/(time_period-1)
    last_treat_thres <- untreated_thres
    for(t in time_period:2){ #unit start getting treated in t = 2
      treat_thres <- quantile(dt_i$treat_latent, untreated_prop + cohort_prop*(time_period - t + 1))
      dt_i[treat_latent <= treat_thres & treat_latent > last_treat_thres, G := t]
      last_treat_thres <- treat_thres
    }
    rm(t)
  } else if (treatment_assign == "uniform"){
    #when treatment is set to 'uniform', untreated propensity is fixed
    dt_i[,G := floor((unit-1)/(sample_size/time_period))]
    dt_i[G < 2, G := Inf]
  }
  
  if(second_cohort){
    setnames(dt_i, "G", "G2")
    if(treatment_assign == "latent"){
      dt_i[, treat_latent := x*0.2 + x2*0.2 + ep1*confound_ratio + rnorm(sample_size)] #unit with larger X tend to be treated and treated earlier
      untreated_thres <- quantile(dt_i$treat_latent, untreated_prop)
      dt_i[treat_latent <= untreated_thres, G := Inf] #unit with low latent is never treated
      
      cohort_prop <- (1-untreated_prop)/(time_period-1)
      last_treat_thres <- untreated_thres
      for(t in time_period:2){ #unit start getting treated in t = 2
        treat_thres <- quantile(dt_i$treat_latent, untreated_prop + cohort_prop*(time_period - t + 1))
        dt_i[treat_latent <= treat_thres & treat_latent > last_treat_thres, G := t]
        last_treat_thres <- treat_thres
      }
      rm(t)
    } else if (treatment_assign == "uniform"){
      #when treatment is set to 'uniform', untreated propensity is fixed
      dt_i[,G := floor((unit-1)/(sample_size/time_period))]
      dt_i[G < 2, G := Inf]
    }
  }

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
  
  #add time_varying covariates
  if(vary_cov){
    dt[, xvar := pmin(G, time_period+4)*time^(1/3)*0.1+rnorm(sample_size*time_period, 0,10)] #should be confounding....?
  } else {
    dt[, xvar := 1]
  }
  
  #untreated potential outcomes
  dt[, y0 := unit_fe + time_fe + (x+x2)*x_trend + xvar + rnorm(sample_size*time_period, sd = epsilon_size)]
  
  dt[, D := as.integer(time >= G)]
  #generate gtatt
  att <- CJ(G = 1:time_period, time = 1:time_period)
  if(hetero == "all"){
    att[, attgt := rnorm(time_period*time_period, mean = 2, sd = 1)]
  } else if (hetero == "dynamic"){
    for(event_t in 0:max(att[,time-G])){
      att[time - G == event_t, attgt := rnorm(1, mean = 2, sd = 1)]
    }
  }
  att[time < G, attgt := 0] #no anticipation
  dt <- dt |> merge(att, by = c("time", "G"), all.x = TRUE, all.y = FALSE)
  dt[is.na(attgt), attgt := 0]
  dt[, tau := attgt*s]
  
  if(second_cohort){
    dt[, D2 := as.integer(time >= G2)]
    #generate gtatt
    att2 <- CJ(G2 = 1:time_period, time = 1:time_period)
    if(second_het == "no"){
      att2[, attgt2 := 10]
    } else {
      if(hetero == "all"){
        att2[, attgt2 := rnorm(time_period*time_period, mean = 2, sd = 1)]
      } else if (hetero == "dynamic"){
        for(event_t in 0:max(att2[,time-G2])){
          att2[time - G2 == event_t, attgt2 := rnorm(1, mean = 2, sd = 1)]
        }
      }
    }
    att2[time < G2, attgt2 := 0] #no anticipation
    
    #add att2 to att
    att[, event := 1]
    att2[, event := 2]
    att <- rbind(att, att2, fill = TRUE)
    
    dt <- dt |> merge(att2, by = c("time", "G2"), all.x = TRUE, all.y = FALSE)
    dt[is.na(attgt2), attgt2 := 0]
    dt[, tau2 := attgt2*s]
    
    #potential outcome
    dt[, y10 := y0 + tau]
    dt[, y01 := y0 + tau2]
    dt[, y11 := y0 + tau + tau2]
    dt[, y := y0*(1-D)*(1-D2)+ y10*D*(1-D2) + y01*(1-D)*D2 + y11*D*D2]
    cols <- c("time", "G", "G2", "unit", "x", "x2", "y", "s", "xvar")
    
  } else {
    #potential outcome
    dt[, y1 := y0 + tau]
    dt[, y := y1*D + y0*(1-D)]
    cols <- c("time", "G", "unit", "x", "x2", "y", "s", "xvar")
  }
  
  dt <- dt[, .SD, .SDcols = cols]
  
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
    dt <- dt[sample(1:size, size*0.99)]
  }
  if(second_outcome == TRUE){  dt[, y2 := y + 1 + rnorm(fnrow(dt), 0, 0.1)]}
  
  return(list(dt = dt, att = att))
  
}




