#' Diagnose Event Confoundedness
#'
#' @param data A data.table containing the panel data.
#' @param timevar The name of the time variable.
#' @param cohortvar The name of the cohort (group) variable.
#' @param control_option The control units used for the DiD estimates. Default is "both".
#' @param weightvar The name of the weight variable, if not specified will cluster on unit level (optional).
#' @param cohortvar2 The name of the second cohort (group) variable.
#' 
#' @return A data.table containing the contamination ratio for each group-time.
#' @export
#'
#' @examples
#' # simulated data
#' simdt <- sim_did(1e+03, 10, cov = "cont", second_cov = TRUE, second_outcome = TRUE)
#'
#' @keywords difference-in-differences fast computation panel data estimation did
diagnose_confound_event <- function(data, timevar, cohortvar, cohortvar2, control_option = "both", weightvar = NA, anticipation = 0){
  
  if(is.na(weightvar)){
    dt <- data[, .(w = .N), by = c(timevar, cohortvar, cohortvar2)]
  } else {
    dt <- data[, .(w = sum(get(weightvar))), by = c(timevar, cohortvar2, cohortvar2)]
  }
  setnames(dt, c("time", "G", "G2", "w"))
  dt[, D2 := as.numeric(time >= G2)]
  
  expo <- data.table() 
  for(g in dt[, unique(G)]){
    for(t in dt[, unique(time)]){ 
      
      base <- g-1-anticipation
      dt[,`:=`(tp = 0, cp = 0, tb = 0, cb = 0)]
      dt[G == g & time == t, tp := w/sum(w)]
      dt[G > max(g,t) & time == t, cp := w/sum(w)]
      dt[G == g & time == base, tb := w/sum(w)]
      dt[G > max(g,t) & time == base, cb := w/sum(w)]
      if(any(dt[, lapply(.SD, sum), .SDcols = c("tp", "cp", "tb", "cb")] == 0)){next} #rule out invalid gt
      
      gamma <- dt[, sum(D2*(tp-cp-tb+cb))] #btw, this is how simple staggered DiD can be without conditional PT
      w = dt[G == g & time == t, sum(w)] #weight for the simple aggregation
      expo <- expo |> rbind(data.table(G = g, time = t, gamma = gamma, w = w))
      
    }
  }
  
  class(expo) <- c("confound_diagnosis", "data.table", "data.frame") #add class for generics
  return(expo)
  
}

#' @export
print.confound_diagnosis <- function(x,...) {
  cat(paste0("Average confoundedness: ", mean(x)))
}


#' @export
plot.confound_diagnosis <- function(expo) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("The ggplot2 package must be installed to use plotting functions")   #Either exit or do something without rgl
    return(NULL)
  }

  expo[time >= G] |> ggplot(aes(x = time-G, y = as.factor(G), fill = gamma)) + geom_tile() + scale_fill_distiller(palette = "RdBu", limits = c(-1, 1)) + 
    labs(subtitle = paste0("average confoundedness: ", mean(expo)), y = "G", x = "event_time") 
}

#' @export
mean.confound_diagnosis <- function(x){
  return(x[time>=G,weighted.mean(gamma)])
}