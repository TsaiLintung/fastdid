
#g11 <- c(1,1,2,3,4)
#g22 <- c(2,1,3,2,4)
#GG <- as.factor(paste0(g11, ".", g22))

g1 <- function(GG){
  if(is.numeric(GG)){return(GG)}
  return(as.numeric(str_split_i(as.character(GG), "-", 1)))
}
g2 <- function(GG){
  return(as.numeric(str_split_i(as.character(GG), "-", 2)))
}
ming <- function(GG){
  if(is.numeric(GG)){return(GG)}
  else {pmin(g1(GG), g2(GG))}
}


#the scheme for getting event-specific effect
get_es_scheme <- function(group_time, aux, p){
  
  es_group_time <- copy(group_time) #group_time with available es effect
  es_weight <- data.table()
  for(ggt in seq_len(nrow(group_time))){
    
    group_time <- get_es_ggt_weight(group_time, ggt, aux, p)
    
    if(group_time[, all(weight == 0)]){ #no available stuff
      t <- group_time[ggt, time]
      gg <- group_time[ggt, G]
      es_group_time <- es_group_time[!(time == t & G == gg)] #remove the ggt from new group time
    } else {
      es_ggt_weights <- group_time[, .(weight)] |> transpose()
      es_weight <- es_weight |> rbind(es_ggt_weights)
    }
    
  }
  
  return(list(group_time = es_group_time, es_weight = es_weight))
  
}

#get the scheme for retriving group-group-time estimates
get_es_ggt_weight <- function(group_time, ggt, aux, p){
  
  group_time[, weight := 0] #reset 
  t <- group_time[ggt, time]
  g1 <- group_time[ggt, g1(G)]
  g2 <- group_time[ggt, g2(G)]
  gg <- group_time[ggt, G]
  
  if(t < g2){ #direct pure effect
    
    group_time[ggt, weight := 1] #just use the observed effect
    
  } else if(g1 < g2) { #imputation = treat-pre + (control-post - control-pre)
    
    base_period <- g2 - 1
    if(base_period == t){return(group_time)}
    
    #get the cohorts
    tb <- group_time[,G == gg & time == base_period]
    c <-  group_time[,g1(G) == g1 & g2(G) > max(t,base_period) & g2(G) != g2]
    cp <- group_time[, c & time == t]
    cb <- group_time[, c & time == base_period]
    
    #if any group have no available cohort, skip
    if(sum(tb) == 0 | sum(cp) == 0 | sum(cb) == 0){return(group_time)}
    
    #assign the weights
    group_time[tb, weight := pg/sum(pg)]
    group_time[cp, weight := pg/sum(pg)]
    group_time[cb, weight := -pg/sum(pg)]
    
    
  } else if (g1 > g2) { #double did = (treat-post - treat-base) - (control-post - control-pre)
    
    base_period <- g1 - 1
    if(base_period == t){return(group_time)}
    
    #get the cohorts
    tp <- group_time[,.I == ggt]
    tb <- group_time[,G == gg & time == base_period]
    c <-  group_time[,g2(G) == g2 & g1(G) > max(t,base_period) & g1(G) != g1]
    cp <- group_time[, c & time == t]
    cb <- group_time[, c & time == base_period]
    
    #if any group have no available cohort, skip
    if(sum(tp) == 0 | sum(tb) == 0 | sum(cp) == 0 | sum(cb) == 0){return(group_time)}
    
    #assign the weights
    group_time[tp, weight := pg/sum(pg)]
    group_time[tb, weight := -pg/sum(pg)]
    group_time[cp, weight := -pg/sum(pg)]
    group_time[cb, weight := pg/sum(pg)]
    
  } 
  
  return(group_time)
  
}
