# estimation --------------------------------------------

event_ATTs<-function(eventdata,
                     outcomes,#vector of variable names
                     covariates=NULL, #vector of covariates on which to check balance
                     clustervar="id", 
                     weights="pweight",
                     keep_trends=TRUE,
                     se="cluster",
                     lean=TRUE,
                     ssc=NULL){
  
  #These regressions should work identically if the fixed effects (after the "|") were replaced with:
  # interaction(time_pair,id,cohort)
  if(se != "cluster") clustervar<-NULL
  if(is.null(ssc)) ssc<-ssc()
  
  construct_event_variables(eventdata)
  eventdata$stratify <- factor(eventdata$stratify, levels = c(levels(eventdata$stratify), "empty"))
  
  results_dynamic <- event_ATTs_dynamic(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  results_pooled <- event_ATTspooled(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  results_means <- event_ATTs_means(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  results_covariates <- event_ATTs_covariates(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  
  return(list(pooled = results_pooled,
              dynamic = results,
              means = results_means,
              covariates = results_covariates))
}


# pre-process ------------------------------------------

process_stratify_ipw <- function(eventdata, stratify_balance_val){
  
  if(!is.na(stratify_balance_val) & stratify_balance_val != "mean"){
    eventdata[,treated_base := treated == 1 & stratify == stratify_balance_val]
    
    stratvals<-unique(eventdata$stratify)
    for(strat in stratvals){
      if(strat != stratify_balance_val){
        stratbalmodel <- feols(treated_base ~ 1 | interaction(cohort,event_time,time_pair,balancevars, drop = TRUE),
                               data = eventdata[(treated == 1 & stratify == stratify_balance_val) | 
                                                  (treated == 1 & stratify == strat),], lean = TRUE, mem.clean=TRUE)
        eventdata[, eval(paste0("pval_1",strat)) := predict(stratbalmodel,
                                                            eventdata)]
        
        eventdata[treated == 1 & stratify == strat, pweight_stratbal := get(paste0("pval_1",strat))/(1-get(paste0("pval_1",strat)))]
      }
      else    eventdata[treated == 1 & stratify == strat, pweight_stratbal := 1]
      
      stratbalmodel <- feols(treated_base ~ 1 | interaction(cohort,event_time,time_pair,balancevars, drop = TRUE),
                             data = eventdata[(treated == 1 & stratify == stratify_balance_val) | 
                                                (treated == 0 & stratify == strat),], lean = TRUE, mem.clean=TRUE)
      
      eventdata[, eval(paste0("pval_0",strat)) := predict(stratbalmodel,
                                                          eventdata)]
      eventdata[treated == 0 & stratify == strat, pweight_stratbal := get(paste0("pval_0",strat))/(1-get(paste0("pval_0",strat)))]
      
      
    }
    
    checkvars <- names(eventdata)[grepl("pval_",names(eventdata))]
    for(var in checkvars){
      eventdata[get(var) == 1 | get(var) == 0 | is.na(get(var)),  pweight_stratbal := NA ]
    }
  }
  if(!is.na(stratify_balance_val) & stratify_balance_val == "mean"){
    
    stratvals<-unique(eventdata$stratify)
    for(strat in stratvals){
      eventdata[,treated_base := treated == 1 & stratify == strat]
      stratbalmodel <- feols(treated_base ~ 1 | interaction(cohort,event_time,time_pair,balancevars, drop = TRUE),
                             data = eventdata, lean = TRUE, mem.clean=TRUE)
      eventdata[, eval(paste0("pval_1",strat)) := predict(stratbalmodel,
                                                          eventdata)]
      
      eventdata[treated == 1 & stratify == strat, pweight_stratbal :=1/(get(paste0("pval_1",strat)))]
      
      eventdata[,treated_base := treated == 0 & stratify == strat]
      stratbalmodel <- feols(treated_base ~ 1 | interaction(cohort,event_time,time_pair,balancevars, drop = TRUE),
                             data = eventdata, lean = TRUE, mem.clean=TRUE)
      
      eventdata[, eval(paste0("pval_0",strat)) := predict(stratbalmodel,
                                                          eventdata)]
      eventdata[treated == 0 & stratify == strat, pweight_stratbal := 1/(get(paste0("pval_0",strat)))]
      
      
      
    }
    #Need to rebalance everyone to look like distribution of TREATED guys:
    stratbalmodel <- feols(treated ~ 1 | interaction(cohort,event_time,time_pair,balancevars, drop = TRUE),
                           data = eventdata, lean = TRUE, mem.clean=TRUE)
    eventdata[, pval_treat := predict(stratbalmodel,
                                      eventdata)]
    eventdata[,pweight_stratbal:=pweight_stratbal*pval_treat]
    
    
    checkvars <- names(eventdata)[grepl("pval_",names(eventdata))]
    for(var in checkvars){
      eventdata[get(var) == 1 | get(var) == 0 | is.na(get(var)),  pweight_stratbal := NA ]
    }
  }
  
  return(eventdata)
}


process_iv <- function(eventdata, instrument, instrument_exposure, covs_instrument_base_balance, saturate){
  eventdata[,instrument_group_now:=max(get(instrument) == 1 & as.numeric(as.character(time_pair)) == as.numeric(as.character(event_time))),by=.(id,cohort,time_pair)]
  eventdata[,instrument_group_base:=max(get(instrument) == 1 & as.numeric(as.character(base_time)) == as.numeric(as.character(event_time))),by=.(id,cohort,time_pair)]
  if(instrument_exposure=="full") eventdata <- eventdata[instrument_group_now == instrument_group_base,]
  if(instrument_exposure=="partial") eventdata <- eventdata[(instrument_group_now == 0 & instrument_group_base == 0) | (instrument_group_base == 0 & instrument_group_now == 1),]
  if(instrument_exposure=="all") {
    eventdata <- eventdata[!is.na(instrument_group_now) & !is.na(instrument_group_base),]
    eventdata[,instrument_group_now := pmax(instrument_group_now,instrument_group_base)]
  }
  if(instrument_exposure=="base") {
    eventdata <- eventdata[!is.na(instrument_group_base),]
    eventdata[,instrument_group_now := instrument_group_base]
  }
  #After restricting on instrument status, need to re-impose balance:
  eventdata[,obscount:=1]
  eventdata[,obscount:=sum(ifelse(event_time==base_time,0,obscount)),by=.(id,cohort)]
  if(balanced_panel==TRUE) {
    eventdata <- eventdata[obscount==numperiods-1,]
  }
  
  eventdata[,instrument_group:=instrument_group_now]
  eventdata[,instrument_group_now:=NULL]
  eventdata[,instrument_group_base:=NULL]
  if(length(unique(eventdata$instrument_group)) == 1) 
    stop("No variation in instrument left in the panel. Maybe because you're trying to impose panel balance?")
  
  
  if(!is.na(covs_instrument_base_balance)) eventdata[,instrument_balancevars:=interaction(eventdata[,covs_instrument_base_balance,with=FALSE], drop=TRUE)]
  else eventdata[,instrument_balancevars:=factor(1,levels=c(1,"OMIT"))]
  
  eventdata[,ivgroup:=interaction(treated,instrument_group)]
  ivlevels<-unique(eventdata$ivgroup)
  for(strat in ivlevels){
    eventdata[,treated_base := ivgroup == strat]
    stratbalmodel <- feols(treated_base ~ 1 | interaction(cohort,event_time,time_pair,instrument_balancevars, drop = TRUE),
                           data = eventdata, lean = TRUE, mem.clean=TRUE,
                           weights = eventdata$pweight)
    eventdata[, eval(paste0("pval_",strat)) := predict(stratbalmodel,
                                                       eventdata)]
    
    eventdata[ivgroup == strat, pweight_new :=pweight/(get(paste0("pval_",strat)))]
    
  }
  #Need to rebalance everyone to look like distribution of TREATED guys:
  stratbalmodel <- feols(treated ~ 1 | interaction(cohort,event_time,time_pair,instrument_balancevars, drop = TRUE),
                         data = eventdata, lean = TRUE, mem.clean=TRUE,
                         weights = eventdata$pweight)
  eventdata[, pval_treat := predict(stratbalmodel,
                                    eventdata)]
  
  eventdata[,pweight_new:=pweight_new*pval_treat]
  eventdata[,pval_treat:=NULL]
  eventdata[,pweight:=pweight_new]
  eventdata[,pweight_new:=NULL]
  
  checkvars <- names(eventdata)[grepl("pval_",names(eventdata))]
  for(var in checkvars){
    eventdata[get(var) == 1 | get(var) == 0 | is.na(get(var)),  pweight := NA ]
    eventdata[,eval(var):=NULL]
  }
  eventdata<-eventdata[!is.na(pweight),]
  if(length(unique(eventdata$instrument_group)) == 1) stop("No variation in instrument left in the panel after imposing common support on observables
                                                             across instrument X treatment cells.")
  return(eventdata)
  
}

construct_event_variables_iv <- function(eventdata,saturate=FALSE,response=NULL){
  
  eventdata[,unitfe := finteraction(time_pair,id,treated,cohort,stratify)]
  
  eventdata[,Z:=instrument_group[event_time != base_time],by=.(unitfe)]
  eventdata[,R:=get(response)[event_time != base_time],by=.(unitfe)]
  eventdata[,response_event_time_stratify := interaction(event_time,stratify, drop = TRUE)]
  eventdata[,instrument_event_time_stratify := interaction(event_time,stratify, drop = TRUE)]
  #Omitting base year for all levels of --stratify--, for response interaction
  eventdata[R == 0 ,response_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  eventdata[event_time==base_time,response_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  eventdata[,response_event_time_stratify:=relevel(response_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  #Omitting base year for all levels of --stratify--, for instrument interaction
  eventdata[Z == 0 ,instrument_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  eventdata[event_time==base_time,instrument_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  eventdata[,instrument_event_time_stratify:=relevel(instrument_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  #Defining response-event-time interactions, setting to omitted value if response == 0:
  eventdata[,treated_response_event_time_stratify := treated_event_time_stratify]
  eventdata[R==0,treated_response_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  
  #Defining pooled response interactions, setting to omitted value if response == 0:
  eventdata[,treated_response_post_stratify := treated_post_stratify]
  eventdata[R==0,treated_response_post_stratify := paste0(c(0,1),collapse=".")]
  eventdata[,treated_response_pre_stratify := treated_pre_stratify]
  eventdata[R==0,treated_response_pre_stratify := 1]
  
  #Defining instrument-event-time interactions, setting to omitted value if instrument == 0:
  eventdata[,treated_instrument_event_time_stratify := treated_event_time_stratify]
  eventdata[Z==0,treated_instrument_event_time_stratify := paste0(c(max(eventdata$base_time),1),collapse=".")]
  
  
  #Defining pooled instrument interactions, setting to omitted value if instrument == 0:
  eventdata[,treated_instrument_post_stratify := treated_post_stratify]
  eventdata[Z==0,treated_instrument_post_stratify := paste0(c(0,1),collapse=".")]
  eventdata[,treated_instrument_pre_stratify := treated_pre_stratify]
  eventdata[Z==0,treated_instrument_pre_stratify := 1]
  
  if(saturate==TRUE) {
    eventdata[,event_time_stratify:=interaction(event_time_stratify,balancevars, drop=TRUE)]
    eventdata[,instrument_event_time_stratify:=interaction(instrument_event_time_stratify,balancevars, drop=TRUE)]
    eventdata[,treated_event_time_stratify:=interaction(treated_event_time_stratify,balancevars, drop=TRUE)]
    eventdata[,treated_pre_stratify:=interaction(treated_pre_stratify,balancevars, drop=TRUE)]
    eventdata[,treated_post_stratify:=interaction(treated_post_stratify,balancevars, drop=TRUE)]
    # if(fullsaturate==TRUE){
    # eventdata[,treated_instrument_event_time_stratify:=interaction(treated_instrument_event_time_stratify,balancevars, drop=TRUE)]
    # eventdata[,treated_instrument_event_time_stratify:=relevel(treated_instrument_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1,0),collapse="."))]
    # eventdata[event_time==base_time | Z == 0 | treated==0,treated_instrument_event_time_stratify:=levels(treated_instrument_event_time_stratify)[1]]
    # eventdata[,treated_instrument_post_stratify:=interaction(treated_instrument_post_stratify,balancevars, drop=TRUE)]
    # eventdata[,treated_instrument_pre_stratify:=interaction(treated_instrument_pre_stratify,balancevars, drop=TRUE)]
    # }
  }
  
  return(eventdata)
  
}



# estimation ----------------------------------

event_IVs<-function(eventdata,
                    outcomes,#vector of variable names
                    response,#binary, endogenous response variable
                    #                     instrument,#instrument we use to estimate causal effect of the response variable
                    saturate=FALSE,
                    #                    fullsaturate=FALSE,
                    clustervar="id", 
                    weights="pweight",
                    keep_trends=TRUE,
                    lean = TRUE,
                    ssc = NULL){
  if(is.null(ssc)) ssc<-ssc()
  
  construct_event_variables(eventdata,saturate=saturate,IV=TRUE,response=response)
  
  eventdata$stratify <- factor(eventdata$stratify, levels = c(levels(eventdata$stratify), "empty"))
  
  if(keep_trends == TRUE){
    results<-feols(as.formula(paste0("c(",
                                     paste0(outcomes,collapse=",") , ") ~ event_time_stratify + instrument_event_time_stratify + treated_event_time_stratify | unitfe |  treated_response_event_time_stratify  ~ treated_instrument_event_time_stratify"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
    
    
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ event_time_stratify + instrument_event_time_stratify + treated_pre_stratify + treated_post_stratify  | unitfe | treated_response_pre_stratify  + treated_response_post_stratify ~ treated_instrument_pre_stratify + treated_instrument_post_stratify"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  }else{
    results<-feols(as.formula(paste0("c(",
                                     paste0(outcomes,collapse=",") , ") ~ treated_event_time_stratify | event_time_stratify + instrument_event_time_stratify + unitfe | treated_response_event_time_stratify ~  treated_instrument_event_time_stratify "
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
    
    
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ treated_pre_stratify + treated_post_stratify | event_time_stratify + instrument_event_time_stratify + unitfe |  treated_response_pre_stratify +  treated_response_post_stratify ~ treated_instrument_pre_stratify + treated_instrument_post_stratify"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE
    )
  }
  
  
  results_means<-feols(as.formula(paste0("c(",
                                         paste0(outcomes,collapse=",") , 
                                         ") ~ interaction(stratify) - 1"
  )),
  data = eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),],
  weights= eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),get(weights)],
  cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  
  return(list(pooled = results_pooled,
              dynamic = results,
              means = results_means))
}

se_adjusted<-function(semats,eventdata,pw="pweight"){
  construct_event_variables(eventdata)
  baset<-max(eventdata$base_time)
  
  V<-model.matrix( ~event_time_stratify -1 
                   ,data =eventdata[event_time!=base_time])
  V<-as.matrix(V[,!grepl(paste0("event_time_stratify",baset),colnames(V))])
  V2<-V
  V2[eventdata[event_time!=base_time,treated==0],]<-0
  V<-cbind(V,V2)
  V<-V*eventdata[event_time!=base_time,get(pw)]
  
  
  D<-diag(1/eventdata[event_time!=base_time,2*get(pw)])
  
  correctingmat<-t(V)%*%D%*%V
  ZZ_inv<-solve(semats$ZZ-correctingmat)
  
  cov<-list()
  for(outcome in names(semats$sandwiches)){
    cov<-c(cov,list(ZZ_inv%*%semats$sandwich[[outcome]]%*%ZZ_inv))
  }
  names(cov)<-names(semats$sandwiches)
  return(cov)
}

se_iv_adjusted<-function(semats,eventdata,response,pw="pweight"){
  construct_event_variables(eventdata,saturate=FALSE,IV=TRUE,response=response)
  eventdata[,Ztreated:=Z*as.numeric(as.character(treated))]
  eventdata[,Rtreated:=R*as.numeric(as.character(treated))]
  baset<-max(eventdata$base_time)
  
  V <- model.matrix( ~event_time_stratify -1 
                     ,data =eventdata[event_time!=base_time])
  V<-as.matrix(V[,!grepl(paste0("event_time_stratify",baset),colnames(V))])
  V<-V*eventdata[event_time!=base_time,get(pw)]
  
  VZ <-V
  VX <-V
  for(vb in c("Z","treated","Ztreated")){
    V2<-V
    V2[eventdata[event_time!=base_time,get(vb)==0],]<-0
    VZ<-cbind(VZ,V2)
  }
  for(vb in c("Z","treated","Rtreated")){
    V2<-V
    V2[eventdata[event_time!=base_time,get(vb)==0],]<-0
    VX<-cbind(VX,V2)
  }
  
  D<-diag(1/eventdata[event_time!=base_time,2*get(pw)])
  
  correctingmatZ<-t(VZ)%*%D%*%VZ
  correctingmatX<-t(VZ)%*%D%*%VX
  ZZ<-semats$ZZ-correctingmatZ
  ZX<-semats$ZX-correctingmatX
  ZZ_inv <- solve(semats$ZZ-correctingmatZ)
  ZX_inv <- solve(semats$ZX-correctingmatX)
  
  #  cov<-list()
  #  for(outcome in names(semats$sandwiches)){
  #  cov<-c(cov,list(ZX_inv%*%semats$sandwiches[[outcome]]%*%t(ZX_inv)))
  #  }
  #  names(cov)<-names(semats$sandwiches)
  cov<-ZX_inv%*%semats$sandwich%*%t(ZX_inv)
  return(cov)
  #  ses<- solve(t(ZZ_inv%*%ZX)%*%ZX)%*%t(ZZ_inv%*%ZX)%*%semats$sandwich%*%ZZ_inv%*%ZX%*%solve(t(ZZ_inv%*%ZX)%*%ZX)
}

event_ATT_ses<-function(outcomes,
                        results,
                        eventdata,
                        includeFEs=FALSE,
                        clustervar="id", 
                        weights="pweight"){
  #Define useful variables:
  construct_event_variables(eventdata)
  
  
  if(includeFEs==TRUE){
    Z<-model.matrix( ~event_time_stratify + treated_event_time_stratify + unitfe
                     ,data =eventdata)
    ndof<-0
  }
  else{
    Z<-model.matrix( ~event_time_stratify + treated_event_time_stratify
                     ,data =eventdata)
    ndof<-nrow(eventdata)/2
  }
  Z<-Z[,which(apply(Z,MARGIN =2,FUN=var)!=0)]
  ZZ <- t(Z)%*%(Z*eventdata[,get(weights)])
  
  sandwiches<-list()
  for(outcome in outcomes){
    if(length(outcomes==1)) result <- results$dynamic
    if(length(outcomes> 1)) {
      outpos<-which(names(results$dynamic)==outcome)
      result <- results$dynamic[[outpos]]
    }
    eventdata[,e:=get(outcome) - predict(result,eventdata)]
    Ze <- Z*eventdata[,e*get(weights)]
    Ze<-as.data.table(Ze)
    Ze[,id:=eventdata[,get(clustervar)]]
    Ze<-as.matrix(Ze[,lapply(.SD,sum),by=id][,id:=NULL])
    
    
    sandwich<-t(Ze)%*%Ze
    sandwiches<-c(sandwiches,list(sandwich))
  }
  dof<-c(obs=nrow(eventdata),
         cl=length(unique(eventdata[,get(clustervar)]))
  )  
  #  dof<- (nrow(eventdata)-1)/(nrow(eventdata)-ncol(Ze)+1-ndof)
  #  dof.cl<-length(unique(eventdata[,get(clustervar)]))/(length(unique(eventdata[,get(clustervar)]))-1)
  names(sandwiches)<-outcomes
  
  return(list(ZZ=ZZ,
              sandwiches=sandwiches,
              dof=dof)
  )
}

event_IV_ses<-function(outcome,
                       response,
                       results,
                       eventdata,
                       includeFEs=FALSE,
                       clustervar="id", 
                       weights="pweight"){
  
  #Define useful variables:
  construct_event_variables(eventdata,saturate=FALSE,IV=TRUE,response=response)
  
  
  if(includeFEs==TRUE) {
    Z<-model.matrix( ~event_time_stratify + instrument_event_time_stratify + treated_event_time_stratify  + treated_instrument_event_time_stratify + unitfe
                     ,data =eventdata)
    
    X<-model.matrix( ~event_time_stratify + instrument_event_time_stratify + treated_event_time_stratify  + treated_response_event_time_stratify + unitfe
                     ,data =eventdata)
    ndof<-0
  }
  else{
    Z<-model.matrix( ~event_time_stratify + instrument_event_time_stratify + treated_event_time_stratify  + treated_instrument_event_time_stratify
                     ,data =eventdata)
    X<-model.matrix( ~event_time_stratify + instrument_event_time_stratify + treated_event_time_stratify  + treated_response_event_time_stratify 
                     ,data =eventdata)
    
    ndof<-nrow(eventdata)/2
  }
  Z<-Z[,which(apply(Z,MARGIN =2,FUN=var)!=0)]
  X<-X[,which(apply(X,MARGIN =2,FUN=var)!=0)]
  ZZ <- t(Z)%*%(Z*eventdata[,get(weights)])
  ZX <- t(Z)%*%(X*eventdata[,get(weights)])
  
  varnames<-rownames(results)
  #correcting name on the IV coefficient(s):
  varnames<-gsub("fit_","",varnames)
  if(!all(colnames(X)%in%varnames)) stop("The variables we need to predict a person's outcome (columns of matrix X) do not all appear in the given coeftable 'results'")
  
  #sandwiches<-list()
  #eventdata[,e:=get(outcome) - predict(result,eventdata)]
  e<-eventdata[,get(outcome)-mean(get(outcome)),by=unitfe]$V1
  Xres<-as.data.table(X)
  Xres[,unitfe:=eventdata[,unitfe]]
  Xres<-Xres[,.SD-lapply(.SD,mean),by=unitfe]
  Xres[,unitfe:=NULL]
  Xres<-as.matrix(Xres)
  e<-e - apply(Xres[,varnames], MARGIN=1,FUN=function(x) sum(x*results[,"Estimate"]))
  Xres<-NULL
  Ze <- Z*e*eventdata[,get(weights)]
  Ze<-as.data.table(Ze)
  Ze[,id:=eventdata[,get(clustervar)]]
  Ze<-as.matrix(Ze[,lapply(.SD,sum),by=id][,id:=NULL])
  
  
  sandwich<-t(Ze)%*%Ze
  # sandwiches<-c(sandwiches,list(sandwich))
  
  dof<-c(obs=nrow(eventdata),
         cl=length(unique(eventdata[,get(clustervar)]))
  )  
  #  dof<- (nrow(eventdata)-1)/(nrow(eventdata)-ncol(Ze)+1-ndof)
  #  dof.cl<-length(unique(eventdata[,get(clustervar)]))/(length(unique(eventdata[,get(clustervar)]))-1)
  #names(sandwiches)<-outcome
  
  return(list(ZZ=ZZ,
              ZX=ZX,
              sandwich=sandwich,
              dof=dof)
  )
}