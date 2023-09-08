
# estimation ---------------------------------------

event_ATTs_dynamic<-function(eventdata,
                             outcomes,#vector of variable names
                             clustervar="id", 
                             weights="pweight",
                             keep_trends=TRUE,
                             se="cluster",
                             lean=TRUE,
                             ssc=NULL,
                             mem.clean = TRUE){
  
  #eventdata[,event_time_stratify:=as.factor(as.character(event_time_stratify))]
  #eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  #eventdata[,treated_event_time_stratify:=as.factor(as.character(treated_event_time_stratify))]
  #eventdata[,treated_event_time_stratify:=relevel(treated_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]

  #construct the call
  outcomes_call <- paste0("c(", paste0(outcomes,collapse=","), ")")
  event_stratify_call <- ifelse(keep_trends, "~ event_time_stratify + treated_event_time_stratify |", 
                                "~ treated_event_time_stratify | event_time_stratify +")
  call <- paste0(outcomes_call, event_stratify_call, " unitfe" )

  results<-feols(as.formula(call),
                 data = eventdata,
                 weights= eventdata[,get(weights)],
                 cluster=clustervar, lean = TRUE, mem.clean = mem.clean)
  
  return(list(dynamic = results))
  
}


event_ATTs_means<-function(eventdata,
                           outcomes,#vector of variable names
                           clustervar="id", 
                           weights="pweight",
                           keep_trends=TRUE,
                           se="cluster",
                           lean=TRUE,
                           ssc=NULL){
  
  eventdata$stratify <- factor(eventdata$stratify, levels = c(levels(eventdata$stratify), "empty"))
  
  results_means<-feols(as.formula(paste0("c(",
                                         paste0(outcomes,collapse=",") , 
                                         ") ~ interaction(stratify) - 1"
  )),
  data = eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),],
  weights= eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),get(weights)],
  cluster=clustervar, lean = TRUE, mem.clean = TRUE)
  
  return(list(means = results_means))
}


event_ATTs_pooled<-function(eventdata,
                            outcomes,#vector of variable names
                            clustervar="id", 
                            weights="pweight",
                            keep_trends=TRUE,
                            se="cluster",
                            lean=TRUE,
                            ssc=NULL){
  
  eventdata[,treated_post := qF((treated == 1) * (post == 1))]
  eventdata[,treated_pre := qF((treated == 1) * (post == 0))]
  
  eventdata[,treated_post_stratify := finteraction(treated_post,stratify)]
  
  eventdata[,treated_pre_stratify := finteraction(treated_pre,stratify)]
  eventdata[treated_pre==0,treated_pre_stratify := 1]
  eventdata[event_time==base_time,treated_pre_stratify := 1]
  
  #Omitting effect for untreated people or observations in pre-period:
  eventdata[treated_post == 0 ,treated_post_stratify := paste0(c(0,1),collapse=".")]
  eventdata[,treated_post_stratify:=relevel(treated_post_stratify,ref = paste0(c(0,1),collapse="."))]
  
  eventdata[,event_time_stratify:=as.factor(as.character(event_time_stratify))]
  eventdata[,event_time_stratify:=relevel(event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  eventdata[,treated_event_time_stratify:=as.factor(as.character(treated_event_time_stratify))]
  eventdata[,treated_event_time_stratify:=relevel(treated_event_time_stratify,ref = paste0(c(max(eventdata$base_time),1),collapse="."))]
  eventdata[,treated_pre_stratify:=as.factor(as.character(treated_pre_stratify))]
  eventdata[,treated_post_stratify:=as.factor(as.character(treated_post_stratify))]
  eventdata[,treated_post_stratify:=relevel(treated_post_stratify,ref = paste0(c(0,1),collapse="."))]
  
  if(keep_trends == TRUE){
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ event_time_stratify + treated_pre_stratify + treated_post_stratify | unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = TRUE, mem.clean = TRUE)
  }
  else{
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ treated_pre_stratify + treated_post_stratify | event_time_stratify + unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    cluster=clustervar, lean = TRUE, mem.clean = TRUE
    )
  }
  return(list(pooled = results_pooled))
}


event_ATTs_covariates<-function(eventdata,
                                outcomes=NULL,#vector of variable names
                                covariates=NULL, #vector of covariates on which to check balance
                                clustervar="id", 
                                weights="pweight",
                                keep_trends=TRUE,
                                se="cluster",
                                lean=TRUE,
                                ssc=NULL){
  
  #eventdata$stratify <- factor(eventdata$stratify, levels = c(levels(eventdata$stratify), "empty"))
  eventdata[,int_event_time_stratify:=interaction(event_time,stratify)]
  
  if(!is.null(covariates)) results_covariates<-feols(as.formula(paste0("c(",
                                                                       paste0(covariates,collapse=",") , ") ~ int_event_time_stratify + int_event_time_stratify:treated- 1"
  )),
  data = eventdata[as.numeric(as.character(event_time)) == as.numeric(as.character(time_pair)) | as.numeric(as.character(time_pair)) == 0,],
  weights= eventdata[as.numeric(as.character(event_time)) == as.numeric(as.character(time_pair)) | as.numeric(as.character(time_pair)) == 0,get(weights)],
  se=se,
  cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  else results_covariates=NULL
  
  return(list(covariates = results_covariates))
}


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
  
  if(keep_trends == TRUE){
    results<-feols(as.formula(paste0("c(",
                                     paste0(outcomes,collapse=",") , ") ~ event_time_stratify + treated_event_time_stratify | unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    se=se,
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
    
    
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ event_time_stratify + treated_pre_stratify + treated_post_stratify | unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    se=se,
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  }else{
    results<-feols(as.formula(paste0("c(",
                                     paste0(outcomes,collapse=",") , ") ~ treated_event_time_stratify | event_time_stratify + unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    se=se,
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
    
    
    results_pooled<-feols(as.formula(paste0("c(",
                                            paste0(outcomes,collapse=",") , 
                                            ") ~ treated_pre_stratify + treated_post_stratify | event_time_stratify + unitfe"
    )),
    data = eventdata,
    weights= eventdata[,get(weights)],
    se=se,
    cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE
    )
  }
  
  results_means<-feols(as.formula(paste0("c(",
                                         paste0(outcomes,collapse=",") , 
                                         ") ~ interaction(stratify) - 1"
  )),
  data = eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),],
  weights= eventdata[treated == 1 & event_time == base_time & !is.na(get(weights)),get(weights)],
  se=se,
  cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  
  if(!is.null(covariates)) results_covariates<-feols(as.formula(paste0("c(",
                                                                       paste0(covariates,collapse=",") , ") ~ interaction(event_time,stratify) + interaction(event_time,stratify):treated- 1"
  )),
  data = eventdata[as.numeric(as.character(event_time)) == as.numeric(as.character(time_pair)) | as.numeric(as.character(time_pair)) == 0,],
  weights= eventdata[as.numeric(as.character(event_time)) == as.numeric(as.character(time_pair)) | as.numeric(as.character(time_pair)) == 0,get(weights)],
  se=se,
  cluster=clustervar, lean = lean, ssc=ssc, mem.clean = TRUE)
  else results_covariates=NULL
  
  return(list(pooled = results_pooled,
              dynamic = results,
              means = results_means,
              covariates = results_covariates))
}

# get result -------------------

get_result_dynamic<-function(eventdata_panel,variable,trends=TRUE, mem.clean = TRUE){
  
  if(nrow(eventdata_panel)==0){ 
    warning("event panel is empty, returning empty data.table")
    return(data.table())
  }
  
  if(any(eventdata_panel[,lapply(.SD, stats::var), .SDcols = variable] == 0)){
    warning("some outcome have no variation, returning empty data.table")
    return(data.table())
  }
    
  results<-event_ATTs_dynamic(eventdata_panel,outcomes = c(variable),keep_trends = trends, mem.clean = mem.clean)
  
  if(length(variable) > 1){
    dt <- data.table()
    for(i in 1:length(variable)){
      
      result <- results$dynamic[[1]]
      dt <- rbind(dt, data.table(variable = row.names(result$coeftable),result$coeftable,obs=result$nobs,outcome = variable[i]))
      
    }
  } else {dt<-data.table(variable = row.names(results$dynamic$coeftable),results$dynamic$coeftable,obs=results$dynamic$nobs)}

  dt[,result:="dynamic"]
  
  get_event_time <- function(x){
    start <- str_locate(x, "stratify")[2]
    end <- str_locate(x, "\\.")[1]
    return(as.numeric(str_sub(x, start + 1, end - 1)))
  }
  
  dt[,event_time:= lapply(variable, get_event_time)]

  return(dt)
  
}

get_result_pooled<-function(eventdata_panel,variable,trends=TRUE){
  
  if(nrow(eventdata_panel)==0){
    dt<-data.table()
    return(dt)
  }else if ((eventdata_panel[,var(get(variable))]==0)==T){
    dt<-data.table()
    return(dt)
  }else{
    
    results<-event_ATTs_pooled(eventdata_panel,outcomes = c(variable),keep_trends = trends)
    dt<-data.table(variable = row.names(results$pooled$coeftable),model=i,results$pooled$coeftable,obs=results$pooled$nobs)
    dt<-dt[,result:="pooled"]
    rm(eventdata_panel)
    rm(results)
    gc()
    return(dt)
    
  }
}

get_result_means<-function(eventdata_panel,variable,trends=TRUE){
  
  results<-event_ATTs_means(eventdata_panel,outcomes = c(variable),keep_trends = trends)
  dt<-data.table(variable = row.names(results$means$coeftable),model=i,results$means$coeftable,obs=results$means$nobs)
  dt<-dt[,result:="means"]
  
  return(dt)
  
}

get_result_covariates<-function(eventdata_panel,covariate,variable=NULL,trends=TRUE){
  
  results<-event_ATTs_covariates(eventdata_panel,covariates=covariate,outcomes = c(variable),keep_trends = trends)
  dt<-data.table(variable = row.names(results$covariates$coeftable),results$covariates$coeftable,obs=results$covariates$nobs)
  dt<-dt[,result:="covariates"][,model:=paste0("additional_covariate--",covariate)]
  
  return(dt)
  
}

# IV related ----------------------------------

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




