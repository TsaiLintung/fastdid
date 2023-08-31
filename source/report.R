# reporting ------------------------------

print_ATT<-function(results,
                    outcomes=NULL,
                    outcome_names=NULL,
                    event_name="Event",
                    pooled_tables=TRUE,
                    dynamic_plots=TRUE,
                    pooled_tables_name=NULL,
                    base_time = NULL,
                    dynamic_pdfname=NULL,
                    stratify_values = NULL,
                    stratify_names = NULL,
                    decimals=0,
                    plot_pval=0.05){
  if(!is.null(stratify_values) & is.null(stratify_names)) stop("Must provide a name for each entry of stratify_values in stratify_names")
  if(pooled_tables==TRUE){
    
    if(is.null(outcome_names)) outcome_names <- names(results$pooled)
    if(is.null(outcomes)) outcomes <- names(results$pooled)
    if(!is.null(stratify_values)){
      tab<-TR(c("","Pre-period Mean","Treatment Effect"),cspan=c(1,length(stratify_values),length(stratify_values)))
      tab<-tab + TR(c("Outcome",stratify_names,stratify_names))
    }
    else{
      tab<-TR(c("Outcome","Pre-period Mean","Treatment Effect"))
      stratify_values<-1
    }
    tab<-tab + midrulep(list(c(2,length(stratify_values)+1),c(length(stratify_values)+2,2*length(stratify_values)+1))) 
    p<-1
    for(outcome in outcomes){
      pos<-which(names(results$pooled)==outcome)
      tab<-tab + TR(outcome_names[p]) %:%
        #Adding in the pre-period means:
        TR(results$means[[pos]]$coefficients[paste0("interaction(stratify)",stratify_values)],
           pvalues = ifelse(is.na(results$means[[pos]]$coeftable[paste0("interaction(stratify)",stratify_values),"Pr(>|t|))"]),
                            1,
                            results$means[[pos]]$coeftable[paste0("interaction(stratify)",stratify_values),"Pr(>|t|))"]
           ), dec = decimals) %:%
        #Adding in the treatment effects:
        TR(results$pooled[[pos]]$coefficients[paste0("treated_post_stratify1.",stratify_values)],
           pvalues = ifelse(is.na(results$pooled[[pos]]$coeftable[paste0("treated_post_stratify1.",stratify_values),"Pr(>|t|))"]),
                            1,
                            results$pooled[[pos]]$coeftable[paste0("treated_post_stratify1.",stratify_values),"Pr(>|t|))"]
           ), dec = decimals)
      #STANDARD ERRORS:
      tab<- tab + TR("") %:%
        #For pre-period means
        TR(results$means[[pos]]$se[paste0("interaction(stratify)",stratify_values)],se=TRUE,dec = decimals)%:%
        #For treatment effects:
        TR(results$pooled[[pos]]$se[paste0("treated_post_stratify1.",stratify_values)],se=TRUE,dec = decimals)
      
      p<-p+1
    }
    TS(tab, file=paste0(paste(c(pooled_tables_name,"pooled_table"),collapse="_")),
       output_path=".",
       pretty_rules=T,
       header=c('r',rep('c',2*length(stratify_values))))
    
  }
  
  if(dynamic_plots==TRUE){
    p<-1
    for(outcome in outcomes){
      pos<-which(names(results$pooled)==outcome)
      
      plot_table<-data.table(
        coef = results$dynamic[[pos]]$coefficients[
          grepl("treated_event_time_stratify",names(results$dynamic[[pos]]$coefficients))
        ],
        se = results$dynamic[[pos]]$se[
          grepl("treated_event_time_stratify",names(results$dynamic[[pos]]$se))
        ],
        event_time = as.numeric(as.character(gsub("\\..*","", 
                                                  gsub("treated_event_time_stratify", "",
                                                       names(results$dynamic[[pos]]$coefficients)[grepl("treated_event_time_stratify",names(results$dynamic[[pos]]$se))])
        ))),
        stratify_value = as.numeric(as.character(gsub(".*\\.","", 
                                                      gsub("treated_event_time_stratify", "",
                                                           names(results$dynamic[[pos]]$coefficients)[grepl("treated_event_time_stratify",names(results$dynamic[[pos]]$se))])
        )))
      )
      #adding in omitted year:
      #Looking for the first missing period in the interval of plotted periods
      #If the interval of plotted periods has no gap, I assume refernce period is first
      #period in the data.
      if(is.null(base_time)){ 
        missingperiods<-(min(plot_table$event_time):max(plot_table$event_time))[!(min(plot_table$event_time):max(plot_table$event_time))%in%unique(plot_table$event_time)]
        if(length(missingperiods)==0) refperiod <- min(plot_table$event_time)-1
        else refperiod <- min(missingperiods)
      }
      else refperiod <- base_time
      plot_table<-rbind(plot_table, data.table(coef=0,
                                               se = 0,
                                               event_time = refperiod,
                                               stratify_value = unique(plot_table$stratify_value)
      ))
      
      plot_table[,upper:= coef + abs(qt(plot_pval/2,
                                        df = results$dynamic[[pos]]$nobs - results$dynamic[[pos]]$nparams))*se]
      plot_table[,lower:= coef - abs(qt(plot_pval/2,
                                        df = results$dynamic[[pos]]$nobs - results$dynamic[[pos]]$nparams))*se]
      pdf(paste0("./",paste(c(dynamic_pdfname,outcome),collapse="_"),".pdf"))
      for(stratval in stratify_values){
        print(ggplot(data=plot_table[stratify_value==stratval,]) + 
                geom_ribbon(aes(ymin = lower, ymax=upper, x=event_time), fill="grey50", alpha=0.5) + 
                geom_hline(yintercept = 0, linetype="dashed") + 
                geom_vline(xintercept = 0, linetype="dashed") + 
                geom_line(aes(x = event_time,
                              y = coef)) +
                scale_x_continuous( breaks = pretty_breaks(12)) +
                ylim(c(min(plot_table$lower),max(plot_table$upper))) +
                labs(y=varnames[p], x = paste0("Years from ",event_name))
        )
      }
      dev.off()
      p<-p+1
    }
  }
  
}

counter_eventdata <- function(eventdata){
  eventdata_counter<-eventdata[treated==1 & pweight ==1 ,.N, by=.(stratify,event_time)]
  return(eventdata_counter)
}


summary_eventdata1 <- function(eventdata,
                               summarylist,
                               summarylevel){
  
  eventdata<-eventdata[,stratify2:=paste(get(summarylevel),post,treatgroup,sep=".")]
  weights<-"pweight"
  summarytable1<-data.table()
  for(h in summarylist){
    
    results_means<-feols(as.formula(paste0(h, "~ interaction(stratify2) - 1")),
                         data = eventdata, cluster = "id", weights = eventdata[,get(weights)],
                         lean = TRUE, mem.clean = TRUE)
    tab<-data.table(variable = row.names(results_means$coeftable),results_means$coeftable,outcome=h,obs=results_means$nobs)
    summarytable1<-rbind(summarytable1,tab)
    console.log(h)
  }
  tab_N<-eventdata[,.N,by=.(stratify2)][,variable:=paste0("interaction(stratify2)",stratify2)]
  summarytable1<-tab_N[summarytable1,on="variable"]
  
  return(summarytable1)
  rm(results_means)
  rm(summarytable1)
  gc()
  
}


summary_eventdata2 <- function(eventdata,
                               summarylist,
                               summarylevel){
  
  eventdata<-eventdata[,stratify2:=paste(get(summarylevel),post,treatgroup,sep=".")]
  weights<-"pweight"
  summarytable2<-data.table()
  for(h in summarylist){
    tab1<-eventdata[,.(weighted.mean(get(h),pweight),wtd.var(get(h),pweight),mean(get(h)),var(get(h))),by=.(stratify2)]%>%rename(weighted_mean=V1,
                                                                                                                                 weighted_var=V2,
                                                                                                                                 mean=V3,
                                                                                                                                 var=V4)
    tab1<-tab1[,outcome:=h]
    summarytable2<-rbind(summarytable2,tab1,fill=T)
    console.log(h)
    gc()
    
  }
  tab_N<-eventdata[,.N,by=.(stratify2)][,variable:=paste0("interaction(stratify2)",stratify2)]
  summarytable2<-tab_N[summarytable2,on="stratify2"]
  return(summarytable2) 
  rm(list=ls(pattern="^tab"))
}

