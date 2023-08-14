setwd("~/github/DiDforBigData")
library(stringr, warn.conflicts = F, quietly = T)
library(didimputation, warn.conflicts = F, quietly = T)
library(did, warn.conflicts = F, quietly = T)
library(DIDmultiplegt, warn.conflicts = F, quietly = T)
library(DiDforBigData, warn.conflicts = F, quietly = T)
library(profmem, warn.conflicts = F, quietly = T)
library(fixest, warn.conflicts = F, quietly = T)
library(parallel, warn.conflicts = F, quietly = T)

# didimputation
Borusyaketal <- function(inputdata, varnames, horizon=c(0,1,2,3)){
  BJS = did_imputation(data=inputdata, yname = varnames$outcome_name, gname=varnames$cohort_name,
                       tname = varnames$time_name, idname = varnames$id_name, horizon=horizon)
  BJS = as.data.table(BJS)[,list(EventTime=term, ATTe=estimate, ATTe_SE=std.error)]
  return(BJS)
}

# did
CallawaySantanna <- function(inputdata, varnames, est_method="reg", bstrap=TRUE){
  CS =  att_gt(yname = varnames$outcome_name,
               gname = varnames$cohort_name,
               idname = varnames$id_name,
               tname = varnames$time_name,
               xformla = ~1,
               data = inputdata,
               est_method = "reg",
               control_group = "notyettreated",
               base_period="universal",
               bstrap=bstrap,
               anticipation=0)
  CS_main = data.table(Cohort=CS$group, EventTime=CS$t-CS$group, CalendarTime=CS$t, ATTge=CS$att, ATTge_SE=CS$se)
  CS_ES = aggte(CS,type="dynamic")
  CS_ES = data.table(EventTime=CS_ES$egt, ATTe=CS_ES$att.egt, ATTe_SE=CS_ES$se.egt)
  return(list(results_cohort=CS, results_average=CS_ES))
}

# DIDmultiplegt
deChaisemartin <- function(inputdata,varnames,dynamic=3,placebo=0,brep=40){
  # set up variable names
  time_name = varnames$time_name
  outcome_name = varnames$outcome_name
  cohort_name = varnames$cohort_name
  id_name = varnames$id_name
  # estimation
  inputdata[, D := as.numeric(year >= get(cohort_name))]
  CH_res = did_multiplegt(df=copy(inputdata), Y=outcome_name, T = time_name, G = id_name, D = "D",
                          dynamic=dynamic, placebo=placebo, brep=brep, cluster=id_name)
  # collect results
  if(dynamic>0){
    CH_res_dyn = data.table( EventTime=0, ATTe=CH_res$effect, ATTe_SE=CH_res$se_effect)
    for(ee in 1:dynamic){
      CH_res_dyn = rbindlist(list(CH_res_dyn, data.table( EventTime=ee, ATTe=CH_res[[paste0("dynamic_",ee)]], ATTe_SE=CH_res[[paste0("se_dynamic_",ee)]])))
    }
    if(placebo>0){
      for(ee in 1:placebo){
        CH_res_dyn = rbindlist(list(CH_res_dyn, data.table( EventTime=(-1*ee), ATTe=CH_res[[paste0("placebo_",ee)]], ATTe_SE=CH_res[[paste0("se_placebo_",ee)]])))
      }
    }
    CH_res_dyn = CH_res_dyn[order(EventTime)]
    return(CH_res_dyn)
  }
}

# Test each DiD estimator's speed and memory usage
TestDiD <- function(estimator = "DiDforBigData", sample_sizes=c(1e3,1e4,1e5), reps=3, output_dir = "inst/speed_tests/"){

  # prepare variable names
  varnames = list()
  varnames$time_name = "year"
  varnames$outcome_name = "Y"
  varnames$cohort_name = "cohort"
  varnames$id_name = "id"

  # define speed tester
  speed_tester <- function(sample_size,reps){
    all_res = data.table()
    for(seed in 1:reps){
      this_res = data.table()
      # simulate
      inputdata = SimDiD(sample_size = sample_size, seed =1+seed, ATTcohortdiff = 2, minyear=2004, maxyear=2013)$simdata
      # DiDforBigData
      if(str_detect(estimator,"DiDforBigData")){
        time0 = proc.time()[3]
        parallel_cores = 1
        if (estimator == "DiDforBigDataMils") {
          parallel_cores = parallel::detectCores() - 2
        }
        My_p <- profmem({ My_res <- DiD(copy(inputdata), varnames = varnames, min_event = -1, max_event = 3) })
        time1 = proc.time()[3]
        My_time = (time1 - time0)/60
        My_mem = sum(My_p$bytes,na.rm=T)/1e9
        My_ATT = My_res$results_average[EventTime==1]$ATTe
        My_ATTse = My_res$results_average[EventTime==1]$ATTe_SE
        this_res = rbindlist(list(this_res, data.table(method="DiDforBigData", ATTe=My_ATT,ATTse=My_ATTse, Mem=My_mem, Time=My_time)))
        gc()
        print(sprintf("finished DiDforBigData in %s",My_time))
      }
      # Callaway Sant'Anna, did reg
      if("CSreg" %in% estimator){
        time0 = proc.time()[3]
        CSreg_p <- profmem({ CSreg_res <- CallawaySantanna(copy(inputdata), varnames, est_method = "reg") })
        time1 = proc.time()[3]
        CSreg_time = (time1 - time0)/60
        CSreg_mem = sum(CSreg_p$bytes,na.rm=T)/1e9
        CSreg_ATT = CSreg_res$results_average[EventTime==1]$ATTe
        CSreg_ATTse = CSreg_res$results_average[EventTime==1]$ATTe_SE
        this_res = rbindlist(list(this_res, data.table(method="CSreg", ATTe=CSreg_ATT, ATTse=CSreg_ATTse, Mem=CSreg_mem, Time=CSreg_time)))
        gc()
        print(sprintf("finished CSreg in %s",CSreg_time))
      }
      # Callaway Sant'Anna, did dr
      if("CSdr" %in% estimator){
        time0 = proc.time()[3]
        CSdr_p <- profmem({ CSdr_res <- CallawaySantanna(copy(inputdata), varnames, est_method = "dr") })
        time1 = proc.time()[3]
        CSdr_time = (time1 - time0)/60
        CSdr_mem = sum(CSdr_p$bytes,na.rm=T)/1e9
        CSdr_ATT = CSdr_res$results_average[EventTime==1]$ATTe
        CSdr_ATTse = CSdr_res$results_average[EventTime==1]$ATTe_SE
        this_res = rbindlist(list(this_res, data.table(method="CSdr", ATTe=CSdr_ATT, ATTse=CSdr_ATTse, Mem=CSdr_mem, Time=CSdr_time)))
        gc()
        print(sprintf("finished CSdr in %s",CSdr_time))
      }
      # Callaway Sant'Anna, did without bstrap
      if("CSbs" %in% estimator){
        time0 = proc.time()[3]
        CSbs_p <- profmem({ CSbs_res <- CallawaySantanna(copy(inputdata), varnames, bstrap=FALSE) })
        time1 = proc.time()[3]
        CSbs_time = (time1 - time0)/60
        CSbs_mem = sum(CSbs_p$bytes,na.rm=T)/1e9
        CSbs_ATT = CSbs_res$results_average[EventTime==1]$ATTe
        CSbs_ATTse = CSbs_res$results_average[EventTime==1]$ATTe_SE
        this_res = rbindlist(list(this_res, data.table(method="CSbs", ATTe=CSbs_ATT, ATTse=CSbs_ATTse, Mem=CSbs_mem, Time=CSbs_time)))
        gc()
        print(sprintf("finished CSbs in %s",CSbs_time))
      }
      # Borusyak et al, didimputation
      if("BJS" %in% estimator){
        time0 = proc.time()[3]
        BJS_p <- profmem({ BJS_res <-  Borusyaketal(copy(inputdata), varnames, horizon=c(0,1,2,3)) })
        time1 = proc.time()[3]
        BJS_time = (time1 - time0)/60
        BJS_mem = sum(BJS_p$bytes,na.rm=T)/1e9
        BJS_ATT = BJS_res[EventTime==1]$ATTe
        BJS_ATTse = BJS_res[EventTime==1]$ATTe_SE
        this_res = rbindlist(list(this_res, data.table(method="BJS", ATTe=BJS_ATT, ATTse=BJS_ATTse, Mem=BJS_mem, Time=BJS_time)))
        gc()
        print(sprintf("finished BJS in %s",BJS_time))
      }
      # de Chaisemartin & D'Haultfoeuille, DIDmultiplegt
      if(str_detect(estimator,"CH")){
        time0 = proc.time()[3]
        CH_p <- profmem({ CH_res <- deChaisemartin(copy(inputdata),varnames,dynamic=3,brep=40) })
        time1 = proc.time()[3]
        CH_time = (time1 - time0)/60
        CH_mem = sum(CH_p$bytes,na.rm=T)/1e9
        CH_ATT = CH_res[EventTime==1]$ATTe
        CH_ATTse = CH_res[EventTime==1]$ATTe_SE
        this_res = rbindlist(list(this_res, data.table(method="CH", ATTe=CH_ATT, ATTse=CH_ATTse, Mem=CH_mem, Time=CH_time)))
        gc()
        print(sprintf("finished CH in %s",CH_time))
      }
      # collect
      all_res = rbindlist(list(all_res,this_res))
    }
    if(reps>1){
      all_res = all_res[,list(ATTe=median(ATTe),ATTse=median(ATTse),Mem=median(Mem),Time=median(Time)),method]
    }
    all_res$sample_size = sample_size
    return(all_res)
  }

  all_speeds = data.table()
  for(nn in sample_sizes){
    print(sprintf("starting Sample Size (Unique Individuals) %s",nn))
    speeds = speed_tester(sample_size=nn,reps=reps)
    all_speeds = rbindlist(list(all_speeds, speeds))
    gc()
    print(all_speeds[])
    print(sprintf("sample %s done",nn))
  }

  file = sprintf("%s/speed_test_%s.csv", output_dir, estimator)
  write.csv(all_speeds, file=file, row.names = FALSE)

  return(all_speeds)

}


plot_results <- function(input_dir="inst/speed_tests", output_dir="docs/articles"){
  suppressMessages(library(ggplot2, warn.conflicts = F, quietly = T))
  suppressMessages(library(scales, warn.conflicts = F, quietly = T))
  suppressMessages(library(data.table, warn.conflicts = F, quietly = T))

  testresults = rbindlist(list(
    setDT(read.csv(file=sprintf("%s/speed_test_DiDforBigData.csv", input_dir))),
    setDT(read.csv(file=sprintf("%s/speed_test_CSdr.csv", input_dir))),
    setDT(read.csv(file=sprintf("%s/speed_test_CSbs.csv", input_dir))),
    setDT(read.csv(file=sprintf("%s/speed_test_CH1.csv", input_dir))),
    setDT(read.csv(file=sprintf("%s/speed_test_CH5.csv", input_dir))),
    setDT(read.csv(file=sprintf("%s/speed_test_CH10.csv", input_dir))),
    setDT(read.csv(file=sprintf("%s/speed_test_CH20.csv", input_dir))),
    setDT(read.csv(file=sprintf("%s/speed_test_BJS.csv", input_dir)))
  ))
  testresults = testresults[sample_size <= 2e4]
  testresults[method=="BJS", variable := "Borusyak Jaravel\n& Spiess\ndidimputation"]
  testresults[method=="CSdr", variable := "Callaway &\nSant'Anna\ndid"]
  testresults[method=="CSbs", variable := "Callaway &\nSant'Anna\ndid bstrap=F"]
  testresults[method=="CH", variable := "Chaisemartin &\nD'Haultfoeuille\nDIDmultiplegt"]
  testresults[method=="DiDforBigData", variable := "DiD for\nBig Data"]
  testresults[, sample_size_char := order(sample_size), variable]
  testresults[, variable := factor(variable,levels=c("Borusyak Jaravel\n& Spiess\ndidimputation", "Chaisemartin &\nD'Haultfoeuille\nDIDmultiplegt", "Callaway &\nSant'Anna\ndid", "Callaway &\nSant'Anna\ndid bstrap=F","DiD for\nBig Data"))]
  testresults = testresults[order(sample_size,variable)]

  gg = ggplot(aes(x=sample_size_char, y=Time, fill=variable),data=testresults)+
    theme_bw(base_size=22) + theme(legend.position = "bottom") +
    labs(x="Sample Size (Unique Individuals)",y="",title="Estimation Time (Minutes)",fill="") +
    scale_y_continuous(breaks=pretty_breaks()) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values=c('red','blue','black','purple','green')) +
    scale_x_continuous(breaks=c(1,2,3,4),labels=c("1"="1,000" , "2"="5,000" , "3"="10,000" , "4"="20,000"))
  ggsave(gg,filename=sprintf("%s/speedtest_small.png",output_dir),width=11,height=6)

  gg = ggplot(aes(x=sample_size_char, y=Mem, fill=variable),data=testresults)+
    theme_bw(base_size=22) + theme(legend.position = "bottom") +
    labs(x="Sample Size (Unique Individuals)",y="",title="Memory Usage (Gigabits)",fill="") +
    scale_y_continuous(breaks=pretty_breaks()) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values=c('red','blue','black','purple','green')) +
    scale_x_continuous(breaks=c(1,2,3,4),labels=c("1"="1,000" , "2"="5,000" , "3"="10,000" , "4"="20,000"))
  ggsave(gg,filename=sprintf("%s/memorytest_small.png",output_dir),width=11,height=6)

  # make ggplot
  variable_vals = sort(testresults[,unique(variable)])
  testresults[, variable := factor(variable, levels=variable_vals)]
  testresults = testresults[order(sample_size,variable)]
  ci_factor = 1.96
  testresults[, sample_size_char := as.character(sample_size)]
  testresults[, sample_size_num := order(as.numeric(sample_size)), method]
  mapping = unique(testresults[,.(sample_size_char,sample_size_num)])
  mapping_lab = c()
  for(ii in 1:nrow(mapping)){
    mapping_lab = c(mapping_lab, paste0(as.character(mapping[ii,sample_size_num]), "=", as.character(mapping[ii,sample_size_char])))
  }
  variables = testresults[, sort(unique(variable)) ]
  max_jitter = 0.25
  jitter_vals = seq(-max_jitter, max_jitter, length.out=length(variables))
  for(ii in 1:length(jitter_vals)){
    testresults[variable==variables[ii], jitter_val := jitter_vals[ii]]
  }
  testresults[, sample_size_num := sample_size_num + jitter_val ]
  gg <- ggplot(aes( x=sample_size_num, y = ATTe, color=variable), data = testresults) +
    labs(x = "Sample Size (Unique Individuals)", y ="", color="", title="DiD Estimates (95% CI)") +
    geom_point() + theme_bw(base_size = 20) +
    geom_errorbar(aes(ymin = ATTe - ci_factor * ATTse, ymax = ATTe + ci_factor * ATTse), width=0.1) +
    scale_y_continuous(breaks = pretty_breaks(), limits=c(3.5,4.5)) +
    geom_hline(yintercept=4, color="black", linetype="dashed") +
    theme(legend.position = "bottom") +
    scale_color_manual(values=c('red','blue','black','purple','green')) +
    scale_x_continuous(breaks=c(1,2,3,4),labels=c("1"="1,000" , "2"="5,000" , "3"="10,000" , "4"="20,000"))
  ggsave(gg,filename=sprintf("%s/estimates_small.png",output_dir),width=10,height=6)



  testresults = rbindlist(list(
    setDT(read.csv(file=sprintf("%s/speed_test_DiDforBigData.csv", input_dir))),
    setDT(read.csv(file=sprintf("%s/speed_test_CSdr.csv", input_dir))),
    setDT(read.csv(file=sprintf("%s/speed_test_CSbs.csv", input_dir)))
  ))
  testresults = testresults[sample_size > 2e4]
  testresults[method=="CSreg", variable := "Callaway &\nSant'Anna\ndid using reg"]
  testresults[method=="CSdr", variable := "Callaway &\nSant'Anna\ndid"]
  testresults[method=="CSbs", variable := "Callaway &\nSant'Anna\ndid bstrap=F"]
  testresults[method=="DiDforBigData", variable := "DiD for\nBig Data"]
  testresults[, sample_size_char := order(sample_size), variable]

  gg = ggplot(aes(x=sample_size_char, y=Time, fill=variable),data=testresults)+
    theme_bw(base_size=22) + theme(legend.position = "bottom") +
    labs(x="Sample Size (Unique Individuals)",y="",title="Estimation Time (Minutes)",fill="") +
    scale_y_continuous(breaks=pretty_breaks()) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values=c('black','purple','green')) +
    scale_x_continuous(breaks=c(1,2,3,4),labels=c("1"="50,000" , "2"="100,000" , "3"="500,000" , "4"="1,000,000"))
  ggsave(gg,filename=sprintf("%s/speedtest_large.png",output_dir),width=11,height=6)

  gg = ggplot(aes(x=sample_size_char, y=Mem, fill=variable),data=testresults)+
    theme_bw(base_size=22) + theme(legend.position = "bottom") +
    labs(x="Sample Size (Unique Individuals)",y="",title="Memory Usage (Gigabits)",fill="") +
    scale_y_continuous(breaks=pretty_breaks()) +
    geom_bar(stat='identity', position='dodge') +
    scale_fill_manual(values=c('black','purple','green')) +
    scale_x_continuous(breaks=c(1,2,3,4),labels=c("1"="50,000" , "2"="100,000" , "3"="500,000" , "4"="1,000,000"))
  ggsave(gg,filename=sprintf("%s/memorytest_large.png",output_dir),width=10,height=6)

}

## run one at a time:
# speedtest = TestDiD(estimator = "DiDforBigDataMils", sample_sizes=c(10e6), reps=3, input_dir=NULL, output_dir=NULL)
# speedtest = TestDiD(estimator = "DiDforBigData", sample_sizes=c(1e3,5e3,1e4,2e4,5e4,1e5,5e5,1e6), reps=3, input_dir=NULL, output_dir=NULL)
# speedtest = TestDiD(estimator = "CSdr", sample_sizes=c(1e3,5e3,1e4,2e4,5e4,1e5,5e5,1e6), reps=3, input_dir=NULL, output_dir=NULL)
# speedtest = TestDiD(estimator = "CSbs", sample_sizes=c(1e3,5e3,1e4,2e4,5e4,1e5,5e5,1e6), reps=3, input_dir=NULL, output_dir=NULL)
# speedtest = TestDiD(estimator = "BJS", sample_sizes=c(1e3,5e3,1e4,2e4), reps=3, input_dir=NULL, output_dir=NULL)
# speedtest = TestDiD(estimator = "CH1", sample_sizes=c(1e3), reps=3, input_dir=NULL, output_dir=NULL)
# speedtest = TestDiD(estimator = "CH5", sample_sizes=c(5e3), reps=3, input_dir=NULL, output_dir=NULL)
# speedtest = TestDiD(estimator = "CH10", sample_sizes=c(1e4), reps=3, input_dir=NULL, output_dir=NULL)
# speedtest = TestDiD(estimator = "CH20", sample_sizes=c(2e4), reps=3, input_dir=NULL, output_dir=NULL)
# plot_results(input_dir=NULL, output_dir=NULL)
