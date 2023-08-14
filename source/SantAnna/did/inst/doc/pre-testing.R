## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE, results="hide", warning=FALSE, message=FALSE----------------
library(did)
# Source the currently version of the did package (based on our Dropbox)
# fldr <- here::here("R/")
# sapply(paste0(fldr,list.files(fldr)), source)
# Source simulation designs
source(here::here("vignettes/setup_sims.R"))

## ----echo=FALSE---------------------------------------------------------------
time.periods <- 4
reset.sim()
bett <- betu <- rep(0, time.periods)
te <- 0
set.seed(1814)

## -----------------------------------------------------------------------------
# generate dataset with 4 time periods
time.periods <- 4

# generate dynamic effects
te.e <- time.periods:1

# generate data set with these parameters
# (main thing: it generates a dataset that satisfies
# parallel trends in all periods...including pre-treatment)
data <- build_sim_dataset()

head(data)

## ---- fig.width=8,fig.height=5, fig.align='center', out.width="90%", dpi = 200----
#-----------------------------------------------------------------------------
# modify the dataset a bit so that we can run an event study
#-----------------------------------------------------------------------------

# generate leads and lags of the treatment
Dtl <- sapply(-(time.periods-1):(time.periods-2), function(l) {
    dtl <- 1*( (data$period == data$G + l) & (data$G > 0) )
    dtl
})
Dtl <- as.data.frame(Dtl)
cnames1 <- paste0("Dtmin",(time.periods-1):1)
colnames(Dtl) <- c(cnames1, paste0("Dt",0:(time.periods-2)))
data <- cbind.data.frame(data, Dtl)
row.names(data) <- NULL

head(data)

#-----------------------------------------------------------------------------
# run the event study regression
#-----------------------------------------------------------------------------

# load plm package
library(plm)

# run event study regression
# normalize effect to be 0 in pre-treatment period
es <- plm(Y ~ Dtmin3 + Dtmin2 + Dt0 + Dt1 + Dt2, 
          data=data, model="within", effect="twoways",
          index=c("id","period"))

summary(es)

#-----------------------------------------------------------------------------
# make an event study plot
#-----------------------------------------------------------------------------

# some housekeeping for making the plot
# add 0 at event time -1
coefs1 <- coef(es)
ses1 <- sqrt(diag(summary(es)$vcov))
idx.pre <- 1:(time.periods-2)
idx.post <- (time.periods-1):length(coefs1)
coefs <- c(coefs1[idx.pre], 0, coefs1[idx.post])
ses <- c(ses1[idx.pre], 0, ses1[idx.post])
exposure <- -(time.periods-1):(time.periods-2)

cmat <- data.frame(coefs=coefs, ses=ses, exposure=exposure)

library(ggplot2)

ggplot(data=cmat, mapping=aes(y=coefs, x=exposure)) +
  geom_line(linetype="dashed") +
  geom_point() + 
  geom_errorbar(aes(ymin=(coefs-1.96*ses), ymax=(coefs+1.96*ses)), width=0.2) +
  ylim(c(-2,5)) +
  theme_bw()

## ---- fig.width=8,fig.height=10, fig.align='center', out.width="90%", dpi = 200----
# estimate group-group time average treatment effects
did_att_gt <- att_gt(yname="Y",
                     tname="period",
                     idname="id",
                     gname="G",
                     data=data,
                     bstrap=FALSE,
                     cband=FALSE)
summary(did_att_gt)

# plot them
ggdid(did_att_gt)


## ---- fig.width=8,fig.height=5, fig.align='center', out.width="90%", dpi = 200----
# aggregate them into event study plot
did_es <- aggte(did_att_gt, type="dynamic")

# plot the event study
ggdid(did_es)

## ---- echo=FALSE--------------------------------------------------------------
reset.sim()
bett <- betu <- rep(0, time.periods)
te <- 0
set.seed(1814)

## -----------------------------------------------------------------------------
# generate dataset with 4 time periods
time.periods <- 4

# generate dynamic effects
te.e <- time.periods:1

# generate selective treatment timing
# (*** this is what is different here ***)
te.bet.ind <- time.periods:1 / (time.periods/2)

# generate data set with these parameters
# (main thing: it generates a dataset that satisfies
# parallel trends in all periods...including pre-treatment)
data <- build_sim_dataset()

## -----------------------------------------------------------------------------
# run through same code as in earlier example...omitted

## ---- echo=FALSE--------------------------------------------------------------
#-----------------------------------------------------------------------------
# modify the dataset a bit so that we can run an event study
#-----------------------------------------------------------------------------

# generate leads and lags of the treatment
Dtl <- sapply(-(time.periods-1):(time.periods-2), function(l) {
    dtl <- 1*( (data$period == data$G + l) & (data$G > 0) )
    dtl
})
Dtl <- as.data.frame(Dtl)
cnames1 <- paste0("Dtmin",(time.periods-1):1)
colnames(Dtl) <- c(cnames1, paste0("Dt",0:(time.periods-2)))
data <- cbind.data.frame(data, Dtl)
row.names(data) <- NULL

#-----------------------------------------------------------------------------
# run the event study regression
#-----------------------------------------------------------------------------

# load plm package
library(plm)

## -----------------------------------------------------------------------------
# run event study regression
# normalize effect to be 0 in pre-treatment period
es <- plm(Y ~ Dtmin3 + Dtmin2 + Dt0 + Dt1 + Dt2, 
          data=data, model="within", effect="twoways", 
          index=c("id","period"))

summary(es)

## ---- echo=FALSE--------------------------------------------------------------
#-----------------------------------------------------------------------------
# make an event study plot
#-----------------------------------------------------------------------------

# some housekeeping for making the plot
# add 0 at event time -1
coefs1 <- coef(es)
ses1 <- sqrt(diag(summary(es)$vcov))
idx.pre <- 1:(time.periods-2)
idx.post <- (time.periods-1):length(coefs1)
coefs <- c(coefs1[idx.pre], 0, coefs1[idx.post])
ses <- c(ses1[idx.pre], 0, ses1[idx.post])
exposure <- -(time.periods-1):(time.periods-2)

cmat <- data.frame(coefs=coefs, ses=ses, exposure=exposure)

## ---- fig.width=8,fig.height=5, fig.align='center', out.width="90%", dpi = 200----
# run through same code as before...omitted

# new event study plot
ggplot(data=cmat, mapping=aes(y=coefs, x=exposure)) +
  geom_line(linetype="dashed") +
  geom_point() + 
  geom_errorbar(aes(ymin=(coefs-1.96*ses), ymax=(coefs+1.96*ses)), width=0.2) +
  ylim(c(-2,5)) +
  theme_bw()

## ---- fig.width=8,fig.height=10, fig.align='center', out.width="90%", dpi = 200----
# estimate group-group time average treatment effects
did.att.gt <- att_gt(yname="Y",
                     tname="period",
                     idnam="id",
                     gname="G",
                     data=data
                     )
summary(did.att.gt)

# plot them
ggdid(did.att.gt)

## ---- fig.width=8,fig.height=5, fig.align='center', out.width="90%", dpi = 200----
# aggregate them into event study plot
did.es <- aggte(did.att.gt, type="dynamic")

# plot the event study
ggdid(did.es)

## ----eval=FALSE---------------------------------------------------------------
#  # not run (this code can be substantially slower)
#  reset.sim()
#  set.seed(1814)
#  nt <- 1000
#  nu <- 1000
#  cdp <- conditional_did_pretest("Y", "period", "id", "G", xformla=~X, data=data)
#  cdp

