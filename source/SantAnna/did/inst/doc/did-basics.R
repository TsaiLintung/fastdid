## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo=FALSE, results="hide", warning=FALSE, message=FALSE----------------
library(did)
# Source the currently version of the did package (based on our Dropbox)
#fldr <- here::here("R/")
#sapply(paste0(fldr,list.files(fldr)), source)
# Source simulation designs
#source(here::here("vignettes/setup_sims.R"))

## ---- echo=FALSE--------------------------------------------------------------
time.periods <- 4
sp <- reset.sim()
sp$te <- 0

## ---- message = FALSE, warning = FALSE----------------------------------------
# set seed so everything is reproducible
set.seed(1814)

# generate dataset with 4 time periods
time.periods <- 4

# add dynamic effects
sp$te.e <- 1:time.periods

# generate data set with these parameters
# here, we dropped all units who are treated in time period 1 as they do not help us recover ATT(g,t)'s.
dta <- build_sim_dataset(sp)

# How many observations remained after dropping the ``always-treated'' units
nrow(dta)
#This is what the data looks like
head(dta)

## -----------------------------------------------------------------------------
# estimate group-time average treatment effects using att_gt method
example_attgt <- att_gt(yname = "Y",
                        tname = "period",
                        idname = "id",
                        gname = "G",
                        xformla = ~X,
                        data = dta
                        )

# summarize the results
summary(example_attgt)

## ---- fig.width=8,fig.height=10, fig.align='center', out.width="90%", dpi = 200----
# plot the results
ggdid(example_attgt)

## -----------------------------------------------------------------------------
agg.simple <- aggte(example_attgt, type = "simple")
summary(agg.simple)

## ---- fig.width=8,fig.height=5, fig.align='center', out.width="90%", dpi = 200----
agg.es <- aggte(example_attgt, type = "dynamic")
summary(agg.es)
ggdid(agg.es)

## ---- fig.width=6, fig.align='center', out.width="90%", dpi = 200-------------
agg.gs <- aggte(example_attgt, type = "group")
summary(agg.gs)
ggdid(agg.gs)

## ---- fig.width=8, fig.height=5, fig.align='center', out.width="90%", dpi = 200----
agg.ct <- aggte(example_attgt, type = "calendar")
summary(agg.ct)
ggdid(agg.ct)

## ----eval=FALSE---------------------------------------------------------------
#  example_attgt_altcontrol <- att_gt(yname = "Y",
#                                     tname = "period",
#                                     idname = "id",
#                                     gname = "G",
#                                     xformla = ~X,
#                                     data = dta,
#                                     control_group = "notyettreated"			
#                                     )
#  summary(example_attgt_altcontrol)

## ----eval=FALSE---------------------------------------------------------------
#  example_attgt_reg <- att_gt(yname = "Y",
#                              tname = "period",
#                              idname = "id",
#                              gname = "G",
#                              xformla = ~X,
#                              data = data,
#                              est_method = "reg"
#                              )
#  summary(example_attgt_reg)

## -----------------------------------------------------------------------------
data(mpdta)

## -----------------------------------------------------------------------------
head(mpdta)

## ---- fig.width=8,fig.height=10, fig.align='center', out.width="90%", dpi = 200----
# estimate group-time average treatment effects without covariates
mw.attgt <- att_gt(yname = "lemp",
                   gname = "first.treat",
                   idname = "countyreal",
                   tname = "year",
                   xformla = ~1,
                   data = mpdta,
                   )

# summarize the results
summary(mw.attgt)

# plot the results
# set ylim so that all plots have the same scale along y-axis
ggdid(mw.attgt, ylim = c(-.3,.3))

## ----  fig.width=8,fig.height=5, fig.align='center', out.width="90%", dpi = 200----
# aggregate the group-time average treatment effects
mw.dyn <- aggte(mw.attgt, type = "dynamic")
summary(mw.dyn)
ggdid(mw.dyn, ylim = c(-.3,.3))

## ----  fig.width=8,fig.height=5, fig.align='center', out.width="90%", dpi = 200----
mw.dyn.balance <- aggte(mw.attgt, type = "dynamic", balance_e=1)
summary(mw.dyn.balance)
ggdid(mw.dyn.balance, ylim = c(-.3,.3))

## ----eval=FALSE---------------------------------------------------------------
#  mw.attgt.X <- att_gt(yname = "lemp",
#                     gname = "first.treat",
#                     idname = "countyreal",
#                     tname = "year",
#                     xformla = ~lpop,
#                     data = mpdta,
#                     )

