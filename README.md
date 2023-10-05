# fastdid - fast Difference-in-Differences

The **fastdid** package is a lightning-fast implementation of the staggered Difference-in-difference estimators proposed by ![https://www.sciencedirect.com/science/article/pii/S0304407620303948](Callaway and Sant'Anna 2021).

# Installation

You can install **fastdid** from GitHub (CRAN release coming soon.)

```
# install.packages("devtools")
devtools::install_github("TsaiLintung/fastdid")
```

# Getting started

A basic call.

```
library(fastdid)

#generate a simulated data
simdt <- sim_did(1e+03, 10, cov = "cont", second_cov = TRUE)
dt <- simdt$dt

result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time")
```

Control for covariates by providing the name of the data columns. 

```
result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                  covaraitesvar = c("x", "x2"))
```

Obtain clustered standard error from multiplier bootstrap. 

```
result <- fastdid(dt, timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",  result_type = "group_time",
                  clustervar = "x", boot = TRUE)
```

# Performance

# Comparison to **did**

# TODO

1. multiple outcome
2. balanced composition
3. DR estimator
4. Larger-than-memory dataset
