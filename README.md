# fastdid - fast Difference-in-Differences

Welcome to fastdid: A lightning-fast implementation of [Callaway and Sant'Anna's (2021)](https://www.sciencedirect.com/science/article/pii/S0304407620303948) staggered Difference-in-Difference estimators.

**fastdid** is not created because we need another DiD estimator, but because existing solutions are often too slow to be used on adminstrative data with limited computing resource. **fastdid** is just trying to be fast **did**. With **fastdid**, you can estimate DiD designs on millions of units in just seconds, not hours.

# Installation

You can install **fastdid** from GitHub (CRAN release coming soon.)

```
# install.packages("devtools")
devtools::install_github("TsaiLintung/fastdid")
```

# Getting started

`fastdid` is the main function provided by **fastdid**. When using `fastdid`, you need to provide the dataset (`dt`), specify names of the relevant columns (`-var`), and the type of target (aggregated) parameters (`result_type` such as `"group_time"`, `"time"`, `"dynamic"`, or `"simple"`). Here is a simple call. 

```
#loading the library
library(fastdid)

#generate a simulated data
simdt <- sim_did(1e+03, 10, cov = "cont", second_cov = TRUE)
dt <- simdt$dt

#calling fastdid
result <- fastdid(dt, #the dataset
                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", #name of the columns
                  result_type = "group_time") #the result type
```

You can control for covariates by providing the name of the data columns. 

```
result <- fastdid(dt, 
                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                  result_type = "group_time",
                  covaraitesvar = c("x", "x2")) #add covariates
```

Clustered standard error can be obtained from multiplier bootstrap. 

```
result <- fastdid(dt,
                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                  result_type = "group_time",
                  clustervar = "x", boot = TRUE) #add clustering by using bootstrap
```

# Performance

**fastdid** is magnitudes faster than **did** and about 15x faster than the fastest alternative **DiDforBigData**. Here is a comparison of run time for **fastdid**, **did**, and **DiDforBigData** (dfbd for short) at different sample size and 10 periods.

![time comparison](https://i.imgur.com/s5v32Rw.png)

The Author's computer fails to execute **did** at 1 million sample. For a rough idea, **DiDforBigData** is about 100x faster than **did** in Bradley Setzler's [benchmark](https://setzler.github.io/DiDforBigData/articles/Background.html). Other staggered DiD implementations are even slower than **did**. 

The time gain is not from time-memory trade-offs, as **fastdid** also uses less memory.

![RAM comparison](https://i.imgur.com/7emkgOz.png)

# **fastdid** and **did**

**fastdid**'s estimators is identical to the ones in **did**. As the performance gains mostly come from efficient data manipulation, the implementations of key steps of estimation are analogous. For example, 2x2 DiD (`estimate_did.R` and `DRDID::std_ipw_did_panel`), influence function from weights (`aggregate_gt.R/get_weight_influence`, `compute.aggte.R/wif`), and multiplier bootstrap (`get_se.R` and `mboot.R`).

Therefore, The estimates are practically identical. The negligible difference (smaller than 1e-12) in point estimates can be atttributed to [floating-point error](https://en.wikipedia.org/wiki/Floating-point_error_mitigation). The standard errors are  slightly different only when clustering (due to randomness from bootstrap) or controlling for covariates (due to different package used for logit estimation), and the difference is always less than 1\% of standard error.

Notably, **fastdid** only offers inverse probability weights estimators for controlling for covariates (OR and DR likely to be added), and can only deal with panel data, but not repeated cross-section (unlikely to be added).

To learn more about staggered Difference-in-differences estimators, visit [[Callaway and Sant'Anna's](https://github.com/bcallaway11/did)'s website.

# TODO

**fastdid** is still in active development, many features are planned to be added.

1. Multiple outcomes
2. Balanced composition
3. DR and OR estimators
4. Larger-than-memory data support
5. User-provided aggregation scheme
6. Frop-in interface to did
7. Further optimization!
