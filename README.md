# fastdid - fast Difference-in-Differences

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/TsaiLintung/fastdid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TsaiLintung/fastdid/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

**fastdid** is a lightning-fast implementation of [Callaway and Sant'Anna's (2021)](https://www.sciencedirect.com/science/article/pii/S0304407620303948) staggered Difference-in-Differences (DiD) estimators. DiD setup for millions of units used to take hours to run. With **fastdid**, it takes seconds. 

To learn more about the staggered Difference-in-differences estimators implemented, visit Callaway and Sant'Anna's [website](https://bcallaway11.github.io/did/articles/did-basics.html).

# Installation

You can install **fastdid** from GitHub.

```
# install.packages("devtools")
devtools::install_github("TsaiLintung/fastdid")
```

# Usage

`fastdid` is the main function provided by **fastdid**. When using `fastdid`, you need to provide the dataset (`data`), specify the names of the relevant columns (`-var`), and the type of target (aggregated) parameters (`result_type` such as `"group_time"`, `"time"`, `"dynamic"`, or `"simple"`.) Here is a simple call. 

```
#loading the package
library(fastdid)

#generate simulated data
simdt <- sim_did(1e+03, 10, cov = "cont", second_cov = TRUE, second_outcome = TRUE)
dt <- simdt$dt

#calling fastdid
result <- fastdid(data = dt, #the dataset
                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y", #name of the columns
                  result_type = "group_time") #the result type
```

You can control for covariates by providing the name of the data columns, and choose the method among doubly-robust (`"dr"`), inverse probability weight (`"ipw"`), and outcome regression (`"or"`). 

```
result <- fastdid(data = dt, 
                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                  result_type = "group_time",
                  control_option = "dr", #choose the control method
                  covaraitesvar = c("x", "x2")) #add covariates
```

While the default is to coerce the data into a balanced panel, you can allow for unbalanced panel. Note that currently only "ipw" is available when `allow_unbalance_panel = TRUE`.

```
result <- fastdid(data = dt, 
                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                  result_type = "group_time",
                  allow_unbalance_panel = TRUE, #allow for unbalanced panel
                  covaraitesvar = c("x", "x2"))
```


Clustered standard error can be obtained from multiplier bootstrap. 

```
result <- fastdid(data = dt,
                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = "y",
                  result_type = "group_time",
                  clustervar = "x", boot = TRUE) #add clustering by using bootstrap
```

Estimation for multiple outcomes can be done in one call by providing a vector of outcome column names (saves a lot of time when controlling for covariates since logit estimates can be recycled across outcomes). 

```
#calling fastdid
result <- fastdid(data = dt, #the dataset
                  timevar = "time", cohortvar = "G", unitvar = "unit", outcomevar = c("y", "y2"), #name of the columns
                  result_type = "group_time") #the result type
```

# Performance

**fastdid** is magnitudes faster than **did**, and 15x faster than the fastest alternative **DiDforBigData** for large dataset. 

Here is a comparison of run time for **fastdid**, **did**, and **DiDforBigData** (dfbd for short) using a panel of 10 periods and varying samples sizes.

![time comparison](https://i.imgur.com/s5v32Rw.png)

Unfortunately, the Author's computer fails to run **did** at 1 million sample. For a rough idea, **DiDforBigData** is about 100x faster than **did** in Bradley Setzler's [benchmark](https://setzler.github.io/DiDforBigData/articles/Background.html). Other staggered DiD implementations are even slower than **did**. 

**fastdid** also uses less memory.

![RAM comparison](https://i.imgur.com/7emkgOz.png)

For the benchmark, a baseline group-time ATT is estimated with no covariates control, no bootstrap, and no explicit parallelization. Computing time is measured by `microbenchmark` and peak RAM by `peakRAM`.

# **fastdid** and **did**

As the name suggests, **fastdid**'s goal is to be fast **did**. Besides performance, here are some comparisons between the two packages.

## Estimates

**fastdid**'s estimators is identical to **did**'s. As the performance gains mostly come from efficient data manipulation, the key estimation implementations are analogous. For example, 2x2 DiD (`estimate_did.R` and `DRDID::std_ipw_did_panel`), influence function from weights (`aggregate_gt.R/get_weight_influence`, `compute.aggte.R/wif`), and multiplier bootstrap (`get_se.R` and `mboot.R`).

Therefore, the estimates are practically identical. For point estimates, the difference is negligible (smaller than 1e-12), and is most likely the result of [floating-point error](https://en.wikipedia.org/wiki/Floating-point_error_mitigation). For standard errors, the estimates can be slightly different sometimes, but the difference never exceeds 1\% of **did**'s standard error estimates. 

## Interface

**fastdid** should feel very similar to `att_gt`. But there are a few differences:

Control group option: 
| fastdid | did | control group used |
|-|-|-|
| both | notyettreated | never-treated + not-yet-but-eventually-treated |
| never| nevertreated  | never-treated |
| notyet | | not-yet-but-eventually-treated |

Aggregated parameters: `fastdid` aggregates in the same function.
| fastdid | did |
|-|-|
| group_time | no aggregation |
|dynamic|dynamic|
|time|calendar|
|group|group|
|simple|simple|

## Other

1. **fastdid** only offers inverse probability weights estimators for controlling for covariates when allowing for unbalanced panels.
2. **fastdid** use universal base periods as default. 
4. **fastdid** only reports the pointwise confidence intervals, instead of the simultaneously valid confidence intervals (check section 4.1 of Callaway and Sant'Anna's (2021) for more detail.)

# Roadmap

**fastdid** is still in active development. Many features are planned to be added:

- Multiple outcomes :white_check_mark:
- Min/max event time and balanced composition :white_check_mark:
- DR and OR estimators :white_check_mark:
- Allowing for unbalanced panels :white_check_mark: (*well, not fully because DR and OR still need to be added*)
- Anticipation :white_check_mark:
- Varying base periods :white_check_mark:
- User-provided aggregation scheme
- User-provided control formula
- simultaneously valid confidence bands
- Further optimization

# Source version

Since **fastdid** is not on CRAN yet, it needs to be converted to R scripts to be used in some restricted environments. This can be done with `development/build_source.R`. After changing the working directory, the script will produce `development/fastdid_VERNAME.R`, which can be sourced to mimic the functionalities of the package.

# Update

## 0.9.3 (2024/5/7)

- add anticipation and varying base period option
- add min and max control cohort difference
- add time-varying control ([reference](https://arxiv.org/abs/2202.02903))
- add filtervar 

0.9.3.1 (2024/5/24): fix the bug with `univar == clustervar` (TODO: address problems with name-changing and collision)

## 0.9.2 (2023/12/20)

- add support to doubly robust and outcome regression estimators
- add support to unbalanced panels (simple and ipw only)
- add support to balanced composition option in dynamics aggregation
- fixed argument checking that was not working properly
- set the default to copying the entire dataset to avoid unexpected modification of the original data (thanks @grantmcdermott for the suggestion.)

## 0.9.1 (2023/10/20)

- now supprts estimation for multiple outcomes in one go! 
- data validation: no longer check missing values for columns not used. 

# Acknowledgments

**fastdid** is created by Lin-Tung Tsai, Maxwell Kellogg, and Kuan-Ju Tseng. 


