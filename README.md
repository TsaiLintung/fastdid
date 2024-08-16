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

You can control for covariates by providing the name of the data columns (`x, x2`), and choose the method via the `control_type` argument among doubly-robust (`"dr"`), inverse probability weight (`"ipw"`), and outcome regression (`"or"`). Ex: `fastdid(..., control_type = "dr", covaraitesvar = c("x", "x2"))`. 

Clustered standard error can be obtained from multiplier bootstrap. Ex: `fastdid(..., boot = TRUE, clustervar = "x")`. 

Estimation for multiple outcomes can be done in one call by providing a vector of outcome column names. This saves a lot of time when using `ipw` since logit estimates can be recycled across outcomes. Ex: `fastdid(...,outcomevar = c("y", "y2"))`.

# Performance

**fastdid** is magnitudes faster than **did**, and 15x faster than the fastest alternative **DiDforBigData** for large dataset. 

Here is a comparison of run time for **fastdid**, **did**, and **DiDforBigData** (dfbd for short) using a panel of 10 periods and varying samples sizes.

![time comparison](https://i.imgur.com/s5v32Rw.png)

Unfortunately, the Author's computer fails to run **did** at 1 million sample. For a rough idea, **DiDforBigData** is about 100x faster than **did** in Bradley Setzler's [benchmark](https://setzler.github.io/DiDforBigData/articles/Background.html). Other staggered DiD implementations are even slower than **did**. 

**fastdid** also uses less memory.

![RAM comparison](https://i.imgur.com/7emkgOz.png)

For the benchmark, a baseline group-time ATT is estimated with no covariates control, no bootstrap, and no explicit parallelization. Computing time is measured by `microbenchmark` and peak RAM by `peakRAM`.

# Validity

Before each release, we conduct tests to ensure the validity of estimates from `fastdid`.

## Basics: comparison with `did`

For features included in CS, `fastdid` maintains a maximum of 1% difference from results from the `did` package. This margin of error is mostly for bootstrapped results due to its inherent randomess. For point estimates, the difference is smaller than 1e-12, and is most likely the result of [floating-point error](https://en.wikipedia.org/wiki/Floating-point_error_mitigation). The relevant test files are [group-time](https://github.com/TsaiLintung/fastdid/blob/main/inst/tinytest/test_2_compare_gt.R) and [aggregates](https://github.com/TsaiLintung/fastdid/blob/main/inst/tinytest/test_3_compare_agg.R). 

## Extensions: direct coverage test

For features not included in CS, `fastdid` maintains that the 95% confidence intervals have a coverage rate between 94% and 96%. The coverage rate is calculated by running 200 iterations. In each iteration, we test whether the confidence interval estimated covers the group-truth values. We then average the rate across iterations. Due to the randomness of coverage, the realized coverage fall outside of the thresholds in about 1% of the time. The relevant test file is [coverage](https://github.com/TsaiLintung/fastdid/blob/main/inst/tinytest/test_5_coverage.R). 

## Experimental: not tested

Experimental features are not tested. The validity of its estimates are not guaranteed. 

# **fastdid** and **did**

As the name suggests, **fastdid**'s goal is to be fast **did**. Besides performance, here are some comparisons between the two packages.

## Estimator

**fastdid**'s estimators is identical to **did**'s. As the performance gains mostly come from efficient data manipulation, the key estimation implementations are analogous. For example, 2x2 DiD (`estimate_did.R` and `DRDID::std_ipw_did_panel`), influence function from weights (`aggregate_gt.R/get_weight_influence`, `compute.aggte.R/wif`), and multiplier bootstrap (`get_se.R` and `mboot.R`).

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

# Roadmap

**fastdid** is still in active development. Many features are planned to be added:

- Multiple outcomes :white_check_mark:
- Min/max event time and balanced composition :white_check_mark:
- DR and OR estimators :white_check_mark:
- Allowing for unbalanced panels :white_check_mark: (*well, not fully because DR and OR still need to be added*)
- Anticipation :white_check_mark:
- Varying base periods :white_check_mark:
- Simultaneously valid confidence bands :white_check_mark:
- User-provided aggregation scheme :white_check_mark:
- Zero-trust between steps for full customization
- Further optimization

# Source version

Since **fastdid** is not on CRAN yet, it needs to be converted to R scripts to be used in some restricted environments. This can be done with `development/build_source.R`. After changing the working directory, the script will produce `development/fastdid_VERNAME.R`, which can be sourced to mimic the functionalities of the package.

# Experimental features

`fastdid` is intended to be a fast implementation of `did` and nothing more. However, as the list of features grows, it has also become a more flexible version of `did`. While such flexibility can be great, these features often lack proof of validity and can easily be misused.

As an attempt to balance the validity and flexibility of `fastdid`, "experimental features" is introduced in version 0.9.4. These features will be less tested and documented, and it is generally advised to not use them unless the user know what they and the package are doing. These experimental features can be accessed via the `exper` argument. For example, to use the `filtervar` feature, call `fastdid(..., exper = list(filtervar = "FF"))`. 

The current list of experimental features are

- `max_control_cohort_diff`: limit the max cohort difference between treated and control group
- `filtervar`: limit the units being used as treated and control group with a potentially-time-varying variable
- `only_balance_2by2`: only require observations to have non-NA values within each 2 by 2 DiD, instead of throughout all time periods. Can be an alternative way of dealing with unbalanced panel by filling the missing periods with NAs. Not recommended as CS only have `allow_unbalance_panel`, which uses a repeated cross-section 2 by 2 DiD estimator.

# Update

## 0.9.9

remove: mincontrol cohort diff, min dynamic max dynamic
add: full, custom scheme


## 0.9.4 (2024/8/2)

> [!WARNING]
> Some BREAKING change is introduced in this update. 

- add uniform confidence interval option with `cband` and significance level `alpha`, confidence interval are now provided in result as column `att_ciub` and `att_cilb`
- BREAKING: `filtervar`, `max_control_cohort_diff`, `min_control_cohort_diff` are moved into the experimental features. See the above section for the explanation.
- add `max_dynamic` and `min_dynamic` as experimental features. 
- more informative error message when estimation fails for a specific `gt`, some internal interface overhaul

## 0.9.3 (2024/5/7)

- add anticipation and varying base period option
- add min and max control cohort difference
- add time-varying control ([reference](https://arxiv.org/abs/2202.02903))
- add filtervar 

0.9.3.1 (2024/5/24): fix the bug with `univar == clustervar` (TODO: address problems with name-changing and collision). 
0.9.3.2 (2024/7/17): fix group_time result when using `control_type = "notyet"` and make the base period in plots adapt to anticipation.
0.9.3.3 (2024/7/22): fix anticipation out of bound problem, more permanent solution for group_time target problem

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

**fastdid** is created and maintained by Lin-Tung Tsai. Many thanks to Maxwell Kellogg and Kuan-Ju Tseng for their contribution. 


