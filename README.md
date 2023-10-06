# fastdid - fast Difference-in-Differences

**fastdid** is a lightning-fast implementation of [Callaway and Sant'Anna's (2021)](https://www.sciencedirect.com/science/article/pii/S0304407620303948) staggered Difference-in-Differences (DiD) estimators. With **fastdid**, you can run DiD setup for millions of units in just seconds, not hours. 

To learn more about the staggered Difference-in-differences estimators implemented, visit Callaway and Sant'Anna's [website](https://bcallaway11.github.io/did/articles/did-basics.html).

# Installation

You can install **fastdid** from GitHub (CRAN release coming soon.)

```
# install.packages("devtools")
devtools::install_github("TsaiLintung/fastdid")
```

# Usage

`fastdid` is the main function provided by **fastdid**. When using `fastdid`, you need to provide the dataset (`dt`), specify the names of the relevant columns (`-var`), and the type of target (aggregated) parameters (`result_type` such as `"group_time"`, `"time"`, `"dynamic"`, or `"simple"`.) Here is a simple call. 

```
#loading the package
library(fastdid)

#generate simulated data
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

**fastdid** is magnitudes faster than **did**, and 15x faster than the fastest alternative **DiDforBigData** for large dataset. 

Here is a comparison of run time for **fastdid**, **did**, and **DiDforBigData** (dfbd for short) using a panel of 10 periods and varying samples sizes.

![time comparison](https://i.imgur.com/s5v32Rw.png)

Unfortunately, the Author's computer fails to run **did** at 1 million sample. For a rough idea, **DiDforBigData** is about 100x faster than **did** in Bradley Setzler's [benchmark](https://setzler.github.io/DiDforBigData/articles/Background.html). Other staggered DiD implementations are even slower than **did**. 

**fastdid** also uses less memory.

![RAM comparison](https://i.imgur.com/7emkgOz.png)

For the benchmark, a baseline group-time ATT is estimated with no covariates control, no bootstrap, no explicit parallelization. Computing time is measured by `microbenchmark` and peak RAM by `peakRAM`.

# **fastdid** and **did**

As the name suggests, **fastdid**'s goal is to be fast **did**. Besides performance, here are some comparisons between the two packages .

## Estimates

**fastdid**'s estimators is identical to **did**'s. As the performance gains mostly come from efficient data manipulation, the key estimation implementation are analogous. For example, 2x2 DiD (`estimate_did.R` and `DRDID::std_ipw_did_panel`), influence function from weights (`aggregate_gt.R/get_weight_influence`, `compute.aggte.R/wif`), and multiplier bootstrap (`get_se.R` and `mboot.R`).

Therefore, the estimates are practically identical. For point estimates, the difference is negligible (smaller than 1e-12), and is most likely the result of [floating-point error](https://en.wikipedia.org/wiki/Floating-point_error_mitigation).

For standard errors, the estimates can be slightly different in certain situations, but the difference never exceeds 1\% of **did**'s standard error estimates. The situations where estimates differ include clustering, due to randomness from bootstrap, and controlling for covariates, due to different package used for logit estimation. 

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

## Feature

Notable differences in feature includes:
1. **fastdid** currently only offer inverse probability weights estimators for controlling for covariates (OR and DR likely to be added soon)
2. **fastdid** only uses the time before event as base periods ("universal" in `attgt`)
3. **fastdid** can only deal with balanced panel, no repeated cross-sections, no missing observations.

# Roadmap

**fastdid** is still in active development, many features are planned to be added:

1. Multiple outcomes
2. Min/max event time and balanced composition
3. DR and OR estimators
4. Larger-than-memory data support
5. User-provided aggregation scheme
6. drop-in interface for did
7. Anticipation
8. Further optimization!

# Acknowledgments

**fastdid** is created by Maxwell Kellogg, Lin-Tung Tsai, and Kuan-Ju Tseng. Contact the authors by either email
