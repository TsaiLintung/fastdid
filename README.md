# fastdid

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/TsaiLintung/fastdid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TsaiLintung/fastdid/actions/workflows/R-CMD-check.yaml)
  [![](https://www.r-pkg.org/badges/version/fastdid?color=blue)](https://cran.r-project.org/package=fastdid)
  [![CRANd](https://cranlogs.r-pkg.org/badges/grand-total/fastdid?color=blue)](https://cran.r-project.org/web/packages/fastdid/index.html)
  <!-- badges: end -->

fastdid implements the Difference-in-Differences (DiD) estimators in [Callaway and Sant'Anna's (2021)](https://doi.org/10.1016/j.jeconom.2020.12.001). fastdid is

  - fast, reducing the computation time with millions of units from hours to [seconds](https://tsailintung.github.io/fastdid/articles/misc.html#performance),
  - flexible, allowing extensions such as time-varying covariates [Caetano and Callaway, 2024](https://arxiv.org/abs/2406.15288) and multiple events ([Tsai 2024](https://arxiv.org/abs/2409.05184)). 

# Getting Started

fastdid can be installed from CRAN,

```
install.packages("fastdid")
```

or the latest developmental version can be installed via GitHub,

```
# install.packages("devtools")
devtools::install_github("TsaiLintung/fastdid")
```

To use `fastdid`, you need to provide the dataset `data`, the column name of time `timevar`, cohort `cohortvar`, unit `unitvar`, and outcome(s) `outcomevar`. Here is a simple call:

```
library(fastdid) #loading the package
did_sim <- sim_did(1e+03, 10) #simulate some data
did_estimate <- fastdid(data = did_sim$dt, timevar = "time",
                  cohortvar = "G", unitvar = "unit", outcomevar = "y")
```
The function returns a `data.table` that includes the estimates. Column `att` is the point estimate, `se` the standard error of the estimate, `att_ciub` and `att_cilb` the confidence interval. The other columns indexes the estimated parameter. 

To create event study plots, use `plot_did_dynamics(did_estimate)`. 

# More

  - [did](https://bcallaway11.github.io/did/articles/did-basics.html): staggered Difference in Difference by Callaway and Sant'Anna
  - [fastdid](https://tsailintung.github.io/fastdid/reference/fastdid.html): full list of arguments and features.
  - [double](https://tsailintung.github.io/fastdid/articles/double.html): introduction to DiD with multiple events.
  - [misc](https://tsailintung.github.io/fastdid/articles/misc.html): comparison with [did](https://github.com/bcallaway11/did), benchmark, tests, and experimental features.

# Acknowledgments

**fastdid** is created and maintained by Lin-Tung Tsai. Many thanks to Maxwell Kellogg and Kuan-Ju Tseng for their contribution. 


