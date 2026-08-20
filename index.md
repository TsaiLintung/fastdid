# fastdid

[![R-CMD-check](https://github.com/TsaiLintung/fastdid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TsaiLintung/fastdid/actions/workflows/R-CMD-check.yaml)
[![](https://www.r-pkg.org/badges/version/fastdid?color=blue)](https://cran.r-project.org/package=fastdid)
[![CRANd](https://cranlogs.r-pkg.org/badges/grand-total/fastdid?color=blue)](https://CRAN.R-project.org/package=fastdid)

fastdid implements the difference-in-differences estimators of [Callaway
and Sant’Anna (2021)](https://doi.org/10.1016/j.jeconom.2020.12.001).
fastdid is:

- **fast**. On millions of units it cuts the computation time from hours
  to
  [seconds](https://tsailintung.github.io/fastdid/articles/misc.html#performance).
- **flexible**. It supports time-varying covariates ([Caetano and
  Callaway, 2024](https://arxiv.org/abs/2406.15288)) and multiple
  events, M \>= 2 ([Tsai, 2026](https://arxiv.org/abs/2409.05184)).

# Getting started

Install fastdid from CRAN:

``` r

install.packages("fastdid")
```

Or install the development version from GitHub:

``` r

# install.packages("devtools")
devtools::install_github("TsaiLintung/fastdid")
```

A call needs five things: the dataset `data`, and the column names for
time (`timevar`), cohort (`cohortvar`), unit (`unitvar`), and the
outcome or outcomes (`outcomevar`).

``` r

library(fastdid)
did_sim <- sim_did(1e+03, 10)                   # simulate some data
did_estimate <- fastdid(data = did_sim$dt, timevar = "time",
                        cohortvar = "G", unitvar = "unit", outcomevar = "y")
```

The function returns a `data.table` of estimates. Column `att` is the
point estimate. Column `se` is its standard error. Columns `att_cilb`
and `att_ciub` give the confidence interval. The remaining columns index
the estimated parameter.

To draw an event-study plot, call `plot_did_dynamics(did_estimate)`.

# More

- [did](https://bcallaway11.github.io/did/articles/did-basics.html) —
  staggered difference-in-differences, by Callaway and Sant’Anna
- [fastdid](https://tsailintung.github.io/fastdid/reference/fastdid.html)
  — the full list of arguments and features
- [double](https://tsailintung.github.io/fastdid/articles/double.html) —
  an introduction to DiD with multiple events. For M \>= 2 confounding
  events, pass a vector to `cohortvar2`, for example
  `cohortvar2 = c("G2", "G3")` for M = 3.
- [misc](https://tsailintung.github.io/fastdid/articles/misc.html) — the
  comparison with [did](https://github.com/bcallaway11/did), the
  benchmark, the tests, and the experimental features

# Acknowledgments

Lin-Tung Tsai created and maintains **fastdid**. Many thanks to Maxwell
Kellogg and Kuan-Ju Tseng for their contribution.
