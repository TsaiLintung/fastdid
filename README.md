# fastdid - fast Difference-in-Differences

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/TsaiLintung/fastdid/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/TsaiLintung/fastdid/actions/workflows/R-CMD-check.yaml)
  [![codecov](https://codecov.io/github/TsaiLintung/fastdid/graph/badge.svg?token=0EYF1DIBRF)](https://codecov.io/github/TsaiLintung/fastdid)
  <!-- badges: end -->

fastdid is a fast implementation of [Callaway and Sant'Anna's (2021)](https://www.sciencedirect.com/science/article/pii/S0304407620303948) staggered Difference-in-Differences (DiD) estimators, it 

  - reduces the computation time with millions of units from hours to seconds,
  - incorporates extensions such as [time-varying covariates](https://arxiv.org/abs/2202.02903) and [double event](not_ready_yet). 

fastdid can be installed from GitHub. 

```
# install.packages("devtools")
devtools::install_github("TsaiLintung/fastdid")
```

Here is a simple call to `fastdid`:

```
library(fastdid) #loading the package
did_sim <- sim_did(1e+03, 10) #generate simulated data 
did_estimate <- fastdid(data = did_sim$dt, timevar = "time", #calling fastdid
  cohortvar = "G", unitvar = "unit", outcomevar = "y")
```

To learn more about fastdid: 

  - [basics](https://bcallaway11.github.io/did/articles/did-basics.html): basics of staggered Difference in Difference by Callaway and Sant'Anna
  - [usage](articles/usage.html): guide to fastdid
  - [fastdid and did](articles/did.html): comparison with [did](https://github.com/bcallaway11/did) package
  - [performance](articles/performance.html): memory and speed benchmark
  - [validity](articles/validity.html): testing conducted to ensure the validity
  - [double](articles/double.html): guide to double event 
  - [experimental](articles/experimental.html): flexible but less tested features

# Acknowledgments

**fastdid** is created and maintained by Lin-Tung Tsai. Many thanks to Maxwell Kellogg and Kuan-Ju Tseng for their contribution. 


