# Fast Staggered DID Estimation

Performs Difference-in-Differences (DID) estimation.

## Usage

``` r
fastdid(
  data,
  timevar,
  cohortvar,
  unitvar,
  outcomevar,
  control_option = "both",
  result_type = "group_time",
  balanced_event_time = NA,
  control_type = "ipw",
  allow_unbalance_panel = FALSE,
  boot = FALSE,
  biters = 1000,
  cband = FALSE,
  alpha = 0.05,
  weightvar = NA,
  clustervar = NA,
  covariatesvar = NA,
  varycovariatesvar = NA,
  copy = TRUE,
  validate = TRUE,
  anticipation = 0,
  anticipation2 = 0,
  base_period = "universal",
  exper = NULL,
  full = FALSE,
  parallel = FALSE,
  cohortvar2 = NA,
  event_specific = TRUE,
  double_control_option = "both"
)
```

## Arguments

- data:

  data.table, the dataset.

- timevar:

  character, name of the time variable.

- cohortvar:

  character, name of the cohort (group) variable.

- unitvar:

  character, name of the unit (id) variable.

- outcomevar:

  character vector, name(s) of the outcome variable(s).

- control_option:

  character, control units used for the DiD estimates, options are
  "both", "never", or "notyet".

- result_type:

  character, type of result to return, options are "group_time", "time",
  "group", "simple", "dynamic" (time since event), "group_group_time",
  or "dynamic_stagger".

- balanced_event_time:

  number, max event time to balance the cohort composition.

- control_type:

  character, estimator for controlling for covariates, options are "ipw"
  (inverse probability weighting), "reg" (outcome regression), or "dr"
  (doubly-robust).

- allow_unbalance_panel:

  logical, allow unbalance panel as input or coerce dataset into one.

- boot:

  logical, whether to use bootstrap standard error.

- biters:

  number, bootstrap iterations. Default is 1000.

- cband:

  logical, whether to use uniform confidence band or point-wise.

- alpha:

  number, the significance level. Default is 0.05.

- weightvar:

  character, name of the weight variable.

- clustervar:

  character, name of the cluster variable.

- covariatesvar:

  character vector, names of time-invariant covariate variables.

- varycovariatesvar:

  character vector, names of time-varying covariate variables.

- copy:

  logical, whether to copy the dataset.

- validate:

  logical, whether to validate the dataset.

- anticipation:

  number, periods with anticipation.

- anticipation2:

  number, periods with anticipation for the second event.

- base_period:

  character, type of base period in pre-preiods, options are
  "universal", or "varying".

- exper:

  list, arguments for experimental features.

- full:

  logical, whether to return the full result (influence function, call,
  weighting scheme, etc,.).

- parallel:

  logical, whether to use parallization on unix system.

- cohortvar2:

  character, name of the second cohort (group) variable.

- event_specific:

  logical, whether to recover target treatment effect or use combined
  effect.

- double_control_option:

  character, control units used for the double DiD, options are "both",
  "never", or "notyet".

## Value

A data.table containing the estimated treatment effects and standard
errors or a list of all results when \`full == TRUE\`.

## Details

\`balanced_event_time\` is only meaningful when \`result_type ==
"dynamic\`.

\`result_type\` as \`group-group-time\` and \`dynamic staggered\` is
only meaningful when using double did.

\`biter\` and \`clustervar\` is only used when \`boot == TRUE\`.

## Examples

``` r
# simulated data
simdt <- sim_did(1e+02, 10, cov = "cont", second_cov = TRUE, second_outcome = TRUE, seed = 1)
dt <- simdt$dt

# basic call
result <- fastdid(
  data = dt, timevar = "time", cohortvar = "G",
  unitvar = "unit", outcomevar = "y",
  result_type = "group_time"
)
```
