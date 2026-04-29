# Changelog

## fastdid 1.0.7

- Extended double DiD to support M\>2 treatment events: `cohortvar2` now
  accepts a character vector of length M-1 (e.g. `c("G2", "G3")` for
  three events)
- Added `add_base_period` parameter: inserts a zero-ATT placeholder at
  the base period in `result_type = "dynamic"` results
- Added experimental options `only_est_min` and `only_est_max` in
  `exper`: restrict estimation to a specific event-time range in dynamic
  mode, skipping g-t pairs outside the window
- Improved input validation: errors on negative
  `anticipation`/`anticipation2`, time-varying weights, and
  `balanced_event_time` exceeding data range
- Restored `parglm` dependency for multi-threaded propensity score
  estimation; thread count automatically matched to
  [`getDTthreads()`](https://rdrr.io/pkg/data.table/man/openmp-utils.html)
  (or 1 when `parallel = TRUE`)
- Various bug fixes and robustness improvements for double DiD
  aggregation and weight handling

## fastdid 1.0.6

CRAN release: 2026-01-14

- Remove archived parglm dependency

## fastdid 1.0.5

CRAN release: 2025-06-23

- Various bug fixes

## fastdid 1.0.4

- Fixed Double DiD inference: account for stolastic part of the double
  DiD weight

## fastdid 1.0.3

CRAN release: 2024-11-04

- Fixed Check for CRAN

## fastdid 1.0.2

CRAN release: 2024-10-28

- Fixed Typo for CRAN

## fastdid 1.0.1

- Fixed Typo for CRAN

## Version 1.0.0

- Release to CRAN!

## Version 0.9.9

- add double did (see the vignette for the introduction)
- add `parallel`, parallization for unix systems, useful if the number
  of g-t is large.
- add `full`, return full result such as influence function, aggregate
  scheme, and such
- add `min`/`max_dynamic`, `custom_scheme` to experimental features

0.9.9.1 (2024/9/13): fix a bug that affects not-yet control with max
treated group != max time

## Version 0.9.4

> Some BREAKING change is introduced in this update.

- add uniform confidence interval option with `cband` and significance
  level `alpha`, confidence interval are now provided in result as
  column `att_ciub` and `att_cilb`
- BREAKING: `filtervar`, `max_control_cohort_diff`,
  `min_control_cohort_diff` are moved into the experimental features.
  See the above section for the explanation.
- add `max_dynamic` and `min_dynamic` as experimental features.
- more informative error message when estimation fails for a specific
  `gt`, some internal interface overhaul

## Version 0.9.3

- add anticipation and varying base period option
- add min and max control cohort difference
- add time-varying control
  ([reference](https://arxiv.org/abs/2202.02903))
- add filtervar

0.9.3.1 (2024/5/24): fix the bug with `univar == clustervar` (TODO:
address problems with name-changing and collision). 0.9.3.2 (2024/7/17):
fix group_time result when using `control_type = "notyet"` and make the
base period in plots adapt to anticipation. 0.9.3.3 (2024/7/22): fix
anticipation out of bound problem, more permanent solution for
group_time target problem

## Version 0.9.2

- add support to doubly robust and outcome regression estimators
- add support to unbalanced panels (simple and ipw only)
- add support to balanced composition option in dynamics aggregation
- fixed argument checking that was not working properly
- set the default to copying the entire dataset to avoid unexpected
  modification of the original data (thanks
  [@grantmcdermott](https://github.com/grantmcdermott) for the
  suggestion.)

## Version 0.9.1

- now supprts estimation for multiple outcomes in one go!
- data validation: no longer check missing values for columns not used.
