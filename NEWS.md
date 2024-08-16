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