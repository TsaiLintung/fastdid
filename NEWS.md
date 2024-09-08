# Version 0.9.9

- add double did (see the vignette for the introduction)
- add `parallel`, parallization for unix systems, useful if the number of g-t is large. 
- add `full`, return full result such as influence function, aggregate scheme, and such
- add `min`/`max_dynamic`, `custom_scheme` to experimental features

# Version  0.9.4

> Some BREAKING change is introduced in this update. 

- add uniform confidence interval option with `cband` and significance level `alpha`, confidence interval are now provided in result as column `att_ciub` and `att_cilb`
- BREAKING: `filtervar`, `max_control_cohort_diff`, `min_control_cohort_diff` are moved into the experimental features. See the above section for the explanation.
- add `max_dynamic` and `min_dynamic` as experimental features. 
- more informative error message when estimation fails for a specific `gt`, some internal interface overhaul

# Version 0.9.3

- add anticipation and varying base period option
- add min and max control cohort difference
- add time-varying control ([reference](https://arxiv.org/abs/2202.02903))
- add filtervar 

0.9.3.1 (2024/5/24): fix the bug with `univar == clustervar` (TODO: address problems with name-changing and collision). 
0.9.3.2 (2024/7/17): fix group_time result when using `control_type = "notyet"` and make the base period in plots adapt to anticipation.
0.9.3.3 (2024/7/22): fix anticipation out of bound problem, more permanent solution for group_time target problem

# Version  0.9.2

- add support to doubly robust and outcome regression estimators
- add support to unbalanced panels (simple and ipw only)
- add support to balanced composition option in dynamics aggregation
- fixed argument checking that was not working properly
- set the default to copying the entire dataset to avoid unexpected modification of the original data (thanks @grantmcdermott for the suggestion.)

# Version  0.9.1

- now supprts estimation for multiple outcomes in one go! 
- data validation: no longer check missing values for columns not used. 