# STATE

What the package is now: the layout, the data flow, and the dependencies.

## Public API

Key `fastdid()` parameters:

| Parameter | Meaning |
|---|---|
| `data` | a data.table, or an object that coerces to one |
| `timevar`, `cohortvar`, `unitvar`, `outcomevar` | column names, as strings |
| `control_option` | `"both"`, `"never"`, or `"notyet"` |
| `result_type` | `"group_time"`, `"time"`, `"group"`, `"simple"`, `"dynamic"`, `"group_group_time"`, `"dynamic_stagger"`, or `"dynamic_event"` |
| `control_type` | `"ipw"`, `"reg"`, `"dr"` (doubly robust), or `"strata"` (the DiD within each value of the covariates, from means) |
| `cohortvar2` | the confounding events, one column name per event, for the two-stage estimator |
| `effect_model` | the second-stage model: `"parallel"`, `"unrestricted"`, or a one-sided formula on the cell features |
| `effect_fit` | how a formula is fit: `"separate"`, `"joint"`, `"ordered"`, or `"state"` |
| `boot` | bootstrap standard errors. The default is 1000 iterations |
| `base_period` | `"universal"` or `"varying"` |
| `covariatesvar`, `varycovariatesvar` | time-invariant and time-varying covariates |
| `anticipation`, `anticipation2` | anticipation periods |
| `exper` | experimental features, such as `filtervar` and `aggregate_scheme` |
| `full` | returns the influence functions and the weights as well |
| `parallel` | Unix parallelization through `mclapply` |

## Internal architecture

### File map

```
R/
├── fastdid.R          entry point: validate, coerce, estimate, aggregate
├── validate.R         input validation through dreamerr
├── aux_funcs.R        data coercion, aux data structures, locked lists
├── estimate_gtatt.R   core loop over cohorts and times
├── estimate_did.R     the 2x2 engine: IPW, OR, doubly robust, influence functions
├── first_stage.R      the first-stage object, double_did(), and the lazy linear map to the targets
├── aggregate_gt.R     the aggregation schemes and the influence of the aggregation weights
├── second_stage.R     the cell table, and the second stage as weights over first-stage cells
├── effect_model.R     the formula-based second stage: design, WLS, estimability, diagnostics
├── strata.R           the strata first stage from means, its influence function, combine_cells()
├── double_did.R       multiple events: the G-string label helpers and the coercion
├── sim_did.R          data simulation
├── generics.R         S3 methods for the fastdid_result class
└── global.R           global-variable suppressions for R CMD check
```

### Data flow

```
fastdid()
  -> validate_argument() + validate_dt()
  -> coerce_dt()        normalize time and cohort to 1,2,3...; store time_offset, time_step
  -> get_auxdata()      fast-access structures, clustered by time period
  -> estimate_gtatt()   per outcome: cohorts x times, then estimate_did(), or
                        estimate_gtatt_outcome_strata() for control_type = "strata"
  -> new_first_stage()  per outcome: the cell table, the shares, and the uncertainty object
  -> aggregate_gt()     per outcome: run_double_did(): second_stage(), the targets, and the SEs
  -> convert_targets()  map back to the original time and cohort scales
```


## Dependencies

- **Imports**: `data.table (>= 1.15.0)`, `stringr`, `BMisc`, `collapse`,
  `dreamerr (>= 1.4.0)`, `parglm`, `ggplot2`
- **Suggests**: `did`, `knitr`, `rmarkdown`, `parallel`, `tinytest`
- **R**: >= 4.1.0
