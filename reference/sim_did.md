# Simulate a Difference-in-Differences (DiD) dataset

Simulates a dataset for a Difference-in-Differences analysis with
various customizable options.

## Usage

``` r
sim_did(
  sample_size,
  time_period,
  untreated_prop = 0.3,
  epsilon_size = 0.001,
  cov = "no",
  hetero = "all",
  second_outcome = FALSE,
  second_cov = FALSE,
  vary_cov = FALSE,
  na = "none",
  balanced = TRUE,
  seed = NA,
  stratify = FALSE,
  treatment_assign = "latent",
  second_cohort = FALSE,
  confound_ratio = 1,
  second_het = "all"
)
```

## Arguments

- sample_size:

  The number of units in the dataset.

- time_period:

  The number of time periods in the dataset.

- untreated_prop:

  The proportion of untreated units.

- epsilon_size:

  The standard deviation for the error term in potential outcomes.

- cov:

  The type of covariate to include ("no", "int", or "cont").

- hetero:

  The type of heterogeneity in treatment effects ("all" or "dynamic").

- second_outcome:

  Whether to include a second outcome variable.

- second_cov:

  Whether to include a second covariate.

- vary_cov:

  include time-varying covariates

- na:

  Whether to generate missing data ("none", "y", "x", or "both").

- balanced:

  Whether to balance the dataset by random sampling.

- seed:

  Seed for random number generation.

- stratify:

  Whether to stratify the dataset based on a binary covariate.

- treatment_assign:

  The method for treatment assignment ("latent" or "uniform").

- second_cohort:

  include confounding events

- confound_ratio:

  extent of event confoundedness

- second_het:

  heterogeneity of the second event

## Value

A list containing the simulated dataset (dt) and the treatment effect
values (att).

## Examples

``` r
# Simulate a DiD dataset with default settings
data <- sim_did(sample_size = 100, time_period = 5)
```
