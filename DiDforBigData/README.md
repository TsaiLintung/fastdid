DiD for Big Data in R
================

This R package provides a big-data-friendly and memory-efficient
difference-in-differences (DiD) estimator for staggered (and
non-staggered) treatment contexts. It supports controlling for
time-varying covariates, heteroskedasticity-robust standard errors, and
(single and multi-way) clustered standard errors. It addresses 4 issues
that arise in the context of large administrative datasets:

1.  **Speed:** In less than 1 minute, `DiDforBigData` will provide
    estimation and inference for staggered DiD with millions of
    observations on a personal laptop. It is orders of magnitude faster
    than other available software if the sample size is large; see the
    demonstration
    [here](https://setzler.github.io/DiDforBigData/articles/Background.html).
2.  **Memory:** Administrative data is often stored on crowded servers
    with limited memory available for researchers to use.
    `DiDforBigData` helps by using much less memory than other software;
    see the demonstration
    [here](https://setzler.github.io/DiDforBigData/articles/Background.html).
3.  **Dependencies:** Administrative servers often do not have outside
    internet access, making it difficult to install dependencies. This
    package has only *two* dependencies, `data.table` for big data
    management and `sandwich` for robust standard error estimation,
    which are already installed with most R distributions. Optionally,
    it will use the `fixest` package to speed up the estimation if it is
    installed. If the `progress` package is installed, it will also
    provide a progress bar so you know how much longer the estimation
    will take.
4.  **Parallelization:** Administrative servers often have a large
    number of available processors, but each processor may be slow, so
    it is important to parallelize. `DiDforBigData` makes
    parallelization easy as long as the `parallel` package is installed.

## Installation

To install the package from CRAN:

``` r
install.packages("DiDforBigData")
```

To install the package from Github:

``` r
devtools::install_github("setzler/DiDforBigData")
```

To use the package after it is installed:

``` r
library(DiDforBigData)
```

It is recommended to also make sure these optional packages have been
installed:

``` r
library("progress")
library("fixest")
library("parallel")
```

## Basic Use

There are only 3 functions in this package:

1.  `SimDiD()`: This function simulates data.
2.  `DiDge()`: This function estimates DiD for a single cohort and a
    single event time.
3.  `DiD()`: This function estimates DiD for all available cohorts and
    event times.

Details for each function are available from the [Function
Documentation](https://setzler.github.io/DiDforBigData/reference/index.html).

Before estimation, set up a variable list with the names of your
variables:

``` r
varnames = list()
varnames$time_name = "year"
varnames$outcome_name = "Y"
varnames$cohort_name = "cohort"
varnames$id_name = "id"
```

To estimate DiD for a single cohort and event time, use the `DiDge`
command. For example:

``` r
DiDge(inputdata = yourdata, varnames = varnames,
             cohort_time = 2010, event_postperiod = 3)
```

A detailed manual explaining the various features available in `DiDge`
is available
[here](https://setzler.github.io/DiDforBigData/reference/index.html) or
by running this command in R:

``` r
?DiDge
```

To estimate DiD for many cohorts and event times, use the `DiD` command.
For example:

``` r
DiD(inputdata = yourdata, varnames = varnames,
    min_event = -3, max_event = 5)
```

A detailed manual explaining the various features available in `DiD` is
available
[here](https://setzler.github.io/DiDforBigData/reference/index.html) or
by running this command in R:

``` r
?DiD
```

## Further Information

For more information, read the following articles:

- [Get
  Started](https://setzler.github.io/DiDforBigData/articles/DiDforBigData.html)
- [Background and
  Demonstration](https://setzler.github.io/DiDforBigData/articles/Background.html)
- [Theory and
  Methods](https://setzler.github.io/DiDforBigData/articles/Theory.html)
- [Function
  Documentation](https://setzler.github.io/DiDforBigData/reference/index.html)
- [Detailed
  Examples](https://setzler.github.io/DiDforBigData/articles/Examples.html)

Acknowledgements: Thanks to Mert Demirer and Kirill Borusyak for helpful
comments.
