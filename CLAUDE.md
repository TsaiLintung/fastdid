# CLAUDE.md

Guidance for Claude Code (claude.ai/code) in this repository.

fastdid is an R package. It implements the staggered difference-in-differences
estimator of Callaway and Sant'Anna (2021). It adds two extensions: time-varying
covariates (Caetano and Callaway, 2024) and multiple treatment events (Tsai,
2026). It is on CRAN. The current version is in `DESCRIPTION`.

The package trades memory for speed. It sorts once, then reads the data through
fast-access structures. On millions of units it cuts estimation from hours to
seconds.

## Public API

fastdid exports seven functions:

| Function | What it does |
|---|---|
| `fastdid()` | the main estimator |
| `plot_did_dynamics()` | an event-study plot through ggplot2 |
| `sim_did()` | simulates a staggered DiD dataset for testing |
| `combine_cells()` | a target with its own weights on the cells of a `control_type = "strata"` result |
| `first_stage()` | the first-stage object of a full result: the cells, the shares, and the uncertainty of the cells |
| `as_first_stage()` | a first-stage object from the cells of another estimator, with the influence functions of its units or a covariance matrix |
| `double_did()` | the second stage and the aggregation on a first-stage object, or custom weights on its cells and means |


## Where the documentation lives

- `agent/STATE.md` — what the package is: the key parameters of `fastdid()`, the file map, the data flow, and the dependencies. Read it before you change anything.
- `agent/ISSUES.md` — decisions with their reasons, and traps to avoid. Read it before you edit the code that an entry names.

## Commands

There is no Makefile. Use the standard devtools workflow:

```r
roxygen2::roxygenise()    # regenerate man/ after a roxygen edit
devtools::check()         # full package check
devtools::install()       # install the development version
```

Run the tests from the package root. The suite must pass before a change is
done:

```r
suppressMessages(devtools::load_all(quiet = TRUE))
data.table::setDTthreads(0)
tinytest::run_test_dir("inst/tinytest")
```

The same suite runs one test file per worker, in about half the time. The mono
deploy script and its CI use this form:

```bash
Rscript development/run_tests.R
```

A test file runs in its own worker, so it must not need another test file.

CI runs R-CMD-check on macOS, Windows, and Ubuntu, plus coverage and pkgdown.
The workflows are in `.github/workflows/`.

## Inbox

`INBOX.md` at the repo root is a drop-box. The root CLAUDE.md of the mono repo gives its rules. Never act on it without a prompt.

## Comments

The code guide `STYLE.md` at the mono root gives the comment rules (COMMENT-1 to COMMENT-7) and the tool rules (TOOL-1 on).

## Prose

Write every document, error message, and commit message in Simplified Technical English. The user-level `~/.claude/CLAUDE.md` gives the rules.
