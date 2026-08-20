# CLAUDE.md

Guidance for Claude Code (claude.ai/code) in this repository.

fastdid is an R package. It implements the staggered
difference-in-differences estimator of Callaway and Sant’Anna (2021). It
adds two extensions: time-varying covariates (Caetano and Callaway,
2024) and multiple treatment events (Tsai, 2026). It is on CRAN. The
current version is in `DESCRIPTION`.

The package trades memory for speed. It sorts once, then reads the data
through fast-access structures. On millions of units it cuts estimation
from hours to seconds.

## Public API

fastdid exports three functions:

| Function | What it does |
|----|----|
| [`fastdid()`](https://tsailintung.github.io/fastdid/reference/fastdid.md) | the main estimator |
| [`plot_did_dynamics()`](https://tsailintung.github.io/fastdid/reference/plot_did_dynamics.md) | an event-study plot through ggplot2 |
| [`sim_did()`](https://tsailintung.github.io/fastdid/reference/sim_did.md) | simulates a staggered DiD dataset for testing |

Key
[`fastdid()`](https://tsailintung.github.io/fastdid/reference/fastdid.md)
parameters:

| Parameter | Meaning |
|----|----|
| `data` | a data.table, or an object that coerces to one |
| `timevar`, `cohortvar`, `unitvar`, `outcomevar` | column names, as strings |
| `control_option` | `"both"`, `"never"`, or `"notyet"` |
| `result_type` | `"group_time"`, `"time"`, `"group"`, `"simple"`, `"dynamic"`, `"group_group_time"`, or `"dynamic_stagger"` |
| `control_type` | `"ipw"`, `"reg"`, or `"dr"` (doubly robust) |
| `cohortvar2` | the second treatment, for double DiD |
| `boot` | bootstrap standard errors. The default is 1000 iterations |
| `base_period` | `"universal"` or `"varying"` |
| `covariatesvar`, `varycovariatesvar` | time-invariant and time-varying covariates |
| `anticipation`, `anticipation2` | anticipation periods |
| `exper` | experimental features, such as `filtervar` and `aggregate_scheme` |
| `full` | returns the influence functions and the weights as well |
| `parallel` | Unix parallelization through `mclapply` |

## Internal architecture

### File map

    R/
    ├── fastdid.R          entry point: validate, coerce, estimate, aggregate
    ├── validate.R         input validation through dreamerr
    ├── aux_funcs.R        data coercion, aux data structures, locked lists
    ├── estimate_gtatt.R   core loop over cohorts and times
    ├── estimate_did.R     the 2x2 engine: IPW, OR, doubly robust, influence functions
    ├── aggregate_gt.R     aggregates g-t ATTs to the target parameters, and their SEs
    ├── double_did.R       multiple events / double DiD
    ├── sim_did.R          data simulation
    ├── generics.R         S3 methods for the fastdid_result class
    └── global.R           global-variable suppressions for R CMD check

### Data flow

    fastdid()
      -> validate_argument() + validate_dt()
      -> coerce_dt()        normalize time and cohort to 1,2,3...; store time_offset, time_step
      -> get_auxdata()      fast-access structures, clustered by time period
      -> estimate_gtatt()   per outcome: cohorts x times, then estimate_did()
      -> aggregate_gt()     weight the g-t estimates to the targets, then the SEs
      -> convert_targets()  map back to the original time and cohort scales

### Key design decisions

1.  **Sort once.** Coercion sorts the data by `(time, G, unit)`. Every
    later access is sequential and fast. Do not reorder the data after
    coercion.
2.  **Two estimation paths.** `estimate_did_bp()` handles a balanced
    panel. `estimate_did_rc()` handles an unbalanced one.
3.  **Influence functions.** The package computes a full semiparametric
    influence function at the 2x2 level, then aggregates it. The
    bootstrap runs through
    [`BMisc::multiplier_bootstrap`](https://bcallaway11.github.io/BMisc/reference/multiplier_bootstrap.html).
4.  **Locked lists.** A parameter object carries the `locked` S3 class.
    `$<-`, `[<-` and `[[<-` all raise an error, so nothing mutates it.
5.  **Time normalization.** User time (for example 2000, 2005, 2010)
    becomes 1, 2, 3 internally. `time_offset` and `time_step` recover
    the original scale.
6.  **Double DiD.** When `cohortvar2` is set, the estimator decomposes
    the effects across two staggered events. It uses three cases, chosen
    by the relative timing of g1 and g2. `double_did.R` implements
    Theorem 3 of Tsai (2026).
7.  **IPW caching.** When several outcomes share the same covariates,
    the propensity scores are computed once and reused.

## Dependencies

- **Imports**: `data.table (>= 1.15.0)`, `stringr`, `BMisc`, `collapse`,
  `dreamerr (>= 1.4.0)`, `parglm`, `ggplot2`
- **Suggests**: `did`, `knitr`, `rmarkdown`, `parallel`, `tinytest`
- **R**: \>= 4.1.0

## Commands

There is no Makefile. Use the standard devtools workflow:

``` r

roxygen2::roxygenise()    # regenerate man/ after a roxygen edit
devtools::check()         # full package check
devtools::install()       # install the development version
```

Run the tests from the package root. The suite must pass before a change
is done:

``` r

suppressMessages(devtools::load_all(quiet = TRUE))
data.table::setDTthreads(0)
tinytest::run_test_dir("inst/tinytest")
```

CI runs R-CMD-check on macOS, Windows, and Ubuntu, plus coverage and
pkgdown. The workflows are in `.github/workflows/`.

## Inbox

`INBOX.md` at the repo root is a drop-box. Humans and agents at work in
other repos leave todo items or information there. Never act on it
without a prompt. Process entries only when the user asks, for example
“resolve the issues in INBOX.md”.

To resolve an entry, move anything durable into this file, then delete
the entry. To leave an entry, append it with a date and a source. Make
it complete enough to act on without the conversation that produced it.

## Comments

Match the density of the code around you. About one line in eight
carries a comment here.

- **A function gets a roxygen block** (`#'`): what it does, what it
  takes, and what it returns. In a package, end an internal function’s
  block with `@noRd`, so that roxygen generates no `man/` page for it.
- **A chunk that does one thing gets a short label above it**, even when
  the chunk is not worth a function of its own. A few words is enough:
  `# filters`, `# create daily panel`.
- **Explain the goal or the behavior where the code does not declare
  it**: a magic code or constant, a quirk of the raw data, a workaround,
  an ordering that matters, or a line that looks wrong and is not.
- Use `# name ----` dividers to structure a long file.
- Do not restate what the line already says.
- Do not leave commented-out code. Git history keeps it.
- Do not record history: no “this used to…”, no dated incident, no
  `@created`.
- Keep an open question out of the code. Put it in `TODO.md` and leave a
  one-line pointer where it applies.

## Prose

Write every document, error message, and commit message in Simplified
Technical English (ASD-STE100).

Classify the text first:

- Procedural text tells the reader what to do. Use the imperative. Use
  one instruction per sentence. Use a maximum of 20 words per sentence.
- Descriptive text explains. Use simple tenses. Use one topic per
  paragraph. Use a maximum of 25 words per sentence, and six sentences
  per paragraph.

Then apply these rules:

- Use these verb forms only: infinitive, imperative, simple present,
  simple past, simple future, and past participle as an adjective.
- Do not use the present perfect. Do not use “-ing” verb forms.
- Use the active voice. Use the passive voice only in descriptive text,
  when the agent is unknown.
- Use these modals only: can, will, must. Do not use should, would, may,
  might, or could.
- Keep the articles and the word “that”. Do not use contractions. Do not
  use semicolons.
- Put the condition before the command: “If the test fails, read the
  log.”
- Use a vertical list for more than two items or steps.
- Keep one meaning for one word through the whole document.
- Delete words that carry no fact: simply, robust, powerful,
  comprehensive, leverage, “in order to”.
- Use American spelling.

Do not apply these rules to code, or to quoted output.
