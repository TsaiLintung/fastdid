# CLAUDE.md — fastdid

## Package Overview

**fastdid** is an R package implementing fast staggered Difference-in-Differences (DiD) estimators. It is a high-performance implementation of Callaway & Sant'Anna (2021) with extensions for time-varying covariates (Caetano & Callaway 2024) and multiple treatment events (Tsai 2024). Published on CRAN (v1.0.6).

**Key selling point**: reduces estimation from hours to seconds for large datasets (millions of units).

---

## Public API

Three exported functions:

- **`fastdid()`** — main estimation function
- **`plot_did_dynamics()`** — event study plot via ggplot2
- **`sim_did()`** — simulate staggered DiD datasets for testing

### `fastdid()` key parameters

| Parameter | Description |
|-----------|-------------|
| `data` | data.table (or coercible) |
| `timevar`, `cohortvar`, `unitvar`, `outcomevar` | column name strings |
| `control_option` | `"both"`, `"never"`, `"notyet"` |
| `result_type` | `"group_time"`, `"time"`, `"group"`, `"simple"`, `"dynamic"`, `"group_group_time"`, `"dynamic_stagger"` |
| `control_type` | `"ipw"`, `"reg"`, `"dr"` (doubly robust) |
| `cohortvar2` | second treatment for double DiD |
| `boot` | bootstrap SE (1000 iterations default) |
| `base_period` | `"universal"` or `"varying"` |
| `covariatesvar`, `varycovariatesvar` | time-invariant / time-varying covariates |
| `anticipation`, `anticipation2` | anticipation periods |
| `exper` | experimental features (filtervar, aggregate_scheme, etc.) |
| `full` | return full results including influence functions and weights |
| `parallel` | Unix parallelization via `mclapply` |

---

## Internal Architecture

### File Map

```
R/
├── fastdid.R             # Entry point: validate → coerce → estimate → aggregate
├── validate.R            # Input validation via dreamerr
├── aux_funcs.R           # Data coercion, aux data structures, locked lists
├── estimate_gtatt.R      # Core loop: expand.grid(cohorts × times), lapply over (g,t)
├── estimate_did.R        # 2×2 DiD engine: IPW, OR, doubly-robust + influence functions
├── aggregate_gt.R        # Aggregate g-t ATTs to target params; SE via influence functions
├── double_did.R          # Multiple events / double DiD support
├── sim_did.R             # Data simulation
├── generics.R            # S3 methods for fastdid_result class
└── global.R              # R CMD check global variable suppressions
```

### Data Flow

```
fastdid()
  → validate_argument() + validate_dt()
  → coerce_dt()          # normalize time/cohort to 1,2,3...; store time_offset, time_step
  → get_auxdata()        # fast-access data structures clustered by time period
  → estimate_gtatt()     # per outcome: expand.grid(cohorts×times) → lapply → estimate_did()
  → aggregate_gt()       # weight g-t estimates → target params → SE
  → convert_targets()    # map back to original time/cohort scales
```

### Key Design Decisions

1. **Sort-once architecture**: data sorted by `(time, G, unit)` at coercion; fast sequential access thereafter.
2. **Two estimation paths**: separate code for balanced (`estimate_did_bp`) and unbalanced (`estimate_did_rc`) panels.
3. **Influence function framework**: full semiparametric efficiency; computed at 2×2 level, aggregated up. Bootstrap via `BMisc::multiplier_bootstrap`.
4. **Locked lists**: parameter objects wrapped in `locked` S3 class — `$<-`, `[<-`, `[[<-` all throw errors to prevent mutation.
5. **Time normalization**: arbitrary user time (e.g., 2000, 2005, 2010) normalized to 1,2,3… internally; `time_offset`/`time_step` stored for recovery.
6. **Double DiD**: when `cohortvar2` is set, decomposes treatment effects across two staggered events using three scenarios based on relative timing of g1 vs g2.
7. **IPW caching**: propensity scores cached when multiple outcomes share the same covariates.

---

## Dependencies

**Imports**: `data.table (>= 1.15.0)`, `stringr`, `BMisc`, `collapse`, `dreamerr (>= 1.4.0)`, `ggplot2`

**Suggests**: `did`, `knitr`, `rmarkdown`, `parallel`, `tinytest`

**R requirement**: >= 4.1.0

---

## Development Workflow

No Makefile. Standard R/devtools workflow:

```r
# Documentation (roxygen2)
roxygen2::roxygenise()

# Check package
devtools::check()

# Run tests
tinytest::test_package("fastdid", remove_side_effects = FALSE)

# Install dev version
devtools::install()
```

**Test framework**: `tinytest` — test files auto-discovered in `tests/`. CI via GitHub Actions (`.github/workflows/`): R-CMD-check on macOS/Windows/Ubuntu, coverage, pkgdown.
