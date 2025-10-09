# Fixes for "seq.default() wrong sign in 'by' argument" Error

## Summary

Added informative error checks throughout the codebase to catch and report the root cause of the `seq.default()` error with clear, actionable messages.

## Changes Made

### 1. `/R/validate.R` - Early Parameter Validation

**Added checks for:**
- `balanced_event_time` must be non-negative (>= 0)
- `anticipation` must be non-negative (>= 0)
- `anticipation2` must be non-negative (>= 0)

**Why:** Catches invalid parameter values before they cause sequence generation errors.

### 2. `/R/aggregate_gt.R` - Balanced Event Time Logic

**Added check in `get_agg_targets()` function:**
- Validates that `min_event_time <= max_event_time` before filtering targets
- Provides clear error message showing:
  - The minimum available event time in the data
  - The specified `balanced_event_time`
  - Suggested corrective action

**Why:** This is the most likely source of the error. When `balanced_event_time` is less than the minimum available event time in the data, the filtering creates an invalid range.

**Example error message:**
```
Invalid balanced_event_time: The minimum available event time (-5) is greater than 
balanced_event_time (2). Please specify a balanced_event_time >= -5 or use a smaller 
value that matches your data structure.
```

### 3. `/R/estimate_gtatt.R` - Cohort Position Functions

**Added checks in:**
- `get_control_pos()`: Validates start <= end before seq()
- `get_treat_pos()`: Validates start <= end before seq()

**Why:** Catches issues with cohort structure that could lead to invalid sequences.

### 4. `/R/aux_funcs.R` - Data Extraction Sequences

**Added checks in `get_auxdata()` function:**
- When extracting outcomes for each time period
- When extracting varying covariates for each time period

**Why:** Catches data structure issues early with clear context about which time period failed.

## Root Cause Analysis

The error `"seq.default() wrong sign in 'by' argument"` occurs when:

```r
seq(from = 10, to = 5, by = 1)  # Error: 'by' should be negative
```

In the fastdid package, this can happen when:

1. **Most Common:** `balanced_event_time` parameter is incompatible with the data
   - User specifies `balanced_event_time = 2`
   - But all cohorts only have event times starting at 5 or higher
   - Filtering `targets[targets <= 2 & targets >= 5]` creates empty/invalid range

2. **Less Common:** 
   - Time periods not properly ordered
   - Cohort structure issues
   - Data structure problems (id_size miscalculation)

## Testing Recommendations

Test with:
1. Negative `balanced_event_time` values
2. `balanced_event_time` values smaller than minimum event time in data
3. Datasets with sparse event time coverage
4. Unbalanced panels with varying time observations

## User Guidance

When users encounter the original error, they should:
1. Check their `balanced_event_time` parameter
2. Examine the event time distribution in their data
3. Ensure `anticipation` and `anticipation2` are non-negative
4. Verify time and cohort variables are properly coded
