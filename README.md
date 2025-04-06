# WAreg
R package for paper "While-alive regression analysis of composite survival endpoints"
**WAreg** is an R package for modeling weighted composite endpoints using *while-alive* frequency-based methods. It is especially useful in analyzing recurrent and terminal events in clinical trials, supporting both individual- and cluster-randomized designs with time-varying effects modeled via splines.

## Key Features

- Joint modeling of recurrent and terminal events using weighted composite scores.
- Flexible spline-based modeling of time-varying effects.
- Support for cluster-randomized trials.
- Cross-validation for optimal spline complexity selection.
- Variance estimation via estimating equations.
- Visualization of time-varying effects with confidence bands.

---

## Installation

Currently, this package is not on CRAN. You can install it manually by sourcing the `.R` files or using `devtools` if it’s hosted on GitHub.

```r
# Example (assuming source is available locally or via GitHub)
# devtools::install_github("yourname/WAreg")
```

---

## Usage

### Basic Workflow

```r
library(WAreg)

# Fit model
fit <- WA_fit(
  model = Surv(time, type) ~ Z1 + Z2,
  data = your_data,
  id = "Subject",
  cluster = "Cluster",         # Optional
  wr = c(1, 1),                # Weights for recurrent events
  wd = 1,                      # Weight for terminal event
  cens_mod = Surv(time, delta) ~ Z2,
  sp_knots = c(5, 10, 15, 20), # Internal spline knots
  degree = 1,                  # Degree of splines
  cv = FALSE                   # Cross-validation (optional)
)

# Summary of model
summary(fit)

# Plot time-varying effects
plot(fit, ts = seq(0, 40, by = 1))
```

---

## Functions

### `WA_fit()`

Fits the weighted composite endpoint model with while-alive frequency integration.

**Arguments:**
- `model`: Survival formula (`Surv(time, type) ~ covariates`)
- `data`: Input data frame
- `id`: Subject ID column
- `cluster`: Cluster ID column (optional)
- `wr`, `wd`: Event weights
- `link`: Link function for loss rate
- `cens_mod`: Formula for censoring model
- `sp_knots`: Internal spline knots
- `degree`: Degree of spline
- `cv`: Logical, use cross-validation
- `K`, `tr`, `nk`: Cross-validation parameters
- `ts`: Time sequence for evaluation

### `summary.WA_obj()`

Prints spline configuration and parameter estimates with standard errors, z-values, and p-values.

### `plot.WA_obj()`

Plots smoothed time-varying covariate effects with confidence bands.

---

## Output

A `WA_obj` object containing:
- Estimated spline coefficients and variances
- Time-varying covariate effect estimates
- Plot-ready result data frames

---

## Example Data Format

Your dataset should include:
- `time`: Follow-up time
- `type`: Event type (0 for censoring, 1+ for recurrent events, max for terminal)
- Covariates: e.g., `Z1`, `Z2`
- `delta`: Censoring indicator (0 = censored, 1 = observed)

---

## Dependencies

- `dplyr`
- `survival`
- `nleqslv`
- `ggplot2`
- `tidyr`

---

## Citation

If you use **WAreg** in published work, please cite the associated methodology or this package appropriately (citation info to be provided).

---

## License

MIT License (or your choice)
