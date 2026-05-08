# Detect Temporal Change Points (ECP)

Runs Energy Change Point detection on time-series matrices (typically
Hellinger distance matrices) over one or more time windows.

## Usage

``` r
detect_changepoints_ecp(
  data_matrix,
  min_window,
  max_window,
  n_timesteps,
  rolling_window = FALSE,
  dynamic_k = FALSE,
  ...
)
```

## Arguments

- data_matrix:

  A numeric matrix (Time x Features). Rows are time points, columns are
  features (e.g. Hellinger distances per site).

- min_window:

  Integer. Starting row index of the initial window.

- max_window:

  Integer. Ending row index of the initial window. Must not exceed
  `nrow(data_matrix)`.

- n_timesteps:

  Integer \>= 0. Number of *additional* detections to run after the
  initial window. The function performs `n_timesteps + 1` detections in
  total. Set to `0` for a single detection on the initial window.

- rolling_window:

  Logical. If `TRUE`, both window endpoints advance per step (sliding
  window of constant length). If `FALSE` (default), only the end
  advances (expanding window).

- dynamic_k:

  Logical. If `TRUE`, sets `K` to `nrow(subset) - 2L` at each step, so
  the algorithm considers the maximum number of change points permitted
  by the window length (subject to `minsize`). Useful for exploratory
  analysis where the true number of change points is unknown a priori.
  Most idiomatic when combined with `n_timesteps > 0`, where the maximum
  admissible K legitimately scales with window size, but also valid for
  single-window exploratory use. If `FALSE` (default), the value of `K`
  supplied via `...` is honoured at every step; if no `K` is supplied,
  the `ks.cp3o` default of 1 is used.

- ...:

  Additional arguments passed to
  [`ks.cp3o`](https://rdrr.io/pkg/ecp/man/ks.cp3o.html). Common choices
  include `K` (maximum number of change points, default 1) and `minsize`
  (minimum segment size, default 30). Note that `K` is overridden when
  `dynamic_k = TRUE`; a warning is emitted in that case.

## Value

A named list:

- Points_List:

  List of integer vectors `c(start, end)` giving the row index window
  used at each step.

- ECP_list:

  List of full `ks.cp3o` result objects, one per step. Inspect `$cpLoc`
  for the optimal change-point locations at each candidate count, and
  `$gofM` for the goodness-of-fit curve, to apply post-hoc filtering.

- ECP_est_list:

  List of change-point estimate vectors (the algorithm's own selection,
  equivalent to `$estimates`).

## Details

The function applies
[`ks.cp3o`](https://rdrr.io/pkg/ecp/man/ks.cp3o.html) repeatedly to
slices of `data_matrix`, advancing the slice forward by one row at each
step. A total of `n_timesteps + 1` detections are performed: one on the
initial window `[min_window, max_window]`, then `n_timesteps` additional
detections on subsequent windows.

Two window-advancement modes are supported:

- **Expanding** (`rolling_window = FALSE`, default): the start index is
  fixed at `min_window`; the end index grows by 1 each step. Each
  iteration uses one more row than the previous. Natural for online
  surveillance, where change-point detection is re-run as new data
  accumulates.

- **Rolling** (`rolling_window = TRUE`): both start and end indices
  advance by 1 each step, keeping the window length constant. Natural
  for retrospective windowed analysis, where local change-point
  structure is examined in fixed-size segments.

## Warnings and errors

- [`stop()`](https://rdrr.io/r/base/stop.html):

  Triggered if `max_window` exceeds `nrow(data_matrix)`: the initial
  window cannot be constructed and no detection can run.

- [`warning()`](https://rdrr.io/r/base/warning.html):

  Triggered (and execution continues) if (a) `max_window + n_timesteps`
  exceeds `nrow(data_matrix)`, in which case the offending later
  iterations are skipped and `NULL` entries appear in the result
  lists; (b) `dynamic_k = TRUE` and `K` is also supplied via `...`, in
  which case the user-supplied `K` is silently overridden by
  `nrow(subset) - 2L`.

## Examples

``` r
set.seed(123)
baseline = matrix(rnorm(50, mean = 0, sd = 0.1), nrow = 10, ncol = 5)
variant  = matrix(rnorm(50, mean = 3, sd = 0.1), nrow = 10, ncol = 5)
data_mat = rbind(baseline, variant)

# Single-window detection: one true change point at row 11.
res = detect_changepoints_ecp(
  data_matrix = data_mat,
  min_window  = 1,
  max_window  = 20,
  n_timesteps = 0,
  minsize     = 5
)
print(res$ECP_est_list[[1]])
#> [1] 11

# Strongest single change point: inspect $cpLoc[[1]] in the full result.
print(res$ECP_list[[1]]$cpLoc[[1]])
#> NULL
```
