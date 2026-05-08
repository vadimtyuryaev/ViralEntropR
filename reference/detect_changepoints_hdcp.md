# Detect Temporal Change Points (HDcpDetect)

Runs high-dimensional change point detection on time-series matrices
(typically Hellinger distance matrices) over one or more time windows,
using either binary segmentation or wild binary segmentation.

## Usage

``` r
detect_changepoints_hdcp(
  data_matrix,
  min_window,
  max_window,
  n_timesteps,
  rolling_window = FALSE,
  wild = FALSE,
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

- wild:

  Logical. If `TRUE`, uses
  [`wild.binary.segmentation`](https://rdrr.io/pkg/HDcpDetect/man/wild.binary.segmentation.html);
  if `FALSE` (default), uses
  [`binary.segmentation`](https://rdrr.io/pkg/HDcpDetect/man/binary.segmentation.html).
  See Details for the methodological distinction.

- ...:

  Additional arguments passed to the chosen segmentation function. Note
  that `M` (number of random intervals) is only accepted by
  `wild.binary.segmentation`; passing it with `wild = FALSE` produces an
  "unused argument" error.

## Value

A named list:

- Points_List:

  List of integer vectors `c(start, end)` giving the row index window
  used at each step.

- HDcp_list:

  List of full segmentation result objects, one per step. Steps where
  `end_idx` exceeds the matrix rows return `NULL` with a warning.

## Details

The function applies
[`binary.segmentation`](https://rdrr.io/pkg/HDcpDetect/man/binary.segmentation.html)
or
[`wild.binary.segmentation`](https://rdrr.io/pkg/HDcpDetect/man/wild.binary.segmentation.html)
repeatedly to slices of `data_matrix`, advancing the slice forward by
one row at each step. A total of `n_timesteps + 1` detections are
performed: one on the initial window `[min_window, max_window]`, then
`n_timesteps` additional detections on subsequent windows.

Two window-advancement modes are supported, mirroring
[`detect_changepoints_ecp`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_ecp.md):

- **Expanding** (`rolling_window = FALSE`, default): the start index is
  fixed at `min_window`; the end index grows by 1 each step. Each
  iteration uses one more row than the previous. Natural for online
  surveillance, where change-point detection is re-run as new data
  accumulates.

- **Rolling** (`rolling_window = TRUE`): both start and end indices
  advance by 1 each step, keeping the window length constant. Natural
  for retrospective windowed analysis.

**Choice of segmentation method.** Binary segmentation (`wild = FALSE`,
default) is the classical recursive method: it finds the most likely
change point, splits the series, and recurses. Wild binary segmentation
(`wild = TRUE`) is the variant of Fryzlewicz (2014) that draws random
subintervals to improve detection of multiple closely-spaced change
points; it accepts an additional `M` argument controlling the number of
random intervals.

## Warnings and errors

- [`stop()`](https://rdrr.io/r/base/stop.html):

  Triggered if `max_window` exceeds `nrow(data_matrix)`: the initial
  window cannot be constructed and no detection can run.

- [`warning()`](https://rdrr.io/r/base/warning.html):

  Triggered (and execution continues) if `max_window + n_timesteps`
  exceeds `nrow(data_matrix)`, in which case the offending later
  iterations are skipped and `NULL` entries appear in `HDcp_list`.

## See also

[`detect_changepoints_ecp`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_ecp.md)
for the energy-statistic alternative;
[`calculate_hellinger_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_hellinger_matrix.md)
for the typical upstream input.

## Examples

``` r
set.seed(42)
baseline = matrix(rnorm(50, mean = 0, sd = 0.1), nrow = 10, ncol = 5)
variant  = matrix(rnorm(50, mean = 3, sd = 0.1), nrow = 10, ncol = 5)
data_mat = rbind(baseline, variant)

# Single-window detection with the default binary segmentation.
res = detect_changepoints_hdcp(
  data_matrix = data_mat,
  min_window  = 1,
  max_window  = 20,
  n_timesteps = 0
)
print(res$HDcp_list[[1]])
#>      FoundList pvalues
#> [1,]        10       0

# Wild binary segmentation with a custom number of random intervals.
res_wild = detect_changepoints_hdcp(
  data_matrix = data_mat,
  min_window  = 1,
  max_window  = 20,
  n_timesteps = 0,
  wild        = TRUE,
  M           = 100
)
print(res_wild$HDcp_list[[1]])
#> [1] 10
```
