# Partition Data into Time Windows

Splits a data frame of time-stamped sequences into discrete time
windows, computes per-site entropy and Mclust clustering for each
window.

## Usage

``` r
partition_time_windows(
  data,
  n_sites,
  window_length = 1,
  window_type = 3,
  start_date = NULL,
  end_date = NULL,
  date_format = "%b-%Y",
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  Data frame. Must contain a column named `Date` coercible to `Date`.
  Columns `1` through `n_sites` must be the numeric site columns; any
  additional columns (such as `Country` or other per-sequence metadata)
  may follow.

- n_sites:

  Integer. Number of site columns. Site columns must occupy positions
  `1` through `n_sites` of `data`; any other columns (`Date`, optionally
  `Country`, etc.) must come after.

- window_length:

  Integer. Window duration in months. Default is `1`.

- window_type:

  Integer (1, 2, or 3). Windowing strategy: `1` = Cumulative, `2` =
  Sliding, `3` = Disjoint. Default is `3`.

- start_date:

  Date or character. Window start. Defaults to earliest date in `data`.

- end_date:

  Date or character. Window end (cutoff). Defaults to latest date in
  `data`.

- date_format:

  Character. Format string for date labels. Default is `"%b-%Y"`.

- verbose:

  Logical. If `TRUE`, prints a progress bar and completion summary.
  Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`cluster_sites_by_entropy`](https://vadimtyuryaev.github.io/ViralEntropR/reference/cluster_sites_by_entropy.md).

## Value

A named list:

- Partitions:

  List of data frame chunks, one per window.

- Entropies:

  List of numeric entropy vectors, one per window.

- Clusters:

  List of clustering results from
  [`cluster_sites_by_entropy`](https://vadimtyuryaev.github.io/ViralEntropR/reference/cluster_sites_by_entropy.md),
  one per window.

- Max_Entropy:

  Numeric vector. Maximum cluster class label per window (equals number
  of clusters found; `NA` if window was empty).

- Dates_Labels:

  Character vector of window label strings.

- N_partitions:

  Integer. Total number of windows.

## Details

Three windowing strategies are supported via `window_type`. Throughout
the descriptions below, `T` denotes the number of whole months between
`start_date` and `end_date`, and `w` denotes `window_length`:

- **1 — Cumulative**: Start is fixed; end expands by `w` months each
  step. Each window includes all prior data plus one more period.
  Produces `floor(T / w)` chunks.

- **2 — Sliding**: A window of `w` months slides one month at a time.
  Produces `T - w + 1` chunks.

- **3 — Disjoint**: Consecutive non-overlapping windows of `w` months.
  Produces `floor(T / w)` chunks. Default.

For example, with `T = 12` months and `w = 2`: cumulative produces 6
chunks (each progressively larger); sliding produces 11 chunks
(overlapping, each 2 months wide); disjoint produces 6 chunks (each
exactly 2 months, non-overlapping).

**Empty windows.** When a window contains no observations, the
corresponding entries in the returned lists carry placeholder values:
`Entropies` is a zero-vector of length `n_sites`; `Clusters` is a
schema-consistent empty result matching `cluster_sites_by_entropy`'s
empty-input return; `Max_Entropy` is `NA_integer_`. Downstream consumers
can guard on `nrow(Partitions[[i]]) > 0` or `!is.na(Max_Entropy[i])`.

**Column layout requirement.** The function extracts site columns as
`data[, seq_len(n_sites), drop = FALSE]`. The site columns must
therefore occupy positions `1` through `n_sites`. `Date` (and any other
metadata columns such as `Country`) must come after the site columns. A
common error is to place `Date` first; this will produce incorrect
entropy values.

## See also

[`calculate_entropy`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_entropy.md)
for the per-site entropy computation;
[`cluster_sites_by_entropy`](https://vadimtyuryaev.github.io/ViralEntropR/reference/cluster_sites_by_entropy.md)
for the GMM clustering applied per window;
[`relabel_entropy_classes`](https://vadimtyuryaev.github.io/ViralEntropR/reference/relabel_entropy_classes.md)
for standardizing the cluster labels in downstream consumers; and
[`calculate_hellinger_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_hellinger_matrix.md)
for typical downstream use of the `Partitions` output.

## Examples

``` r
dates <- seq(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
df <- data.frame(
  s1   = 1L,
  s2   = c(rep(1L, 15), rep(2L, 15)),
  Date = rep(dates, each = 5)
)
res <- partition_time_windows(df, n_sites = 2, window_length = 2,
                              window_type = 3, verbose = TRUE)
#> Partitioning 2 disjoint windows ...
#>   |                                                          |                                                  |   0%  |                                                          |=========================                         |  50%  |                                                          |==================================================| 100%
#> Partitioning complete: 2 partitions generated (Jan-2020 to Jun-2020).
print(res$Dates_Labels)
#> [1] "Jan-2020 - Feb-2020" "Mar-2020 - Apr-2020"
print(res$Max_Entropy)
#> [1] NA  1
```
