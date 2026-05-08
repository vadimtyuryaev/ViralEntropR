# Relabel Entropy Classes

Relabels GMM cluster labels so that class 1 is always the
highest-entropy group, class 2 the second highest, and so on.

## Usage

``` r
relabel_entropy_classes(df)
```

## Arguments

- df:

  A data frame containing clustering results. Must have columns `class`
  (integer cluster labels from
  [`cluster_sites_by_entropy`](https://vadimtyuryaev.github.io/ViralEntropR/reference/cluster_sites_by_entropy.md))
  and `entropies` (numeric Shannon entropy values per row).

## Value

The data frame with a relabeled `class` column (integer) and an added
`num_classes` column reporting the count of distinct labels present
after relabeling.

## Details

For univariate input (the package's typical use case),
[`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.html)
orders cluster components by increasing mean. Applied to per-site
Shannon entropy values, this places the highest-entropy group at the
highest label number (label `G` for `G` fitted components). This
function flips that convention so that label `1` always denotes the
highest-entropy group, which is more natural for downstream filtering
("class 1 = top") and visual labeling. Cluster identities and means are
unchanged; only the integer labels are remapped.

**Sentinel preservation.** The class label `999L` marks undifferentiated
groups (all entropies equal in
[`cluster_sites_by_entropy`](https://vadimtyuryaev.github.io/ViralEntropR/reference/cluster_sites_by_entropy.md))
and is never relabeled. If a data frame contains a mixture of real GMM
labels and one or more `999` entries, only the real labels are ranked
and overwritten; the `999` rows pass through unchanged.

**No-op return paths.** The input is returned unchanged (with at most a
`num_classes` column added) in four situations:

- Missing required columns (`class` or `entropies`) — warns and returns
  input.

- Zero rows.

- All rows already carry the sentinel `999L`.

- Only one class label is present (relabeling is a no-op).

## See also

[`cluster_sites_by_entropy`](https://vadimtyuryaev.github.io/ViralEntropR/reference/cluster_sites_by_entropy.md)
for the upstream clustering step;
[`plot_entropy_trajectories`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_entropy_trajectories.md)
and
[`plot_site_class_trajectory`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_site_class_trajectory.md)
for downstream consumers that rely on the relabeled class convention;
[`partition_time_windows`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md)
which calls `cluster_sites_by_entropy` per window.

## Examples

``` r
# Standard case: three classes, ranked by mean entropy.
df <- data.frame(
  sites     = 1:6,
  entropies = c(0.1, 0.1, 0.5, 0.5, 1.2, 1.3),
  class     = c(1L, 1L, 2L, 2L, 3L, 3L)
)
relabel_entropy_classes(df)
#>   sites entropies class num_classes
#> 1     1       0.1     3           3
#> 2     2       0.1     3           3
#> 3     3       0.5     2           3
#> 4     4       0.5     2           3
#> 5     5       1.2     1           3
#> 6     6       1.3     1           3
# Class 3 (highest mean entropy 1.25) -> 1; class 2 (0.5) -> 2;
# class 1 (0.1) -> 3.

# Sentinel preservation: 999 rows pass through unchanged even when
# mixed with real classes.
df_mixed <- data.frame(
  sites     = 1:4,
  entropies = c(0.5, 1.2, 0.3, 0.4),
  class     = c(1L, 2L, 999L, 1L)
)
relabel_entropy_classes(df_mixed)
#>   sites entropies class num_classes
#> 1     1       0.5     2           3
#> 2     2       1.2     1           3
#> 3     3       0.3   999           3
#> 4     4       0.4     2           3
# Class 2 (highest) -> 1; class 1 -> 2; class 999 stays 999.
```
