# Cluster a Univariate Numeric Vector by Gaussian Mixture Model

Wraps
[`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.html)
for unsupervised clustering of a univariate numeric vector, with
preprocessing rules and edge-case handling tailored to per-site Shannon
entropy values from viral sequence data, which is the package's primary
use case, but applicable to any univariate data the user wishes to
cluster by GMM.

## Usage

``` r
cluster_sites_by_entropy(
  entropies,
  nr,
  nsites = length(entropies),
  precision = 6L,
  removez = TRUE,
  removesngl = TRUE,
  transfr = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- entropies:

  Numeric vector to cluster. In the package's primary use case these are
  per-site Shannon entropy values, but any univariate numeric vector is
  accepted.

- nr:

  Integer. Total number of sequences from which the entropies were
  computed. Required only when `removesngl = TRUE`.

- nsites:

  Integer. Expected number of sites. If it mismatches
  `length(entropies)`, the actual length is used with a warning. Default
  is `length(entropies)`.

- precision:

  Integer. Decimal places for rounding during singleton threshold
  comparison and the all-identical uniqueness check. Default is `6`.

- removez:

  Logical. If `TRUE`, removes sites with entropy = 0 (invariant sites),
  using a small tolerance `1e-9` to absorb floating-point near-zeros.
  Default is `TRUE`.

- removesngl:

  Logical. If `TRUE`, removes sites whose entropy equals the singleton
  value (one differing sequence out of `nr`). Uses tolerance-based
  comparison. Default is `TRUE`.

- transfr:

  A function, or an object of class `transform` with a `$transform()`
  method, applied to entropies before clustering. Default is `NULL` (no
  transformation).

- verbose:

  Logical. If `TRUE`, emits diagnostic warnings for non-fatal events
  (empty partitions, Mclust failures, etc.). Default is `FALSE`.

- ...:

  Additional arguments passed to
  [`Mclust`](https://mclust-org.github.io/mclust/reference/Mclust.html).

## Value

A named list with two elements:

- FitObject:

  The raw `Mclust` result, or a minimal
  `list(classification = integer(0L))` when clustering was bypassed or
  failed.

- DataFrame:

  A data frame with columns `sites` (original site indices), `entropies`
  (values after any transformation), and `class` (GMM cluster label).
  The `class` column is **always** present in every return path,
  including zero-row DataFrames. Downstream consumers need only guard on
  `nrow(df) > 0` before accessing class values. Raw Mclust labels are
  returned as-is; call
  [`relabel_entropy_classes`](https://vadimtyuryaev.github.io/ViralEntropR/reference/relabel_entropy_classes.md)
  to standardise label ordering.

## Details

In the package's typical use, sites are clustered by their Shannon
entropy to identify groups of residue positions with similar variability
across a sequence collection. Two preprocessing rules apply when
clustering entropies: `removez = TRUE` drops invariant sites (entropy =
0), and `removesngl = TRUE` drops singleton sites whose entropy
corresponds to exactly one differing sequence across `nr` rows.

Class assignment rules (applied in priority order):

- **No rows remaining after filtering**: empty DataFrame returned with a
  zero-length `class` column (consistent schema).

- **Single row remaining**: class `1` assigned directly; Mclust is not
  called (undefined on 1 observation).

- **All entropies identical**: class `999` for all sites (sentinel — one
  undifferentiated group).

- **Normal Mclust result**: raw class labels `1, 2, ..., G`. These are
  Mclust's own integer labels, ordered by increasing component mean
  (univariate Mclust orders components by mean) — call
  [`relabel_entropy_classes()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/relabel_entropy_classes.md)
  on the returned data frame to obtain application-friendly class labels
  (highest-entropy class `= 1`, lowest-entropy class `= G`).
  [`relabel_entropy_classes`](https://vadimtyuryaev.github.io/ViralEntropR/reference/relabel_entropy_classes.md)
  on the returned DataFrame to standardise so that class 1 =
  highest-entropy group.

- **Mclust failure**: empty DataFrame returned (same schema as the
  no-rows case), treating the partition as uninformative.

## See also

[`calculate_entropy`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_entropy.md)
for computing per-site entropy values,
[`relabel_entropy_classes`](https://vadimtyuryaev.github.io/ViralEntropR/reference/relabel_entropy_classes.md)
for standardising the returned class labels, and
[`partition_time_windows`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md),
which calls this function on each temporal partition.

## Examples

``` r
# Clear bimodal structure: 5 low-entropy + 5 high-entropy sites.
set.seed(42)
entropies <- c(rnorm(5, mean = 0.1, sd = 0.01),
               rnorm(5, mean = 1.5, sd = 0.1))
result <- cluster_sites_by_entropy(entropies, removez = FALSE,
                                   removesngl = FALSE)
print(result$DataFrame)
#>    sites  entropies class
#> 2      2 0.09435302     1
#> 3      3 0.10363128     2
#> 5      5 0.10404268     2
#> 4      4 0.10632863     2
#> 1      1 0.11370958     1
#> 6      6 1.48938755     4
#> 8      8 1.49053410     4
#> 10    10 1.49372859     4
#> 7      7 1.65115220     5
#> 9      9 1.70184237     5

# Single-row edge case: class = 1 assigned directly.
res1 <- cluster_sites_by_entropy(0.35, removesngl = FALSE)
print(res1$DataFrame)
#>   sites entropies class
#> 1     1      0.35     1

# All-identical edge case: class = 999 (sentinel, one undifferentiated group).
res2 <- cluster_sites_by_entropy(c(0.35, 0.35, 0.35), removesngl = FALSE)
print(res2$DataFrame)
#>   sites entropies class
#> 1     1      0.35   999
#> 2     2      0.35   999
#> 3     3      0.35   999
```
