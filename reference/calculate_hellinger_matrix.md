# Calculate Hellinger Distance Matrix

Computes the Hellinger distance between the amino acid distribution of a
reference time point (first partition) and all subsequent time points,
for each requested site.

## Usage

``` r
calculate_hellinger_matrix(
  partitions,
  sites = seq_len(ncol(partitions[[1]])),
  aa_levels = 25L,
  normalized = FALSE,
  labels = paste0("T", seq_along(partitions)),
  include_freq_tables = FALSE
)
```

## Arguments

- partitions:

  A list of data frames, one per time window. Each data frame must have
  numeric-encoded amino acid sequences as columns (integers 1 to
  `aa_levels`). This is typically the `$Partitions` element returned by
  [`partition_time_windows`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md).

- sites:

  Integer vector. Indices of the sites to analyse. Defaults to all sites
  (`seq_len(ncol(partitions[[1]]))`).

- aa_levels:

  Integer. Alphabet size. Must match the encoding used when the
  partitions were created:
  [`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md)
  produces values 1–25 by default (20 standard residues, three ambiguous
  codes, `*`, `-`). Default is `25L`.

- normalized:

  Logical. If `TRUE`, scales distances by \\1/\sqrt{2}\\ to bound the
  result in \\\[0, 1\]\\. Default is `FALSE`.

- labels:

  Character vector of partition labels. Length must equal
  `length(partitions)`. The first label names the reference partition;
  subsequent labels become column names of the returned matrix. Defaults
  to `"T1"`, `"T2"`, …

- include_freq_tables:

  Logical. If `TRUE`, the return value also includes raw count tables
  and proportion tables for each site. Default is `FALSE`.

## Value

When `include_freq_tables = FALSE` (default): a numeric matrix with rows
corresponding to `sites` and columns corresponding to `labels[-1]`, each
entry being the Hellinger distance from the reference partition
(`labels[1]`) at that site. Row names are the site indices; column names
are taken from `labels[-1]`. With default labels, the reference
partition is `"T1"` and the matrix has columns `"T2"`, `"T3"`, ….

When `include_freq_tables = TRUE`: a named list with elements `Sites`,
`Hellinger_Distances`, `Frequency_Tables`, and `Proportions_Tables`.

## Details

The Hellinger distance between two discrete distributions \\P\\ and
\\Q\\ is: \$\$H(P, Q) = \sqrt{\sum\_{i=1}^{k} (\sqrt{p_i} -
\sqrt{q_i})^2}\$\$ When `normalized = TRUE` the result is scaled by
\\1/\sqrt{2}\\, bounding the distance to \\\[0, 1\]\\. Otherwise the
range is \\\[0, \sqrt{2}\]\\.

Internally, amino acid counts per site per partition are tabulated using
[`get_site_counts`](https://vadimtyuryaev.github.io/ViralEntropR/reference/get_site_counts.md)
(built on [`tabulate`](https://rdrr.io/r/base/tabulate.html)).
Per-partition proportions and Hellinger distances are then computed by
fully vectorised matrix operations — no inner loop over partitions.

## See also

[`partition_time_windows`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md)
for producing temporally partitioned data,
[`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md)
for integer encoding consistent with `aa_levels`, and
[`detect_changepoints_ecp`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_ecp.md)
or
[`detect_changepoints_hdcp`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_hdcp.md)
for downstream change-point detection on the returned matrix.

## Examples

``` r
# Toy 3-partition data: site 1 starts homogeneous (all Alanine, 1),
# acquires Valine (20) across two later partitions; site 2 stays
# constant (all Valine throughout).
p1 = data.frame(s1 = c(1, 1, 1, 1, 1),  s2 = c(20, 20, 20, 20, 20))
p2 = data.frame(s1 = c(20, 20, 20, 20, 20), s2 = c(20, 20, 20, 20, 20))
p3 = data.frame(s1 = c(1, 1, 20, 20, 20), s2 = c(20, 20, 20, 20, 20))
parts = list(T1 = p1, T2 = p2, T3 = p3)

result = calculate_hellinger_matrix(parts, sites = 1:2)
print(result)
#>         T2        T3
#> 1 1.414214 0.8573733
#> 2 0.000000 0.0000000
# Site 1 has nonzero distance in T2 and T3 (composition shift).
# Site 2 has zero distance throughout (constant).

# With raw frequency tables for inspection.
result2 = calculate_hellinger_matrix(parts, sites = 1:2,
                                      include_freq_tables = TRUE)
print(result2$Frequency_Tables[[1]])
#>       T1 T2 T3
#>  [1,]  5  0  2
#>  [2,]  0  0  0
#>  [3,]  0  0  0
#>  [4,]  0  0  0
#>  [5,]  0  0  0
#>  [6,]  0  0  0
#>  [7,]  0  0  0
#>  [8,]  0  0  0
#>  [9,]  0  0  0
#> [10,]  0  0  0
#> [11,]  0  0  0
#> [12,]  0  0  0
#> [13,]  0  0  0
#> [14,]  0  0  0
#> [15,]  0  0  0
#> [16,]  0  0  0
#> [17,]  0  0  0
#> [18,]  0  0  0
#> [19,]  0  0  0
#> [20,]  0  5  3
#> [21,]  0  0  0
#> [22,]  0  0  0
#> [23,]  0  0  0
#> [24,]  0  0  0
#> [25,]  0  0  0
```
