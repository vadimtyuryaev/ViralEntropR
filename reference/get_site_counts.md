# Get Site Counts Matrix (Internal)

Extract and tabulate amino acid counts for a specific site across all
partitions.

## Usage

``` r
get_site_counts(partitions, site_index, alphabet_size)
```

## Arguments

- partitions:

  List of data frames, one per time window. Each must have
  integer-encoded amino acid columns.

- site_index:

  Integer. Index of the site (column) to count.

- alphabet_size:

  Integer. Size of the amino acid alphabet (typically `25` for the
  package's default encoding).

## Value

Integer matrix with `alphabet_size` rows (one per amino acid code, 1
through `alphabet_size`) and one column per partition. Each cell is the
count of that amino acid at `site_index` in that partition.

## Details

Built on [`tabulate`](https://rdrr.io/r/base/tabulate.html) for fast
integer counting. Consumed internally by
[`calculate_hellinger_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_hellinger_matrix.md).

## Examples

``` r
# 1. Create dummy partitions.
# Partition 1: site 1 has three 'A's (code 1).
p1 = data.frame(s1 = c(1, 1, 1))
# Partition 2: site 1 has one 'A' (1) and two 'R's (2).
p2 = data.frame(s1 = c(1, 2, 2))
parts = list(p1, p2)

# 2. Get counts for site 1.
# Internal function — accessed via the triple-colon operator.
counts = ViralEntropR:::get_site_counts(parts, site_index = 1,
                                         alphabet_size = 25)

# Row 1 (A) is [3, 1]; row 2 (R) is [0, 2].
print(counts[1:5, ])
#>      [,1] [,2]
#> [1,]    3    1
#> [2,]    0    2
#> [3,]    0    0
#> [4,]    0    0
#> [5,]    0    0
```
