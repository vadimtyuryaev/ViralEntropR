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
