# Build Distance/Entropy Matrix

Converts the output of
[`calculate_hellinger_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_hellinger_matrix.md)
or the legacy `get_hellinger_dist` / entropy list into a numeric
`sites × time_steps` matrix suitable for change point detection.

## Usage

``` r
build_distance_matrix(data_input)
```

## Arguments

- data_input:

  Either:

  - A numeric matrix (`sites × time_steps`) as returned by
    [`calculate_hellinger_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_hellinger_matrix.md).

  - A named list with elements `Sites` and one of `Hellinger_Distances`
    or `Entropies`, as returned by the legacy `get_hellinger_dist`
    function.

## Value

A numeric matrix with:

- Rows:

  Sites (1 to `max(sites)`, sparse if sites are non-contiguous).

- Columns:

  Time steps (partitions T2, T3, … or entropy periods).

Row names are character site indices. Column names are preserved from
the input where available.

## Details

In the current pipeline
[`calculate_hellinger_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_hellinger_matrix.md)
already returns a `sites × time_steps` matrix, so calling this function
is usually unnecessary — transposing that matrix directly gives `dat_t`
for [`e.agglo`](https://rdrr.io/pkg/ecp/man/e.agglo.html) or
[`ks.cp3o`](https://rdrr.io/pkg/ecp/man/ks.cp3o.html):


      hell_mat = calculate_hellinger_matrix(partitions, sites = seq_len(n_sites))
      dat_t    = t(hell_mat)   # time_steps × sites — ready for ECP

This function is retained for **backward compatibility** with code that
was written against the old `get_hellinger_dist` list format, and as a
convenience wrapper that also accepts entropy lists.

**Input formats accepted:**

1.  A numeric matrix (e.g. from `calculate_hellinger_matrix`) — returned
    as-is with site rownames normalised.

2.  The legacy named list from `get_hellinger_dist` containing `$Sites`
    and `$Hellinger_Distances` (list of vectors).

3.  An entropy list containing `$Sites` and `$Entropies` (list of
    numeric vectors), used when no Hellinger distances are present.

## See also

[`calculate_hellinger_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_hellinger_matrix.md),
[`detect_changepoints_ecp`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_ecp.md)
