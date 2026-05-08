# Simulate Viral Variant Evolution

Generates a synthetic, ground-truth-traceable time series of viral
amino-acid sequences in which a reference strain is progressively
challenged by emerging variants under stochastic growth and pairwise
competition. Designed as a controllable signal generator for
benchmarking the entropy, Hellinger, and change-point components of the
ViralEntropR pipeline; the dynamics are deliberately simplified rather
than a faithful biological model of evolutionary fitness.

## Usage

``` r
simulate_variant_evolution(
  ref_sequences,
  n_ref_months = 3L,
  start_date,
  end_date,
  variants_config,
  variant_intervals,
  n_new_mutations = 1L,
  mutation_rate = 2,
  mutation_rate_variability = 0.25,
  deleterious_rate = 0,
  n_deleterious_limit = 1L,
  n_sequences_total = 140L,
  ref_variability = FALSE,
  n_ref_sequences = 100L,
  prob_deletion_event = 0,
  n_rows_to_delete = 0L,
  seed = NULL
)
```

## Arguments

- ref_sequences:

  Either a single character string giving the reference amino-acid
  sequence, or a data.frame with a `Date` column and site columns named
  with integer position labels (`"1"`, `"2"`, ..., `as.character(L)`;
  one column per site). The sequence length `L` must be at least 3,
  since variant mutations and random substitution events are sampled
  from sites \\\>2\\ (sites 1 and 2 are reserved for reference
  variability when `ref_variability = TRUE`).

- n_ref_months:

  Integer. Duration of the initial reference phase (months). Default
  `3`.

- start_date:

  Character or Date. Simulation start.

- end_date:

  Character or Date. Simulation end (inclusive month).

- variants_config:

  Integer vector. Number of mutations per variant, e.g., `c(2, 6)`.

- variant_intervals:

  Integer vector. Months between consecutive variant emergences (length
  `length(variants_config) - 1`).

- n_new_mutations:

  Integer. Number of de-novo sequence *copies* introduced at a variant's
  first appearance. Despite the name, this controls the count of seeded
  sequences, not the number of mutations carried by the variant (which
  is set per-variant by `variants_config`). Default `1`.

- mutation_rate:

  Numeric scalar or vector. Monthly growth multiplier per variant.
  Default `2`.

- mutation_rate_variability:

  Numeric in `[0, 1)`. Fractional spread around `mutation_rate`: growth
  is drawn from `Uniform(mult*(1-v), mult*(1+v))`. Set to `0` for
  deterministic growth. Default `0.25`.

- deleterious_rate:

  Numeric in `[0, 1]`. Per-mutation probability of being flagged as
  "deleterious". Default `0`.

- n_deleterious_limit:

  Integer. Per-site cap on the number of variant rows allowed to retain
  a flagged ("deleterious") allele in a given month; overflow is
  reverted to the reference allele at the flagged site only. See
  **Capped "Deleterious" Sites** in Details. Default `1`.

- n_sequences_total:

  Integer. Total sequences generated over the full mutation phase
  (spread evenly per month). Default `140`.

- ref_variability:

  Logical. If `TRUE`, introduces low-level variability at sites 1 and 2
  of the reference pool. Requires `n_ref_sequences >= 6 * n_ref_months`
  for the variability rows to be correctly captured by the multi-variant
  fill pool; smaller values trigger a warning. Default `FALSE`.

- n_ref_sequences:

  Integer. When `ref_sequences` is a single string, the number of
  reference rows to create. Default `100`.

- prob_deletion_event:

  Numeric in `[0, 1]`. Per-month probability of triggering a
  substitution event affecting `n_rows_to_delete` rows. (Parameter name
  retained for backward compatibility; the operation is substitution,
  not literal deletion. See Details.) Default `0`.

- n_rows_to_delete:

  Integer. Number of rows affected when a substitution event fires (see
  `prob_deletion_event`). Default `0`.

- seed:

  Integer or `NULL`. Random seed for reproducibility. Default `NULL`.

## Value

An object of class `"viralSim"`, a named list:

- Simulation_Output:

  Data frame of all sequences, one column per site (named `"1"`, ...,
  `as.character(L)`) plus `Variant` (per-row strain label), `Phase`
  (`"R"` or `"RV"` for reference rows depending on whether the row
  matches the wild-type sequence; `"M1"`, `"M2"`, ... in mutation-phase
  months — note that `"Mn"` is a *month label* and applies to every row
  in that month regardless of strain), `Date`, `Period` (1-indexed
  month), and `Delet` (`"Yes"` for rows affected by substitution events,
  `"No"` otherwise).

- Variant_Details:

  List of per-variant metadata. Each element has `em` (emergence-month
  index), `pos` (mutated site indices), `flags` (logical vector of
  "deleterious" flags per mutated site), `vseq` (full mutated reference
  sequence as character vector), `mult` (growth multiplier), `last`
  (count in the most recent month), `cum` (cumulative count), and `repl`
  (logical flag set when the variant has been displaced by its
  successor).

- Simulation_Dates:

  Date vector of all simulated months.

- Baseline_Ref_Sequence:

  Character. The wild-type reference string.

- Delet_Records:

  Named list of substitution events keyed by date (`"YYYY-MM-DD"`). Each
  entry holds `period`, `date`, `site`, `old_aa`, `new_aa`, and `rows`
  (absolute row indices in `Simulation_Output` that were modified).

- Pool:

  Fill-pool data frame from the final multi-variant competition period
  (includes `._weight` column); `NULL` if at most one variant was ever
  active.

## Details

Simulates the trajectory of a viral population over discrete monthly
time steps. All stochastic components are reproducible under a fixed
`seed` and are designed to expose a known emergence pattern that
downstream detection methods can be evaluated against.

**Simulation Phases:**

1.  **Reference Phase:** A stable period in which only the reference
    strain exists. Reference months are resampled with replacement to
    assemble `n_ref_months` monthly batches.

2.  **Variant Emergence:** New variants are introduced at the monthly
    intervals specified by `variant_intervals`.

3.  **Pairwise Growth:** Each month the newest variant grows alongside
    its immediate predecessor. Both grow stochastically: each variant's
    next-month count is its previous-month count times
    `Uniform(mult*(1-variability), mult*(1+variability))`. The two are
    not coupled by a shared frequency constraint — only by the hard
    monthly quota
    `ceiling(n_sequences_total / (n_periods - n_ref_months))` that trims
    any overshoot. Older, non-paired variants enter the fill pool with
    sampling weights `2^(i+1)`.

4.  **Capped "Deleterious" Sites:** A Bernoulli-random subset of each
    variant's mutated positions is flagged with probability
    `deleterious_rate`. At each post-emergence month the count of
    variant rows still carrying the flagged allele is capped at
    `n_deleterious_limit` by reverting overflow sequences to the
    reference allele *at the flagged site only*. This is a per-site cap,
    not a genotype-level fitness penalty.

5.  **Random Substitution Events:** Each month independently fires a
    Bernoulli coin flip with probability `prob_deletion_event`. If it
    fires, the last `n_rows_to_delete` rows of that month's batch
    receive a single non-reference amino-acid substitution at a randomly
    chosen site (sites 1 and 2 excluded), and those rows are flagged
    with `Delet = "Yes"`. Parameter names retain "deletion" / "delete"
    for backward compatibility; the operation itself is a substitution.

## See also

[`partition_time_windows`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md),
[`calculate_hellinger_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_hellinger_matrix.md),
[`detect_changepoints_ecp`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_ecp.md),
[`detect_changepoints_hdcp`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_hdcp.md)

## Examples

``` r
ref_seq <- "MKTIIALSYIFCLVFADYKDDDDK"

sim <- simulate_variant_evolution(
  ref_sequences   = ref_seq,
  n_ref_months    = 3,
  start_date      = "2021-01-01",
  end_date        = "2021-12-01",
  variants_config = c(3, 5),
  variant_intervals = c(4),
  n_sequences_total = 50,
  mutation_rate   = 1.5,
  seed            = 123
)
head(sim$Simulation_Output)
#>   1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24   Variant
#> 1 M K T I I A L S Y  I  F  C  L  V  F  A  D  Y  K  D  D  D  D  K Reference
#> 2 M K T I I A L S Y  I  F  C  L  V  F  A  D  Y  K  D  D  D  D  K Reference
#> 3 M K T I I A L S Y  I  F  C  L  V  F  A  D  Y  K  D  D  D  D  K Reference
#> 4 M K T I I A L S Y  I  F  C  L  V  F  A  D  Y  K  D  D  D  D  K Reference
#> 5 M K T I I A L S Y  I  F  C  L  V  F  A  D  Y  K  D  D  D  D  K Reference
#> 6 M K T I I A L S Y  I  F  C  L  V  F  A  D  Y  K  D  D  D  D  K Reference
#>   Phase       Date Period Delet
#> 1     R 2021-01-01      1    No
#> 2     R 2021-01-01      1    No
#> 3     R 2021-01-01      1    No
#> 4     R 2021-01-01      1    No
#> 5     R 2021-01-01      1    No
#> 6     R 2021-01-01      1    No
```
