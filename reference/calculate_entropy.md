# Calculate Shannon Entropy

Computes the Shannon entropy of a categorical vector.

## Usage

``` r
calculate_entropy(vctr, base = 2, precision = 6)
```

## Arguments

- vctr:

  A vector (character, factor, or integer) representing categorical
  data.

- base:

  A numeric scalar. The base of the logarithm. Default is 2.

- precision:

  Integer. The number of decimal places to round the result to. Default
  is 6.

## Value

A numeric scalar representing the entropy. Returns 0 if the vector
contains only one unique value or has length 0.

## Details

Entropy is calculated as \\H(X) = -\sum p(x) \log_b p(x)\\, where
\\p(x)\\ is the proportion of observations belonging to category \\x\\.

## Examples

``` r
seq_vec = c("A", "A", "T", "G", "C", "A")
calculate_entropy(seq_vec)
#> [1] 1.792481

# Pure homogeneity
calculate_entropy(rep("A", 10))
#> [1] 0
```
