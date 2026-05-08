# Decode Amino Acid Sequences

Converts an integer-encoded matrix of amino acids back to its character
representation under the package's 25-symbol alphabet. Inverse of
[`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md).

## Usage

``` r
decode_aa_sequence(matrix_input)
```

## Arguments

- matrix_input:

  Numeric matrix of integer-encoded amino acids, typically the output of
  [`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md).
  Values outside `1:25` (including `0` and `NA`) decode to the sentinel
  string `"0"`. A non-matrix input is coerced via
  [`as.matrix`](https://rdrr.io/r/base/matrix.html).

## Value

A character matrix of the same dimensions as `matrix_input`, with the
same `dimnames`. Each cell is either a one-character amino acid code
from the 25-symbol alphabet, or the sentinel `"0"` for out-of-range and
missing values.

## Details

Decoding is a single vectorised lookup against the fixed alphabet
(`A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V` for the
twenty standard residues, then `B, Z, X, *, -` for ambiguous codes, stop
codons, and gaps). The input matrix is flattened, indexed against the
alphabet vector in one operation, and reshaped — there is no per-row or
per-element loop.

Values outside the valid range `1:25` (including `0`, which
[`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md)
produces for unrecognised characters) are returned as the sentinel
string `"0"`. This preserves round-trip consistency with the encoder:
encoding then decoding any character originally outside the alphabet
yields `"0"` rather than throwing an error. `NA` values are also mapped
to `"0"`.

Row and column names of the input matrix are preserved on the output.

## See also

[`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md)
for the inverse operation;
[`fasta_to_char_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/fasta_to_char_matrix.md)
for the FASTA-to-character-matrix step that typically precedes encoding.

## Examples

``` r
# 1. Decode a numeric matrix.
num_mat = matrix(c(1, 2, 25, 10), nrow = 2, byrow = TRUE)
decoded = decode_aa_sequence(num_mat)
print(decoded)
#>      [,1] [,2]
#> [1,] "A"  "R" 
#> [2,] "-"  "I" 
# 1 -> "A", 2 -> "R", 25 -> "-", 10 -> "I"

# 2. Round-trip consistency check (excluding unknowns).
orig = matrix(c("A", "C", "W", "G"), nrow = 2)
enc  = encode_aa_sequence(orig)
dec  = decode_aa_sequence(enc)
all.equal(orig, dec)
#> [1] TRUE

# 3. Out-of-range and NA values both decode to the sentinel "0".
decode_aa_sequence(matrix(c(0, NA, 30, 5), nrow = 2))
#>      [,1] [,2]
#> [1,] "0"  "0" 
#> [2,] "0"  "C" 
```
