# Encode Amino Acid Sequences

Converts a character matrix of amino acid codes to its integer-encoded
representation under the package's 25-symbol alphabet. Inverse of
[`decode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/decode_aa_sequence.md).

## Usage

``` r
encode_aa_sequence(matrix_input)
```

## Arguments

- matrix_input:

  A character matrix of amino acid codes, typically produced by
  [`fasta_to_char_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/fasta_to_char_matrix.md).
  A non-matrix input is coerced via
  [`as.matrix`](https://rdrr.io/r/base/matrix.html). Lowercase input is
  accepted; characters are uppercased before lookup.

## Value

A numeric matrix of the same dimensions as `matrix_input`, with the same
`dimnames`. Each cell is either an integer in `1:25` from the alphabet
table above, or `0` for any input not in the alphabet (including `NA`
and the empty string).

## Details

Encoding is a single vectorised lookup against a fixed named map. The
input matrix is flattened, normalised to uppercase, indexed against the
alphabet's name vector in one operation, and reshaped — there is no
per-row or per-element loop.

The 25-symbol alphabet covers the twenty standard residues
(`A=1, R=2, N=3, D=4, C=5, Q=6, E=7, G=8, H=9, I=10, L=11, K=12, M=13, F=14, P=15, S=16, T=17, W=18, Y=19, V=20`),
the three IUPAC ambiguous codes `B=21, Z=22, X=23`, the stop codon
`*=24`, and the alignment gap `-=25`.

Any character outside this alphabet — including `NA`, the empty string,
lowercase letters not in the standard set after uppercasing (e.g. `J`
for leucine/isoleucine), and non-letter symbols — maps to the sentinel
`0`. This sentinel is recognised by downstream functions:
[`filter_ambiguous_sequences`](https://vadimtyuryaev.github.io/ViralEntropR/reference/filter_ambiguous_sequences.md)
treats `0` as ambiguous (alongside `B`, `X`, `Z`) and removes affected
sequences, and
[`decode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/decode_aa_sequence.md)
round-trips it back to the string `"0"`.

Row and column names of the input matrix are preserved on the output.

## See also

[`decode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/decode_aa_sequence.md)
for the inverse operation;
[`fasta_to_char_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/fasta_to_char_matrix.md)
for the FASTA-to-character-matrix step that typically precedes encoding;
[`filter_ambiguous_sequences`](https://vadimtyuryaev.github.io/ViralEntropR/reference/filter_ambiguous_sequences.md)
for downstream removal of sequences containing ambiguous or unrecognised
residues.

## Examples

``` r
# 1. Encode a simple matrix of sequences.
seq_mat = matrix(c("A", "R", "N", "D"), nrow = 2, byrow = TRUE)
encoded = encode_aa_sequence(seq_mat)
print(encoded)
#>      [,1] [,2]
#> [1,]    1    2
#> [2,]    3    4

# 2. Gaps and unknown characters.
# '-' maps to 25; '?' is not in the alphabet and maps to the sentinel 0.
gapped_mat = matrix(c("A", "-", "G", "?"), nrow = 1)
encode_aa_sequence(gapped_mat)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1   25    8    0

# 3. Lowercase input is accepted (uppercased before lookup).
encode_aa_sequence(matrix(c("a", "r", "n", "d"), nrow = 2))
#>      [,1] [,2]
#> [1,]    1    3
#> [2,]    2    4

# 4. NA and empty string both map to 0.
encode_aa_sequence(matrix(c("A", NA, "G", ""), nrow = 2))
#>      [,1] [,2]
#> [1,]    1    8
#> [2,]    0    0
```
