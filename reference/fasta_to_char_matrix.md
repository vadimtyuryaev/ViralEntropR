# Convert FASTA Object to Character Matrix

Converts an `AAStringSet` object loaded via
[`readAAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html)
into a character matrix where rows are sequences and columns are residue
positions (sites). Inverse structural transformation of
[`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md)'s
expected input shape.

## Usage

``` r
fasta_to_char_matrix(fsta)
```

## Arguments

- fsta:

  An `AAStringSet` object, typically the output of
  [`readAAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html).
  May be aligned or unaligned (see Details).

## Value

A character matrix with `length(fsta)` rows and
`max(nchar(as.character(fsta)))` columns. Each cell contains a
single-character amino acid code from the input sequences (or the gap
character `"-"` for padded positions in unaligned input). The matrix has
no row or column names; sequence names from the `AAStringSet` are not
carried over. An empty input (`length(fsta) == 0`) returns an empty
0-by-0 character matrix.

## Details

**Alignment.** The function expects an aligned `AAStringSet` — all
sequences of equal width. Unaligned input is accepted and shorter
sequences are right-padded with the gap character `"-"` to match the
longest sequence, but downstream entropy-based analysis assumes
positional homology across rows; if sequences in the input are not
biologically aligned, results from per-site computations will not be
meaningful. For unaligned input, run a multiple-sequence alignment (e.g.
`msa::msa()` or `DECIPHER::AlignSeqs()`) before calling this function.

**Performance.** Conversion is fully vectorised: all sequences are
coerced to a single character string vector, split simultaneously, and
reshaped into a matrix in one operation. No per-row loop, no
intermediate list of split sequences kept alive — substantially faster
than per-row `strsplit` on large inputs (100k+ sequences).

## See also

[`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md)
for converting the resulting character matrix to an integer-encoded
matrix;
[`filter_ambiguous_sequences`](https://vadimtyuryaev.github.io/ViralEntropR/reference/filter_ambiguous_sequences.md)
for removing rows containing ambiguous residues;
[`readAAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html)
for loading FASTA files into the input format.

## Examples

``` r
# \donttest{
# Convert the bundled sample to a character matrix.
path  <- system.file("extdata", "sarscov2_sample.fasta.gz",
                      package = "ViralEntropR")
fasta <- Biostrings::readAAStringSet(path)
mat   <- fasta_to_char_matrix(fasta)
dim(mat)
#> [1]  100 1273
mat[1:3, 1:10]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,] "M"  "F"  "V"  "F"  "L"  "V"  "L"  "L"  "P"  "L"  
#> [2,] "M"  "F"  "V"  "F"  "L"  "V"  "L"  "L"  "P"  "L"  
#> [3,] "M"  "F"  "V"  "F"  "L"  "V"  "L"  "L"  "P"  "L"  
# }
```
