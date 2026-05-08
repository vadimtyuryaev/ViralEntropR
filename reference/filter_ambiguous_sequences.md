# Remove Sequences Containing Ambiguous Residues

Removes rows (sequences) that contain at least one ambiguous amino acid
residue (B, J, X, or Z) — and, under integer-encoded input, any
unrecognised character — from a sequence matrix. Accepts both
integer-encoded matrices and character matrices.

## Usage

``` r
filter_ambiguous_sequences(NumMatrix, option = 1)
```

## Arguments

- NumMatrix:

  A matrix. Rows are sequences, columns are sites. Either
  integer-encoded (`option = 1`) or character (`option = 2`). Despite
  the name, character matrices are also accepted under `option = 2`.

- option:

  Integer. `1` (default) for integer-encoded matrices produced by
  [`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md);
  `2` for character matrices produced by
  [`fasta_to_char_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/fasta_to_char_matrix.md).

## Value

A named list:

- OriginalDim:

  Character string reporting the number of input sequences.

- NewDim:

  Character string reporting the number of sequences remaining after
  filtering.

- NumberAmbiguous:

  Character string reporting the number of sequences that contained at
  least one ambiguous residue.

- RangeAmbiguous:

  Character string reporting the min and max count of ambiguous residues
  per removed sequence, or `"No ambiguous sequences found"` when none
  were removed.

- DeletedSeqId:

  Integer vector of row indices that were removed. Empty integer vector
  if nothing was removed.

- FilteredMatrix:

  The filtered matrix with ambiguous rows removed, preserving the
  original column structure and storage mode.

## Details

**What is removed.** Sequences are flagged for removal if any of their
residue positions contain one of the four IUPAC ambiguous codes:

- `B` — Aspartate / Asparagine.

- `J` — Leucine / Isoleucine.

- `X` — any residue.

- `Z` — Glutamate / Glutamine.

**What is NOT removed.** Standard alignment gaps (`-`, integer code
`25`) are retained — gaps represent known absences rather than uncertain
identities and are typically positionally meaningful in aligned data.
Sequences containing only canonical 20 amino acids and gaps are kept.

**How input mode is handled.**

- `option = 1` (integer-encoded input from
  [`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md)):
  rows are removed if any cell equals `0` (unrecognised — including J,
  NA, empty, lowercase mismatches, byte-order marks, and other
  characters that fell outside the encoding alphabet), `21` (B), `22`
  (Z), or `23` (X). The `0` sentinel acts as a catch-all for anything
  not in the 25-symbol alphabet.

- `option = 2` (character input from
  [`fasta_to_char_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/fasta_to_char_matrix.md)):
  rows are removed if any cell is exactly `"B"`, `"J"`, `"X"`, or `"Z"`.
  Unrecognised characters in character input (e.g. lowercase letters,
  `NA`, empty strings) are NOT caught at this stage; encode first if you
  need that catch-all behaviour.

**Performance.** Detection is fully vectorised: a single logical matrix
comparison followed by [`rowSums`](https://rdrr.io/r/base/colSums.html)
counts ambiguous residues per sequence in one C-level call, replacing
the original row-by-row loop for a substantial speed improvement on
large matrices (100k+ rows).

## See also

[`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md)
and
[`fasta_to_char_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/fasta_to_char_matrix.md)
for producing the typical input;
[`decode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/decode_aa_sequence.md)
for inspecting the surviving sequences in character form.

## Examples

``` r
# Synthetic example: 50 sequences, 10 sites, drawn from canonical residues.
set.seed(1)
m <- matrix(sample(1:20, 500, replace = TRUE), nrow = 50, ncol = 10)
# Inject ambiguous codes into 3 specific rows: 21 (B), 23 (X), 0 (unrecognised).
m[c(3, 17, 42), sample(1:10, 3)] <- c(21, 23, 0)
result <- filter_ambiguous_sequences(m, option = 1)
cat(result$NumberAmbiguous, "\n")
#> Number of sequences containing at least one of B, X, Z or J characters is 3 
cat(result$RangeAmbiguous, "\n")
#> Number of ambiguous protein characters per sequence varies between 3 and 3 
dim(result$FilteredMatrix)
#> [1] 47 10

# Character-mode example.
chr <- matrix(c("M", "K", "T", "I", "I", "X", "K", "T", "I", "I"),
               nrow = 2, byrow = TRUE)
filter_ambiguous_sequences(chr, option = 2)$DeletedSeqId
#> [1] 2
```
