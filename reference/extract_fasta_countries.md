# Extract Countries from FASTA Sequence Names

Extracts country names from the sequence name strings of an
`AAStringSet` object loaded via
[`readAAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html).
Handles single-word (e.g. `UK`), hyphenated (e.g. `Timor-Leste`), and
multi-word (e.g. `United States of America`) country names.

## Usage

``` r
extract_fasta_countries(sequence, position, problematic_characters = FALSE)
```

## Arguments

- sequence:

  An `AAStringSet` object.

- position:

  Integer (1–4). Location of the country field within the sequence name
  string:

  - **1** — text before the first `|` (e.g. `SouthKorea|...`).

  - **2** — text between the first and second `|`.

  - **3** — text between the first and second `/`.

  - **4** — text after the last `|`.

- problematic_characters:

  Logical. If `TRUE`, sequence names are re-encoded to UTF-8, replacing
  non-representable bytes with their escaped form. Useful for FASTA
  files with non-ASCII characters in headers. Default is `FALSE`.

## Value

A named list with three elements:

- countries:

  Character vector of extracted country strings, one per sequence. `NA`
  where extraction failed (no match against the chosen pattern).

- message:

  A single character string summarising extraction success.

- missing_id:

  Integer vector of indices where extraction failed, or `NA` if all
  extractions succeeded.

## Details

The function selects one of four regex patterns based on `position` and
applies it to each sequence name via
[`str_extract`](https://stringr.tidyverse.org/reference/str_extract.html).
**Only the first match per header is returned.** If a header contains
multiple delimited fields, the country must be in the first such field
for the corresponding `position` value to extract it correctly. For
example, with a GISAID-style header
`Spike|hCoV-19/USA/OH/.../2021|2021-05-15|EPI_ISL_...|`, `position = 3`
(between slashes) returns `USA`, but `position = 2` (between pipes)
returns `hCoV-19/USA/OH/...`, not `USA`. Inspect representative headers
with `names(sequence)[1]` before choosing `position`.

**Encoding.** FASTA files with non-ASCII characters in headers (accented
characters, byte-order marks, etc.) can break regex extraction. Setting
`problematic_characters = TRUE` re-encodes headers to UTF-8 with
non-representable bytes escaped, allowing the regex to proceed.

## See also

[`extract_fasta_dates`](https://vadimtyuryaev.github.io/ViralEntropR/reference/extract_fasta_dates.md)
for the date-extraction companion;
[`readAAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html)
for loading the input `AAStringSet`.

## Examples

``` r
# \donttest{
path_sample  <- system.file("extdata", "sarscov2_sample.fasta.gz",
                             package = "ViralEntropR")
fasta_sample <- Biostrings::readAAStringSet(path_sample)

# Inspect header structure to confirm field positions before extraction.
sample(names(fasta_sample), 1)
#> [1] "UAB29556.1 |India|2021-07-18"

# Extract countries (position 2 = between first and second pipe).
result <- extract_fasta_countries(fasta_sample, position = 2)
result$message
#> [1] "All countries have been extracted"
sort(table(result$countries), decreasing = TRUE)
#> 
#>          USA    Australia        India  New Zealand   Bangladesh        Chile 
#>           84            7            3            2            1            1 
#>      Germany Saudi Arabia 
#>            1            1 
# }
```
