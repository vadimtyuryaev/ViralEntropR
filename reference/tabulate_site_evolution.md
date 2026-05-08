# Tabulate Site Frequency Evolution

Generates a frequency or proportion table showing the amino acid
distribution at a specific site across multiple time partitions,
optionally styled with kableExtra and saved as standalone HTML.

## Usage

``` r
tabulate_site_evolution(
  partitions,
  site_index,
  labels = NULL,
  alphabet_size = 25L,
  zeros = TRUE,
  use_letters = TRUE,
  relative = FALSE,
  digits = 2L,
  col_width = "100px",
  highlight_col = NULL,
  background = "#f0f8ff",
  wrap_length = 10L,
  save = FALSE,
  save_extension = ".html",
  save_path = getwd(),
  return_table = TRUE
)
```

## Arguments

- partitions:

  A list of data frames, typically produced by
  [`partition_time_windows`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md).
  Each data frame must contain integer-encoded amino acid sequences as
  columns (values 1 to `alphabet_size`).

- site_index:

  Integer. The column index (site) to analyse.

- labels:

  Character vector. Column labels for each partition. Defaults to
  `names(partitions)`, or `"P1"`, `"P2"`, ... if unnamed.

- alphabet_size:

  Integer. Total number of possible amino acid codes. Must match the
  encoding used during integer encoding
  ([`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md)
  produces values 1-25 by default). Default is `25L`.

- zeros:

  Logical. If `TRUE` (default), fills missing counts with `0`. If
  `FALSE`, replaces zeros with `""`; see Details for the type-coercion
  implication.

- use_letters:

  Logical. If `TRUE` (default), uses decoded single-character codes from
  the 25-symbol alphabet (A, R, N, D, C, ..., V, B, Z, X, \*, -) as row
  names. If `FALSE`, uses numeric codes 1 to `alphabet_size`.

- relative:

  Logical. If `TRUE`, converts counts to proportions (column-wise
  division by partition size). Default is `FALSE`.

- digits:

  Integer. Decimal places for rounding when `relative = TRUE`. Default
  is `2L`.

- col_width:

  Character. CSS width string applied to all columns (e.g. `"100px"`).
  Default is `"100px"`.

- highlight_col:

  Integer or `NULL`. 1-based column index (relative to the data columns,
  not counting the row-name column) of a partition to highlight with
  `background`. Out-of-range values trigger a warning and no highlight
  is applied. `NULL` (default) means no highlight.

- background:

  Character. CSS background colour for the highlighted column. Default
  is `"#f0f8ff"` (light blue).

- wrap_length:

  Integer. Character width at which to wrap long column labels using
  HTML line breaks. Default is `10L`.

- save:

  Logical. If `TRUE`, saves the rendered HTML table to disk via
  [`save_kable`](https://rdrr.io/pkg/kableExtra/man/save_kable.html).
  Default is `FALSE`.

- save_extension:

  Character. File extension for the saved file (including leading dot).
  Default is `".html"`.

- save_path:

  Character. Directory in which to save. Default is
  [`getwd()`](https://rdrr.io/r/base/getwd.html).

- return_table:

  Logical. If `TRUE` (default), returns a named list with both the raw
  data frame and the styled kable object. If `FALSE`, returns only the
  styled kable.

## Value

If `return_table = TRUE`, a named list:

- table:

  The raw count (or proportion) data frame, with row names corresponding
  to amino acid codes and column names corresponding to partition
  labels.

- styled:

  The kableExtra HTML kable object.

If `return_table = FALSE`, returns only the styled kable object.

## Details

Aggregates amino acid counts per partition using
[`get_site_counts`](https://vadimtyuryaev.github.io/ViralEntropR/reference/get_site_counts.md),
optionally converts to relative frequencies, applies kableExtra styling
(column width, column highlighting, striped rows), and optionally saves
to disk. Row names are decoded amino acid codes via
[`decode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/decode_aa_sequence.md)
when `use_letters = TRUE`.

**Empty partitions.** Partitions containing no observations contribute
all-zero columns. When `relative = TRUE`, division by zero is avoided by
treating empty-column sums as 1, leaving proportions at zero. Inspect
partition sizes via `sapply(partitions, nrow)` before interpreting the
table.

**Note on `zeros = FALSE`.** Setting `zeros = FALSE` replaces numeric
zeros with empty strings (`""`) for visual clarity in the kable. This
conversion forces the underlying data frame to character storage;
numeric operations (`sum`, `mean`, etc.) will not work on the returned
`table` element. Use `zeros = TRUE` (default) if downstream numerical
use is intended.

## See also

[`partition_time_windows`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md)
for producing the typical input list of partitions;
[`get_site_counts`](https://vadimtyuryaev.github.io/ViralEntropR/reference/get_site_counts.md)
for the count-tabulation primitive;
[`decode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/decode_aa_sequence.md)
for the alphabet code mapping;
[`calculate_hellinger_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_hellinger_matrix.md)
for the related cross-partition distance calculation on the same data
shape.

## Examples

``` r
p1 = data.frame(s1 = c(1L, 1L, 1L, 1L, 2L))
p2 = data.frame(s1 = c(1L, 1L, 2L, 2L, 2L))
parts = list(T1 = p1, T2 = p2)

# Default: counts, letters, no save
tbl = tabulate_site_evolution(parts, site_index = 1)
tbl$styled
#> <table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
#>  <thead>
#>   <tr>
#>    <th style="text-align:left;">   </th>
#>    <th style="text-align:right;"> T1 </th>
#>    <th style="text-align:right;"> T2 </th>
#>   </tr>
#>  </thead>
#> <tbody>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> A </td>
#>    <td style="text-align:right;width: 100px; "> 4 </td>
#>    <td style="text-align:right;"> 2 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> R </td>
#>    <td style="text-align:right;width: 100px; "> 1 </td>
#>    <td style="text-align:right;"> 3 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> N </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> D </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> C </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> Q </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> E </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> G </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> H </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> I </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> L </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> K </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> M </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> F </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> P </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> S </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> T </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> W </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> Y </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> V </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> B </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> Z </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> X </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> * </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> - </td>
#>    <td style="text-align:right;width: 100px; "> 0 </td>
#>    <td style="text-align:right;"> 0 </td>
#>   </tr>
#> </tbody>
#> </table>

# Relative frequencies, highlight second partition
tbl2 = tabulate_site_evolution(parts, site_index = 1,
                                relative = TRUE, highlight_col = 2)
tbl2$styled
#> <table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
#>  <thead>
#>   <tr>
#>    <th style="text-align:left;">   </th>
#>    <th style="text-align:right;"> T1 </th>
#>    <th style="text-align:right;"> T2 </th>
#>   </tr>
#>  </thead>
#> <tbody>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> A </td>
#>    <td style="text-align:right;width: 100px; "> 0.8 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.4 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> R </td>
#>    <td style="text-align:right;width: 100px; "> 0.2 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.6 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> N </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> D </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> C </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> Q </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> E </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> G </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> H </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> I </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> L </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> K </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> M </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> F </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> P </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> S </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> T </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> W </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> Y </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> V </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> B </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> Z </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> X </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> * </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#>   <tr>
#>    <td style="text-align:left;width: 100px; "> - </td>
#>    <td style="text-align:right;width: 100px; "> 0.0 </td>
#>    <td style="text-align:right;background-color: rgba(240, 248, 255, 255) !important;"> 0.0 </td>
#>   </tr>
#> </tbody>
#> </table>

# Numeric codes (skip alphabet decoding)
tbl3 = tabulate_site_evolution(parts, site_index = 1, use_letters = FALSE)
rownames(tbl3$table)
#>  [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15"
#> [16] "16" "17" "18" "19" "20" "21" "22" "23" "24" "25"
```
