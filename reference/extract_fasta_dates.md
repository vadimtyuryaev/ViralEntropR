# Extract Dates from FASTA Sequence Names

Extracts date strings from the sequence name strings of an `AAStringSet`
object loaded via
[`readAAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html).
Several built-in date patterns are provided for the common date
conventions used in NCBI and GISAID exports; a fully custom regex can
also be supplied.

## Usage

``` r
extract_fasta_dates(
  sequence,
  option = 1,
  date_format = "%Y-%m-%d",
  custom_pattern = NULL
)
```

## Arguments

- sequence:

  An `AAStringSet` object.

- option:

  Integer (1, 2, 3, or 4). Selects the built-in pattern when
  `custom_pattern` is not supplied. See Details.

- date_format:

  Character. `strptime`-style format string used to coerce extracted
  strings to `Date`. Default is `"%Y-%m-%d"`. Use `"%Y-%m"` together
  with `option = 3`.

- custom_pattern:

  Character or `NULL`. A custom regex passed directly to
  [`str_extract`](https://stringr.tidyverse.org/reference/str_extract.html).
  When supplied, `option` is ignored. Default is `NULL`.

## Value

A named list of six elements, each aligned with the input `sequence`:

- raw_date_strings:

  Character vector of extracted date strings before any correction. `NA`
  where extraction failed.

- corrected_date_strings:

  Character vector with `-00` substrings replaced by `-01`. `NA` where
  extraction failed.

- raw_dates:

  `Date` vector coerced from `raw_date_strings`. `NA` for unparseable or
  missing strings (including any record where day = `00`).

- corrected_dates:

  `Date` vector coerced from `corrected_date_strings`. `NA` for
  unparseable or missing strings.

- message:

  Character string summarising extraction success.

- missing_id:

  Integer vector of indices where extraction failed, or `NA` if all
  extractions succeeded.

## Details

Date strings of the form `yyyy-mm-dd` are matched between pipe
characters (`|...|`) by default. Day value `00` (a common GISAID
convention indicating unknown collection day) is accepted in the raw
string and corrected to `01` before coercion to `Date`. Both raw and
corrected versions are returned, so the caller can decide how to treat
unknown-day records downstream.

**Choosing a built-in pattern.** The four options correspond to the four
most common date conventions in viral sequence repositories:

- `option = 1`: `yyyy-mm-dd` between pipes — GISAID export format, where
  the date is followed by additional pipe-delimited fields.

- `option = 2`: `yyyy-dd-mm` between pipes — some European data sources
  reverse day and month.

- `option = 3`: `yyyy-mm` between pipes — month-level resolution, useful
  when the source omits or hides the day. Pair with
  `date_format = "%Y-%m"`.

- `option = 4`: `yyyy-mm-dd` at end of header — NCBI Virus export
  format, where the collection date is the final field with no trailing
  delimiter. This is the format of the bundled `sarscov2_sample`.

For datasets where the date does not lie between pipes, supply a
`custom_pattern` matching whatever surrounding context the headers
provide.

**Coercion to Date.** When `date_format = "%Y-%m"` the function uses
[`as.yearmon`](https://rdrr.io/pkg/zoo/man/yearmon.html) for coercion so
that year-month strings are handled correctly (base `as.Date` cannot
parse `"2021-05"` alone). For all other formats, base `as.Date` with the
supplied `date_format` is used.

**Output alignment.** All six elements of the return list are the same
length as the input `sequence`. Where extraction fails for a record, the
corresponding entries are `NA`; `missing_id` lists the affected indices.

## See also

[`extract_fasta_countries`](https://vadimtyuryaev.github.io/ViralEntropR/reference/extract_fasta_countries.md)
for the country-extraction companion;
[`readAAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html)
for loading the input `AAStringSet`;
[`as.yearmon`](https://rdrr.io/pkg/zoo/man/yearmon.html) for the
year-month coercion path.

## Examples

``` r
# \donttest{
path_sample  <- system.file("extdata", "sarscov2_sample.fasta.gz",
                             package = "ViralEntropR")
fasta_sample <- Biostrings::readAAStringSet(path_sample)

# Inspect header structure to confirm date field position.
sample(names(fasta_sample), 1)
#> [1] "QXT18620.1 |USA|2021-01-21"
# The bundled sample uses NCBI Virus format: date is at end of header.

# Default usage on bundled sample: option = 4 for end-of-header dates.
dates <- extract_fasta_dates(fasta_sample, option = 4)
dates$message
#> [1] "There are date strings that have not been recognized"
head(dates$corrected_dates)
#> [1] "2021-02-12" "2020-04-03" "2021-03-05" "2021-02-09" "2020-07-18"
#> [6] "2020-03-23"
range(dates$corrected_dates, na.rm = TRUE)
#> [1] "2020-03-09" "2021-09-10"

# Custom regex for non-standard headers:
dates_custom <- extract_fasta_dates(
  fasta_sample,
  custom_pattern = "[0-9]{4}-(0?[1-9]|1[0-2])-(0?[1-9]|[12][0-9]|3[01]|00)"
)
# }
```
