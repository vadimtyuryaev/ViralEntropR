# Parse SARS-CoV-2 VOC/VOI Variant Metadata from Excel

Reads a structured Excel workbook of SARS-CoV-2 Variants of Concern
(VOC) and Variants of Interest (VOI) and returns a named list containing
mutation profiles, nomenclature, temporal detection metadata, defining
SNPs, and a fully citable reference table with interactive display
support.

## Usage

``` r
get_variants(tibble, check = FALSE)
```

## Arguments

- tibble:

  A `tibble` or `data.frame` produced by
  [`read_excel`](https://readxl.tidyverse.org/reference/read_excel.html)
  from `SARS_CoV_2_VOC_VOI.xlsx`. Column 1 must be the `Variant` field
  column; columns 2 onward are one column per variant in the order
  listed above.

- check:

  Logical. If `TRUE`, validates that all internal vectors have equal
  length and stops with an informative message on mismatch. Default
  `FALSE`.

## Value

A named list with 13 elements; see Details for full descriptions.

## Details

Intended to be called once during data preparation via
`data-raw/sarscov2_variants.R`:


      variants_dat  <- readxl::read_excel("SARS_CoV_2_VOC_VOI.xlsx")
      variants_list <- get_variants(variants_dat)
      saveRDS(variants_list,
              file = "inst/extdata/sarscov2_variants.rds")

The saved object can then be loaded anywhere in the package or
vignettes:


      voc_data <- readRDS(system.file("extdata", "sarscov2_variants.rds",
                                      package = "ViralEntropR"))

**Column order of variants in the Excel workbook:** Alpha, Beta,
Epsilon, Eta, Iota, Kappa, Delta, Lambda, Gamma, Zeta, Theta, Omicron.

**Returned list elements:**

- `WHO_Label`:

  List. WHO variant label strings (e.g. `"Alpha"`).

- `Pango_Lineage`:

  List. Pango lineage designations.

- `GISAID_Clade`:

  Character vector. GISAID clade strings.

- `Nextstrain_Clade`:

  Character vector. Nextstrain clade strings.

- `Country_First_Detected`:

  List. Country of first documented detection per variant.

- `Date_Earliest_Sample`:

  Character vector. Month-Year of the earliest documented sample per
  variant.

- `Date_First_Detected`:

  Character vector. Month-Year of world-level first detection.

- `Date_First_Detected_US`:

  Character vector. Month-Year of first US detection.

- `Spike_Mutations`:

  List. Spike protein mutation strings read from the Excel workbook.

- `Mutation_Sites`:

  List. Integer site positions extracted from `Spike_Mutations`.

- `Defining_SNPs`:

  List. Canonical defining SNP strings per variant (`NA` where not
  characterised).

- `Defining_SNP_Sites`:

  List. Integer positions extracted from `Defining_SNPs`.

- `References`:

  Named list with three elements: `$data` — data frame of 21 verified
  references; `$display(variant = NULL)` — interactive
  [`datatable`](https://rdrr.io/pkg/DT/man/datatable.html), optionally
  filtered by WHO label; `$cite(variant)` — character vector of
  formatted citation strings suitable for manuscript use.

## Examples

``` r
if (FALSE) { # \dontrun{
variants_dat  <- readxl::read_excel("SARS_CoV_2_VOC_VOI.xlsx")
voc_data      <- get_variants(variants_dat, check = TRUE)

# Access specific fields
voc_data$Pango_Lineage[[which(voc_data$WHO_Label == "Alpha")]]
voc_data$Defining_SNPs[[which(voc_data$WHO_Label == "Gamma")]]

# Reference table
voc_data$References$data                        # full data frame (20 refs)
voc_data$References$display()                   # interactive DT (all refs)
voc_data$References$display(variant = "Alpha")  # filtered to Alpha refs
voc_data$References$cite("Omicron")             # formatted citation strings
} # }
```
