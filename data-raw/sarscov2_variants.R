## data-raw/sarscov2_variants.R
##
## Builds the curated SARS-CoV-2 VOC/VOI metadata object and saves it to
## data/sarscov2_variants.rda for lazy loading inside the package.
##
## HOW TO RUN (from package root):
##   source("data-raw/sarscov2_variants.R")
##
## Prerequisites:
##   - SARS_CoV_2_VOC_VOI.xlsx is in data-raw/
##   - Package readxl is installed
##   - get_variants.R is in R/  (or package is loaded via devtools::load_all())
##
## Output:
##   data/sarscov2_variants.rda
##
## Users access the object simply by loading the package:
##   library(ViralEntropR)
##   sarscov2_variants        # available immediately via lazy loading
##   ?sarscov2_variants       # help page
## ---------------------------------------------------------------------------

# -- 0. Resolve package root --------------------------------------------------

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(0)$ofile)),
  error = function(e) getwd()
)

pkg_root <- if (basename(script_dir) == "data-raw") {
  dirname(script_dir)
} else if (file.exists(file.path(script_dir, "DESCRIPTION"))) {
  script_dir
} else if (file.exists(file.path(getwd(), "DESCRIPTION"))) {
  getwd()
} else {
  stop(
    "Cannot locate the package root (no DESCRIPTION file found).\n",
    "Please source this script from the package root directory:\n",
    "  source('data-raw/sarscov2_variants.R')",
    call. = FALSE
  )
}

message("Package root: ", pkg_root)

# -- 1. Read Excel ------------------------------------------------------------

if (!requireNamespace("readxl", quietly = TRUE))
  stop("Package 'readxl' is required: install.packages('readxl')", call. = FALSE)

excel_path <- file.path(pkg_root, "data-raw", "SARS_CoV_2_VOC_VOI.xlsx")
if (!file.exists(excel_path))
  stop(
    "Excel workbook not found at:\n  ", excel_path,
    "\nPlace SARS_CoV_2_VOC_VOI.xlsx in data-raw/ and re-run.",
    call. = FALSE
  )

variants_dat <- readxl::read_excel(excel_path)
message("Read Excel: ", nrow(variants_dat), " rows, ",
        ncol(variants_dat) - 1L, " variants.")

# -- 2. Load get_variants() ---------------------------------------------------

if (!exists("get_variants", mode = "function")) {
  r_path <- file.path(pkg_root, "R", "get_variants.R")
  if (file.exists(r_path)) {
    source(r_path)
    message("Sourced: ", r_path)
  } else {
    stop(
      "get_variants() not found. Run devtools::load_all() first.", 
      call. = FALSE
    )
  }
}

# -- 3. Build -----------------------------------------------------------------

sarscov2_variants <- get_variants(variants_dat, check = TRUE)

# -- 4. Verify ----------------------------------------------------------------

expected_names <- c(
  "WHO_Label", "Pango_Lineage", "GISAID_Clade", "Nextstrain_Clade",
  "Country_First_Detected",
  "Date_Earliest_Sample", "Date_First_Detected", "Date_First_Detected_US",
  "Spike_Mutations", "Mutation_Sites",
  "Defining_SNPs", "Defining_SNP_Sites",
  "References"
)

stopifnot(
  "Result is not a list"                = is.list(sarscov2_variants),
  "Expected 12 variants"                = length(sarscov2_variants$WHO_Label) == 12L,
  "Missing list elements"               = all(expected_names %in% names(sarscov2_variants)),
  "References$data is not a data frame" = is.data.frame(sarscov2_variants$References$data),
  "Expected 21 references"              = nrow(sarscov2_variants$References$data) == 21L
)

message("Integrity check passed: 12 variants, 21 references.")

# -- 5. Save as .rda ----------------------------------------------------------

out_dir <- file.path(pkg_root, "data")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_path <- file.path(out_dir, "sarscov2_variants.rda")
save(sarscov2_variants, file = out_path, compress = "bzip2")

message("Saved -> ", out_path)
message(sprintf("File size: %.1f KB", file.size(out_path) / 1024))
message("\nObject is now available via lazy loading after library(ViralEntropR):")
message("  sarscov2_variants")
message("  ?sarscov2_variants")
