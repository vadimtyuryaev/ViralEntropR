## data-raw/sarscov2_ncbi.R
##
## Extracts a random sample of 100 sequences from the full NCBI FASTA and
## saves them to inst/extdata/sarscov2_sample.fasta.gz unchanged.
##
## HOW TO RUN (from package root):
##   source("data-raw/sarscov2_ncbi.R")
## ---------------------------------------------------------------------------

if (!requireNamespace("Biostrings", quietly = TRUE))
  stop("Package 'Biostrings' is required: ",
       "BiocManager::install('Biostrings')", call. = FALSE)
if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required: install.packages('here')", call. = FALSE)

set.seed(2021)

# -- 0. Paths -----------------------------------------------------------------

pkg_root <- tryCatch(here::here(), error = function(e) getwd())

if (!file.exists(file.path(pkg_root, "DESCRIPTION")))
  stop("Cannot locate package root (no DESCRIPTION at: ", pkg_root, ").",
       call. = FALSE)

message("Package root: ", pkg_root)

fasta_path <- file.path(pkg_root, "data-raw", "sequences.fasta")
out_path   <- file.path(pkg_root, "inst", "extdata", "sarscov2_sample.fasta.gz")
n_sample   <- 100L

if (!file.exists(fasta_path))
  stop("FASTA not found at: ", fasta_path,
       "\nUpdate fasta_path at the top of the script.", call. = FALSE)

# -- 1. Read ------------------------------------------------------------------

message("Reading FASTA...")
fasta <- Biostrings::readAAStringSet(fasta_path)
message(sprintf("  Loaded %d sequences.", length(fasta)))

# -- 2. Sample ----------------------------------------------------------------

idx          <- sample(length(fasta), min(n_sample, length(fasta)))
sample_fasta <- fasta[idx]
message(sprintf("  Sampled %d sequences.", length(sample_fasta)))

# -- 3. Save ------------------------------------------------------------------

if (!dir.exists(dirname(out_path)))
  dir.create(dirname(out_path), recursive = TRUE)

Biostrings::writeXStringSet(sample_fasta, filepath = out_path, compress = TRUE)

message("Saved -> ", out_path)
message(sprintf("File size: %.1f KB", file.size(out_path) / 1024))
