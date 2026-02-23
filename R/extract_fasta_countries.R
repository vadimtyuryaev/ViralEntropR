#' @title Extract Countries from FASTA Sequence Names
#' @description Extracts country names from the sequence name strings of an
#'   \code{AAStringSet} object loaded via
#'   \code{\link[Biostrings]{readAAStringSet}}. Handles single-word (e.g.
#'   \code{UK}), hyphenated (e.g. \code{Timor-Leste}), and multi-word (e.g.
#'   \code{United States of America}) country names.
#'
#' @param sequence An \code{AAStringSet} object.
#' @param position Integer (1–4). Location of the country field within the
#'   sequence name string:
#'   \itemize{
#'     \item \strong{1} — before the first \code{|} (e.g.
#'       \code{SouthKorea|...}).
#'     \item \strong{2} — between pipe characters \code{|...|}.
#'     \item \strong{3} — between forward slashes \code{/.../}.
#'     \item \strong{4} — after the last \code{|}.
#'   }
#' @param problematic_characters Logical. If \code{TRUE}, sequence names are
#'   re-encoded to UTF-8, replacing non-representable bytes with their escaped
#'   form. Useful for FASTA files with non-ASCII characters in headers.
#'   Default is \code{FALSE}.
#'
#' @return A named list with three elements:
#' \item{countries}{Character vector of extracted country strings, one per
#'   sequence. \code{NA} where extraction failed.}
#' \item{message}{A single character string summarising extraction success.}
#' \item{missing_id}{Integer vector of indices where extraction failed, or
#'   \code{NA} if all extractions succeeded.}
#'
#' @importFrom stringr str_extract
#'
#' @examples
#' \dontrun{
#' fasta <- Biostrings::readAAStringSet("sequences.fasta")
#' result <- extract_fasta_countries(fasta, position = 2)
#' print(result$countries)
#' print(result$message)
#' }
#'
#' @export
extract_fasta_countries <- function(sequence,
                                       position,
                                       problematic_characters = FALSE) {
  
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required. Install it with:\n",
         "  install.packages('BiocManager')\n",
         "  BiocManager::install('Biostrings')",
         call. = FALSE)
  }
  
  if (!position %in% 1:4) {
    stop("`position` must be an integer between 1 and 4.")
  }

  # Use names() — the public accessor — instead of internal slot access.
  sqnce <- names(sequence)

  if (isTRUE(problematic_characters)) {
    sqnce <- iconv(sqnce, "UTF-8", "UTF-8", sub = "byte")
  }

  # Regex patterns — first occurrence of delimiter sought for positions 2 and 3
  pattern <- switch(as.character(position),
    "1" = "^[^|]*",               # text before the first |
    "2" = "(?<=\\|)[^|]+(?=\\|)", # text between pipes |...|
    "3" = "(?<=\\/)[^\\/]+",      # text between slashes /.../
    "4" = "(?<=\\|)[^|]+$"        # text after the last |
  )

  countries  <- stringr::str_extract(sqnce, pattern)
  is_missing <- is.na(countries)

  if (any(is_missing)) {
    msg        <- "Several countries have not been extracted or are missing"
    missing_id <- which(is_missing)
  } else {
    msg        <- "All countries have been extracted"
    missing_id <- NA
  }

  list(countries  = countries,
       message    = msg,
       missing_id = missing_id)
}
