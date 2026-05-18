#' @title Extract Countries from FASTA Sequence Names
#' @description Extracts country names from the sequence name strings of an
#'   \code{AAStringSet} object loaded via
#'   \code{\link[Biostrings]{readAAStringSet}}. Handles single-word (e.g.
#'   \code{UK}), hyphenated (e.g. \code{Timor-Leste}), and multi-word (e.g.
#'   \code{United States of America}) country names.
#'
#' @details
#' The function selects one of four regex patterns based on \code{position}
#' and applies it to each sequence name via
#' \code{\link[stringr]{str_extract}}. \strong{Only the first match per
#' header is returned.} If a header contains multiple delimited fields,
#' the country must be in the first such field for the corresponding
#' \code{position} value to extract it correctly. For example, with a
#' GISAID-style header
#' \code{Spike|hCoV-19/USA/OH/.../2021|2021-05-15|EPI_ISL_...|},
#' \code{position = 3} (between slashes) returns \code{USA}, but
#' \code{position = 2} (between pipes) returns \code{hCoV-19/USA/OH/...},
#' not \code{USA}. Inspect representative headers with
#' \code{names(sequence)[1]} before choosing \code{position}.
#'
#' \strong{Encoding.} FASTA files with non-ASCII characters in headers
#' (accented characters, byte-order marks, etc.) can break regex
#' extraction. Setting \code{problematic_characters = TRUE} re-encodes
#' headers to UTF-8 with non-representable bytes escaped, allowing the
#' regex to proceed.
#'
#' @param sequence An \code{AAStringSet} object.
#' @param position Integer (1–4). Location of the country field within the
#'   sequence name string:
#'   \itemize{
#'     \item \strong{1} — text before the first \code{|} (e.g.
#'       \code{SouthKorea|...}).
#'     \item \strong{2} — text between the first and second \code{|}.
#'     \item \strong{3} — text between the first and second \code{/}.
#'     \item \strong{4} — text after the last \code{|}.
#'   }
#' @param problematic_characters Logical. If \code{TRUE}, sequence names are
#'   re-encoded to UTF-8, replacing non-representable bytes with their
#'   escaped form. Useful for FASTA files with non-ASCII characters in
#'   headers. Default is \code{FALSE}.
#'
#' @return A named list with three elements:
#' \item{countries}{Character vector of extracted country strings, one per
#'   sequence. \code{NA} where extraction failed (no match against the
#'   chosen pattern).}
#' \item{message}{A single character string summarising extraction
#'   success.}
#' \item{missing_id}{Integer vector of indices where extraction failed, or
#'   \code{NA} if all extractions succeeded.}
#'
#' @seealso \code{\link{extract_fasta_dates}} for the date-extraction
#'   companion; \code{\link[Biostrings]{readAAStringSet}} for loading the
#'   input \code{AAStringSet}.
#'
#' @importFrom stringr str_extract
#' @export
#'
#' @examples
#' \donttest{
#' path_sample  <- system.file("extdata", "sarscov2_sample.fasta.gz",
#'                              package = "ViralEntropR")
#' fasta_sample <- Biostrings::readAAStringSet(path_sample)
#'
#' # Inspect header structure to confirm field positions before extraction.
#' sample(names(fasta_sample), 1)
#'
#' # Extract countries (position 2 = between first and second pipe).
#' result <- extract_fasta_countries(fasta_sample, position = 2)
#' result$message
#' sort(table(result$countries), decreasing = TRUE)
#' }
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
    stop("`position` must be an integer between 1 and 4.", call. = FALSE)
  }
  
  # Use names() — the public accessor — instead of internal slot access.
  sqnce <- names(sequence)
  
  if (isTRUE(problematic_characters)) {
    sqnce <- iconv(sqnce, "UTF-8", "UTF-8", sub = "byte")
  }
  
  # Regex patterns. Each pattern matches the FIRST occurrence of the field;
  # see Details for header-format implications.
  pattern <- switch(as.character(position),
                    "1" = "^[^|]+",                # text before the first | (non-empty)
                    "2" = "(?<=\\|)[^|]+(?=\\|)",  # text between pipes |...|
                    "3" = "(?<=\\/)[^\\/]+",       # text between slashes /.../
                    "4" = "(?<=\\|)[^|]+$"         # text after the last |
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