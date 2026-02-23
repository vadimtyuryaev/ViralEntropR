#' @title Extract Dates from FASTA Sequence Names
#' @description Extracts date strings from the sequence name strings of an
#'   \code{AAStringSet} object loaded via
#'   \code{\link[Biostrings]{readAAStringSet}}. Several built-in date patterns
#'   are provided; a fully custom regex can also be supplied.
#'
#' @details
#' Date strings of the form \code{yyyy-mm-dd} are matched between pipe
#' characters (\code{|...|}) by default. Day value \code{00} is accepted in
#' the raw string and corrected to \code{01} before coercion to \code{Date}.
#'
#' When \code{date_format = "\%Y-\%m"} the function uses
#' \code{\link[zoo]{as.yearmon}} for coercion so that year-month strings are
#' handled correctly.
#'
#' @param sequence An \code{AAStringSet} object.
#' @param option Integer. Selects the built-in pattern when \code{custom_pattern}
#'   is not supplied:
#'   \itemize{
#'     \item \strong{1} (default) — \code{yyyy-mm-dd} between pipes.
#'     \item \strong{2} — \code{yyyy-dd-mm} between pipes.
#'     \item \strong{3} — \code{yyyy-mm} between pipes.
#'   }
#' @param date_format Character. \code{strptime}-style format string used to
#'   coerce extracted strings to \code{Date}. Default is \code{"\%Y-\%m-\%d"}.
#' @param custom_pattern Character or \code{NULL}. A custom regex passed
#'   directly to \code{\link[stringr]{str_extract}}. When supplied,
#'   \code{option} is ignored. Default is \code{NULL}.
#'
#' @return A named list:
#' \item{raw_date_strings}{Character vector of extracted date strings before
#'   any correction.}
#' \item{corrected_date_strings}{Character vector with \code{-00} replaced by
#'   \code{-01}.}
#' \item{raw_dates}{\code{Date} vector coerced from \code{raw_date_strings}.}
#' \item{corrected_dates}{\code{Date} vector coerced from
#'   \code{corrected_date_strings}.}
#' \item{message}{Character string summarising extraction success.}
#' \item{missing_id}{Integer vector of indices where extraction failed, or
#'   \code{NA} if all extractions succeeded.}
#'
#' @importFrom stringr str_extract
#' @importFrom zoo as.yearmon
#'
#' @examples
#' \dontrun{
#' fasta  <- Biostrings::readAAStringSet("sequences.fasta")
#' result <- extract_fasta_dates(fasta, option = 1)
#' print(head(result$corrected_dates))
#' print(result$message)
#'
#' # Custom pattern
#' result2 <- extract_fasta_dates(fasta,
#'                                   custom_pattern = "[0-9]{4}-[0-9]{2}-[0-9]{2}")
#' }
#'
#' @export
extract_fasta_dates <- function(sequence,
                                   option         = 1,
                                   date_format    = "%Y-%m-%d",
                                   custom_pattern = NULL) {
  
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required. Install it with:\n",
         "  install.packages('BiocManager')\n",
         "  BiocManager::install('Biostrings')",
         call. = FALSE)
  }

  # Use names() — public accessor — instead of internal slot access
  seq_names <- names(sequence)

  # Select regex pattern
  if (!is.null(custom_pattern)) {
    pattern <- custom_pattern
  } else {
    pattern <- switch(as.character(option),
      "1" = "(?<=\\|)[0-9]{4}-(0?[1-9]|1[0-2])-(0?[1-9]|[12][0-9]|3[01]|00)(?=\\|)",
      "2" = "(?<=\\|)[0-9]{4}-(0?[1-9]|[12][0-9]|3[01]|00)-(0?[1-9]|1[0-2])(?=\\|)",
      "3" = "(?<=\\|)[0-9]{4}-(0?[1-9]|1[0-2])(?=\\|)",
      stop("`option` must be 1, 2, or 3.")
    )
  }

  date_strings <- stringr::str_extract(seq_names, pattern)
  is_missing   <- is.na(date_strings)

  if (!any(is_missing)) {
    msg      <- "All date strings have been successfully extracted"
    miss_id  <- NA
  } else {
    msg      <- "There are date strings that have not been recognized"
    miss_id  <- which(is_missing)
  }

  raw_date_strings       <- date_strings
  corrected_date_strings <- gsub("-00", "-01", raw_date_strings)

  # Coerce to Date — yearmon path for yyyy-mm format
  case_yyyy_mm <- identical(date_format, "%Y-%m")

  if (case_yyyy_mm) {
    raw_dates       <- as.Date(zoo::as.yearmon(raw_date_strings,       format = "%Y-%m"))
    corrected_dates <- as.Date(zoo::as.yearmon(corrected_date_strings, format = "%Y-%m"))
  } else {
    raw_dates       <- as.Date(raw_date_strings,       format = date_format)
    corrected_dates <- as.Date(corrected_date_strings, format = date_format)
  }

  list(raw_date_strings       = raw_date_strings,
       corrected_date_strings = corrected_date_strings,
       raw_dates              = raw_dates,
       corrected_dates        = corrected_dates,
       message                = msg,
       missing_id             = miss_id)
}
