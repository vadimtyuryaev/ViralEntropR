#' @title Tabulate Site Frequency Evolution
#' @description Generates a frequency or proportion table showing the amino
#'   acid distribution at a specific site across multiple time partitions,
#'   optionally styled with \pkg{kableExtra} and saved as standalone HTML.
#'
#' @details
#' Aggregates amino acid counts per partition using
#' \code{\link{get_site_counts}}, optionally converts to relative
#' frequencies, applies kableExtra styling (column width, column
#' highlighting, striped rows), and optionally saves to disk. Row names
#' are decoded amino acid codes via \code{\link{decode_aa_sequence}} when
#' \code{use_letters = TRUE}.
#'
#' \strong{Empty partitions.} Partitions containing no observations
#' contribute all-zero columns. When \code{relative = TRUE}, division by
#' zero is avoided by treating empty-column sums as 1, leaving
#' proportions at zero. Inspect partition sizes via
#' \code{sapply(partitions, nrow)} before interpreting the table.
#'
#' \strong{Note on \code{zeros = FALSE}.} Setting \code{zeros = FALSE}
#' replaces numeric zeros with empty strings (\code{""}) for visual
#' clarity in the kable. This conversion forces the underlying data
#' frame to character storage; numeric operations (\code{sum},
#' \code{mean}, etc.) will not work on the returned \code{table}
#' element. Use \code{zeros = TRUE} (default) if downstream numerical
#' use is intended.
#'
#' @param partitions A list of data frames, typically produced by
#'   \code{\link{partition_time_windows}}. Each data frame must contain
#'   integer-encoded amino acid sequences as columns (values 1 to
#'   \code{alphabet_size}).
#' @param site_index Integer. The column index (site) to analyse.
#' @param labels Character vector. Column labels for each partition.
#'   Defaults to \code{names(partitions)}, or \code{"P1"}, \code{"P2"},
#'   ... if unnamed.
#' @param alphabet_size Integer. Total number of possible amino acid
#'   codes. Must match the encoding used during integer encoding
#'   (\code{\link{encode_aa_sequence}} produces values 1-25 by default).
#'   Default is \code{25L}.
#' @param zeros Logical. If \code{TRUE} (default), fills missing counts
#'   with \code{0}. If \code{FALSE}, replaces zeros with \code{""}; see
#'   Details for the type-coercion implication.
#' @param use_letters Logical. If \code{TRUE} (default), uses decoded
#'   single-character codes from the 25-symbol alphabet (A, R, N, D,
#'   C, ..., V, B, Z, X, *, -) as row names. If \code{FALSE}, uses
#'   numeric codes 1 to \code{alphabet_size}.
#' @param relative Logical. If \code{TRUE}, converts counts to
#'   proportions (column-wise division by partition size). Default is
#'   \code{FALSE}.
#' @param digits Integer. Decimal places for rounding when
#'   \code{relative = TRUE}. Default is \code{2L}.
#' @param col_width Character. CSS width string applied to all columns
#'   (e.g. \code{"100px"}). Default is \code{"100px"}.
#' @param highlight_col Integer or \code{NULL}. 1-based column index
#'   (relative to the data columns, not counting the row-name column)
#'   of a partition to highlight with \code{background}. Out-of-range
#'   values trigger a warning and no highlight is applied. \code{NULL}
#'   (default) means no highlight.
#' @param background Character. CSS background colour for the
#'   highlighted column. Default is \code{"#f0f8ff"} (light blue).
#' @param wrap_length Integer. Character width at which to wrap long
#'   column labels using HTML line breaks. Default is \code{10L}.
#' @param save Logical. If \code{TRUE}, saves the rendered HTML table
#'   to disk via \code{\link[kableExtra]{save_kable}}. Default is
#'   \code{FALSE}.
#' @param save_extension Character. File extension for the saved file
#'   (including leading dot). Default is \code{".html"}.
#' @param save_path Character or \code{NULL}.  Directory in which to save the
#'   file.  ...  Must be supplied when \code{save = TRUE}.  Default is
#'   \code{NULL}.
#' @param return_table Logical. If \code{TRUE} (default), returns a
#'   named list with both the raw data frame and the styled kable
#'   object. If \code{FALSE}, returns only the styled kable.
#'
#' @return If \code{return_table = TRUE}, a named list:
#'   \item{table}{The raw count (or proportion) data frame, with row
#'     names corresponding to amino acid codes and column names
#'     corresponding to partition labels.}
#'   \item{styled}{The \pkg{kableExtra} HTML kable object.}
#'   If \code{return_table = FALSE}, returns only the styled kable
#'   object.
#'
#' @seealso \code{\link{partition_time_windows}} for producing the
#'   typical input list of partitions; \code{\link{get_site_counts}}
#'   for the count-tabulation primitive; \code{\link{decode_aa_sequence}}
#'   for the alphabet code mapping; \code{\link{calculate_hellinger_matrix}}
#'   for the related cross-partition distance calculation on the same
#'   data shape.
#'
#' @importFrom kableExtra kable kable_styling column_spec save_kable
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' p1 = data.frame(s1 = c(1L, 1L, 1L, 1L, 2L))
#' p2 = data.frame(s1 = c(1L, 1L, 2L, 2L, 2L))
#' parts = list(T1 = p1, T2 = p2)
#'
#' # Default: counts, letters, no save
#' tbl = tabulate_site_evolution(parts, site_index = 1)
#' tbl$styled
#'
#' # Relative frequencies, highlight second partition
#' tbl2 = tabulate_site_evolution(parts, site_index = 1,
#'                                 relative = TRUE, highlight_col = 2)
#' tbl2$styled
#'
#' # Numeric codes (skip alphabet decoding)
#' tbl3 = tabulate_site_evolution(parts, site_index = 1, use_letters = FALSE)
#' rownames(tbl3$table)
tabulate_site_evolution = function(partitions,
                                   site_index,
                                   labels         = NULL,
                                   alphabet_size  = 25L,
                                   zeros          = TRUE,
                                   use_letters    = TRUE,
                                   relative       = FALSE,
                                   digits         = 2L,
                                   col_width      = "100px",
                                   highlight_col  = NULL,
                                   background     = "#f0f8ff",
                                   wrap_length    = 10L,
                                   save           = FALSE,
                                   save_extension = ".html",
                                   save_path      = NULL,
                                   return_table   = TRUE) {
  
  # --- 1. Input validation ---------------------------------------------------
  if (!is.list(partitions) || length(partitions) == 0L)
    stop("`partitions` must be a non-empty list of data frames.",
         call. = FALSE)
  
  if (!is.numeric(site_index) || length(site_index) != 1L || site_index < 1L)
    stop("`site_index` must be a single positive integer.",
         call. = FALSE)
  
  if (!is.numeric(alphabet_size) || alphabet_size < 1L)
    stop("`alphabet_size` must be a positive integer.",
         call. = FALSE)
  
  n_part = length(partitions)
  
  # --- 2. Column labels ------------------------------------------------------
  if (is.null(labels)) {
    p_names = names(partitions)
    labels  = if (!is.null(p_names)) p_names else paste0("P", seq_len(n_part))
  }
  
  # --- 3. Pre-compute row names ----------------------------------------------
  # Built once and reused; avoids the duplicate rowname assignment that
  # the original code performed inside the relative-frequency branch.
  row_labels = if (isTRUE(use_letters)) {
    as.vector(decode_aa_sequence(matrix(seq_len(alphabet_size), ncol = 1L)))
  } else {
    as.character(seq_len(alphabet_size))
  }
  
  # --- 4. Counts via get_site_counts -----------------------------------------
  counts_matrix = get_site_counts(partitions, site_index, alphabet_size)
  
  # --- 5. Optional conversion to relative frequencies ------------------------
  if (isTRUE(relative)) {
    col_sums = colSums(counts_matrix)
    col_sums[col_sums == 0] = 1L              # guard against empty partitions
    prot_df = as.data.frame(
      round(sweep(counts_matrix, 2L, col_sums, "/"), digits)
    )
  } else {
    prot_df = as.data.frame(counts_matrix)
  }
  
  rownames(prot_df) = row_labels
  
  # --- 6. Optional zero -> empty-string conversion ---------------------------
  # Note: this forces the data frame to character storage. See @details.
  if (!isTRUE(zeros)) {
    prot_df[prot_df == 0] = ""
  }
  
  # --- 7. Wrap column labels for HTML ----------------------------------------
  wrap_text = function(text, width = wrap_length) {
    paste(strwrap(text, width = width), collapse = "<br>")
  }
  formatted_labels = vapply(labels, wrap_text, character(1L))
  
  # --- 8. Build kable --------------------------------------------------------
  styled = prot_df %>%
    kableExtra::kable("html",
                      row.names = TRUE,
                      col.names = formatted_labels,
                      escape    = FALSE) %>%
    kableExtra::kable_styling("striped", full_width = FALSE) %>%
    kableExtra::column_spec(seq_len(n_part), width = col_width)
  
  # --- 9. Optional column highlight ------------------------------------------
  if (!is.null(highlight_col)) {
    if (highlight_col < 1L || highlight_col > n_part) {
      warning(sprintf(
        "`highlight_col` (%d) is out of range [1, %d]; no highlight applied.",
        highlight_col, n_part), call. = FALSE)
    } else {
      styled = styled %>%
        kableExtra::column_spec(highlight_col + 1L,   # +1: row-name column comes first
                                background = background)
    }
  }
  
  # --- 10. Optional save -----------------------------------------------------
  if (isTRUE(save)) {
    if (is.null(save_path))
      stop("`save_path` must be supplied when `save = TRUE`.", call. = FALSE)
    file_path = file.path(save_path,
                          paste0("Site_", site_index, save_extension))
    styled %>% kableExtra::save_kable(file_path)
  }
  
  # --- 11. Return ------------------------------------------------------------
  if (isTRUE(return_table)) {
    list(table  = prot_df,
         styled = styled)
  } else {
    styled
  }
}