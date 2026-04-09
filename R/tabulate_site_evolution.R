#' @title Tabulate Site Frequency Evolution
#' @description Generates a frequency table showing the amino acid distribution
#'   at a specific site across multiple time partitions, optionally styled and
#'   saved as HTML.
#'
#' @details
#' Aggregates amino acid counts per partition using \code{\link{get_site_counts}},
#' then optionally converts to relative frequencies, applies kableExtra styling
#' (column width, column highlighting), and/or saves to disk. Row names are
#' decoded amino acid letters via \code{\link{decode_aa_sequence}}.
#'
#' @param partitions A list of data frames, typically produced by
#'   \code{\link{partition_time_windows}}.
#' @param site_index Integer. The column index (site) to analyse.
#' @param labels Character vector. Column labels for each partition. Defaults
#'   to \code{names(partitions)}, or \code{"P1", "P2", ...} if unnamed.
#' @param alphabet_size Integer. Total number of possible amino acid codes.
#'   Default is 25.
#' @param zeros Logical. If \code{TRUE} (default), fills missing counts with
#'   \code{0}. If \code{FALSE}, leaves them as \code{""}.
#' @param letters Logical. If \code{TRUE} (default), uses decoded amino acid
#'   letters as row names. If \code{FALSE}, uses numeric codes.
#' @param relative Logical. If \code{TRUE}, converts counts to proportions.
#'   Default is \code{FALSE}.
#' @param digits Integer. Decimal places for rounding when \code{relative =
#'   TRUE}. Default is \code{2}.
#' @param col_width Character. CSS width string applied to all columns (e.g.
#'   \code{"100px"}). Default is \code{"100px"}.
#' @param highlight_col Integer. Column index (1-based, relative to data
#'   columns) to highlight. \code{NULL} for no highlight. Default is
#'   \code{NULL}.
#' @param background Character. CSS background colour for the highlighted
#'   column. Default is \code{"#f0f8ff"}.
#' @param wrap_length Integer. Character width at which to wrap long column
#'   labels using HTML line breaks. Default is \code{10}.
#' @param save Logical. If \code{TRUE}, saves the HTML table to disk. Default
#'   is \code{FALSE}.
#' @param save_extension Character. File extension for saved file. Default is
#'   \code{".html"}.
#' @param save_path Character. Directory for saved file. Default is
#'   \code{getwd()}.
#' @param return_table Logical. If \code{TRUE} (default), returns the raw data
#'   frame in addition to the styled kable as a list. If \code{FALSE}, returns
#'   only the kable object.
#'
#' @return If \code{return_table = TRUE}, a named list:
#'   \item{table}{The raw count (or proportion) data frame.}
#'   \item{styled}{The \code{kableExtra} HTML kable object.}
#'   If \code{return_table = FALSE}, returns only the styled kable object.
#'
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
tabulate_site_evolution = function(partitions,
                                    site_index,
                                    labels         = NULL,
                                    alphabet_size  = 25L,
                                    zeros          = TRUE,
                                    letters        = TRUE,
                                    relative       = FALSE,
                                    digits         = 2L,
                                    col_width      = "100px",
                                    highlight_col  = NULL,
                                    background     = "#f0f8ff",
                                    wrap_length    = 10L,
                                    save           = FALSE,
                                    save_extension = ".html",
                                    save_path      = getwd(),
                                    return_table   = TRUE) {
  
  n_part = length(partitions)
  
  # --- 1. Column labels ------------------------------------------------------
  if (is.null(labels)) {
    p_names = names(partitions)
    labels  = if (!is.null(p_names)) p_names else paste0("P", seq_len(n_part))
  }
  
  # --- 2. Counts via get_site_counts -----------------------------------------
  counts_matrix = get_site_counts(partitions, site_index, alphabet_size)
  prot_df       = as.data.frame(counts_matrix)
  
  # --- 3. Row names ----------------------------------------------------------
  rownames(prot_df) = if (isTRUE(letters)) {
    as.vector(decode_aa_sequence(matrix(seq_len(alphabet_size), ncol = 1L)))
  } else {
    as.character(seq_len(alphabet_size))
  }
  
  # --- 4. Zero / empty fill --------------------------------------------------
  if (!isTRUE(zeros)) {
    prot_df[prot_df == 0] = ""
  }
  
  # --- 5. Relative frequencies -----------------------------------------------
  if (isTRUE(relative)) {
    col_sums = colSums(counts_matrix)
    col_sums[col_sums == 0] = 1               # guard against empty partitions
    prot_df  = as.data.frame(
      round(sweep(counts_matrix, 2L, col_sums, "/"), digits)
    )
    rownames(prot_df) = if (isTRUE(letters)) {
      as.vector(decode_aa_sequence(matrix(seq_len(alphabet_size), ncol = 1L)))
    } else {
      as.character(seq_len(alphabet_size))
    }
    if (!isTRUE(zeros)) prot_df[prot_df == 0] = ""
  }
  
  # --- 6. Wrap column labels for HTML ----------------------------------------
  wrap_text = function(text, width = wrap_length) {
    paste(strwrap(text, width = width), collapse = "<br>")
  }
  formatted_labels = vapply(labels, wrap_text, character(1L))
  
  # --- 7. Build kable --------------------------------------------------------
  styled = prot_df %>%
    kableExtra::kable("html",
                      row.names = TRUE,
                      col.names = formatted_labels,
                      escape    = FALSE) %>%
    kableExtra::kable_styling("striped", full_width = FALSE) %>%
    kableExtra::column_spec(seq_len(n_part), width = col_width)
  
  # --- 8. Optional column highlight ------------------------------------------
  if (!is.null(highlight_col) &&
      highlight_col > 0L &&
      highlight_col <= n_part) {
    styled = styled %>%
      kableExtra::column_spec(highlight_col + 1L,   # +1: kable counts row-name col first
                              background = background)
  }
  
  # --- 9. Optional save ------------------------------------------------------
  if (isTRUE(save)) {
    file_path = file.path(save_path,
                           paste0("Site_", site_index, save_extension))
    styled %>% kableExtra::save_kable(file_path)
  }
  
  # --- 10. Return ------------------------------------------------------------
  if (isTRUE(return_table)) {
    list(table  = prot_df,
         styled = styled)
  } else {
    styled
  }
}