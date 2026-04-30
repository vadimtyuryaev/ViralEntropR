# =============================================================================
# Comparison Table Functions
# =============================================================================
#
# Table design — wide format
# --------------------------
# Rows  = metrics (N, N Overlap, % Overlap, N Missed, Selected Sites,
#                  Missed Sites)
# Columns = GMM subsets (Important reference + one column per candidate)
#
# Highlighting design
# -------------------
# Sites are processed from integer vectors in a single pass — no string
# splitting after HTML injection.  Each integer site is checked individually
# with %in%; matching sites receive a compact yellow background + bold via a
# raw HTML span (no padding, no border-radius).  This avoids the pill-badge
# spacing artefact produced by kableExtra::cell_spec().
# =============================================================================


# ---------------------------------------------------------------------------
# format_site_vector  (internal)
# ---------------------------------------------------------------------------
# Formats an integer site vector as a single HTML string.
# Sites in ref_sites receive yellow background + bold (compact inline span).
# Lines are wrapped every items_per_line tokens.
format_site_vector <- function(site_vec, ref_sites, items_per_line = 10L) {
  if (length(site_vec) == 0L) return("")
  site_vec <- as.integer(site_vec)

  formatted <- vapply(site_vec, function(s) {
    if (!is.na(s) && s %in% ref_sites)
      sprintf('<span style="background-color:#FFFF00;font-weight:bold;">%s</span>',
              as.character(s))
    else
      as.character(s)
  }, character(1L))

  chunks <- split(formatted, ceiling(seq_along(formatted) / items_per_line))
  paste(vapply(chunks, paste, character(1L), collapse = ", "),
        collapse = "<br>")
}


# ---------------------------------------------------------------------------
# extract_sites_from_result  (internal)
# ---------------------------------------------------------------------------
# Extracts the named list of site vectors from a gmm_site_selection() result.
# Handles both original-scale (GMM_orig_sites) and transformed-scale
# (GMM_transfr_sites) results automatically.
extract_sites_from_result <- function(result) {
  os <- result$Sites$GMM_orig_sites
  ts <- result$Sites$GMM_transfr_sites

  if (!is.null(os)) {
    list(
      complete = os$GMM_orig_sites,
      q1plus   = os$GMM_75_sites,
      median   = os$GMM_50_sites,
      q3plus   = os$GMM_25_sites,
      mean     = os$GMM_mean_sites,
      quant    = os$GMM_quant_sites
    )
  } else if (!is.null(ts)) {
    list(
      complete = ts$GMM_transfr_orig_sites,
      q1plus   = ts$GMM_transfr_75_sites,
      median   = ts$GMM_transfr_50_sites,
      q3plus   = ts$GMM_transfr_25_sites,
      mean     = ts$GMM_transfr_mean_sites,
      quant    = ts$GMM_transfr_quant_sites
    )
  } else {
    stop("result contains neither GMM_orig_sites nor GMM_transfr_sites.")
  }
}


# ---------------------------------------------------------------------------
# generate_comparison_table
# ---------------------------------------------------------------------------
#' @title Generate GMM Site Selection Comparison Table (Wide Format)
#'
#' @description
#' Produces a styled HTML comparison table in wide format:
#' \itemize{
#'   \item \strong{Rows} — N, N Overlap, \% Overlap, N Missed,
#'     Selected Sites, Missed Sites.
#'   \item \strong{Columns} — reference (Important sites) plus one column
#'     per candidate GMM subset.
#' }
#' Sites present in the reference set are highlighted yellow and bold in the
#' Selected Sites row.
#'
#' @param vectors_list List of integer vectors.  First = reference sites;
#'   remainder = candidate site sets in column order.
#' @param table_caption Character.  Caption rendered bold and centred.
#' @param col_names Character vector of length \code{length(vectors_list)}.
#'   Column header labels.
#' @param items_per_line Integer.  Site indices per display line.  Default 10.
#'
#' @return A \code{kableExtra} HTML table object.

generate_comparison_table <- function(vectors_list,
                                      table_caption  = NULL,
                                      col_names      = NULL,
                                      items_per_line = 10L) {

  stopifnot(is.list(vectors_list), length(vectors_list) >= 2L)

  reference_sites <- as.integer(vectors_list[[1L]])
  candidate_sets  <- vectors_list[-1L]
  n_candidates    <- length(candidate_sets)
  n_important     <- length(reference_sites)
  n_cols          <- 1L + n_candidates

  if (is.null(col_names))
    col_names <- c("Important", paste0("Vector", seq_len(n_candidates)))

  # ── Per-column statistics ─────────────────────────────────────────────────
  n_selected <- c(n_important,
                  vapply(candidate_sets, function(v) length(as.integer(v)),
                         integer(1L)))

  n_overlap  <- c(NA_integer_,
                  vapply(candidate_sets, function(v)
                    length(intersect(reference_sites, as.integer(v))),
                    integer(1L)))

  pct_overlap <- c("\u2014",
                   vapply(n_overlap[-1L], function(n)
                     sprintf("%.1f%%", n / n_important * 100),
                     character(1L)))

  n_missed   <- c(NA_integer_,
                  vapply(candidate_sets, function(v)
                    length(setdiff(reference_sites, as.integer(v))),
                    integer(1L)))

  sites_row <- c(
    format_site_vector(reference_sites, integer(0L), items_per_line),
    vapply(candidate_sets, function(v)
      format_site_vector(as.integer(v), reference_sites, items_per_line),
      character(1L))
  )

  missed_row <- c(
    "\u2014",
    vapply(candidate_sets, function(v) {
      missed <- setdiff(reference_sites, as.integer(v))
      if (length(missed) == 0L) "None"
      else format_site_vector(missed, integer(0L), items_per_line)
    }, character(1L))
  )

  # ── Assemble wide data frame ──────────────────────────────────────────────
  metric_labels <- c("N", "N Overlap", "% Overlap",
                     "N Missed", "Selected Sites", "Missed Sites")

  values_mat <- rbind(
    ifelse(is.na(n_selected), "\u2014", as.character(n_selected)),
    ifelse(is.na(n_overlap),  "\u2014", as.character(n_overlap)),
    pct_overlap,
    ifelse(is.na(n_missed),   "\u2014", as.character(n_missed)),
    sites_row,
    missed_row
  )

  df <- as.data.frame(values_mat, stringsAsFactors = FALSE)
  colnames(df) <- col_names
  df <- cbind(Metric = metric_labels, df)

  caption_html <- if (!is.null(table_caption))
    paste0("<center><strong>", table_caption, "</strong><br>",
           "<small>Reference sites highlighted in ",
           "<span style='background:#FFFF00;font-weight:bold;padding:0 3px;'>",
           "yellow</span>. ",
           "N\u00a0Overlap and %\u00a0Overlap relative to N\u00a0Important.",
           "</small></center>")
  else
    paste0("<center><small>Reference sites highlighted in ",
           "<span style='background:#FFFF00;font-weight:bold;padding:0 3px;'>",
           "yellow</span>.</small></center>")

  kableExtra::kable(
    df,
    format    = "html",
    escape    = FALSE,
    caption   = caption_html,
    row.names = FALSE
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width        = TRUE,
      font_size         = 12L
    ) %>%
    kableExtra::column_spec(1L, bold = TRUE, width = "9em",
                            extra_css = "white-space:nowrap;") %>%
    kableExtra::column_spec(2L, width = "14em") %>%
    kableExtra::column_spec(seq(3L, n_cols + 1L),
                            width = sprintf("%.1fem",
                                            (100 - 9 - 14) / n_candidates)) %>%
    kableExtra::row_spec(5L, extra_css = "vertical-align:top;") %>%
    kableExtra::row_spec(6L, extra_css = "vertical-align:top;")
}


# ---------------------------------------------------------------------------
# generate_summary_table
# ---------------------------------------------------------------------------
#' @title Generate Cross-Transformation Summary Table
#'
#' @description
#' Aggregates \% Overlap with the reference VOC sites across all
#' transformations and GMM subsets into a single compact summary table.
#'
#' Rows correspond to GMM subsets (Complete, Top 75\%, Top 50\%, Top 25\%,
#' Mean+, 90th percentile+); columns correspond to transformations.  Cell
#' values are \% Overlap with the reference site set.
#'
#' @param results_list Named list of \code{gmm_site_selection()} result
#'   objects.  Names become column headers (transformation labels).  Each
#'   element may contain either \code{GMM_orig_sites} (untransformed) or
#'   \code{GMM_transfr_sites} (transformed); the correct slot is detected
#'   automatically.
#' @param important_sites Integer vector.  Reference VOC-defining sites.
#' @param table_caption Character.  Table caption.  Default \code{NULL}.
#' @param save_path Character.  Full file path (\code{.html}) for saving.
#'   If \code{NULL} the table is returned but not saved.  Default \code{NULL}.
#'
#' @return Invisibly returns the \code{kableExtra} HTML table object.

generate_summary_table <- function(results_list,
                                   important_sites,
                                   table_caption = NULL,
                                   save_path     = NULL) {

  stopifnot(is.list(results_list), length(results_list) >= 1L,
            !is.null(base::names(results_list)))

  important_sites <- as.integer(important_sites)
  n_important     <- length(important_sites)

  subset_labels <- c("Complete", "Top 75%", "Top 50%",
                     "Top 25%", "Mean+", "90th pct+")
  subset_keys   <- c("complete", "q1plus", "median",
                     "q3plus", "mean", "quant")

  # ── Compute % Overlap for every (transformation × subset) cell ───────────
  pct_mat <- vapply(results_list, function(res) {
    sites <- extract_sites_from_result(res)
    vapply(subset_keys, function(key) {
      sv <- as.integer(sites[[key]])
      if (length(sv) == 0L) return(NA_real_)
      length(intersect(important_sites, sv)) / n_important * 100
    }, numeric(1L))
  }, numeric(length(subset_keys)))

  # pct_mat: rows = subsets, columns = transformations
  df <- as.data.frame(pct_mat, stringsAsFactors = FALSE)
  colnames(df) <- base::names(results_list)

  # Format cells: "XX.X%" with NA shown as "—"
  df_fmt <- as.data.frame(
    lapply(df, function(col)
      ifelse(is.na(col), "\u2014", sprintf("%.1f%%", col))),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  df_fmt <- cbind(Subset = subset_labels, df_fmt)

  # ── Colour scale: higher % = darker green background ─────────────────────
  # Applied to the numeric matrix before formatting for correct ordering
  pct_vec    <- as.vector(pct_mat)
  pct_finite <- pct_vec[is.finite(pct_vec)]
  pct_min    <- if (length(pct_finite)) min(pct_finite) else 0
  pct_max    <- if (length(pct_finite)) max(pct_finite) else 100

  # Map each numeric value to a green shade (#FFFFFF at min → #1a7a4a at max)
  green_shade <- function(pct) {
    if (!is.finite(pct)) return(NULL)
    t   <- (pct - pct_min) / max(pct_max - pct_min, 1)
    r   <- as.integer(255 - t * (255 - 26))
    g   <- as.integer(255 - t * (255 - 122))
    b   <- as.integer(255 - t * (255 - 74))
    sprintf("background-color:rgb(%d,%d,%d);", r, g, b)
  }

  caption_html <- if (!is.null(table_caption))
    paste0("<center><strong>", table_caption, "</strong><br>",
           "<small>Cell values: \u0025 overlap with VOC reference sites. ",
           "Colour intensity proportional to overlap percentage.",
           "</small></center>")
  else
    paste0("<center><small>Cell values: \u0025 overlap with VOC reference ",
           "sites. Colour intensity proportional to overlap percentage.",
           "</small></center>")

  n_transf <- length(results_list)

  tbl <- kableExtra::kable(
    df_fmt,
    format    = "html",
    escape    = FALSE,
    caption   = caption_html,
    row.names = FALSE,
    align     = c("l", rep("c", n_transf))
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width        = TRUE,
      font_size         = 12L
    ) %>%
    kableExtra::column_spec(1L, bold = TRUE, width = "9em",
                            extra_css = "white-space:nowrap;")

  # Apply per-cell colour using column_spec with cell-level background
  for (col_i in seq_len(n_transf)) {
    col_vals   <- pct_mat[, col_i]
    col_css    <- vapply(col_vals, function(p) {
      shade <- green_shade(p)
      if (is.null(shade)) "" else shade
    }, character(1L))

    # Bold the maximum value in each column
    is_max <- is.finite(col_vals) & col_vals == max(col_vals, na.rm = TRUE)
    col_css[is_max] <- paste0(col_css[is_max], "font-weight:bold;")

    tbl <- kableExtra::column_spec(
      tbl,
      column    = col_i + 1L,   # +1 for Subset label column
      extra_css = col_css
    )
  }

  if (!is.null(save_path))
    kableExtra::save_kable(tbl, file = save_path)

  invisible(tbl)
}


# ---------------------------------------------------------------------------
# save_orig_comparison_table
# ---------------------------------------------------------------------------
save_orig_comparison_table <- function(result_list,
                                       table_caption,
                                       save_path,
                                       important_sites) {

  os <- result_list$Sites$GMM_orig_sites
  if (is.null(os))
    stop("calc_gmm_orig must be TRUE to call save_orig_comparison_table().")

  col_names_vec <- c("Important", "Complete", "Top 75%", "Top 50%",
                     "Top 25%", "Mean+", "90th pct+")

  tbl <- generate_comparison_table(
    list(as.integer(important_sites),
         as.integer(os$GMM_orig_sites),
         as.integer(os$GMM_75_sites),
         as.integer(os$GMM_50_sites),
         as.integer(os$GMM_25_sites),
         as.integer(os$GMM_mean_sites),
         as.integer(os$GMM_quant_sites)),
    table_caption  = table_caption,
    col_names      = col_names_vec,
    items_per_line = 8L
  )

  kableExtra::save_kable(tbl, file = save_path)
  invisible(tbl)
}


# ---------------------------------------------------------------------------
# save_transformed_comparison_table
# ---------------------------------------------------------------------------
save_transformed_comparison_table <- function(result_list,
                                              table_caption,
                                              save_path,
                                              important_sites) {

  ts <- result_list$Sites$GMM_transfr_sites
  if (is.null(ts))
    stop("calc_gmm_transfr must be TRUE to call save_transformed_comparison_table().")

  col_names_vec <- c("Important", "Complete", "Top 75%", "Top 50%",
                     "Top 25%", "Mean+", "90th pct+")

  tbl <- generate_comparison_table(
    list(as.integer(important_sites),
         as.integer(ts$GMM_transfr_orig_sites),
         as.integer(ts$GMM_transfr_75_sites),
         as.integer(ts$GMM_transfr_50_sites),
         as.integer(ts$GMM_transfr_25_sites),
         as.integer(ts$GMM_transfr_mean_sites),
         as.integer(ts$GMM_transfr_quant_sites)),
    table_caption  = table_caption,
    col_names      = col_names_vec,
    items_per_line = 8L
  )

  kableExtra::save_kable(tbl, file = save_path)
  invisible(tbl)
}
