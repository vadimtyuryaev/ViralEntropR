# =============================================================================
# empirical_pdel.R
# =============================================================================
#
# Empirical estimator of the per-(sequence, position) deleterious-mutation
# rate p_del, computed from a surveillance feature matrix. Used to anchor
# the simulator's noise model in the Sample-Size Simulation Study.
#
# Definition (per-site noise rate at a frequency threshold):
#   For one position s with observed residue counts c_1, ..., c_K
#   summing to total_s = sum_k c_k, a residue is "real" if its observed
#   proportion at s is at or above the threshold, and "noise" otherwise:
#     noise_count_s = sum_{k : c_k / total_s < threshold} c_k
#     rate_s        = noise_count_s / total_s
#   The pooled rate across multiple sites is
#     pooled_rate   = sum_s noise_count_s / sum_s total_s.
#
# Empirically across 21 SARS-CoV-2 variant-defining Spike positions
# (Alpha, Beta, Gamma, Delta, Omicron defining sites) with the GISAID
# cumulative feature matrix through May 2022, the 0.1%-threshold pooled
# rate is 9.56e-4 with a per-site range of 1.70e-4 to 1.74e-3, supporting
# the simulator's default p_del = 1e-3 as a conservative round value.
#
# The function is purely computational. The README chunk in section 6.2
# is the consumer and is responsible for rendering the returned tables.
#
# Author : Vadim Tyuryaev
# =============================================================================

# -----------------------------------------------------------------------------
# empirical_pdel
# -----------------------------------------------------------------------------
# Computes the per-site empirical noise rate at a given proportion threshold
# from a surveillance feature matrix.
#
# Inputs
#   feature_matrix    : data frame or integer matrix. Rows = sequences,
#                       columns = positions (and optionally metadata
#                       columns named in `drop_columns`). Integer-encoded
#                       amino acids in 1..aa_levels following the package's
#                       encode_aa_sequence convention.
#   sites             : integer vector of position-column indices, or NULL.
#                       NULL (default) evaluates every position column
#                       (i.e., the entire data after dropping metadata).
#                       Non-NULL evaluates only the specified positions
#                       and is the typical mode for variant-defining-site
#                       audits (e.g., the 21-site VOC list for the README).
#   threshold         : numeric scalar in [0, 1). Per-site proportion
#                       threshold below which a residue is classified as
#                       noise. Default 0.001 (0.1%).
#   aa_levels         : integer alphabet size, passed to tabulate's nbins
#                       and to decode_aa_sequence. Default 25L, matching
#                       the package's encode_aa_sequence.
#   position_columns  : integer vector of column indices that contain
#                       position data, or NULL. NULL triggers auto-detect:
#                       all columns of an input matrix, or all columns
#                       of a data frame except those named in drop_columns.
#                       Pass explicitly when your matrix has a non-standard
#                       layout (e.g., metadata columns with non-default
#                       names).
#   drop_columns      : character vector of metadata column names to
#                       exclude during auto-detect. Default c("Date",
#                       "Country") matches the NCBI and GISAID Spike
#                       feature matrices bundled with this study. Ignored
#                       when `position_columns` is supplied.
#
# Returns: a list with
#   $per_site : data frame with one row per evaluated site. Columns:
#                 Site            integer position index
#                 Total           integer count of observations at the site
#                 Real_residues   character string of "+"-concatenated
#                                 single-letter codes for residues whose
#                                 proportion at the site is >= threshold
#                 Noise_count     integer sum of counts at residues below
#                                 the threshold
#                 Per_site_rate   numeric noise_count / total
#   $summary  : list with
#                 pooled_rate           pooled noise count / total obs
#                 mean_rate             mean of per-site rates
#                 median_rate           median of per-site rates
#                 min_rate, min_site    minimum per-site rate and its site
#                 max_rate, max_site    maximum per-site rate and its site
#                 n_sites               number of sites evaluated
#                 threshold             input threshold (for audit trail)
#                 total_noise_events    sum of noise counts across sites
#                 total_observations    sum of totals across sites
#   $params   : list with the input threshold, aa_levels, the resolved
#               position_columns vector, and the resolved sites vector.
#
# Notes
#   - The function uses base::tabulate on integer column vectors directly,
#     matching the primitive that ViralEntropR:::get_site_counts wraps for
#     the partitioned case. Residue letters come from
#     ViralEntropR::decode_aa_sequence so the alphabet mapping is
#     guaranteed consistent with the rest of the package.
#   - Sites with total_s = 0 (a fully-NA column, in principle) return
#     Per_site_rate = NA and are dropped from the summary statistics.
#   - The 25-symbol alphabet includes three ambiguous codes (B/Z/X), the
#     stop codon (*), and the gap (-). At biologically meaningful sites
#     these typically appear at very low frequency and are classified as
#     noise; that is a feature, not a bug, since the simulator's noise
#     model intentionally lumps all transient low-frequency variation
#     (de novo deleterious mutations + sequencing artifacts + ambiguous
#     calls) into a single rate.
# -----------------------------------------------------------------------------
empirical_pdel <- function(feature_matrix,
                           sites             = NULL,
                           threshold         = 0.001,
                           aa_levels         = 25L,
                           position_columns  = NULL,
                           drop_columns      = c("Date", "Country")) {

  # --- 1. Input validation -------------------------------------------------
  if (!is.data.frame(feature_matrix) && !is.matrix(feature_matrix))
    stop("`feature_matrix` must be a data frame or an integer matrix.",
         call. = FALSE)

  threshold <- as.numeric(threshold)
  if (length(threshold) != 1L || is.na(threshold) ||
      threshold < 0 || threshold >= 1)
    stop("`threshold` must be a single numeric value in [0, 1).",
         call. = FALSE)

  aa_levels <- as.integer(aa_levels)
  if (length(aa_levels) != 1L || is.na(aa_levels) || aa_levels < 2L)
    stop("`aa_levels` must be a single integer >= 2.", call. = FALSE)

  # --- 2. Determine position columns ---------------------------------------
  n_cols <- ncol(feature_matrix)
  if (is.null(position_columns)) {
    if (is.data.frame(feature_matrix)) {
      all_names <- colnames(feature_matrix)
      drop_idx  <- stats::na.omit(match(drop_columns, all_names))
      position_columns <- setdiff(seq_len(n_cols), as.integer(drop_idx))
    } else {
      position_columns <- seq_len(n_cols)
    }
  } else {
    position_columns <- as.integer(position_columns)
    if (any(is.na(position_columns)) ||
        any(position_columns < 1L) || any(position_columns > n_cols))
      stop("`position_columns` contains invalid indices for this matrix.",
           call. = FALSE)
  }
  if (length(position_columns) == 0L)
    stop("No position columns identified. Pass `position_columns` ",
         "explicitly or adjust `drop_columns`.", call. = FALSE)

  # --- 3. Determine sites --------------------------------------------------
  if (is.null(sites)) {
    sites <- position_columns
  } else {
    sites <- as.integer(sites)
    if (any(is.na(sites)))
      stop("`sites` contains NA values.", call. = FALSE)
    if (!all(sites %in% position_columns))
      stop("`sites` contains indices outside the position-column range. ",
           "If your matrix has a non-standard layout, pass ",
           "`position_columns` explicitly.", call. = FALSE)
  }
  if (length(sites) == 0L)
    stop("`sites` resolved to an empty vector.", call. = FALSE)

  # --- 4. Letter lookup via the package's decode_aa_sequence --------------
  letters_lookup <- as.vector(ViralEntropR::decode_aa_sequence(
    matrix(seq_len(aa_levels), ncol = 1L)
  ))

  # --- 5. Per-site computation --------------------------------------------
  # Same primitive as ViralEntropR:::get_site_counts (base::tabulate on
  # integer column vectors), applied directly to the feature matrix
  # without wrapping in a single-element partitions list.
  is_df <- is.data.frame(feature_matrix)

  per_site_rows <- lapply(sites, function(s) {
    col_raw <- if (is_df) feature_matrix[[s]] else feature_matrix[, s]
    col     <- as.integer(col_raw)
    counts  <- tabulate(col, nbins = aa_levels)
    total   <- sum(counts)

    if (total == 0L) {
      return(data.frame(
        Site          = s,
        Total         = 0L,
        Real_residues = NA_character_,
        Noise_count   = 0L,
        Per_site_rate = NA_real_,
        stringsAsFactors = FALSE,
        check.names      = FALSE
      ))
    }

    cutoff       <- total * threshold
    is_real      <- counts >= cutoff
    real_codes   <- which(is_real)
    real_letters <- letters_lookup[real_codes]
    noise_count  <- sum(counts[!is_real])

    data.frame(
      Site          = s,
      Total         = total,
      Real_residues = paste(real_letters, collapse = "+"),
      Noise_count   = noise_count,
      Per_site_rate = noise_count / total,
      stringsAsFactors = FALSE,
      check.names      = FALSE
    )
  })
  per_site <- do.call(rbind, per_site_rows)
  rownames(per_site) <- NULL

  # --- 6. Summary across evaluated sites ----------------------------------
  rates       <- per_site$Per_site_rate
  total_noise <- sum(per_site$Noise_count, na.rm = TRUE)
  total_obs   <- sum(per_site$Total,       na.rm = TRUE)
  pooled_rate <- if (total_obs > 0L) total_noise / total_obs else NA_real_

  valid       <- !is.na(rates)
  rates_v     <- rates[valid]
  sites_v     <- per_site$Site[valid]

  summary_list <- list(
    pooled_rate         = pooled_rate,
    mean_rate           = if (length(rates_v) > 0L) mean(rates_v)            else NA_real_,
    median_rate         = if (length(rates_v) > 0L) stats::median(rates_v)   else NA_real_,
    min_rate            = if (length(rates_v) > 0L) min(rates_v)             else NA_real_,
    max_rate            = if (length(rates_v) > 0L) max(rates_v)             else NA_real_,
    min_site            = if (length(rates_v) > 0L) as.integer(sites_v[which.min(rates_v)]) else NA_integer_,
    max_site            = if (length(rates_v) > 0L) as.integer(sites_v[which.max(rates_v)]) else NA_integer_,
    n_sites             = nrow(per_site),
    n_sites_valid       = sum(valid),
    threshold           = threshold,
    total_noise_events  = as.integer(total_noise),
    total_observations  = as.integer(total_obs)
  )

  list(
    per_site = per_site,
    summary  = summary_list,
    params   = list(
      threshold        = threshold,
      aa_levels        = aa_levels,
      position_columns = as.integer(position_columns),
      sites_evaluated  = as.integer(sites),
      drop_columns     = drop_columns
    )
  )
}

# -----------------------------------------------------------------------------
# summarize_empirical_pdel
# -----------------------------------------------------------------------------
# Convenience helper that turns the summary list into a 2-column data frame
# suitable for kable display. Used by the README chunk in section 6.2.
#
# Inputs
#   summary_list : the $summary slot from empirical_pdel().
#   sci_digits   : integer. Number of significant digits in the displayed
#                  scientific-notation rate values. Default 3L.
#
# Returns: data frame with columns Statistic and Value, both character.
# -----------------------------------------------------------------------------
summarize_empirical_pdel <- function(summary_list, sci_digits = 3L) {
  fmt <- function(x) {
    if (is.na(x)) return("NA")
    formatC(x, format = "e", digits = sci_digits - 1L)
  }
  data.frame(
    Statistic = c(
      "Pooled rate (total noise / total observations)",
      "Mean per-site rate",
      "Median per-site rate",
      sprintf("Minimum per-site rate (site %s)",
              as.character(summary_list$min_site)),
      sprintf("Maximum per-site rate (site %s)",
              as.character(summary_list$max_site))
    ),
    Value = c(
      fmt(summary_list$pooled_rate),
      fmt(summary_list$mean_rate),
      fmt(summary_list$median_rate),
      fmt(summary_list$min_rate),
      fmt(summary_list$max_rate)
    ),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
}
