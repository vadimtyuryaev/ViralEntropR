################################################################################
## helpers_hellinger.R — Changepoint_Detection_Study
##
## Per-window Hellinger trajectory and per-window Gaussian-mixture fit.
##
## All functions assume integer-encoded feature matrices with the 25-residue
## alphabet (codes 1..25). aa_levels = 25L is a package invariant.
################################################################################

# ── 1. Partition a window's sequences into bins ──────────────────────────────
#
# Returns a named list of data.frames, one per bin, retaining only site
# columns (Date and Country dropped). Bin names T1, T2, … follow the
# convention used by ViralEntropR::calculate_hellinger_matrix().

partition_window <- function(feature_matrix, window_start, window_end,
                             site_cols, bin_months = 2L) {

  ranges <- bin_row_ranges(feature_matrix$Date, window_start, window_end,
                           bin_months)
  n_bins <- length(ranges)
  parts  <- vector("list", n_bins)
  for (i in seq_len(n_bins)) {
    idx <- ranges[[i]]
    parts[[i]] <- if (length(idx) == 0L) {
      feature_matrix[integer(0L), site_cols, drop = FALSE]
    } else {
      feature_matrix[idx, site_cols, drop = FALSE]
    }
  }
  names(parts) <- paste0("T", seq_len(n_bins))
  parts
}

# ── 2. Hellinger trajectory ──────────────────────────────────────────────────
#
# Returns the matrix as-is from calculate_hellinger_matrix: rows = sites,
# columns = T2, T3, …, T_n. Downstream code transposes to T × d before
# passing to CP methods.

compute_hellinger_full <- function(partitions, sites, aa_levels = 25L,
                                   normalized = FALSE) {
  if (length(partitions) < 2L)
    stop("Hellinger trajectory needs at least 2 partitions.", call. = FALSE)
  ViralEntropR::calculate_hellinger_matrix(
    partitions = partitions,
    sites      = sites,
    aa_levels  = aa_levels,
    normalized = normalized,
    labels     = paste0("T", seq_along(partitions))
  )
}

# ── 3. Per-window GMM fit and class-1 site extraction ────────────────────────
#
# Fits a GMM on per-site Shannon entropy over the WINDOW'S sequences only
# (no oracle leakage). Returns:
#
#   $class1_sites      : integer vector of site indices in class 1 after
#                        relabel_entropy_classes. May be integer(0L).
#   $reduced_skipped   : TRUE iff GMM returned the 999 sentinel.
#   $G                 : chosen number of components (or NA).
#   $modelName         : chosen Mclust model name (or NA).
#   $n_class1_sites    : length(class1_sites).
#
# Edge cases:
#   - 999 sentinel (all-identical entropy)     → reduced_skipped = TRUE.
#   - G = 1                                    → class 1 is the single class;
#                                                class1_sites = all sites that
#                                                passed the entropy filters.
#   - Empty filtered dataframe                 → reduced_skipped = TRUE.

fit_window_gmm <- function(window_seq_matrix, config) {

  if (nrow(window_seq_matrix) == 0L)
    stop("Cannot fit GMM: window has zero sequences.", call. = FALSE)

  entropies <- apply(window_seq_matrix, 2L, ViralEntropR::calculate_entropy)

  args <- list(
    entropies  = entropies,
    nr         = nrow(window_seq_matrix),
    G          = config$MCLUST_G,
    removez    = TRUE,
    removesngl = TRUE,
    verbose    = FALSE
  )
  if (!is.null(config$MCLUST_MODELS))
    args$modelNames <- config$MCLUST_MODELS

  fit <- do.call(ViralEntropR::cluster_sites_by_entropy, args)
  df  <- fit$DataFrame

  if (nrow(df) == 0L)
    return(list(class1_sites    = integer(0L),
                reduced_skipped = TRUE,
                G               = NA_integer_,
                modelName       = NA_character_,
                n_class1_sites  = 0L))

  if (any(df$class == 999L))
    return(list(class1_sites    = integer(0L),
                reduced_skipped = TRUE,
                G               = NA_integer_,
                modelName       = NA_character_,
                n_class1_sites  = 0L))

  df_rl        <- ViralEntropR::relabel_entropy_classes(df)
  class1_sites <- sort(as.integer(df_rl$sites[df_rl$class == 1L]))

  mn <- if (is.list(fit$FitObject) && !is.null(fit$FitObject$modelName))
          as.character(fit$FitObject$modelName) else NA_character_
  g  <- if (is.list(fit$FitObject) && !is.null(fit$FitObject$G))
          as.integer(fit$FitObject$G) else NA_integer_

  list(
    class1_sites    = class1_sites,
    reduced_skipped = length(class1_sites) == 0L,
    G               = g,
    modelName       = mn,
    n_class1_sites  = length(class1_sites)
  )
}
