# =============================================================================
# detect_in_snapshot.R
# =============================================================================
#
# Detection logic for one GMM fit on one population snapshot.
#
# Contains the inner-loop helpers used by run_one_replicate.R during the
# sweep over candidate n_emerge values. The entropy computation is
# vectorised at the column level via base::tabulate(), avoiding the slower
# table()-based path inside ViralEntropR::calculate_entropy when called via
# apply(., 2, calculate_entropy).
#
# Author : Vadim Tyuryaev
# =============================================================================

# -----------------------------------------------------------------------------
# compute_entropy_per_site
# -----------------------------------------------------------------------------
# Per-site Shannon entropy across rows of an integer matrix.
#
# Inputs
#   mat       : integer matrix with rows = sequences, columns = sites
#               and values in 1..aa_levels.
#   n_rows    : integer scalar. Number of leading rows to include in the
#               computation. Use values < nrow(mat) to compute entropy over
#               a row prefix without copying the matrix; the function only
#               subsets one column at a time (a cheap copy of length n_rows).
#   aa_levels : integer. Alphabet size (default 25, matching the package
#               convention). All codes outside 1..aa_levels contribute zero
#               to tabulate() so they neither inflate nor zero out entropy.
#   base      : numeric. Log base (default 2). Same default as
#               ViralEntropR::calculate_entropy.
#
# Returns: numeric vector of length ncol(mat). Sites with fewer than two
# distinct residues across the n_rows rows have entropy 0.
#
# Implementation notes
#   - tabulate() is implemented in C and is significantly faster than
#     table() on integer vectors.
#   - The per-column loop is unavoidable in pure R without inflating memory
#     to an n_rows x n_sites x aa_levels indicator array. The work per
#     column is O(n_rows) in C, so total cost is O(n_rows * n_sites) with a
#     small constant.
# -----------------------------------------------------------------------------
compute_entropy_per_site <- function(mat,
                                     n_rows    = nrow(mat),
                                     aa_levels = 25L,
                                     base      = 2) {
  if (!is.matrix(mat) || !is.integer(mat))
    stop("`mat` must be an integer matrix.", call. = FALSE)

  n_rows    <- as.integer(n_rows)
  aa_levels <- as.integer(aa_levels)
  if (n_rows < 1L || n_rows > nrow(mat))
    stop("`n_rows` out of bounds: ", n_rows, call. = FALSE)

  ncols <- ncol(mat)
  out   <- numeric(ncols)
  idx   <- if (n_rows == nrow(mat)) NULL else seq_len(n_rows)
  inv_n <- 1 / n_rows
  log_base <- log(base)

  for (j in seq_len(ncols)) {
    col <- if (is.null(idx)) mat[, j] else mat[idx, j]
    counts <- tabulate(col, nbins = aa_levels)
    nz <- counts > 0L
    if (sum(nz) <= 1L) {
      out[j] <- 0
    } else {
      probs <- counts[nz] * inv_n
      out[j] <- -sum(probs * log(probs)) / log_base
    }
  }

  out
}

# -----------------------------------------------------------------------------
# .call_cluster_sites
# -----------------------------------------------------------------------------
# Thin wrapper around ViralEntropR::cluster_sites_by_entropy that omits the
# modelNames argument when it is NULL, so that mclust's default model
# selection across both equal-variance ("E") and variable-variance ("V")
# univariate Gaussian mixtures is invoked.
#
# Forcing modelNames = "V" was found to cause BIC degeneracy when the input
# entropy vector contains exact ties: the V model's per-component variance
# is unidentified at zero within-component variance, leaving G = 1 as the
# only well-defined BIC value, which is then selected by default. The E
# model is well-defined under the same conditions because its single
# pooled-variance parameter is positive whenever any component has spread.
# The simulator's discrete population structure produces exact ties at
# every mutation site of any single variant, so the default search across
# both E and V is essential to obtain biologically correct cluster
# assignments.
# -----------------------------------------------------------------------------
.call_cluster_sites <- function(entropies, n_rows, mclust_models, mclust_G) {
  args <- list(
    entropies, nr = n_rows,
    G          = mclust_G,
    removez    = TRUE,
    removesngl = TRUE,
    verbose    = FALSE
  )
  if (!is.null(mclust_models)) {
    args$modelNames <- mclust_models
  }
  do.call(ViralEntropR::cluster_sites_by_entropy, args)
}

# -----------------------------------------------------------------------------
# .extract_fit_metadata
# -----------------------------------------------------------------------------
# Pulls the (modelName, G) chosen by BIC from a cluster_sites_by_entropy
# return value. Returns NA for both fields if the FitObject is missing or
# malformed (e.g., the all-identical sentinel path returns without an
# Mclust fit). Used for audit trail in the per-replicate RDS.
# -----------------------------------------------------------------------------
.extract_fit_metadata <- function(cl) {
  fit <- cl$FitObject
  if (is.null(fit) || !is.list(fit)) {
    return(list(modelName_chosen = NA_character_, G_chosen = NA_integer_))
  }
  mn <- fit$modelName
  g  <- fit$G
  list(
    modelName_chosen = if (is.null(mn)) NA_character_ else as.character(mn),
    G_chosen         = if (is.null(g))  NA_integer_   else as.integer(g)
  )
}

# -----------------------------------------------------------------------------
# test_detection
# -----------------------------------------------------------------------------
# Performs the study's detection test on one population snapshot.
#
# Inputs
#   mat              : integer matrix from simulate_population_snapshot()
#                      (post deleterious-noise application if any).
#   n_rows           : integer. Leading rows to consider.
#   emerge_positions : integer vector. Mutation positions of the emerging
#                      variant whose detection is being tested.
#   mclust_models    : character or NULL. modelNames forwarded to
#                      ViralEntropR::cluster_sites_by_entropy. When NULL the
#                      package default model search (both "E" and "V") is
#                      used; see .call_cluster_sites for rationale.
#   mclust_G         : integer vector. G forwarded to
#                      ViralEntropR::cluster_sites_by_entropy.
#
# Returns: a list with
#   $detected         : logical. TRUE iff every emerge_position is in
#                       class 1 of the relabelled GMM, with the two
#                       degenerate-path tightenings described below.
#   $entropy_at_sites : numeric vector. Per-site entropy at each
#                       emerge_position (from compute_entropy_per_site, which
#                       sees all positions regardless of GMM filtering).
#   $class_at_sites   : integer vector. Relabelled class at each
#                       emerge_position, or NA if the site is not in the
#                       GMM-returned DataFrame (e.g. filtered as invariant
#                       or singleton).
#   $n_class_1_sites  : integer. Number of sites in class 1 across the
#                       entire DataFrame (i.e., genome-wide class-1 count).
#   $modelName_chosen : character. BIC-selected model name ("E" or "V"),
#                       or NA when the all-identical sentinel path was
#                       taken or the fit object is missing.
#   $G_chosen         : integer. BIC-selected number of components, or NA
#                       under the same conditions as above.
#
# Detection semantics
#   - Multi-class path (G >= 2 after relabel): detection iff every
#     emerging position is in the highest-entropy class (class 1).
#   - Single-class path (G = 1): detection iff the surviving sites are
#     EXACTLY the emerging positions. Strict equality prevents G = 1
#     collapse from claiming detection when the single cluster is
#     contaminated with established-variant sites.
#   - All-identical sentinel path (class = 999L for every surviving site):
#     detection iff the surviving sites are EXACTLY the emerging
#     positions. Same defensive reasoning as the G = 1 path.
#
# Notes
#   - cluster_sites_by_entropy is called with the package defaults
#     removez = TRUE and removesngl = TRUE; positions with zero or
#     singleton entropy are dropped from the DataFrame and treated as
#     "not detected" for those sites.
#   - When the DataFrame is empty (all sites filtered), detection is
#     reported as FALSE.
# -----------------------------------------------------------------------------
test_detection <- function(mat,
                           n_rows,
                           emerge_positions,
                           mclust_models = NULL,
                           mclust_G      = 1:15) {

  emerge_positions <- as.integer(emerge_positions)

  # ----- 1. Per-site entropy on the row prefix -----------------------------
  entropies <- compute_entropy_per_site(mat, n_rows = n_rows, aa_levels = 25L)

  # Trivial failure: emerging-variant sites all carry zero entropy
  # (cannot occur under the simulator's invariants but guarded for safety).
  if (all(entropies[emerge_positions] == 0)) {
    return(list(
      detected         = FALSE,
      entropy_at_sites = entropies[emerge_positions],
      class_at_sites   = rep(NA_integer_, length(emerge_positions)),
      n_class_1_sites  = 0L,
      modelName_chosen = NA_character_,
      G_chosen         = NA_integer_
    ))
  }

  # ----- 2. GMM clustering -------------------------------------------------
  cl  <- .call_cluster_sites(entropies, n_rows, mclust_models, mclust_G)
  df  <- cl$DataFrame
  fit_meta <- .extract_fit_metadata(cl)

  if (nrow(df) == 0L) {
    return(list(
      detected         = FALSE,
      entropy_at_sites = entropies[emerge_positions],
      class_at_sites   = rep(NA_integer_, length(emerge_positions)),
      n_class_1_sites  = 0L,
      modelName_chosen = fit_meta$modelName_chosen,
      G_chosen         = fit_meta$G_chosen
    ))
  }

  # All-identical entropies sentinel. cluster_sites_by_entropy returns
  # class = 999L for every surviving site when entropies are bit-identical
  # across all retained positions. Detection succeeds iff the surviving
  # sites are EXACTLY the emerging positions. The strict equality prevents
  # false positives when the sentinel is reached but background-variant
  # sites also survived alongside the emerging ones.
  if (all(df$class == 999L)) {
    surviving_sites <- df$sites
    match_idx       <- match(emerge_positions, surviving_sites)
    detected_999    <- setequal(surviving_sites, emerge_positions)
    return(list(
      detected         = detected_999,
      entropy_at_sites = entropies[emerge_positions],
      class_at_sites   = ifelse(is.na(match_idx), NA_integer_, 999L),
      n_class_1_sites  = length(surviving_sites),
      modelName_chosen = fit_meta$modelName_chosen,
      G_chosen         = fit_meta$G_chosen
    ))
  }

  # Single-component result (G = 1 selected by BIC). Every variable site
  # is in one undifferentiated class. Detection succeeds iff the surviving
  # sites are EXACTLY the emerging positions. The strict equality prevents
  # G = 1 collapse from claiming detection when V1, V2, or V3 sites also
  # contaminate the single cluster.
  unique_classes <- unique(df$class)
  if (length(unique_classes) == 1L) {
    surviving_sites <- df$sites
    detected_g1     <- setequal(surviving_sites, emerge_positions)
    return(list(
      detected         = detected_g1,
      entropy_at_sites = entropies[emerge_positions],
      class_at_sites   = ifelse(emerge_positions %in% surviving_sites,
                                unique_classes,
                                NA_integer_),
      n_class_1_sites  = length(surviving_sites),
      modelName_chosen = fit_meta$modelName_chosen,
      G_chosen         = fit_meta$G_chosen
    ))
  }

  # ----- 3. Relabel so class 1 = highest entropy ---------------------------
  df <- ViralEntropR::relabel_entropy_classes(df)

  # ----- 4. Detection test -------------------------------------------------
  class_1_sites <- df$sites[df$class == 1L]
  detected      <- all(emerge_positions %in% class_1_sites)

  match_idx <- match(emerge_positions, df$sites)
  class_at_sites <- ifelse(is.na(match_idx), NA_integer_, df$class[match_idx])

  list(
    detected         = detected,
    entropy_at_sites = entropies[emerge_positions],
    class_at_sites   = class_at_sites,
    n_class_1_sites  = length(class_1_sites),
    modelName_chosen = fit_meta$modelName_chosen,
    G_chosen         = fit_meta$G_chosen
  )
}
