################################################################################
## helpers_fda.R — FDA_Analysis
##
## Core science layer:
##   • Per-cell window construction with coverage clipping
##   • Window slicing with site-first column reordering
##   • Custom strategy-aware partitioner (sliding / disjoint / cumulative)
##   • Cell-level GMM with class-1 site extraction
##   • Per-site entropy curves and Hellinger trajectories on class-1 sites
##   • log1p transformation
##   • Safe fdahclust wrapper
##   • Mean silhouette via L2 dissimilarities
##   • Silhouette-driven n_clusters selection over a grid
##   • Hennig (2007) bootstrap Jaccard stability
##   • Fisher exact SNP enrichment with Bonferroni correction
##
## Sourced by run_one_cell.R and (read-only helpers) by aggregate_tables.R.
################################################################################

# Pre-requisite: setup.R has been sourced in this R session.


# ── 1. Per-cell window construction ──────────────────────────────────────────
#
# Window = [detection − K_pre·m, detection + K_post·m] where m = REF_BIN_MONTHS,
# clipped to dataset coverage. Returns the clipped window and the realised
# K_pre / K_post (an integer number of REF_BIN_MONTHS-bins actually contained
# on each side).
#
# Out-of-coverage logic: a detection is out-of-coverage iff its month is
# strictly before the dataset's first month or strictly after the dataset's
# last month. Edge months are kept (clipping handles them).

build_cell_window <- function(detection_date, fm_start, fm_end,
                              K_pre, K_post, ref_bin_months) {

  if (any(is.na(c(detection_date, fm_start, fm_end))))
    stop("NA dates passed to build_cell_window().", call. = FALSE)

  d_snap     <- snap_to_month(detection_date)
  fm_s_snap  <- snap_to_month(fm_start)
  fm_e_snap  <- snap_to_month(fm_end)
  fm_e_excl  <- add_months(fm_e_snap, 1L)            # exclusive upper bound

  if (d_snap < fm_s_snap || d_snap >= fm_e_excl) {
    return(list(
      window_start    = as.Date(NA),
      window_end      = as.Date(NA),
      K_pre_realised  = NA_integer_,
      K_post_realised = NA_integer_,
      status          = "out_of_coverage",
      reason          = sprintf(
        "Detection %s outside dataset coverage [%s, %s).",
        format(d_snap), format(fm_s_snap), format(fm_e_excl))
    ))
  }

  desired_start <- add_months(d_snap, -K_pre  * ref_bin_months)
  desired_end   <- add_months(d_snap,  K_post * ref_bin_months)

  win_s <- max(desired_start, fm_s_snap)
  win_e <- min(desired_end,   fm_e_excl)

  if (win_e <= win_s)
    return(list(
      window_start    = as.Date(NA),
      window_end      = as.Date(NA),
      K_pre_realised  = NA_integer_,
      K_post_realised = NA_integer_,
      status          = "out_of_coverage",
      reason          = "Clipped window has zero or negative width."
    ))

  K_pre_realised  <- diff_months(d_snap,  win_s) %/% ref_bin_months
  K_post_realised <- diff_months(win_e,   d_snap) %/% ref_bin_months

  list(
    window_start    = win_s,
    window_end      = win_e,
    K_pre_realised  = as.integer(K_pre_realised),
    K_post_realised = as.integer(K_post_realised),
    status          = "ok",
    reason          = NA_character_
  )
}


# ── 2. Window slicing with site-first column reordering ─────────────────────
#
# Returns a data.frame containing only rows in [window_start, window_end), with
# columns reordered to (sites_1..n_sites, Date, [Country if present]). This is
# the canonical layout for both partition_time_windows and our custom
# partitioner.

slice_window <- function(feature_matrix, window_start, window_end,
                         n_sites = 1273L) {
  in_window <- feature_matrix$Date >= window_start &
               feature_matrix$Date <  window_end
  fm_slice  <- feature_matrix[in_window, , drop = FALSE]

  site_cols  <- setdiff(colnames(fm_slice), c("Date", "Country"))
  if (length(site_cols) != n_sites)
    stop("Window slice has ", length(site_cols),
         " site columns; expected ", n_sites, ".", call. = FALSE)
  other_cols <- intersect(c("Date", "Country"), colnames(fm_slice))

  out <- fm_slice[, c(site_cols, other_cols), drop = FALSE]
  rownames(out) <- NULL
  out
}


# ── 3. Custom strategy-aware partitioner ────────────────────────────────────
#
# Returns the list of per-bin site-only data frames, bin date metadata, and
# sequence counts. Skips the per-bin GMM that partition_time_windows performs
# (which we don't need here — we have the cell-level GMM).

fda_partition <- function(fm_slice, window_start, window_end, strategy_def,
                          n_sites = 1273L) {

  wl <- as.integer(strategy_def$window_length)
  wt <- as.integer(strategy_def$window_type)

  total_months <- diff_months(window_end, window_start)
  if (total_months < wl)
    return(list(Partitions = list(), Dates_Labels = character(0L),
                bin_starts = as.Date(character(0L)),
                bin_ends   = as.Date(character(0L)),
                bin_counts = integer(0L), n_bins = 0L))

  n_chunks <- switch(as.character(wt),
                     "1" = total_months %/% wl,                # cumulative
                     "2" = total_months -  wl + 1L,            # sliding
                     "3" = total_months %/% wl)                # disjoint
  if (is.null(n_chunks) || n_chunks < 2L)
    return(list(Partitions = list(), Dates_Labels = character(0L),
                bin_starts = as.Date(character(0L)),
                bin_ends   = as.Date(character(0L)),
                bin_counts = integer(0L), n_bins = 0L))

  partitions   <- vector("list", n_chunks)
  bin_starts   <- as.Date(rep(NA_real_, n_chunks), origin = "1970-01-01")
  bin_ends     <- as.Date(rep(NA_real_, n_chunks), origin = "1970-01-01")
  dates_labels <- character(n_chunks)
  bin_counts   <- integer(n_chunks)

  for (i in seq_len(n_chunks)) {
    if (wt == 1L) {
      win_s <- window_start
      win_e <- add_months(window_start, wl * i)
    } else if (wt == 2L) {
      win_s <- add_months(window_start, i - 1L)
      win_e <- add_months(win_s, wl)
    } else {
      win_s <- add_months(window_start, wl * (i - 1L))
      win_e <- add_months(win_s, wl)
    }

    bin_starts[i] <- win_s
    bin_ends[i]   <- add_months(win_e, -1L)   # inclusive last month

    in_bin <- fm_slice$Date >= win_s &
              fm_slice$Date <  win_e &
              fm_slice$Date <  window_end
    chunk <- fm_slice[in_bin, , drop = FALSE]
    partitions[[i]] <- chunk[, seq_len(n_sites), drop = FALSE]
    bin_counts[i]   <- nrow(chunk)

    dates_labels[i] <- if (wt == 1L)
      sprintf("%s\u2197%s", format(window_start, "%b-%y"),
              format(bin_ends[i], "%b-%y"))
    else
      sprintf("%s\u2013%s", format(win_s,         "%b-%y"),
              format(bin_ends[i], "%b-%y"))
  }

  names(partitions) <- paste0("T", seq_len(n_chunks))

  list(
    Partitions   = partitions,
    Dates_Labels = dates_labels,
    bin_starts   = bin_starts,
    bin_ends     = bin_ends,
    bin_counts   = bin_counts,
    n_bins       = as.integer(n_chunks)
  )
}


# ── 4. Cell-level GMM (per-cell window; shared semantics across strategies) ──
#
# Fits a Gaussian Mixture Model on the per-site entropies computed on the
# entire cell window (not per-strategy bins). Extracts class-1 sites (the
# highest-entropy class after descending-mean relabelling).
#
# Returns either a list with non-empty $class1_sites, or a list with
# $reduced_skipped = TRUE explaining why FDA cannot proceed.

fit_cell_gmm <- function(fm_slice, n_sites, config) {

  if (nrow(fm_slice) < config$MIN_SEQUENCES_PER_BIN)
    return(list(class1_sites    = integer(0L),
                reduced_skipped = TRUE,
                G               = NA_integer_,
                modelName       = NA_character_,
                n_class1_sites  = 0L,
                reason          = sprintf("Window has %d < %d sequences.",
                                          nrow(fm_slice),
                                          config$MIN_SEQUENCES_PER_BIN)))

  seq_mat <- as.matrix(fm_slice[, seq_len(n_sites), drop = FALSE])
  entropies <- apply(seq_mat, 2L, ViralEntropR::calculate_entropy)

  args <- list(
    entropies  = entropies,
    nr         = nrow(seq_mat),
    nsites     = n_sites,
    G          = config$MCLUST_G,
    removez    = TRUE,
    removesngl = TRUE,
    verbose    = FALSE
  )
  if (!is.null(config$MCLUST_MODELS))
    args$modelNames <- config$MCLUST_MODELS

  fit <- tryCatch(do.call(ViralEntropR::cluster_sites_by_entropy, args),
                  error = function(e)
                    structure(list(error = conditionMessage(e)),
                              class = "fit_error"))

  if (inherits(fit, "fit_error"))
    return(list(class1_sites    = integer(0L),
                reduced_skipped = TRUE,
                G               = NA_integer_,
                modelName       = NA_character_,
                n_class1_sites  = 0L,
                reason          = paste0("cluster_sites_by_entropy() failed: ",
                                         fit$error)))

  df <- fit$DataFrame
  if (is.null(df) || nrow(df) == 0L)
    return(list(class1_sites    = integer(0L),
                reduced_skipped = TRUE,
                G               = NA_integer_,
                modelName       = NA_character_,
                n_class1_sites  = 0L,
                reason          = "GMM returned no sites after filtering."))

  if (any(df$class == 999L))
    return(list(class1_sites    = integer(0L),
                reduced_skipped = TRUE,
                G               = NA_integer_,
                modelName       = NA_character_,
                n_class1_sites  = 0L,
                reason          = "GMM 999 sentinel: entropies all identical."))

  df_rl <- ViralEntropR::relabel_entropy_classes(df)
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
    n_class1_sites  = length(class1_sites),
    reason          = if (length(class1_sites) == 0L)
                        "Class-1 set is empty after relabelling." else NA_character_
  )
}


# ── 5. Per-site entropy matrix on class-1 sites ─────────────────────────────
#
# Returns a numeric matrix with rows = class-1 sites, cols = bins (T1..T_n).
# Row names = site indices (as character), col names = bin labels.

compute_entropy_matrix <- function(partition_result, class1_sites) {

  parts  <- partition_result$Partitions
  n_bins <- partition_result$n_bins
  if (n_bins == 0L)
    stop("compute_entropy_matrix: no partitions.", call. = FALSE)
  if (length(class1_sites) == 0L)
    stop("compute_entropy_matrix: empty class1_sites.", call. = FALSE)

  ent_mat <- matrix(0,
                    nrow = length(class1_sites), ncol = n_bins,
                    dimnames = list(as.character(class1_sites),
                                    names(parts)))

  for (k in seq_len(n_bins)) {
    chunk <- parts[[k]]
    if (nrow(chunk) == 0L) {
      ent_mat[, k] <- 0
      next
    }
    sub <- chunk[, class1_sites, drop = FALSE]
    ent_mat[, k] <- apply(sub, 2L, ViralEntropR::calculate_entropy)
  }

  ent_mat
}


# ── 6. Per-site Hellinger matrix on class-1 sites ───────────────────────────
#
# Wraps ViralEntropR::calculate_hellinger_matrix and returns a numeric matrix
# of shape (n_class1_sites × (n_bins - 1)). The reference bin is partitions
# [[reference_idx]]; the column for bin T_k is D_H(T_k, T_ref).

compute_per_site_hellinger <- function(partition_result, class1_sites,
                                       aa_levels = 25L, normalized = FALSE,
                                       reference_idx = 1L) {

  parts <- partition_result$Partitions
  n     <- length(parts)
  if (n < 2L)
    stop("compute_per_site_hellinger: need >= 2 partitions.", call. = FALSE)
  if (reference_idx < 1L || reference_idx > n)
    stop("reference_idx ", reference_idx, " out of range [1, ", n, "].",
         call. = FALSE)

  # Reorder so the reference is first (calculate_hellinger_matrix takes T1 as ref).
  if (reference_idx != 1L) {
    new_order <- c(reference_idx, setdiff(seq_len(n), reference_idx))
    parts    <- parts[new_order]
    new_labs <- paste0("T", new_order)
    names(parts) <- new_labs
  } else {
    new_labs <- names(parts)
  }

  hel <- ViralEntropR::calculate_hellinger_matrix(
    partitions = parts,
    sites      = class1_sites,
    aa_levels  = aa_levels,
    normalized = normalized,
    labels     = new_labs
  )

  if (!is.matrix(hel))
    stop("calculate_hellinger_matrix did not return a matrix.", call. = FALSE)
  if (nrow(hel) != length(class1_sites))
    stop("Hellinger matrix has ", nrow(hel),
         " rows; expected ", length(class1_sites), ".", call. = FALSE)
  if (ncol(hel) != n - 1L)
    stop("Hellinger matrix has ", ncol(hel),
         " cols; expected ", n - 1L, ".", call. = FALSE)

  rownames(hel) <- as.character(class1_sites)
  hel
}


# ── 7. log1p transformation ─────────────────────────────────────────────────

log1p_transform <- function(M) {
  if (!is.matrix(M))
    stop("log1p_transform: input must be a matrix.", call. = FALSE)
  if (any(M < 0, na.rm = TRUE))
    stop("log1p_transform: negative values encountered.", call. = FALSE)
  log1p(M)
}


# ── 8. Safe fdahclust wrapper ───────────────────────────────────────────────
#
# fdahclust signature (verified against fdacluster docs):
#   x = numeric grid vector of length M
#   y = numeric matrix of shape N × M
#   returns a `caps` object with $memberships (integer length N)
#
# We pin: warping_class="none", metric="l2", centroid_type="mean",
# linkage_criterion="complete", use_verbose=FALSE.

fit_fdahclust_safe <- function(curve_matrix, grid_vec, n_clusters, config) {

  n_curves <- nrow(curve_matrix)
  n_points <- ncol(curve_matrix)

  if (n_curves < n_clusters + 1L)
    stop("fdahclust: need >= n_clusters+1 curves; got n=", n_curves,
         " for k=", n_clusters, ".", call. = FALSE)
  if (n_points < 2L)
    stop("fdahclust: need >= 2 grid points; got M=", n_points, ".",
         call. = FALSE)
  if (length(grid_vec) != n_points)
    stop("fdahclust: grid_vec length (", length(grid_vec),
         ") does not match curve_matrix ncol (", n_points, ").",
         call. = FALSE)
  if (anyNA(curve_matrix) || anyNA(grid_vec))
    stop("fdahclust: NA in curve_matrix or grid_vec.", call. = FALSE)

  fit <- tryCatch(
    fdacluster::fdahclust(
      x                 = as.numeric(grid_vec),
      y                 = curve_matrix,
      n_clusters        = as.integer(n_clusters),
      warping_class     = config$FDA_WARPING_CLASS,
      metric            = config$FDA_METRIC,
      centroid_type     = config$FDA_CENTROID_TYPE,
      linkage_criterion = config$FDA_LINKAGE,
      use_verbose       = FALSE
    ),
    error = function(e)
      stop("fdahclust failed for k=", n_clusters, ": ",
           conditionMessage(e), call. = FALSE)
  )

  if (is.null(fit$memberships) || length(fit$memberships) != n_curves)
    stop("fdahclust returned malformed memberships.", call. = FALSE)

  fit
}


# ── 9. Mean silhouette via L2 dissimilarity ─────────────────────────────────

mean_silhouette <- function(curve_matrix, memberships) {

  K <- length(unique(memberships))
  if (K < 2L) return(NA_real_)

  D <- as.matrix(stats::dist(curve_matrix, method = "euclidean"))
  sil <- tryCatch(cluster::silhouette(as.integer(memberships), dmatrix = D),
                  error = function(e) NULL)
  if (is.null(sil) || !is.matrix(sil) || !("sil_width" %in% colnames(sil)))
    return(NA_real_)
  mean(sil[, "sil_width"], na.rm = TRUE)
}


# ── 10. Silhouette-driven n_clusters selection ──────────────────────────────
#
# Fits fdahclust at each feasible k ∈ k_grid and returns the k that maximises
# mean silhouette width. A k is feasible iff n_curves > k.

select_n_clusters <- function(curve_matrix, grid_vec, k_grid, config) {

  n_curves <- nrow(curve_matrix)
  feas     <- k_grid[k_grid + 1L <= n_curves]
  if (length(feas) == 0L)
    stop("select_n_clusters: no feasible k for n_curves=", n_curves, ".",
         call. = FALSE)

  fits <- vector("list", length(feas))
  sils <- numeric(length(feas))

  for (i in seq_along(feas)) {
    k   <- feas[i]
    fit <- tryCatch(fit_fdahclust_safe(curve_matrix, grid_vec, k, config),
                    error = function(e) NULL)
    fits[[i]] <- fit
    sils[i]   <- if (is.null(fit)) NA_real_ else
                 mean_silhouette(curve_matrix, as.integer(fit$memberships))
  }

  if (all(is.na(sils)))
    stop("select_n_clusters: all fdahclust fits failed across k grid.",
         call. = FALSE)

  best <- which.max(sils)
  list(
    k                = as.integer(feas[best]),
    fit              = fits[[best]],
    silhouette_mean  = unname(sils[best]),
    silhouette_per_k = setNames(sils, paste0("k=", feas)),
    k_grid_feasible  = as.integer(feas)
  )
}


# ── 11. Hennig (2007) bootstrap Jaccard stability ───────────────────────────
#
# For each of B bootstrap resamples of the N curves (with replacement):
#   1. Refit fdahclust at the SAME n_clusters as the original fit.
#   2. Map each original curve to a bootstrap cluster via its FIRST occurrence
#      in the bootstrap sample (Hennig's prescription). Curves not sampled
#      contribute NA.
#   3. For each original cluster, record the max Jaccard against any
#      bootstrap cluster.
# Per-cluster stability = mean over the B max-Jaccards.
# A cluster is "stable" iff its mean Jaccard >= STABLE_JACCARD_THRESHOLD.

bootstrap_jaccard_stability <- function(curve_matrix, grid_vec, n_clusters,
                                        B, seed, config) {

  set.seed(seed)
  N <- nrow(curve_matrix)

  original_fit <- fit_fdahclust_safe(curve_matrix, grid_vec, n_clusters, config)
  orig_mem <- as.integer(original_fit$memberships)

  per_cluster_jaccard <- matrix(NA_real_, nrow = B, ncol = n_clusters)
  n_failed <- 0L

  for (b in seq_len(B)) {
    idx    <- sample.int(N, size = N, replace = TRUE)
    boot_y <- curve_matrix[idx, , drop = FALSE]

    boot_fit <- tryCatch(
      fit_fdahclust_safe(boot_y, grid_vec, n_clusters, config),
      error = function(e) NULL
    )
    if (is.null(boot_fit)) { n_failed <- n_failed + 1L; next }

    boot_mem <- as.integer(boot_fit$memberships)

    # First-occurrence mapping: for each original site s, find the position of
    # its first appearance in `idx` (NA if s never sampled).
    first_occ <- match(seq_len(N), idx)
    orig_to_boot <- ifelse(is.na(first_occ),
                            NA_integer_,
                            boot_mem[first_occ])

    for (k in seq_len(n_clusters)) {
      orig_in_k <- which(orig_mem == k)
      if (length(orig_in_k) == 0L) next

      best_jacc <- 0
      for (l in seq_len(n_clusters)) {
        boot_in_l <- which(orig_to_boot == l)
        if (length(boot_in_l) == 0L) next
        n_int   <- length(intersect(orig_in_k, boot_in_l))
        n_union <- length(union(orig_in_k, boot_in_l))
        if (n_union > 0L) best_jacc <- max(best_jacc, n_int / n_union)
      }
      per_cluster_jaccard[b, k] <- best_jacc
    }
  }

  mean_jacc <- apply(per_cluster_jaccard, 2L, function(v) {
    if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  })

  list(
    n_bootstrap         = as.integer(B),
    n_failed_bootstrap  = as.integer(n_failed),
    per_cluster_jaccard = mean_jacc,
    n_stable_clusters   = sum(mean_jacc >= config$STABLE_JACCARD_THRESHOLD,
                               na.rm = TRUE),
    original_memberships = orig_mem
  )
}


# ── 12. Fisher exact SNP enrichment on class-1 sites ────────────────────────
#
# Tests whether SNP sites are differentially distributed across clusters.
# Contingency table is 2 (snp / non-snp) × n_clusters; Fisher exact in R
# handles arbitrary 2×K via simulated p when K > 2. Bonferroni denominator
# is fixed at N_TOTAL_ENRICHMENT_TESTS (29 cells × 2 metrics + spare = 60).

test_snp_enrichment <- function(class1_sites, cluster_assignments, snp_sites,
                                n_total_tests) {

  if (length(class1_sites) != length(cluster_assignments))
    stop("class1_sites / cluster_assignments length mismatch.", call. = FALSE)

  is_snp <- class1_sites %in% snp_sites
  clust  <- as.integer(cluster_assignments)

  ctbl <- table(
    snp     = factor(is_snp, levels = c("FALSE", "TRUE")),
    cluster = factor(clust,   levels = sort(unique(clust)))
  )

  blank <- list(
    fisher_p                = NA_real_,
    fisher_p_bonferroni     = NA_real_,
    n_snp_in_class1         = sum(is_snp),
    n_class1_total          = length(class1_sites),
    n_snp_total_truth       = length(snp_sites),
    max_overlap_cluster     = NA_integer_,
    max_overlap_OR          = NA_real_,
    snp_overlap_per_cluster = setNames(integer(0L), character(0L)),
    contingency_table       = ctbl,
    status                  = "insufficient_variation"
  )

  if (nrow(ctbl) < 2L || ncol(ctbl) < 2L) return(blank)

  ft <- tryCatch(
    if (ncol(ctbl) == 2L) fisher.test(ctbl)
    else                  fisher.test(ctbl, simulate.p.value = TRUE,
                                       B = 10000L),
    error = function(e) NULL
  )

  if (is.null(ft)) {
    blank$status <- "fisher_failed"
    return(blank)
  }

  per_cluster_snp_n <- as.integer(ctbl["TRUE", , drop = TRUE])
  names(per_cluster_snp_n) <- colnames(ctbl)
  per_cluster_total <- as.integer(colSums(ctbl))
  per_cluster_rate  <- per_cluster_snp_n / pmax(per_cluster_total, 1L)
  max_clust         <- as.integer(colnames(ctbl)[which.max(per_cluster_rate)])

  # 2×2 odds ratio for max_clust vs rest.
  in_max <- clust == max_clust
  ctbl_2x2 <- table(
    snp    = factor(is_snp, levels = c("FALSE", "TRUE")),
    in_max = factor(in_max, levels = c("FALSE", "TRUE"))
  )
  or_estimate <- if (all(dim(ctbl_2x2) == 2L)) {
    ft2 <- tryCatch(fisher.test(ctbl_2x2), error = function(e) NULL)
    if (is.null(ft2)) NA_real_ else unname(ft2$estimate)
  } else NA_real_

  list(
    fisher_p                = unname(ft$p.value),
    fisher_p_bonferroni     = min(1, unname(ft$p.value) * n_total_tests),
    n_snp_in_class1         = sum(is_snp),
    n_class1_total          = length(class1_sites),
    n_snp_total_truth       = length(snp_sites),
    max_overlap_cluster     = max_clust,
    max_overlap_OR          = or_estimate,
    snp_overlap_per_cluster = per_cluster_snp_n,
    contingency_table       = ctbl,
    status                  = "ok"
  )
}


# ── 13. Detection-bin lookup (for plotting) ─────────────────────────────────

detection_bin_index <- function(detection_date, bin_starts, bin_ends) {
  d_snap <- snap_to_month(detection_date)
  idx <- which(bin_starts <= d_snap & d_snap <= bin_ends)
  if (length(idx) == 0L) return(NA_integer_)
  as.integer(idx[1L])
}


# ── 14. Anchor-bin resolution for dual Hellinger anchors ────────────────────
#
# Given an anchor Date (e.g., a predecessor variant's US detection), return the
# partition index in the cell whose [bin_starts, bin_ends] contains the anchor,
# or — if the anchor is outside the window — the nearest in-window bin (with a
# flag indicating the snap).
#
# Returns: list(bin_idx, snapped, distance_months).
#   bin_idx       integer 1..n_bins
#   snapped       TRUE if the anchor fell outside [bin_starts[1], bin_ends[n]]
#   distance_months  abs(months) from anchor to bin centre (0 if direct hit)

anchor_bin_idx <- function(anchor_date, bin_starts, bin_ends) {

  if (length(bin_starts) == 0L || length(bin_ends) != length(bin_starts))
    stop("anchor_bin_idx: bad bin metadata.", call. = FALSE)
  if (is.na(anchor_date))
    return(list(bin_idx = NA_integer_, snapped = NA, distance_months = NA_real_))

  a_snap <- snap_to_month(anchor_date)
  n_bins <- length(bin_starts)

  # Direct hit?
  hit <- which(bin_starts <= a_snap & a_snap <= bin_ends)
  if (length(hit) > 0L)
    return(list(bin_idx        = as.integer(hit[1L]),
                snapped        = FALSE,
                distance_months = 0))

  # Outside window: snap to nearest bin by month distance to bin start.
  d <- vapply(seq_len(n_bins),
              function(i) abs(diff_months(a_snap, bin_starts[i])),
              numeric(1L))
  best <- which.min(d)
  list(bin_idx        = as.integer(best),
       snapped        = TRUE,
       distance_months = as.numeric(d[best]))
}


# ── 15. compute_per_site_hellinger_anchor — arbitrary anchor index ──────────
#
# Variant of compute_per_site_hellinger() that takes the anchor bin index
# explicitly rather than a hard-coded reference_idx of 1. Returns a matrix of
# shape (n_class1_sites × (n_bins - 1)), with rows the class-1 sites and
# columns the non-anchor bins in original order (T1..T_(anchor-1), T_(anchor+1)..T_n).
#
# The column labels are paste0("T", non_anchor_idx) so callers can map back to
# real bin indices for plotting.

compute_per_site_hellinger_anchor <- function(partition_result, class1_sites,
                                              anchor_idx,
                                              aa_levels = 25L,
                                              normalized = FALSE) {

  parts <- partition_result$Partitions
  n     <- length(parts)
  if (n < 2L)
    stop("compute_per_site_hellinger_anchor: need >= 2 partitions.", call. = FALSE)
  if (anchor_idx < 1L || anchor_idx > n)
    stop("anchor_idx ", anchor_idx, " out of range [1, ", n, "].", call. = FALSE)

  # Reorder so anchor is first; remember original positions for relabelling.
  new_order <- c(anchor_idx, setdiff(seq_len(n), anchor_idx))
  parts_re  <- parts[new_order]
  labs_re   <- paste0("T", new_order)
  names(parts_re) <- labs_re

  hel <- ViralEntropR::calculate_hellinger_matrix(
    partitions = parts_re,
    sites      = class1_sites,
    aa_levels  = aa_levels,
    normalized = normalized,
    labels     = labs_re
  )

  if (!is.matrix(hel))
    stop("calculate_hellinger_matrix did not return a matrix.", call. = FALSE)
  if (nrow(hel) != length(class1_sites))
    stop("Hellinger matrix has ", nrow(hel),
         " rows; expected ", length(class1_sites), ".", call. = FALSE)
  if (ncol(hel) != n - 1L)
    stop("Hellinger matrix has ", ncol(hel),
         " cols; expected ", n - 1L, ".", call. = FALSE)

  # Column order from calculate_hellinger_matrix matches labs_re[-1]. Re-order
  # columns to follow original bin index (skipping the anchor).
  non_anchor <- setdiff(seq_len(n), anchor_idx)
  reorder    <- order(non_anchor)
  hel        <- hel[, reorder, drop = FALSE]
  colnames(hel) <- paste0("T", non_anchor[reorder])
  rownames(hel) <- as.character(class1_sites)
  hel
}


# ── 16. Per-site silhouette via L2 dissimilarity ────────────────────────────
#
# Returns the per-site silhouette widths as a named numeric vector aligned with
# rows of curve_matrix. NA-handling: sites with NA membership get NA silhouette.

per_site_silhouette <- function(curve_matrix, memberships) {

  K <- length(unique(memberships[!is.na(memberships)]))
  if (K < 2L)
    return(setNames(rep(NA_real_, nrow(curve_matrix)),
                    rownames(curve_matrix)))

  D <- as.matrix(stats::dist(curve_matrix, method = "euclidean"))
  sil <- tryCatch(cluster::silhouette(as.integer(memberships), dmatrix = D),
                  error = function(e) NULL)
  if (is.null(sil) || !is.matrix(sil) || !("sil_width" %in% colnames(sil)))
    return(setNames(rep(NA_real_, nrow(curve_matrix)),
                    rownames(curve_matrix)))

  setNames(as.numeric(sil[, "sil_width"]), rownames(curve_matrix))
}
