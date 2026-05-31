################################################################################
## helpers_frames.R — FDA_Analysis (v2)
##
## Per-frame logic for the v2 animation pipeline. A "frame" is a slice of
## FITTING_WINDOW_HALFWIDTH partitions on each side of a centre partition.
## Each frame runs its own GMM on the slice's sequences, picks class-1 sites,
## fits fdahclust at the silhouette-best k over N_CLUSTERS_GRID, and records
## per-site silhouettes plus convergence flags.
##
## Three guards from the methodology (Modes A/B/C):
##   • Mode A: per-frame GMM picks G = 1 → no class-1, frame null.
##   • Mode B: best mean silhouette over the k-grid < MIN_MEAN_SILHOUETTE →
##     frame is structurally null (cluster IDs still recorded but flagged).
##   • Mode C: per-site silhouette < 0 → that site is "marginal", to be
##     rendered red in the bottom panel.
##
## Plus the v2-specific visual cap (G2):
##   • A cluster larger than MAX_CLUSTER_SIZE (= 33, Omicron's mutation count)
##     is rendered as background light grey in the panels (state recorded as
##     "oversized" per-site). Underlying tab_02 cluster_membership.csv is NOT
##     truncated — the cap is purely visualisation.
##
## Loaded by run_one_cell.R (per-frame loop) and by animate_results.R / 
## plot_results.R (frame replay for rendering).
################################################################################

# Pre-requisite: setup.R + helpers_fda.R sourced.


# ── 1. Frame-centre enumeration ─────────────────────────────────────────────
#
# Returns the set of valid frame centre indices for a cell with n_bins
# partitions and fitting half-width h. Each centre c yields a fitting slice
# (c-h)..(c+h). Centres are restricted to [h+1, n_bins-h].

frame_centre_indices <- function(n_bins, halfwidth) {
  if (n_bins < 2L * halfwidth + 1L) return(integer(0L))
  as.integer(seq.int(halfwidth + 1L, n_bins - halfwidth))
}


# ── 2. Per-frame GMM on a partition slice ───────────────────────────────────
#
# Pools all sequences in the slice's partitions, computes per-site entropy,
# fits Mclust over G = MCLUST_G, and extracts class-1 sites. Returns a
# structured list with explicit status.

fit_frame_gmm <- function(partition_result, frame_bins, n_sites, config,
                          frame_seed) {

  parts  <- partition_result$Partitions
  if (length(parts) == 0L)
    return(list(status = "no_partitions", class1_sites = integer(0L),
                G = NA_integer_, modelName = NA_character_,
                n_sequences = 0L, n_class1_sites = 0L,
                reason = "Partition result empty."))

  # Pool sequences across the frame's partitions.
  pooled <- do.call(rbind, parts[frame_bins])
  n_seq  <- nrow(pooled)

  if (n_seq < config$MIN_SEQUENCES_PER_BIN)
    return(list(status = "underpowered", class1_sites = integer(0L),
                G = NA_integer_, modelName = NA_character_,
                n_sequences = n_seq, n_class1_sites = 0L,
                reason = sprintf("Frame slice has %d < %d sequences.",
                                 n_seq, config$MIN_SEQUENCES_PER_BIN)))

  set.seed(frame_seed)
  ent <- apply(as.matrix(pooled[, seq_len(n_sites), drop = FALSE]), 2L,
               ViralEntropR::calculate_entropy)

  extra_args <- list(G          = config$MCLUST_G,
                     removez    = TRUE,
                     removesngl = TRUE,
                     verbose    = FALSE)
  if (!is.null(config$MCLUST_MODELS))
    extra_args$modelNames <- config$MCLUST_MODELS

  fit <- tryCatch(
    do.call(
      ViralEntropR::cluster_sites_by_entropy,
      c(list(entropies = ent, nr = n_seq, nsites = n_sites), extra_args)
    ),
    error = function(e)
      structure(list(error = conditionMessage(e)), class = "fit_error")
  )

  if (inherits(fit, "fit_error"))
    return(list(status = "gmm_error", class1_sites = integer(0L),
                G = NA_integer_, modelName = NA_character_,
                n_sequences = n_seq, n_class1_sites = 0L,
                reason = paste0("cluster_sites_by_entropy() failed: ",
                                fit$error)))

  df <- fit$DataFrame
  if (is.null(df) || nrow(df) == 0L)
    return(list(status = "no_sites_post_filter", class1_sites = integer(0L),
                G = NA_integer_, modelName = NA_character_,
                n_sequences = n_seq, n_class1_sites = 0L,
                reason = "GMM returned no sites after filtering."))

  if (any(df$class == 999L))
    return(list(status = "sentinel_999", class1_sites = integer(0L),
                G = NA_integer_, modelName = NA_character_,
                n_sequences = n_seq, n_class1_sites = 0L,
                reason = "GMM 999 sentinel: entropies all identical."))

  g_used <- if (is.list(fit$FitObject) && !is.null(fit$FitObject$G))
              as.integer(fit$FitObject$G) else NA_integer_
  mn     <- if (is.list(fit$FitObject) && !is.null(fit$FitObject$modelName))
              as.character(fit$FitObject$modelName) else NA_character_

  if (!is.na(g_used) && g_used == 1L)
    return(list(status = "G_eq_1", class1_sites = integer(0L),
                G = g_used, modelName = mn,
                n_sequences = n_seq, n_class1_sites = 0L,
                reason = "GMM collapsed to G=1: window non-informative."))

  df_rl <- ViralEntropR::relabel_entropy_classes(df)
  class1_sites <- sort(as.integer(df_rl$sites[df_rl$class == 1L]))

  if (length(class1_sites) < config$MIN_CLASS1_PER_FRAME)
    return(list(status = "too_few_class1", class1_sites = class1_sites,
                G = g_used, modelName = mn,
                n_sequences = n_seq, n_class1_sites = length(class1_sites),
                reason = sprintf("Class-1 has %d < %d sites.",
                                 length(class1_sites),
                                 config$MIN_CLASS1_PER_FRAME)))

  list(status = "ok", class1_sites = class1_sites,
       G = g_used, modelName = mn,
       n_sequences = n_seq, n_class1_sites = length(class1_sites),
       reason = NA_character_)
}


# ── 3. Per-frame fdahclust fit + per-site silhouette ────────────────────────
#
# Given the frame's class-1 curve_matrix (n_class1 × n_grid_points) and the
# grid vector, picks k by silhouette over N_CLUSTERS_GRID, runs the fit, and
# computes per-site silhouettes.
#
# Returns:
#   status              "ok" | "no_structure" | "fit_failed"
#   n_clusters          integer or NA
#   silhouette_mean     numeric or NA
#   silhouette_per_k    named numeric vector
#   cluster_assignments integer vector aligned with rows of curve_matrix
#   per_site_sil        numeric vector aligned with rows of curve_matrix
#   structurally_null   TRUE iff mean silhouette < MIN_MEAN_SILHOUETTE
#   oversized           logical vector of length n_clusters: cluster too big?

fit_frame_clustering <- function(curve_matrix, grid_vec, config) {

  n_curves <- nrow(curve_matrix)
  if (n_curves < min(config$N_CLUSTERS_GRID) + 1L)
    return(list(status = "fit_failed",
                n_clusters = NA_integer_,
                silhouette_mean = NA_real_,
                silhouette_per_k = setNames(numeric(0L), character(0L)),
                cluster_assignments = integer(0L),
                per_site_sil = numeric(0L),
                structurally_null = NA,
                oversized = logical(0L)))

  sel <- tryCatch(
    select_n_clusters(curve_matrix, grid_vec,
                       k_grid = config$N_CLUSTERS_GRID, config = config),
    error = function(e) NULL
  )
  if (is.null(sel))
    return(list(status = "fit_failed",
                n_clusters = NA_integer_,
                silhouette_mean = NA_real_,
                silhouette_per_k = setNames(numeric(0L), character(0L)),
                cluster_assignments = integer(0L),
                per_site_sil = numeric(0L),
                structurally_null = NA,
                oversized = logical(0L)))

  mem <- as.integer(sel$fit$memberships)
  sils_per_site <- per_site_silhouette(curve_matrix, mem)

  cluster_sizes <- tabulate(mem, nbins = sel$k)
  oversized <- cluster_sizes > config$MAX_CLUSTER_SIZE

  null_flag <- is.na(sel$silhouette_mean) ||
               sel$silhouette_mean < config$MIN_MEAN_SILHOUETTE

  list(status              = "ok",
       n_clusters          = as.integer(sel$k),
       silhouette_mean     = sel$silhouette_mean,
       silhouette_per_k    = sel$silhouette_per_k,
       cluster_assignments = mem,
       per_site_sil        = sils_per_site,
       structurally_null   = null_flag,
       oversized           = oversized,
       cluster_sizes       = as.integer(cluster_sizes))
}


# ── 4. Per-frame convergence flags ──────────────────────────────────────────
#
# For a frame's clustering, determines whether (a) all defining-SNP sites in
# class-1 share a single cluster, (b) same for Mutation_Sites, (c) same for
# all class-1 sites. Each test honours the three guards:
#   - Frame must not be structurally null (mean silhouette ≥ threshold).
#   - Sites with per-site silhouette < 0 are excluded ("marginal").
#   - Cluster must be ≤ MAX_CLUSTER_SIZE.
#   - At least one defining/mutation site must be present in class-1 for the
#     respective predicate to be evaluable (else returns FALSE: not yet).
#
# Returns three logical flags + the qualifying-cluster id (or NA if none).

evaluate_convergence <- function(class1_sites, cluster_assignments,
                                  per_site_sil, structurally_null, oversized,
                                  snp_sites, mutation_sites,
                                  config) {

  result <- list(
    converged_snp        = FALSE,
    converged_mut        = FALSE,
    converged_class1     = FALSE,
    qualifying_cluster_snp    = NA_integer_,
    qualifying_cluster_mut    = NA_integer_,
    qualifying_cluster_class1 = NA_integer_,
    n_snp_in_class1      = sum(class1_sites %in% snp_sites),
    n_mut_in_class1      = sum(class1_sites %in% mutation_sites)
  )

  if (isTRUE(structurally_null)) return(result)
  if (length(class1_sites) == 0L) return(result)
  if (length(cluster_assignments) != length(class1_sites)) return(result)

  # Identify well-clustered (sil >= 0) sites in usable (non-oversized) clusters.
  good_mask <- !is.na(per_site_sil) &
               per_site_sil >= 0 &
               !oversized[cluster_assignments]
  good_sites    <- class1_sites[good_mask]
  good_clusters <- cluster_assignments[good_mask]

  one_cluster <- function(target_sites) {
    target_in_good <- which(good_sites %in% target_sites)
    if (length(target_in_good) == 0L) return(c(FALSE, NA_integer_))
    cls <- good_clusters[target_in_good]
    if (length(unique(cls)) == 1L)
      return(c(TRUE, as.integer(cls[1L])))
    c(FALSE, NA_integer_)
  }

  snp_eval <- one_cluster(snp_sites)
  mut_eval <- one_cluster(mutation_sites)
  cls_eval <- one_cluster(good_sites)

  result$converged_snp    <- as.logical(snp_eval[1L])
  result$converged_mut    <- as.logical(mut_eval[1L])
  result$converged_class1 <- as.logical(cls_eval[1L])
  result$qualifying_cluster_snp    <- as.integer(snp_eval[2L])
  result$qualifying_cluster_mut    <- as.integer(mut_eval[2L])
  result$qualifying_cluster_class1 <- as.integer(cls_eval[2L])
  result
}


# ── 5. First-frame convergence summary across a cell ────────────────────────
#
# Given a list of per-frame convergence flags (in centre order), returns the
# first centre index at which each predicate is TRUE (NA if never).

first_convergence_frame <- function(per_frame_conv, frame_centres) {

  if (length(per_frame_conv) == 0L)
    return(list(first_snp = NA_integer_,
                first_mut = NA_integer_,
                first_class1 = NA_integer_))

  vec_snp <- vapply(per_frame_conv, function(x) isTRUE(x$converged_snp),
                     logical(1L))
  vec_mut <- vapply(per_frame_conv, function(x) isTRUE(x$converged_mut),
                     logical(1L))
  vec_cls <- vapply(per_frame_conv, function(x) isTRUE(x$converged_class1),
                     logical(1L))

  first_or_na <- function(v) {
    i <- which(v)[1L]
    if (length(i) == 0L || is.na(i)) NA_integer_ else as.integer(frame_centres[i])
  }

  list(first_snp    = first_or_na(vec_snp),
       first_mut    = first_or_na(vec_mut),
       first_class1 = first_or_na(vec_cls))
}


# ── 6. Per-frame entropy-curve matrix ───────────────────────────────────────
#
# Computes entropy of the frame's class-1 sites across ALL bins of the cell
# (not just the fitting window). The fitting window is used only for the GMM
# pre-selection; the curves are plotted across the full cell window so the
# reader can see what's happening outside the current fitting slice.

compute_frame_entropy_curves <- function(partition_result, class1_sites) {

  parts  <- partition_result$Partitions
  n_bins <- partition_result$n_bins
  if (length(class1_sites) == 0L)
    return(matrix(numeric(0L), nrow = 0L, ncol = n_bins))

  ent_mat <- matrix(0,
                    nrow = length(class1_sites), ncol = n_bins,
                    dimnames = list(as.character(class1_sites),
                                    names(parts)))
  for (k in seq_len(n_bins)) {
    chunk <- parts[[k]]
    if (nrow(chunk) == 0L) { ent_mat[, k] <- 0; next }
    sub <- chunk[, class1_sites, drop = FALSE]
    ent_mat[, k] <- apply(sub, 2L, ViralEntropR::calculate_entropy)
  }
  ent_mat
}


# ── 7. Per-frame Hellinger trajectories ─────────────────────────────────────
#
# Wraps compute_per_site_hellinger_anchor for a specific anchor index, using
# the per-frame class-1 set.

compute_frame_hellinger <- function(partition_result, class1_sites,
                                     anchor_idx, aa_levels = 25L,
                                     normalized = FALSE) {

  if (length(class1_sites) == 0L)
    return(matrix(numeric(0L), nrow = 0L,
                  ncol = max(0L, partition_result$n_bins - 1L)))

  compute_per_site_hellinger_anchor(partition_result, class1_sites,
                                     anchor_idx = anchor_idx,
                                     aa_levels  = aa_levels,
                                     normalized = normalized)
}


# ── 8. Frame-driver wrapper ─────────────────────────────────────────────────
#
# Single-frame orchestrator: GMM → curves → clustering → convergence. Returns
# the full per-frame record for one (metric × frame). Used by run_one_cell.R
# inside its per-frame loop.

run_one_frame <- function(partition_result, frame_bins, centre_bin,
                          metric, n_sites,
                          snp_sites, mutation_sites,
                          anchor_idx_T1, anchor_idx_Tpred,
                          frame_seed, config) {

  gmm <- fit_frame_gmm(partition_result, frame_bins, n_sites, config,
                       frame_seed = frame_seed)
  if (gmm$status != "ok") {
    return(list(
      centre_bin        = as.integer(centre_bin),
      frame_bins        = as.integer(frame_bins),
      metric            = metric,
      gmm_status        = gmm$status,
      gmm_G             = gmm$G,
      gmm_modelName     = gmm$modelName,
      gmm_n_sequences   = gmm$n_sequences,
      class1_sites      = integer(0L),
      n_class1_sites    = 0L,
      reason            = gmm$reason,
      curve_matrix      = matrix(numeric(0L), 0L, 0L),
      curve_matrix_transformed = matrix(numeric(0L), 0L, 0L),
      grid_vec          = integer(0L),
      n_clusters        = NA_integer_,
      silhouette_mean   = NA_real_,
      silhouette_per_k  = setNames(numeric(0L), character(0L)),
      cluster_assignments = integer(0L),
      per_site_sil      = numeric(0L),
      cluster_sizes     = integer(0L),
      oversized         = logical(0L),
      structurally_null = NA,
      converged_snp     = FALSE,
      converged_mut     = FALSE,
      converged_class1  = FALSE,
      qualifying_cluster_snp    = NA_integer_,
      qualifying_cluster_mut    = NA_integer_,
      qualifying_cluster_class1 = NA_integer_,
      n_snp_in_class1   = 0L,
      n_mut_in_class1   = 0L,
      anchor_idx        = NA_integer_,
      fit_status        = "skipped_no_class1"
    ))
  }

  class1_sites <- gmm$class1_sites

  # Build the curve matrix and grid for this metric.
  anchor_idx <- NA_integer_
  if (metric == "entropy") {
    curve <- compute_frame_entropy_curves(partition_result, class1_sites)
    grid_vec <- seq_len(ncol(curve))
  } else if (metric == "hellinger_T1") {
    anchor_idx <- anchor_idx_T1
    curve <- compute_frame_hellinger(partition_result, class1_sites,
                                       anchor_idx = anchor_idx_T1,
                                       aa_levels  = config$AA_LEVELS,
                                       normalized = config$HELLINGER_NORMALIZED)
    # grid: non-anchor bins by original order
    grid_vec <- setdiff(seq_len(partition_result$n_bins), anchor_idx_T1)
  } else if (metric == "hellinger_Tpred") {
    anchor_idx <- anchor_idx_Tpred
    if (is.na(anchor_idx_Tpred)) {
      return(list(
        centre_bin = as.integer(centre_bin),
        frame_bins = as.integer(frame_bins),
        metric     = metric,
        gmm_status = gmm$status, gmm_G = gmm$G, gmm_modelName = gmm$modelName,
        gmm_n_sequences = gmm$n_sequences,
        class1_sites = class1_sites, n_class1_sites = length(class1_sites),
        reason     = "No Tpred anchor available (no predecessor).",
        curve_matrix = matrix(numeric(0L), 0L, 0L),
        curve_matrix_transformed = matrix(numeric(0L), 0L, 0L),
        grid_vec = integer(0L),
        n_clusters = NA_integer_, silhouette_mean = NA_real_,
        silhouette_per_k = setNames(numeric(0L), character(0L)),
        cluster_assignments = integer(0L), per_site_sil = numeric(0L),
        cluster_sizes = integer(0L), oversized = logical(0L),
        structurally_null = NA,
        converged_snp = FALSE, converged_mut = FALSE, converged_class1 = FALSE,
        qualifying_cluster_snp = NA_integer_,
        qualifying_cluster_mut = NA_integer_,
        qualifying_cluster_class1 = NA_integer_,
        n_snp_in_class1 = 0L, n_mut_in_class1 = 0L,
        anchor_idx = NA_integer_,
        fit_status = "no_predecessor_anchor"
      ))
    }
    curve <- compute_frame_hellinger(partition_result, class1_sites,
                                       anchor_idx = anchor_idx_Tpred,
                                       aa_levels  = config$AA_LEVELS,
                                       normalized = config$HELLINGER_NORMALIZED)
    grid_vec <- setdiff(seq_len(partition_result$n_bins), anchor_idx_Tpred)
  } else {
    stop("Unknown metric: ", metric, call. = FALSE)
  }

  if (nrow(curve) == 0L || ncol(curve) < 2L)
    return(list(
      centre_bin = as.integer(centre_bin),
      frame_bins = as.integer(frame_bins),
      metric     = metric,
      gmm_status = gmm$status, gmm_G = gmm$G, gmm_modelName = gmm$modelName,
      gmm_n_sequences = gmm$n_sequences,
      class1_sites = class1_sites, n_class1_sites = length(class1_sites),
      reason     = "Curve matrix has zero rows or < 2 columns.",
      curve_matrix = curve,
      curve_matrix_transformed = curve,
      grid_vec   = grid_vec,
      n_clusters = NA_integer_, silhouette_mean = NA_real_,
      silhouette_per_k = setNames(numeric(0L), character(0L)),
      cluster_assignments = integer(0L), per_site_sil = numeric(0L),
      cluster_sizes = integer(0L), oversized = logical(0L),
      structurally_null = NA,
      converged_snp = FALSE, converged_mut = FALSE, converged_class1 = FALSE,
      qualifying_cluster_snp = NA_integer_,
      qualifying_cluster_mut = NA_integer_,
      qualifying_cluster_class1 = NA_integer_,
      n_snp_in_class1 = 0L, n_mut_in_class1 = 0L,
      anchor_idx = as.integer(anchor_idx),
      fit_status = "curve_too_small"
    ))

  curve_t <- if (isTRUE(config$LOG1P_TRANSFORM)) log1p_transform(curve) else curve

  cl <- fit_frame_clustering(curve_t, grid_vec, config)

  if (cl$status != "ok")
    return(list(
      centre_bin = as.integer(centre_bin),
      frame_bins = as.integer(frame_bins),
      metric     = metric,
      gmm_status = gmm$status, gmm_G = gmm$G, gmm_modelName = gmm$modelName,
      gmm_n_sequences = gmm$n_sequences,
      class1_sites = class1_sites, n_class1_sites = length(class1_sites),
      reason     = "fit_frame_clustering failed.",
      curve_matrix = curve,
      curve_matrix_transformed = curve_t,
      grid_vec   = grid_vec,
      n_clusters = NA_integer_, silhouette_mean = NA_real_,
      silhouette_per_k = setNames(numeric(0L), character(0L)),
      cluster_assignments = integer(0L), per_site_sil = numeric(0L),
      cluster_sizes = integer(0L), oversized = logical(0L),
      structurally_null = NA,
      converged_snp = FALSE, converged_mut = FALSE, converged_class1 = FALSE,
      qualifying_cluster_snp = NA_integer_,
      qualifying_cluster_mut = NA_integer_,
      qualifying_cluster_class1 = NA_integer_,
      n_snp_in_class1 = 0L, n_mut_in_class1 = 0L,
      anchor_idx = as.integer(anchor_idx),
      fit_status = "clustering_failed"
    ))

  conv <- evaluate_convergence(
    class1_sites        = class1_sites,
    cluster_assignments = cl$cluster_assignments,
    per_site_sil        = cl$per_site_sil,
    structurally_null   = cl$structurally_null,
    oversized           = cl$oversized,
    snp_sites           = snp_sites,
    mutation_sites      = mutation_sites,
    config              = config
  )

  list(
    centre_bin          = as.integer(centre_bin),
    frame_bins          = as.integer(frame_bins),
    metric              = metric,
    gmm_status          = gmm$status,
    gmm_G               = gmm$G,
    gmm_modelName       = gmm$modelName,
    gmm_n_sequences     = gmm$n_sequences,
    class1_sites        = class1_sites,
    n_class1_sites      = length(class1_sites),
    reason              = NA_character_,
    curve_matrix        = curve,
    curve_matrix_transformed = curve_t,
    grid_vec            = grid_vec,
    n_clusters          = cl$n_clusters,
    silhouette_mean     = cl$silhouette_mean,
    silhouette_per_k    = cl$silhouette_per_k,
    cluster_assignments = cl$cluster_assignments,
    per_site_sil        = cl$per_site_sil,
    cluster_sizes       = cl$cluster_sizes,
    oversized           = cl$oversized,
    structurally_null   = cl$structurally_null,
    converged_snp       = conv$converged_snp,
    converged_mut       = conv$converged_mut,
    converged_class1    = conv$converged_class1,
    qualifying_cluster_snp    = conv$qualifying_cluster_snp,
    qualifying_cluster_mut    = conv$qualifying_cluster_mut,
    qualifying_cluster_class1 = conv$qualifying_cluster_class1,
    n_snp_in_class1     = conv$n_snp_in_class1,
    n_mut_in_class1     = conv$n_mut_in_class1,
    anchor_idx          = as.integer(anchor_idx),
    fit_status          = "ok"
  )
}
