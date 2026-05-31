################################################################################
## run_one_cell.R — FDA_Analysis (v2)
##
## Per-cell driver. A "cell" is a triple (dataset, variant, strategy).
##
## v2 changes from v1:
##   • Three metrics: entropy, hellinger_T1, hellinger_Tpred.
##   • The Hellinger anchors are: T1 = window-start partition; Tpred =
##     partition containing the predecessor variant's US detection (resolved
##     via config$PREDECESSOR_MAP). Tpred is clipped to the cell window if
##     the predecessor date is outside.
##   • Full (1273-site) entropy and Hellinger matrices are computed and
##     stored — needed for the "background curves" in the new dual-panel
##     visualization.
##   • A per-frame loop runs after the cell-level fits. Frames are centred at
##     bins (h+1)..(n_bins-h) with h = FITTING_WINDOW_HALFWIDTH. Each frame:
##         - Re-fits a GMM on the pooled sequences of the (centre ± h)
##           partitions of the cell's window.
##         - Computes per-frame entropy / Hellinger curves on the frame
##           class-1 set.
##         - Runs fdahclust at the silhouette-best k over the grid.
##         - Records cluster assignments, per-site silhouettes, oversized
##           flags, convergence flags (SNP / Mutation / class-1).
##   • Bootstrap stability and SNP enrichment continue to be computed once
##     per cell on the FULL-WINDOW fit, as in v1 (Option γ design).
##   • New convergence summary: `convergence_first` records the first frame
##     centre index at which each predicate becomes TRUE.
##
## Workflow (v2):
##   1. Resolve IDs and seeds; resolve predecessor.
##   2. Build the cell window; emit stub on out_of_coverage.
##   3. Slice the feature matrix.
##   4. Fit the cell-level GMM. Stub on reduced_skipped / too_few_sites.
##   5. Partition; stub on underpowered_bin.
##   6. Build the FULL 1273-site entropy matrix and FULL 1273-site Hellinger
##      matrices (T1 and Tpred). These are needed for the background panel.
##   7. Build class-1 entropy / Hellinger matrices. log1p-transform.
##   8. Cell-level fits + bootstrap + enrichment for each of the 3 metrics.
##   9. Per-frame loop: build the per-frame records (one list per metric).
##  10. Compute the first-convergence summary per metric.
##  11. Assemble the RDS payload and write atomically.
##
## Resilience: every step that can fail is wrapped in tryCatch; non-fatal
## failures populate stub fields in the payload and the cell-level status is
## downgraded but the file is always written.
################################################################################

# Pre-requisite: setup.R + helpers_fda.R + helpers_frames.R have been sourced.


run_one_cell <- function(dataset_name, variant_name, strategy_name,
                          feature_matrix, truth_df, config) {

  t_start  <- Sys.time()
  out_path <- cell_path(config, dataset_name, variant_name, strategy_name)

  # ── (0) IDs and seeds ─────────────────────────────────────────────────────
  ds_id  <- dataset_id_for (dataset_name,  config)
  var_id <- variant_id_for (variant_name,  config)
  str_id <- strategy_id_for(strategy_name, config)

  cell_seed     <- seed_for_cell(config$BASE_SEED, ds_id, var_id)
  strategy_seed <- seed_for_strategy(cell_seed, str_id)

  base_meta <- list(
    schema_version  = config$SCHEMA_VERSION %||% 2L,
    dataset         = dataset_name,
    variant         = variant_name,
    strategy        = strategy_name,
    dataset_id      = ds_id,
    variant_id      = var_id,
    strategy_id     = str_id,
    seed_cell       = cell_seed,
    seed_strategy   = strategy_seed,
    run_timestamp   = format(t_start, "%Y-%m-%dT%H:%M:%S%z"),
    r_version       = paste(R.version$major, R.version$minor, sep = "."),
    package_version = as.character(utils::packageVersion("ViralEntropR"))
  )

  emit_stub <- function(status, reason, extra = list()) {
    payload <- c(
      base_meta,
      list(
        status        = status,
        status_reason = as.character(reason),
        walltime_s    = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      ),
      extra
    )
    save_rds_atomic(payload, out_path)
    invisible(payload)
  }

  # ── (1) Truth lookup + predecessor resolution ────────────────────────────
  truth_row <- truth_df[truth_df$WHO_Label == variant_name, , drop = FALSE]
  if (nrow(truth_row) == 0L)
    return(emit_stub("error",
                      sprintf("Variant %s not in truth catalogue.",
                              variant_name)))
  if (nrow(truth_row) > 1L) truth_row <- truth_row[1L, , drop = FALSE]

  detection_date <- truth_row$Date_First_Detected_US[1L]
  snp_sites      <- truth_row$Defining_SNP_Sites[[1L]]
  mutation_sites <- truth_row$Mutation_Sites[[1L]]
  if (is.na(detection_date))
    return(emit_stub("error",
                      sprintf("No US detection date for %s.", variant_name)))

  pred <- resolve_predecessor(variant_name, truth_df, config)
  base_meta$predecessor_name     <- pred$predecessor_name
  base_meta$predecessor_us_date  <- pred$predecessor_us_date
  base_meta$predecessor_status   <- pred$status

  # ── (2) Window construction ───────────────────────────────────────────────
  fm_start <- min(feature_matrix$Date)
  fm_end   <- max(feature_matrix$Date)
  win      <- build_cell_window(detection_date, fm_start, fm_end,
                                 K_pre  = config$K_BINS_BEFORE,
                                 K_post = config$K_BINS_AFTER,
                                 ref_bin_months = config$REF_BIN_MONTHS)

  if (win$status != "ok")
    return(emit_stub(win$status, win$reason,
                      extra = list(detection_date = detection_date,
                                    fm_start       = fm_start,
                                    fm_end         = fm_end)))

  # ── (3) Window slicing ────────────────────────────────────────────────────
  fm_slice <- tryCatch(
    slice_window(feature_matrix, win$window_start, win$window_end,
                  n_sites = config$N_SITES),
    error = function(e) NULL
  )
  if (is.null(fm_slice) || nrow(fm_slice) == 0L)
    return(emit_stub("error", "Window slice produced zero sequences.",
                      extra = list(detection_date = detection_date,
                                    window_start   = win$window_start,
                                    window_end     = win$window_end)))

  # ── (4) Cell-level GMM ────────────────────────────────────────────────────
  set.seed(cell_seed)
  gmm <- tryCatch(
    fit_cell_gmm(fm_slice, n_sites = config$N_SITES, config = config),
    error = function(e)
      list(reduced_skipped = TRUE, class1_sites = integer(0L),
           reason = paste0("fit_cell_gmm error: ", conditionMessage(e)),
           G = NA_integer_, modelName = NA_character_, n_class1_sites = 0L)
  )

  gmm_meta <- list(
    G              = gmm$G,
    modelName      = gmm$modelName,
    n_class1_sites = gmm$n_class1_sites,
    class1_sites   = gmm$class1_sites,
    reduced_skipped = isTRUE(gmm$reduced_skipped),
    reason         = gmm$reason
  )

  if (gmm$reduced_skipped || length(gmm$class1_sites) == 0L)
    return(emit_stub("reduced_skipped", gmm$reason,
                      extra = list(detection_date = detection_date,
                                    window_start   = win$window_start,
                                    window_end     = win$window_end,
                                    K_pre_realised = win$K_pre_realised,
                                    K_post_realised = win$K_post_realised,
                                    n_seqs_window  = nrow(fm_slice),
                                    gmm_meta       = gmm_meta)))

  if (length(gmm$class1_sites) < config$MIN_CLASS1_SITES_FOR_FDA)
    return(emit_stub("too_few_sites",
                      sprintf("Class-1 has %d < %d sites.",
                              length(gmm$class1_sites),
                              config$MIN_CLASS1_SITES_FOR_FDA),
                      extra = list(detection_date = detection_date,
                                    window_start   = win$window_start,
                                    window_end     = win$window_end,
                                    K_pre_realised = win$K_pre_realised,
                                    K_post_realised = win$K_post_realised,
                                    n_seqs_window  = nrow(fm_slice),
                                    gmm_meta       = gmm_meta)))

  # ── (5) Strategy-specific partitioning ────────────────────────────────────
  strategy_def <- config$STRATEGY_DEFS[[strategy_name]]
  part <- tryCatch(
    fda_partition(fm_slice, win$window_start, win$window_end, strategy_def,
                   n_sites = config$N_SITES),
    error = function(e) NULL
  )
  if (is.null(part) || part$n_bins < 2L)
    return(emit_stub("error",
                      sprintf("fda_partition produced %s bins.",
                              if (is.null(part)) "NA" else part$n_bins),
                      extra = list(detection_date = detection_date,
                                    window_start   = win$window_start,
                                    window_end     = win$window_end,
                                    gmm_meta       = gmm_meta)))

  underpowered <- which(part$bin_counts < config$MIN_SEQUENCES_PER_BIN)
  if (length(underpowered) > 0L)
    return(emit_stub("underpowered_bin",
                      sprintf("%d / %d bins have < %d seqs (worst = %d).",
                              length(underpowered), part$n_bins,
                              config$MIN_SEQUENCES_PER_BIN,
                              min(part$bin_counts)),
                      extra = list(detection_date = detection_date,
                                    window_start   = win$window_start,
                                    window_end     = win$window_end,
                                    K_pre_realised = win$K_pre_realised,
                                    K_post_realised = win$K_post_realised,
                                    n_seqs_window  = nrow(fm_slice),
                                    gmm_meta       = gmm_meta,
                                    bin_starts     = part$bin_starts,
                                    bin_ends       = part$bin_ends,
                                    bin_counts     = part$bin_counts,
                                    bin_labels     = part$Dates_Labels,
                                    n_bins         = part$n_bins)))

  # ── (6) Hellinger anchor indices ─────────────────────────────────────────
  anchor_T1   <- list(bin_idx = 1L, snapped = FALSE, distance_months = 0)
  anchor_Tpred <- if (!is.na(pred$predecessor_us_date))
                    anchor_bin_idx(pred$predecessor_us_date,
                                    part$bin_starts, part$bin_ends)
                  else list(bin_idx = NA_integer_, snapped = NA,
                            distance_months = NA_real_)

  # ── (7) FULL 1273-site matrices (for background panel) ───────────────────
  all_sites <- seq_len(config$N_SITES)

  full_entropy <- tryCatch(
    compute_entropy_matrix(part, all_sites),
    error = function(e) NULL
  )
  full_hellinger_T1 <- tryCatch(
    compute_per_site_hellinger_anchor(part, all_sites,
                                       anchor_idx = anchor_T1$bin_idx,
                                       aa_levels  = config$AA_LEVELS,
                                       normalized = config$HELLINGER_NORMALIZED),
    error = function(e) NULL
  )
  full_hellinger_Tpred <- if (!is.na(anchor_Tpred$bin_idx))
    tryCatch(
      compute_per_site_hellinger_anchor(part, all_sites,
                                         anchor_idx = anchor_Tpred$bin_idx,
                                         aa_levels  = config$AA_LEVELS,
                                         normalized = config$HELLINGER_NORMALIZED),
      error = function(e) NULL
    )
  else NULL

  if (is.null(full_entropy))
    return(emit_stub("error",
                      "compute_entropy_matrix on all sites failed.",
                      extra = list(gmm_meta = gmm_meta)))
  if (is.null(full_hellinger_T1))
    return(emit_stub("error",
                      "compute_per_site_hellinger on T1 anchor failed.",
                      extra = list(gmm_meta = gmm_meta)))

  # ── (8) Class-1 matrices ─────────────────────────────────────────────────
  class1 <- gmm$class1_sites
  ent_class1   <- full_entropy[as.character(class1), , drop = FALSE]
  hel_T1_class1   <- full_hellinger_T1[as.character(class1), , drop = FALSE]
  hel_Tpred_class1 <- if (!is.null(full_hellinger_Tpred))
                       full_hellinger_Tpred[as.character(class1), , drop = FALSE]
                     else NULL

  log1p_t <- function(M) if (isTRUE(config$LOG1P_TRANSFORM)) log1p_transform(M) else M
  ent_class1_t      <- log1p_t(ent_class1)
  hel_T1_class1_t   <- log1p_t(hel_T1_class1)
  hel_Tpred_class1_t <- if (!is.null(hel_Tpred_class1)) log1p_t(hel_Tpred_class1)
                      else NULL

  transformation <- if (isTRUE(config$LOG1P_TRANSFORM)) "log1p" else "identity"

  # ── (9) Cell-level fits + stability + enrichment per metric ──────────────
  detection_bin <- detection_bin_index(detection_date, part$bin_starts,
                                        part$bin_ends)

  fit_one_metric <- function(M_t, grid_vec, label) {
    if (is.null(M_t) || nrow(M_t) == 0L || ncol(M_t) < 2L)
      return(list(status = "skipped", status_reason = "Empty / too-small matrix.",
                  fit = NULL, n_clusters = NA_integer_,
                  silhouette_mean = NA_real_,
                  silhouette_per_k = setNames(numeric(0L), character(0L)),
                  cluster_assignments = integer(0L),
                  transformation = transformation))
    sel <- tryCatch(
      select_n_clusters(M_t, grid_vec, k_grid = config$N_CLUSTERS_GRID,
                         config = config),
      error = function(e) NULL
    )
    if (is.null(sel))
      return(list(status = "error",
                  status_reason = paste0("select_n_clusters failed for ",
                                          label, "."),
                  fit = NULL, n_clusters = NA_integer_,
                  silhouette_mean = NA_real_,
                  silhouette_per_k = setNames(numeric(0L), character(0L)),
                  cluster_assignments = integer(0L),
                  transformation = transformation))

    list(status              = "ok",
         status_reason       = NA_character_,
         fit                 = sel$fit,
         n_clusters          = sel$k,
         silhouette_mean     = sel$silhouette_mean,
         silhouette_per_k    = sel$silhouette_per_k,
         cluster_assignments = as.integer(sel$fit$memberships),
         transformation      = transformation)
  }

  ent_grid    <- seq_len(ncol(ent_class1_t))
  hel_T1_grid <- as.integer(sub("^T", "", colnames(hel_T1_class1_t)))
  hel_Tpred_grid <- if (!is.null(hel_Tpred_class1_t))
                      as.integer(sub("^T", "", colnames(hel_Tpred_class1_t)))
                    else integer(0L)

  fda_entropy        <- fit_one_metric(ent_class1_t,     ent_grid,    "entropy")
  fda_hellinger_T1   <- fit_one_metric(hel_T1_class1_t,  hel_T1_grid, "hellinger_T1")
  fda_hellinger_Tpred <- fit_one_metric(hel_Tpred_class1_t, hel_Tpred_grid,
                                          "hellinger_Tpred")

  stability_for <- function(M_t, grid_vec, fit_block, seed_offset) {
    if (fit_block$status != "ok") return(NULL)
    tryCatch(
      bootstrap_jaccard_stability(M_t, grid_vec,
                                   n_clusters = fit_block$n_clusters,
                                   B          = config$N_BOOTSTRAP,
                                   seed       = strategy_seed + seed_offset,
                                   config     = config),
      error = function(e)
        list(n_bootstrap = config$N_BOOTSTRAP, n_failed_bootstrap = NA_integer_,
             per_cluster_jaccard = rep(NA_real_, fit_block$n_clusters),
             n_stable_clusters = NA_integer_, error = conditionMessage(e))
    )
  }

  enrichment_for <- function(fit_block) {
    if (fit_block$status != "ok") return(NULL)
    tryCatch(
      test_snp_enrichment(class1_sites = class1,
                           cluster_assignments = fit_block$cluster_assignments,
                           snp_sites = snp_sites,
                           n_total_tests = config$N_TOTAL_ENRICHMENT_TESTS),
      error = function(e)
        list(fisher_p = NA_real_, fisher_p_bonferroni = NA_real_,
             n_snp_in_class1 = NA_integer_, n_class1_total = NA_integer_,
             n_snp_total_truth = length(snp_sites),
             max_overlap_cluster = NA_integer_, max_overlap_OR = NA_real_,
             snp_overlap_per_cluster = setNames(integer(0L), character(0L)),
             contingency_table = NULL, status = "error",
             error = conditionMessage(e))
    )
  }

  stability_entropy        <- stability_for(ent_class1_t, ent_grid,
                                              fda_entropy, 1L)
  stability_hellinger_T1   <- stability_for(hel_T1_class1_t, hel_T1_grid,
                                              fda_hellinger_T1, 2L)
  stability_hellinger_Tpred <- stability_for(hel_Tpred_class1_t, hel_Tpred_grid,
                                               fda_hellinger_Tpred, 3L)

  snp_entropy        <- enrichment_for(fda_entropy)
  snp_hellinger_T1   <- enrichment_for(fda_hellinger_T1)
  snp_hellinger_Tpred <- enrichment_for(fda_hellinger_Tpred)

  # ── (10) Per-frame loop ───────────────────────────────────────────────────
  halfwidth     <- config$FITTING_WINDOW_HALFWIDTH
  frame_centres <- frame_centre_indices(part$n_bins, halfwidth)

  per_frame_entropy        <- list()
  per_frame_hellinger_T1   <- list()
  per_frame_hellinger_Tpred <- list()

  for (i in seq_along(frame_centres)) {
    c_bin <- frame_centres[i]
    fb    <- seq.int(c_bin - halfwidth, c_bin + halfwidth)
    seed_frame <- strategy_seed + 100L + as.integer(c_bin)

    per_frame_entropy[[i]] <- run_one_frame(
      partition_result = part,  frame_bins = fb,  centre_bin = c_bin,
      metric           = "entropy",
      n_sites          = config$N_SITES,
      snp_sites        = snp_sites,
      mutation_sites   = mutation_sites,
      anchor_idx_T1    = anchor_T1$bin_idx,
      anchor_idx_Tpred = anchor_Tpred$bin_idx,
      frame_seed       = seed_frame,
      config           = config
    )
    per_frame_hellinger_T1[[i]] <- run_one_frame(
      partition_result = part,  frame_bins = fb,  centre_bin = c_bin,
      metric           = "hellinger_T1",
      n_sites          = config$N_SITES,
      snp_sites        = snp_sites,
      mutation_sites   = mutation_sites,
      anchor_idx_T1    = anchor_T1$bin_idx,
      anchor_idx_Tpred = anchor_Tpred$bin_idx,
      frame_seed       = seed_frame + 1L,
      config           = config
    )
    per_frame_hellinger_Tpred[[i]] <- run_one_frame(
      partition_result = part,  frame_bins = fb,  centre_bin = c_bin,
      metric           = "hellinger_Tpred",
      n_sites          = config$N_SITES,
      snp_sites        = snp_sites,
      mutation_sites   = mutation_sites,
      anchor_idx_T1    = anchor_T1$bin_idx,
      anchor_idx_Tpred = anchor_Tpred$bin_idx,
      frame_seed       = seed_frame + 2L,
      config           = config
    )
  }

  # ── (11) Convergence summaries (first frame at which each predicate holds)
  conv_entropy        <- first_convergence_frame(per_frame_entropy,        frame_centres)
  conv_hellinger_T1   <- first_convergence_frame(per_frame_hellinger_T1,   frame_centres)
  conv_hellinger_Tpred <- first_convergence_frame(per_frame_hellinger_Tpred, frame_centres)

  # Overall cell-level status: ok iff at least one metric succeeded.
  metric_statuses <- c(fda_entropy$status,
                        fda_hellinger_T1$status,
                        fda_hellinger_Tpred$status)
  overall_status <- if (any(metric_statuses == "ok")) "ok" else "error"
  overall_reason <- if (overall_status == "ok") NA_character_ else
                    paste(c(fda_entropy$status_reason,
                            fda_hellinger_T1$status_reason,
                            fda_hellinger_Tpred$status_reason),
                          collapse = " | ")

  # ── (12) Assemble RDS payload ─────────────────────────────────────────────
  payload <- c(
    base_meta,
    list(
      detection_date    = detection_date,
      window_start      = win$window_start,
      window_end        = win$window_end,
      K_pre_realised    = win$K_pre_realised,
      K_post_realised   = win$K_post_realised,
      n_seqs_window     = nrow(fm_slice),

      bin_starts        = part$bin_starts,
      bin_ends          = part$bin_ends,
      bin_labels        = part$Dates_Labels,
      bin_counts        = part$bin_counts,
      n_bins            = part$n_bins,
      detection_bin     = detection_bin,

      # Anchor metadata
      anchor_T1         = anchor_T1,
      anchor_Tpred      = anchor_Tpred,

      gmm_meta              = gmm_meta,
      snp_sites_truth       = snp_sites,
      mutation_sites_truth  = mutation_sites,

      # Full 1273-site matrices for background curves
      entropy_full_matrix         = full_entropy,
      hellinger_T1_full_matrix    = full_hellinger_T1,
      hellinger_Tpred_full_matrix = full_hellinger_Tpred,

      # Class-1 matrices (raw and transformed)
      entropy_class1_matrix         = ent_class1,
      hellinger_T1_class1_matrix    = hel_T1_class1,
      hellinger_Tpred_class1_matrix = hel_Tpred_class1,
      entropy_class1_matrix_transformed         = ent_class1_t,
      hellinger_T1_class1_matrix_transformed    = hel_T1_class1_t,
      hellinger_Tpred_class1_matrix_transformed = hel_Tpred_class1_t,
      transformation    = transformation,

      # Cell-level fits / stability / enrichment per metric
      fda_entropy             = fda_entropy,
      fda_hellinger_T1        = fda_hellinger_T1,
      fda_hellinger_Tpred     = fda_hellinger_Tpred,
      stability_entropy        = stability_entropy,
      stability_hellinger_T1   = stability_hellinger_T1,
      stability_hellinger_Tpred = stability_hellinger_Tpred,
      snp_enrichment_entropy        = snp_entropy,
      snp_enrichment_hellinger_T1   = snp_hellinger_T1,
      snp_enrichment_hellinger_Tpred = snp_hellinger_Tpred,

      # Per-frame records
      frame_centres                  = frame_centres,
      frames_entropy                 = per_frame_entropy,
      frames_hellinger_T1            = per_frame_hellinger_T1,
      frames_hellinger_Tpred         = per_frame_hellinger_Tpred,

      # Convergence summaries
      convergence_entropy            = conv_entropy,
      convergence_hellinger_T1       = conv_hellinger_T1,
      convergence_hellinger_Tpred    = conv_hellinger_Tpred,

      status         = overall_status,
      status_reason  = overall_reason,
      walltime_s     = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    )
  )

  save_rds_atomic(payload, out_path)
  invisible(payload)
}
