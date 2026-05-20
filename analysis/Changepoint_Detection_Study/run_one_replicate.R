################################################################################
## run_one_replicate.R — Changepoint_Detection_Study
##
## Per-window driver. Inputs:
##   - cell_id, dataset name (resolves the row of the manifest)
##   - the manifest (loaded once at subprocess startup)
##   - the feature matrix (loaded once at subprocess startup)
##   - the config
##
## Per replicate the driver:
##   1. Slices sequences in the window from the feature matrix.
##   2. Partitions into 2-month bins.
##   3. Computes full-1273-site Hellinger trajectory and runs all four
##      methods (plus the K-sweep for ks.cp3o) on it.
##   4. Fits a per-window GMM on the window's sequences; extracts class-1
##      sites; computes the reduced-site Hellinger trajectory and runs all
##      four methods on it (unless the GMM degenerates to the 999 sentinel,
##      in which case the reduced-site analysis is recorded as skipped).
##   5. Computes per-method metrics (P, R, F1, TLE) at the primary tolerance
##      and across the truth-shift sweep.
##   6. Saves a per-replicate RDS via save_rds_atomic.
##
## Schema of the saved RDS is documented in the README §"Per-replicate output
## schema".
################################################################################

# Required upstream sources: setup.R, helpers_windows.R, helpers_hellinger.R,
# cp_methods.R, metrics.R.

run_one_replicate <- function(dataset_name, cell_id, manifest, feature_matrix,
                              config, output_base_dir = NULL) {

  t_start_total <- Sys.time()
  require_runtime_packages()

  # ── 1. Resolve manifest row ───────────────────────────────────────────────
  row_idx <- which(manifest$cell_id == cell_id)
  if (length(row_idx) != 1L)
    stop("cell_id ", cell_id, " not found in manifest (or duplicated).",
         call. = FALSE)
  row <- manifest[row_idx, , drop = FALSE]

  dataset_id  <- as.integer(row$dataset_id)
  tier_id     <- as.integer(row$tier_id)
  run_in_tier <- as.integer(row$run_in_tier)
  seed_rep    <- as.integer(row$seed_replicate)

  # Re-derive seed defensively and assert agreement (catches manifest drift).
  seed_check <- seed_for_replicate(config$BASE_SEED, dataset_id, tier_id,
                                   run_in_tier)
  if (seed_rep != seed_check)
    stop(sprintf(
      "Seed drift detected: manifest stored %d, recomputed %d.",
      seed_rep, seed_check), call. = FALSE)
  set.seed(seed_rep)

  window_start <- row$window_start[[1L]]
  window_end   <- row$window_end[[1L]]
  n_bins       <- as.integer(row$n_bins)

  truth_hellinger_idx <- row$truth_hellinger_idx[[1L]]
  if (length(truth_hellinger_idx) == 0L)
    stop("Manifest row has empty truth_hellinger_idx; window-acceptance bug.",
         call. = FALSE)

  # ── 2. Slice window sequences ─────────────────────────────────────────────
  site_cols <- setdiff(colnames(feature_matrix), c("Date", "Country"))
  in_window <- feature_matrix$Date >= window_start &
               feature_matrix$Date <  window_end
  win_mat <- feature_matrix[in_window, , drop = FALSE]
  if (nrow(win_mat) == 0L)
    stop("Empty window slice; window-acceptance contract violated.",
         call. = FALSE)

  # ── 3. Partition into 2-month bins ────────────────────────────────────────
  partitions <- partition_window(
    feature_matrix = win_mat,
    window_start   = window_start,
    window_end     = window_end,
    site_cols      = site_cols,
    bin_months     = config$BIN_MONTHS
  )

  # ── 4. Full-site Hellinger + four methods ─────────────────────────────────
  hellinger_full <- compute_hellinger_full(
    partitions  = partitions,
    sites       = seq_len(length(site_cols)),
    aa_levels   = config$AA_LEVELS,
    normalized  = config$HELLINGER_NORMALIZED
  )
  dat_t_full   <- t(hellinger_full)
  n_hell_rows  <- nrow(dat_t_full)

  methods_full <- run_all_methods(dat_t_full, config)

  metrics_full <- evaluate_methods(
    method_results    = methods_full,
    truth_idx         = truth_hellinger_idx,
    n_hellinger_rows  = n_hell_rows,
    tolerance         = 1L,
    truth_shifts      = if (isTRUE(config$RUN_TRUTH_SHIFT_SWEEP))
                          config$TRUTH_SHIFT_BINS else NULL
  )

  # ── 5. Per-window GMM + reduced-site Hellinger + four methods ─────────────
  win_seq_matrix <- as.matrix(win_mat[, site_cols, drop = FALSE])
  gmm_meta <- fit_window_gmm(win_seq_matrix, config)

  if (isTRUE(gmm_meta$reduced_skipped) ||
      length(gmm_meta$class1_sites) < 2L) {
    # < 2 class-1 sites is also degenerate: a Hellinger trajectory on a
    # single column is trivially constant — no CP can be detected.
    methods_reduced <- NULL
    metrics_reduced <- NULL
    reduced_status  <- "skipped"
    reduced_reason  <- if (isTRUE(gmm_meta$reduced_skipped))
                         "GMM returned 999 sentinel or zero class-1 sites"
                       else
                         sprintf("only %d class-1 sites (< 2)",
                                 length(gmm_meta$class1_sites))
  } else {
    hellinger_red <- compute_hellinger_full(
      partitions  = partitions,
      sites       = gmm_meta$class1_sites,
      aa_levels   = config$AA_LEVELS,
      normalized  = config$HELLINGER_NORMALIZED
    )
    dat_t_red   <- t(hellinger_red)
    methods_reduced <- run_all_methods(dat_t_red, config)
    metrics_reduced <- evaluate_methods(
      method_results    = methods_reduced,
      truth_idx         = truth_hellinger_idx,
      n_hellinger_rows  = nrow(dat_t_red),
      tolerance         = 1L,
      truth_shifts      = if (isTRUE(config$RUN_TRUTH_SHIFT_SWEEP))
                            config$TRUTH_SHIFT_BINS else NULL
    )
    reduced_status <- "ok"
    reduced_reason <- NA_character_
  }

  # ── 6. Assemble result list ──────────────────────────────────────────────
  result <- list(
    # Identifying metadata
    dataset             = dataset_name,
    dataset_id          = dataset_id,
    cell_id             = as.integer(cell_id),
    tier                = as.character(row$tier),
    tier_id             = tier_id,
    run_in_tier         = run_in_tier,
    seed                = seed_rep,

    # Window metadata
    window_start        = window_start,
    window_end          = window_end,
    length_months       = as.integer(row$length_months),
    n_bins              = n_bins,
    n_seqs              = nrow(win_mat),
    bin_counts          = row$bin_counts[[1L]],

    # Truth metadata
    truths_raw          = as.integer(row$truths_raw),
    truths_dropped_T1   = as.integer(row$truths_dropped_T1),
    truths_effective    = as.integer(row$truths_effective),
    truth_dates         = row$truth_dates[[1L]],
    truth_labels        = row$truth_labels[[1L]],
    truth_hellinger_idx = truth_hellinger_idx,

    # GMM diagnostics
    gmm_meta = list(
      G                = gmm_meta$G,
      modelName        = gmm_meta$modelName,
      n_class1_sites   = gmm_meta$n_class1_sites,
      class1_sites     = gmm_meta$class1_sites,
      reduced_skipped  = gmm_meta$reduced_skipped
    ),
    reduced_status      = reduced_status,
    reduced_reason      = reduced_reason,

    # Method results — full site set
    # Note: metrics_full contains each method's detection set and metadata
    # PLUS the per-method metrics_primary and metrics_truth_shift fields,
    # so we do not separately store methods_full (would be a strict subset).
    metrics_full        = metrics_full,

    # Method results — reduced site set (NULL when skipped)
    metrics_reduced     = metrics_reduced,

    # Versioning + walltime
    walltime_s          = as.numeric(difftime(Sys.time(), t_start_total,
                                              units = "secs")),
    r_version           = R.version.string,
    package_version     = as.character(utils::packageVersion("ViralEntropR"))
  )

  # ── 7. Atomic write ──────────────────────────────────────────────────────
  base_dir <- if (is.null(output_base_dir)) output_dir(config) else output_base_dir
  out_path <- replicate_path(config, dataset_name, cell_id, base_dir = base_dir)
  save_rds_atomic(result, out_path)

  invisible(out_path)
}
