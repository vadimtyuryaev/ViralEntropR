################################################################################
## run_full_dataset_detection.R — Changepoint_Detection_Study
##
## For each dataset (NCBI_US, GISAID_US):
##   1. Load the feature matrix.
##   2. Partition the ENTIRE dataset's date range into 2-month non-overlapping
##      bins (no sub-windowing).
##   3. Compute the full-1,273-site Hellinger trajectory T2..T_n vs T1.
##   4. Run all 9 method configurations on the trajectory.
##   5. Save the harmonised detection sets, the trajectory, and the bin grid
##      to outputs/full_dataset_detection_<dataset>.rds.
##
## This is a once-off computation (not part of the 5,001-window benchmark).
## Used by plot_results.R for the per-method timeline figures.
##
## Usage:
##   Rscript run_full_dataset_detection.R
##
## Env overrides (same as simulation_study.R):
##   VIRAL_CP_STUDY_DIR=/path
##   FEATURE_RDS_NCBI_US=/path
##   FEATURE_RDS_GISAID_US=/path
################################################################################

suppressPackageStartupMessages({
  source("setup.R")
  source("helpers_windows.R")
  source("helpers_hellinger.R")
  source("cp_methods.R")
  source("metrics.R")
})

config <- build_config()
dir.create(output_dir(config), recursive = TRUE, showWarnings = FALSE)

log_msg("==============================================================")
log_msg("Full-dataset detection — single run per dataset, per method")
log_msg("==============================================================")

truth_df <- load_truth_catalogue()
log_msg(sprintf("Truth catalogue: %d VOC/VOI variants with US dates.",
                nrow(truth_df)))

# ── Per-dataset driver ───────────────────────────────────────────────────────

run_one_full_dataset <- function(dataset_name, config, truth_df) {

  out_path <- file.path(
    output_dir(config),
    sprintf("full_dataset_detection_%s.rds", dataset_name)
  )

  if (file.exists(out_path)) {
    log_msg(sprintf("[%s] Existing detection file found: %s. Skipping.",
                    dataset_name, out_path))
    return(invisible(readRDS(out_path)))
  }

  require_runtime_packages()

  fm <- load_feature_matrix(config$DATASETS[[dataset_name]]$feature_rds)
  data_start <- min(fm$Date)
  data_end   <- max(fm$Date)
  log_msg(sprintf(
    "[%s] Feature matrix: %d sequences, dates %s → %s.",
    dataset_name, nrow(fm), format(data_start), format(data_end)
  ))

  # Snap the window to first-of-month boundaries and ensure the total month
  # count is divisible by BIN_MONTHS. If not, drop the trailing partial month.
  total_months <- diff_months(add_months(data_end, 1L), data_start)
  full_months  <- (total_months %/% config$BIN_MONTHS) * config$BIN_MONTHS
  window_start <- data_start
  window_end   <- add_months(window_start, full_months)
  n_bins       <- full_months %/% config$BIN_MONTHS

  log_msg(sprintf(
    "[%s] Full-dataset window: %s → %s (%d months, %d bins of %d months).",
    dataset_name, format(window_start), format(window_end),
    full_months, n_bins, config$BIN_MONTHS
  ))

  # Slice the feature matrix to the full-dataset window [start, end).
  in_window <- fm$Date >= window_start & fm$Date < window_end
  win_mat   <- fm[in_window, , drop = FALSE]
  n_seqs    <- nrow(win_mat)
  log_msg(sprintf("[%s] %d sequences retained after snap-to-month-boundary.",
                  dataset_name, n_seqs))

  site_cols <- setdiff(colnames(fm), c("Date", "Country"))

  # Bin sequence counts and date endpoints.
  bin_counts <- count_sequences_per_bin(win_mat$Date, window_start, window_end,
                                        config$BIN_MONTHS)
  bin_starts <- bin_endpoints(window_start, window_end, config$BIN_MONTHS)
  bin_labels <- format(bin_starts[-length(bin_starts)], "%Y-%m")

  log_msg(sprintf("[%s] Bin sequence count: median=%d, min=%d, max=%d.",
                  dataset_name, median(bin_counts),
                  min(bin_counts), max(bin_counts)))

  # Partition + Hellinger.
  partitions <- partition_window(
    feature_matrix = win_mat,
    window_start   = window_start,
    window_end     = window_end,
    site_cols      = site_cols,
    bin_months     = config$BIN_MONTHS
  )

  log_msg(sprintf("[%s] Computing full-1,273-site Hellinger trajectory...",
                  dataset_name))
  t_hell <- Sys.time()
  hellinger <- compute_hellinger_full(
    partitions  = partitions,
    sites       = seq_len(length(site_cols)),
    aa_levels   = config$AA_LEVELS,
    normalized  = config$HELLINGER_NORMALIZED
  )
  log_msg(sprintf("[%s] Hellinger computed in %.1f s; matrix shape: %d × %d.",
                  dataset_name,
                  as.numeric(difftime(Sys.time(), t_hell, units = "secs")),
                  nrow(hellinger), ncol(hellinger)))

  dat_t       <- t(hellinger)
  n_hell_rows <- nrow(dat_t)

  # Truths for the full window (no sub-windowing; truth set is the entire
  # catalogue intersected with the dataset's bin grid).
  truths  <- truths_in_window(truth_df, window_start, window_end)
  mapping <- map_truth_to_hellinger_idx(truths$Date_First_Detected_US,
                                        window_start, config$BIN_MONTHS)

  log_msg(sprintf(
    "[%s] Truth catalogue in range: %d events → %d unique non-T1 bins.",
    dataset_name, mapping$events_raw, mapping$bins_unique
  ))

  # Run all 9 methods.
  log_msg(sprintf("[%s] Running all 9 method configurations...", dataset_name))
  methods <- run_all_methods(dat_t, config)

  # Compute metrics against the truth set (single-call, no tier).
  metrics_eval <- evaluate_methods(
    method_results   = methods,
    truth_idx        = mapping$hellinger_idx,
    n_hellinger_rows = n_hell_rows,
    tolerance        = 1L,
    truth_shifts     = config$TRUTH_SHIFT_BINS
  )

  for (k in names(methods)) {
    m <- metrics_eval[[k]]
    log_msg(sprintf(
      "[%s] %-18s status=%s n_detected=%d P=%.2f R=%.2f F1=%.2f TLE=%s",
      dataset_name, k, m$status, length(m$detected_cps),
      m$metrics_primary$P, m$metrics_primary$R, m$metrics_primary$F1,
      if (is.na(m$metrics_primary$TLE)) "NA"
      else sprintf("%.2f", m$metrics_primary$TLE)
    ))
  }
  
  # ── GMM-based reduced-site analysis ────────────────────────────────────
  log_msg(sprintf("[%s] Fitting full-dataset GMM for site reduction...",
                  dataset_name))
  win_seq_matrix <- as.matrix(win_mat[, site_cols, drop = FALSE])
  gmm_meta <- fit_window_gmm(win_seq_matrix, config)
  log_msg(sprintf(
    "[%s] GMM: G=%s, modelName=%s, class-1 sites=%d, skipped=%s",
    dataset_name,
    if (is.na(gmm_meta$G))         "NA" else as.character(gmm_meta$G),
    if (is.na(gmm_meta$modelName)) "NA" else gmm_meta$modelName,
    gmm_meta$n_class1_sites, isTRUE(gmm_meta$reduced_skipped)
  ))
  
  if (isTRUE(gmm_meta$reduced_skipped) ||
      length(gmm_meta$class1_sites) < 2L) {
    hellinger_reduced <- NULL
    methods_reduced   <- NULL
    metrics_reduced   <- NULL
    reduced_status    <- "skipped"
    reduced_reason    <- if (isTRUE(gmm_meta$reduced_skipped))
      "GMM returned 999 sentinel or zero class-1 sites"
    else
      sprintf("only %d class-1 sites (< 2)",
              length(gmm_meta$class1_sites))
    log_msg(sprintf("[%s] Reduced-site analysis SKIPPED: %s",
                    dataset_name, reduced_reason))
  } else {
    log_msg(sprintf(
      "[%s] Computing reduced-site Hellinger trajectory on %d class-1 sites...",
      dataset_name, gmm_meta$n_class1_sites))
    t_red <- Sys.time()
    hellinger_reduced <- compute_hellinger_full(
      partitions  = partitions,
      sites       = gmm_meta$class1_sites,
      aa_levels   = config$AA_LEVELS,
      normalized  = config$HELLINGER_NORMALIZED
    )
    log_msg(sprintf(
      "[%s] Reduced Hellinger computed in %.1f s; matrix shape: %d × %d.",
      dataset_name,
      as.numeric(difftime(Sys.time(), t_red, units = "secs")),
      nrow(hellinger_reduced), ncol(hellinger_reduced)
    ))
    dat_t_red <- t(hellinger_reduced)
    methods_reduced <- run_all_methods(dat_t_red, config)
    metrics_reduced <- evaluate_methods(
      method_results   = methods_reduced,
      truth_idx        = mapping$hellinger_idx,
      n_hellinger_rows = nrow(dat_t_red),
      tolerance        = 1L,
      truth_shifts     = config$TRUTH_SHIFT_BINS
    )
    for (k in names(methods_reduced)) {
      m <- metrics_reduced[[k]]
      log_msg(sprintf(
        "[%s] %-18s status=%s n_detected=%d P=%.2f R=%.2f F1=%.2f TLE=%s (reduced)",
        dataset_name, k, m$status, length(m$detected_cps),
        m$metrics_primary$P, m$metrics_primary$R, m$metrics_primary$F1,
        if (is.na(m$metrics_primary$TLE)) "NA"
        else sprintf("%.2f", m$metrics_primary$TLE)
      ))
    }
    reduced_status <- "ok"
    reduced_reason <- NA_character_
  }
  
  result <- list(
    dataset             = dataset_name,
    window_start        = window_start,
    window_end          = window_end,
    n_bins              = n_bins,
    n_seqs              = n_seqs,
    bin_starts          = bin_starts,           # length n_bins + 1
    bin_labels          = bin_labels,           # length n_bins
    bin_counts          = bin_counts,           # length n_bins
    hellinger           = hellinger,            # sites × time_steps (full)
    n_hellinger_rows    = n_hell_rows,
    
    # GMM-based reduced-site analysis
    gmm_meta            = list(
      G                 = gmm_meta$G,
      modelName         = gmm_meta$modelName,
      n_class1_sites    = gmm_meta$n_class1_sites,
      class1_sites      = gmm_meta$class1_sites,
      reduced_skipped   = gmm_meta$reduced_skipped
    ),
    hellinger_reduced   = hellinger_reduced,    # class-1 sites × time_steps (or NULL)
    metrics_reduced     = metrics_reduced,      # list per method (or NULL)
    reduced_status      = reduced_status,
    reduced_reason      = reduced_reason,
    
    truth_dates         = truths$Date_First_Detected_US,
    truth_labels        = truths$WHO_Label,
    truth_hellinger_idx = mapping$hellinger_idx,
    truths_dropped_T1   = mapping$dropped_T1,
    truths_raw          = mapping$events_raw,
    truths_unique_bins  = mapping$bins_unique,

    metrics             = metrics_eval,

    r_version           = R.version.string,
    package_version     = as.character(utils::packageVersion("ViralEntropR")),
    run_timestamp       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )

  save_rds_atomic(result, out_path)
  log_msg(sprintf("[%s] Detection saved: %s", dataset_name, out_path))

  rm(fm); gc()
  invisible(result)
}

# ── Main loop ────────────────────────────────────────────────────────────────

for (ds_name in names(config$DATASETS)) {
  log_msg("==============================================================")
  log_msg(sprintf("Dataset: %s", ds_name))
  log_msg("==============================================================")
  tryCatch(
    run_one_full_dataset(ds_name, config, truth_df),
    error = function(e) {
      log_error(config, sprintf("[%s] FULL-DATASET DETECTION ERROR: %s",
                                ds_name, conditionMessage(e)))
      log_msg(sprintf("[%s] ERROR: %s", ds_name, conditionMessage(e)))
    }
  )
}

log_msg("==============================================================")
log_msg("Full-dataset detection complete.")
log_msg("==============================================================")
