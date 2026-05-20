################################################################################
## precompute_windows.R — Changepoint_Detection_Study
##
## Builds the per-dataset window manifest: 1,667 accepted windows per tier
## (short / medium / long), 5,001 windows per dataset. Deterministic in
## BASE_SEED. Saved to outputs/windows_<dataset>.rds for downstream workers.
##
## Manifest schema (one row per accepted window):
##   dataset           : character (e.g. "NCBI_US")
##   dataset_id        : integer (1 or 2)
##   cell_id           : integer 1..5001, unique within dataset
##   tier              : character ("short" / "medium" / "long")
##   tier_id           : integer 1..3
##   run_in_tier       : integer 1..1667
##   window_start      : Date (first of month)
##   window_end        : Date (first of month, exclusive)
##   length_months     : integer
##   n_bins            : integer
##   n_seqs            : integer (total sequences across the window)
##   bin_counts        : list-col, integer vector of length n_bins
##   truths_raw        : integer (truths in [start, end) before T1-drop)
##   truths_dropped_T1 : integer
##   truths_effective  : integer (truths after T1-drop; ≥ MIN_TRUTHS_PER_WINDOW)
##   truth_dates       : list-col, Date vector
##   truth_labels      : list-col, character vector (WHO_Label)
##   truth_hellinger_idx : list-col, integer vector (Hellinger row indices)
##   seed_replicate    : integer (per-replicate RNG seed)
################################################################################

# precompute_windows.R requires: setup.R, helpers_windows.R sourced upstream.

precompute_windows_for_dataset <- function(dataset_name, config,
                                           truth_df, dates,
                                           verbose = TRUE) {

  ds  <- config$DATASETS[[dataset_name]]
  if (is.null(ds))
    stop("Unknown dataset: ", dataset_name, call. = FALSE)

  dataset_id <- ds$id
  data_start <- min(dates)
  data_end   <- max(dates)
  max_window_end <- add_months(data_end, 1L)   # see helpers_windows.R rationale

  if (verbose) {
    log_msg(sprintf(
      "[%s] Window pre-computation: data_start = %s, data_end = %s, max_window_end = %s",
      dataset_name, format(data_start), format(data_end), format(max_window_end)
    ))
  }

  tier_names <- c("short", "medium", "long")
  tier_results <- vector("list", length(tier_names))
  names(tier_results) <- tier_names

  for (t_idx in seq_along(tier_names)) {
    tier_name <- tier_names[t_idx]
    tier_id   <- t_idx
    lengths   <- config$TIER_LENGTHS[[tier_name]]
    if (is.null(lengths))
      stop("No tier lengths defined for: ", tier_name, call. = FALSE)

    seed <- seed_for_window_pool(config$BASE_SEED, dataset_id, tier_id)

    if (verbose) {
      log_msg(sprintf(
        "[%s] Tier '%s' (id %d): target %d windows, valid lengths = {%s}, seed = %d",
        dataset_name, tier_name, tier_id, config$N_RUNS_PER_TIER,
        paste(lengths, collapse = ", "), seed
      ))
    }

    res <- draw_tier(
      tier_name      = tier_name,
      tier_lengths   = lengths,
      data_start     = data_start,
      max_window_end = max_window_end,
      dates          = dates,
      truth_df       = truth_df,
      n_target       = config$N_RUNS_PER_TIER,
      seed           = seed,
      config         = config
    )

    tier_results[[tier_name]] <- res

    if (verbose) {
      reason_tbl <- if (length(res$reasons) > 0L) {
        sort(table(res$reasons), decreasing = TRUE)
      } else {
        integer(0L)
      }
      log_msg(sprintf(
        "[%s] Tier '%s': accepted %d / %d after %d attempts (%.1f%% acceptance).",
        dataset_name, tier_name, res$n_kept, res$n_target,
        res$attempts, 100 * res$n_kept / max(1L, res$attempts)
      ))
      if (length(reason_tbl) > 0L) {
        log_msg(sprintf(
          "[%s] Tier '%s' rejection reasons (top 5):", dataset_name, tier_name))
        top5 <- head(reason_tbl, 5L)
        for (i in seq_along(top5)) {
          log_msg(sprintf("    %s : %d", names(top5)[i], top5[[i]]))
        }
      }
    }
  }

  # ── Assemble manifest data.frame ───────────────────────────────────────────

  rows <- list()
  cell_id <- 0L
  for (t_idx in seq_along(tier_names)) {
    tier_name <- tier_names[t_idx]
    tier_id   <- t_idx
    for (run_in_tier in seq_along(tier_results[[tier_name]]$accepted)) {
      cell_id <- cell_id + 1L
      a <- tier_results[[tier_name]]$accepted[[run_in_tier]]

      seed_rep <- seed_for_replicate(
        config$BASE_SEED, dataset_id, tier_id, run_in_tier
      )

      rows[[cell_id]] <- data.frame(
        dataset           = dataset_name,
        dataset_id        = dataset_id,
        cell_id           = cell_id,
        tier              = tier_name,
        tier_id           = tier_id,
        run_in_tier       = run_in_tier,
        window_start      = a$window_start,
        window_end        = a$window_end,
        length_months     = a$length_months,
        n_bins            = a$n_bins,
        n_seqs            = sum(a$bin_counts),
        truths_raw        = a$truths_raw,
        truths_dropped_T1 = a$truths_dropped_T1,
        truths_effective  = a$truths_effective,
        seed_replicate    = seed_rep,
        stringsAsFactors  = FALSE
      )
      # Attach list-column payloads as attributes; they'll be reconstructed
      # after rbind via a side-channel list. (rbind on data.frames drops
      # list columns silently; safest to store separately.)
      attr(rows[[cell_id]], "bin_counts")          <- a$bin_counts
      attr(rows[[cell_id]], "truth_dates")         <- a$truth_dates
      attr(rows[[cell_id]], "truth_labels")        <- a$truth_labels
      attr(rows[[cell_id]], "truth_hellinger_idx") <- a$truth_hellinger_idx
    }
  }

  manifest <- do.call(rbind, rows)
  manifest$bin_counts          <- lapply(rows, attr, "bin_counts")
  manifest$truth_dates         <- lapply(rows, attr, "truth_dates")
  manifest$truth_labels        <- lapply(rows, attr, "truth_labels")
  manifest$truth_hellinger_idx <- lapply(rows, attr, "truth_hellinger_idx")
  rownames(manifest) <- NULL

  # ── Save to disk ──────────────────────────────────────────────────────────
  out_path <- windows_path(config, dataset_name)
  save_rds_atomic(manifest, out_path)
  if (verbose)
    log_msg(sprintf("[%s] Manifest saved (%d rows): %s",
                    dataset_name, nrow(manifest), out_path))

  manifest
}

# Convenience wrapper that loads the feature matrix and truth catalogue.
precompute_all_windows <- function(config, verbose = TRUE) {
  truth_df <- load_truth_catalogue()
  if (verbose) {
    log_msg(sprintf("Truth catalogue: %d VOC/VOI variants with US dates.",
                    nrow(truth_df)))
  }

  manifests <- list()
  for (ds_name in names(config$DATASETS)) {
    rds_path <- config$DATASETS[[ds_name]]$feature_rds
    if (verbose)
      log_msg(sprintf("[%s] Loading feature matrix: %s", ds_name, rds_path))
    fm <- load_feature_matrix(rds_path)
    if (verbose)
      log_msg(sprintf("[%s] Feature matrix: %d sequences × %d columns.",
                      ds_name, nrow(fm), ncol(fm)))

    manifests[[ds_name]] <- precompute_windows_for_dataset(
      dataset_name = ds_name,
      config       = config,
      truth_df     = truth_df,
      dates        = fm$Date,
      verbose      = verbose
    )
    rm(fm); gc()
  }
  manifests
}
