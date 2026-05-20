################################################################################
## simulation_study.R — Changepoint_Detection_Study
##
## Orchestrator. Iterates over (dataset, cell_id), dispatching each replicate
## either sequentially (in-process for-loop) or in parallel via callr::r_bg
## subprocess pool. Mirrors the architecture of
## Sample_Size_Simulation_Study/simulation_study.R.
##
## Usage:
##   Rscript simulation_study.R
##
## Env-var overrides:
##   N_WORKERS=24                 number of parallel workers (default: config)
##   RUN_TRUTH_SHIFT_SWEEP=1      enable truth-shift sensitivity (default: 1)
##   VIRAL_CP_STUDY_DIR=/path     override study directory (read by setup.R)
##   FEATURE_RDS_NCBI_US=/path    override NCBI feature matrix path
##   FEATURE_RDS_GISAID_US=/path  override GISAID feature matrix path
##
## Resumability: on startup, the orchestrator scans
##   outputs/replicates_<dataset>/run_<NNNN>.rds
## and skips replicates already on disk. Interrupted runs can be resumed by
## re-executing the script.
################################################################################

suppressPackageStartupMessages({
  source("setup.R")
  source("helpers_windows.R")
  source("helpers_hellinger.R")
  source("cp_methods.R")
  source("metrics.R")
  source("precompute_windows.R")
  source("run_one_replicate.R")
})

# ── 1. Configuration + env overrides ─────────────────────────────────────────

config <- build_config()

env_n_workers <- Sys.getenv("N_WORKERS", unset = "")
if (nzchar(env_n_workers)) {
  config$N_WORKERS <- suppressWarnings(as.integer(env_n_workers))
  if (is.na(config$N_WORKERS) || config$N_WORKERS < 1L)
    stop("N_WORKERS env override must be a positive integer; got: ",
         env_n_workers, call. = FALSE)
}

env_rob <- Sys.getenv("RUN_TRUTH_SHIFT_SWEEP", unset = "")
if (nzchar(env_rob))
  config$RUN_TRUTH_SHIFT_SWEEP <- tolower(env_rob) %in% c("1", "true", "yes")

# Auto-detect N_WORKERS if unset.
if (is.null(config$N_WORKERS) || is.na(config$N_WORKERS))
  config$N_WORKERS <- auto_detect_n_workers()
config$N_WORKERS <- as.integer(config$N_WORKERS)

# Validate callr availability for parallel mode.
if (config$N_WORKERS > 1L && !requireNamespace("callr", quietly = TRUE))
  stop("Parallel mode (N_WORKERS > 1) requires the `callr` package.",
       call. = FALSE)

# Ensure output directory exists.
dir.create(output_dir(config), recursive = TRUE, showWarnings = FALSE)

# ── 2. Startup banner ────────────────────────────────────────────────────────

log_msg("==============================================================")
log_msg("Changepoint_Detection_Study — orchestrator startup")
log_msg("==============================================================")
log_msg("Configuration:")
log_msg(sprintf("  STUDY_DIR              = %s", config$STUDY_DIR))
log_msg(sprintf("  BASE_SEED              = %d", config$BASE_SEED))
log_msg(sprintf("  N_RUNS_PER_TIER        = %d", config$N_RUNS_PER_TIER))
log_msg(sprintf("  N_WORKERS              = %d (%s)",
                config$N_WORKERS,
                if (config$N_WORKERS == 1L) "sequential" else "parallel via callr"))
log_msg(sprintf("  RUN_TRUTH_SHIFT_SWEEP  = %s", config$RUN_TRUTH_SHIFT_SWEEP))
log_msg(sprintf("  MIN_SEQUENCES_PER_BIN  = %d", config$MIN_SEQUENCES_PER_BIN))
log_msg(sprintf("  MIN_TRUTHS_PER_WINDOW  = %d", config$MIN_TRUTHS_PER_WINDOW))
log_msg(sprintf("  BIN_MONTHS             = %d", config$BIN_MONTHS))
log_msg(sprintf("  Datasets               = %s",
                paste(names(config$DATASETS), collapse = ", ")))

# ── 3. Build / load window manifests ─────────────────────────────────────────

manifests <- list()
for (ds_name in names(config$DATASETS)) {
  wp <- windows_path(config, ds_name)
  if (file.exists(wp)) {
    log_msg(sprintf("[%s] Loading existing manifest: %s", ds_name, wp))
    manifests[[ds_name]] <- readRDS(wp)
    log_msg(sprintf("[%s] Manifest has %d windows.",
                    ds_name, nrow(manifests[[ds_name]])))
  } else {
    log_msg(sprintf("[%s] Manifest not found; building...", ds_name))
    truth_df <- load_truth_catalogue()
    rds_path <- config$DATASETS[[ds_name]]$feature_rds
    log_msg(sprintf("[%s] Loading feature matrix: %s", ds_name, rds_path))
    fm <- load_feature_matrix(rds_path)
    log_msg(sprintf("[%s] Feature matrix: %d sequences × %d columns.",
                    ds_name, nrow(fm), ncol(fm)))
    manifests[[ds_name]] <- precompute_windows_for_dataset(
      dataset_name = ds_name,
      config       = config,
      truth_df     = truth_df,
      dates        = fm$Date,
      verbose      = TRUE
    )
    rm(fm); gc()
  }
}

# ── 4. Sequential execution path ─────────────────────────────────────────────

run_dataset_sequential <- function(dataset_name, manifest, config) {

  log_msg(sprintf("[%s] Sequential mode: loading feature matrix once...",
                  dataset_name))
  fm <- load_feature_matrix(config$DATASETS[[dataset_name]]$feature_rds)

  cell_ids <- manifest$cell_id
  rds_dir  <- replicates_dir(config, dataset_name)
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

  # Resume-scan.
  existing <- file.exists(file.path(
    rds_dir, sprintf("run_%04d.rds", cell_ids)
  ))
  todo <- cell_ids[!existing]
  log_msg(sprintf("[%s] Resume scan: %d / %d already complete; %d to run.",
                  dataset_name, sum(existing), length(cell_ids), length(todo)))

  for (i in seq_along(todo)) {
    cid <- todo[i]
    tryCatch({
      run_one_replicate(dataset_name, cid, manifest, fm, config)
    }, error = function(e) {
      log_error(config, sprintf("[%s cell %d] ERROR: %s",
                                dataset_name, cid, conditionMessage(e)))
    })
    if (i %% 100L == 0L)
      log_msg(sprintf("[%s] Progress: %d / %d", dataset_name, i, length(todo)))
  }
  rm(fm); gc()
  log_msg(sprintf("[%s] Sequential pass complete.", dataset_name))
}

# ── 5. Parallel execution path (callr) ───────────────────────────────────────
#
# Each subprocess sources the orchestrator's helper files, loads the feature
# matrix and the manifest, and runs one replicate. The master polls for
# completion and dispatches new subprocesses to maintain a pool of size
# config$N_WORKERS. Failures are logged but do not interrupt the run.

launch_replicate_subprocess <- function(dataset_name, cell_id, config) {

  payload <- list(
    config       = config,
    dataset_name = dataset_name,
    cell_id      = as.integer(cell_id)
  )

  callr::r_bg(
    func = function(payload) {
      setwd(payload$config$STUDY_DIR)
      source("setup.R")
      source("helpers_windows.R")
      source("helpers_hellinger.R")
      source("cp_methods.R")
      source("metrics.R")
      source("run_one_replicate.R")

      ds_name  <- payload$dataset_name
      cid      <- payload$cell_id

      manifest <- readRDS(windows_path(payload$config, ds_name))
      fm <- load_feature_matrix(
        payload$config$DATASETS[[ds_name]]$feature_rds
      )
      run_one_replicate(ds_name, cid, manifest, fm, payload$config)
    },
    args        = list(payload = payload),
    stderr      = file.path(tempdir(), sprintf("cp_stderr_%s_%04d.log",
                                               dataset_name, cell_id)),
    stdout      = file.path(tempdir(), sprintf("cp_stdout_%s_%04d.log",
                                               dataset_name, cell_id)),
    supervise   = TRUE
  )
}

run_dataset_parallel <- function(dataset_name, manifest, config) {

  cell_ids <- manifest$cell_id
  rds_dir  <- replicates_dir(config, dataset_name)
  dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

  # Resume-scan.
  existing <- file.exists(file.path(
    rds_dir, sprintf("run_%04d.rds", cell_ids)
  ))
  todo <- cell_ids[!existing]
  log_msg(sprintf("[%s] Resume scan: %d / %d already complete; %d to run.",
                  dataset_name, sum(existing), length(cell_ids), length(todo)))
  if (length(todo) == 0L) {
    log_msg(sprintf("[%s] All replicates already complete; skipping dispatch.",
                    dataset_name))
    return(invisible(NULL))
  }

  # Worker pool dispatch.
  poll_s <- config$SUBPROCESS_POLL_S
  n_workers <- config$N_WORKERS
  active   <- list()   # named list(cell_id = process)
  cursor   <- 1L
  n_done   <- 0L
  n_failed <- 0L
  t0       <- Sys.time()

  while (cursor <= length(todo) || length(active) > 0L) {

    # Fill the pool.
    while (length(active) < n_workers && cursor <= length(todo)) {
      cid <- todo[cursor]
      proc <- tryCatch(
        launch_replicate_subprocess(dataset_name, cid, config),
        error = function(e) {
          log_error(config, sprintf("[%s cell %d] LAUNCH ERROR: %s",
                                    dataset_name, cid, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(proc)) {
        active[[as.character(cid)]] <- proc
      }
      cursor <- cursor + 1L
    }

    # Poll active subprocesses.
    Sys.sleep(poll_s)
    finished_keys <- character(0L)
    for (key in names(active)) {
      proc <- active[[key]]
      if (!proc$is_alive()) {
        finished_keys <- c(finished_keys, key)
        exit_status <- tryCatch(proc$get_exit_status(), error = function(e) NA)
        if (is.na(exit_status) || exit_status != 0L) {
          n_failed <- n_failed + 1L
          err_out <- tryCatch(proc$read_all_error(), error = function(e) "")
          log_error(config, sprintf(
            "[%s cell %s] SUBPROCESS EXIT %s: %s",
            dataset_name, key, as.character(exit_status), substr(err_out, 1L, 500L)
          ))
        } else {
          n_done <- n_done + 1L
        }
      }
    }
    for (key in finished_keys) active[[key]] <- NULL

    # Progress.
    if (n_done > 0L && (n_done + n_failed) %% 100L == 0L) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      log_msg(sprintf(
        "[%s] Progress: %d done, %d failed of %d (%.1f min elapsed)",
        dataset_name, n_done, n_failed, length(todo), elapsed
      ))
    }
  }

  elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  log_msg(sprintf(
    "[%s] Parallel pass complete: %d done, %d failed, %.1f min total.",
    dataset_name, n_done, n_failed, elapsed_min))
}

# ── 6. Aggregation: per-dataset summary table ────────────────────────────────
#
# Aggregates per-replicate RDS files into a flat per-(method, site_set,
# replicate) data.frame for downstream plotting and statistical analysis.

aggregate_dataset <- function(dataset_name, config) {

  rds_dir <- replicates_dir(config, dataset_name)
  files   <- list.files(rds_dir, pattern = "^run_\\d{4}\\.rds$",
                        full.names = TRUE)
  if (length(files) == 0L) {
    log_msg(sprintf("[%s] Aggregation: no per-replicate files; skipping.",
                    dataset_name))
    return(invisible(NULL))
  }

  rows <- vector("list", length(files))
  for (i in seq_along(files)) {
    r <- tryCatch(readRDS(files[i]), error = function(e) NULL)
    if (is.null(r)) next

    base_meta <- data.frame(
      dataset           = r$dataset,
      cell_id           = r$cell_id,
      tier              = r$tier,
      tier_id           = r$tier_id,
      run_in_tier       = r$run_in_tier,
      window_start      = r$window_start,
      window_end        = r$window_end,
      length_months     = r$length_months,
      n_bins            = r$n_bins,
      n_seqs            = r$n_seqs,
      truths_effective  = r$truths_effective,
      gmm_G             = r$gmm_meta$G,
      gmm_modelName     = r$gmm_meta$modelName,
      n_class1_sites    = r$gmm_meta$n_class1_sites,
      reduced_status    = r$reduced_status,
      walltime_s        = r$walltime_s,
      seed              = r$seed,
      stringsAsFactors  = FALSE
    )

    # Long-format: one row per (method, site_set).
    for (site_set in c("full", "reduced")) {
      metrics_blk <- if (site_set == "full") r$metrics_full else r$metrics_reduced
      if (is.null(metrics_blk)) next

      for (method_key in names(metrics_blk)) {
        m    <- metrics_blk[[method_key]]
        prim <- m$metrics_primary
        rows[[length(rows) + 1L]] <- cbind(
          base_meta,
          data.frame(
            site_set     = site_set,
            method       = method_key,
            status       = m$status,
            n_detected   = length(m$detected_cps),
            P            = prim$P,
            R            = prim$R,
            F1           = prim$F1,
            TLE          = prim$TLE,
            TP           = prim$TP,
            FP           = prim$FP,
            FN           = prim$FN,
            method_walltime_s = m$walltime_s,
            stringsAsFactors  = FALSE
          )
        )
      }
    }
  }
  # Strip pre-allocated NULL entries.
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L)
    return(invisible(NULL))

  summary_df <- do.call(rbind, rows)
  rownames(summary_df) <- NULL
  save_rds_atomic(summary_df, summary_path(config, dataset_name))
  log_msg(sprintf("[%s] Summary written: %s (%d rows).",
                  dataset_name, summary_path(config, dataset_name),
                  nrow(summary_df)))
  invisible(summary_df)
}

# ── 7. Main loop ─────────────────────────────────────────────────────────────

for (ds_name in names(manifests)) {
  log_msg("==============================================================")
  log_msg(sprintf("Dataset: %s (%d windows)", ds_name, nrow(manifests[[ds_name]])))
  log_msg("==============================================================")

  if (config$N_WORKERS == 1L) {
    run_dataset_sequential(ds_name, manifests[[ds_name]], config)
  } else {
    run_dataset_parallel(ds_name, manifests[[ds_name]], config)
  }

  aggregate_dataset(ds_name, config)
  gc()
}

log_msg("==============================================================")
log_msg("All datasets complete.")
log_msg("==============================================================")
