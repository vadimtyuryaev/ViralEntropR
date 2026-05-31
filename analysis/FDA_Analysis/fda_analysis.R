################################################################################
## fda_analysis.R — FDA_Analysis
##
## Top-level orchestrator. Dispatches per-cell jobs across the cross product of
## (dataset × variant × strategy) cells, sequentially or in parallel via
## callr::r_bg. Each completed cell is saved to outputs/cells/cell_*.rds.
##
## Invocation:
##   Rscript fda_analysis.R
##
## Environment overrides (optional):
##   N_WORKERS              integer, parallel workers (default: auto-detect)
##   VIRAL_FDA_STUDY_DIR    path to FDA_Analysis/  (default: auto-detect)
##   FEATURE_RDS_NCBI_US    override NCBI feature-matrix RDS path
##   FEATURE_RDS_GISAID_US  override GISAID feature-matrix RDS path
##
## Resume safety: cells with an existing RDS at outputs/cells/cell_*.rds are
## skipped. Stub RDS files (for out_of_coverage, reduced_skipped, etc.) count
## as completed.
################################################################################

suppressPackageStartupMessages({
  source("setup.R")
  source("helpers_fda.R")
  source("helpers_frames.R")
  source("run_one_cell.R")
})

require_runtime_packages(animation = FALSE)


# ── 1. Configuration & banner ────────────────────────────────────────────────

config <- build_config()
setwd(config$STUDY_DIR)

dir.create(output_dir(config), recursive = TRUE, showWarnings = FALSE)
dir.create(cells_dir (config), recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir  (config), recursive = TRUE, showWarnings = FALSE)

log_msg("================================================================")
log_msg("FDA_Analysis — orchestrator starting")
log_msg(sprintf("  STUDY_DIR : %s", config$STUDY_DIR))
log_msg(sprintf("  N_WORKERS : %d", config$N_WORKERS))
log_msg(sprintf("  BASE_SEED : %d", config$BASE_SEED))
log_msg(sprintf("  Variants  : %s", paste(config$VARIANTS_TO_RUN,
                                            collapse = ", ")))
log_msg(sprintf("  Strategies: %s", paste(config$STRATEGIES_TO_RUN,
                                            collapse = ", ")))
log_msg("================================================================")


# ── 2. Truth catalogue (loaded once, in every subprocess too) ───────────────

truth_df <- load_truth_catalogue(config$VARIANTS_TO_RUN)
log_msg(sprintf("Loaded truth catalogue: %d / %d variants resolved.",
                 nrow(truth_df), length(config$VARIANTS_TO_RUN)))
absent <- setdiff(config$VARIANTS_TO_RUN, truth_df$WHO_Label)
if (length(absent) > 0L)
  log_msg(sprintf("  Variants absent from truth: %s",
                   paste(absent, collapse = ", ")))


# ── 3. Cell enumeration & resume scan ───────────────────────────────────────

cells_df <- enumerate_cells(config)
cells_df$out_path <- vapply(seq_len(nrow(cells_df)), function(i)
  cell_path(config, cells_df$dataset[i], cells_df$variant[i],
             cells_df$strategy[i]),
  character(1L))

# Schema-aware resume: only treat cells as complete if they exist AND match
# the current SCHEMA_VERSION. Older RDS files are flagged stale.
schema_ok <- function(path) {
  if (!file.exists(path)) return(FALSE)
  tryCatch({
    r <- readRDS(path)
    isTRUE(identical(r$schema_version, config$SCHEMA_VERSION %||% 2L))
  }, error = function(e) FALSE)
}
cells_df$exists <- vapply(cells_df$out_path, schema_ok, logical(1L))

n_total    <- nrow(cells_df)
n_existing <- sum(cells_df$exists)
n_todo     <- n_total - n_existing

log_msg(sprintf("Resume scan: %d / %d cells already complete; %d to run.",
                 n_existing, n_total, n_todo))

if (n_todo == 0L) {
  log_msg("Nothing to do; exiting.")
  quit(save = "no", status = 0L)
}

todo <- cells_df[!cells_df$exists, , drop = FALSE]
rownames(todo) <- NULL


# ── 4. Per-dataset feature-matrix cache (for sequential mode) ───────────────

.fm_cache <- new.env(parent = emptyenv())

get_feature_matrix <- function(ds_name) {
  if (is.null(.fm_cache[[ds_name]])) {
    rds <- config$DATASETS[[ds_name]]$feature_rds
    log_msg(sprintf("Loading feature matrix [%s] from %s ...", ds_name, rds))
    .fm_cache[[ds_name]] <- load_feature_matrix(rds,
                                                  expected_n_sites = config$N_SITES)
    log_msg(sprintf("  → %d sequences × %d sites loaded.",
                     nrow(.fm_cache[[ds_name]]),
                     config$N_SITES))
  }
  .fm_cache[[ds_name]]
}


# ── 5. Sequential mode ──────────────────────────────────────────────────────

run_sequential <- function(todo, truth_df, config) {

  # Sort by dataset so we only swap the feature matrix when needed.
  todo <- todo[order(todo$dataset, todo$variant, todo$strategy), , drop = FALSE]
  rownames(todo) <- NULL

  t0  <- Sys.time()
  done <- 0L; err <- 0L; skip <- 0L

  for (i in seq_len(nrow(todo))) {
    row <- todo[i, ]
    cid <- sprintf("[%s | %s | %s]", row$dataset, row$variant, row$strategy)

    fm <- tryCatch(get_feature_matrix(row$dataset),
                   error = function(e) {
                     log_error(config, sprintf("%s feature-matrix load: %s",
                                                 cid, conditionMessage(e)))
                     NULL
                   })
    if (is.null(fm)) { err <- err + 1L; next }

    result <- tryCatch(
      run_one_cell(row$dataset, row$variant, row$strategy,
                    feature_matrix = fm,
                    truth_df       = truth_df,
                    config         = config),
      error = function(e) {
        log_error(config, sprintf("%s run_one_cell uncaught: %s",
                                    cid, conditionMessage(e)))
        NULL
      }
    )

    if (is.null(result)) {
      err <- err + 1L
    } else if (identical(result$status, "ok")) {
      done <- done + 1L
    } else {
      skip <- skip + 1L
    }

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    log_msg(sprintf("%s status=%s walltime=%.1fs  | totals: ok=%d skip=%d err=%d  (%.1f min elapsed)",
                     cid,
                     if (is.null(result)) "ERROR" else result$status,
                     if (is.null(result)) NA_real_ else result$walltime_s,
                     done, skip, err, elapsed))
  }

  list(done = done, skipped = skip, errored = err)
}


# ── 6. Parallel mode (callr::r_bg) ──────────────────────────────────────────
#
# Each subprocess:
#   • setwd(STUDY_DIR), sources setup.R / helpers_fda.R / run_one_cell.R
#   • Loads ITS dataset's feature matrix (≈ 5 s startup)
#   • Loads truth catalogue
#   • Runs one cell
#   • Writes the per-cell RDS via save_rds_atomic
#   • Exits

launch_cell_subprocess <- function(dataset_name, variant_name, strategy_name,
                                     config) {

  payload <- list(
    dataset_name  = dataset_name,
    variant_name  = variant_name,
    strategy_name = strategy_name,
    config        = config
  )
  tag <- sprintf("%s__%s__%s", dataset_name, variant_name, strategy_name)

  callr::r_bg(
    func = function(payload) {
      setwd(payload$config$STUDY_DIR)
      source("setup.R")
      source("helpers_fda.R")
      source("helpers_frames.R")
      source("run_one_cell.R")

      truth_df <- load_truth_catalogue(payload$config$VARIANTS_TO_RUN)
      fm <- load_feature_matrix(
        payload$config$DATASETS[[payload$dataset_name]]$feature_rds,
        expected_n_sites = payload$config$N_SITES
      )

      run_one_cell(payload$dataset_name,
                    payload$variant_name,
                    payload$strategy_name,
                    feature_matrix = fm,
                    truth_df       = truth_df,
                    config         = payload$config)
    },
    args      = list(payload = payload),
    stderr    = file.path(logs_dir(config),
                           sprintf("subproc_stderr_%s.log", tag)),
    stdout    = file.path(logs_dir(config),
                           sprintf("subproc_stdout_%s.log", tag)),
    supervise = TRUE
  )
}

run_parallel <- function(todo, config) {

  if (!.have_callr)
    stop("Parallel mode requires the 'callr' package.", call. = FALSE)

  todo <- todo[order(todo$dataset, todo$variant, todo$strategy), , drop = FALSE]
  rownames(todo) <- NULL

  poll_s    <- config$SUBPROCESS_POLL_S
  n_workers <- config$N_WORKERS
  active    <- list()
  cursor    <- 1L
  done <- 0L; err <- 0L
  t0   <- Sys.time()

  while (cursor <= nrow(todo) || length(active) > 0L) {

    # Fill pool.
    while (length(active) < n_workers && cursor <= nrow(todo)) {
      row <- todo[cursor, ]
      tag <- sprintf("%s__%s__%s", row$dataset, row$variant, row$strategy)

      proc <- tryCatch(
        launch_cell_subprocess(row$dataset, row$variant, row$strategy, config),
        error = function(e) {
          log_error(config, sprintf("[%s] LAUNCH ERROR: %s",
                                      tag, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(proc)) active[[tag]] <- proc
      cursor <- cursor + 1L
    }

    Sys.sleep(poll_s)

    finished_keys <- character(0L)
    for (key in names(active)) {
      proc <- active[[key]]
      if (!proc$is_alive()) {
        finished_keys <- c(finished_keys, key)
        exit <- tryCatch(proc$get_exit_status(), error = function(e) NA)
        if (is.na(exit) || exit != 0L) {
          err <- err + 1L
          stderr_tail <- tryCatch(proc$read_all_error(),
                                   error = function(e) "")
          log_error(config, sprintf("[%s] SUBPROCESS EXIT %s: %s",
                                      key, as.character(exit),
                                      substr(stderr_tail, 1L, 800L)))
        } else {
          done <- done + 1L
        }
      }
    }
    for (key in finished_keys) active[[key]] <- NULL

    if ((done + err) > 0L && (done + err) %% 5L == 0L) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      log_msg(sprintf("Progress: %d done, %d failed of %d  (%.1f min, pool=%d)",
                       done, err, nrow(todo), elapsed, length(active)))
    }
  }

  list(done = done, skipped = 0L, errored = err)
}


# ── 7. Dispatch ─────────────────────────────────────────────────────────────

t_dispatch <- Sys.time()
totals <- if (config$N_WORKERS == 1L) {
  run_sequential(todo, truth_df, config)
} else {
  run_parallel(todo, config)
}
elapsed_min <- as.numeric(difftime(Sys.time(), t_dispatch, units = "mins"))

log_msg("================================================================")
log_msg(sprintf("Dispatch complete: ok+skip=%d errored=%d  (%.1f min total)",
                 totals$done + totals$skipped, totals$errored, elapsed_min))
log_msg(sprintf("  RDS files at: %s", cells_dir(config)))
log_msg(sprintf("  Logs at:      %s", logs_dir(config)))
log_msg("================================================================")
