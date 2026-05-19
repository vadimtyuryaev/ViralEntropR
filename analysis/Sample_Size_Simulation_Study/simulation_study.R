# =============================================================================
# simulation_study.R - Orchestrator
# =============================================================================
#
# Entry point for the Sample-Size Simulation Study.
#
# Responsibilities:
#   1. Resolve the study directory, locate the reference FASTA, build the
#      configuration list, set N_WORKERS (auto-detect or env override).
#   2. Load the reference sequence and the empirical constants from
#      sarscov2_variants.
#   3. Build or load the four cells tables.
#   4. Per scenario: enumerate the (cell, rep) work queue, skip already-
#      completed replicates (resume from disk), and run remaining
#      replicates either sequentially (N_WORKERS = 1) or in parallel via
#      callr::r_bg with at most N_WORKERS concurrent processes.
#   5. After each scenario completes, aggregate the per-replicate RDS files
#      into a tidy summary data frame and save.
#   6. Optionally run the p_del robustness sweep on a stratified-random
#      subset of cells across {0, 1e-4, 1e-3, 1e-2}, writing outputs under
#      outputs/robustness_pdel/.
#
# Concurrency switch (single field):
#   N_WORKERS = 1L    -> sequential in-process execution
#   N_WORKERS > 1L    -> parallel callr pool with N_WORKERS slots
#
# Env overrides recognized at startup:
#   N_WORKERS=<int>             override config$N_WORKERS
#   RUN_ROBUSTNESS_SWEEP=1      enable robustness sweep
#   REF_SEQ_FASTA=<path>        override reference-FASTA location
#   VIRAL_STUDY_DIR=<path>      override study-directory resolution
#
# Invocation examples:
#   # Linux/macOS background run, 16 parallel workers, log to file:
#   N_WORKERS=16 nohup Rscript simulation_study.R > run.log 2>&1 &
#
#   # Sequential debug run on a laptop:
#   N_WORKERS=1 Rscript simulation_study.R
#
#   # Sequential debug + robustness sweep:
#   N_WORKERS=1 RUN_ROBUSTNESS_SWEEP=1 Rscript simulation_study.R
#
# Resume behaviour:
#   Implicit. At startup, the orchestrator scans each scenario's output
#   directory for existing per-replicate RDS files and only re-runs missing
#   ones. To force a clean restart, delete the relevant replicates_sc<N>/
#   directory before invocation.
#
# Author : Vadim Tyuryaev
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Source helpers
# -----------------------------------------------------------------------------
.script_path <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  fn_arg <- args[grep("^--file=", args)]
  if (length(fn_arg) > 0L) {
    normalizePath(sub("^--file=", "", fn_arg[1L]))
  } else if (sys.nframe() > 0L) {
    normalizePath(sys.frame(1)$ofile)
  } else {
    file.path(getwd(), "simulation_study.R")
  }
}, error = function(e) file.path(getwd(), "simulation_study.R"))

setwd(dirname(.script_path))

source("setup.R")
source("simulator.R")
source("detect_in_snapshot.R")
source("precompute_cells.R")
source("run_one_replicate.R")

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------
config <- build_config()
config$STUDY_DIR <- resolve_study_dir()

# Reference FASTA: defaults to spike_reference_1273aa.fasta inside the study
# directory; override via REF_SEQ_FASTA env var.
config$REF_SEQ_FASTA <- {
  override <- Sys.getenv("REF_SEQ_FASTA", unset = "")
  if (nzchar(override)) {
    override
  } else {
    file.path(config$STUDY_DIR, "spike_reference_1273aa.fasta")
  }
}

# N_WORKERS: env override > config build > auto-detect.
env_n_workers <- Sys.getenv("N_WORKERS", unset = "")
if (nzchar(env_n_workers)) {
  config$N_WORKERS <- as.integer(env_n_workers)
}
if (is.null(config$N_WORKERS) || is.na(config$N_WORKERS)) {
  config$N_WORKERS <- auto_detect_n_workers()
}
config$N_WORKERS <- as.integer(max(1L, config$N_WORKERS))

# Robustness-sweep flag.
env_rob <- tolower(Sys.getenv("RUN_ROBUSTNESS_SWEEP", unset = ""))
if (env_rob %in% c("1", "true", "yes")) config$RUN_ROBUSTNESS_SWEEP <- TRUE

# callr only required when running in parallel.
if (config$N_WORKERS > 1L && !requireNamespace("callr", quietly = TRUE))
  stop("Package 'callr' is required for parallel subprocess execution ",
       "(N_WORKERS > 1). Install via install.packages('callr'), or set ",
       "N_WORKERS=1 for sequential mode.", call. = FALSE)

log_msg("Configuration:")
log_msg(sprintf("  STUDY_DIR              = %s", config$STUDY_DIR))
log_msg(sprintf("  REF_SEQ_FASTA          = %s", config$REF_SEQ_FASTA))
log_msg(sprintf("  N_WORKERS              = %d (%s mode)",
                config$N_WORKERS,
                if (config$N_WORKERS == 1L) "sequential" else "parallel"))
log_msg(sprintf("  BASE_SEED              = %d", config$BASE_SEED))
log_msg(sprintf("  N_REPS_PER_CELL        = %d", config$N_REPS_PER_CELL))
log_msg(sprintf("  GRID_SIZE              = %d", config$GRID_SIZE))
log_msg(sprintf("  P_DEL                  = %g", config$P_DEL))
log_msg(sprintf("  MCLUST_MODELS          = %s",
                if (is.null(config$MCLUST_MODELS)) "NULL (default E and V via BIC)"
                else paste(config$MCLUST_MODELS, collapse = ",")))
log_msg(sprintf("  MCLUST_G               = %d..%d",
                min(config$MCLUST_G), max(config$MCLUST_G)))
log_msg(sprintf("  RUN_ROBUSTNESS_SWEEP   = %s",
                as.character(config$RUN_ROBUSTNESS_SWEEP)))

# -----------------------------------------------------------------------------
# 2. Reference sequence and empirical constants
# -----------------------------------------------------------------------------
log_msg("Loading reference sequence...")
ref_seq_int <- load_reference_sequence(config$REF_SEQ_FASTA, expected_L = 1273L)
log_msg(sprintf("Reference loaded: length = %d (codes %d:%d)",
                length(ref_seq_int), min(ref_seq_int), max(ref_seq_int)))

empirical <- build_empirical_constants()
log_msg(sprintf(
  "Empirical constants: |POOL_11|=%d, |POOL_12|=%d, MUT_COUNT_SET_11=[%s], OMICRON=%d",
  length(empirical$POOL_11), length(empirical$POOL_12),
  paste(empirical$MUT_COUNT_SET_11, collapse = ","),
  empirical$MUT_COUNT_OMICRON))

# -----------------------------------------------------------------------------
# 3. Build / load cell tables
# -----------------------------------------------------------------------------
dir.create(output_dir(config), recursive = TRUE, showWarnings = FALSE)
precompute_cells_for_all_scenarios(config)

# -----------------------------------------------------------------------------
# 4. Work-queue execution (sequential and parallel branches)
# -----------------------------------------------------------------------------

# Sequential executor: runs each replicate in-process. tryCatch isolates
# replicate-level failures so a single bad cell does not abort the run.
run_replicate_sequential <- function(scenario, cell_id, rep_id, cells, config,
                                     p_del_override = NULL,
                                     output_base_dir = NULL) {
  tryCatch({
    run_one_replicate(
      scenario        = scenario,
      cell_id         = cell_id,
      rep_id          = rep_id,
      cells           = cells,
      ref_seq_int     = ref_seq_int,
      empirical       = empirical,
      config          = config,
      p_del_override  = p_del_override,
      output_base_dir = output_base_dir
    )
    list(ok = TRUE, err = NULL)
  }, error = function(e) {
    list(ok = FALSE, err = conditionMessage(e))
  })
}

# Parallel executor: launches a fresh R session via callr::r_bg per
# replicate. The subprocess sources the study files, loads its own copies
# of the reference and empirical constants, reads the cells table, and
# runs one replicate.
launch_replicate_callr <- function(scenario, cell_id, rep_id, config,
                                   p_del_override  = NULL,
                                   output_base_dir = NULL) {
  callr::r_bg(
    function(scenario, cell_id, rep_id, config_payload,
             p_del_override, output_base_dir) {
      setwd(config_payload$STUDY_DIR)
      source("setup.R")
      source("simulator.R")
      source("detect_in_snapshot.R")
      source("run_one_replicate.R")

      ref_seq_int <- load_reference_sequence(
        config_payload$REF_SEQ_FASTA, expected_L = 1273L
      )
      empirical <- build_empirical_constants()
      cells     <- readRDS(cells_path(config_payload, scenario))

      run_one_replicate(
        scenario        = scenario,
        cell_id         = cell_id,
        rep_id          = rep_id,
        cells           = cells,
        ref_seq_int     = ref_seq_int,
        empirical       = empirical,
        config          = config_payload,
        p_del_override  = p_del_override,
        output_base_dir = output_base_dir
      )
    },
    args = list(scenario = scenario, cell_id = cell_id, rep_id = rep_id,
                config_payload = config,
                p_del_override = p_del_override,
                output_base_dir = output_base_dir),
    supervise = TRUE,
    stderr = tempfile(pattern = sprintf("sc%d_c%d_r%d_err_",
                                        scenario, cell_id, rep_id),
                      fileext = ".log"),
    stdout = tempfile(pattern = sprintf("sc%d_c%d_r%d_out_",
                                        scenario, cell_id, rep_id),
                      fileext = ".log")
  )
}

# Shared work-queue executor used by both the main study and the
# robustness sweep. Branches on config$N_WORKERS.
execute_work_queue <- function(scenario, work_queue, cells, config,
                               p_del_override  = NULL,
                               output_base_dir = NULL,
                               label = "main") {

  total_work <- nrow(work_queue)
  if (total_work == 0L) return(list(completed = 0L, failed = 0L))

  log_every <- max(50L, total_work %/% 100L)
  n_completed <- 0L
  n_failed    <- 0L

  if (config$N_WORKERS == 1L) {
    # ---- Sequential branch ----------------------------------------------
    for (i in seq_len(total_work)) {
      res <- run_replicate_sequential(
        scenario        = scenario,
        cell_id         = work_queue$cell_id[i],
        rep_id          = work_queue$rep_id[i],
        cells           = cells,
        config          = config,
        p_del_override  = p_del_override,
        output_base_dir = output_base_dir
      )
      if (res$ok) {
        n_completed <- n_completed + 1L
      } else {
        n_failed <- n_failed + 1L
        log_err(config, sprintf(
          "Replicate FAILED [%s] sc=%d cell=%d rep=%d: %s",
          label, scenario, work_queue$cell_id[i],
          work_queue$rep_id[i], res$err))
      }
      if ((n_completed + n_failed) %% log_every == 0L ||
          (n_completed + n_failed) == total_work) {
        log_msg(sprintf("    [%s] progress: %d done (%d ok, %d failed) / %d (%.1f%%)",
                        label, n_completed + n_failed, n_completed, n_failed,
                        total_work,
                        100 * (n_completed + n_failed) / total_work))
      }
    }
  } else {
    # ---- Parallel callr pool branch -------------------------------------
    next_idx     <- 1L
    active_procs <- list()

    while (n_completed + n_failed < total_work) {
      # Top up active slots
      while (length(active_procs) < config$N_WORKERS &&
             next_idx <= total_work) {
        proc <- launch_replicate_callr(
          scenario        = scenario,
          cell_id         = work_queue$cell_id[next_idx],
          rep_id          = work_queue$rep_id[next_idx],
          config          = config,
          p_del_override  = p_del_override,
          output_base_dir = output_base_dir
        )
        active_procs[[length(active_procs) + 1L]] <- list(
          proc    = proc,
          cell_id = work_queue$cell_id[next_idx],
          rep_id  = work_queue$rep_id[next_idx]
        )
        next_idx <- next_idx + 1L
      }

      # Poll for completion
      Sys.sleep(config$SUBPROCESS_POLL_S)
      done_idx <- integer(0L)
      for (i in seq_along(active_procs)) {
        if (!active_procs[[i]]$proc$is_alive()) {
          done_idx <- c(done_idx, i)
        }
      }

      if (length(done_idx) > 0L) {
        for (i in rev(done_idx)) {
          done <- active_procs[[i]]
          exit_code <- tryCatch(done$proc$get_exit_status(),
                                error = function(e) -1L)
          if (is.null(exit_code) || exit_code != 0L) {
            err_text <- tryCatch(readLines(done$proc$get_error_file(),
                                           warn = FALSE),
                                 error = function(e) character(0))
            out_text <- tryCatch(readLines(done$proc$get_output_file(),
                                           warn = FALSE),
                                 error = function(e) character(0))
            crash_log <- file.path(
              output_dir(config),
              sprintf("crash_%s_sc%d_cell%04d_rep%02d.log",
                      label, scenario, done$cell_id, done$rep_id))
            cat(sprintf("=== exit_code = %s ===\n", as.character(exit_code)),
                "=== stderr ===\n", err_text,
                "\n=== stdout ===\n", out_text,
                "\n", sep = "", file = crash_log)
            log_err(config, sprintf(
              "Subprocess FAILED [%s] sc=%d cell=%d rep=%d exit=%s -> %s",
              label, scenario, done$cell_id, done$rep_id,
              as.character(exit_code), crash_log))
            n_failed <- n_failed + 1L
          } else {
            n_completed <- n_completed + 1L
          }
          active_procs[[i]] <- NULL
        }

        if ((n_completed + n_failed) %% log_every == 0L ||
            (n_completed + n_failed) == total_work) {
          log_msg(sprintf("    [%s] progress: %d done (%d ok, %d failed) / %d (%.1f%%)",
                          label, n_completed + n_failed, n_completed, n_failed,
                          total_work,
                          100 * (n_completed + n_failed) / total_work))
        }
      }
    }
  }

  list(completed = n_completed, failed = n_failed)
}

# -----------------------------------------------------------------------------
# 5. Per-scenario aggregation
# -----------------------------------------------------------------------------
aggregate_scenario <- function(scenario, cells, config,
                               base_dir = output_dir(config)) {
  rds_dir   <- replicates_dir(config, scenario, base_dir)
  rds_files <- list.files(rds_dir,
                          pattern = sprintf("^sc%d_cell\\d+_rep\\d+\\.rds$",
                                            scenario),
                          full.names = TRUE)
  if (length(rds_files) == 0L) {
    log_msg(sprintf("  No replicate RDS files found in %s", rds_dir))
    return(data.frame())
  }

  rows <- vector("list", length(rds_files))
  for (i in seq_along(rds_files)) {
    r <- tryCatch(readRDS(rds_files[i]), error = function(e) NULL)
    if (is.null(r)) {
      log_err(config, sprintf("Failed to read %s; skipping in aggregation.",
                              rds_files[i]))
      next
    }
    rows[[i]] <- data.frame(
      scenario        = r$scenario,
      cell_id         = r$cell_id,
      rep_id          = r$rep_id,
      band            = r$band,
      N_ref           = r$N_ref,
      n_V1            = r$n_V1,
      n_V2            = r$n_V2,
      n_V3            = r$n_V3,
      ratio_V1        = r$ratio_V1,
      ratio_V2        = r$ratio_V2,
      ratio_V3        = r$ratio_V3,
      ceiling         = r$ceiling,
      p_del           = r$p_del,
      mutation_count  = r$mutation_count,
      n_emerge_needed = r$n_emerge_needed,
      detected        = r$detected,
      hit_modelName   = r$hit_modelName,
      hit_G           = r$hit_G,
      walltime_s      = r$walltime_s,
      seed            = r$seed,
      stringsAsFactors = FALSE
    )
  }
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
# 6. Main per-scenario loop
# -----------------------------------------------------------------------------
for (scenario in 1:4) {
  log_msg(sprintf("=== Scenario %d ===", scenario))

  cells <- readRDS(cells_path(config, scenario))
  log_msg(sprintf("  %d cells, %d reps/cell, %d total replicates",
                  nrow(cells), config$N_REPS_PER_CELL,
                  nrow(cells) * config$N_REPS_PER_CELL))

  # Build work queue and filter out completed replicates
  work_queue <- expand.grid(
    cell_id = cells$cell_id,
    rep_id  = seq_len(config$N_REPS_PER_CELL),
    KEEP.OUT.ATTRS = FALSE
  )

  exists_vec <- file.exists(file.path(
    replicates_dir(config, scenario),
    sprintf("sc%d_cell%04d_rep%02d.rds",
            scenario, work_queue$cell_id, work_queue$rep_id)
  ))
  work_queue <- work_queue[!exists_vec, , drop = FALSE]

  if (nrow(work_queue) == 0L) {
    log_msg("  All replicates already complete; aggregating only.")
  } else {
    log_msg(sprintf("  %d replicates pending (resume mode)",
                    nrow(work_queue)))
    res <- execute_work_queue(
      scenario   = scenario,
      work_queue = work_queue,
      cells      = cells,
      config     = config,
      label      = sprintf("sc%d", scenario)
    )
    log_msg(sprintf("  Scenario %d: %d completed, %d failed.",
                    scenario, res$completed, res$failed))
  }

  # Aggregate to summary
  log_msg("  Aggregating per-replicate RDS files...")
  summary_df <- aggregate_scenario(scenario, cells, config)
  if (nrow(summary_df) > 0L) {
    save_rds_atomic(summary_df, summary_path(config, scenario))
    log_msg(sprintf("  Summary: %s (%d rows)",
                    summary_path(config, scenario), nrow(summary_df)))
  } else {
    log_err(config, sprintf("Scenario %d produced no aggregable replicates.",
                            scenario))
  }
}

# -----------------------------------------------------------------------------
# 7. Optional robustness sweep over p_del values
# -----------------------------------------------------------------------------
if (isTRUE(config$RUN_ROBUSTNESS_SWEEP)) {
  log_msg("=== Robustness sweep over p_del ===")

  rob_cells_df <- stratified_robustness_cells(config)
  log_msg(sprintf("Robustness cells: %d total across scenarios/bands.",
                  nrow(rob_cells_df)))

  for (p_del_value in config$P_DEL_ROBUSTNESS_VALUES) {
    log_msg(sprintf("  --- p_del = %g ---", p_del_value))
    rob_base <- robustness_output_dir(config, p_del_value)
    dir.create(rob_base, recursive = TRUE, showWarnings = FALSE)

    for (scenario in unique(rob_cells_df$scenario)) {
      cells <- readRDS(cells_path(config, scenario))
      sc_cells <- rob_cells_df[rob_cells_df$scenario == scenario, ,
                               drop = FALSE]
      if (nrow(sc_cells) == 0L) next

      work_queue <- expand.grid(
        cell_id = sc_cells$cell_id,
        rep_id  = seq_len(config$N_REPS_PER_CELL),
        KEEP.OUT.ATTRS = FALSE
      )
      exists_vec <- file.exists(replicate_path(
        config, scenario, work_queue$cell_id, work_queue$rep_id,
        base_dir = rob_base))
      work_queue <- work_queue[!exists_vec, , drop = FALSE]

      if (nrow(work_queue) == 0L) {
        log_msg(sprintf("    sc%d: all robustness replicates complete.",
                        scenario))
        next
      }
      log_msg(sprintf("    sc%d: %d robustness replicates pending.",
                      scenario, nrow(work_queue)))
      res <- execute_work_queue(
        scenario        = scenario,
        work_queue      = work_queue,
        cells           = cells,
        config          = config,
        p_del_override  = p_del_value,
        output_base_dir = rob_base,
        label           = sprintf("rob_p%g_sc%d", p_del_value, scenario)
      )
      log_msg(sprintf("    sc%d: %d completed, %d failed.",
                      scenario, res$completed, res$failed))

      summary_df <- aggregate_scenario(scenario, cells, config,
                                       base_dir = rob_base)
      if (nrow(summary_df) > 0L) {
        save_rds_atomic(summary_df,
                        file.path(rob_base,
                                  sprintf("summary_sc%d.rds", scenario)))
      }
    }
  }
}

log_msg("All scenarios complete.")

# ============================================================================
# Final-results summary block
# ============================================================================
# Reads summary_sc<N>.rds for every scenario that completed, prints a
# band-by-detection-status table per scenario plus an overall version, and
# adds median/mean of n_emerge_needed per scenario and overall. Designed to
# run silently if any scenario's summary is missing.

.print_detection_table <- function(df, label) {
  if (is.null(df) || nrow(df) == 0L) {
    cat(sprintf("\n=== %s ===\n  (no replicates)\n", label))
    return(invisible(NULL))
  }
  df$status <- ifelse(is.na(df$n_emerge_needed), "NA", "detected")
  tbl <- table(df$band, df$status, useNA = "no")
  
  # Per-band median and mean (over detected replicates only)
  detected <- df[!is.na(df$n_emerge_needed), , drop = FALSE]
  if (nrow(detected) == 0L) {
    med <- mean_ <- rep(NA_real_, nrow(tbl))
    names(med) <- names(mean_) <- rownames(tbl)
  } else {
    med   <- tapply(detected$n_emerge_needed, detected$band,
                    stats::median, na.rm = TRUE)
    mean_ <- tapply(detected$n_emerge_needed, detected$band,
                    mean, na.rm = TRUE)
    # Align with table row order (some bands may have no detected rows)
    med   <- med[rownames(tbl)]
    mean_ <- mean_[rownames(tbl)]
  }
  
  cat(sprintf("\n=== %s ===\n", label))
  print(tbl)
  cat(sprintf("\n%-10s  %12s  %12s\n",
              "band", "median(det)", "mean(det)"))
  for (b in rownames(tbl)) {
    cat(sprintf("%-10s  %12s  %12s\n",
                b,
                if (is.na(med[b]))   "NA" else format(med[b],   digits = 4),
                if (is.na(mean_[b])) "NA" else format(mean_[b], digits = 4)))
  }
  
  # Overall (all bands collapsed)
  total_med  <- if (nrow(detected) > 0L)
    stats::median(detected$n_emerge_needed, na.rm = TRUE) else NA_real_
  total_mean <- if (nrow(detected) > 0L)
    mean(detected$n_emerge_needed, na.rm = TRUE) else NA_real_
  cat(sprintf("%-10s  %12s  %12s\n",
              "(all)",
              if (is.na(total_med))  "NA" else format(total_med,  digits = 4),
              if (is.na(total_mean)) "NA" else format(total_mean, digits = 4)))
  cat(sprintf("\nDetection rate: %d / %d (%.1f%%)\n",
              nrow(detected), nrow(df),
              100 * nrow(detected) / nrow(df)))
  invisible(NULL)
}

log_msg("=== Final detection summary ===")

scenario_summaries <- vector("list", 4L)
for (s in 1:4) {
  p <- summary_path(config, s)
  if (file.exists(p)) {
    scenario_summaries[[s]] <- readRDS(p)
    .print_detection_table(scenario_summaries[[s]],
                           sprintf("Scenario %d", s))
  } else {
    log_msg(sprintf("  Sc%d summary missing: %s", s, p))
  }
}

# Overall across all scenarios
non_null <- scenario_summaries[!vapply(scenario_summaries, is.null,
                                       logical(1L))]
if (length(non_null) > 0L) {
  overall_df <- do.call(rbind, lapply(non_null, function(d) {
    # Restrict to columns present in every summary to avoid rbind errors
    common <- Reduce(intersect,
                     lapply(non_null, colnames))
    d[, common, drop = FALSE]
  }))
  .print_detection_table(overall_df, "Overall (all scenarios)")
}
