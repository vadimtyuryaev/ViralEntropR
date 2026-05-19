# =============================================================================
# setup.R
# =============================================================================
#
# Libraries, configuration constants, paths, and helper functions used by
# every component of the Sample-Size Simulation Study.
#
# Sourced once at the top of each callr subprocess and once by the
# orchestrator simulation_study.R. Does not depend on any other file in the
# study folder. Loading is idempotent; sourcing twice is safe.
#
# Author : Vadim Tyuryaev
# =============================================================================

suppressPackageStartupMessages({
  library(ViralEntropR)
})

# -----------------------------------------------------------------------------
# resolve_study_dir
# -----------------------------------------------------------------------------
# Returns the absolute path to the Sample_Size_Simulation_Study folder.
# Honours the VIRAL_STUDY_DIR environment variable if set; otherwise tries
# the user's standard Windows and Linux paths in order.
# -----------------------------------------------------------------------------
resolve_study_dir <- function() {
  candidates <- c(
    Sys.getenv("VIRAL_STUDY_DIR", unset = ""),
    "C:/YORK_PhD/RESEARCH/PAPERS/GitHub/ViralEntropR/analysis/Sample_Size_Simulation_Study",
    "/home/vadimtyu/Clean_Code_Running/Dissertation/Simulation_Study",
    file.path(getwd(), "analysis", "Sample_Size_Simulation_Study"),
    getwd()
  )
  for (p in candidates) {
    if (nzchar(p) && dir.exists(p)) {
      return(normalizePath(p, mustWork = FALSE))
    }
  }
  stop("Cannot locate Sample_Size_Simulation_Study directory. ",
       "Set the VIRAL_STUDY_DIR environment variable to the absolute path.",
       call. = FALSE)
}

# -----------------------------------------------------------------------------
# auto_detect_n_workers
# -----------------------------------------------------------------------------
# Heuristic worker-count chooser. Used by build_config() when N_WORKERS is
# NULL and by the orchestrator when the env override is not set.
#
# Returns 1L when CPU or RAM headroom would force single-process operation.
# A return value of 1L is honoured by the orchestrator as sequential mode.
# -----------------------------------------------------------------------------
auto_detect_n_workers <- function(per_subproc_gb = 2,
                                  headroom_gb    = 16,
                                  cpu_reserve    = 4L,
                                  hard_cap       = 16L) {
  total_mem_gb <- tryCatch({
    if (.Platform$OS.type == "unix" && file.exists("/proc/meminfo")) {
      mem_kb <- as.numeric(gsub("\\D", "",
                                grep("^MemTotal:", readLines("/proc/meminfo"),
                                     value = TRUE)[1L]))
      mem_kb / 1024^2
    } else {
      tryCatch(utils::memory.limit() / 1024, error = function(e) 16)
    }
  }, error = function(e) 16)

  ram_cap <- max(1L, floor((total_mem_gb - headroom_gb) / per_subproc_gb))
  cpu_cap <- max(1L, parallel::detectCores() - cpu_reserve)
  as.integer(max(1L, min(ram_cap, cpu_cap, hard_cap)))
}

# -----------------------------------------------------------------------------
# build_config
# -----------------------------------------------------------------------------
# Returns the configuration list. Every quantity that varies between
# development and production lives here. The orchestrator may overwrite
# STUDY_DIR, REF_SEQ_FASTA, and N_WORKERS after calling this function.
#
# Concurrency switch:
#   N_WORKERS = 1L    -> sequential in-process execution
#   N_WORKERS > 1L    -> parallel callr::r_bg pool with N_WORKERS slots
#   N_WORKERS = NULL  -> auto-detected by the orchestrator from RAM and CPU
#
# The orchestrator additionally honours the env var N_WORKERS (integer)
# which overrides whatever build_config() returns. The env var RUN_ROBUSTNESS_SWEEP
# enables the optional p_del sensitivity sweep when set to "1" or "true".
# -----------------------------------------------------------------------------
build_config <- function() {
  list(
    # Reproducibility
    BASE_SEED         = 2025L,

    # Design grid
    N_REF_GRID_SIZE   = 10L,    
    RATIOS            = c(1.25, 1.5, 1.75, 2.5, 5),
    N_REPS_PER_CELL   = 30L,
    BAND_RANGES       = list(
      small  = c(10L,    100L),
      medium = c(100L,   1000L),
      large  = c(1000L,  10000L)
    ),

    # Sweep: 50-point log-spaced grid in every band; no refinement
    GRID_SIZE         = 50L,

    # GMM. NULL = let mclust evaluate both "E" and "V" via BIC.
    MCLUST_MODELS     = NULL,
    MCLUST_G          = 1:15,

    # Deleterious-mutation noise.
    # P_DEL is the per-(sequence, position) Bernoulli probability that a
    # given residue is replaced with a uniformly drawn non-reference standard
    # amino acid. Empirically anchored at 9.6e-4 via the 0.1%-frequency
    # threshold analysis on 21 SARS-CoV-2 variant-defining Spike positions;
    # rounded upward to 1e-3 as a conservative default.
    P_DEL             = 1e-3,

    # Robustness sweep over p_del values (optional).
    RUN_ROBUSTNESS_SWEEP    = TRUE,
    P_DEL_ROBUSTNESS_VALUES = c(0, 1e-4, 1e-3, 1e-2),
    ROBUSTNESS_CELLS_PER_BAND = 5L,

    # Reference sequence FASTA (set by orchestrator at startup)
    REF_SEQ_FASTA     = NULL,

    # Output (STUDY_DIR set by orchestrator at startup)
    STUDY_DIR         = NULL,
    OUTPUT_SUBDIR     = "outputs",
    ROBUSTNESS_SUBDIR = "robustness_pdel",

    # Concurrency
    N_WORKERS         = 24,   # NULL = auto-detect at orchestrator startup
    SUBPROCESS_POLL_S = 0.2   # subprocess completion poll interval
  )
}

# -----------------------------------------------------------------------------
# Amino-acid integer encoding
# -----------------------------------------------------------------------------
# Matches ViralEntropR::encode_aa_sequence's 25-symbol alphabet exactly.
# We only need the 20 standard residues (codes 1..20) for the simulator,
# since substitute amino acids are always standard.
# -----------------------------------------------------------------------------
.AA_TO_INT <- c(
  A = 1L,  R = 2L,  N = 3L,  D = 4L,  C = 5L,
  Q = 6L,  E = 7L,  G = 8L,  H = 9L,  I = 10L,
  L = 11L, K = 12L, M = 13L, F = 14L, P = 15L,
  S = 16L, T = 17L, W = 18L, Y = 19L, V = 20L
)

# Substitute lookup matrix: row r contains the 19 standard codes NOT equal to r.
# Used to draw substitute amino acids at mutated positions, and to draw
# replacement residues for deleterious-noise events, without an inner
# setdiff() call per position. Construction is one-time at file load.
.SUBSTITUTE_MATRIX <- local({
  m <- matrix(0L, nrow = 20L, ncol = 19L)
  for (r in seq_len(20L)) {
    m[r, ] <- setdiff(seq_len(20L), r)
  }
  m
})

# -----------------------------------------------------------------------------
# build_empirical_constants
# -----------------------------------------------------------------------------
# Computes the empirical sets and pools used by the simulator, from the
# sarscov2_variants catalogue bundled with the ViralEntropR package.
# Returns a named list. Called once per subprocess at startup.
# -----------------------------------------------------------------------------
build_empirical_constants <- function() {
  e <- new.env()
  data("sarscov2_variants", package = "ViralEntropR", envir = e)
  sv <- e$sarscov2_variants

  labels    <- unlist(sv$WHO_Label)
  mut_sites <- sv$Mutation_Sites

  mut_counts <- vapply(mut_sites, length, integer(1L))
  names(mut_counts) <- labels

  pool_11 <- sort(unique(unlist(mut_sites[labels != "Omicron"])))
  pool_12 <- sort(unique(unlist(mut_sites)))

  if (!(614L %in% pool_11) || !(614L %in% pool_12)) {
    stop("Empirical invariant violated: position 614 missing from POOL_11/POOL_12. ",
         "Check sarscov2_variants$Mutation_Sites integrity.", call. = FALSE)
  }

  list(
    MUT_COUNT_SET_11  = as.integer(sort(unique(mut_counts[labels != "Omicron"]))),
    MUT_COUNT_OMICRON = as.integer(unname(mut_counts["Omicron"])),
    POOL_11           = as.integer(pool_11),
    POOL_12           = as.integer(pool_12),
    POOL_11_DRAW      = as.integer(setdiff(pool_11, 614L)),
    POOL_12_DRAW      = as.integer(setdiff(pool_12, 614L))
  )
}

# -----------------------------------------------------------------------------
# RNG seed scheme
# -----------------------------------------------------------------------------
# Each (scenario, cell_id, rep_id) maps to a unique integer seed. With
# scenario in 1..4, cell_id up to 3,750 (Sc 3/4), and rep_id up to 30,
# the value stays well within 32-bit integer range.
# -----------------------------------------------------------------------------
seed_for_replicate <- function(base_seed, scenario, cell_id, rep_id) {
  as.integer(base_seed + scenario * 100000000L + cell_id * 1000L + rep_id)
}

# -----------------------------------------------------------------------------
# Path helpers
# -----------------------------------------------------------------------------
output_dir <- function(config) {
  file.path(config$STUDY_DIR, config$OUTPUT_SUBDIR)
}

robustness_output_dir <- function(config, p_del_value) {
  file.path(output_dir(config),
            config$ROBUSTNESS_SUBDIR,
            sprintf("p_del_%g", p_del_value))
}

replicates_dir <- function(config, scenario, base_dir = output_dir(config)) {
  file.path(base_dir, sprintf("replicates_sc%d", scenario))
}

cells_path <- function(config, scenario) {
  file.path(output_dir(config), sprintf("cells_sc%d.rds", scenario))
}

# Note: 4-digit cell id to accommodate up to 3,750 cells in scenarios 3 and 4.
replicate_path <- function(config, scenario, cell_id, rep_id,
                           base_dir = output_dir(config)) {
  file.path(replicates_dir(config, scenario, base_dir),
            sprintf("sc%d_cell%04d_rep%02d.rds",
                    scenario, cell_id, rep_id))
}

summary_path <- function(config, scenario) {
  file.path(output_dir(config), sprintf("summary_sc%d.rds", scenario))
}

error_log_path <- function(config) {
  file.path(output_dir(config), "error_log.txt")
}

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------
# Timestamped message emitter. Goes to stderr (R's message channel) so that
# stdout stays clean for tabular output if needed.
# -----------------------------------------------------------------------------
log_msg <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", ts, paste(..., collapse = " ")))
}

log_err <- function(config, ...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", ts, paste(..., collapse = " "))
  message(line)
  if (!is.null(config$STUDY_DIR)) {
    cat(line, "\n", file = error_log_path(config), append = TRUE)
  }
}

# -----------------------------------------------------------------------------
# Atomic RDS write
# -----------------------------------------------------------------------------
# Writes to a .tmp file then renames. Guarantees the destination file
# either does not exist or contains a fully-written object. Avoids the
# resume-from-corrupt-RDS edge case where a crashed subprocess leaves a
# partial file at the destination path.
# -----------------------------------------------------------------------------
save_rds_atomic <- function(object, path) {
  tmp_path <- paste0(path, ".tmp")
  saveRDS(object, file = tmp_path)
  file.rename(tmp_path, path)
  invisible(path)
}
