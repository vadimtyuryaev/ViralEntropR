# =============================================================================
# ViralEntropR — Simulation Study
# =============================================================================
#
# Purpose
# -------
# Benchmarks the entropy-GMM-PAM detection pipeline across biological
# scenarios by sweeping over n_new_mut_seq_vec and recording site-detection
# and change-point accuracy metrics.  Three detection criteria are evaluated
# in parallel: strict on-time, any-time eventual, and partial (threshold).
#
# Phased execution
# ----------------
#   Phase 1  Single-variant scenarios (four noise conditions).
#            Intended for dissertation figures on minimal detectable
#            sample size under varying background noise.
#   Phase 2  Two-variant scenarios (timing and rate conditions).
#            Tests pipeline sensitivity to competitive variant dynamics.
#   Phase 3  Full sweep (3–10 variants, all noise conditions).
#            Production-scale benchmark; requires HPC Linux.
#
# Parallelism
# -----------
#   Windows (local dev) : future::multisession  (background R sessions)
#   Linux HPC           : future::multicore     (fork-based, zero overhead)
#   Switching is automatic via .Platform$OS.type.
#
# Output
# ------
#   Per-iteration : RDS  (results + params + metrics) + optional HTML table
#   Summary RDS   : FD/CP matrices, detection-rate arrays, scenario metadata
#   Plots         : PDF (vector) + PNG (300 dpi) per figure
#
# Author:  Vadim Tyuryaev, York University
# =============================================================================


# ===========================================================================
# 0.  LIBRARIES
# ===========================================================================

required_pkgs <- c("ViralEntropR", "future", "future.apply",
                   "ggplot2", "tidyr", "dplyr", "scales",
                   "kableExtra", "lubridate", "tictoc")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1L), quietly = TRUE)]
if (length(missing_pkgs) > 0L)
  stop("Install missing packages before running:\n  ",
       paste(missing_pkgs, collapse = ", "))

library(ViralEntropR)
library(future)
library(future.apply)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(kableExtra)
library(lubridate)
library(tictoc)


# ===========================================================================
# 1.  CONTROL PANEL  — edit only here
# ===========================================================================

# ── Phase selection ──────────────────────────────────────────────────────────
# 1 = single-variant | 2 = two-variant | 3 = full sweep (3–10 variants)
PHASE <- 1L

# ── Iteration count ─────────────────────────────────────────────────────────
# Development / quick test : 10–20
# Dissertation figures      : 100–250
# Production benchmark      : 500+
N_ITER    <- 50L
BASE_SEED <- 2025L   # scenario-offset added per scenario for reproducibility

# ── Reference sequence ──────────────────────────────────────────────────────
# Preliminary phases: 94-aa synthetic sequence (fast).
# Full-scale runs: replace with SARS-CoV-2 Spike reference, e.g.:
#   spike_fasta <- readLines("data-raw/spike_ref.fasta")
#   REF_SEQ     <- paste(spike_fasta[!grepl("^>", spike_fasta)], collapse = "")
REF_SEQ <- paste0(
  "MKTIIALSYIFCLVFADYKDDDDKKDDDDKYDAFVLCFIYSLAIITKM",
  "MKTIIALSYIFCLVFADYKDDDDKKDDDDKYDAFVLCFIYSLAIITKM"
)

# ── De-novo sequence count sweep ─────────────────────────────────────────────
# Coarser grid (e.g. seq(2, 51, 5)) for development; 2:51 for publication.
# N_NEW_MUT_SEQ_VEC <- 2:51

N_NEW_MUT_SEQ_VEC <- seq(2, 51, 5)

# ── Simulation timeline ──────────────────────────────────────────────────────
# END_DATE must be set one month beyond the intended last simulation month.
# Phase 1/2 : "2020-01-01" → "2022-01-01"  (24-month simulation, 12 windows)
# Phase 3   : "2020-01-01" → "2023-01-01"  (36-month simulation, 18 windows)
#             Needed so 10-variant configurations fit within MAX_EMERG_ALLOWED.
START_DATE <- "2020-01-01"
END_DATE   <- "2022-01-01"

# ── Sequence generation ──────────────────────────────────────────────────────
N_SEQ_PER_MONTH <- 1000L
MUT_RATE_VAR    <- 0        # 0 = deterministic growth in simulate_variant_evolution
REF_VARIABILITY <- FALSE

# ── Partition settings ───────────────────────────────────────────────────────
SLIDING_WINDOW_LEN <- 2L
WINDOW_OPTION      <- 3L    # 3 = non-overlapping (jumping windows)

# ── Detection criteria ───────────────────────────────────────────────────────
# FD_strict  : First_Detect_Sample from run_detection_study() — all sites in
#              top cluster at the expected partition; site sets mutually distinct.
# FD_anytime : First n where all variants are eventually detected (dp != "-").
# FD_partial : First n where every variant has >= DETECTION_THRESHOLD fraction
#              of its sites in the top cluster at the expected partition.
DETECTION_THRESHOLD <- 0.75   # used only for FD_partial

# ── Output paths ─────────────────────────────────────────────────────────────
# Set both; the correct one is selected automatically by OS.
OUT_DIR_WIN   <- "C:/YORK_PhD/RESEARCH/PAPERS/GitHub/ViralEntropR/output/simulation_study"
OUT_DIR_LINUX <- "/home/vadimtyu/Simulate_AA_mutations"
OUT_DIR <- if (.Platform$OS.type == "windows") OUT_DIR_WIN else OUT_DIR_LINUX

# ── Per-iteration HTML tables ─────────────────────────────────────────────────
# Produces one .html file per (scenario × iteration).  Leave FALSE during
# development — each file is ~100 KB; 2000 files = ~200 MB.
SAVE_ITER_HTML <- FALSE

# ── Parallelism ───────────────────────────────────────────────────────────────
# Capped at detectCores() - 1 to leave one core for OS and I/O.
# Override manually for HPC (e.g. 32L or 64L).
N_CORES <- max(1L, parallel::detectCores() - 1L)

# ===========================================================================
# 2.  DERIVED CONSTANTS  — nothing below this block is hardcoded
# ===========================================================================

# Simulation time axis — matches run_detection_study() internals
SIM_DATES  <- seq.Date(as.Date(START_DATE), as.Date(END_DATE), by = "month")
N_PERIODS  <- length(SIM_DATES)

# partition_time_windows uses lubridate::interval %/% months(1) for n_chunks,
# which gives one fewer month than length(seq.Date(...)).
# MAX_EMERG_ALLOWED is therefore derived from the same formula to ensure every
# variant's expected detection partition (ceiling(em / wl)) is strictly within
# the available partition range AND has at least one subsequent partition for
# the forward search.
TOTAL_MONTHS_PTW <- as.integer(
  lubridate::interval(as.Date(START_DATE), as.Date(END_DATE)) %/%
    lubridate::period(1, "months")
)
N_CHUNKS_PTW     <- TOTAL_MONTHS_PTW %/% SLIDING_WINDOW_LEN
MAX_EMERG_ALLOWED <- N_CHUNKS_PTW * SLIDING_WINDOW_LEN - 1L

stopifnot(
  "END_DATE too close to START_DATE: need at least 2 non-overlapping windows" =
    N_CHUNKS_PTW >= 2L,
  "MAX_EMERG_ALLOWED must be positive" =
    MAX_EMERG_ALLOWED >= 1L
)

N_NEW_LEN <- length(N_NEW_MUT_SEQ_VEC)

message(sprintf(
  "Derived constants: N_PERIODS=%d | N_CHUNKS_PTW=%d | MAX_EMERG_ALLOWED=%d",
  N_PERIODS, N_CHUNKS_PTW, MAX_EMERG_ALLOWED
))


# ===========================================================================
# 3.  SCENARIO DEFINITIONS
# ===========================================================================
#
# Each scenario is a named list consumed by sample_run_params().
#
# Fields
# ------
#   name             Character.  Scenario identifier (used in plots and filenames).
#   nv_range         Integer vector.  Range from which n_variants is drawn.
#   mut_first_range  Integer vector.  Mutation count for Variant 1.
#   mut_second_range Integer vector.  Mutation count for Variant 2.
#   mut_later_range  Integer vector.  Mutation counts for Variants 3+.
#   rate_range_v1    Numeric vector.  Growth rate pool for Variant 1 (and all
#                    variants when rate_range_later is NULL).
#   rate_range_later Numeric vector or NULL.  When non-NULL, Variants 2+ draw
#                    from this pool independently (differential competition).
#   ref_months_range Integer vector.  Reference phase duration (months).
#   vi_range         Integer vector.  Pool for inter-variant interval draws.
#   prob_delet_range Numeric or vector.  Scalar 0 for none; vector to sample from.
#   drate_range      Numeric or vector.  Deleterious rate pool.
#   n_delet_range    Integer or vector.  Deletion count pool.
#
# Phase 1 — single variant; four noise conditions
PHASE1_SCENARIOS <- list(

  list(name             = "NoDel_NoSiteDel",
       nv_range         = 1L,
       mut_first_range  = 1:4,
       mut_second_range = 4:8,
       mut_later_range  = 8:33,
       rate_range_v1    = seq(1.125, 2.000, 0.125),
       rate_range_later = NULL,
       ref_months_range = 1:4,
       vi_range         = 2:12,
       prob_delet_range = 0,
       drate_range      = 0,
       n_delet_range    = 0L),

  list(name             = "Del_NoSiteDel",
       nv_range         = 1L,
       mut_first_range  = 1:4,
       mut_second_range = 4:8,
       mut_later_range  = 8:33,
       rate_range_v1    = seq(1.125, 2.000, 0.125),
       rate_range_later = NULL,
       ref_months_range = 1:4,
       vi_range         = 2:12,
       prob_delet_range = seq(0.05, 0.10, 0.01),
       drate_range      = 0,
       n_delet_range    = 1:9),

  list(name             = "NoDel_SiteDel",
       nv_range         = 1L,
       mut_first_range  = 1:4,
       mut_second_range = 4:8,
       mut_later_range  = 8:33,
       rate_range_v1    = seq(1.125, 2.000, 0.125),
       rate_range_later = NULL,
       ref_months_range = 1:4,
       vi_range         = 2:12,
       prob_delet_range = 0,
       drate_range      = seq(0.05, 0.10, 0.01),
       n_delet_range    = 1:9),

  list(name             = "Del_And_SiteDel",
       nv_range         = 1L,
       mut_first_range  = 1:4,
       mut_second_range = 4:8,
       mut_later_range  = 8:33,
       rate_range_v1    = seq(1.125, 2.000, 0.125),
       rate_range_later = NULL,
       ref_months_range = 1:4,
       vi_range         = 2:12,
       prob_delet_range = seq(0.05, 0.10, 0.01),
       drate_range      = seq(0.05, 0.10, 0.01),
       n_delet_range    = 1:9)
)

# Phase 2 — two variants; timing and rate conditions
PHASE2_SCENARIOS <- list(

  list(name             = "SameRate_ShortInterval",
       nv_range         = 2L,
       mut_first_range  = 1:4,
       mut_second_range = 4:8,
       mut_later_range  = 8:33,
       rate_range_v1    = seq(1.500, 2.000, 0.125),
       rate_range_later = NULL,
       ref_months_range = 1:3,
       vi_range         = 2:4,
       prob_delet_range = 0,
       drate_range      = 0,
       n_delet_range    = 0L),

  list(name             = "SameRate_LongInterval",
       nv_range         = 2L,
       mut_first_range  = 1:4,
       mut_second_range = 4:8,
       mut_later_range  = 8:33,
       rate_range_v1    = seq(1.500, 2.000, 0.125),
       rate_range_later = NULL,
       ref_months_range = 1:3,
       vi_range         = 8:12,
       prob_delet_range = 0,
       drate_range      = 0,
       n_delet_range    = 0L),

  list(name             = "DiffRate_Compete",
       nv_range         = 2L,
       mut_first_range  = 1:4,
       mut_second_range = 4:8,
       mut_later_range  = 8:33,
       rate_range_v1    = seq(1.125, 1.500, 0.125),   # V1: slow competitor
       rate_range_later = seq(1.625, 2.000, 0.125),   # V2+: fast competitor
       ref_months_range = 1:3,
       vi_range         = 3:6,
       prob_delet_range = 0,
       drate_range      = 0,
       n_delet_range    = 0L),

  list(name             = "DiffRate_WithNoise",
       nv_range         = 2L,
       mut_first_range  = 1:4,
       mut_second_range = 4:8,
       mut_later_range  = 8:33,
       rate_range_v1    = seq(1.125, 1.500, 0.125),
       rate_range_later = seq(1.625, 2.000, 0.125),
       ref_months_range = 1:3,
       vi_range         = 3:6,
       prob_delet_range = seq(0.05, 0.10, 0.01),
       drate_range      = seq(0.05, 0.10, 0.01),
       n_delet_range    = 1:9)
)

# Phase 3 — full sweep: inherit Phase 1 noise structure, expand nv_range.
# Note: END_DATE = "2023-01-01" strongly recommended for Phase 3 to give
# sufficient temporal range for up to 10 variants with 2–12 month intervals.
PHASE3_SCENARIOS <- lapply(PHASE1_SCENARIOS, function(sc) {
  sc$nv_range <- 3:10
  sc
})

SCENARIOS <- switch(
  as.character(PHASE),
  "1" = PHASE1_SCENARIOS,
  "2" = PHASE2_SCENARIOS,
  "3" = PHASE3_SCENARIOS,
  stop(sprintf("PHASE must be 1, 2, or 3. Received: %s", PHASE))
)

N_SCEN   <- length(SCENARIOS)
SC_NAMES <- vapply(SCENARIOS, `[[`, character(1L), "name")


# ===========================================================================
# 4.  OUTPUT INFRASTRUCTURE
# ===========================================================================

timestamp  <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_SUBDIR <- file.path(OUT_DIR, sprintf("Phase%d_%s", PHASE, timestamp))
dir.create(OUT_SUBDIR, recursive = TRUE, showWarnings = FALSE)

error_log_path <- file.path(OUT_SUBDIR,
                             sprintf("error_log_Phase%d.txt", PHASE))
cat(sprintf(
  "ViralEntropR Simulation Study\nPhase:   %d\nStarted: %s\nN_ITER:  %d\nN_CORES: %d\nREF_SEQ length: %d aa\nScenarios: %s\n\n",
  PHASE, timestamp, N_ITER, N_CORES, nchar(REF_SEQ),
  paste(SC_NAMES, collapse = ", ")
), file = error_log_path)


# ===========================================================================
# 5.  PARALLEL BACKEND
# ===========================================================================

if (.Platform$OS.type == "windows") {
  future::plan(future::multisession, workers = N_CORES)
  message(sprintf("Windows: multisession plan, %d workers.", N_CORES))
} else {
  future::plan(future::multicore, workers = N_CORES)
  message(sprintf("Linux/Mac: multicore plan, %d workers.", N_CORES))
}


# ===========================================================================
# 6.  HELPER FUNCTIONS
# ===========================================================================

# ---------------------------------------------------------------------------
# sample_run_params
# ---------------------------------------------------------------------------
# Draws all stochastic parameters for one iteration of scenario sc.
# The rejection loop for variant intervals uses max_emerg_allowed (derived
# in Section 2), guaranteeing every variant's expected detection partition
# is within the available partition range.
# ---------------------------------------------------------------------------
sample_run_params <- function(sc, max_emerg_allowed) {

  nv <- if (length(sc$nv_range) == 1L) sc$nv_range
        else sample(sc$nv_range, 1L)

  # Mutation counts — biologically increasing across variant order
  mut_counts     <- integer(nv)
  mut_counts[1L] <- sample(sc$mut_first_range, 1L)
  if (nv >= 2L) mut_counts[2L]    <- sample(sc$mut_second_range, 1L)
  if (nv >= 3L) mut_counts[3L:nv] <- sample(sc$mut_later_range,
                                              nv - 2L, replace = TRUE)

  # Growth rates (per variant)
  if (!is.null(sc$rate_range_later) && nv > 1L) {
    rates <- c(sample(sc$rate_range_v1,    1L),
               sample(sc$rate_range_later, nv - 1L, replace = TRUE))
  } else {
    rates <- sample(sc$rate_range_v1, nv, replace = TRUE)
  }

  # Reference phase length
  num_months_ref <- sample(sc$ref_months_range, 1L)

  # Variant intervals — rejection-sampled so every variant's expected
  # detection partition falls within the available partition range.
  if (nv == 1L) {
    vi <- integer(0L)
    # Guard: single-variant case; adjust ref months if necessary
    if ((num_months_ref + 1L) > max_emerg_allowed)
      num_months_ref <- max(1L, max_emerg_allowed - 1L)

  } else {
    max_tries <- 10000L
    found     <- FALSE
    for (try_i in seq_len(max_tries)) {
      vi    <- sample(sc$vi_range, nv - 1L, replace = TRUE)
      emerg <- num_months_ref + 1L + cumsum(c(0L, vi))
      if (max(emerg) <= max_emerg_allowed) {
        found <- TRUE
        break
      }
      # After 200 failed draws reduce ref months to widen the window
      if (try_i == 200L && num_months_ref > 1L)
        num_months_ref <- 1L
    }
    if (!found)
      stop(sprintf(paste0(
        "sample_run_params: could not draw valid variant intervals after %d tries. ",
        "Extend END_DATE or reduce vi_range / nv_range for this PHASE."
      ), max_tries))
  }

  # Noise parameters
  prob_delet <- if (length(sc$prob_delet_range) == 1L) sc$prob_delet_range
                else sample(sc$prob_delet_range, 1L)
  d_rate     <- if (length(sc$drate_range) == 1L)      sc$drate_range
                else sample(sc$drate_range, 1L)
  n_delet    <- if (length(sc$n_delet_range)    == 1L)  sc$n_delet_range
                else sample(sc$n_delet_range, 1L)

  list(nv = nv, mut_counts = mut_counts, rates = rates,
       num_months_ref = num_months_ref, vi = vi,
       prob_delet = prob_delet, d_rate = d_rate, n_delet = n_delet)
}

# ---------------------------------------------------------------------------
# compute_fd_anytime
# ---------------------------------------------------------------------------
# First n_new_seq at which ALL variants are eventually detected
# (detection_partition != "-"), regardless of timing.
# ---------------------------------------------------------------------------
compute_fd_anytime <- function(results, n_new_mut_seq_vec) {
  for (n in n_new_mut_seq_vec) {
    sub <- results[results$num_denovo_sequences == n, ]
    if (nrow(sub) == 0L) next
    if (all(!is.na(sub$detection_partition) & sub$detection_partition != "-"))
      return(n)
  }
  NA_real_
}

# ---------------------------------------------------------------------------
# compute_fd_partial
# ---------------------------------------------------------------------------
# First n_new_seq at which every variant has >= threshold fraction of its
# mutation sites appearing in the top-entropy cluster at the expected
# partition.  detected_sites_raw always holds the top-cluster sites at
# part_em, making this metric independent of the forward-search result.
# ---------------------------------------------------------------------------
compute_fd_partial <- function(results, threshold, n_new_mut_seq_vec) {
  for (n in n_new_mut_seq_vec) {
    sub <- results[results$num_denovo_sequences == n, ]
    if (nrow(sub) == 0L) next
    fracs <- mapply(function(a_raw, d_raw) {
      a <- unlist(strsplit(a_raw, ","))
      d <- unlist(strsplit(d_raw, ","))
      if (length(a) == 0L) return(1.0)
      length(intersect(a, d)) / length(a)
    }, sub$sites_raw, sub$detected_sites_raw)
    if (all(fracs >= threshold, na.rm = TRUE)) return(n)
  }
  NA_real_
}

# ---------------------------------------------------------------------------
# compute_detect_curves
# ---------------------------------------------------------------------------
# Returns three logical/numeric vectors of length(n_new_mut_seq_vec):
#   strict    — all variants detected on time (dp == actual_partition)
#   anytime   — all variants eventually detected (dp != "-")
#   mean_frac — mean fraction of mutation sites in top cluster at part_em
# Used to build per-iteration detection curves; these are averaged across
# iterations in the aggregation step.
# ---------------------------------------------------------------------------
compute_detect_curves <- function(results, n_new_mut_seq_vec) {
  strict    <- logical(length(n_new_mut_seq_vec))
  anytime   <- logical(length(n_new_mut_seq_vec))
  mean_frac <- numeric(length(n_new_mut_seq_vec))

  for (k in seq_along(n_new_mut_seq_vec)) {
    n   <- n_new_mut_seq_vec[k]
    sub <- results[results$num_denovo_sequences == n, ]
    if (nrow(sub) == 0L) next

    # detection_partition is character when any value is "-";
    # convert to numeric for the strict comparison, treating "-" as NA.
    dp_num  <- suppressWarnings(as.numeric(sub$detection_partition))
    ap_num  <- sub$actual_partition

    strict[k]  <- all(!is.na(dp_num) & dp_num == ap_num)
    anytime[k] <- all(!is.na(sub$detection_partition) &
                        sub$detection_partition != "-")

    fracs <- mapply(function(a_raw, d_raw) {
      a <- unlist(strsplit(a_raw, ","))
      d <- unlist(strsplit(d_raw, ","))
      if (length(a) == 0L) return(1.0)
      length(intersect(a, d)) / length(a)
    }, sub$sites_raw, sub$detected_sites_raw)
    mean_frac[k] <- mean(fracs, na.rm = TRUE)
  }

  list(strict = strict, anytime = anytime, mean_frac = mean_frac)
}


# ===========================================================================
# 7.  WORKER FUNCTION
# ===========================================================================
# Executes one full run_detection_study() call and extracts all metrics.
# Wrapped in tryCatch so that a single failed iteration (GMM degeneracy,
# OOB partition, etc.) logs an error without crashing the outer loop.
#
# All UPPER_CASE objects are sourced from the parent (global) environment;
# future_lapply exports them automatically to each worker process.
# ===========================================================================

run_one_iter <- function(j, s_idx, params) {

  tryCatch({

    run_res <- run_detection_study(
      ref_seq               = REF_SEQ,
      variants_list         = list(params$mut_counts),
      mutation_rate_list    = list(params$rates),
      n_new_mut_seq_vec     = N_NEW_MUT_SEQ_VEC,
      num_months_ref_seq    = params$num_months_ref,
      variant_interval      = params$vi,
      start_date            = START_DATE,
      end_date              = END_DATE,
      mutation_rate_variability = MUT_RATE_VAR,
      deleterious_rate      = params$d_rate,
      ref_variability       = REF_VARIABILITY,
      n_seq_per_month       = N_SEQ_PER_MONTH,
      prob_delet            = params$prob_delet,
      n_delet               = params$n_delet,
      sliding_window_length = SLIDING_WINDOW_LEN,
      window_option         = WINDOW_OPTION,
      save_html             = FALSE
    )

    results    <- run_res$Results

    # Scalar metrics — run_detection_study is called with a single-element
    # variants_list, so First_Detect_Sample and CP_Accuracy are length-1.
    fd_strict  <- run_res$First_Detect_Sample[[1L]]
    cp         <- run_res$CP_Accuracy[[1L]]

    # Post-hoc detection metrics computed from the Results data frame
    fd_anytime <- compute_fd_anytime(results, N_NEW_MUT_SEQ_VEC)
    fd_partial <- compute_fd_partial(results, DETECTION_THRESHOLD,
                                     N_NEW_MUT_SEQ_VEC)
    dc         <- compute_detect_curves(results, N_NEW_MUT_SEQ_VEC)

    # Save per-iteration RDS (results + metrics; heavy Sim_List discarded)
    rds_path <- file.path(OUT_SUBDIR,
                          sprintf("s%02d_iter%03d.rds", s_idx, j))
    saveRDS(
      list(params     = params,
           results    = results,
           fd_strict  = fd_strict,
           fd_anytime = fd_anytime,
           fd_partial = fd_partial,
           cp         = cp,
           dc         = dc),
      file = rds_path
    )

    # Optionally save HTML table
    html_path <- NULL
    if (SAVE_ITER_HTML) {
      html_path <- sub("\\.rds$", ".html", rds_path)
      tryCatch(
        kableExtra::save_kable(run_res$Table, file = html_path),
        error = function(e) NULL
      )
    }

    list(fd_strict  = fd_strict,
         fd_anytime = fd_anytime,
         fd_partial = fd_partial,
         cp         = cp,
         dc         = dc,
         rds_path   = rds_path,
         html_path  = html_path,
         error      = NULL)

  }, error = function(e) {
    list(fd_strict  = NA_real_,
         fd_anytime = NA_real_,
         fd_partial = NA_real_,
         cp         = NA_real_,
         dc         = NULL,
         rds_path   = NULL,
         html_path  = NULL,
         error      = conditionMessage(e))
  })
}


# ===========================================================================
# 8.  PRE-ALLOCATE RESULT STRUCTURES
# ===========================================================================

# Scalar metric matrices  [iteration × scenario]
FD_strict_mat  <- matrix(NA_real_, nrow = N_ITER, ncol = N_SCEN,
                         dimnames = list(NULL, SC_NAMES))
FD_anytime_mat <- matrix(NA_real_, nrow = N_ITER, ncol = N_SCEN,
                         dimnames = list(NULL, SC_NAMES))
FD_partial_mat <- matrix(NA_real_, nrow = N_ITER, ncol = N_SCEN,
                         dimnames = list(NULL, SC_NAMES))
CP_mat         <- matrix(NA_real_, nrow = N_ITER, ncol = N_SCEN,
                         dimnames = list(NULL, SC_NAMES))

# Detection curve arrays  [n_idx × iteration × scenario]
detect_strict_arr  <- array(NA_real_,
                             dim      = c(N_NEW_LEN, N_ITER, N_SCEN),
                             dimnames = list(N_NEW_MUT_SEQ_VEC, NULL, SC_NAMES))
detect_anytime_arr <- array(NA_real_,
                             dim      = c(N_NEW_LEN, N_ITER, N_SCEN),
                             dimnames = list(N_NEW_MUT_SEQ_VEC, NULL, SC_NAMES))
mean_frac_arr      <- array(NA_real_,
                             dim      = c(N_NEW_LEN, N_ITER, N_SCEN),
                             dimnames = list(N_NEW_MUT_SEQ_VEC, NULL, SC_NAMES))

# Lightweight simulation index (parameters + file paths per run)
Simulation <- vector("list", N_SCEN)
for (.s in seq_len(N_SCEN)) Simulation[[.s]] <- vector("list", N_ITER)
rm(.s)


# ===========================================================================
# 9.  MAIN LOOP
# ===========================================================================

tic("Simulation study total")

for (s_idx in seq_len(N_SCEN)) {

  sc <- SCENARIOS[[s_idx]]
  message(sprintf("\n=== Scenario %d / %d : %s ===", s_idx, N_SCEN, sc$name))

  # Pre-draw all iteration parameters before parallel dispatch.
  # set.seed here makes parameter draws reproducible across re-runs and
  # independent of N_CORES (fork/multisession differences do not affect this).
  set.seed(BASE_SEED + s_idx * 1000L)
  iter_params <- lapply(seq_len(N_ITER),
                        function(j) sample_run_params(sc, MAX_EMERG_ALLOWED))

  # Parallel dispatch — future_lapply handles seed propagation to workers
  tic(sprintf("Scenario %d (%s)", s_idx, sc$name))
  iter_results <- future.apply::future_lapply(
    seq_len(N_ITER),
    function(j) run_one_iter(j, s_idx, iter_params[[j]]),
    future.seed = TRUE
  )
  toc()

  # Unpack and store results
  n_errors <- 0L
  for (j in seq_len(N_ITER)) {
    r <- iter_results[[j]]

    if (!is.null(r$error)) {
      n_errors <- n_errors + 1L
      cat(sprintf("[s%02d iter%03d] %s\n", s_idx, j, r$error),
          file = error_log_path, append = TRUE)
    }

    FD_strict_mat[j,  s_idx] <- r$fd_strict
    FD_anytime_mat[j, s_idx] <- r$fd_anytime
    FD_partial_mat[j, s_idx] <- r$fd_partial
    CP_mat[j,         s_idx] <- r$cp

    if (!is.null(r$dc)) {
      detect_strict_arr[,  j, s_idx] <- as.numeric(r$dc$strict)
      detect_anytime_arr[, j, s_idx] <- as.numeric(r$dc$anytime)
      mean_frac_arr[,      j, s_idx] <- r$dc$mean_frac
    }

    Simulation[[s_idx]][[j]] <- list(params    = iter_params[[j]],
                                     rds_path  = r$rds_path,
                                     html_path = r$html_path,
                                     error     = r$error)
  }

  message(sprintf("  Done: %d / %d successful, %d errors.",
                  N_ITER - n_errors, N_ITER, n_errors))
}

toc()


# ===========================================================================
# 10.  AGGREGATION AND SUMMARY
# ===========================================================================

# Detection rate curves averaged across iterations
detect_rate_df <- do.call(rbind, lapply(seq_len(N_SCEN), function(s) {
  data.frame(
    n         = as.integer(dimnames(detect_strict_arr)[[1L]]),
    strict    = rowMeans(detect_strict_arr[,,  s, drop = FALSE], na.rm = TRUE),
    anytime   = rowMeans(detect_anytime_arr[,, s, drop = FALSE], na.rm = TRUE),
    mean_frac = rowMeans(mean_frac_arr[,,      s, drop = FALSE], na.rm = TRUE),
    scenario  = SC_NAMES[s],
    stringsAsFactors = FALSE
  )
}))

# Per-scenario summary statistics
iqr_fn    <- function(x) diff(quantile(x, c(0.25, 0.75), na.rm = TRUE))
na_pct_fn <- function(x) mean(is.na(x)) * 100

summary_df <- data.frame(
  scenario           = SC_NAMES,

  mean_fd_strict     = colMeans(FD_strict_mat,  na.rm = TRUE),
  median_fd_strict   = apply(FD_strict_mat,  2, median, na.rm = TRUE),
  iqr_fd_strict      = apply(FD_strict_mat,  2, iqr_fn),
  pct_NA_fd_strict   = apply(FD_strict_mat,  2, na_pct_fn),

  mean_fd_anytime    = colMeans(FD_anytime_mat, na.rm = TRUE),
  median_fd_anytime  = apply(FD_anytime_mat, 2, median, na.rm = TRUE),
  iqr_fd_anytime     = apply(FD_anytime_mat, 2, iqr_fn),

  mean_fd_partial    = colMeans(FD_partial_mat, na.rm = TRUE),
  median_fd_partial  = apply(FD_partial_mat, 2, median, na.rm = TRUE),
  iqr_fd_partial     = apply(FD_partial_mat, 2, iqr_fn),

  mean_cp_accuracy   = colMeans(CP_mat, na.rm = TRUE),
  median_cp_accuracy = apply(CP_mat, 2, median, na.rm = TRUE),

  row.names = NULL,
  stringsAsFactors = FALSE
)
print(summary_df, digits = 3)

# Save full summary RDS
summary_path <- file.path(OUT_SUBDIR,
                           sprintf("summary_Phase%d.rds", PHASE))
saveRDS(
  list(
    FD_strict_mat      = FD_strict_mat,
    FD_anytime_mat     = FD_anytime_mat,
    FD_partial_mat     = FD_partial_mat,
    CP_mat             = CP_mat,
    detect_rate_df     = detect_rate_df,
    detect_strict_arr  = detect_strict_arr,
    detect_anytime_arr = detect_anytime_arr,
    mean_frac_arr      = mean_frac_arr,
    summary_df         = summary_df,
    Simulation         = Simulation,
    SCENARIOS          = SCENARIOS,
    run_config         = list(
      PHASE              = PHASE,
      N_ITER             = N_ITER,
      BASE_SEED          = BASE_SEED,
      REF_SEQ            = REF_SEQ,
      START_DATE         = START_DATE,
      END_DATE           = END_DATE,
      N_SEQ_PER_MONTH    = N_SEQ_PER_MONTH,
      MUT_RATE_VAR       = MUT_RATE_VAR,
      REF_VARIABILITY    = REF_VARIABILITY,
      SLIDING_WINDOW_LEN = SLIDING_WINDOW_LEN,
      WINDOW_OPTION      = WINDOW_OPTION,
      DETECTION_THRESHOLD = DETECTION_THRESHOLD,
      N_NEW_MUT_SEQ_VEC  = N_NEW_MUT_SEQ_VEC
    )
  ),
  file = summary_path
)
message(sprintf("Summary saved: %s", summary_path))


# ===========================================================================
# 11.  PUBLICATION-QUALITY PLOTS
# ===========================================================================

# ---------------------------------------------------------------------------
# Shared theme
# ---------------------------------------------------------------------------
theme_pub <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_blank(),
    strip.background    = element_rect(fill = "grey92", colour = "grey60"),
    strip.text          = element_text(face = "bold", size = 10),
    legend.position     = "bottom",
    legend.key.width    = unit(1.8, "cm"),
    legend.text         = element_text(size = 9),
    plot.title          = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle       = element_text(hjust = 0.5, size = 10, colour = "grey35"),
    plot.caption        = element_text(hjust = 1,   size  = 8,  colour = "grey45",
                                       lineheight = 1.3),
    axis.title          = element_text(size = 11),
    axis.text           = element_text(size = 10)
  )

# ---------------------------------------------------------------------------
# Utility: save as PDF (vector) and PNG (raster)
# ---------------------------------------------------------------------------
save_plot <- function(p, name, width = 8, height = 5.5) {
  base <- file.path(OUT_SUBDIR, sprintf("%s_Phase%d", name, PHASE))
  # cairo_pdf provides correct Unicode (e.g. >= symbol); fall back to pdf()
  pdf_dev <- if (capabilities("cairo")) cairo_pdf else pdf
  ggplot2::ggsave(paste0(base, ".pdf"), p,
                  width = width, height = height, device = pdf_dev)
  ggplot2::ggsave(paste0(base, ".png"), p,
                  width = width, height = height, dpi = 300)
  invisible(p)
}

# ---------------------------------------------------------------------------
# NA counts annotation helper
# ---------------------------------------------------------------------------
na_label <- function(mat) {
  n_na <- sum(is.na(mat))
  if (n_na == 0L) return(NULL)
  sprintf("%d iteration%s not detected within sweep (NA, omitted)",
          n_na, if (n_na == 1L) "" else "s")
}

# ---------------------------------------------------------------------------
# PLOT 1 — Distribution of First Detection Sample Size
# Three criteria shown side-by-side per scenario via dodged violin + boxplot.
# ---------------------------------------------------------------------------
crit_levels <- c(
  "Strict (on-time, all sites)",
  sprintf("Partial (\u2265%d%% of sites, on-time)",
          round(DETECTION_THRESHOLD * 100)),
  "Any-time (all sites)"
)

fd_long <- rbind(
  data.frame(scenario  = factor(rep(SC_NAMES, each = N_ITER), levels = SC_NAMES),
             value     = as.vector(FD_strict_mat),
             criterion = crit_levels[1L]),
  data.frame(scenario  = factor(rep(SC_NAMES, each = N_ITER), levels = SC_NAMES),
             value     = as.vector(FD_partial_mat),
             criterion = crit_levels[2L]),
  data.frame(scenario  = factor(rep(SC_NAMES, each = N_ITER), levels = SC_NAMES),
             value     = as.vector(FD_anytime_mat),
             criterion = crit_levels[3L])
)
fd_long$criterion <- factor(fd_long$criterion, levels = crit_levels)

y_max_fd <- max(N_NEW_MUT_SEQ_VEC) + 3
na_note  <- na_label(FD_strict_mat)

p1 <- ggplot(fd_long, aes(x = scenario, y = value, fill = criterion)) +
  geom_violin(alpha        = 0.22,
              scale        = "width",
              trim         = TRUE,
              colour       = NA,
              position     = position_dodge(0.85),
              na.rm        = TRUE) +
  geom_boxplot(width        = 0.12,
               outlier.size  = 0.9,
               outlier.alpha = 0.45,
               colour        = "grey25",
               linewidth     = 0.35,
               position      = position_dodge(0.85),
               na.rm         = TRUE) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  scale_y_continuous(
    limits = c(NA, y_max_fd),
    breaks = pretty(N_NEW_MUT_SEQ_VEC, n = 6),
    name   = expression(italic(n)[new~sequences])
  ) +
  labs(
    title    = "First Detection Sample Size by Scenario",
    subtitle = sprintf(
      "Phase %d  \u2014  %d iterations per scenario  \u2014  %d-aa reference",
      PHASE, N_ITER, nchar(REF_SEQ)
    ),
    x       = NULL,
    caption = paste0(
      "Strict: all variant sites in top-entropy cluster at expected partition; ",
      "site sets mutually distinct.\n",
      sprintf("Partial: \u2265%d%% of sites in top cluster at expected partition.  ",
              round(DETECTION_THRESHOLD * 100)),
      "Any-time: eventual detection permitted.\n",
      if (!is.null(na_note)) na_note else ""
    )
  ) +
  coord_flip() +
  theme_pub

save_plot(p1, "plot_fd_distribution")

# ---------------------------------------------------------------------------
# PLOT 2 — Change-Point Accuracy Distribution
# ---------------------------------------------------------------------------
cp_long <- data.frame(
  scenario = factor(rep(SC_NAMES, each = N_ITER), levels = SC_NAMES),
  value    = as.vector(CP_mat)
)

p2 <- ggplot(cp_long, aes(x = scenario, y = value, fill = scenario)) +
  geom_violin(alpha  = 0.28, scale = "width", trim = TRUE,
              colour = NA, na.rm = TRUE) +
  geom_boxplot(width         = 0.15,
               outlier.size  = 0.9,
               outlier.alpha = 0.45,
               colour        = "grey25",
               linewidth     = 0.35,
               fill          = "white",
               na.rm         = TRUE) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_y_continuous(
    limits = c(0, 1.04),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(accuracy = 1),
    name   = "Change-point accuracy"
  ) +
  labs(
    title    = "Change-Point Detection Accuracy by Scenario",
    subtitle = sprintf("Phase %d  \u2014  %d iterations per scenario", PHASE, N_ITER),
    x        = NULL,
    caption  = paste0(
      "CP accuracy: proportion of variants whose expected change point ",
      "was correctly identified by e.agglo on the local (up to part_em - 1) ",
      "Hellinger distance sequence."
    )
  ) +
  coord_flip() +
  theme_pub

save_plot(p2, "plot_cp_accuracy")

# ---------------------------------------------------------------------------
# PLOT 3 — Detection Rate Curves
# Detection probability as a function of n_new_seq, averaged across
# iterations.  Three criteria per scenario via linetype.
# ---------------------------------------------------------------------------
detect_long <- tidyr::pivot_longer(
  detect_rate_df,
  cols      = c("strict", "anytime", "mean_frac"),
  names_to  = "criterion",
  values_to = "rate"
) %>%
  dplyr::mutate(
    criterion = dplyr::recode(
      criterion,
      "strict"    = "Strict (on-time)",
      "anytime"   = "Any-time",
      "mean_frac" = "Mean site fraction"
    ),
    criterion = factor(criterion,
                       levels = c("Strict (on-time)",
                                  "Any-time",
                                  "Mean site fraction")),
    scenario  = factor(scenario, levels = SC_NAMES)
  )

p3 <- ggplot(detect_long, aes(x = n, y = rate,
                               colour   = scenario,
                               linetype = criterion)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  geom_hline(yintercept = 0.80, linetype = "dotdash",
             colour = "grey55", linewidth = 0.45) +
  annotate("text",
           x = max(N_NEW_MUT_SEQ_VEC),
           y = 0.83,
           label    = "80%",
           colour   = "grey45",
           size     = 3.2,
           hjust    = 1) +
  scale_colour_brewer(palette = "Set2", name = "Scenario") +
  scale_linetype_manual(
    values = c("Strict (on-time)"    = "solid",
               "Any-time"            = "dashed",
               "Mean site fraction"  = "dotted"),
    name   = "Detection criterion"
  ) +
  scale_x_continuous(
    name   = expression(italic(n)[new~sequences]),
    breaks = pretty(N_NEW_MUT_SEQ_VEC, n = 8)
  ) +
  scale_y_continuous(
    name   = "Detection probability",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title    = "Variant Detection Rate vs. De-Novo Sequence Count",
    subtitle = sprintf("Phase %d  \u2014  %d iterations per scenario", PHASE, N_ITER),
    caption  = paste0(
      "Each curve is the mean detection rate across iterations at a given ",
      "n_new_seq value.\nHorizontal reference at 80% detection rate."
    )
  ) +
  guides(colour   = guide_legend(nrow = 2L, override.aes = list(linewidth = 1.2)),
         linetype = guide_legend(nrow = 2L)) +
  theme_pub

save_plot(p3, "plot_detect_rate_curves", width = 9, height = 5.5)

# ---------------------------------------------------------------------------
# PLOT 4 — Time-to-Detection Distribution  (Phase 2 and 3 only)
# Requires extracting time_to_detection from per-iteration RDS files.
# ---------------------------------------------------------------------------
if (PHASE >= 2L) {

  ttd_rows <- do.call(rbind, lapply(seq_len(N_SCEN), function(s) {
    do.call(rbind, lapply(seq_len(N_ITER), function(j) {
      meta <- Simulation[[s]][[j]]
      if (is.null(meta$rds_path) || !file.exists(meta$rds_path)) return(NULL)
      r <- tryCatch(readRDS(meta$rds_path), error = function(e) NULL)
      if (is.null(r) || is.null(r$results)) return(NULL)
      res <- r$results
      # One representative row per (variant, n_new_seq): use the row at
      # the first n_new_seq where the variant was detected on time.
      n_ref <- unique(res$num_months_ref_seq)
      n_ref <- if (length(n_ref) == 1L) n_ref else NA_integer_
      data.frame(
        scenario          = SC_NAMES[s],
        iteration         = j,
        variant           = res$variant,
        num_denovo        = res$num_denovo_sequences,
        time_to_detection = as.numeric(res$time_to_detection),
        stringsAsFactors  = FALSE
      )
    }))
  }))

  if (!is.null(ttd_rows) && nrow(ttd_rows) > 0L) {
    ttd_rows$scenario <- factor(ttd_rows$scenario, levels = SC_NAMES)

    p4 <- ggplot(ttd_rows[!is.na(ttd_rows$time_to_detection), ],
                 aes(x = scenario, y = time_to_detection, fill = scenario)) +
      geom_violin(alpha  = 0.28, scale = "width", trim = FALSE,
                  colour = NA, na.rm = TRUE) +
      geom_boxplot(width         = 0.15,
                   outlier.size  = 0.9,
                   outlier.alpha = 0.45,
                   colour        = "grey25",
                   linewidth     = 0.35,
                   fill          = "white",
                   na.rm         = TRUE) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 colour = "grey50", linewidth = 0.4) +
      scale_fill_brewer(palette = "Set2", guide = "none") +
      scale_y_continuous(name = "Time to detection (months)") +
      labs(
        title    = "Time to Detection Relative to Emergence",
        subtitle = sprintf("Phase %d  \u2014  %d iterations per scenario", PHASE, N_ITER),
        x        = NULL,
        caption  = paste0(
          "Time to detection = detection_month \u2212 actual_emergence_month. ",
          "Zero indicates detection within the expected 2-month window. ",
          "Positive values indicate delayed detection."
        )
      ) +
      coord_flip() +
      theme_pub

    save_plot(p4, "plot_time_to_detection")
  }
}

# ---------------------------------------------------------------------------
# Print file locations
# ---------------------------------------------------------------------------
message(sprintf("\nAll outputs written to:\n  %s", OUT_SUBDIR))
message(sprintf("Plots: plot_fd_distribution, plot_cp_accuracy, plot_detect_rate_curves%s",
                if (PHASE >= 2L) ", plot_time_to_detection" else ""))
message(sprintf("Summary RDS: %s", basename(summary_path)))

# Shutdown parallel workers cleanly
future::plan(future::sequential)
