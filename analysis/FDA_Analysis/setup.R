################################################################################
## setup.R — FDA_Analysis
##
## Configuration, paths, RNG seeds, atomic I/O, logging, date arithmetic,
## truth-catalogue loader, and feature-matrix loader.
##
## Mirrors Changepoint_Detection_Study/setup.R conventions exactly.
## Sourced by every entry point of FDA_Analysis.
################################################################################

# ── 1. Packages ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ViralEntropR)
})

.have_fdacluster <- requireNamespace("fdacluster", quietly = TRUE)
.have_mclust     <- requireNamespace("mclust",     quietly = TRUE)
.have_cluster    <- requireNamespace("cluster",    quietly = TRUE)
.have_magick     <- requireNamespace("magick",     quietly = TRUE)
.have_callr      <- requireNamespace("callr",      quietly = TRUE)
.have_ggplot2    <- requireNamespace("ggplot2",    quietly = TRUE)


# ── 2. Study-directory resolution ────────────────────────────────────────────
#
# Resolution order:
#   (1) $VIRAL_FDA_STUDY_DIR (env override)
#   (2) Known absolute paths on the dissertation server / workstation
#   (3) ./analysis/FDA_Analysis under getwd() (package working-copy)
#   (4) getwd() as a last resort

resolve_study_dir <- function() {
  candidates <- c(
    Sys.getenv("VIRAL_FDA_STUDY_DIR", unset = ""),
    "/home/vadimtyu/Clean_Code_Running/Dissertation/FDA_Analysis",
    "/home/vadimtyu/Clean_Code_Running/Dissertation/ViralEntropR/analysis/FDA_Analysis",
    "C:/YORK_PhD/RESEARCH/PAPERS/GitHub/ViralEntropR/analysis/FDA_Analysis",
    file.path(getwd(), "analysis", "FDA_Analysis"),
    getwd()
  )
  candidates <- candidates[nzchar(candidates)]
  for (cand in candidates) {
    if (dir.exists(cand) && file.exists(file.path(cand, "setup.R")))
      return(normalizePath(cand, winslash = "/"))
  }
  normalizePath(getwd(), winslash = "/")
}


# ── 3. Worker auto-detection ─────────────────────────────────────────────────
#
# FDA cells are memory-bound from fdacluster's pairwise curve distances *and*
# from the bootstrap (B = 100 fdahclust refits per cell × 2 metrics). Hard cap
# is 8 — lower than Changepoint_Detection_Study's 24 — to leave headroom.

auto_detect_n_workers <- function(worker_ram_gb = 2.0, hard_cap = 64L) {
  cpu_cap <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  ram_cap <- tryCatch({
    if (Sys.info()[["sysname"]] == "Linux") {
      ram_kb <- as.numeric(system(
        "awk '/MemAvailable/ {print $2}' /proc/meminfo",
        intern = TRUE, ignore.stderr = TRUE))
      floor((ram_kb / 1e6) / worker_ram_gb)
    } else {
      hard_cap
    }
  }, error = function(e) hard_cap)
  max(1L, min(as.integer(cpu_cap), as.integer(ram_cap), as.integer(hard_cap)))
}


# ── 4. Configuration builder ─────────────────────────────────────────────────

build_config <- function() {
  config <- list(

    # — Schema (bump when RDS payload shape changes) —
    SCHEMA_VERSION           = 2L,

    # — Reproducibility —
    BASE_SEED                = 2025L,

    # — Cell definition (full WHO catalogue) —
    # 12 WHO-labelled variants. `NULL` here means "all variants with a
    # parseable US detection date in sarscov2_variants" — resolved at runtime.
    VARIANTS_TO_RUN          = NULL,
    STRATEGIES_TO_RUN        = c("sliding_2m", "disjoint_2m", "cumulative_1m"),

    # — Per-strategy partitioning specs (window_type per partition_time_windows) —
    STRATEGY_DEFS = list(
      sliding_2m    = list(strategy_id = 1L, window_type = 2L, window_length = 2L),
      disjoint_2m   = list(strategy_id = 2L, window_type = 3L, window_length = 2L),
      cumulative_1m = list(strategy_id = 3L, window_type = 1L, window_length = 1L)
    ),

    # — Window construction (per cell, around detection date) —
    K_BINS_BEFORE            = 5L,
    K_BINS_AFTER             = 5L,
    REF_BIN_MONTHS           = 2L,   # used only for window sizing

    # — Encoding / metric —
    N_SITES                  = 1273L,
    AA_LEVELS                = 25L,
    HELLINGER_NORMALIZED     = FALSE,
    REFERENCE_IDX_HELLINGER  = 1L,

    # — GMM (per-cell window, shared across strategies) —
    MCLUST_G                 = 1:15,
    MCLUST_MODELS            = NULL,

    # — Quality thresholds —
    MIN_SEQUENCES_PER_BIN    = 30L,
    MIN_CLASS1_SITES_FOR_FDA = 4L,

    # — fdahclust —
    N_CLUSTERS_GRID          = 2:6,
    FDA_WARPING_CLASS        = "none",
    FDA_METRIC               = "l2",
    FDA_CENTROID_TYPE        = "mean",
    FDA_LINKAGE              = "complete",
    LOG1P_TRANSFORM          = TRUE,

    # — Per-frame design (v2) —
    # Fitting sub-window for the animation is `centre ± FITTING_WINDOW_HALFWIDTH`
    # partitions; frames slide one partition per step.
    FITTING_WINDOW_HALFWIDTH = 2L,
    MIN_CLASS1_PER_FRAME     = 3L,        # need ≥ k+1 ≥ 3 for k=2
    MIN_MEAN_SILHOUETTE      = 0.25,      # Kaufman & Rousseeuw 1990 floor
    MAX_CLUSTER_SIZE         = 33L,       # Omicron's mutation count

    # — Hellinger anchor scheme (v2 dual anchor) —
    # T1   = window-start partition (always used)
    # Tpred = partition containing the predecessor variant's US detection
    # The predecessor mapping is hard-coded below.
    HELLINGER_ANCHORS        = c("T1", "Tpred"),
    D614G_US_DETECTION       = as.Date("2020-04-01"),   # Korber et al. 2020

    # Predecessor variant per WHO label. Variants whose predecessor is
    # "D614G" use the canonical Apr-2020 anchor directly.
    PREDECESSOR_MAP = c(
      Alpha   = "D614G",
      Beta    = "Alpha",
      Gamma   = "Alpha",
      Delta   = "Beta",
      Omicron = "Delta",
      Epsilon = "D614G",
      Zeta    = "D614G",
      Eta     = "D614G",
      Iota    = "D614G",
      Theta   = "Alpha",
      Kappa   = "Alpha",
      Lambda  = "Alpha"
    ),

    # — Bootstrap stability (Hennig 2007) —
    N_BOOTSTRAP              = 100L,
    STABLE_JACCARD_THRESHOLD = 0.75,

    # — SNP enrichment Bonferroni denominator —
    # With 12 variants × 2 datasets × 3 strategies × 3 metrics, minus
    # out-of-coverage cells (Omicron × NCBI_US: 3 cells). Cap ≈ 213; we use
    # 220 for a small safety margin.
    N_TOTAL_ENRICHMENT_TESTS = 220L,

    # — Animation / static plot dimensions (dual-panel layout) —
    ANIMATION_FPS            = 1L, #2L
    GIF_WIDTH_PX             = 2040L,
    GIF_HEIGHT_PX            = 1800L,
    GIF_RES                  = 300L,
    PANEL_TOP_FRAC           = 0.60,     # top panel ≈ 60% of vertical real estate
    PNG_WIDTH_CM             = 17,
    PNG_HEIGHT_CM            = 15,

    # — Datasets —
    DATASETS = list(
      NCBI_US = list(
        id          = 1L,
        label       = "NCBI US",
        feature_rds = "/home/vadimtyu/Clean_Code_Running/Dissertation/NCBI_US_unaligned_feature_matrix_1273aa.rds"
      ),
      GISAID_US = list(
        id          = 2L,
        label       = "GISAID US",
        feature_rds = "/home/vadimtyu/Clean_Code_Running/Dissertation/GISAID_US_aligned_feature_matrix_1273aa.rds"
      )
    ),

    # — Paths —
    STUDY_DIR                = resolve_study_dir(),
    OUTPUT_SUBDIR            = "outputs",

    # — Concurrency —
    N_WORKERS                = NULL,       # auto-detect when NULL
    SUBPROCESS_POLL_S        = 0.25
  )

  # ENV overrides for feature-matrix paths.
  for (ds_name in names(config$DATASETS)) {
    env_key  <- sprintf("FEATURE_RDS_%s", ds_name)
    override <- Sys.getenv(env_key, unset = "")
    if (nzchar(override))
      config$DATASETS[[ds_name]]$feature_rds <- override
  }

  # ENV override for worker count.
  env_workers <- Sys.getenv("N_WORKERS", unset = "")
  if (nzchar(env_workers)) {
    val <- suppressWarnings(as.integer(env_workers))
    if (!is.na(val) && val >= 1L) config$N_WORKERS <- val
  }
  if (is.null(config$N_WORKERS)) config$N_WORKERS <- auto_detect_n_workers()

  config
}


# ── 5. ID lookups (deterministic, used in seed derivation) ───────────────────

variant_id_for <- function(variant_name, config) {
  vars <- config$VARIANTS_TO_RUN
  if (is.null(vars)) {
    truth <- load_truth_catalogue(variants_to_run = NULL)
    vars  <- sort(unique(truth$WHO_Label))
  }
  vidx <- match(variant_name, vars)
  if (is.na(vidx))
    stop("Unknown variant: ", variant_name, call. = FALSE)
  as.integer(vidx)
}

strategy_id_for <- function(strategy_name, config) {
  sdef <- config$STRATEGY_DEFS[[strategy_name]]
  if (is.null(sdef))
    stop("Unknown strategy: ", strategy_name, call. = FALSE)
  as.integer(sdef$strategy_id)
}

dataset_id_for <- function(dataset_name, config) {
  ds <- config$DATASETS[[dataset_name]]
  if (is.null(ds))
    stop("Unknown dataset: ", dataset_name, call. = FALSE)
  as.integer(ds$id)
}


# ── 6. Path helpers ──────────────────────────────────────────────────────────

output_dir <- function(config)
  file.path(config$STUDY_DIR, config$OUTPUT_SUBDIR)

cells_dir <- function(config)
  file.path(output_dir(config), "cells")

cell_path <- function(config, dataset_name, variant_name, strategy_name)
  file.path(cells_dir(config),
            sprintf("cell_%s__%s__%s.rds",
                    dataset_name, variant_name, strategy_name))

plots_dir   <- function(config) file.path(output_dir(config), "plots")
gifs_dir    <- function(config) file.path(output_dir(config), "gifs")
tables_dir  <- function(config) file.path(output_dir(config), "tables")
logs_dir    <- function(config) file.path(output_dir(config), "logs")

error_log_path <- function(config) file.path(logs_dir(config), "error_log.txt")
run_log_path   <- function(config) file.path(logs_dir(config), "run_log.txt")


# ── 7. RNG seeds ─────────────────────────────────────────────────────────────
#
# Cell-level:     BASE_SEED + dataset_id*1e8 + variant_id*1e6
# Strategy-level: cell_seed + strategy_id*1e3
#
# All seeds fit in 32-bit signed int (max value ~2.05e8 ≪ 2^31).

seed_for_cell <- function(base_seed, dataset_id, variant_id) {
  as.integer(base_seed +
             as.numeric(dataset_id) * 1e8 +
             as.numeric(variant_id) * 1e6)
}

seed_for_strategy <- function(cell_seed, strategy_id) {
  as.integer(cell_seed + as.numeric(strategy_id) * 1e3)
}


# ── 8. Null-coalescing operator (generic utility) ────────────────────────────
#
# Used throughout aggregate_tables.R to provide safe defaults when fields are
# absent from older or partial cell RDS files.

`%||%` <- function(a, b) if (is.null(a)) b else a


# ── 9. Atomic RDS write (parallel-safe) ──────────────────────────────────────
#
# Writes to a unique temp file, then renames. file.rename on POSIX is atomic
# within a filesystem; on Windows we fall back to copy+remove.

save_rds_atomic <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid(), "-",
                sample.int(.Machine$integer.max, 1L))
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    file.remove(tmp)
  }
  invisible(path)
}


# ── 10. Logging ───────────────────────────────────────────────────────────────

log_msg <- function(..., con = stderr()) {
  msg <- sprintf("[%s] %s",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, "\n", sep = "", file = con)
  if (inherits(con, "connection")) try(flush(con), silent = TRUE)
  invisible(msg)
}

log_error <- function(config, ...) {
  path <- error_log_path(config)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  msg <- sprintf("[%s] %s\n",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, file = path, append = TRUE)
  invisible(msg)
}


# ── 11. Date arithmetic (timezone-independent, month-aligned) ────────────────

add_months <- function(date, n_months) {
  if (length(date) == 0L) return(date)
  y     <- as.integer(format(date, "%Y"))
  m0    <- as.integer(format(date, "%m")) - 1L
  total <- (y * 12L + m0) + as.integer(n_months)
  new_y <- total %/% 12L
  new_m <- (total %% 12L) + 1L
  as.Date(sprintf("%04d-%02d-01", new_y, new_m))
}

diff_months <- function(date_to, date_from) {
  y_to   <- as.integer(format(date_to,   "%Y"))
  m_to   <- as.integer(format(date_to,   "%m"))
  y_from <- as.integer(format(date_from, "%Y"))
  m_from <- as.integer(format(date_from, "%m"))
  as.integer(12L * (y_to - y_from) + (m_to - m_from))
}

snap_to_month <- function(date) {
  if (length(date) == 0L) return(date)
  as.Date(format(date, "%Y-%m-01"))
}


# ── 12. Month-Year date parser (Sep-2020 → 2020-09-01) ───────────────────────

.MONTH_ABBR <- c(Jan = 1L,  Feb = 2L,  Mar = 3L,  Apr = 4L,
                 May = 5L,  Jun = 6L,  Jul = 7L,  Aug = 8L,
                 Sep = 9L,  Oct = 10L, Nov = 11L, Dec = 12L)

parse_month_year <- function(s) {
  if (length(s) == 0L) return(as.Date(character(0L)))
  out <- as.Date(rep(NA_character_, length(s)))
  for (i in seq_along(s)) {
    si <- s[i]
    if (is.na(si) || identical(si, "") || identical(si, "NA")) next
    parts <- strsplit(si, "-", fixed = TRUE)[[1L]]
    if (length(parts) != 2L) next
    mn <- .MONTH_ABBR[parts[1L]]
    if (is.na(mn)) next
    yr <- suppressWarnings(as.integer(parts[2L]))
    if (is.na(yr)) next
    out[i] <- as.Date(sprintf("%04d-%02d-01", yr, mn))
  }
  out
}


# ── 13. Truth-catalogue loader (ground-truth SNP and mutation sites) ─────────

load_truth_catalogue <- function(variants_to_run = NULL) {
  v <- ViralEntropR::sarscov2_variants
  expected <- c("WHO_Label", "Pango_Lineage", "Date_First_Detected_US",
                "Defining_SNP_Sites", "Mutation_Sites")
  missing  <- setdiff(expected, names(v))
  if (length(missing) > 0L)
    stop("sarscov2_variants is missing expected fields: ",
         paste(missing, collapse = ", "), call. = FALSE)

  who   <- unlist(v$WHO_Label, use.names = FALSE)
  pango <- vapply(v$Pango_Lineage,
                  function(x) paste(unlist(x, use.names = FALSE), collapse = "/"),
                  character(1L))
  us_date <- parse_month_year(as.character(v$Date_First_Detected_US))

  country <- if ("Country_First_Detected" %in% names(v))
               vapply(v$Country_First_Detected,
                      function(x) paste(unlist(x, use.names = FALSE),
                                        collapse = "/"),
                      character(1L))
             else
               rep(NA_character_, length(who))

  to_int_vec <- function(x) {
    out <- suppressWarnings(as.integer(unlist(x, use.names = FALSE)))
    sort(unique(out[!is.na(out)]))
  }

  df <- data.frame(
    WHO_Label              = who,
    Pango_Lineage          = pango,
    Country_First_Detected = country,
    Date_First_Detected_US = us_date,
    stringsAsFactors       = FALSE
  )
  df$Defining_SNP_Sites <- I(lapply(v$Defining_SNP_Sites, to_int_vec))
  df$Mutation_Sites     <- I(lapply(v$Mutation_Sites,     to_int_vec))

  if (!is.null(variants_to_run))
    df <- df[df$WHO_Label %in% variants_to_run, , drop = FALSE]
  df <- df[!is.na(df$Date_First_Detected_US), , drop = FALSE]
  rownames(df) <- NULL
  df
}


# ── 14. Feature-matrix loader ───────────────────────────────────────────────

load_feature_matrix <- function(rds_path, expected_n_sites = 1273L) {
  if (!file.exists(rds_path))
    stop("Feature matrix RDS not found: ", rds_path, call. = FALSE)

  fm <- readRDS(rds_path)
  if (!is.data.frame(fm))
    stop("Feature matrix must be a data.frame; got ",
         paste(class(fm), collapse = "/"), ".", call. = FALSE)
  if (!"Date" %in% colnames(fm))
    stop("Feature matrix lacks a Date column.", call. = FALSE)
  if (!inherits(fm$Date, "Date"))
    stop("Feature matrix Date column must be class Date.", call. = FALSE)

  site_cols <- setdiff(colnames(fm), c("Date", "Country"))
  if (length(site_cols) != expected_n_sites)
    stop("Feature matrix has ", length(site_cols),
         " site columns; expected ", expected_n_sites, ".", call. = FALSE)

  if (is.unsorted(fm$Date)) {
    fm <- fm[order(fm$Date), , drop = FALSE]
    rownames(fm) <- NULL
  }

  fm
}


# ── 15. Runtime-package guard ───────────────────────────────────────────────

require_runtime_packages <- function(animation = FALSE) {
  if (!.have_fdacluster)
    stop("Package 'fdacluster' is required: install.packages('fdacluster').",
         call. = FALSE)
  if (!.have_mclust)
    stop("Package 'mclust' is required.", call. = FALSE)
  if (!.have_cluster)
    stop("Package 'cluster' is required.", call. = FALSE)
  if (!.have_ggplot2)
    stop("Package 'ggplot2' is required.", call. = FALSE)
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required for dual-panel plots.",
         call. = FALSE)
  if (animation && !.have_magick)
    stop("Package 'magick' is required for animation.",
         " On Linux: apt install libmagick++-dev, then install.packages('magick').",
         call. = FALSE)
  invisible(TRUE)
}


# ── 16. Cell enumeration ────────────────────────────────────────────────────

enumerate_cells <- function(config) {
  ds   <- names(config$DATASETS)
  vars <- config$VARIANTS_TO_RUN
  # Resolve NULL → all variants in the catalogue with a parseable US date.
  if (is.null(vars)) {
    truth <- load_truth_catalogue(variants_to_run = NULL)
    vars  <- sort(unique(truth$WHO_Label))
  }
  strs <- config$STRATEGIES_TO_RUN
  expand.grid(
    dataset  = ds,
    variant  = vars,
    strategy = strs,
    KEEP.OUT.ATTRS   = FALSE,
    stringsAsFactors = FALSE
  )
}


# ── 17. Predecessor anchor resolution ───────────────────────────────────────
#
# Maps a variant to the partition index in the cell's window that contains the
# variant's "closest predecessor" detection date (or the D614G canonical
# anchor for variants without a prior named predecessor).
#
# Returns:
#   list(
#     predecessor_name      = char,   "D614G" or another variant's WHO_Label
#     predecessor_us_date   = Date    "first detected in US" of predecessor
#     bin_idx               = int     index in part$bin_starts / NA if outside
#                                      (then snapped to nearest bin)
#     snapped               = bool    TRUE if outside the window and clipped
#   )

resolve_predecessor <- function(variant_name, truth_df, config) {
  pred_name <- config$PREDECESSOR_MAP[[variant_name]]
  if (is.null(pred_name) || !nzchar(pred_name))
    return(list(predecessor_name   = NA_character_,
                predecessor_us_date = as.Date(NA),
                status             = "no_predecessor"))

  if (identical(pred_name, "D614G"))
    return(list(predecessor_name   = "D614G",
                predecessor_us_date = config$D614G_US_DETECTION,
                status             = "ok"))

  row <- truth_df[truth_df$WHO_Label == pred_name, , drop = FALSE]
  if (nrow(row) == 0L || is.na(row$Date_First_Detected_US[1L]))
    return(list(predecessor_name   = pred_name,
                predecessor_us_date = as.Date(NA),
                status             = "predecessor_missing"))

  list(predecessor_name   = pred_name,
       predecessor_us_date = row$Date_First_Detected_US[1L],
       status             = "ok")
}
