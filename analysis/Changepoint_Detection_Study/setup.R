################################################################################
## setup.R — Changepoint_Detection_Study
##
## Configuration, paths, RNG seeds, utility helpers, and reference-data loading.
## Sourced by every entry point.
##
## Mirrors Sample_Size_Simulation_Study/setup.R conventions.
################################################################################

# ── 1. Packages ──────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ViralEntropR)
})

.have_ecp   <- requireNamespace("ecp",        quietly = TRUE)
.have_hdcp  <- requireNamespace("HDcpDetect", quietly = TRUE)
.have_callr <- requireNamespace("callr",      quietly = TRUE)

# ── 2. Study-directory resolution ────────────────────────────────────────────

resolve_study_dir <- function() {
  candidates <- c(
    Sys.getenv("VIRAL_CP_STUDY_DIR", unset = ""),
    "/home/vadimtyu/Clean_Code_Running/Dissertation/Changepoint_Detection_Study",
    "/home/vadimtyu/Clean_Code_Running/Dissertation/ViralEntropR/analysis/Changepoint_Detection_Study",
    "C:/YORK_PhD/RESEARCH/PAPERS/GitHub/ViralEntropR/analysis/Changepoint_Detection_Study",
    file.path(getwd(), "analysis", "Changepoint_Detection_Study"),
    getwd()
  )
  candidates <- candidates[nzchar(candidates)]
  for (cand in candidates) {
    if (dir.exists(cand) && file.exists(file.path(cand, "setup.R")))
      return(normalizePath(cand))
  }
  normalizePath(getwd())
}

# ── 3. Worker auto-detection ─────────────────────────────────────────────────

auto_detect_n_workers <- function(worker_ram_gb = 2.5, hard_cap = 24L) {
  cpu_cap <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  ram_cap <- tryCatch({
    if (Sys.info()[["sysname"]] == "Linux") {
      ram_kb <- as.numeric(system(
        "awk '/MemAvailable/ {print $2}' /proc/meminfo",
        intern = TRUE))
      floor((ram_kb / 1e6) / worker_ram_gb)
    } else {
      hard_cap
    }
  }, error = function(e) hard_cap)
  max(1L, min(as.integer(cpu_cap), as.integer(ram_cap), hard_cap))
}

# ── 4. Configuration builder ─────────────────────────────────────────────────

build_config <- function() {

  config <- list(

    # — Reproducibility —
    BASE_SEED               = 2025L,

    # — Design constants —
    # Even-only lengths so bin count = length / 2 is integer.
    TIER_LENGTHS = list(
      short  = c(6L, 8L, 10L),
      medium = c(12L, 14L, 16L),
      long   = c(18L, 20L, 22L, 24L, 26L, 28L)
    ),
    N_RUNS_PER_TIER         = 1667L, 
    BIN_MONTHS              = 2L,
    MIN_SEQUENCES_PER_BIN   = 30L,
    MIN_TRUTHS_PER_WINDOW   = 2L,            # in UNIQUE Hellinger bins

    # — Hellinger —
    AA_LEVELS               = 25L,
    HELLINGER_NORMALIZED    = FALSE,

    # — Change-point methods —
    KS_CP3O_MINSIZE         = 2L,
    KS_CP3O_K_SWEEP         = c(1L, 2L, 3L, 5L, 10L),
    E_AGGLO_ALPHA           = 1,
    HDCP_M                  = 100L,

    # — GMM (per-window) —
    MCLUST_MODELS           = NULL,
    MCLUST_G                = 1:15,

    # — Sensitivity sweeps —
    TRUTH_SHIFT_BINS        = c(-2L, -1L, 0L, 1L, 2L),
    RUN_TRUTH_SHIFT_SWEEP   = TRUE,

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
    STUDY_DIR               = resolve_study_dir(),
    OUTPUT_SUBDIR           = "outputs", 

    # — Concurrency —
    N_WORKERS               = 24L,
    SUBPROCESS_POLL_S       = 0.20
  )

  for (ds_name in names(config$DATASETS)) {
    env_key <- sprintf("FEATURE_RDS_%s", ds_name)
    override <- Sys.getenv(env_key, unset = "")
    if (nzchar(override))
      config$DATASETS[[ds_name]]$feature_rds <- override
  }

  config
}

# ── 5. Path helpers ──────────────────────────────────────────────────────────

output_dir <- function(config) {
  file.path(config$STUDY_DIR, config$OUTPUT_SUBDIR)
}

windows_path <- function(config, dataset_name) {
  file.path(output_dir(config), sprintf("windows_%s.rds", dataset_name))
}

replicates_dir <- function(config, dataset_name, base_dir = output_dir(config)) {
  file.path(base_dir, sprintf("replicates_%s", dataset_name))
}

replicate_path <- function(config, dataset_name, cell_id,
                           base_dir = output_dir(config)) {
  file.path(
    replicates_dir(config, dataset_name, base_dir = base_dir),
    sprintf("run_%04d.rds", cell_id)
  )
}

summary_path <- function(config, dataset_name) {
  file.path(output_dir(config), sprintf("summary_%s.rds", dataset_name))
}

error_log_path <- function(config) {
  file.path(output_dir(config), "error_log.txt")
}

run_log_path <- function(config) {
  file.path(output_dir(config), "run_log.txt")
}

# ── 6. RNG seeds ─────────────────────────────────────────────────────────────
#
# seed = BASE_SEED + dataset_id*1e8 + tier_id*1e6 + run_id*1e3
# Max with dataset_id ≤ 2, tier_id ≤ 3, run_id ≤ 1667:
#   2025 + 2e8 + 3e6 + 1.667e6 ≈ 2.05e8, well within int32 range.

seed_for_replicate <- function(base_seed, dataset_id, tier_id, run_id) {
  as.integer(base_seed +
             dataset_id * 1e8 +
             tier_id    * 1e6 +
             run_id     * 1e3)
}

seed_for_window_pool <- function(base_seed, dataset_id, tier_id) {
  as.integer(base_seed + dataset_id * 1e8 + tier_id * 1e6)
}

# ── 7. I/O helpers ───────────────────────────────────────────────────────────

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

# ── 8. Logging ───────────────────────────────────────────────────────────────

log_msg <- function(..., con = stderr()) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, "\n", sep = "", file = con)
  flush(con)
  invisible(msg)
}

log_error <- function(config, ..., con = error_log_path(config)) {
  dir.create(dirname(con), recursive = TRUE, showWarnings = FALSE)
  msg <- sprintf("[%s] %s\n",
                 format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                 paste0(..., collapse = ""))
  cat(msg, file = con, append = TRUE)
  invisible(msg)
}

# ── 9. "Month-Year" date parser ──────────────────────────────────────────────
#
# `sarscov2_variants$Date_First_Detected_US` is a character vector of strings
# like "Sep-2020". We standardise to first-of-month Date class. NA strings,
# empty strings, and the literal "NA" all map to NA Date.

.MONTH_ABBR <- c(Jan = 1L, Feb = 2L, Mar = 3L, Apr = 4L, May  = 5L,  Jun = 6L,
                 Jul = 7L, Aug = 8L, Sep = 9L, Oct = 10L, Nov = 11L, Dec = 12L)

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

# ── 10. Ground-truth catalogue loader ────────────────────────────────────────
#
# `sarscov2_variants` is a named LIST of 12 WHO-labelled variants built from
# `SARS_CoV_2_VOC_VOI.xlsx` (every entry is by construction either a VOC or
# a VOI — there is no separate status column to filter on). This loader:
#   1. Reads the list elements relevant for truth matching.
#   2. Parses Date_First_Detected_US from "Sep-2020" → as.Date("2020-09-01").
#   3. Drops entries with NA US-detection dates.
#   4. Returns a data.frame sorted by US date, ready for window filtering.

load_truth_catalogue <- function() {

  v <- ViralEntropR::sarscov2_variants

  expected <- c("WHO_Label", "Pango_Lineage", "Date_First_Detected_US")
  missing  <- setdiff(expected, names(v))
  if (length(missing) > 0L)
    stop("sarscov2_variants is missing expected fields: ",
         paste(missing, collapse = ", "), call. = FALSE)

  who <- unlist(v$WHO_Label,    use.names = FALSE)
  # Pango_Lineage may be a list of character vectors (some variants have
  # multiple sublineages). Collapse to slash-separated string.
  pango <- vapply(v$Pango_Lineage,
                  function(x) paste(unlist(x, use.names = FALSE),
                                    collapse = "/"),
                  character(1L))
  us_str  <- as.character(v$Date_First_Detected_US)
  us_date <- parse_month_year(us_str)

  df <- data.frame(
    WHO_Label              = who,
    Pango_Lineage          = pango,
    Date_First_Detected_US = us_date,
    stringsAsFactors       = FALSE
  )

  # Drop rows with NA US-detection dates (variants never sequenced in US).
  df <- df[!is.na(df$Date_First_Detected_US), , drop = FALSE]
  df <- df[order(df$Date_First_Detected_US), , drop = FALSE]
  rownames(df) <- NULL

  df
}

# ── 11. Feature-matrix loader ────────────────────────────────────────────────

load_feature_matrix <- function(rds_path, expected_n_sites = 1273L) {
  if (!file.exists(rds_path))
    stop("Feature matrix RDS not found: ", rds_path, call. = FALSE)

  fm <- readRDS(rds_path)

  if (!is.data.frame(fm))
    stop("Feature matrix must be a data.frame; got class ",
         paste(class(fm), collapse = "/"), ".", call. = FALSE)
  if (!"Date" %in% colnames(fm))
    stop("Feature matrix lacks a Date column.", call. = FALSE)
  if (!inherits(fm$Date, "Date"))
    stop("Feature matrix Date column must be of class Date.", call. = FALSE)

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

# ── 12. Worker package banner (called inside subprocesses) ───────────────────

require_runtime_packages <- function() {
  if (!.have_ecp)
    stop("Package 'ecp' is required but not installed.", call. = FALSE)
  if (!.have_hdcp)
    stop("Package 'HDcpDetect' is required but not installed.", call. = FALSE)
  invisible(TRUE)
}
