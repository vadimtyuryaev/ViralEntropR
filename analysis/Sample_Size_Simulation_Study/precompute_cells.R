# =============================================================================
# precompute_cells.R
# =============================================================================
#
# Deterministic cell-table builder for the full factorial design.
#
# Produces, for each of the four scenarios, a data frame in which each row
# is one cell of the (band x N_ref x ratio_axes) grid. Cells are
# deterministic; no randomness is involved. The tables are serialised once
# at study startup and reused across runs.
#
# Scenarios (in this study; all detect a newly arising variant against a
# non-trivial established background):
#   Sc1: ref + V1                  -> detect V2  (POOL_11_DRAW, free count)
#   Sc2: ref + V1 + V2             -> detect V3  (POOL_11_DRAW, free count)
#   Sc3: ref + V1 + V2 + V3        -> detect V4  (POOL_11_DRAW, free count)
#   Sc4: ref + V1 + V2 + V3        -> detect V4  (POOL_12_DRAW, 33 mutations)
#
# Factorial structure (ratios = config$RATIOS, 5 values by default):
#   Sc1: band x N_ref x ratio_V1                          ->     150 cells
#   Sc2: band x N_ref x ratio_V1 x ratio_V2               ->     750 cells
#   Sc3: band x N_ref x ratio_V1 x ratio_V2 x ratio_V3    ->  3,750 cells
#   Sc4: same structure as Sc3                            ->  3,750 cells
#   Total: 8,400 cells x 30 replicates = 252,000 simulations.
#
# Ceiling per scenario:
#   Sc1: n_V1
#   Sc2: max(n_V1, n_V2)
#   Sc3, Sc4: max(n_V1, n_V2, n_V3)
#
# Author : Vadim Tyuryaev
# =============================================================================

# -----------------------------------------------------------------------------
# n_ref_grid
# -----------------------------------------------------------------------------
# Log-spaced N_ref integer grid within a single band.
#
# Inputs
#   band_range : integer vector of length 2, c(lo, hi).
#   n_points   : integer scalar. Number of grid points.
#
# Returns: integer vector of unique grid values. The endpoints lo and hi
# are always included. Duplicates that arise from rounding are removed via
# unique(), so the returned length may be slightly less than n_points in
# very narrow bands (not an issue for the bands used in this study).
# -----------------------------------------------------------------------------
n_ref_grid <- function(band_range, n_points) {
  unique(as.integer(round(exp(
    seq(log(band_range[1L]), log(band_range[2L]),
        length.out = n_points)
  ))))
}

# -----------------------------------------------------------------------------
# precompute_cells_for_scenario
# -----------------------------------------------------------------------------
# Builds the cells data frame for one scenario.
#
# Inputs
#   scenario : integer scalar in 1..4.
#   config   : output of build_config().
#
# Returns: data frame with columns
#   scenario, cell_id, band, N_ref,
#   ratio_V1, n_V1, ratio_V2, n_V2, ratio_V3, n_V3,
#   ceiling.
#
# Columns absent from a scenario's design carry NA. Construction is
# vectorised via expand.grid() inside a per-band lapply. No loops over
# cells.
# -----------------------------------------------------------------------------
precompute_cells_for_scenario <- function(scenario, config) {
  stopifnot(scenario %in% 1:4)
  scenario <- as.integer(scenario)

  bands  <- names(config$BAND_RANGES)
  ratios <- as.numeric(config$RATIOS)
  ngrid  <- as.integer(config$N_REF_GRID_SIZE)

  per_band <- lapply(bands, function(b) {
    grid_n <- n_ref_grid(config$BAND_RANGES[[b]], ngrid)

    if (scenario == 1L) {
      # band x N_ref x ratio_V1
      df <- expand.grid(
        N_ref    = as.integer(grid_n),
        ratio_V1 = ratios,
        KEEP.OUT.ATTRS   = FALSE,
        stringsAsFactors = FALSE
      )
      df$band     <- b
      df$n_V1     <- as.integer(round(df$N_ref * df$ratio_V1))
      df$ratio_V2 <- NA_real_
      df$n_V2     <- NA_integer_
      df$ratio_V3 <- NA_real_
      df$n_V3     <- NA_integer_
    } else if (scenario == 2L) {
      # band x N_ref x ratio_V1 x ratio_V2
      df <- expand.grid(
        N_ref    = as.integer(grid_n),
        ratio_V1 = ratios,
        ratio_V2 = ratios,
        KEEP.OUT.ATTRS   = FALSE,
        stringsAsFactors = FALSE
      )
      df$band     <- b
      df$n_V1     <- as.integer(round(df$N_ref * df$ratio_V1))
      df$n_V2     <- as.integer(round(df$N_ref * df$ratio_V2))
      df$ratio_V3 <- NA_real_
      df$n_V3     <- NA_integer_
    } else {
      # Sc3 or Sc4: band x N_ref x ratio_V1 x ratio_V2 x ratio_V3
      df <- expand.grid(
        N_ref    = as.integer(grid_n),
        ratio_V1 = ratios,
        ratio_V2 = ratios,
        ratio_V3 = ratios,
        KEEP.OUT.ATTRS   = FALSE,
        stringsAsFactors = FALSE
      )
      df$band <- b
      df$n_V1 <- as.integer(round(df$N_ref * df$ratio_V1))
      df$n_V2 <- as.integer(round(df$N_ref * df$ratio_V2))
      df$n_V3 <- as.integer(round(df$N_ref * df$ratio_V3))
    }

    df[, c("band", "N_ref", "ratio_V1", "n_V1",
           "ratio_V2", "n_V2", "ratio_V3", "n_V3")]
  })

  cells <- do.call(rbind, per_band)

  # Ceiling
  cells$ceiling <- if (scenario == 1L) {
    cells$n_V1
  } else if (scenario == 2L) {
    pmax(cells$n_V1, cells$n_V2)
  } else {
    pmax(cells$n_V1, cells$n_V2, cells$n_V3)
  }

  cells$scenario <- scenario
  cells$cell_id  <- seq_len(nrow(cells))

  cells[, c("scenario", "cell_id", "band", "N_ref",
            "ratio_V1", "n_V1", "ratio_V2", "n_V2",
            "ratio_V3", "n_V3", "ceiling")]
}

# -----------------------------------------------------------------------------
# precompute_cells_for_all_scenarios
# -----------------------------------------------------------------------------
# Builds the cells table for every scenario and saves to disk. Skips any
# scenario whose RDS already exists, supporting resume-from-disk behaviour
# without re-derivation.
#
# Inputs
#   config : output of build_config().
#
# Side effects
#   Creates output_dir(config) if missing. Writes cells_sc<N>.rds for any
#   missing scenario.
# -----------------------------------------------------------------------------
precompute_cells_for_all_scenarios <- function(config) {
  dir.create(output_dir(config), recursive = TRUE, showWarnings = FALSE)

  for (s in 1:4) {
    cp <- cells_path(config, s)
    if (file.exists(cp)) {
      log_msg(sprintf("Cells for scenario %d already exist: %s", s, cp))
      next
    }
    cells <- precompute_cells_for_scenario(s, config)
    save_rds_atomic(cells, cp)
    log_msg(sprintf("Scenario %d: %d cells written to %s",
                    s, nrow(cells), cp))
  }
}

# -----------------------------------------------------------------------------
# stratified_robustness_cells
# -----------------------------------------------------------------------------
# Stratified-random selection of cells for the optional p_del robustness
# sweep. Selects up to config$ROBUSTNESS_CELLS_PER_BAND cells per (scenario,
# band) stratum, drawn uniformly at random from the scenario's cells table
# with a deterministic seed for reproducibility.
#
# Inputs
#   config : output of build_config().
#
# Returns: data frame with columns scenario, cell_id, band. Total row count
# is at most 4 scenarios x 3 bands x config$ROBUSTNESS_CELLS_PER_BAND.
# -----------------------------------------------------------------------------
stratified_robustness_cells <- function(config) {
  per_band <- as.integer(config$ROBUSTNESS_CELLS_PER_BAND)
  out_list <- list()
  rng_seed <- config$BASE_SEED + 999L
  i <- 0L
  for (s in 1:4) {
    cells <- readRDS(cells_path(config, s))
    for (b in names(config$BAND_RANGES)) {
      i <- i + 1L
      band_cells <- cells[cells$band == b, , drop = FALSE]
      if (nrow(band_cells) == 0L) next
      set.seed(rng_seed + i)
      k <- min(per_band, nrow(band_cells))
      pick <- sort(sample.int(nrow(band_cells), k))
      out_list[[length(out_list) + 1L]] <- data.frame(
        scenario = s,
        cell_id  = as.integer(band_cells$cell_id[pick]),
        band     = b,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out_list)
}
