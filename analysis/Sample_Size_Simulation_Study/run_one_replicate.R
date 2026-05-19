# =============================================================================
# run_one_replicate.R
# =============================================================================
#
# Per-(scenario, cell_id, rep_id) driver.
#
# Workflow:
#   1. Look up the cell's structural parameters.
#   2. Seed the RNG deterministically.
#   3. Draw the per-replicate biological parameters: mutation counts and
#      positions for V2 (Sc 2-4), V3 (Sc 3-4), and the emerging variant.
#   4. Build the population matrix at the ceiling row count.
#   5. Apply deleterious-mutation noise once at p_del to the full-ceiling
#      matrix. The noisy matrix is then reused across all sweep points
#      (the noise is a property of the simulated dataset, not regenerated
#      per detection test).
#   6. Run the 50-point log-spaced sweep over n_emerge in [2, ceiling].
#      Break at the first grid point that detects; no bisection refinement.
#   7. Write the result RDS atomically with full audit metadata.
#
# Returns the result invisibly. Designed to be called once per subprocess
# in parallel mode (run inside a fresh R session via callr::r_bg) or once
# per worker iteration in sequential mode.
#
# A separate entry point run_one_replicate_robustness() redirects the
# output to the p_del robustness sub-tree.
#
# Author : Vadim Tyuryaev
# =============================================================================

# Null-coalescing helper used during result construction. Defined at file
# scope so it is available throughout the sourced file regardless of where
# run_one_replicate() is invoked from.
`%||%` <- function(a, b) if (is.null(a)) b else a

# -----------------------------------------------------------------------------
# .draw_variants_for_scenario
# -----------------------------------------------------------------------------
# Internal helper. Builds the variants_spec list for a given scenario and
# replicate, using the cell's n_V_i values for established variants and
# cell$ceiling for the emerging variant.
#
# All variants force D614G as their first mutation. V2, V3, and the
# non-Omicron emerging variant draw their remaining positions from
# POOL_11_DRAW. The Omicron-like emerging variant in scenario 4 draws
# its 32 remaining positions from POOL_12_DRAW.
# -----------------------------------------------------------------------------
.draw_variants_for_scenario <- function(scenario, cell, empirical, ref_seq_int) {

  variants_spec <- list()

  # V1: D614G only. Always present (every scenario).
  variants_spec[["V1"]] <- list(
    positions   = 614L,
    substitutes = .AA_TO_INT[["G"]],
    n_sequences = as.integer(cell$n_V1)
  )

  # V2: present in Sc 2, 3, 4.
  if (scenario >= 2L) {
    n_v2_muts <- sample(empirical$MUT_COUNT_SET_11, 1L)
    variants_spec[["V2"]] <- draw_variant_spec(
      n_mutations  = as.integer(n_v2_muts),
      force_d614g  = TRUE,
      draw_pool    = empirical$POOL_11_DRAW,
      n_sequences  = as.integer(cell$n_V2),
      ref_seq_int  = ref_seq_int
    )
  }

  # V3: present in Sc 3, 4.
  if (scenario >= 3L) {
    n_v3_muts <- sample(empirical$MUT_COUNT_SET_11, 1L)
    variants_spec[["V3"]] <- draw_variant_spec(
      n_mutations  = as.integer(n_v3_muts),
      force_d614g  = TRUE,
      draw_pool    = empirical$POOL_11_DRAW,
      n_sequences  = as.integer(cell$n_V3),
      ref_seq_int  = ref_seq_int
    )
  }

  # Emerging variant. Sc 1, 2, 3 draw count from MUT_COUNT_SET_11 and
  # positions from POOL_11_DRAW. Sc 4 uses fixed Omicron-like 33 mutations
  # drawn from POOL_12_DRAW.
  if (scenario %in% 1:3) {
    n_emerge_muts <- sample(empirical$MUT_COUNT_SET_11, 1L)
    emerge_spec <- draw_variant_spec(
      n_mutations  = as.integer(n_emerge_muts),
      force_d614g  = TRUE,
      draw_pool    = empirical$POOL_11_DRAW,
      n_sequences  = as.integer(cell$ceiling),
      ref_seq_int  = ref_seq_int
    )
  } else {
    n_emerge_muts <- empirical$MUT_COUNT_OMICRON  # fixed at 33
    emerge_spec <- draw_variant_spec(
      n_mutations  = as.integer(n_emerge_muts),
      force_d614g  = TRUE,
      draw_pool    = empirical$POOL_12_DRAW,
      n_sequences  = as.integer(cell$ceiling),
      ref_seq_int  = ref_seq_int
    )
  }
  variants_spec[["Emerging"]] <- emerge_spec

  variants_spec
}

# -----------------------------------------------------------------------------
# run_one_replicate
# -----------------------------------------------------------------------------
run_one_replicate <- function(scenario,
                              cell_id,
                              rep_id,
                              cells,
                              ref_seq_int,
                              empirical,
                              config,
                              p_del_override = NULL,
                              output_base_dir = NULL) {

  t0 <- Sys.time()

  scenario <- as.integer(scenario)
  cell_id  <- as.integer(cell_id)
  rep_id   <- as.integer(rep_id)

  # ----- 1. Cell lookup ----------------------------------------------------
  cell <- cells[cells$cell_id == cell_id, ]
  if (nrow(cell) != 1L)
    stop(sprintf("Cell %d not found in cells table for scenario %d.",
                 cell_id, scenario), call. = FALSE)

  # ----- 2. Seed -----------------------------------------------------------
  rep_seed <- seed_for_replicate(config$BASE_SEED, scenario, cell_id, rep_id)
  set.seed(rep_seed)

  # ----- 3. Build variants_spec -------------------------------------------
  variants_spec    <- .draw_variants_for_scenario(scenario, cell,
                                                  empirical, ref_seq_int)
  emerge_positions <- variants_spec[["Emerging"]]$positions
  mutation_count   <- length(emerge_positions)

  # ----- 4. Build population matrix at the ceiling ------------------------
  pop <- simulate_population_snapshot(
    ref_seq_int   = ref_seq_int,
    n_ref         = as.integer(cell$N_ref),
    variants_spec = variants_spec
  )
  mat <- pop$matrix
  rm(pop)

  # ----- 5. Apply deleterious-mutation noise once -------------------------
  # The noise is a property of the simulated dataset. We apply it once on
  # the full-ceiling matrix; subsequent sweep points reuse the same noisy
  # matrix via a row prefix. Seed offset distinguishes noise RNG from
  # variant-draw RNG so the same replicate at different p_del values gets
  # the same variant placement but different (independent) noise patterns.
  p_del_value <- if (is.null(p_del_override)) config$P_DEL else p_del_override
  mat <- apply_deleterious_noise(mat,
                                 p_del = p_del_value,
                                 seed  = rep_seed + 1L)

  n_base <- as.integer(cell$N_ref) +
            as.integer(cell$n_V1) +
            (if (scenario >= 2L) as.integer(cell$n_V2) else 0L) +
            (if (scenario >= 3L) as.integer(cell$n_V3) else 0L)

  ceiling_val <- as.integer(cell$ceiling)

  # ----- 6. 50-point log-spaced sweep, no refinement ---------------------
  grid <- unique(as.integer(round(exp(
    seq(log(2), log(ceiling_val), length.out = config$GRID_SIZE)
  ))))

  sweep_outcomes <- logical(length(grid))
  sweep_modelName <- character(length(grid))
  sweep_G         <- integer(length(grid))
  sweep_n_class_1 <- integer(length(grid))
  first_hit_idx   <- NA_integer_

  for (k in seq_along(grid)) {
    n_emerge <- grid[k]
    fit <- test_detection(
      mat              = mat,
      n_rows           = n_base + n_emerge,
      emerge_positions = emerge_positions,
      mclust_models    = config$MCLUST_MODELS,
      mclust_G         = config$MCLUST_G
    )
    sweep_outcomes[k]  <- fit$detected
    sweep_modelName[k] <- fit$modelName_chosen %||% NA_character_
    sweep_G[k]         <- fit$G_chosen         %||% NA_integer_
    sweep_n_class_1[k] <- fit$n_class_1_sites
    if (fit$detected) {
      first_hit_idx <- k
      break
    }
  }

  # Trim sweep arrays to executed points only
  executed_len <- if (is.na(first_hit_idx)) length(grid) else first_hit_idx
  sweep_outcomes  <- sweep_outcomes[seq_len(executed_len)]
  sweep_modelName <- sweep_modelName[seq_len(executed_len)]
  sweep_G         <- sweep_G[seq_len(executed_len)]
  sweep_n_class_1 <- sweep_n_class_1[seq_len(executed_len)]
  sweep_grid_exec <- grid[seq_len(executed_len)]

  if (is.na(first_hit_idx)) {
    n_emerge_needed <- NA_integer_
    detected        <- FALSE
    hit_modelName   <- NA_character_
    hit_G           <- NA_integer_
  } else {
    n_emerge_needed <- grid[first_hit_idx]
    detected        <- TRUE
    hit_modelName   <- sweep_modelName[first_hit_idx]
    hit_G           <- sweep_G[first_hit_idx]
  }

  walltime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # ----- 7. Build result record -------------------------------------------
  result <- list(
    scenario           = scenario,
    cell_id            = cell_id,
    rep_id             = rep_id,
    band               = as.character(cell$band),
    N_ref              = as.integer(cell$N_ref),
    n_V1               = as.integer(cell$n_V1),
    n_V2               = if (scenario >= 2L) as.integer(cell$n_V2)
                         else NA_integer_,
    n_V3               = if (scenario >= 3L) as.integer(cell$n_V3)
                         else NA_integer_,
    ratio_V1           = as.numeric(cell$ratio_V1),
    ratio_V2           = if (scenario >= 2L) as.numeric(cell$ratio_V2)
                         else NA_real_,
    ratio_V3           = if (scenario >= 3L) as.numeric(cell$ratio_V3)
                         else NA_real_,
    ceiling            = ceiling_val,
    p_del              = p_del_value,
    mutation_count     = mutation_count,
    mutation_positions = as.integer(emerge_positions),
    n_emerge_needed    = n_emerge_needed,
    detected           = detected,
    hit_modelName      = hit_modelName,
    hit_G              = hit_G,
    sweep_grid         = sweep_grid_exec,
    sweep_outcomes     = sweep_outcomes,
    sweep_modelName    = sweep_modelName,
    sweep_G            = sweep_G,
    sweep_n_class_1    = sweep_n_class_1,
    seed               = rep_seed,
    walltime_s         = walltime_s,
    r_version          = R.version.string,
    package_version    = as.character(utils::packageVersion("ViralEntropR"))
  )

  # ----- 8. Atomic write --------------------------------------------------
  base_dir <- if (is.null(output_base_dir)) output_dir(config) else output_base_dir
  out_path <- replicate_path(config, scenario, cell_id, rep_id,
                             base_dir = base_dir)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  save_rds_atomic(result, out_path)

  # ----- 9. Free memory ---------------------------------------------------
  rm(mat)
  gc(verbose = FALSE)

  invisible(result)
}

# -----------------------------------------------------------------------------
# run_one_replicate_robustness
# -----------------------------------------------------------------------------
# Convenience wrapper for the optional p_del robustness sweep. Writes the
# replicate's RDS under outputs/robustness_pdel/p_del_<value>/replicates_sc<N>/
# rather than the main outputs/replicates_sc<N>/.
# -----------------------------------------------------------------------------
run_one_replicate_robustness <- function(scenario, cell_id, rep_id,
                                         cells, ref_seq_int, empirical,
                                         config, p_del_value) {
  run_one_replicate(
    scenario        = scenario,
    cell_id         = cell_id,
    rep_id          = rep_id,
    cells           = cells,
    ref_seq_int     = ref_seq_int,
    empirical       = empirical,
    config          = config,
    p_del_override  = p_del_value,
    output_base_dir = robustness_output_dir(config, p_del_value)
  )
}
