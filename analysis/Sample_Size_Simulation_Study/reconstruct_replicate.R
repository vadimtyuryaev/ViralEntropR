# =============================================================================
# reconstruct_replicate.R
# =============================================================================
#
# Deterministic reconstruction of per-replicate state from the on-disk
# summary plus the reference FASTA and the empirical constants.
#
# Public entry points:
#   reconstruct_replicate_matrix()    rebuild one replicate's matrix
#   reconstruct_replicate_pipeline()  replay the full detection pipeline
#                                     on one replicate
#   list_na_replicates()              enumerate replicates with NA
#                                     n_emerge_needed across scenarios
#   reconstruct_na_replicates()       batch driver: loop the per-replicate
#                                     reconstruction over a list of
#                                     replicates (e.g. the output of
#                                     list_na_replicates())
#
# Internal helper:
#   .resolve_replicate_record()       locate and load the on-disk record
#                                     for one replicate, with fallback
#                                     from the per-replicate RDS to the
#                                     scenario summary RDS
#
# Requires (source first, in this order):
#   setup.R, simulator.R, detect_in_snapshot.R, run_one_replicate.R
#
# Author : Vadim Tyuryaev
# =============================================================================

# -----------------------------------------------------------------------------
# .resolve_replicate_record
# -----------------------------------------------------------------------------
# Locates the per-replicate RDS for (scenario, cell_id, rep_id), loads it,
# and returns the result object. Falls back to the summary RDS if the
# per-replicate file is missing.
# -----------------------------------------------------------------------------
.resolve_replicate_record <- function(scenario, cell_id, rep_id, config,
                                      base_dir = output_dir(config)) {
  rep_path <- replicate_path(config, scenario, cell_id, rep_id,
                             base_dir = base_dir)
  if (file.exists(rep_path)) {
    return(list(record = readRDS(rep_path), source = "per_replicate_rds"))
  }
  summ_path <- summary_path(config, scenario)
  if (file.exists(summ_path)) {
    df <- readRDS(summ_path)
    row <- df[df$cell_id == cell_id & df$rep_id == rep_id, , drop = FALSE]
    if (nrow(row) == 1L) {
      return(list(record = as.list(row), source = "summary_rds"))
    }
  }
  stop(sprintf(
    "No record found for scenario=%d, cell_id=%d, rep_id=%d.",
    scenario, cell_id, rep_id), call. = FALSE)
}

# -----------------------------------------------------------------------------
# reconstruct_replicate_matrix
# -----------------------------------------------------------------------------
# Rebuilds the integer feature matrix for one replicate at any chosen
# n_emerge in 0..ceiling. Variant draws and deleterious-noise pattern
# are reseeded with the same scheme the orchestrator uses, so the
# returned matrix is bit-identical to the production run at the recorded
# detection point.
#
# Returns a list with $matrix, $variants_spec, $emerge_positions,
# $base_rows, $n_emerge, $ceiling, $seed, $noise_seed, $p_del,
# $noise_applied, $record_source.
# -----------------------------------------------------------------------------
reconstruct_replicate_matrix <- function(scenario, cell_id, rep_id,
                                         cells       = NULL,
                                         ref_seq_int = NULL,
                                         empirical   = NULL,
                                         config      = NULL,
                                         n_emerge    = NULL,
                                         apply_noise = TRUE) {
  
  if (is.null(config)) {
    config <- build_config()
    config$STUDY_DIR <- resolve_study_dir()
  }
  if (is.null(ref_seq_int)) {
    ref_seq_int <- load_reference_sequence(config$REF_SEQ_FASTA)
  }
  if (is.null(empirical)) {
    empirical <- build_empirical_constants()
  }
  if (is.null(cells)) {
    cells <- readRDS(cells_path(config, scenario))
  }
  
  cell <- cells[cells$cell_id == cell_id, ]
  if (nrow(cell) != 1L)
    stop(sprintf("Cell %d not found in cells table for scenario %d.",
                 cell_id, scenario), call. = FALSE)
  
  rec_info  <- .resolve_replicate_record(scenario, cell_id, rep_id, config)
  rec       <- rec_info$record
  rep_seed  <- as.integer(rec$seed)
  p_del     <- as.numeric(rec$p_del)
  ceiling_v <- as.integer(cell$ceiling)
  
  if (is.null(n_emerge)) {
    n_emerge <- rec$n_emerge_needed
    if (is.na(n_emerge)) n_emerge <- ceiling_v
  }
  n_emerge <- as.integer(n_emerge)
  if (n_emerge < 0L || n_emerge > ceiling_v)
    stop(sprintf("`n_emerge` must be in 0..ceiling (%d); got %d.",
                 ceiling_v, n_emerge), call. = FALSE)
  
  set.seed(rep_seed)
  variants_spec    <- .draw_variants_for_scenario(scenario, cell,
                                                  empirical, ref_seq_int)
  emerge_positions <- variants_spec[["Emerging"]]$positions
  
  pop <- simulate_population_snapshot(
    ref_seq_int   = ref_seq_int,
    n_ref         = as.integer(cell$N_ref),
    variants_spec = variants_spec
  )
  mat <- pop$matrix
  
  noise_seed <- rep_seed + 1L
  if (apply_noise) {
    mat <- apply_deleterious_noise(mat, p_del = p_del, seed = noise_seed)
  }
  
  n_base <- as.integer(cell$N_ref) +
    as.integer(cell$n_V1) +
    (if (scenario >= 2L) as.integer(cell$n_V2) else 0L) +
    (if (scenario >= 3L) as.integer(cell$n_V3) else 0L)
  total_rows <- n_base + n_emerge
  
  list(
    matrix           = mat[seq_len(total_rows), , drop = FALSE],
    variants_spec    = variants_spec,
    emerge_positions = emerge_positions,
    base_rows        = n_base,
    n_emerge         = n_emerge,
    ceiling          = ceiling_v,
    seed             = rep_seed,
    noise_seed       = noise_seed,
    p_del            = p_del,
    noise_applied    = apply_noise,
    record_source    = rec_info$source
  )
}

# -----------------------------------------------------------------------------
# reconstruct_replicate_pipeline
# -----------------------------------------------------------------------------
# Replays the full detection pipeline at one sweep step. Returns the
# matrix plus the per-site entropy vector, the cluster_sites_by_entropy
# fit, the relabelled DataFrame (when the multi-class path was taken),
# the detection-path label, and the test_detection() output.
# -----------------------------------------------------------------------------
reconstruct_replicate_pipeline <- function(scenario, cell_id, rep_id,
                                           cells       = NULL,
                                           ref_seq_int = NULL,
                                           empirical   = NULL,
                                           config      = NULL,
                                           n_emerge    = NULL) {
  
  if (is.null(config)) {
    config <- build_config()
    config$STUDY_DIR <- resolve_study_dir()
  }
  
  recon <- reconstruct_replicate_matrix(
    scenario    = scenario, cell_id = cell_id, rep_id = rep_id,
    cells       = cells,    ref_seq_int = ref_seq_int,
    empirical   = empirical, config = config,
    n_emerge    = n_emerge,  apply_noise = TRUE
  )
  
  mat        <- recon$matrix
  n_rows     <- nrow(mat)
  entropies  <- compute_entropy_per_site(mat, n_rows = n_rows,
                                         aa_levels = 25L)
  
  cluster_fit <- .call_cluster_sites(
    entropies     = entropies,
    n_rows        = n_rows,
    mclust_models = config$MCLUST_MODELS,
    mclust_G      = config$MCLUST_G
  )
  
  detection <- test_detection(
    mat              = mat,
    n_rows           = n_rows,
    emerge_positions = recon$emerge_positions,
    mclust_models    = config$MCLUST_MODELS,
    mclust_G         = config$MCLUST_G
  )
  
  df <- cluster_fit$DataFrame
  if (all(entropies[recon$emerge_positions] == 0)) {
    detection_path <- "all_zero_entropy"
    relabelled_df  <- NULL
  } else if (nrow(df) == 0L) {
    detection_path <- "empty_dataframe"
    relabelled_df  <- NULL
  } else if (all(df$class == 999L)) {
    detection_path <- "sentinel_999"
    relabelled_df  <- NULL
  } else if (length(unique(df$class)) == 1L) {
    detection_path <- "G1_collapse"
    relabelled_df  <- NULL
  } else {
    detection_path <- "multi_class"
    relabelled_df  <- ViralEntropR::relabel_entropy_classes(df)
  }
  
  c(recon, list(
    entropies      = entropies,
    cluster_fit    = cluster_fit,
    relabelled_df  = relabelled_df,
    detection_path = detection_path,
    detection      = detection
  ))
}

# -----------------------------------------------------------------------------
# list_na_replicates
# -----------------------------------------------------------------------------
# Returns a data frame of all replicates whose sweep failed to detect
# (n_emerge_needed is NA) across the requested scenarios. Useful as the
# work list for batch reconstruction.
#
# Inputs
#   config    : output of build_config() with STUDY_DIR set. Pass NULL to
#               auto-build.
#   scenarios : integer vector of scenarios to scan. Default 1:4.
#   bands     : optional character vector to filter by band
#               ("small", "medium", "large"). Default NULL (all bands).
#   base_dir  : output directory to read summaries from. Defaults to
#               output_dir(config); pass robustness_output_dir(config, p)
#               to scan a robustness-sweep sub-tree instead.
# -----------------------------------------------------------------------------
list_na_replicates <- function(config    = NULL,
                               scenarios = 1:4,
                               bands     = NULL,
                               base_dir  = NULL) {
  
  if (is.null(config)) {
    config <- build_config()
    config$STUDY_DIR <- resolve_study_dir()
  }
  if (is.null(base_dir)) base_dir <- output_dir(config)
  
  out <- vector("list", length(scenarios))
  for (k in seq_along(scenarios)) {
    s  <- as.integer(scenarios[k])
    sp <- file.path(base_dir, sprintf("summary_sc%d.rds", s))
    if (!file.exists(sp)) next
    df <- readRDS(sp)
    df <- df[is.na(df$n_emerge_needed), , drop = FALSE]
    if (!is.null(bands)) df <- df[df$band %in% bands, , drop = FALSE]
    if (nrow(df) == 0L) next
    
    keep <- intersect(
      c("scenario", "cell_id", "rep_id", "band", "N_ref",
        "n_V1", "n_V2", "n_V3", "ceiling", "mutation_count",
        "seed", "hit_modelName", "hit_G"),
      colnames(df)
    )
    out[[k]] <- df[, keep, drop = FALSE]
  }
  out <- out[!vapply(out, is.null, logical(1L))]
  if (length(out) == 0L) {
    return(data.frame(
      scenario = integer(0), cell_id = integer(0), rep_id = integer(0),
      band = character(0), stringsAsFactors = FALSE
    ))
  }
  
  res <- do.call(rbind, out)
  band_order <- factor(res$band, levels = c("small", "medium", "large"))
  res <- res[order(res$scenario, band_order, res$cell_id, res$rep_id), ,
             drop = FALSE]
  rownames(res) <- NULL
  res
}

# -----------------------------------------------------------------------------
# reconstruct_na_replicates
# -----------------------------------------------------------------------------
# Reconstructs the matrix-and-pipeline state for every row of `na_table`.
#
# `what`:
#   "pipeline" : full reconstruct_replicate_pipeline() per row (heavy)
#   "matrix"   : reconstruct_replicate_matrix() per row (no GMM replay)
#   "summary"  : tidy data frame summarising the detection path and the
#                emerging-site class assignments per row (lightweight)
#
# `workers > 1L` uses parallel::mclapply (forking; Linux/macOS only).
# Set workers = 1L on Windows.
# -----------------------------------------------------------------------------
reconstruct_na_replicates <- function(na_table,
                                      config  = NULL,
                                      what    = c("pipeline", "matrix",
                                                  "summary"),
                                      workers = 1L,
                                      verbose = TRUE) {
  what    <- match.arg(what)
  workers <- as.integer(workers)
  
  if (nrow(na_table) == 0L) {
    if (what == "summary") return(data.frame())
    return(list())
  }
  
  if (is.null(config)) {
    config <- build_config()
    config$STUDY_DIR <- resolve_study_dir()
  }
  
  ref_seq_int <- load_reference_sequence(config$REF_SEQ_FASTA)
  empirical   <- build_empirical_constants()
  cells_cache <- new.env(parent = emptyenv())
  get_cells <- function(s) {
    key <- as.character(s)
    if (!exists(key, envir = cells_cache, inherits = FALSE)) {
      assign(key, readRDS(cells_path(config, s)), envir = cells_cache)
    }
    get(key, envir = cells_cache, inherits = FALSE)
  }
  
  worker_fn <- function(i) {
    row <- na_table[i, , drop = FALSE]
    s   <- as.integer(row$scenario)
    out <- if (what == "matrix") {
      reconstruct_replicate_matrix(
        scenario    = s, cell_id = as.integer(row$cell_id),
        rep_id      = as.integer(row$rep_id),
        cells       = get_cells(s),
        ref_seq_int = ref_seq_int, empirical = empirical,
        config      = config
      )
    } else {
      reconstruct_replicate_pipeline(
        scenario    = s, cell_id = as.integer(row$cell_id),
        rep_id      = as.integer(row$rep_id),
        cells       = get_cells(s),
        ref_seq_int = ref_seq_int, empirical = empirical,
        config      = config
      )
    }
    if (verbose && i %% 10L == 0L) {
      message(sprintf("  reconstructed %d / %d", i, nrow(na_table)))
    }
    out
  }
  
  if (workers > 1L) {
    if (!requireNamespace("parallel", quietly = TRUE))
      stop("Package 'parallel' is required for workers > 1L.",
           call. = FALSE)
    results <- parallel::mclapply(seq_len(nrow(na_table)),
                                  worker_fn, mc.cores = workers)
  } else {
    results <- lapply(seq_len(nrow(na_table)), worker_fn)
  }
  
  names(results) <- sprintf("sc%d_cell%04d_rep%02d",
                            na_table$scenario,
                            na_table$cell_id, na_table$rep_id)
  
  if (what == "summary") {
    rows <- lapply(seq_along(results), function(i) {
      r <- results[[i]]
      classes_emerge <- r$detection$class_at_sites
      data.frame(
        scenario                     = na_table$scenario[i],
        cell_id                      = na_table$cell_id[i],
        rep_id                       = na_table$rep_id[i],
        band                         = na_table$band[i],
        ceiling                      = na_table$ceiling[i],
        mutation_count               = length(r$emerge_positions),
        n_emerge_used                = r$n_emerge,
        detection_path               = r$detection_path,
        modelName_chosen             = r$detection$modelName_chosen %||% NA,
        G_chosen                     = r$detection$G_chosen         %||% NA,
        n_class_1_sites              = r$detection$n_class_1_sites,
        entropy_at_emerge_min        = min(r$entropies[r$emerge_positions]),
        entropy_at_emerge_max        = max(r$entropies[r$emerge_positions]),
        entropy_at_emerge_median     = stats::median(
          r$entropies[r$emerge_positions]),
        n_emerging_in_class_1        = sum(classes_emerge == 1L, na.rm = TRUE),
        n_emerging_in_higher_classes = sum(classes_emerge > 1L,
                                           na.rm = TRUE),
        n_emerging_missing           = sum(is.na(classes_emerge)),
        stringsAsFactors             = FALSE,
        check.names                  = FALSE
      )
    })
    return(do.call(rbind, rows))
  }
  
  results
}

# -----------------------------------------------------------------------------
# reconstruct_noise_pattern
# -----------------------------------------------------------------------------
# Returns the set of (row, col) cells overwritten by the deleterious-noise
# step of one replicate, together with the original residue (after variant
# substitution, before noise) and the noise-substituted residue.
#
# Strategy: build the matrix twice via reconstruct_replicate_matrix(),
# once with apply_noise = FALSE (the clean variant matrix) and once with
# apply_noise = TRUE (the production matrix). Cell-wise difference gives
# the noise mask exactly.
#
# Inputs match reconstruct_replicate_matrix(). The default n_emerge is
# the cell's full ceiling (so the returned mask covers every noise event
# the replicate ever encountered, not just the rows below the detection
# point).
#
# Returns a list with
#   $mask         : logical matrix of nrow × 1273, TRUE where noise hit.
#   $events       : data frame of noise events, one row per affected cell.
#                   Columns: row (matrix row index), col (Spike position),
#                   from (original code), to (post-noise code),
#                   from_aa (single-letter), to_aa (single-letter),
#                   row_origin ("reference", "V1", "V2", "V3", "emerging").
#   $by_row       : integer vector of length nrow, count of noise events
#                   per row.
#   $by_position  : integer vector of length 1273, count of noise events
#                   per Spike position.
#   $by_row_origin: named integer vector summarising events per variant
#                   block (reference, V1, V2, V3, emerging).
#   $matrix_clean : the variant-only matrix (before noise) for reference.
#   $matrix_noisy : the production matrix (variant + noise).
# -----------------------------------------------------------------------------
reconstruct_noise_pattern <- function(scenario, cell_id, rep_id,
                                      cells       = NULL,
                                      ref_seq_int = NULL,
                                      empirical   = NULL,
                                      config      = NULL,
                                      n_emerge    = NULL) {
  
  if (is.null(config)) {
    config <- build_config()
    config$STUDY_DIR <- resolve_study_dir()
  }
  if (is.null(ref_seq_int)) {
    ref_seq_int <- load_reference_sequence(config$REF_SEQ_FASTA)
  }
  if (is.null(empirical)) {
    empirical <- build_empirical_constants()
  }
  if (is.null(cells)) {
    cells <- readRDS(cells_path(config, scenario))
  }
  
  # Default to the cell's full ceiling so the returned mask covers the
  # whole matrix the replicate ever held in memory.
  if (is.null(n_emerge)) {
    cell <- cells[cells$cell_id == cell_id, ]
    if (nrow(cell) != 1L)
      stop(sprintf("Cell %d not found in scenario %d.", cell_id, scenario),
           call. = FALSE)
    n_emerge <- as.integer(cell$ceiling)
  }
  
  recon_clean <- reconstruct_replicate_matrix(
    scenario = scenario, cell_id = cell_id, rep_id = rep_id,
    cells = cells, ref_seq_int = ref_seq_int, empirical = empirical,
    config = config, n_emerge = n_emerge, apply_noise = FALSE
  )
  recon_noisy <- reconstruct_replicate_matrix(
    scenario = scenario, cell_id = cell_id, rep_id = rep_id,
    cells = cells, ref_seq_int = ref_seq_int, empirical = empirical,
    config = config, n_emerge = n_emerge, apply_noise = TRUE
  )
  
  clean <- recon_clean$matrix
  noisy <- recon_noisy$matrix
  mask  <- clean != noisy
  
  # Row origin labels
  cell <- cells[cells$cell_id == cell_id, ]
  origin_breaks <- c(
    as.integer(cell$N_ref),
    as.integer(cell$N_ref) + as.integer(cell$n_V1),
    as.integer(cell$N_ref) + as.integer(cell$n_V1) +
      (if (scenario >= 2L) as.integer(cell$n_V2) else 0L),
    as.integer(cell$N_ref) + as.integer(cell$n_V1) +
      (if (scenario >= 2L) as.integer(cell$n_V2) else 0L) +
      (if (scenario >= 3L) as.integer(cell$n_V3) else 0L)
  )
  origin_labels <- c("reference", "V1", "V2", "V3", "emerging")
  
  row_origin <- character(nrow(clean))
  for (i in seq_len(nrow(clean))) {
    if      (i <= origin_breaks[1]) row_origin[i] <- origin_labels[1]
    else if (i <= origin_breaks[2]) row_origin[i] <- origin_labels[2]
    else if (i <= origin_breaks[3]) row_origin[i] <- origin_labels[3]
    else if (i <= origin_breaks[4]) row_origin[i] <- origin_labels[4]
    else                             row_origin[i] <- origin_labels[5]
  }
  
  # Letter map: integer code -> single-letter via decode_aa_sequence
  aa_letters <- as.vector(ViralEntropR::decode_aa_sequence(
    matrix(seq_len(25L), ncol = 1L)
  ))
  
  hit_rc <- which(mask, arr.ind = TRUE)
  if (nrow(hit_rc) > 0L) {
    from_code <- clean[mask]
    to_code   <- noisy[mask]
    events <- data.frame(
      row        = as.integer(hit_rc[, "row"]),
      col        = as.integer(hit_rc[, "col"]),
      from       = as.integer(from_code),
      to         = as.integer(to_code),
      from_aa    = aa_letters[from_code],
      to_aa      = aa_letters[to_code],
      row_origin = row_origin[hit_rc[, "row"]],
      stringsAsFactors = FALSE,
      check.names      = FALSE
    )
    events <- events[order(events$row, events$col), , drop = FALSE]
    rownames(events) <- NULL
  } else {
    events <- data.frame(
      row = integer(0), col = integer(0),
      from = integer(0), to = integer(0),
      from_aa = character(0), to_aa = character(0),
      row_origin = character(0),
      stringsAsFactors = FALSE
    )
  }
  
  list(
    mask           = mask,
    events         = events,
    by_row         = rowSums(mask),
    by_position    = colSums(mask),
    by_row_origin  = table(factor(row_origin[hit_rc[, "row"]],
                                  levels = origin_labels)),
    matrix_clean   = clean,
    matrix_noisy   = noisy
  )
}

# -----------------------------------------------------------------------------
# reconstruct_noise_patterns_batch
# -----------------------------------------------------------------------------
# Batch driver around reconstruct_noise_pattern(). Consumes the output of
# list_na_replicates() (or any equivalent data frame with scenario,
# cell_id, rep_id columns) and runs the per-replicate noise-pattern
# reconstruction over every row.
#
# Two output modes:
#   what = "events"   : returns one tidy data frame with all events from
#                       all replicates, augmented with the (scenario,
#                       cell_id, rep_id) triple per row. Best for plotting
#                       and group-by analysis across many replicates.
#   what = "full"     : returns a named list of full reconstruct_noise_pattern()
#                       outputs (mask, events, by_row, by_position,
#                       by_row_origin, matrix_clean, matrix_noisy) keyed
#                       by "sc<S>_cell<CCCC>_rep<RR>". Heavy — each
#                       large-band replicate carries two ~1 GB matrices.
#                       Use only on small subsets.
#   what = "summary"  : returns a data frame with one row per replicate
#                       containing only the aggregate event counts
#                       (n_events, events_per_origin), discarding the
#                       per-event detail. Cheapest for cross-replicate
#                       overview.
#
# `workers > 1L` uses parallel::mclapply (Linux/macOS only). Each
# reconstruction in the large band materialises two matrices of up to
# ~2 GB each, so worker count is constrained by RAM, not just CPU. On a
# 128 GB server, workers = 8L is comfortable for large-band cells; drop
# to 4L if you expect many large-band reconstructions in flight.
# -----------------------------------------------------------------------------
reconstruct_noise_patterns_batch <- function(replicates_table,
                                             config  = NULL,
                                             what    = c("events", "full",
                                                         "summary"),
                                             workers = 1L,
                                             verbose = TRUE) {
  what    <- match.arg(what)
  workers <- as.integer(workers)
  
  if (nrow(replicates_table) == 0L) {
    if (what == "events" || what == "summary") return(data.frame())
    return(list())
  }
  
  if (is.null(config)) {
    config <- build_config()
    config$STUDY_DIR <- resolve_study_dir()
  }
  
  ref_seq_int <- load_reference_sequence(config$REF_SEQ_FASTA)
  empirical   <- build_empirical_constants()
  cells_cache <- new.env(parent = emptyenv())
  get_cells <- function(s) {
    key <- as.character(s)
    if (!exists(key, envir = cells_cache, inherits = FALSE)) {
      assign(key, readRDS(cells_path(config, s)), envir = cells_cache)
    }
    get(key, envir = cells_cache, inherits = FALSE)
  }
  
  worker_fn <- function(i) {
    row <- replicates_table[i, , drop = FALSE]
    s   <- as.integer(row$scenario)
    np  <- reconstruct_noise_pattern(
      scenario    = s,
      cell_id     = as.integer(row$cell_id),
      rep_id      = as.integer(row$rep_id),
      cells       = get_cells(s),
      ref_seq_int = ref_seq_int,
      empirical   = empirical,
      config      = config
    )
    if (verbose && i %% 10L == 0L) {
      message(sprintf("  noise pattern %d / %d  (%s)",
                      i, nrow(replicates_table),
                      sprintf("sc%d_cell%04d_rep%02d",
                              s, row$cell_id, row$rep_id)))
    }
    
    if (what == "events") {
      if (nrow(np$events) == 0L) return(NULL)
      ev <- np$events
      ev$scenario <- s
      ev$cell_id  <- as.integer(row$cell_id)
      ev$rep_id   <- as.integer(row$rep_id)
      ev[, c("scenario", "cell_id", "rep_id",
             "row", "col", "from", "to",
             "from_aa", "to_aa", "row_origin")]
    } else if (what == "summary") {
      ot <- np$by_row_origin
      data.frame(
        scenario             = s,
        cell_id              = as.integer(row$cell_id),
        rep_id               = as.integer(row$rep_id),
        band                 = row$band,
        ceiling              = row$ceiling,
        n_events_total       = nrow(np$events),
        n_events_reference   = as.integer(ot["reference"] %||% 0L),
        n_events_V1          = as.integer(ot["V1"]        %||% 0L),
        n_events_V2          = as.integer(ot["V2"]        %||% 0L),
        n_events_V3          = as.integer(ot["V3"]        %||% 0L),
        n_events_emerging    = as.integer(ot["emerging"]  %||% 0L),
        n_positions_hit      = sum(np$by_position > 0L),
        max_events_at_one_pos = max(np$by_position),
        stringsAsFactors     = FALSE,
        check.names          = FALSE
      )
    } else {
      np  # full
    }
  }
  
  if (workers > 1L) {
    if (!requireNamespace("parallel", quietly = TRUE))
      stop("Package 'parallel' is required for workers > 1L.",
           call. = FALSE)
    results <- parallel::mclapply(seq_len(nrow(replicates_table)),
                                  worker_fn, mc.cores = workers)
  } else {
    results <- lapply(seq_len(nrow(replicates_table)), worker_fn)
  }
  
  if (what == "events") {
    results <- results[!vapply(results, is.null, logical(1L))]
    if (length(results) == 0L) return(data.frame())
    return(do.call(rbind, results))
  }
  
  if (what == "summary") {
    return(do.call(rbind, results))
  }
  
  # Full mode: keyed list
  names(results) <- sprintf("sc%d_cell%04d_rep%02d",
                            replicates_table$scenario,
                            replicates_table$cell_id,
                            replicates_table$rep_id)
  results
}