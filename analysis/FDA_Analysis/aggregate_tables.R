################################################################################
## aggregate_tables.R — FDA_Analysis (v2)
##
## Reads all outputs/cells/cell_*.rds and produces six supplementary CSV
## tables under outputs/tables/:
##
##   tab_01_internal_quality.csv  per (cell × metric):
##                                 n_class1, gmm_G, k_selected, mean si,
##                                 metric_status, convergence_*_frame
##   tab_02_cluster_membership.csv per (cell × metric × site)
##   tab_03_stability.csv         per (cell × metric × cluster_id) Hennig
##                                 bootstrap Jaccard
##   tab_04_cross_strategy_ari.csv per (dataset × variant × metric × strategy
##                                 pair) — ARI on common class-1 sites
##   tab_05_snp_enrichment.csv    per (cell × metric) Fisher exact result
##   tab_06_per_frame.csv         per (cell × metric × frame) — frame
##                                 diagnostics for the animations
##
## Three metrics in v2: entropy, hellinger_T1, hellinger_Tpred.
##
## Invocation: Rscript aggregate_tables.R
################################################################################

suppressPackageStartupMessages({
  source("setup.R")
  source("helpers_fda.R")
  .have_mclust_ari <- requireNamespace("mclust", quietly = TRUE)
})

if (!.have_mclust_ari)
  stop("Package 'mclust' is required for adjustedRandIndex.", call. = FALSE)

config <- build_config()
setwd(config$STUDY_DIR)
dir.create(tables_dir(config), recursive = TRUE, showWarnings = FALSE)

METRICS <- c("entropy", "hellinger_T1", "hellinger_Tpred")


# ── Atomic CSV writer ────────────────────────────────────────────────────────

write_csv_atomic <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid(), "-",
                 sample.int(.Machine$integer.max, 1L))
  utils::write.csv(df, tmp, row.names = FALSE, fileEncoding = "UTF-8")
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    file.remove(tmp)
  }
  invisible(path)
}


# ── Helpers ──────────────────────────────────────────────────────────────────

fda_block <- function(cell, metric) {
  switch(metric,
         entropy         = cell$fda_entropy,
         hellinger_T1    = cell$fda_hellinger_T1,
         hellinger_Tpred = cell$fda_hellinger_Tpred)
}

stability_block <- function(cell, metric) {
  switch(metric,
         entropy         = cell$stability_entropy,
         hellinger_T1    = cell$stability_hellinger_T1,
         hellinger_Tpred = cell$stability_hellinger_Tpred)
}

enrichment_block <- function(cell, metric) {
  switch(metric,
         entropy         = cell$snp_enrichment_entropy,
         hellinger_T1    = cell$snp_enrichment_hellinger_T1,
         hellinger_Tpred = cell$snp_enrichment_hellinger_Tpred)
}

frames_block <- function(cell, metric) {
  switch(metric,
         entropy         = cell$frames_entropy,
         hellinger_T1    = cell$frames_hellinger_T1,
         hellinger_Tpred = cell$frames_hellinger_Tpred)
}

convergence_block <- function(cell, metric) {
  switch(metric,
         entropy         = cell$convergence_entropy,
         hellinger_T1    = cell$convergence_hellinger_T1,
         hellinger_Tpred = cell$convergence_hellinger_Tpred)
}


# ── Load cells ───────────────────────────────────────────────────────────────

files <- list.files(cells_dir(config), pattern = "^cell_.+\\.rds$",
                     full.names = TRUE)
log_msg(sprintf("aggregate_tables: %d cell RDS files.", length(files)))

cells <- lapply(files, function(f)
  tryCatch(readRDS(f), error = function(e) {
    log_error(config, sprintf("aggregate_tables: read %s failed: %s",
                                basename(f), conditionMessage(e)))
    NULL
  }))
cells <- cells[!vapply(cells, is.null, logical(1L))]
log_msg(sprintf("Loaded %d cells.", length(cells)))


# ── tab_01: internal quality + convergence ──────────────────────────────────

build_internal_quality <- function(cells) {
  rows <- list(); k <- 0L
  for (cell in cells) {
    for (metric in METRICS) {
      k <- k + 1L
      fb <- fda_block(cell, metric)
      cv <- convergence_block(cell, metric)
      rows[[k]] <- data.frame(
        dataset           = cell$dataset,
        variant           = cell$variant,
        strategy          = cell$strategy,
        metric            = metric,
        cell_status       = cell$status %||% NA_character_,
        cell_status_reason = cell$status_reason %||% NA_character_,
        predecessor_name  = cell$predecessor_name %||% NA_character_,
        n_class1_sites    = cell$gmm_meta$n_class1_sites %||% NA_integer_,
        gmm_G             = cell$gmm_meta$G %||% NA_integer_,
        gmm_modelName     = cell$gmm_meta$modelName %||% NA_character_,
        n_bins            = cell$n_bins %||% NA_integer_,
        n_frames          = length(frames_block(cell, metric)),
        k_selected        = (if (!is.null(fb)) fb$n_clusters else NA_integer_),
        silhouette_mean   = (if (!is.null(fb)) fb$silhouette_mean else NA_real_),
        metric_status     = (if (!is.null(fb)) fb$status else NA_character_),
        convergence_first_snp    = (if (!is.null(cv)) cv$first_snp else NA_integer_),
        convergence_first_mut    = (if (!is.null(cv)) cv$first_mut else NA_integer_),
        convergence_first_class1 = (if (!is.null(cv)) cv$first_class1 else NA_integer_),
        walltime_s        = cell$walltime_s %||% NA_real_,
        stringsAsFactors  = FALSE
      )
    }
  }
  do.call(rbind, rows)
}


# ── tab_02: cluster membership ──────────────────────────────────────────────

build_cluster_membership <- function(cells) {
  rows <- list(); k <- 0L
  for (cell in cells) {
    if (!identical(cell$status, "ok")) next
    snp_set <- cell$snp_sites_truth
    mut_set <- cell$mutation_sites_truth

    for (metric in METRICS) {
      fb <- fda_block(cell, metric)
      if (is.null(fb) || !identical(fb$status, "ok")) next
      sites <- cell$gmm_meta$class1_sites
      mem   <- fb$cluster_assignments
      if (length(sites) != length(mem)) next
      for (j in seq_along(sites)) {
        k <- k + 1L
        rows[[k]] <- data.frame(
          dataset    = cell$dataset,
          variant    = cell$variant,
          strategy   = cell$strategy,
          metric     = metric,
          cluster_id = as.integer(mem[j]),
          site       = as.integer(sites[j]),
          is_snp     = sites[j] %in% snp_set,
          is_mutation_site = sites[j] %in% mut_set,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L)
    return(data.frame(dataset = character(0), variant = character(0),
                       strategy = character(0), metric = character(0),
                       cluster_id = integer(0), site = integer(0),
                       is_snp = logical(0), is_mutation_site = logical(0)))
  do.call(rbind, rows)
}


# ── tab_03: stability ───────────────────────────────────────────────────────

build_stability <- function(cells, jaccard_threshold) {
  rows <- list(); k <- 0L
  for (cell in cells) {
    if (!identical(cell$status, "ok")) next
    for (metric in METRICS) {
      stab <- stability_block(cell, metric)
      if (is.null(stab)) next
      pcj <- stab$per_cluster_jaccard
      if (length(pcj) == 0L) next
      for (cluster_id in seq_along(pcj)) {
        k <- k + 1L
        rows[[k]] <- data.frame(
          dataset      = cell$dataset,
          variant      = cell$variant,
          strategy     = cell$strategy,
          metric       = metric,
          cluster_id   = as.integer(cluster_id),
          mean_jaccard = unname(pcj[cluster_id]),
          is_stable    = (!is.na(pcj[cluster_id])) &&
                          pcj[cluster_id] >= jaccard_threshold,
          n_bootstrap  = stab$n_bootstrap,
          n_failed_bootstrap = stab$n_failed_bootstrap %||% NA_integer_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L)
    return(data.frame(dataset = character(0), variant = character(0),
                       strategy = character(0), metric = character(0),
                       cluster_id = integer(0), mean_jaccard = numeric(0),
                       is_stable = logical(0), n_bootstrap = integer(0),
                       n_failed_bootstrap = integer(0)))
  do.call(rbind, rows)
}


# ── tab_04: cross-strategy ARI ──────────────────────────────────────────────

build_cross_strategy_ari <- function(cells) {

  grp_keys <- vapply(cells, function(c) sprintf("%s__%s", c$dataset, c$variant),
                      character(1L))
  groups <- split(seq_along(cells), grp_keys)

  pairs <- list(
    c("sliding_2m",  "disjoint_2m"),
    c("sliding_2m",  "cumulative_1m"),
    c("disjoint_2m", "cumulative_1m")
  )

  rows <- list(); k <- 0L
  for (gk in names(groups)) {
    idxs <- groups[[gk]]
    by_strategy <- setNames(idxs,
                              vapply(idxs, function(i) cells[[i]]$strategy,
                                      character(1L)))
    parts <- strsplit(gk, "__", fixed = TRUE)[[1L]]
    ds_name  <- parts[1L]; var_name <- parts[2L]

    for (metric in METRICS) {
      for (pr in pairs) {
        a_name <- pr[1L]; b_name <- pr[2L]
        if (!(a_name %in% names(by_strategy)) ||
            !(b_name %in% names(by_strategy))) next

        cell_a <- cells[[by_strategy[[a_name]]]]
        cell_b <- cells[[by_strategy[[b_name]]]]
        if (!identical(cell_a$status, "ok") ||
            !identical(cell_b$status, "ok")) next

        fb_a <- fda_block(cell_a, metric)
        fb_b <- fda_block(cell_b, metric)
        if (is.null(fb_a) || is.null(fb_b)) next
        if (!identical(fb_a$status, "ok") || !identical(fb_b$status, "ok")) next

        sites_a <- cell_a$gmm_meta$class1_sites
        sites_b <- cell_b$gmm_meta$class1_sites
        common <- intersect(sites_a, sites_b)
        if (length(common) < 3L) next

        mem_a <- fb_a$cluster_assignments[match(common, sites_a)]
        mem_b <- fb_b$cluster_assignments[match(common, sites_b)]

        ari <- tryCatch(
          mclust::adjustedRandIndex(mem_a, mem_b),
          error = function(e) NA_real_
        )

        k <- k + 1L
        rows[[k]] <- data.frame(
          dataset        = ds_name,
          variant        = var_name,
          metric         = metric,
          strategy_a     = a_name,
          strategy_b     = b_name,
          n_sites_common = length(common),
          k_a            = fb_a$n_clusters,
          k_b            = fb_b$n_clusters,
          ARI            = ari,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L)
    return(data.frame(dataset = character(0), variant = character(0),
                       metric = character(0), strategy_a = character(0),
                       strategy_b = character(0), n_sites_common = integer(0),
                       k_a = integer(0), k_b = integer(0), ARI = numeric(0)))
  do.call(rbind, rows)
}


# ── tab_05: SNP enrichment ──────────────────────────────────────────────────

build_snp_enrichment <- function(cells) {
  rows <- list(); k <- 0L
  for (cell in cells) {
    for (metric in METRICS) {
      snp <- enrichment_block(cell, metric)
      k <- k + 1L
      rows[[k]] <- data.frame(
        dataset             = cell$dataset,
        variant             = cell$variant,
        strategy            = cell$strategy,
        metric              = metric,
        cell_status         = cell$status %||% NA_character_,
        enrichment_status   = (if (is.null(snp)) NA_character_ else
                                 snp$status %||% NA_character_),
        n_snp_in_class1     = (if (is.null(snp)) NA_integer_ else
                                 snp$n_snp_in_class1 %||% NA_integer_),
        n_class1_total      = (if (is.null(snp)) NA_integer_ else
                                 snp$n_class1_total %||% NA_integer_),
        n_snp_total_truth   = (if (is.null(snp)) NA_integer_ else
                                 snp$n_snp_total_truth %||% NA_integer_),
        fisher_p            = (if (is.null(snp)) NA_real_ else
                                 snp$fisher_p %||% NA_real_),
        fisher_p_bonferroni = (if (is.null(snp)) NA_real_ else
                                 snp$fisher_p_bonferroni %||% NA_real_),
        max_overlap_cluster = (if (is.null(snp)) NA_integer_ else
                                 snp$max_overlap_cluster %||% NA_integer_),
        max_overlap_OR      = (if (is.null(snp)) NA_real_ else
                                 snp$max_overlap_OR %||% NA_real_),
        stringsAsFactors    = FALSE
      )
    }
  }
  do.call(rbind, rows)
}


# ── tab_06: per-frame summary ───────────────────────────────────────────────

build_per_frame <- function(cells) {
  rows <- list(); k <- 0L
  for (cell in cells) {
    if (!identical(cell$status, "ok")) next
    for (metric in METRICS) {
      frames <- frames_block(cell, metric)
      if (length(frames) == 0L) next
      for (i in seq_along(frames)) {
        f <- frames[[i]]
        k <- k + 1L
        rows[[k]] <- data.frame(
          dataset            = cell$dataset,
          variant            = cell$variant,
          strategy           = cell$strategy,
          metric             = metric,
          frame_idx          = i,
          centre_bin         = f$centre_bin,
          fitting_start_bin  = min(f$frame_bins),
          fitting_end_bin    = max(f$frame_bins),
          gmm_status         = f$gmm_status,
          gmm_G              = f$gmm_G %||% NA_integer_,
          gmm_n_sequences    = f$gmm_n_sequences %||% NA_integer_,
          fit_status         = f$fit_status,
          n_class1_sites     = f$n_class1_sites,
          n_clusters         = f$n_clusters %||% NA_integer_,
          silhouette_mean    = f$silhouette_mean %||% NA_real_,
          structurally_null  = f$structurally_null %||% NA,
          n_snp_in_class1    = f$n_snp_in_class1,
          n_mut_in_class1    = f$n_mut_in_class1,
          converged_snp      = isTRUE(f$converged_snp),
          converged_mut      = isTRUE(f$converged_mut),
          converged_class1   = isTRUE(f$converged_class1),
          qualifying_cluster_snp    = f$qualifying_cluster_snp,
          qualifying_cluster_mut    = f$qualifying_cluster_mut,
          qualifying_cluster_class1 = f$qualifying_cluster_class1,
          anchor_idx         = f$anchor_idx,
          n_oversized_clusters = (if (is.null(f$oversized)) NA_integer_
                                    else sum(f$oversized)),
          stringsAsFactors   = FALSE
        )
      }
    }
  }
  if (length(rows) == 0L)
    return(data.frame())
  do.call(rbind, rows)
}


# ── Build & write all tables ────────────────────────────────────────────────

tabs <- list(
  tab_01_internal_quality   = build_internal_quality(cells),
  tab_02_cluster_membership = build_cluster_membership(cells),
  tab_03_stability          = build_stability(cells,
                                                config$STABLE_JACCARD_THRESHOLD),
  tab_04_cross_strategy_ari = build_cross_strategy_ari(cells),
  tab_05_snp_enrichment     = build_snp_enrichment(cells),
  tab_06_per_frame          = build_per_frame(cells)
)

for (nm in names(tabs)) {
  out <- file.path(tables_dir(config), paste0(nm, ".csv"))
  write_csv_atomic(tabs[[nm]], out)
  log_msg(sprintf("Wrote %s  (%d rows)", basename(out), nrow(tabs[[nm]])))
}

log_msg("aggregate_tables complete.")
log_msg(sprintf("  → %s", tables_dir(config)))
