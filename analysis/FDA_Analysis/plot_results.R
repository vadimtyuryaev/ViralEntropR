################################################################################
## plot_results.R — FDA_Analysis (v2)
##
## Reads outputs/cells/cell_*.rds and produces:
##
##   1. One PNG per (cell × metric × frame) — `frame_<metric>_<centre:03d>.png`
##      in outputs/plots/<dataset>/<variant>/<strategy>/. These are the same
##      images that compose the GIFs from animate_results.R; storing them as
##      individual PNGs lets the user cite specific frames in a paper.
##
##   2. One "summary" PNG per (cell × metric) — `summary_<metric>.png` in the
##      same directory. The summary plot uses the FULL-WINDOW cell-level fit
##      (not a per-frame fit), rendered in the same dual-panel format. It's
##      the "headline result" image for each cell.
##
##   3. A manifest CSV per (cell × strategy) directory listing all frames and
##      their fit_status — `frames.manifest.csv`. Lets downstream readers
##      find frames without parsing filenames.
##
## Cells with overall status != "ok" are silently skipped at the cell level;
## individual frames with fit_status != "ok" are still rendered (showing the
## null state) so the animation timeline is complete.
##
## Invocation: Rscript plot_results.R
################################################################################

suppressPackageStartupMessages({
  source("setup.R")
  source("helpers_fda.R")
  source("helpers_frames.R")
  source("helpers_panels.R")
  library(ggplot2)
})

require_runtime_packages(animation = FALSE)

config <- build_config()
setwd(config$STUDY_DIR)
dir.create(plots_dir(config), recursive = TRUE, showWarnings = FALSE)


# ── 1. Per-cell directory layout ─────────────────────────────────────────────

cell_plots_dir <- function(cell, config) {
  d <- file.path(plots_dir(config), cell$dataset, cell$variant, cell$strategy)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}


# ── 2. Build a "summary" frame record from the FULL-WINDOW cell-level fit ───
#
# Wraps the cell-level fda_* and stability_* into a frame-record-shaped list
# so the same panel-building functions can render it.

cell_to_summary_frame <- function(cell, metric) {

  fda  <- switch(
    metric,
    entropy         = cell$fda_entropy,
    hellinger_T1    = cell$fda_hellinger_T1,
    hellinger_Tpred = cell$fda_hellinger_Tpred
  )
  if (is.null(fda) || !identical(fda$status, "ok")) return(NULL)

  M_t <- switch(
    metric,
    entropy         = cell$entropy_class1_matrix_transformed,
    hellinger_T1    = cell$hellinger_T1_class1_matrix_transformed,
    hellinger_Tpred = cell$hellinger_Tpred_class1_matrix_transformed
  )
  M_raw <- switch(
    metric,
    entropy         = cell$entropy_class1_matrix,
    hellinger_T1    = cell$hellinger_T1_class1_matrix,
    hellinger_Tpred = cell$hellinger_Tpred_class1_matrix
  )
  if (is.null(M_t) || nrow(M_t) == 0L) return(NULL)

  grid_vec <- if (metric == "entropy") seq_len(ncol(M_t))
              else as.integer(sub("^T", "", colnames(M_t)))

  mem <- fda$cluster_assignments
  sil <- tryCatch(per_site_silhouette(M_t, mem),
                   error = function(e) setNames(rep(NA_real_, nrow(M_t)),
                                                  rownames(M_t)))
  cluster_sizes <- tabulate(mem, nbins = fda$n_clusters)
  oversized     <- cluster_sizes > config$MAX_CLUSTER_SIZE
  structurally_null <- is.na(fda$silhouette_mean) ||
                       fda$silhouette_mean < config$MIN_MEAN_SILHOUETTE

  list(
    centre_bin        = NA_integer_,
    frame_bins        = seq_len(cell$n_bins),
    metric            = metric,
    gmm_status        = "ok",
    gmm_G             = cell$gmm_meta$G,
    gmm_modelName     = cell$gmm_meta$modelName,
    gmm_n_sequences   = cell$n_seqs_window,
    class1_sites      = cell$gmm_meta$class1_sites,
    n_class1_sites    = cell$gmm_meta$n_class1_sites,
    curve_matrix      = M_raw,
    curve_matrix_transformed = M_t,
    grid_vec          = grid_vec,
    n_clusters        = fda$n_clusters,
    silhouette_mean   = fda$silhouette_mean,
    silhouette_per_k  = fda$silhouette_per_k,
    cluster_assignments = mem,
    per_site_sil      = sil,
    cluster_sizes     = cluster_sizes,
    oversized         = oversized,
    structurally_null = structurally_null,
    converged_snp     = NA, converged_mut = NA, converged_class1 = NA,
    qualifying_cluster_snp    = NA_integer_,
    qualifying_cluster_mut    = NA_integer_,
    qualifying_cluster_class1 = NA_integer_,
    n_snp_in_class1   = sum(cell$gmm_meta$class1_sites %in% cell$snp_sites_truth),
    n_mut_in_class1   = sum(cell$gmm_meta$class1_sites %in% cell$mutation_sites_truth),
    anchor_idx        = NA_integer_,
    fit_status        = "ok",
    reason            = NA_character_
  )
}


# ── 3. Manifest builder ──────────────────────────────────────────────────────

build_manifest <- function(cell) {
  rows <- list()
  for (metric in c("entropy", "hellinger_T1", "hellinger_Tpred")) {
    frames <- switch(
      metric,
      entropy         = cell$frames_entropy,
      hellinger_T1    = cell$frames_hellinger_T1,
      hellinger_Tpred = cell$frames_hellinger_Tpred
    )
    if (length(frames) == 0L) next
    for (i in seq_along(frames)) {
      f <- frames[[i]]
      rows[[length(rows) + 1L]] <- data.frame(
        dataset          = cell$dataset,
        variant          = cell$variant,
        strategy         = cell$strategy,
        metric           = metric,
        frame_idx        = i,
        centre_bin       = f$centre_bin,
        fitting_start    = min(f$frame_bins),
        fitting_end      = max(f$frame_bins),
        gmm_status       = f$gmm_status,
        fit_status       = f$fit_status,
        n_class1_sites   = f$n_class1_sites,
        n_clusters       = f$n_clusters %||% NA_integer_,
        silhouette_mean  = f$silhouette_mean %||% NA_real_,
        structurally_null = f$structurally_null %||% NA,
        converged_snp    = f$converged_snp,
        converged_mut    = f$converged_mut,
        converged_class1 = f$converged_class1,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}


# ── 4. Cell renderer ─────────────────────────────────────────────────────────

render_cell <- function(cell, config) {

  if (!identical(cell$status, "ok")) return(c(0L, 0L, 1L))
  outd <- cell_plots_dir(cell, config)

  png_done <- 0L; errored <- 0L

  for (metric in c("entropy", "hellinger_T1", "hellinger_Tpred")) {

    frames <- switch(
      metric,
      entropy         = cell$frames_entropy,
      hellinger_T1    = cell$frames_hellinger_T1,
      hellinger_Tpred = cell$frames_hellinger_Tpred
    )

    # ── Per-frame PNGs ──
    if (length(frames) > 0L) {
      for (i in seq_along(frames)) {
        f <- frames[[i]]
        out_png <- file.path(outd,
                              sprintf("frame_%s_%03d.png", metric,
                                      as.integer(f$centre_bin)))
        ok <- tryCatch({
          p <- build_frame_plot(cell, f, metric, i, length(frames), config)
          ggplot2::ggsave(out_png, p,
                           width = config$PNG_WIDTH_CM,
                           height = config$PNG_HEIGHT_CM,
                           units = "cm", dpi = 300, bg = "white")
          TRUE
        }, error = function(e) {
          log_error(config,
                    sprintf("[plot_results] %s/%s/%s metric=%s frame=%d: %s",
                            cell$dataset, cell$variant, cell$strategy,
                            metric, i, conditionMessage(e)))
          FALSE
        })
        if (ok) png_done <- png_done + 1L else errored <- errored + 1L
      }
    }

    # ── Summary PNG (full-window cell-level fit) ──
    summ <- tryCatch(cell_to_summary_frame(cell, metric),
                      error = function(e) NULL)
    if (!is.null(summ)) {
      out_png <- file.path(outd, sprintf("summary_%s.png", metric))
      ok <- tryCatch({
        p <- build_frame_plot(cell, summ, metric,
                              frame_idx = 0L, n_frames = 0L, config)
        ggplot2::ggsave(out_png, p,
                         width = config$PNG_WIDTH_CM,
                         height = config$PNG_HEIGHT_CM,
                         units = "cm", dpi = 300, bg = "white")
        TRUE
      }, error = function(e) {
        log_error(config,
                  sprintf("[plot_results] %s/%s/%s summary=%s: %s",
                          cell$dataset, cell$variant, cell$strategy,
                          metric, conditionMessage(e)))
        FALSE
      })
      if (ok) png_done <- png_done + 1L else errored <- errored + 1L
    }
  }

  # ── Manifest CSV ──
  man <- build_manifest(cell)
  if (!is.null(man)) {
    out_csv <- file.path(outd, "frames.manifest.csv")
    tryCatch({
      utils::write.csv(man, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
    }, error = function(e) {
      log_error(config,
                sprintf("[plot_results] manifest %s/%s/%s: %s",
                        cell$dataset, cell$variant, cell$strategy,
                        conditionMessage(e)))
    })
  }

  c(png_done, errored, 0L)
}


# ── 5. Main loop ─────────────────────────────────────────────────────────────

files <- list.files(cells_dir(config), pattern = "^cell_.+\\.rds$",
                     full.names = TRUE)
log_msg(sprintf("plot_results: %d cell RDS files.", length(files)))

png_done <- 0L; skipped <- 0L; errored <- 0L

for (f in files) {
  cell <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(cell)) { errored <- errored + 1L; next }

  res <- render_cell(cell, config)
  png_done <- png_done + res[1L]
  errored  <- errored  + res[2L]
  skipped  <- skipped  + res[3L]
}

log_msg(sprintf("plot_results complete: %d PNGs written, %d cells skipped, %d errored.",
                 png_done, skipped, errored))
log_msg(sprintf("  → %s", plots_dir(config)))
