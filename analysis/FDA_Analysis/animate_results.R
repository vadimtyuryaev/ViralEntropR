################################################################################
## animate_results.R — FDA_Analysis (v2)
##
## Reads outputs/cells/cell_*.rds and produces one GIF per (cell, metric) for
## cells with status == "ok". Animations are built from the per-frame records
## already computed by run_one_cell.R (no re-fitting at render time).
##
## Output layout:
##   outputs/gifs/<dataset>/<variant>/<strategy>/anim_<metric>.gif
##
## Invocation: Rscript animate_results.R
################################################################################

suppressPackageStartupMessages({
  source("setup.R")
  source("helpers_fda.R")
  source("helpers_frames.R")
  source("helpers_panels.R")
  library(ggplot2)
  library(magick)
})

require_runtime_packages(animation = TRUE)

config <- build_config()
setwd(config$STUDY_DIR)
dir.create(gifs_dir(config), recursive = TRUE, showWarnings = FALSE)


# ── Per-cell GIF rendering ──────────────────────────────────────────────────

render_cell_gifs <- function(cell, config) {

  if (!identical(cell$status, "ok")) return(c(0L, 0L, 1L))

  outd <- file.path(gifs_dir(config), cell$dataset, cell$variant, cell$strategy)
  dir.create(outd, recursive = TRUE, showWarnings = FALSE)

  done <- 0L; err <- 0L

  for (metric in c("entropy", "hellinger_T1", "hellinger_Tpred")) {
    frames <- switch(
      metric,
      entropy         = cell$frames_entropy,
      hellinger_T1    = cell$frames_hellinger_T1,
      hellinger_Tpred = cell$frames_hellinger_Tpred
    )
    if (length(frames) == 0L) next

    # Skip metrics where the cell-level fit was a no-go (no usable matrix).
    full_mat <- switch(
      metric,
      entropy         = cell$entropy_full_matrix,
      hellinger_T1    = cell$hellinger_T1_full_matrix,
      hellinger_Tpred = cell$hellinger_Tpred_full_matrix
    )
    if (is.null(full_mat)) next

    out_gif <- file.path(outd, sprintf("anim_%s.gif", metric))
    res <- tryCatch(
      assemble_gif(cell, frames, metric, out_gif, config),
      error = function(e) {
        log_error(config,
                  sprintf("[animate_results] %s/%s/%s metric=%s: %s",
                          cell$dataset, cell$variant, cell$strategy,
                          metric, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res)) err <- err + 1L else done <- done + 1L

    # Defensive device cleanup between GIFs
    while (length(dev.list()) > 0L) try(dev.off(), silent = TRUE)
  }

  c(done, err, 0L)
}


# ── Main loop ────────────────────────────────────────────────────────────────

files <- list.files(cells_dir(config), pattern = "^cell_.+\\.rds$",
                     full.names = TRUE)
log_msg(sprintf("animate_results: %d cell RDS files.", length(files)))

done <- 0L; err <- 0L; skip <- 0L
for (f in files) {
  cell <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(cell)) { err <- err + 1L; next }
  res <- render_cell_gifs(cell, config)
  done <- done + res[1L]
  err  <- err  + res[2L]
  skip <- skip + res[3L]
}

log_msg(sprintf("animate_results complete: %d GIFs written, %d cells skipped, %d errored.",
                 done, skip, err))
log_msg(sprintf("  → %s", gifs_dir(config)))
