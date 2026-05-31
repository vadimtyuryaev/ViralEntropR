################################################################################
## render_existing.R — FDA_Analysis (v2.1)
##
## Fast re-render of GIFs from an EXISTING cell RDS. Use this after editing
## helpers_panels.R to preview the new design without re-running the 20-40
## minute per-cell fit. Total time per cell: ~60-90 seconds (3 metrics x
## ~25 s of ggplot rendering + GIF assembly).
##
## Usage:
##   cd ~/Clean_Code_Running/Dissertation/FDA_Analysis
##   Rscript render_existing.R 2>&1 | tee outputs/render_log.txt
##
## Env overrides (same as smoke_test.R):
##   SMOKE_DATASET / SMOKE_VARIANT / SMOKE_STRATEGY
##
## The script reads outputs/cells/cell_<ds>__<vrn>__<strg>.rds, deletes the
## previous GIF files for that cell (so file timestamps clearly reflect the
## re-render), and writes fresh ones to outputs/gifs/<ds>/<vrn>/<strg>/.
################################################################################

t_total <- Sys.time()
step <- function(n, msg) {
  cat(sprintf("[%s] STEP %s: %s\n",
              format(Sys.time(), "%H:%M:%S"), n, msg))
  flush.console()
}

ds   <- Sys.getenv("SMOKE_DATASET",  "NCBI_US")
vrn  <- Sys.getenv("SMOKE_VARIANT",  "Gamma")
strg <- Sys.getenv("SMOKE_STRATEGY", "disjoint_2m")

# SAVE_FRAMES=TRUE (default) writes each frame as a standalone PNG to
# outputs/plots/<ds>/<vrn>/<strg>/frame_<metric>_<centre:03d>.png so the
# user can inspect frames individually. Set SAVE_FRAMES=FALSE to skip.
SAVE_FRAMES_ENV <- Sys.getenv("SAVE_FRAMES", "TRUE")
SAVE_FRAMES     <- !identical(toupper(SAVE_FRAMES_ENV), "FALSE")

step("0/5", sprintf("re-rendering cell = %s x %s x %s  (SAVE_FRAMES=%s)",
                     ds, vrn, strg, SAVE_FRAMES))

# ── 1 ───────────────────────────────────────────────────────────────────────
step("1/5", "Sourcing helpers (no fitting code needed)")
t1 <- Sys.time()
source("setup.R")
source("helpers_fda.R")
source("helpers_frames.R")
source("helpers_panels.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(magick)
  library(patchwork)
})
cat(sprintf("        ...done in %.1fs\n",
            as.numeric(difftime(Sys.time(), t1, units = "secs"))))
flush.console()

# ── 2 ───────────────────────────────────────────────────────────────────────
step("2/5", "build_config() + locate cell RDS")
config <- build_config()
cell_rds <- cell_path(config, ds, vrn, strg)
if (!file.exists(cell_rds)) {
  stop("Cell RDS not found: ", cell_rds,
       "\n  Run smoke_test.R or fda_analysis.R first to produce it.",
       call. = FALSE)
}
cat(sprintf("        cell RDS: %s (%.1f MB)\n",
            cell_rds, file.info(cell_rds)$size / 1e6))

# ── 3 ───────────────────────────────────────────────────────────────────────
step("3/5", "readRDS cell payload")
t1 <- Sys.time()
res <- readRDS(cell_rds)
cat(sprintf("        ...done in %.1fs\n",
            as.numeric(difftime(Sys.time(), t1, units = "secs"))))
cat(sprintf("        schema_version=%s, status=%s, n_bins=%s\n",
            res$schema_version, res$status, res$n_bins))
cat(sprintf("        |class1|=%s, predecessor=%s, detection_bin=%s\n",
            res$gmm_meta$n_class1_sites,
            res$predecessor_name %||% "NA",
            res$detection_bin %||% "NA"))

if (!identical(res$status, "ok"))
  stop("cell status = ", res$status, "; cannot render.", call. = FALSE)

# Build the variants catalogue ONCE for the whole render — used by
# write_cluster_membership_csv() to attach per-site variant-membership
# lists and per-cluster majority-voting columns. The catalogue makes the
# CSV variant-agnostic (identical content across variants of the same
# base cell).
variants_catalogue <- tryCatch(
  build_variants_catalogue(sarscov2_variants),
  error = function(e) {
    cat(sprintf("        WARN: variants catalogue build failed: %s\n",
                conditionMessage(e)))
    NULL
  }
)
if (!is.null(variants_catalogue))
  cat(sprintf("        variants catalogue: %d variants\n",
              length(variants_catalogue)))

# ── 4 ───────────────────────────────────────────────────────────────────────
step("4/5", "Rendering 3 GIFs (entropy, hellinger_T1, hellinger_Tpred)")

gif_dir <- file.path(gifs_dir(config), ds, vrn, strg)
dir.create(gif_dir, recursive = TRUE, showWarnings = FALSE)

# Remove old GIFs so file timestamps unambiguously reflect this re-render
old_gifs <- list.files(gif_dir, pattern = "^anim_.*\\.gif$", full.names = TRUE)
if (length(old_gifs) > 0L) {
  cat(sprintf("        removing %d old GIF(s):\n", length(old_gifs)))
  for (g in old_gifs) {
    cat(sprintf("          %s\n", basename(g)))
    file.remove(g)
  }
}

gif_paths <- character(0L)
frames_dir <- if (SAVE_FRAMES) {
  file.path(plots_dir(config), ds, vrn, strg)
} else {
  NULL
}
if (!is.null(frames_dir)) {
  dir.create(frames_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("        frame PNGs -> %s\n", frames_dir))
}

for (metric in c("entropy", "hellinger_T1", "hellinger_Tpred")) {
  frames <- switch(metric,
                   entropy         = res$frames_entropy,
                   hellinger_T1    = res$frames_hellinger_T1,
                   hellinger_Tpred = res$frames_hellinger_Tpred)
  out_gif <- file.path(gif_dir, sprintf("anim_%s.gif", metric))
  cat(sprintf("        [%s] %d frames -> %s\n",
              metric, length(frames), basename(out_gif)))
  flush.console()

  t1 <- Sys.time()
  err <- tryCatch({
    assemble_gif(cell = res, frame_records = frames,
                  metric = metric, output_path = out_gif,
                  config = config,
                  frames_dir = frames_dir)
    NULL
  }, error = function(e) conditionMessage(e))

  while (length(dev.list()) > 0L) try(dev.off(), silent = TRUE)

  dt <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  if (is.null(err)) {
    cat(sprintf("          rendered in %.1fs\n", dt))
    gif_paths <- c(gif_paths, out_gif)
  } else {
    cat(sprintf("          FAILED in %.1fs: %s\n", dt, err))
  }

  # Per-frame cluster membership CSV (one row per class-1 site per frame)
  if (!is.null(frames_dir)) {
    csv_path <- file.path(frames_dir,
                           sprintf("cluster_membership_%s.csv", metric))
    csv_err <- tryCatch({
      write_cluster_membership_csv(res, frames, metric, csv_path,
                                    variants_catalogue = variants_catalogue)
      NULL
    }, error = function(e) conditionMessage(e))
    if (is.null(csv_err))
      cat(sprintf("          cluster_membership CSV -> %s\n",
                  basename(csv_path)))
    else
      cat(sprintf("          CSV write FAILED: %s\n", csv_err))
  }
  flush.console()
}

# ── 5 ───────────────────────────────────────────────────────────────────────
step("5/5", "Verifying GIF outputs + total elapsed")
for (gp in gif_paths) {
  sz <- file.info(gp)$size
  info <- tryCatch(magick::image_info(magick::image_read(gp)),
                   error = function(e) NULL)
  if (is.null(info))
    cat(sprintf("        %s  (%d bytes, magick read FAILED)\n",
                basename(gp), sz))
  else
    cat(sprintf("        %s  (%d bytes, %d frames, %dx%d px)\n",
                basename(gp), sz, nrow(info),
                info$width[1L], info$height[1L]))
}
cat(sprintf("\nRENDER COMPLETE — wall time %.1fs.\n",
            as.numeric(difftime(Sys.time(), t_total, units = "secs"))))
cat(sprintf("GIFs at: %s\n", gif_dir))
if (!is.null(frames_dir)) {
  n_pngs <- length(list.files(frames_dir, pattern = "^frame_.*\\.png$"))
  cat(sprintf("Frame PNGs at: %s  (%d files)\n", frames_dir, n_pngs))
}
