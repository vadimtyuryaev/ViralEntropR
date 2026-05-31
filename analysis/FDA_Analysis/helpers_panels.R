################################################################################
## helpers_panels.R — FDA_Analysis (v2.1, "clean trajectory" design)
##
## Dual-panel rendering primitives, shared by plot_results.R and
## animate_results.R. v2.1 ports the original `custom_clust_gif_plot`
## aesthetic from FDA_Clustering.Rmd:
##
##   ┌─────────────────────────────────────────────────────────────┐
##   │ TOP PANEL — class-1 trajectories coloured by cluster        │
##   │   • 1273 background curves in light grey                    │
##   │   • Class-1 site curves in cluster colour (Set3 palette)    │
##   │   • Solid inside the fitting window, dashed outside         │
##   │   • Triangles ▲ start, ▼ end at top edge of the panel       │
##   │   • Filled red circle at the detection bin on the x-axis    │
##   │   • Legend at top: "Clusters" with member sites listed      │
##   ├─────────────────────────────────────────────────────────────┤
##   │ BOTTOM PANEL — same trajectory plot, defining-SNP highlight │
##   │   • 1273 background curves in light grey                    │
##   │   • Variant-defining SNP curves in RED                      │
##   │   • Same solid/dashed convention                            │
##   │   • Same window-marker and detection-bin annotations        │
##   │   • Subtitle: "<Variant> SNPs: <site1, site2, ...>"          │
##   └─────────────────────────────────────────────────────────────┘
##
## Y-axis values are displayed RAW (not log1p). The log1p transform is
## still used inside the fitting (clustering uses transformed curves),
## but readers see raw entropy / Hellinger values, matching the original
## design's "Hellinger Distance 0–1.2" axis.
##
## Public entry points (signatures preserved from v2):
##   build_top_panel(cell, frame_record, metric, config)
##   build_bottom_panel(cell, frame_record, metric, config)
##   build_dual_panel(top, bottom, title_str, subtitle_str, config)
##   build_frame_plot(cell, frame_record, metric, frame_idx, n_frames, config)
##   assemble_gif(cell, frame_records, metric, output_path, config)
##
## Pre-requisite: setup.R + helpers_fda.R + helpers_frames.R sourced.
################################################################################


# ── 1. Palettes and constants ───────────────────────────────────────────────
#
# Set3 (RColorBrewer) — 12 pastel hues. Cluster k receives palette[k]
# (cycled if k > 12). Re-cycling beyond 12 is a non-issue since
# N_CLUSTERS_GRID caps k at 6.

SET3_PALETTE <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
                  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                  "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

GREY_BG       <- "#D9D9D9"    # background curves (all 1273 sites)
GREY_OVERSIZE <- "#BBBBBB"    # oversized cluster (visualisation gate)
RED_HIGHLIGHT <- "#E41A1C"    # SNP highlight + marginal sites (sil<0)
DETECTION_COL <- "#A50026"    # filled red circle at the detection bin
WINDOW_COL    <- "grey20"     # fitting-window triangle markers


# ── 2. Metric -> matrix / axis-label dispatch ────────────────────────────────

metric_full_matrix <- function(cell, metric) {
  switch(metric,
         entropy         = cell$entropy_full_matrix,
         hellinger_T1    = cell$hellinger_T1_full_matrix,
         hellinger_Tpred = cell$hellinger_Tpred_full_matrix,
         stop("Unknown metric: ", metric, call. = FALSE))
}

metric_y_label <- function(metric) {
  switch(metric,
         entropy         = "Entropy (bits)",
         hellinger_T1    = "Hellinger Distance",
         hellinger_Tpred = "Hellinger Distance",
         "value")
}

metric_anchor_label <- function(cell, metric) {
  switch(metric,
         entropy         = NULL,
         hellinger_T1    = "T1 (window start)",
         hellinger_Tpred = sprintf("%s (predecessor)",
                                   cell$predecessor_name %||% "T_pred"))
}


# ── 3. Bin-index mapping for the x-axis ─────────────────────────────────────
#
# Entropy matrices use column index = bin index. Hellinger matrices use
# columns "T2", "T3", ... (anchor dropped). Returns integer vector of
# original bin indices that correspond to each matrix column.

matrix_to_bin_indices <- function(M, metric) {
  if (metric == "entropy") return(seq_len(ncol(M)))
  as.integer(sub("^T", "", colnames(M)))
}


# ── 4. Long-form trajectory builder ─────────────────────────────────────────
#
# Returns a data frame with columns site (chr), x (int = bin index),
# y (num), in_window (logical). `bin_indices` maps columns of M to bins.

make_traj_df <- function(M, bin_indices, frame_bins) {
  if (is.null(M) || nrow(M) == 0L || ncol(M) == 0L)
    return(data.frame(site = character(0L), x = integer(0L),
                      y = numeric(0L), in_window = logical(0L),
                      stringsAsFactors = FALSE))

  n_sites <- nrow(M); n_cols <- ncol(M)
  data.frame(
    site      = rep(rownames(M),  times = n_cols),
    x         = rep(bin_indices,  each  = n_sites),
    y         = as.numeric(M),
    in_window = rep(bin_indices %in% frame_bins, each = n_sites),
    stringsAsFactors = FALSE
  )
}


# ── 5. Per-site colour resolution (top panel only) ──────────────────────────
#
# Returns a named character vector aligned with frame_record$class1_sites.
# Three-state rule:
#   - well-clustered, non-oversized: Set3 palette colour for the cluster
#   - per-site silhouette < 0:        RED_HIGHLIGHT
#   - assigned to oversized cluster:  GREY_OVERSIZE

resolve_class1_colours <- function(frame_record) {
  n <- length(frame_record$class1_sites)
  if (n == 0L) return(setNames(character(0L), character(0L)))

  cls  <- frame_record$cluster_assignments
  sil  <- frame_record$per_site_sil
  over <- frame_record$oversized

  # Top-panel rule: only well-clustered sites get a Set3 colour. Marginal
  # (sil < 0), unassigned, and oversized-cluster members all fall back to
  # GREY_OVERSIZE. Red is reserved exclusively for the BOTTOM panel as a
  # "this SNP did not cluster" signal — keeping the two panels semantically
  # disjoint.
  cols <- character(n)
  for (i in seq_len(n)) {
    if (is.na(cls[i])) { cols[i] <- GREY_OVERSIZE; next }
    if (isTRUE(over[cls[i]])) { cols[i] <- GREY_OVERSIZE; next }
    if (is.na(sil[i]) || sil[i] < 0) { cols[i] <- GREY_OVERSIZE; next }
    cols[i] <- SET3_PALETTE[((cls[i] - 1L) %% length(SET3_PALETTE)) + 1L]
  }
  setNames(cols, as.character(frame_record$class1_sites))
}


# ── 6. Cluster-label legend builder ─────────────────────────────────────────
#
# Returns a list with character vectors: labels (one per cluster that
# contains AT LEAST ONE variant-defining site) and colours (matched 1:1).
# Oversized clusters and clusters with no defining members are excluded.
# Label format: "c<k>: site1, site2, ..." for compact display in the
# subtitle band of the top panel.

build_cluster_legend <- function(frame_record, defining_sites = integer(0L)) {
  if (frame_record$n_class1_sites == 0L ||
      !identical(frame_record$fit_status, "ok"))
    return(list(labels = character(0L), colours = character(0L)))

  cls    <- frame_record$cluster_assignments
  sites  <- frame_record$class1_sites
  over   <- frame_record$oversized
  sil    <- frame_record$per_site_sil
  K      <- max(cls, na.rm = TRUE)

  labels  <- character(0L)
  colours <- character(0L)
  for (k in seq_len(K)) {
    if (isTRUE(over[k])) next         # oversized: rendered grey, no legend entry
    in_k    <- which(cls == k & !is.na(sil) & sil >= 0)
    sites_k <- sites[in_k]
    if (length(sites_k) == 0L) next   # all members marginal

    # Filter to variant-defining sites for a compact, on-topic legend.
    if (length(defining_sites) > 0L) {
      defining_in_k <- intersect(sites_k, defining_sites)
      if (length(defining_in_k) == 0L) next  # cluster has no defining members
      sites_for_legend <- defining_in_k
    } else {
      sites_for_legend <- sites_k       # fall-back: show all class-1 members
    }

    if (length(sites_for_legend) <= 8L)
      lab <- sprintf("c%d: %s", k, paste(sites_for_legend, collapse = ", "))
    else
      lab <- sprintf("c%d: %s ... (%d defining)",
                      k,
                      paste(sites_for_legend[1:6], collapse = ", "),
                      length(sites_for_legend))

    labels[length(labels) + 1L]  <- lab
    colours[length(colours) + 1L] <- SET3_PALETTE[((k - 1L) %% length(SET3_PALETTE)) + 1L]
  }
  list(labels = labels, colours = colours)
}


# ── 6b. SNP colour resolution (bottom panel) ────────────────────────────────
#
# Maps each variant-defining SNP site to a colour based on its cluster
# assignment in the current frame. Linking principle: a SNP in the bottom
# panel inherits the colour of the cluster (top panel) it ended up in.
#
# Rules (per SNP):
#   - Not in class-1 (the GMM dropped it):       RED_HIGHLIGHT
#   - In class-1 but per-site silhouette < 0:    RED_HIGHLIGHT  (marginal)
#   - In class-1 in an oversized cluster:        GREY_OVERSIZE  (visual gate)
#   - Otherwise:                                 SET3 cluster colour
#
# If the frame's fit_status is not "ok" (e.g. GMM collapsed to G=1), every
# SNP returns RED_HIGHLIGHT.
#
# Returns a named character vector aligned with as.character(snp_sites).

resolve_snp_colours <- function(snp_sites, frame_record) {
  n <- length(snp_sites)
  if (n == 0L) return(setNames(character(0L), character(0L)))

  cols <- setNames(rep(RED_HIGHLIGHT, n), as.character(snp_sites))
  if (!identical(frame_record$fit_status, "ok")) return(cols)
  if (frame_record$n_class1_sites == 0L)         return(cols)

  cls_sites <- frame_record$class1_sites
  cls_mem   <- frame_record$cluster_assignments
  sil       <- frame_record$per_site_sil
  over      <- frame_record$oversized

  for (i in seq_along(snp_sites)) {
    s <- snp_sites[i]
    j <- match(s, cls_sites)
    if (is.na(j)) next                       # not in class-1: stay red
    k <- cls_mem[j]
    if (is.na(k)) next
    if (isTRUE(over[k])) {
      cols[as.character(s)] <- GREY_OVERSIZE # oversized cluster -> grey
      next
    }
    if (is.na(sil[j]) || sil[j] < 0) next    # marginal site: stay red
    cols[as.character(s)] <-
      SET3_PALETTE[((k - 1L) %% length(SET3_PALETTE)) + 1L]
  }
  cols
}


# ── 7. Window-marker + detection-circle annotation layers ───────────────────
#
# Adds the ▲ start, ▼ end triangles at the top of the panel and the filled
# red circle at the bottom (y = 0) of the detection bin.

add_window_and_detection <- function(p, cell, frame_record, y_max) {

  fb <- frame_record$frame_bins
  if (length(fb) >= 1L) {
    pad <- y_max * 0.06
    tri_df <- data.frame(
      x     = c(min(fb), max(fb)),
      y     = y_max + pad,
      label = c("\u25B2", "\u25BC")   # ▲ filled-up, ▼ filled-down
    )
    p <- p + ggplot2::geom_text(
      data = tri_df,
      ggplot2::aes(x = x, y = y, label = label),
      colour = WINDOW_COL, size = 4.5, inherit.aes = FALSE
    )
  }

  if (!is.null(cell$detection_bin) && !is.na(cell$detection_bin) &&
      cell$detection_bin >= 1L && cell$detection_bin <= cell$n_bins) {
    p <- p + ggplot2::annotate(
      "point",
      x      = cell$detection_bin,
      y      = 0,
      shape  = 16, size = 4,
      colour = DETECTION_COL
    )
  }
  p
}


# ── 8. Top panel — class-1 trajectories coloured by cluster ─────────────────

build_top_panel <- function(cell, frame_record, metric, config) {

  bg_mat <- metric_full_matrix(cell, metric)
  if (is.null(bg_mat))
    stop("No full-cell matrix for metric '", metric, "'.", call. = FALSE)

  bin_idx_vec <- matrix_to_bin_indices(bg_mat, metric)
  fb          <- frame_record$frame_bins

  bg_df <- make_traj_df(bg_mat, bin_idx_vec, fb)

  fg_df    <- NULL
  leg_list <- list(labels = character(0L), colours = character(0L))
  if (frame_record$n_class1_sites > 0L &&
      identical(frame_record$fit_status, "ok")) {

    cls_sites <- as.character(frame_record$class1_sites)
    cls_mat   <- bg_mat[cls_sites, , drop = FALSE]
    site_cols <- resolve_class1_colours(frame_record)

    # Display rule: top-panel legend lists Defining_SNP_Sites only.
    # The Mutation_Sites set is still tracked in the per-frame CSV
    # (write_cluster_membership_csv) so users can filter downstream.
    defining <- cell$snp_sites_truth %||% integer(0L)
    leg_list <- build_cluster_legend(frame_record, defining_sites = defining)

    fg_df         <- make_traj_df(cls_mat, bin_idx_vec, fb)
    # Map each long-form row to its site's colour via named-vector lookup
    fg_df$colour  <- unname(site_cols[fg_df$site])
  }

  y_max <- max(c(bg_df$y, fg_df$y, 0.01), na.rm = TRUE)

  p <- ggplot2::ggplot()

  # Background: 1273 grey curves
  p <- p + ggplot2::geom_line(
    data = bg_df,
    ggplot2::aes(x = x, y = y, group = site),
    colour = GREY_BG, alpha = 0.55, linewidth = 0.18
  )

  if (!is.null(fg_df) && nrow(fg_df) > 0L) {
    # Boundary semantics: solid up to and including the last data point that
    # falls inside (or before) the fitting window; dashed from that same
    # point onwards. For Hellinger anchored at a bin inside the fitting
    # window (i.e. T_pred, where max(fb) often IS the anchor bin and thus
    # has no value), `boundary_x` is the last bin in the data with x <=
    # max(fb), not max(fb) itself. That keeps the line visually continuous
    # across the missing anchor bin.
    end_fb     <- max(fb)
    data_x     <- sort(unique(fg_df$x))
    boundary_x <- if (any(data_x <= end_fb)) max(data_x[data_x <= end_fb])
                  else min(data_x)

    in_df  <- fg_df[fg_df$x <= boundary_x, , drop = FALSE]   # solid
    out_df <- fg_df[fg_df$x >= boundary_x, , drop = FALSE]   # dashed (shares boundary_x)

    if (nrow(in_df) > 0L)
      p <- p + ggplot2::geom_line(
        data = in_df,
        ggplot2::aes(x = x, y = y, group = site, colour = I(colour)),
        linewidth = 0.8, alpha = 0.95
      )
    if (nrow(out_df) > 0L)
      p <- p + ggplot2::geom_line(
        data = out_df,
        ggplot2::aes(x = x, y = y, group = site, colour = I(colour)),
        linewidth = 0.55, alpha = 0.80, linetype = "22"
      )
  }

  p <- add_window_and_detection(p, cell, frame_record, y_max)

  # Legend subtitle: variant-defining sites grouped by cluster, e.g.
  # "Defining sites by cluster:  c1: 501, 681 | c3: 26, 138"
  if (length(leg_list$labels) > 0L) {
    leg_subtitle <- paste0(
      "Defining sites by cluster:  ",
      paste(leg_list$labels, collapse = "   |   ")
    )
  } else {
    leg_subtitle <- "Defining sites: none clustered"
  }

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = seq_len(cell$n_bins),
      labels = cell$bin_labels,
      limits = c(0.5, cell$n_bins + 0.5),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, y_max * 1.13),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::labs(x = NULL,
                   y = metric_y_label(metric),
                   subtitle = leg_subtitle) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y      = ggplot2::element_text(size = 7),
      axis.title.y     = ggplot2::element_text(size = 9),
      plot.subtitle    = ggplot2::element_text(size = 8, hjust = 0.5,
                                                face = "bold",
                                                colour = "grey15"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(2, 6, 0, 6, "pt")
    )
  p
}


# ── 9. Bottom panel — same plot, SNP sites highlighted red ──────────────────

build_bottom_panel <- function(cell, frame_record, metric, config) {

  bg_mat <- metric_full_matrix(cell, metric)
  if (is.null(bg_mat))
    stop("No full-cell matrix for metric '", metric, "'.", call. = FALSE)

  bin_idx_vec <- matrix_to_bin_indices(bg_mat, metric)
  fb          <- frame_record$frame_bins

  snp_sites <- cell$snp_sites_truth
  if (is.null(snp_sites)) snp_sites <- integer(0L)

  bg_df <- make_traj_df(bg_mat, bin_idx_vec, fb)

  fg_df <- NULL
  snp_in_data <- intersect(snp_sites, as.integer(rownames(bg_mat)))
  if (length(snp_in_data) > 0L) {
    snp_mat <- bg_mat[as.character(snp_in_data), , drop = FALSE]
    fg_df   <- make_traj_df(snp_mat, bin_idx_vec, fb)
    # Per-SNP colour: inherits the cluster colour from the top panel if the
    # SNP ended up in a usable cluster; otherwise red (or grey if oversized).
    snp_cols      <- resolve_snp_colours(snp_in_data, frame_record)
    fg_df$colour  <- unname(snp_cols[fg_df$site])
  }

  y_max <- max(c(bg_df$y, fg_df$y, 0.01), na.rm = TRUE)

  p <- ggplot2::ggplot()

  # Background: 1273 grey curves (includes the SNPs, drawn under the highlights)
  p <- p + ggplot2::geom_line(
    data = bg_df,
    ggplot2::aes(x = x, y = y, group = site),
    colour = GREY_BG, alpha = 0.55, linewidth = 0.18
  )

  if (!is.null(fg_df) && nrow(fg_df) > 0L) {
    # Same data-aware boundary as the top panel (see comment there). Avoids
    # the visible gap when the fitting-window upper bound is the anchor bin
    # of hellinger_Tpred (no value at that bin).
    end_fb     <- max(fb)
    data_x     <- sort(unique(fg_df$x))
    boundary_x <- if (any(data_x <= end_fb)) max(data_x[data_x <= end_fb])
                  else min(data_x)

    in_df  <- fg_df[fg_df$x <= boundary_x, , drop = FALSE]   # solid
    out_df <- fg_df[fg_df$x >= boundary_x, , drop = FALSE]   # dashed (shares boundary_x)

    if (nrow(in_df) > 0L)
      p <- p + ggplot2::geom_line(
        data = in_df,
        ggplot2::aes(x = x, y = y, group = site, colour = I(colour)),
        linewidth = 0.85, alpha = 0.95
      )
    if (nrow(out_df) > 0L)
      p <- p + ggplot2::geom_line(
        data = out_df,
        ggplot2::aes(x = x, y = y, group = site, colour = I(colour)),
        linewidth = 0.55, alpha = 0.80, linetype = "22"
      )
  }

  p <- add_window_and_detection(p, cell, frame_record, y_max)

  # Subtitle stays compact — just the SNP list. The per-SNP cluster mapping
  # (e.g. "18 -> c3") is emitted by write_cluster_membership_csv() to the
  # frame PNG directory so it can be inspected without crowding the plot.
  snp_subtitle <- if (length(snp_in_data) > 0L) {
    snp_str <- if (length(snp_in_data) <= 16L)
                  paste(snp_in_data, collapse = ", ")
                else
                  sprintf("%s ... (%d sites)",
                          paste(snp_in_data[1:12], collapse = ", "),
                          length(snp_in_data))
    sprintf("%s SNPs:  %s", cell$variant, snp_str)
  } else {
    sprintf("%s SNPs: none in matrix", cell$variant)
  }

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = seq_len(cell$n_bins),
      labels = cell$bin_labels,
      limits = c(0.5, cell$n_bins + 0.5),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, y_max * 1.13),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::labs(x = "Periods",
                   y = metric_y_label(metric),
                   subtitle = snp_subtitle) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y      = ggplot2::element_text(size = 7),
      axis.title       = ggplot2::element_text(size = 9),
      plot.subtitle    = ggplot2::element_text(size = 8, hjust = 0.5,
                                                face = "bold",
                                                colour = "grey15"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(0, 6, 4, 6, "pt")
    )
  p
}


# ── 10. Title / subtitle helpers ────────────────────────────────────────────

frame_title <- function(cell, metric) {
  metric_lab <- switch(
    metric,
    entropy         = "Entropy",
    hellinger_T1    = "Hellinger (anchor: T1)",
    hellinger_Tpred = sprintf("Hellinger (anchor: %s)",
                               cell$predecessor_name %||% "T_pred"),
    metric
  )
  sprintf("%s  \u00B7  %s  \u00B7  %s  \u00B7  %s",
          cell$variant, cell$dataset, cell$strategy, metric_lab)
}

frame_subtitle <- function(cell, frame_record, frame_idx, n_frames) {

  fb <- frame_record$frame_bins
  win_lab <- if (length(fb) > 0L)
               sprintf("fitting window bins %d\u2013%d (%s \u2192 %s)",
                       min(fb), max(fb),
                       cell$bin_labels[min(fb)], cell$bin_labels[max(fb)])
             else "fitting window: none"

  # Trimmed subtitle: only frame counter, fitting window, and k.
  # class-1 count, mean silhouette, and oversized-cluster flags live in
  # the per-frame CSV (write_cluster_membership_csv).
  status_lab <- if (identical(frame_record$fit_status, "ok"))
                  sprintf("k = %d", frame_record$n_clusters)
                else
                  sprintf("status: %s", frame_record$fit_status)

  sprintf("frame %d/%d  \u00B7  %s  \u00B7  %s",
          frame_idx, n_frames, win_lab, status_lab)
}


# ── 11. Combine into the dual-panel figure ──────────────────────────────────

build_dual_panel <- function(top, bottom, title_str, subtitle_str = NULL,
                              config) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("patchwork required.", call. = FALSE)

  tf <- config$PANEL_TOP_FRAC %||% 0.5
  out <- patchwork::wrap_plots(top, bottom, ncol = 1L,
                                 heights = c(tf, 1 - tf)) +
         patchwork::plot_annotation(
           title    = title_str,
           subtitle = subtitle_str,
           theme    = ggplot2::theme(
             plot.title    = ggplot2::element_text(face = "bold", size = 11,
                                                    hjust = 0.5),
             plot.subtitle = ggplot2::element_text(size = 9, colour = "grey30",
                                                    hjust = 0.5),
             plot.margin   = ggplot2::margin(6, 6, 4, 6, "pt")
           )
         )
  out
}


# ── 12. Single-frame composite ──────────────────────────────────────────────

build_frame_plot <- function(cell, frame_record, metric, frame_idx, n_frames,
                              config) {
  top    <- build_top_panel   (cell, frame_record, metric, config)
  bottom <- build_bottom_panel(cell, frame_record, metric, config)
  build_dual_panel(
    top, bottom,
    title_str    = frame_title(cell, metric),
    subtitle_str = frame_subtitle(cell, frame_record, frame_idx, n_frames),
    config       = config
  )
}


# ── 12b. Per-frame cluster membership CSV writer ────────────────────────────
# ── 12b. Variant catalogue helpers (used by the CSV writer) ─────────────────
#
# build_variants_catalogue() turns the package's `sarscov2_variants` data
# frame (one row per variant, list-columns of integer site positions) into
# a flat named list:
#
#   list(
#     Alpha   = list(snp = c(...), mutation = c(...)),
#     Beta    = list(snp = c(...), mutation = c(...)),
#     ...
#   )
#
# Robust to both list-column and atomic representations. Pass the result
# to write_cluster_membership_csv() via the `variants_catalogue` argument
# to populate per-site and per-cluster variant-attribution columns.

# build_variants_catalogue() handles BOTH possible structures of the input:
#
#   (a) A data.frame with one row per variant and list-columns
#       Defining_SNP_Sites / Mutation_Sites.
#   (b) The package's actual structure: a NAMED LIST where each element is a
#       column-vector of length n_variants (e.g. sarscov2_variants$WHO_Label
#       is chr(12), sarscov2_variants$Defining_SNP_Sites is a list of 12
#       integer vectors). This is what gets exported as
#       `data(sarscov2_variants, package = "ViralEntropR")`.
#
# Both paths converge to the same output: a named list keyed by variant
# label, each entry list(snp = integer, mutation = integer).

build_variants_catalogue <- function(variants_data) {
  if (is.null(variants_data)) return(NULL)

  # Both list and data.frame support `$col` extraction the same way.
  who_lab   <- variants_data$WHO_Label %||% variants_data$Variant_Name %||%
                 variants_data$Name
  snp_col   <- variants_data$Defining_SNP_Sites %||% variants_data$SNP_Sites %||%
                 variants_data$snp_sites
  mut_col   <- variants_data$Mutation_Sites %||% variants_data$mutation_sites

  if (is.null(who_lab) || length(who_lab) == 0L) {
    warning("build_variants_catalogue: no WHO_Label / Variant_Name column.",
            call. = FALSE)
    return(NULL)
  }

  extract_sites <- function(val) {
    if (is.null(val) || length(val) == 0L) return(integer(0L))
    if (is.list(val)) val <- unlist(val, use.names = FALSE)
    val <- val[!is.na(val)]
    if (length(val) == 0L) return(integer(0L))
    as.integer(val)
  }

  n <- length(who_lab)
  catalogue <- vector("list", n)
  names(catalogue) <- as.character(who_lab)
  for (i in seq_len(n)) {
    snp_i <- if (!is.null(snp_col)) extract_sites(snp_col[[i]]) else integer(0L)
    mut_i <- if (!is.null(mut_col)) extract_sites(mut_col[[i]]) else integer(0L)
    catalogue[[i]] <- list(snp = snp_i, mutation = mut_i)
  }
  catalogue
}

# Per-site lookup: for a vector of sites, return a character vector where
# entry i is a comma-separated list of variant names whose `key` set
# (either "snp" or "mutation") contains sites[i]. Empty string if none.

.variants_containing_site <- function(sites, catalogue, key) {
  if (is.null(catalogue) || length(sites) == 0L)
    return(rep("", length(sites)))
  vapply(sites, function(s) {
    hits <- character(0L)
    for (v in names(catalogue)) {
      if (s %in% catalogue[[v]][[key]]) hits <- c(hits, v)
    }
    paste(hits, collapse = ",")
  }, character(1L))
}

# Per-cluster majority: count, across all variants, how many sites in the
# given cluster belong to each variant's `key` set. Return the variant(s)
# with the maximum count as a comma-joined string. Returns "" if no
# variant has ANY sites in this cluster (cluster is variant-unrelated).

.cluster_majority_variant <- function(sites_in_cluster, catalogue, key) {
  if (is.null(catalogue) || length(sites_in_cluster) == 0L)
    return("")
  counts <- vapply(catalogue, function(vsets)
                   length(intersect(sites_in_cluster, vsets[[key]])),
                   integer(1L))
  if (max(counts) == 0L) return("")
  paste(names(counts)[counts == max(counts)], collapse = ",")
}


# ── 12c. Per-frame cluster membership CSV writer ────────────────────────────
#
# Writes one row per (frame, class-1 site) for the given metric. The CSV
# is INTENTIONALLY VARIANT-AGNOSTIC: the same content is valid for every
# variant of the same (dataset, strategy) base cell (and, for
# hellinger_Tpred, every variant sharing the same predecessor). Columns:
#
#   frame_idx                     integer
#   centre_bin                    integer
#   metric                        character
#   n_class1_sites                integer  (frame-level summary)
#   mean_silhouette               numeric  (frame-level summary)
#   cluster_id                    integer
#   cluster_size                  integer
#   site                          integer
#   per_site_silhouette           numeric
#   is_oversized_cluster          logical
#   is_snp                        logical  (in focal-variant SNP set)
#   is_mutation                   logical  (in focal-variant Mutation set)
#   variants_with_snp             character (comma-joined variant names)
#   variants_with_mutation        character (comma-joined variant names)
#   majority_variant_by_snp       character (per-cluster, broadcast)
#   majority_variant_by_mutation  character (per-cluster, broadcast)
#
# Frames with fit_status != "ok" are written with NA cluster_id / NA
# silhouette to preserve the timeline.

write_cluster_membership_csv <- function(cell, frame_records, metric,
                                          output_path,
                                          variants_catalogue = NULL) {
  if (length(frame_records) == 0L)
    return(invisible(NULL))

  snp_truth <- cell$snp_sites_truth      %||% integer(0L)
  mut_truth <- cell$mutation_sites_truth %||% integer(0L)

  empty_row <- function(i, f) data.frame(
    frame_idx                    = i,
    centre_bin                   = f$centre_bin %||% NA_integer_,
    metric                       = metric,
    n_class1_sites               = NA_integer_,
    mean_silhouette              = NA_real_,
    cluster_id                   = NA_integer_,
    cluster_size                 = NA_integer_,
    site                         = NA_integer_,
    per_site_silhouette          = NA_real_,
    is_oversized_cluster         = NA,
    is_snp                       = NA,
    is_mutation                  = NA,
    variants_with_snp            = NA_character_,
    variants_with_mutation       = NA_character_,
    majority_variant_by_snp      = NA_character_,
    majority_variant_by_mutation = NA_character_,
    stringsAsFactors             = FALSE
  )

  rows <- vector("list", length(frame_records))
  for (i in seq_along(frame_records)) {
    f <- frame_records[[i]]

    if (!identical(f$fit_status, "ok") || f$n_class1_sites == 0L) {
      rows[[i]] <- empty_row(i, f)
      next
    }

    sites <- f$class1_sites
    cls   <- f$cluster_assignments
    sil   <- f$per_site_sil
    over  <- f$oversized

    n     <- length(sites)
    # Per-cluster aggregates broadcast to row level
    cluster_size_vec  <- integer(n)
    majority_snp_vec  <- character(n)
    majority_mut_vec  <- character(n)
    for (k in unique(stats::na.omit(cls))) {
      idx_k    <- which(cls == k)
      sites_k  <- sites[idx_k]
      cluster_size_vec[idx_k] <- length(sites_k)
      majority_snp_vec[idx_k] <-
        .cluster_majority_variant(sites_k, variants_catalogue, "snp")
      majority_mut_vec[idx_k] <-
        .cluster_majority_variant(sites_k, variants_catalogue, "mutation")
    }

    rows[[i]] <- data.frame(
      frame_idx                    = i,
      centre_bin                   = f$centre_bin %||% NA_integer_,
      metric                       = metric,
      n_class1_sites               = f$n_class1_sites,
      mean_silhouette              = f$silhouette_mean,
      cluster_id                   = cls,
      cluster_size                 = cluster_size_vec,
      site                         = sites,
      per_site_silhouette          = sil,
      is_oversized_cluster         = vapply(
        cls,
        function(k) if (is.na(k)) NA else isTRUE(over[k]),
        logical(1L)
      ),
      is_snp                       = sites %in% snp_truth,
      is_mutation                  = sites %in% mut_truth,
      variants_with_snp            =
        .variants_containing_site(sites, variants_catalogue, "snp"),
      variants_with_mutation       =
        .variants_containing_site(sites, variants_catalogue, "mutation"),
      majority_variant_by_snp      = majority_snp_vec,
      majority_variant_by_mutation = majority_mut_vec,
      stringsAsFactors             = FALSE
    )
  }

  df <- do.call(rbind, rows)
  # Order rows: frame_idx ASC, then cluster_id ASC (NAs last), then site ASC.
  # Groups all members of cluster 1 together, then cluster 2, etc., within
  # each frame — easier to read by cluster.
  ord <- order(
    df$frame_idx,
    is.na(df$cluster_id), df$cluster_id,
    df$site,
    na.last = TRUE
  )
  df <- df[ord, , drop = FALSE]
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, output_path, row.names = FALSE)
  invisible(output_path)
}


# ── 13. GIF assembler ───────────────────────────────────────────────────────
#
# Streams all frames through a single magick::image_graph virtual device
# and writes the animation.

assemble_gif <- function(cell, frame_records, metric, output_path, config,
                          frames_dir = NULL) {

  if (!requireNamespace("magick", quietly = TRUE))
    stop("magick package required.", call. = FALSE)
  if (length(frame_records) == 0L) {
    log_msg(sprintf("assemble_gif: no frames for %s / %s.",
                    cell$variant, metric))
    return(invisible(NULL))
  }

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  if (!is.null(frames_dir))
    dir.create(frames_dir, recursive = TRUE, showWarnings = FALSE)

  img <- magick::image_graph(
    width  = config$GIF_WIDTH_PX,
    height = config$GIF_HEIGHT_PX,
    res    = config$GIF_RES,
    bg     = "white"
  )
  on.exit({
    if (length(dev.list()) > 0L) try(dev.off(), silent = TRUE)
  }, add = TRUE)

  n <- length(frame_records)
  for (i in seq_len(n)) {
    p <- tryCatch(
      build_frame_plot(cell, frame_records[[i]], metric, i, n, config),
      error = function(e) {
        ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                           label = sprintf("frame %d failed: %s",
                                            i, conditionMessage(e))) +
        ggplot2::theme_void()
      }
    )

    # Optionally save the frame as a standalone PNG before printing it
    # into the GIF stream.
    if (!is.null(frames_dir)) {
      centre <- frame_records[[i]]$centre_bin
      tag    <- if (is.null(centre) || is.na(centre))
                  sprintf("%02d", i)
                else
                  sprintf("%03d", as.integer(centre))
      png_path <- file.path(frames_dir,
                             sprintf("frame_%s_%s.png", metric, tag))
      tryCatch(
        ggplot2::ggsave(png_path, p,
                         width  = config$PNG_WIDTH_CM,
                         height = config$PNG_HEIGHT_CM,
                         units  = "cm", dpi = 300, bg = "white"),
        error = function(e)
          log_msg(sprintf("assemble_gif: frame PNG write failed (%s): %s",
                          png_path, conditionMessage(e)))
      )
    }

    print(p)
  }

  dev.off()
  on.exit()

  ani <- magick::image_animate(img, fps = config$ANIMATION_FPS, loop = 0L)
  magick::image_write(ani, path = output_path, format = "gif")
  invisible(output_path)
}
