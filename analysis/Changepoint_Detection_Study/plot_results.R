################################################################################
## plot_results.R — Changepoint_Detection_Study
##
## Produces eight publication-quality figures:
##
##   fig_01_f1_by_method_tier_dataset    : 6-panel master (3 tiers × 2 datasets)
##   fig_02_operating_curve              : F1 vs truths_effective, crossover plot
##   fig_03_k_sweep_saturation           : K vs F1, faceted by dataset × tier
##   fig_04_full_vs_reduced_scatter      : per-replicate full-vs-reduced F1
##   fig_05_failure_heatmap              : method × tier failure rates
##   fig_06_cross_dataset_rank_scatter   : NCBI rank vs GISAID rank, Spearman ρ
##   fig_07_timeline_NCBI_US             : detected vs true CPs, faceted by method
##   fig_08_timeline_GISAID_US           : detected vs true CPs, faceted by method
##
## All figures: PNG only (300 DPI), 17 cm wide × variable height,
## journal-ready typography, captionless (captions managed in LaTeX).
##
## Inputs:
##   outputs/summary_<dataset>.rds                  (5,001-window benchmark)
##   outputs/full_dataset_detection_<dataset>.rds   (single full-dataset run)
##
## Outputs:
##   outputs/plots/fig_NN_<name>.png
################################################################################

suppressPackageStartupMessages({
  source("setup.R")
  source("helpers_windows.R")   # add_months, diff_months, bin_midpoint helpers
  library(ggplot2)
  library(scales)
})

config   <- build_config()
plot_dir <- file.path(output_dir(config), "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Aesthetic constants ───────────────────────────────────────────────────

PALETTE_METHOD <- c(
  e_agglo          = "#1F77B4",   # blue   — K-free ECP
  ks_cp3o_dynamic  = "#D62728",   # red    — primary K-required ECP
  ks_cp3o_K1       = "#FF9896",   # red    — K-sweep variants (lighter)
  ks_cp3o_K2       = "#E377C2",
  ks_cp3o_K3       = "#FFBB78",
  ks_cp3o_K5       = "#FF7F0E",
  ks_cp3o_K10      = "#8C564B",
  hdcp_binseg      = "#2CA02C",   # green  — HDcp BS
  hdcp_wbs         = "#7F7F7F"    # grey   — HDcp WBS
)

METHOD_LABELS <- c(
  e_agglo          = "e_agglo",
  ks_cp3o_dynamic  = "ks_cp3o (K=dyn)",
  ks_cp3o_K1       = "ks_cp3o (K=1)",
  ks_cp3o_K2       = "ks_cp3o (K=2)",
  ks_cp3o_K3       = "ks_cp3o (K=3)",
  ks_cp3o_K5       = "ks_cp3o (K=5)",
  ks_cp3o_K10      = "ks_cp3o (K=10)",
  hdcp_binseg      = "hdcp_binseg",
  hdcp_wbs         = "hdcp_wbs"
)

METHOD_ORDER <- c("e_agglo", "ks_cp3o_dynamic",
                  "ks_cp3o_K1", "ks_cp3o_K2", "ks_cp3o_K3",
                  "ks_cp3o_K5", "ks_cp3o_K10",
                  "hdcp_binseg", "hdcp_wbs")

TIER_ORDER  <- c("short", "medium", "long")
TIER_LABELS <- c(short = "Short (6–10 mo)",
                 medium = "Medium (12–16 mo)",
                 long   = "Long (18+ mo)")

DATASET_LABELS <- c(NCBI_US = "NCBI US", GISAID_US = "GISAID US")

# Two-letter glyphs for the 12 WHO VOC/VOI variants (Q2 design choice B).
VARIANT_GLYPH <- c(
  Alpha   = "Al", Beta  = "Be", Gamma   = "Ga", Delta = "De",
  Epsilon = "Ep", Eta   = "Et", Iota    = "Io", Kappa = "Ka",
  Lambda  = "La", Omicron = "Om", Theta = "Th", Zeta  = "Ze"
)

# Publication theme: serif, no grid, minimal chrome, ggplot-standard but tight.
theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = "") %+replace%
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey85", linewidth = 0.3),
      panel.border       = element_rect(colour = "grey20",
                                        fill = NA, linewidth = 0.4),
      strip.background   = element_rect(fill = "grey95", colour = NA),
      strip.text         = element_text(size = base_size, face = "bold"),
      axis.text          = element_text(size = base_size - 1L),
      axis.title         = element_text(size = base_size),
      legend.title       = element_text(size = base_size - 1L),
      legend.text        = element_text(size = base_size - 1L),
      legend.key.size    = unit(0.8, "lines"),
      plot.title         = element_blank(),
      plot.subtitle      = element_blank()
    )
}

# Save helper: write PNG (300 DPI) at journal width.
save_plot <- function(p, name, width_cm = 17, height_cm = 11) {
  out <- file.path(plot_dir, sprintf("%s.png", name))
  ggsave(filename = out, plot = p,
         width = width_cm, height = height_cm, units = "cm",
         dpi = 300, device = "png", bg = "white")
  log_msg(sprintf("Saved: %s", out))
  invisible(NULL)
}

# ── 2. Load benchmark summaries ──────────────────────────────────────────────

log_msg("Loading benchmark summary tables...")
sN <- readRDS(summary_path(config, "NCBI_US"))
sG <- readRDS(summary_path(config, "GISAID_US"))

# Factor ordering for plots.
prep_summary <- function(df, dataset_label) {
  df$dataset <- dataset_label
  df$method  <- factor(df$method, levels = METHOD_ORDER)
  df$tier    <- factor(df$tier,   levels = TIER_ORDER)
  df
}
sN <- prep_summary(sN, "NCBI US")
sG <- prep_summary(sG, "GISAID US")
sAll <- rbind(sN, sG)
sAll$dataset <- factor(sAll$dataset, levels = c("NCBI US", "GISAID US"))

# ── 3. Figure 1: F1 by method × tier × dataset (the master figure) ───────────

log_msg("Building figure 01: F1 by method × tier × dataset...")
df_f1 <- subset(sAll, site_set == "full")

fig_01 <- ggplot(df_f1, aes(x = method, y = F1, fill = method)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3,
               linewidth = 0.3, na.rm = TRUE) +
  facet_grid(dataset ~ tier,
             labeller = labeller(tier = TIER_LABELS)) +
  scale_fill_manual(values = PALETTE_METHOD, labels = METHOD_LABELS,
                    guide  = "none") +
  scale_x_discrete(labels = METHOD_LABELS) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = expression(F[1])) +
  theme_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

save_plot(fig_01, "fig_01_f1_by_method_tier_dataset",
          width_cm = 17, height_cm = 12)

# ── 4. Figure 2: F1 vs truths_effective (operating curve) ────────────────────

log_msg("Building figure 02: operating curve F1 vs truths_effective...")
focal_methods <- c("e_agglo", "ks_cp3o_dynamic")
df_op <- subset(sAll, site_set == "full" & method %in% focal_methods)
df_op$method <- droplevels(df_op$method)

agg_op <- aggregate(F1 ~ method + truths_effective + dataset, df_op,
                    function(x) c(mean = mean(x, na.rm = TRUE),
                                  se   = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))),
                                  n    = sum(!is.na(x))))
agg_op <- do.call(data.frame,
                  c(agg_op[, 1:3], list(F1 = agg_op$F1)))
names(agg_op) <- c("method", "truths_effective", "dataset",
                   "F1_mean", "F1_se", "F1_n")

fig_02 <- ggplot(agg_op,
                 aes(x = truths_effective, y = F1_mean,
                     colour = method, fill = method)) +
  geom_ribbon(aes(ymin = F1_mean - 1.96 * F1_se,
                  ymax = F1_mean + 1.96 * F1_se),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2, shape = 21, colour = "white", stroke = 0.5) +
  facet_wrap(~ dataset) +
  scale_colour_manual(values = PALETTE_METHOD, labels = METHOD_LABELS,
                      name = NULL) +
  scale_fill_manual(values   = PALETTE_METHOD, guide = "none") +
  scale_x_continuous(breaks = function(lim) seq(ceiling(lim[1]),
                                                floor(lim[2]), 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Unique truth bins per window", y = expression(F[1])) +
  theme_pub() +
  theme(legend.position = "top")

save_plot(fig_02, "fig_02_operating_curve",
          width_cm = 17, height_cm = 9)

# ── 5. Figure 3: K-sweep saturation ──────────────────────────────────────────

log_msg("Building figure 03: K-sweep saturation...")
ks_methods <- c("ks_cp3o_K1", "ks_cp3o_K2", "ks_cp3o_K3",
                "ks_cp3o_K5", "ks_cp3o_K10")
df_ks <- subset(sAll, site_set == "full" & method %in% ks_methods)
df_ks$method <- droplevels(df_ks$method)
df_ks$K <- as.integer(sub("ks_cp3o_K", "", as.character(df_ks$method)))

agg_ks <- aggregate(F1 ~ K + tier + dataset, df_ks,
                    function(x) c(mean = mean(x, na.rm = TRUE),
                                  se   = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))))
agg_ks <- do.call(data.frame, c(agg_ks[, 1:3], list(F1 = agg_ks$F1)))
names(agg_ks) <- c("K", "tier", "dataset", "F1_mean", "F1_se")

fig_03 <- ggplot(agg_ks, aes(x = K, y = F1_mean)) +
  geom_ribbon(aes(ymin = F1_mean - 1.96 * F1_se,
                  ymax = F1_mean + 1.96 * F1_se),
              alpha = 0.18, fill = "#D62728", colour = NA) +
  geom_line(linewidth = 0.7, colour = "#D62728") +
  geom_point(size = 2, shape = 21, fill = "#D62728",
             colour = "white", stroke = 0.5) +
  facet_grid(dataset ~ tier,
             labeller = labeller(tier = TIER_LABELS)) +
  scale_x_continuous(breaks = c(1, 2, 3, 5, 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = expression(italic(K)~"(max change points)"),
       y = expression(F[1])) +
  theme_pub()

save_plot(fig_03, "fig_03_k_sweep_saturation",
          width_cm = 17, height_cm = 11)

# ── 6. Figure 4: Full vs reduced scatter ─────────────────────────────────────

log_msg("Building figure 04: full vs reduced site-set scatter...")
# Pair full and reduced F1 per (cell_id, dataset, method).
pair_full_reduced <- function(df) {
  full <- df[df$site_set == "full",
             c("dataset", "cell_id", "method", "F1")]
  red  <- df[df$site_set == "reduced",
             c("dataset", "cell_id", "method", "F1")]
  m <- merge(full, red, by = c("dataset", "cell_id", "method"),
             suffixes = c("_full", "_reduced"))
  m
}
df_pair <- pair_full_reduced(sAll)
df_pair$method  <- factor(df_pair$method, levels = METHOD_ORDER)
df_pair$dataset <- factor(df_pair$dataset, levels = c("NCBI US", "GISAID US"))

fig_04 <- ggplot(df_pair, aes(x = F1_full, y = F1_reduced)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "grey50", linetype = "dashed", linewidth = 0.3) +
  geom_point(aes(colour = method), alpha = 0.15, size = 0.4,
             show.legend = FALSE) +
  facet_grid(dataset ~ method,
             labeller = labeller(method = METHOD_LABELS)) +
  scale_colour_manual(values = PALETTE_METHOD) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  coord_fixed() +
  labs(x = expression(F[1]~"(full 1,273 sites)"),
       y = expression(F[1]~"(GMM class-1 reduced)")) +
  theme_pub(base_size = 8) +
  theme(strip.text.x = element_text(size = 7),
        axis.text    = element_text(size = 7))

save_plot(fig_04, "fig_04_full_vs_reduced_scatter",
          width_cm = 22, height_cm = 7)

# ── 7. Figure 5: Failure-rate heatmap ────────────────────────────────────────

log_msg("Building figure 05: failure-rate heatmap...")
df_fail <- aggregate(status ~ method + tier + dataset,
                     subset(sAll, site_set == "full"),
                     function(x) mean(x == "failed"))
names(df_fail)[4L] <- "failure_rate"
df_fail$method  <- factor(df_fail$method, levels = rev(METHOD_ORDER))
df_fail$tier    <- factor(df_fail$tier,   levels = TIER_ORDER)
df_fail$dataset <- factor(df_fail$dataset, levels = c("NCBI US", "GISAID US"))

fig_05 <- ggplot(df_fail,
                 aes(x = tier, y = method, fill = failure_rate)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * failure_rate)),
            colour = "black", size = 3) +
  facet_wrap(~ dataset) +
  scale_fill_gradient(low = "#FFFFFF", high = "#D62728",
                      limits = c(0, 1),
                      labels = scales::percent_format(accuracy = 1),
                      name   = "Failure rate") +
  scale_x_discrete(labels = TIER_LABELS) +
  scale_y_discrete(labels = METHOD_LABELS) +
  labs(x = NULL, y = NULL) +
  theme_pub() +
  theme(legend.position = "right",
        panel.grid      = element_blank(),
        axis.text.x     = element_text(angle = 30, hjust = 1, vjust = 1))

save_plot(fig_05, "fig_05_failure_heatmap",
          width_cm = 18, height_cm = 9)

# ── 8. Figure 6: Cross-dataset rank scatter ──────────────────────────────────

log_msg("Building figure 06: cross-dataset rank scatter...")
agg_method_dataset <- aggregate(F1 ~ method + dataset,
                                subset(sAll, site_set == "full"),
                                function(x) mean(x, na.rm = TRUE))
agg_wide <- merge(
  subset(agg_method_dataset, dataset == "NCBI US",   select = c(method, F1)),
  subset(agg_method_dataset, dataset == "GISAID US", select = c(method, F1)),
  by = "method", suffixes = c("_NCBI", "_GISAID")
)
agg_wide$rank_NCBI   <- rank(-agg_wide$F1_NCBI)
agg_wide$rank_GISAID <- rank(-agg_wide$F1_GISAID)

spearman_rho <- cor(agg_wide$rank_NCBI, agg_wide$rank_GISAID,
                    method = "spearman")

.have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
label_layer <- if (.have_ggrepel) {
  ggrepel::geom_text_repel(aes(label = METHOD_LABELS[as.character(method)]),
                           size = 3, colour = "black",
                           box.padding = 0.4,
                           max.overlaps = Inf,
                           segment.color = "grey60",
                           segment.size = 0.3)
} else {
  log_msg("ggrepel not installed; falling back to geom_text for figure 06.")
  geom_text(aes(label = METHOD_LABELS[as.character(method)]),
            size = 3, colour = "black",
            hjust = -0.15, vjust = -0.5)
}

fig_06 <- ggplot(agg_wide,
                 aes(x = rank_NCBI, y = rank_GISAID, colour = method)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "grey50", linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 4) +
  label_layer +
  annotate("text",
           x = 1.2, y = 9,
           label = sprintf("Spearman~rho == %.3f", spearman_rho),
           parse = TRUE, hjust = 0, size = 3.6,
           family = "") +
  scale_colour_manual(values = PALETTE_METHOD,
                      labels = METHOD_LABELS,
                      guide  = "none") +
  scale_x_continuous(breaks = 1:9, limits = c(0.5, 9.5)) +
  scale_y_continuous(breaks = 1:9, limits = c(0.5, 9.5)) +
  coord_fixed() +
  labs(x = "Rank on NCBI US (1 = best)",
       y = "Rank on GISAID US (1 = best)") +
  theme_pub()

save_plot(fig_06, "fig_06_cross_dataset_rank_scatter",
          width_cm = 14, height_cm = 14)

# ── 9. Figures 7 & 8: per-method timeline (NCBI / GISAID) ────────────────────
#
# For each dataset:
#   - Load outputs/full_dataset_detection_<dataset>.rds
#   - Map detected CPs to mid-bin dates (bin k → midpoint of bin k+1, since
#     Hellinger row k corresponds to bin T_{k+1}).
#   - Map truth dates: each variant glyph at its bin's midpoint date.
#   - Faceted by method (1 panel per method, vertical stack).
#   - Stars (detected) + letter glyphs (truth), stacked vertically when
#     multiple truths share a month.

bin_midpoint <- function(bin_starts, k) {
  # k is a 1-indexed Hellinger row → maps to bin (k + 1) of the n-bin grid.
  # Bin k+1 spans [bin_starts[k+1], bin_starts[k+2]); midpoint is the date
  # one month after bin_starts[k+1] for 2-month bins.
  start_k <- bin_starts[k + 1L]
  add_months(start_k, 1L)
}

# Build the x-axis tick locations and partition-range labels for a timeline.
# bin_starts has length n_bins + 1 (the partition boundaries). Each bin k
# spans [bin_starts[k], bin_starts[k+1]); we place the tick at the bin's
# midpoint and label it with the partition range, e.g. "Mar 2020 - Apr 2020".
partition_breaks_and_labels <- function(bin_starts) {
  n_bins <- length(bin_starts) - 1L
  starts <- bin_starts[seq_len(n_bins)]
  ends   <- add_months(starts, 1L)              # second month of each 2-month bin
  mids   <- add_months(starts, 1L)              # midpoint = start of second month
  labels <- sprintf("%s - %s",
                    format(starts, "%b %Y"),
                    format(ends,   "%b %Y"))
  list(breaks = mids, labels = labels)
}

build_timeline_data <- function(detection_rds,
                                site_set = c("full", "reduced")) {
  
  site_set     <- match.arg(site_set)
  metrics_blk  <- if (site_set == "full") detection_rds$metrics
  else                    detection_rds$metrics_reduced
  
  # If the reduced metrics block is absent (GMM skip), return empty frames.
  if (is.null(metrics_blk)) {
    return(list(
      detected = data.frame(method = character(0L),
                            date   = as.Date(character(0L)),
                            stringsAsFactors = FALSE),
      truth    = data.frame(date = as.Date(character(0L)),
                            label = character(0L),
                            glyph = character(0L),
                            stack_level = integer(0L),
                            stringsAsFactors = FALSE)
    ))
  }
  
  # Detected CPs per method.
  detected_rows <- list()
  for (method_key in names(metrics_blk)) {
    m <- metrics_blk[[method_key]]
    if (length(m$detected_cps) == 0L) next
    detected_rows[[length(detected_rows) + 1L]] <- data.frame(
      method = method_key,
      date   = bin_midpoint(detection_rds$bin_starts, m$detected_cps),
      stringsAsFactors = FALSE
    )
  }
  detected_df <- if (length(detected_rows) == 0L) {
    data.frame(method = character(0L), date = as.Date(character(0L)))
  } else {
    do.call(rbind, detected_rows)
  }
  
  # True CPs: one row per variant. Map each variant's date to its Hellinger
  # bin's midpoint (NOT the raw date), so detections and truths share a
  # common temporal frame. Variants in T1 are dropped.
  truth_dates_in <- detection_rds$truth_dates
  truth_labels   <- detection_rds$truth_labels
  if (length(truth_dates_in) == 0L) {
    truth_df <- data.frame(date = as.Date(character(0L)),
                           label = character(0L),
                           glyph = character(0L),
                           stringsAsFactors = FALSE)
  } else {
    # Drop T1 truths (no Hellinger row → no detectable CP).
    months_from_start <- diff_months(truth_dates_in,
                                     detection_rds$window_start)
    bin_idx <- (months_from_start %/% config$BIN_MONTHS) + 1L
    drop_T1 <- bin_idx == 1L
    kept_bin_idx <- bin_idx[!drop_T1]
    kept_labels  <- truth_labels[!drop_T1]
    
    # Each truth → midpoint of its bin (bin_idx-th bin out of n_bins).
    truth_date_mapped <- add_months(detection_rds$bin_starts[kept_bin_idx], 1L)
    
    # Stack vertically when multiple variants share a month — assign stack
    # level by within-date sorted order.
    truth_df <- data.frame(date  = truth_date_mapped,
                           label = kept_labels,
                           glyph = VARIANT_GLYPH[kept_labels],
                           stringsAsFactors = FALSE)
    truth_df <- truth_df[order(truth_df$date, truth_df$label), , drop = FALSE]
    truth_df$stack_level <- ave(seq_len(nrow(truth_df)),
                                truth_df$date,
                                FUN = seq_along)
  }
  
  list(detected = detected_df, truth = truth_df)
}

build_timeline_figure <- function(detection_rds, dataset_label,
                                  site_set = c("full", "reduced")) {
  
  site_set <- match.arg(site_set)
  td       <- build_timeline_data(detection_rds, site_set = site_set)
  
  # Restrict to the 4 primary methods + dynamic-K only (drop K-sweep variants
  # to keep the figure readable). The K-sweep already has its own panel (fig 3).
  primary_methods <- c("e_agglo", "ks_cp3o_dynamic", "hdcp_binseg", "hdcp_wbs")
  td$detected <- td$detected[td$detected$method %in% primary_methods, ,
                             drop = FALSE]
  td$detected$method <- factor(td$detected$method, levels = primary_methods)
  
  # X-axis range: the dataset's actual window.
  x_min <- detection_rds$window_start
  x_max <- detection_rds$window_end
  
  # One row per method panel. Truth glyphs in a strip at the BOTTOM
  # (replicated across panels for context), stars within each panel.
  # Stack level → y-offset for the truth strip.
  
  # Compute glyph y-positions: just below y = 0, stacked downward.
  glyph_base_y <- -0.1
  glyph_step_y <- -0.18
  truth_plot <- td$truth
  if (nrow(truth_plot) > 0L) {
    truth_plot$y <- glyph_base_y + (truth_plot$stack_level - 1L) * glyph_step_y
  }
  
  # Star plot data: stars at y = 1 within each method panel.
  detected_plot <- td$detected
  if (nrow(detected_plot) > 0L) {
    detected_plot$y <- 1
  }
  
  # X-axis breaks at bin midpoints, labels as partition ranges
  # (e.g. "Mar 2020 - Apr 2020"). Rotated vertically (90 deg) so 25 bins
  # fit comfortably even on the GISAID timeline.
  x_axis_spec <- partition_breaks_and_labels(detection_rds$bin_starts)
  
  # Build the plot.
  p <- ggplot() +
    # Star markers for detected CPs.
    {
      if (nrow(detected_plot) > 0L) {
        geom_point(data = detected_plot,
                   aes(x = date, y = y, colour = method),
                   shape = 8, size = 3, stroke = 0.9,
                   show.legend = FALSE)
      } else NULL
    } +
    # Bottom truth strip (replicated across panels).
    {
      if (nrow(truth_plot) > 0L) {
        geom_text(data = truth_plot,
                  aes(x = date, y = y, label = glyph),
                  size = 2.6, family = "mono", colour = "black")
      } else NULL
    } +
    # Vertical gridlines at truth dates (subtle).
    {
      if (nrow(truth_plot) > 0L) {
        geom_vline(data = unique(truth_plot[, "date", drop = FALSE]),
                   aes(xintercept = date),
                   linetype = "dotted", linewidth = 0.25, colour = "grey70")
      } else NULL
    } +
    facet_wrap(~ method, ncol = 1,
               labeller = labeller(method = METHOD_LABELS),
               strip.position = "left",
               drop = FALSE) +
    scale_colour_manual(values = PALETTE_METHOD) +
    scale_x_date(limits = c(x_min, x_max),
                 breaks = x_axis_spec$breaks,
                 labels = x_axis_spec$labels,
                 expand = expansion(mult = c(0.01, 0.01))) +
    scale_y_continuous(limits = c(glyph_base_y +
                                    glyph_step_y *
                                    max(1L, max(truth_plot$stack_level,
                                                0L, na.rm = TRUE)), 1.15),
                       breaks = NULL) +
    labs(x = NULL, y = NULL) +
    theme_pub(base_size = 9) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1, vjust = 0.5,
                                     size = 7),
          strip.placement = "outside",
          strip.background.y = element_rect(fill = "grey95", colour = NA),
          strip.text.y.left = element_text(angle = 0,
                                           face = "bold",
                                           size = 8,
                                           hjust = 1))
  
  p
}

# Figure 7: NCBI timeline.
log_msg("Building figure 07: NCBI US timeline...")
det_path_N <- file.path(output_dir(config),
                        "full_dataset_detection_NCBI_US.rds")
if (file.exists(det_path_N)) {
  det_N <- readRDS(det_path_N)
  fig_07 <- build_timeline_figure(det_N, "NCBI US")
  save_plot(fig_07, "fig_07_timeline_NCBI_US",
            width_cm = 19, height_cm = 13)
} else {
  log_msg(sprintf("SKIP figure 07: %s not found. Run run_full_dataset_detection.R first.",
                  det_path_N))
}

# Figure 8: GISAID timeline.
log_msg("Building figure 08: GISAID US timeline (full sites)...")
det_path_G <- file.path(output_dir(config),
                        "full_dataset_detection_GISAID_US.rds")
if (file.exists(det_path_G)) {
  det_G <- readRDS(det_path_G)
  fig_08 <- build_timeline_figure(det_G, "GISAID US", site_set = "full")
  save_plot(fig_08, "fig_08_timeline_GISAID_US",
            width_cm = 19, height_cm = 13)
} else {
  log_msg(sprintf("SKIP figure 08: %s not found. Run run_full_dataset_detection.R first.",
                  det_path_G))
}

# Figure 9: NCBI timeline on GMM-class-1 reduced sites.
log_msg("Building figure 09: NCBI US timeline (GMM-reduced sites)...")
if (file.exists(det_path_N)) {
  if (!is.null(det_N$metrics_reduced)) {
    fig_09 <- build_timeline_figure(det_N, "NCBI US", site_set = "reduced")
    save_plot(fig_09, "fig_09_timeline_NCBI_US_reduced",
              width_cm = 19, height_cm = 13)
    log_msg(sprintf("  Reduced-site context: %d GMM class-1 sites (G=%s, model=%s).",
                    det_N$gmm_meta$n_class1_sites,
                    as.character(det_N$gmm_meta$G),
                    as.character(det_N$gmm_meta$modelName)))
  } else {
    log_msg(sprintf("SKIP figure 09: NCBI reduced-site analysis was skipped (%s).",
                    det_N$reduced_reason))
  }
} else {
  log_msg(sprintf("SKIP figure 09: %s not found.", det_path_N))
}

# Figure 10: GISAID timeline on GMM-class-1 reduced sites.
log_msg("Building figure 10: GISAID US timeline (GMM-reduced sites)...")
if (file.exists(det_path_G)) {
  if (!is.null(det_G$metrics_reduced)) {
    fig_10 <- build_timeline_figure(det_G, "GISAID US", site_set = "reduced")
    save_plot(fig_10, "fig_10_timeline_GISAID_US_reduced",
              width_cm = 19, height_cm = 13)
    log_msg(sprintf("  Reduced-site context: %d GMM class-1 sites (G=%s, model=%s).",
                    det_G$gmm_meta$n_class1_sites,
                    as.character(det_G$gmm_meta$G),
                    as.character(det_G$gmm_meta$modelName)))
  } else {
    log_msg(sprintf("SKIP figure 10: GISAID reduced-site analysis was skipped (%s).",
                    det_G$reduced_reason))
  }
} else {
  log_msg(sprintf("SKIP figure 10: %s not found.", det_path_G))
}

log_msg("==============================================================")
log_msg(sprintf("All figures saved to: %s", plot_dir))
log_msg("==============================================================")