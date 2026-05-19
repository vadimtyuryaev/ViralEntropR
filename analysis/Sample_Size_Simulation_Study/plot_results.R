# =============================================================================
# plot_results.R
# =============================================================================
#
# Publication-quality plot and table generators for the Sample-Size
# Simulation Study. Reads summary_sc<N>.rds from disk, produces
# - per-scenario and per-band tidy summary tables
# - histograms of n_emerge_needed (per scenario, per band, overall)
# - boxplots with notched medians (per scenario, per band, overall)
# - mean-with-CI plots (BCa bootstrap; robust to grid quantization)
# - heatmaps of detection rate across the (N_ref, ratio) grid per scenario
# - faceted detection-rate-vs-n_emerge curves
#
# All plots saved as 300 DPI PNG. Captions intentionally omitted; figure
# captions are managed at the manuscript level in LaTeX.
#
# Statistical notes:
# - Confidence intervals around the mean use BCa bootstrap (1e4 replicates
#   default), not normal-theory CIs. Justification: n_emerge_needed is
#   integer-valued, right-skewed, and grid-quantized; the central limit
#   theorem is unreliable at the sample sizes available per cell (30
#   reps).
# - Notch confidence intervals around the median follow McGill, Tukey &
#   Larsen (1978): median +/- 1.58 * IQR / sqrt(n). Notches that do not
#   overlap suggest a 95%-confidence difference in medians (informal but
#   widely accepted). With 30 reps the notch is occasionally wider than
#   the box, which ggplot flags with a warning and which we suppress.
# - Grid-quantization caveat: reported confidence intervals are about the
#   true (unobservable) detection-threshold distribution. Because the
#   sweep grid has approximately 12% relative resolution in the small
#   band, 17% in the medium band, and 23% in the large band, intervals
#   narrower than this resolution should be interpreted as integer-grid
#   precision limited; narrower than this resolution should be interpreted 
#   as integer-grid precision limited; this caveat should be reported 
#   in the figure caption.
#
# Author : Vadim Tyuryaev
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# -----------------------------------------------------------------------------
# Shared style
# -----------------------------------------------------------------------------
.SCENARIO_LABELS <- c(
  `1` = "Scenario 1",
  `2` = "Scenario 2",
  `3` = "Scenario 3",
  `4` = "Scenario 4"
)

.BAND_ORDER  <- c("small", "medium", "large")
.BAND_COLORS <- c(small = "#1b9e77", medium = "#d95f02", large = "#7570b3")

# Q1-journal-style theme: white background, thin axes, sans-serif
.theme_q1 <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.border     = element_rect(color = "grey20", linewidth = 0.4),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(hjust = 0.5),
      legend.position  = "right"
    )
}

.save_png <- function(plot, path, width_in, height_in, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = path, plot = plot,
         width = width_in, height = height_in,
         dpi = dpi, units = "in", device = "png")
  invisible(path)
}

# -----------------------------------------------------------------------------
# load_all_summaries
# -----------------------------------------------------------------------------
# Reads summary_sc<N>.rds for every scenario that completed, binds them
# into one data frame, and returns the result. Adds a factor band column
# with the canonical level order.
# -----------------------------------------------------------------------------
load_all_summaries <- function(config = NULL, scenarios = 1:4,
                               base_dir = NULL) {
  if (is.null(config)) {
    config <- build_config()
    config$STUDY_DIR <- resolve_study_dir()
  }
  if (is.null(base_dir)) base_dir <- output_dir(config)
  
  dfs <- lapply(scenarios, function(s) {
    p <- file.path(base_dir, sprintf("summary_sc%d.rds", s))
    if (file.exists(p)) readRDS(p) else NULL
  })
  dfs <- dfs[!vapply(dfs, is.null, logical(1L))]
  if (length(dfs) == 0L) stop("No summary RDS files found.", call. = FALSE)
  
  common <- Reduce(intersect, lapply(dfs, colnames))
  out    <- do.call(rbind, lapply(dfs, function(d) d[, common, drop = FALSE]))
  out$band     <- factor(out$band, levels = .BAND_ORDER)
  out$scenario <- as.integer(out$scenario)
  out
}

# -----------------------------------------------------------------------------
# bca_ci_mean
# -----------------------------------------------------------------------------
# BCa bootstrap CI for the mean of an integer-valued, grid-quantized
# sample. Excludes NA values (treats them as right-censored; the user
# should also report the NA rate separately).
# Returns c(lo, mean, hi).
# -----------------------------------------------------------------------------
bca_ci_mean <- function(x, n_boot = 1e4L, conf = 0.95) {
  x <- x[!is.na(x)]
  if (length(x) < 3L) return(c(NA_real_, NA_real_, NA_real_))
  if (!requireNamespace("boot", quietly = TRUE))
    stop("Package 'boot' is required for BCa CIs.", call. = FALSE)
  b <- boot::boot(x, statistic = function(d, i) mean(d[i]), R = n_boot)
  ci <- tryCatch(
    boot::boot.ci(b, type = "bca", conf = conf),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(boot::boot.ci(b, type = "bca",
                                                         conf = conf))
  )
  if (is.null(ci) || is.null(ci$bca)) return(c(NA_real_, mean(x), NA_real_))
  c(lo = ci$bca[4L], mean = mean(x), hi = ci$bca[5L])
}

# -----------------------------------------------------------------------------
# build_results_table
# -----------------------------------------------------------------------------
# Produces a tidy summary table: one row per (scenario, band) plus
# overall rows. Columns: n_total, n_detected, detection_rate, median,
# IQR_low, IQR_high, mean, mean_CI_low, mean_CI_high. Means and CIs
# computed over detected replicates only.
#
# Returns a data frame; pass to knitr::kable, gt, or write to CSV.
# -----------------------------------------------------------------------------
build_results_table <- function(df, conf = 0.95, n_boot = 1e4L,
                                include_km = TRUE) {
  agg <- function(sub) {
    n   <- nrow(sub)
    det <- sub[!is.na(sub$n_emerge_needed), , drop = FALSE]
    x   <- det$n_emerge_needed
    ci  <- if (length(x) >= 3L) bca_ci_mean(x, n_boot = n_boot, conf = conf)
    else c(NA, mean(x, na.rm = TRUE), NA)
    out <- data.frame(
      n_total          = n,
      n_detected       = length(x),
      detect_rate      = length(x) / n,
      median_detected  = if (length(x) > 0L) stats::median(x) else NA_real_,
      mean_detected    = ci[["mean"]],
      mean_ci_low      = ci[["lo"]],
      mean_ci_high     = ci[["hi"]],
      stringsAsFactors = FALSE, check.names = FALSE
    )
    if (include_km && requireNamespace("survival", quietly = TRUE)) {
      sd  <- prepare_survival_data(sub)
      fit <- tryCatch(
        survival::survfit(survival::Surv(time_to_detect, status) ~ 1,
                          data = sd),
        error = function(e) NULL
      )
      if (!is.null(fit)) {
        tbl <- summary(fit)$table
        out$km_median     <- unname(tbl["median"])
        out$km_ci_low     <- unname(tbl["0.95LCL"])
        out$km_ci_high    <- unname(tbl["0.95UCL"])
      } else {
        out$km_median <- out$km_ci_low <- out$km_ci_high <- NA_real_
      }
    }
    out
  }
  
  rows <- list()
  for (s in sort(unique(df$scenario))) {
    for (b in .BAND_ORDER) {
      sub <- df[df$scenario == s & df$band == b, , drop = FALSE]
      if (nrow(sub) == 0L) next
      rows[[length(rows) + 1L]] <- cbind(
        data.frame(scenario = s, band = b, stringsAsFactors = FALSE),
        agg(sub)
      )
    }
    sub <- df[df$scenario == s, , drop = FALSE]
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(scenario = s, band = "(all)", stringsAsFactors = FALSE),
      agg(sub)
    )
  }
  rows[[length(rows) + 1L]] <- cbind(
    data.frame(scenario = NA_integer_, band = "(all)",
               stringsAsFactors = FALSE),
    agg(df)
  )
  
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# -----------------------------------------------------------------------------
# plot_histograms
# -----------------------------------------------------------------------------
# Faceted histogram of n_emerge_needed: rows = scenario, columns = band.
# Bins are integer-respecting; the y-axis is count, not density, because
# replicates are the unit of analysis.
# -----------------------------------------------------------------------------
plot_histograms <- function(df, out_path = "outputs/plots/histograms.png",
                            width_in = 9, height_in = 7) {
  det <- df[!is.na(df$n_emerge_needed), , drop = FALSE]
  p <- ggplot(det, aes(x = n_emerge_needed, fill = band)) +
    geom_histogram(bins = 30, color = "white", linewidth = 0.2) +
    scale_x_log10(breaks = c(2, 10, 100, 1000, 10000),
                  labels = scales::label_comma()) +
    scale_fill_manual(values = .BAND_COLORS, name = "Band") +
    facet_grid(scenario ~ band,
               scales   = "free_y",
               labeller = labeller(
                 scenario = .SCENARIO_LABELS,
                 band     = function(x) paste(toupper(substr(x, 1, 1)),
                                              substr(x, 2, nchar(x)), sep = ""))) +
    labs(
      title = expression("Distribution of "*n[emerge]^{"needed"}*" by scenario and band"),
      x     = expression(n[emerge]^{"needed"}*" (log scale)"),
      y     = "Replicate count"
    ) +
    .theme_q1()
  .save_png(p, out_path, width_in, height_in)
  invisible(out_path)
}

# -----------------------------------------------------------------------------
# plot_boxplots
# -----------------------------------------------------------------------------
# Notched boxplots of n_emerge_needed, faceted by scenario, coloured by
# band. Notch width is median +/- 1.58 * IQR / sqrt(n) (McGill-Tukey-Larsen
# 1978); non-overlapping notches suggest a 95%-confidence difference in
# medians. Overlaid points show every replicate (alpha 0.15) for honest
# disclosure of underlying distribution.
# -----------------------------------------------------------------------------
plot_boxplots <- function(df, out_path = "outputs/plots/boxplots.png",
                          width_in = 9, height_in = 6) {
  det <- df[!is.na(df$n_emerge_needed), , drop = FALSE]
  det$scenario_lab <- factor(.SCENARIO_LABELS[as.character(det$scenario)],
                             levels = .SCENARIO_LABELS)
  
  suppressWarnings({  # notch-wider-than-box warning is expected for small n
    p <- ggplot(det, aes(x = band, y = n_emerge_needed, fill = band)) +
      stat_boxplot(geom = "errorbar", width = 0.3, linewidth = 0.35) +
      geom_boxplot(notch = TRUE, outlier.shape = NA, alpha = 0.7,
                   linewidth = 0.35) +
      geom_jitter(width = 0.18, alpha = 0.15, size = 0.6) +
      scale_y_log10(breaks = c(2, 10, 100, 1000, 10000),
                    labels = scales::label_comma()) +
      scale_fill_manual(values = .BAND_COLORS, guide = "none") +
      facet_wrap(~ scenario_lab, nrow = 2) +
      labs(
        title = expression("Distribution of "*n[emerge]^{"needed"}*" by band"),
        x     = "Band",
        y     = expression(n[emerge]^{"needed"}*" (log scale)")
      ) +
      .theme_q1()
  })
  .save_png(p, out_path, width_in, height_in)
  invisible(out_path)
}

# -----------------------------------------------------------------------------
# plot_mean_ci
# -----------------------------------------------------------------------------
# Point-and-error-bar plot of mean n_emerge_needed by (scenario, band),
# with BCa bootstrap 95% CIs.
# -----------------------------------------------------------------------------
plot_mean_ci <- function(df, out_path = "outputs/plots/mean_ci.png",
                         width_in = 8, height_in = 5, n_boot = 1e4L) {
  if (!requireNamespace("boot", quietly = TRUE))
    stop("Package 'boot' is required for mean-CI plots.", call. = FALSE)
  
  tbl <- do.call(rbind, lapply(sort(unique(df$scenario)), function(s) {
    do.call(rbind, lapply(.BAND_ORDER, function(b) {
      x <- df$n_emerge_needed[df$scenario == s & df$band == b]
      x <- x[!is.na(x)]
      if (length(x) < 3L) return(NULL)
      ci <- bca_ci_mean(x, n_boot = n_boot)
      data.frame(scenario = s, band = b, n = length(x),
                 mean = ci[["mean"]], lo = ci[["lo"]], hi = ci[["hi"]])
    }))
  }))
  tbl$scenario_lab <- factor(.SCENARIO_LABELS[as.character(tbl$scenario)],
                             levels = .SCENARIO_LABELS)
  tbl$band <- factor(tbl$band, levels = .BAND_ORDER)
  
  p <- ggplot(tbl, aes(x = band, y = mean, color = band)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.6) +
    geom_point(size = 2.4) +
    geom_text(aes(label = scales::comma(n)),
              vjust = -1.1, hjust = 0.5,
              size  = 3, color = "grey30") +
    scale_y_log10(labels = scales::label_comma(),
                  expand = expansion(mult = c(0.05, 0.15))) +
    scale_color_manual(values = .BAND_COLORS, guide = "none") +
    facet_wrap(~ scenario_lab, nrow = 2, scales = "free_y") +
    labs(
      title = expression("Mean "*n[emerge]^{"needed"}*" with 95% BCa bootstrap CIs"),
      x     = "Band",
      y     = expression("Mean "*n[emerge]^{"needed"}*" (log scale)")
    ) +
    .theme_q1()
  .save_png(p, out_path, width_in, height_in)
  invisible(out_path)
}

# -----------------------------------------------------------------------------
# plot_detection_rate_heatmaps
# -----------------------------------------------------------------------------
# Per-scenario heatmap of detection rate across the (N_ref, ratio_V1)
# grid (Sc1), the (N_ref x mean(ratio_V1, ratio_V2)) grid (Sc2), or
# (N_ref x mean(ratio_Vi)) for Sc3/Sc4. Detection rate is the fraction
# of replicates in the cell that detected (vs NA).
#
# Generates one PNG per scenario.
# -----------------------------------------------------------------------------
plot_detection_rate_heatmaps <- function(df,
                                         out_dir = "outputs/plots",
                                         width_in = 7, height_in = 5) {
  paths <- character()
  for (s in sort(unique(df$scenario))) {
    sub <- df[df$scenario == s, , drop = FALSE]
    if (nrow(sub) == 0L) next
    
    # Effective ratio for the heatmap: average over present ratio columns.
    ratio_cols <- intersect(c("ratio_V1", "ratio_V2", "ratio_V3"),
                            colnames(sub))
    sub$ratio_eff <- rowMeans(sub[, ratio_cols, drop = FALSE], na.rm = TRUE)
    
    # agg <- aggregate(
    #   cbind(n = !is.na(sub$n_emerge_needed),
    #         detected = !is.na(sub$n_emerge_needed)) ~ N_ref + ratio_eff + band,
    #   data = sub,
    #   FUN = function(x) c(n = length(x), detected = sum(x))
    # )
    # aggregate with two-column FUN is fiddly; simpler with split-by:
    agg_list <- split(sub, list(sub$N_ref, sub$ratio_eff, sub$band),
                      drop = TRUE)
    agg <- do.call(rbind, lapply(agg_list, function(g) {
      data.frame(
        N_ref       = g$N_ref[1L],
        ratio_eff   = g$ratio_eff[1L],
        band        = g$band[1L],
        n           = nrow(g),
        detected    = sum(!is.na(g$n_emerge_needed)),
        detect_rate = sum(!is.na(g$n_emerge_needed)) / nrow(g),
        stringsAsFactors = FALSE
      )
    }))
    
    p <- ggplot(agg, aes(x = factor(round(ratio_eff, 2)),
                         y = factor(N_ref), fill = detect_rate)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradient2(low = "#b2182b", mid = "#fddbc7", high = "#1b7837",
                           midpoint = 0.5, limits = c(0, 1),
                           labels = scales::percent_format(accuracy = 1),
                           name = "Detection rate") +
      scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
      facet_wrap(~ band, scales = "free_y") +
      labs(
        title = sprintf("Detection rate across the design grid: %s",
                        .SCENARIO_LABELS[as.character(s)]),
        x     = "Mean dominance ratio",
        y     = expression(N[ref])
      ) +
      .theme_q1() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                       size = rel(0.75)))
    
    pth <- file.path(out_dir, sprintf("heatmap_sc%d.png", s))
    .save_png(p, pth, width_in, height_in)
    paths <- c(paths, pth)
  }
  invisible(paths)
}

# -----------------------------------------------------------------------------
# plot_detection_curves
# -----------------------------------------------------------------------------
# Empirical cumulative detection rate as a function of n_emerge: for each
# (scenario, band), the fraction of cell-replicates that have detected
# by each candidate n_emerge value, treating NA as "never detected".
#
# This is a one-graph operationalisation of "how many sequences do we
# need to see this variant?" — exactly the study's primary question.
# -----------------------------------------------------------------------------
plot_detection_curves <- function(df,
                                  out_path = "outputs/plots/detection_curves.png",
                                  width_in = 9, height_in = 5) {
  # Empirical CDF per (scenario, band), treating NA as +Inf
  df$ne <- ifelse(is.na(df$n_emerge_needed), Inf, df$n_emerge_needed)
  ecdf_df <- do.call(rbind, lapply(sort(unique(df$scenario)), function(s) {
    do.call(rbind, lapply(.BAND_ORDER, function(b) {
      x <- df$ne[df$scenario == s & df$band == b]
      if (length(x) == 0L) return(NULL)
      x_sorted <- sort(x)
      data.frame(
        scenario = s,
        band     = b,
        n_emerge = x_sorted,
        cum_rate = seq_along(x_sorted) / length(x_sorted),
        stringsAsFactors = FALSE
      )
    }))
  }))
  ecdf_df <- ecdf_df[is.finite(ecdf_df$n_emerge), , drop = FALSE]
  ecdf_df$scenario_lab <- factor(.SCENARIO_LABELS[as.character(ecdf_df$scenario)],
                                 levels = .SCENARIO_LABELS)
  ecdf_df$band <- factor(ecdf_df$band, levels = .BAND_ORDER)
  
  p <- ggplot(ecdf_df, aes(x = n_emerge, y = cum_rate, color = band)) +
    geom_step(linewidth = 0.7) +
    scale_x_log10(breaks = c(2, 10, 100, 1000, 10000),
                  labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    scale_color_manual(values = .BAND_COLORS, name = "Band") +
    facet_wrap(~ scenario_lab, nrow = 2) +
    labs(
      title = expression("Cumulative detection rate as a function of "*n[emerge]),
      x     = expression(n[emerge]*" (log scale)"),
      y     = "Cumulative fraction of replicates detected"
    ) +
    .theme_q1()
  .save_png(p, out_path, width_in, height_in)
  invisible(out_path)
}

# -----------------------------------------------------------------------------
# generate_all_plots
# -----------------------------------------------------------------------------
# Convenience wrapper: runs every plot generator above, plus
# build_results_table(), and writes everything to outputs/plots/ and
# outputs/tables/. Returns a named list of file paths.
# -----------------------------------------------------------------------------
generate_all_plots <- function(config = NULL, n_boot = 1e4L) {
  if (is.null(config)) {
    config <- build_config()
    config$STUDY_DIR <- resolve_study_dir()
  }
  df <- load_all_summaries(config)
  
  plots_dir  <- file.path(output_dir(config), "plots")
  tables_dir <- file.path(output_dir(config), "tables")
  dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
  
  log_msg("Generating summary table...")
  tbl <- build_results_table(df, n_boot = n_boot)
  saveRDS(tbl, file.path(tables_dir, "results_table.rds"))
  write.csv(tbl, file.path(tables_dir, "results_table.csv"),
            row.names = FALSE)
  
  log_msg("Generating histograms...")
  hist_path <- plot_histograms(df,
                               out_path = file.path(plots_dir,
                                                    "histograms.png"))
  
  log_msg("Generating boxplots...")
  box_path <- plot_boxplots(df,
                            out_path = file.path(plots_dir,
                                                 "boxplots.png"))
  
  log_msg("Generating mean-CI plot...")
  mean_path <- plot_mean_ci(df,
                            out_path = file.path(plots_dir,
                                                 "mean_ci.png"),
                            n_boot = n_boot)
  
  log_msg("Generating Kaplan-Meier curves...")
  km_path <- plot_km_curves(df,
                            out_path = file.path(plots_dir, "km_curves.png"))
  
  log_msg("Generating K-M summary table...")
  km_tbl <- km_summary_table(df)
  saveRDS(km_tbl, file.path(tables_dir, "km_summary.rds"))
  write.csv(km_tbl, file.path(tables_dir, "km_summary.csv"),
            row.names = FALSE)
  
  log_msg("Generating heatmaps...")
  heat_paths <- plot_detection_rate_heatmaps(df, out_dir = plots_dir)
  
  log_msg("Generating detection curves...")
  curves_path <- plot_detection_curves(df,
                                       out_path = file.path(plots_dir,
                                                            "detection_curves.png"))
  
  paths <- list(
    table_rds       = file.path(tables_dir, "results_table.rds"),
    table_csv       = file.path(tables_dir, "results_table.csv"),
    histograms      = hist_path,
    boxplots        = box_path,
    mean_ci         = mean_path,
    heatmaps        = heat_paths,
    detection_curves = curves_path,
    km_curves = km_path,
    km_table  = file.path(tables_dir, "km_summary.csv")
  )
  
  log_msg("Done. Outputs written to:")
  for (nm in names(paths)) {
    for (p in paths[[nm]]) log_msg(sprintf("  %-18s %s", nm, p))
  }
  invisible(paths)
}

# -----------------------------------------------------------------------------
# Survival-analysis helpers
# -----------------------------------------------------------------------------
# Frame the sample-size sweep as a time-to-event problem to handle the
# right-censored NA replicates (sweep reached the cell ceiling without
# detecting). "Time" is n_emerge_needed; the event is detection; the
# censoring time is the cell's ceiling.
#
# Used to complement the bootstrap-CI summary in build_results_table()
# with Kaplan-Meier medians that are unbiased under right-censoring.
# Requires the survival package; survminer is optional for publication-
# quality plots and is auto-detected.
# -----------------------------------------------------------------------------

# Adds two columns: time_to_detect and status (0 = censored, 1 = event).
prepare_survival_data <- function(df) {
  df$time_to_detect <- ifelse(is.na(df$n_emerge_needed),
                              df$ceiling, df$n_emerge_needed)
  df$status         <- ifelse(is.na(df$n_emerge_needed), 0L, 1L)
  df
}

# Per-group Kaplan-Meier fit. Returns survfit object.
fit_km_per_group <- function(df, group = c("band", "scenario", "both")) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Package 'survival' is required for K-M analysis.", call. = FALSE)
  group <- match.arg(group)
  sd <- prepare_survival_data(df)
  fml <- switch(group,
                band     = survival::Surv(time_to_detect, status) ~ band,
                scenario = survival::Surv(time_to_detect, status) ~ scenario,
                both     = survival::Surv(time_to_detect, status) ~ scenario + band
  )
  survival::survfit(fml, data = sd)
}

# Tidy summary table of K-M medians and 95% CIs per (scenario, band) plus
# log-rank p-values for between-group comparisons within each scenario.
# Returns a data frame with one row per (scenario, band) plus a "logrank"
# attribute carrying scenario-level test results.
km_summary_table <- function(df) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Package 'survival' is required for K-M analysis.", call. = FALSE)
  sd <- prepare_survival_data(df)
  
  rows <- list()
  logrank <- list()
  for (s in sort(unique(sd$scenario))) {
    sub  <- sd[sd$scenario == s, , drop = FALSE]
    if (nrow(sub) == 0L) next
    fit  <- survival::survfit(survival::Surv(time_to_detect, status) ~ band,
                              data = sub)
    fsum <- summary(fit)$table
    if (is.matrix(fsum)) {
      # rows are "band=small", "band=medium", "band=large"
      band_keys <- sub("^band=", "", rownames(fsum))
      for (b in .BAND_ORDER) {
        idx <- match(b, band_keys)
        if (is.na(idx)) next
        rows[[length(rows) + 1L]] <- data.frame(
          scenario   = s, band = b,
          n          = unname(fsum[idx, "records"]),
          events     = unname(fsum[idx, "events"]),
          censored   = unname(fsum[idx, "records"] - fsum[idx, "events"]),
          km_median  = unname(fsum[idx, "median"]),
          km_ci_low  = unname(fsum[idx, "0.95LCL"]),
          km_ci_high = unname(fsum[idx, "0.95UCL"]),
          stringsAsFactors = FALSE, check.names = FALSE
        )
      }
    }
    # Log-rank test across bands within this scenario
    if (length(unique(sub$band)) > 1L) {
      lr <- survival::survdiff(survival::Surv(time_to_detect, status) ~ band,
                               data = sub)
      logrank[[as.character(s)]] <- list(
        chisq = lr$chisq,
        df    = length(lr$n) - 1L,
        p     = pchisq(lr$chisq, df = length(lr$n) - 1L, lower.tail = FALSE)
      )
    }
  }
  out <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()
  attr(out, "logrank") <- logrank
  rownames(out) <- NULL
  out
}

# Publication-quality K-M curves per scenario, with CI bands and the
# cumulative-event (detection) framing rather than the survival framing
# (so the y-axis ascends with sequence depth, matching biologists'
# reading direction).
plot_km_curves <- function(df, out_path = "outputs/plots/km_curves.png",
                           width_in = 9, height_in = 6) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("Package 'survival' is required.", call. = FALSE)
  
  sd <- prepare_survival_data(df)
  
  # Build tidy data per (scenario, band) from survfit objects manually
  # so we have one consistent code path with or without survminer.
  ribbons <- list()
  for (s in sort(unique(sd$scenario))) {
    sub <- sd[sd$scenario == s, , drop = FALSE]
    fit <- survival::survfit(survival::Surv(time_to_detect, status) ~ band,
                             data = sub)
    strata_levels <- if (is.null(fit$strata)) "(all)"
    else rep(sub("^band=", "", names(fit$strata)),
             fit$strata)
    ribbons[[length(ribbons) + 1L]] <- data.frame(
      scenario = s,
      band     = factor(strata_levels, levels = .BAND_ORDER),
      time     = fit$time,
      cum_det  = 1 - fit$surv,
      cum_lo   = 1 - fit$upper,    # K-M upper bound of S -> lower of CDF
      cum_hi   = 1 - fit$lower,
      stringsAsFactors = FALSE
    )
  }
  km_df <- do.call(rbind, ribbons)
  km_df$scenario_lab <- factor(.SCENARIO_LABELS[as.character(km_df$scenario)],
                               levels = .SCENARIO_LABELS)
  
  p <- ggplot(km_df, aes(x = time, y = cum_det, color = band, fill = band)) +
    geom_ribbon(aes(ymin = cum_lo, ymax = cum_hi),
                alpha = 0.18, color = NA) +
    geom_step(linewidth = 0.7) +
    scale_x_log10(breaks = c(2, 10, 100, 1000, 10000),
                  labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, 1)) +
    scale_color_manual(values = .BAND_COLORS, name  = "Band") +
    scale_fill_manual( values = .BAND_COLORS, guide = "none") +
    facet_wrap(~ scenario_lab, nrow = 2) +
    labs(
      title = expression("Kaplan-Meier detection curves with 95% CI"),
      x     = expression(n[emerge]*" (log scale)"),
      y     = "Cumulative probability of detection"
    ) +
    .theme_q1()
  .save_png(p, out_path, width_in, height_in)
  invisible(out_path)
}