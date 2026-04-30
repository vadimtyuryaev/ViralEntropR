# =============================================================================
# GMM Site Selection — Helper Functions
# =============================================================================
#
#   plot_entropy_distributions()
#     Renders a 2×3 panel histogram of the complete entropy vector and five
#     trimmed subsets.  No GMM is fitted.  Used for the complete (unfiltered)
#     dataset visualisation only.
#
#   gmm_site_selection()
#     Fits GMM on zero/singleton-removed entropy vectors (original and/or
#     transformed scale) and renders input and selected-site histograms.
#
# Scientific note on cluster extraction
# --------------------------------------
# Mclust's classification labels (1 to G) are NOT ordered by component mean.
# The highest-entropy cluster is identified via which.max(parameters$mean),
# not max(classification).  Using max(classification) would silently select
# the wrong cluster whenever Mclust assigns a lower label index to the
# highest-mean component.
# =============================================================================


# ---------------------------------------------------------------------------
# fit_gmm_quiet — suppress Mclust progress bar output
# ---------------------------------------------------------------------------
# Mclust writes its "fitting ... |====| 100%" bar to stdout via cat().
# capture.output() with type = "output" intercepts and discards this without
# affecting warnings or errors.
fit_gmm_quiet <- function(x, ...) {
  fit <- NULL
  invisible(utils::capture.output(
    fit <- mclust::Mclust(x, ...),
    type = "output"
  ))
  fit
}


# ---------------------------------------------------------------------------
# extract_max_cl_sites — correct highest-entropy cluster extraction
# ---------------------------------------------------------------------------
extract_max_cl_sites <- function(gmm_fit) {
  max_cl      <- max(gmm_fit$classification)
  whichmax_cl <- which.max(gmm_fit$parameters$mean)
  if (max_cl != whichmax_cl)
    warning(sprintf(
      "Mclust cluster label mismatch: max(classification) = %d but which.max(parameters$mean) = %d. Components may not be ordered by mean. Verify GMM fit.",
      max_cl, whichmax_cl
    ))
  members <- gmm_fit$classification[gmm_fit$classification == max_cl]
  as.integer(names(members))
}

# ---------------------------------------------------------------------------
# plot_entropy_distributions
# ---------------------------------------------------------------------------
#' @title Plot Entropy Distribution Panels (No GMM)
#'
#' @description
#' Renders a 2×3 patchwork panel of histograms for the complete entropy vector
#' and five trimmed subsets.  No GMM is fitted.  Intended for visualising the
#' full dataset distribution prior to zero/singleton removal.
#'
#' @param entrop_vec Named numeric vector of per-site Shannon entropies.
#' @param prob Numeric in (0, 1).  Upper quantile threshold.  Default 0.9.
#' @param transfr A \code{scales} transformation object or \code{NULL}.
#' @param dataset_label Character.  Appended to the panel title.
#' @param ... Additional arguments forwarded to \code{geom_histogram}.

plot_entropy_distributions <- function(entrop_vec,
                                       prob          = 0.9,
                                       transfr       = NULL,
                                       dataset_label = "",
                                       ...) {

  stopifnot(is.numeric(entrop_vec), length(entrop_vec) > 0L,
            is.numeric(prob), prob > 0, prob < 1)

  ss          <- as.numeric(summary(entrop_vec))
  quant_label <- sprintf("Trimmed: %.0f%%+ percentile", prob * 100)

  ev_list <- vector("list", 6L)
  ev_list[[1L]] <- entrop_vec
  ev_list[[2L]] <- entrop_vec[entrop_vec >= ss[2L]]
  ev_list[[3L]] <- entrop_vec[entrop_vec >= ss[3L]]
  ev_list[[4L]] <- entrop_vec[entrop_vec >= ss[4L]]
  ev_list[[5L]] <- entrop_vec[entrop_vec >= ss[5L]]
  ev_list[[6L]] <- entrop_vec[entrop_vec >= quantile(entrop_vec, prob)]
  names(ev_list) <- c("Complete data",
                      "Trimmed: Q1+ (top 75%)",
                      "Trimmed: median+ (top 50%)",
                      "Trimmed: mean+",
                      "Trimmed: Q3+ (top 25%)",
                      quant_label)

  if (!is.null(transfr)) {
    ev_list  <- lapply(ev_list, transfr$transform)
    x_label  <- sprintf("%s-transformed Shannon entropy", transfr$name)
    main_ttl <- sprintf("Entropy distributions (%s) \u2014 %s",
                        transfr$name, dataset_label)
  } else {
    x_label  <- "Shannon entropy"
    main_ttl <- if (nzchar(dataset_label))
      sprintf("Entropy distributions \u2014 %s", dataset_label)
    else "Entropy distributions"
  }

  theme_hist <- ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 9,
                                               hjust = 0.5),
      axis.title       = ggplot2::element_text(size = 8),
      axis.text        = ggplot2::element_text(size = 7),
      panel.grid.minor = ggplot2::element_blank()
    )

  panels <- lapply(names(ev_list), function(nm) {
    ggplot2::ggplot(data.frame(x = ev_list[[nm]]), ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(fill = "steelblue", colour = "white",
                              alpha = 0.85, ...) +
      ggplot2::labs(title = nm, x = x_label, y = "Count") +
      theme_hist
  })

  print(
    patchwork::wrap_plots(panels, ncol = 3L) +
      patchwork::plot_annotation(
        title = main_ttl,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", hjust = 0.5,
                                             size = 11)))
  )
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# gmm_site_selection
# ---------------------------------------------------------------------------
#' @title GMM Site Selection with Optional Entropy Transformation
#'
#' @description
#' Fits GMM via \code{mclust::Mclust} on an entropy vector (original and/or
#' transformed scale) and renders 2×3 panel histograms of input distributions
#' and selected-site distributions.  Call only on zero/singleton-removed
#' entropy vectors; use \code{plot_entropy_distributions()} for complete data.
#'
#' @param entrop_vec Named numeric vector of per-site Shannon entropies
#'   (zeros and singletons already removed).
#' @param prob Numeric in (0, 1).  Upper quantile threshold.  Default 0.9.
#' @param transfr A \code{scales} transformation object or \code{NULL}.
#' @param hist Logical.  Render histograms.  Default \code{TRUE}.
#' @param calc_gmm_orig Logical.  Fit GMM on the original scale.
#' @param calc_gmm_transfr Logical.  Fit GMM on the transformed scale.
#'   Ignored when \code{transfr = NULL}.
#' @param dataset_label Character.  Appended to histogram titles.
#' @param ... Additional arguments forwarded to \code{geom_histogram}.
#'
#' @return Named list with elements \code{Objects} and \code{Sites}.

gmm_site_selection <- function(entrop_vec,
                               prob             = 0.9,
                               transfr          = NULL,
                               hist             = TRUE,
                               calc_gmm_orig    = TRUE,
                               calc_gmm_transfr = TRUE,
                               dataset_label    = "",
                               ...) {

  stopifnot(is.numeric(entrop_vec), length(entrop_vec) > 0L,
            is.numeric(prob), prob > 0, prob < 1)

  # ── Trimmed sub-vectors ───────────────────────────────────────────────────
  ss <- as.numeric(summary(entrop_vec))

  ev_complete <- entrop_vec
  ev_q1plus   <- entrop_vec[entrop_vec >= ss[2L]]
  ev_median   <- entrop_vec[entrop_vec >= ss[3L]]
  ev_mean     <- entrop_vec[entrop_vec >= ss[4L]]
  ev_q3plus   <- entrop_vec[entrop_vec >= ss[5L]]
  ev_quant    <- entrop_vec[entrop_vec >= quantile(entrop_vec, prob)]

  panel_labels <- c(
    "Complete data",
    "Trimmed: Q1+ (top 75%)",
    "Trimmed: median+ (top 50%)",
    "Trimmed: mean+",
    "Trimmed: Q3+ (top 25%)",
    sprintf("Trimmed: %.0f%%+ percentile", prob * 100)
  )

  ttl_sfx <- if (nzchar(dataset_label))
    sprintf(" \u2014 %s", dataset_label) else ""

  # ── Shared themes ─────────────────────────────────────────────────────────
  theme_hist <- ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = 9,
                                               hjust = 0.5),
      axis.title       = ggplot2::element_text(size = 8),
      axis.text        = ggplot2::element_text(size = 7),
      panel.grid.minor = ggplot2::element_blank()
    )

  make_hist <- function(data, title, x_label = "Shannon entropy") {
    ggplot2::ggplot(data.frame(x = data), ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(fill = "steelblue", colour = "white",
                              alpha = 0.85, ...) +
      ggplot2::labs(title = title, x = x_label, y = "Count") +
      theme_hist
  }

  make_hist_sel <- function(data, title,
                            x_label = "Shannon entropy (original scale)") {
    ggplot2::ggplot(data.frame(x = data), ggplot2::aes(x = x)) +
      ggplot2::geom_histogram(fill = "#2ca25f", colour = "white",
                              alpha = 0.85, ...) +
      ggplot2::labs(title = title, x = x_label, y = "Count") +
      theme_hist
  }

  print_panel_pair <- function(in_data_list, sel_data_list,
                               in_labels, sel_labels,
                               in_title, sel_title,
                               x_label_in  = "Shannon entropy",
                               x_label_sel = "Shannon entropy (original scale)") {
    panels_in <- mapply(
      function(d, lbl) make_hist(d, lbl, x_label = x_label_in),
      d = in_data_list, lbl = in_labels, SIMPLIFY = FALSE)
    panels_sl <- mapply(
      function(d, lbl) make_hist_sel(d, lbl, x_label = x_label_sel),
      d = sel_data_list, lbl = sel_labels, SIMPLIFY = FALSE)
    print(patchwork::wrap_plots(panels_in, ncol = 3L) +
            patchwork::plot_annotation(
              title = in_title,
              theme = ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                   size = 11))))
    print(patchwork::wrap_plots(panels_sl, ncol = 3L) +
            patchwork::plot_annotation(
              title = sel_title,
              theme = ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                   size = 11))))
  }

  # ── GMM on original scale ─────────────────────────────────────────────────
  GMM_orig       <- NULL
  GMM_orig_sites <- NULL

  if (calc_gmm_orig) {
    gmm_complete <- fit_gmm_quiet(ev_complete)
    gmm_q1plus   <- fit_gmm_quiet(ev_q1plus)
    gmm_median   <- fit_gmm_quiet(ev_median)
    gmm_q3plus   <- fit_gmm_quiet(ev_q3plus)
    gmm_mean     <- fit_gmm_quiet(ev_mean)
    gmm_quant    <- fit_gmm_quiet(ev_quant)

    GMM_orig <- list(GMM_orig  = gmm_complete, GMM_75   = gmm_q1plus,
                     GMM_50    = gmm_median,   GMM_25   = gmm_q3plus,
                     GMM_mean  = gmm_mean,     GMM_quant = gmm_quant)

    GMM_orig_sites <- list(
      GMM_orig_sites  = extract_max_cl_sites(gmm_complete),
      GMM_75_sites    = extract_max_cl_sites(gmm_q1plus),
      GMM_50_sites    = extract_max_cl_sites(gmm_median),
      GMM_25_sites    = extract_max_cl_sites(gmm_q3plus),
      GMM_mean_sites  = extract_max_cl_sites(gmm_mean),
      GMM_quant_sites = extract_max_cl_sites(gmm_quant)
    )

    if (hist) {
      print_panel_pair(
        in_data_list  = list(ev_complete, ev_q1plus, ev_median,
                             ev_mean, ev_q3plus, ev_quant),
        sel_data_list = list(
          ev_complete[as.character(GMM_orig_sites$GMM_orig_sites)],
          ev_complete[as.character(GMM_orig_sites$GMM_75_sites)],
          ev_complete[as.character(GMM_orig_sites$GMM_50_sites)],
          ev_complete[as.character(GMM_orig_sites$GMM_mean_sites)],
          ev_complete[as.character(GMM_orig_sites$GMM_25_sites)],
          ev_complete[as.character(GMM_orig_sites$GMM_quant_sites)]
        ),
        in_labels  = panel_labels,
        sel_labels = paste("GMM selected:", panel_labels),
        in_title   = sprintf("Entropy distributions (original scale)%s", ttl_sfx),
        sel_title  = sprintf("GMM-selected site entropies (original scale)%s",
                             ttl_sfx)
      )
    }
  }

  # ── GMM on transformed scale ──────────────────────────────────────────────
  GMM_transfr       <- NULL
  GMM_transfr_sites <- NULL

  if (calc_gmm_transfr && !is.null(transfr)) {
    ev_complete_t <- transfr$transform(ev_complete)
    ev_q1plus_t   <- transfr$transform(ev_q1plus)
    ev_median_t   <- transfr$transform(ev_median)
    ev_mean_t     <- transfr$transform(ev_mean)
    ev_q3plus_t   <- transfr$transform(ev_q3plus)
    ev_quant_t    <- transfr$transform(ev_quant)

    gmm_t_complete <- fit_gmm_quiet(ev_complete_t)
    gmm_t_q1plus   <- fit_gmm_quiet(ev_q1plus_t)
    gmm_t_median   <- fit_gmm_quiet(ev_median_t)
    gmm_t_q3plus   <- fit_gmm_quiet(ev_q3plus_t)
    gmm_t_mean     <- fit_gmm_quiet(ev_mean_t)
    gmm_t_quant    <- fit_gmm_quiet(ev_quant_t)

    GMM_transfr <- list(
      GMM_transfr_orig = gmm_t_complete, GMM_75   = gmm_t_q1plus,
      GMM_50           = gmm_t_median,   GMM_25   = gmm_t_q3plus,
      GMM_mean         = gmm_t_mean,     GMM_quant = gmm_t_quant)

    GMM_transfr_sites <- list(
      GMM_transfr_orig_sites  = extract_max_cl_sites(gmm_t_complete),
      GMM_transfr_75_sites    = extract_max_cl_sites(gmm_t_q1plus),
      GMM_transfr_50_sites    = extract_max_cl_sites(gmm_t_median),
      GMM_transfr_25_sites    = extract_max_cl_sites(gmm_t_q3plus),
      GMM_transfr_mean_sites  = extract_max_cl_sites(gmm_t_mean),
      GMM_transfr_quant_sites = extract_max_cl_sites(gmm_t_quant)
    )

    if (hist) {
      x_lab_t  <- sprintf("%s-transformed Shannon entropy", transfr$name)
      t_labels <- sprintf("%s (%s)", panel_labels, transfr$name)

      print_panel_pair(
        in_data_list  = list(ev_complete_t, ev_q1plus_t, ev_median_t,
                             ev_mean_t, ev_q3plus_t, ev_quant_t),
        sel_data_list = list(
          ev_complete[as.character(GMM_transfr_sites$GMM_transfr_orig_sites)],
          ev_complete[as.character(GMM_transfr_sites$GMM_transfr_75_sites)],
          ev_complete[as.character(GMM_transfr_sites$GMM_transfr_50_sites)],
          ev_complete[as.character(GMM_transfr_sites$GMM_transfr_mean_sites)],
          ev_complete[as.character(GMM_transfr_sites$GMM_transfr_25_sites)],
          ev_complete[as.character(GMM_transfr_sites$GMM_transfr_quant_sites)]
        ),
        in_labels   = t_labels,
        sel_labels  = sprintf("GMM selected: %s (%s)",
                              panel_labels, transfr$name),
        in_title    = sprintf("Entropy distributions (%s)%s",
                              transfr$name, ttl_sfx),
        sel_title   = sprintf(
          "GMM-selected site entropies (%s, original scale)%s",
          transfr$name, ttl_sfx),
        x_label_in  = x_lab_t
      )
    }
  }

  list(
    Objects = list(GMM_orig    = GMM_orig,
                   GMM_transfr = GMM_transfr),
    Sites   = list(GMM_orig_sites    = GMM_orig_sites,
                   GMM_transfr_sites = GMM_transfr_sites)
  )
}
