#' @title Plot GMM Entropy Class Trajectory for a Single Site
#'
#' @description Plots the Shannon entropy trajectory for a single sequence site
#'   across time partitions, with geom-label overlays at each partition
#'   recording the site's current GMM class (above the line, in green) and the
#'   partition's highest-entropy class label (below the line, in red).  The
#'   visual contrast between the two labels tracks the moment a site enters the
#'   highest-entropy cluster — the entropy-based analogue of variant emergence
#'   detection.
#'
#' @details
#' The function operates on the \code{$Data_Frame} element returned by
#' \code{\link{plot_entropy_trajectories}}, which contains columns
#' \code{sites}, \code{entropies}, \code{class}, \code{max_class},
#' \code{period}, and \code{coverage}.
#'
#' \strong{Class label interpretation.}  The green label (upper) records the
#' GMM class assigned to the site in that partition; the red label (lower)
#' records \code{max_class} — the label of the highest-entropy component for
#' that partition.  When the two labels are equal the site has entered the
#' highest-entropy class.  Without relabeling, \code{max_class} is the raw
#' Mclust label carrying the highest mean entropy (i.e.
#' \code{max(classification)} at clustering time).  When the user has
#' pre-relabeled the partitions via \code{\link{relabel_entropy_classes}}
#' before calling \code{\link{plot_entropy_trajectories}}, \code{max_class}
#' is \code{1L} throughout (class 1 is the highest-entropy group by
#' definition), and the sentinel \code{999L} is preserved unchanged.
#'
#' \strong{Axis padding.}  Both axes carry generous expansion margins so that
#' geom-label boxes at extreme entropy values or at the first and last
#' partitions do not clip the panel border.  The returned \code{ggplot} object
#' can be further adjusted by appending standard \code{ggplot2} layers with
#' \code{+} before printing or saving.
#'
#' @param data_frame Data frame.  The \code{$Data_Frame} element from
#'   \code{\link{plot_entropy_trajectories}}.  Must contain columns
#'   \code{sites}, \code{entropies}, \code{class}, \code{max_class},
#'   \code{period}, and \code{coverage}.
#' @param site Integer (length 1).  The site index to plot.  Must be present
#'   in \code{data_frame$sites}.
#' @param site_color Character.  Colour of the entropy trajectory line.
#'   Default is \code{"steelblue"}.  Pass
#'   \code{plot_entropy_trajectories()$Colors[as.character(site)]} for
#'   cross-plot colour consistency.
#' @param xbreaks Integer vector.  Partition period indices for x-axis breaks.
#'   Typically \code{plot_entropy_trajectories()$XBreaks}.
#' @param xlabels Character vector.  Partition label strings aligned with
#'   \code{xbreaks}.  Typically \code{plot_entropy_trajectories()$XLabels}.
#' @param col_current Character.  Fill colour of the upper geom-label
#'   (current GMM class at each partition).  Default is \code{"springgreen2"}.
#' @param col_max_class Character.  Fill colour of the lower geom-label
#'   (highest-entropy class label at each partition).  Default is
#'   \code{"red2"}.
#' @param label_size Numeric.  Font size of the class integer text inside the
#'   geom-labels.  Default is \code{3}.
#' @param x_angle Numeric.  Rotation angle of x-axis tick labels in degrees.
#'   Default is \code{45}.
#' @param line_size Numeric.  Width of the entropy trajectory line.  Default
#'   is \code{1.5}.
#' @param plot_title Character or \code{NULL}.  Plot title.  If \code{NULL}
#'   (default), the title is auto-generated as
#'   \code{"GMM Entropy Class Trajectory \u2014 Site <site>"}.
#' @param save Logical.  If \code{TRUE}, the plot is saved to disk via
#'   \code{\link[ggplot2]{ggsave}}.  Default is \code{FALSE}.
#' @param save_path Character or \code{NULL}.  Directory in which to save the
#'   file.  Created recursively if it does not exist.  Must be supplied when
#'   \code{save = TRUE}.  Default is \code{NULL}.
#' @param save_extension Character.  File extension including the leading dot
#'   (e.g. \code{".jpeg"}, \code{".pdf"}, \code{".png"}).  Default is
#'   \code{".jpeg"}.
#' @param width Numeric.  Saved figure width in inches.  Default is \code{20}.
#' @param height Numeric.  Saved figure height in inches.  Default is
#'   \code{15}.
#' @param dpi Numeric.  Resolution of the saved raster output in dots per
#'   inch.  Default is \code{600}.
#'
#' @return A \code{ggplot} object (returned invisibly).  Additional
#'   \code{ggplot2} layers can be appended with \code{+} before printing or
#'   saving, for example:
#'   \preformatted{
#' p <- plot_site_class_trajectory(traj$Data_Frame, site = 681L, ...)
#' p + ggplot2::geom_vline(xintercept = 14, colour = "darkorange",
#'                          linetype = "dashed", linewidth = 1.5)
#'   }
#'
#' @seealso
#' \code{\link{plot_entropy_trajectories}},
#' \code{\link{relabel_entropy_classes}},
#' \code{\link{partition_time_windows}}
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_label scale_x_continuous
#'   scale_y_continuous expansion labs theme_bw theme element_text unit ggsave
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # Shared synthetic dataset used across all three examples.
#' # Partition entropies (bits):
#' #   Jan  s1: 0.000  s2: 0.000  — both removed (zero entropy)
#' #   Feb  s1: 0.000  s2: 0.722  — s1 removed, s2 retained
#' #   Mar  s1: 0.971  s2: 1.000  — both retained
#' #   Apr  s1: 0.881  s2: 0.722  — both retained
#' df <- data.frame(
#'   s1 = c(rep(1L, 10L),
#'          rep(1L, 10L),
#'          c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
#'          c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)),
#'   s2 = c(rep(1L, 10L),
#'          c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L),
#'          c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L),
#'          c(1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)),
#'   Date = rep(
#'     seq(as.Date("2020-01-01"), by = "month", length.out = 4L),
#'     each = 10L
#'   )
#' )
#'
#' part_data <- partition_time_windows(
#'   data          = df,
#'   n_sites       = 2L,
#'   window_length = 1L,
#'   window_type   = 3L,
#'   start_date    = "2020-01-01",
#'   end_date      = "2020-04-01"
#' )
#'
#' # Example 1: without relabeling — class labels are as assigned by Mclust.
#' traj <- plot_entropy_trajectories(part_data)
#' p1 <- plot_site_class_trajectory(
#'   data_frame = traj$Data_Frame,
#'   site       = 2L,
#'   site_color = traj$Colors["2"],
#'   xbreaks    = traj$XBreaks,
#'   xlabels    = traj$XLabels
#' )
#' print(p1)

plot_site_class_trajectory <- function(data_frame,
                                       site,
                                       site_color     = "steelblue",
                                       xbreaks,
                                       xlabels,
                                       col_current    = "springgreen2",
                                       col_max_class  = "red2",
                                       label_size     = 3,
                                       x_angle        = 45,
                                       line_size      = 1.5,
                                       plot_title     = NULL,
                                       save           = FALSE,
                                       save_path      = NULL,
                                       save_extension = ".png",
                                       width          = 20,
                                       height         = 15,
                                       dpi            = 300) {

  # ---------------------------------------------------------------------------
  # 1. Input validation
  # ---------------------------------------------------------------------------
  required_cols <- c("sites", "entropies", "class", "max_class",
                     "period", "coverage")
  missing_cols  <- setdiff(required_cols, names(data_frame))
  if (length(missing_cols) > 0L)
    stop("`data_frame` is missing required columns: ",
         paste(missing_cols, collapse = ", "),
         ". Ensure `data_frame` is the `$Data_Frame` element returned by ",
         "plot_entropy_trajectories().", call. = FALSE)

  if (length(site) != 1L || !is.numeric(site))
    stop("`site` must be a single numeric site index.", call. = FALSE)
  site <- as.integer(site)

  # Subset to the requested site.
  # data_frame$sites is a factor (from plot_entropy_trajectories): convert
  # factor labels to integer for comparison with site.
  df_site <- data_frame[as.integer(as.character(data_frame$sites)) == site, ,
                        drop = FALSE]

  if (nrow(df_site) == 0L)
    stop("Site ", site, " is not present in `data_frame`. ",
         "Verify the site index and ensure it was included when calling ",
         "plot_entropy_trajectories().", call. = FALSE)

  # Coerce class and max_class to integer for unambiguous label display.
  # data_frame$class is a factor; as.character() recovers the actual label
  # string before as.integer() converts it — avoids factor-code confusion.
  df_site$class     <- as.integer(as.character(df_site$class))
  df_site$max_class <- as.integer(df_site$max_class)
  df_site$period    <- as.integer(df_site$period)

  # ---------------------------------------------------------------------------
  # 2. Validate xbreaks / xlabels
  # ---------------------------------------------------------------------------
  if (missing(xbreaks) || is.null(xbreaks))
    stop("`xbreaks` is required. Pass `plot_entropy_trajectories()$XBreaks`.", 
         call. = FALSE)
  if (missing(xlabels) || is.null(xlabels))
    stop("`xlabels` is required. Pass `plot_entropy_trajectories()$XLabels`.", 
         call. = FALSE)
  if (length(xbreaks) != length(xlabels))
    stop("`xbreaks` and `xlabels` must have equal length.", call. = FALSE)

  # ---------------------------------------------------------------------------
  # 3. Plot title
  # ---------------------------------------------------------------------------
  if (is.null(plot_title))
    plot_title <- paste0("GMM Entropy Class Trajectory \u2014 Site ", site)

  # ---------------------------------------------------------------------------
  # 4. Build ggplot
  # ---------------------------------------------------------------------------
  # The .data pronoun from rlang is used throughout aes() to make column
  # references unambiguous to R CMD check's static analyser.
  p <- ggplot2::ggplot(df_site,
                       ggplot2::aes(x = .data$period,
                                    y = .data$entropies)) +

    # Entropy trajectory line
    ggplot2::geom_line(colour    = site_color,
                       linewidth = line_size) +

    # Upper label: current GMM class at each partition (green, above line)
    ggplot2::geom_label(
      ggplot2::aes(label = .data$class),
      hjust         = 0.5,
      vjust         = -0.5,
      fill          = col_current,
      size          = label_size,
      fontface      = "bold",
      label.padding = ggplot2::unit(1, "lines")
    ) +

    # Lower label: highest-entropy class at each partition (red, below line)
    ggplot2::geom_label(
      ggplot2::aes(label = .data$max_class),
      hjust         = 0.5,
      vjust         = 0.5,
      fill          = col_max_class,
      size          = label_size,
      fontface      = "bold",
      label.padding = ggplot2::unit(1, "lines")
    ) +

    # Axis scales: generous expansion on both axes prevents geom-label boxes
    # from clipping the panel border at extreme entropy values or at the first
    # and last partitions.
    ggplot2::scale_x_continuous(
      breaks = xbreaks,
      labels = xlabels,
      expand = ggplot2::expansion(add = 0.6)
    ) +
    ggplot2::scale_y_continuous(
      name   = "Shannon Entropy",
      expand = ggplot2::expansion(mult = c(0.4, 0.4))
    ) +

    ggplot2::labs(
      title = plot_title,
      x     = "Partition"
    ) +

    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = x_angle, hjust = 1)
    )

  # ---------------------------------------------------------------------------
  # 5. Optional save
  # ---------------------------------------------------------------------------
  if (isTRUE(save)) {
    if (is.null(save_path))
      stop("`save_path` must be supplied when `save = TRUE`.", call. = FALSE)
    if (!dir.exists(save_path))
      dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    file_name <- paste0("Site_", site, "_class_trajectory", save_extension)
    ggplot2::ggsave(
      filename = file.path(save_path, file_name),
      plot     = p,
      width    = width,
      height   = height,
      dpi      = dpi
    )
  }

  # ---------------------------------------------------------------------------
  # 6. Return
  # ---------------------------------------------------------------------------
  invisible(p)
}
