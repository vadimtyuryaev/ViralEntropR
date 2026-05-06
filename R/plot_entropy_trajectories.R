#' @title Plot Shannon Entropy Trajectories
#'
#' @description Plots per-site Shannon entropy as continuous trajectories across
#'   time partitions for a selected set of sequence sites, using the output of
#'   \code{\link{partition_time_windows}}.
#'
#' @details
#' For each partition the function extracts the GMM clustering result from
#' \code{part_data$Clusters} and assembles a long-format data frame spanning
#' all selected sites across all partitions.  Sites absent from a given
#' partition (removed by zero-entropy or singleton filtering, or because the
#' partition window was empty) are silently omitted from that partition's
#' trajectory and do not interrupt adjacent observations.
#'
#' \strong{Class relabeling.}  This function does not perform any relabeling of
#' GMM class labels.  If class 1 must denote the highest-entropy group
#' throughout the returned \code{$Data_Frame} (e.g. before passing it to
#' \code{\link{plot_site_class_trajectory}}), the user should call
#' \code{\link{relabel_entropy_classes}} on each partition's
#' \code{Clusters[[i]]$DataFrame} and update \code{Max_Entropy[i]} to
#' \code{1L} prior to calling this function.
#'
#' \strong{Colour scheme.}  Site colours are specified through
#' \code{site_colors}, a named character vector whose names are site indices
#' (as character strings) and whose values are valid R colour strings.  Any
#' site not listed in \code{site_colors} receives an automatically assigned
#' colour from the HCL \code{"Dark 2"} qualitative palette.  The final colour
#' mapping is returned as \code{$Colors} so the same scheme can be passed to
#' subsequent calls for cross-plot consistency.
#'
#' \strong{Group-stratified trajectories (\code{by_group = TRUE}).}  When
#' biological groupings must be distinguished visually (e.g. defining SNP
#' sites vs. other mutation sites), \code{groups_list} partitions
#' \code{sites} into explicitly named groups.  Any site not assigned to an
#' explicit group is automatically collected into a remainder group appended
#' as the final element of \code{groups_list}.  Line type and line width are
#' mapped to group membership via \code{line_type_groups} and
#' \code{line_size_groups}, both of which must have length equal to the
#' total number of groups (explicit plus the automatic remainder).  At most
#' six groups are supported.
#'
#' \strong{\code{max_class} column.}  The returned \code{$Data_Frame} carries
#' a \code{max_class} column recording the label of the highest-entropy GMM
#' component for each partition, taken directly from
#' \code{part_data$Max_Entropy[i]}.  This column is consumed by
#' \code{\link{plot_site_class_trajectory}} (red labels) and by downstream
#' class-assignment tables.
#'
#' @param part_data Named list.  Output of \code{\link{partition_time_windows}},
#'   optionally with per-partition \code{Clusters[[i]]$DataFrame} already
#'   relabeled by the user via \code{\link{relabel_entropy_classes}}.  Must
#'   contain elements \code{Clusters}, \code{Max_Entropy}, \code{Dates_Labels},
#'   and \code{N_partitions}.
#' @param sites Integer vector.  Site indices to include.  Defaults to the
#'   union of all sites observed across all partitions (i.e. every site that
#'   has non-zero, non-singleton entropy in at least one partition window).
#' @param labels Character vector of length \code{N_partitions}.  Partition
#'   labels used on the x-axis.  Defaults to \code{part_data$Dates_Labels}.
#' @param site_colors Named character vector.  Names are site indices as
#'   character strings (e.g. \code{"681"}); values are colour strings
#'   (e.g. \code{"#FB8072"}).  Sites absent from \code{site_colors} receive
#'   automatically assigned colours.  Default is \code{NULL} (all colours
#'   auto-assigned).
#' @param by_group Logical.  If \code{TRUE}, maps line type and line width to
#'   site groups defined by \code{groups_list}.  Default is \code{FALSE}.
#' @param groups_list List of integer vectors.  Each element specifies the
#'   site indices belonging to one explicit group.  Sites in \code{sites} not
#'   covered by any explicit group are automatically assigned to a remainder
#'   group appended as the final element.  Total group count (explicit plus
#'   remainder) must not exceed 6.  Required when \code{by_group = TRUE}.
#' @param line_type_groups Character vector.  One line-type string per group
#'   (in order, including the automatic remainder group).  Must have length
#'   equal to the total number of groups.  Defaults to \code{"solid"} for the
#'   first group and \code{"dashed"} for all remaining groups.
#' @param line_size_groups Numeric vector.  One line-width value per group
#'   (in order, including the automatic remainder group).  Must have length
#'   equal to the total number of groups.  Defaults to \code{2} for the first
#'   group and \code{1} for all remaining groups.
#' @param transformation Object of class \code{"transform"} or \code{"trans"}
#'   as returned by \code{\link[scales]{trans_new}}, or \code{NULL} (identity,
#'   no transformation).  Applied to the y-axis via
#'   \code{\link[ggplot2]{scale_y_continuous}}.  Default is \code{NULL}.
#' @param line_size Numeric.  Line width used when \code{by_group = FALSE}.
#'   Default is \code{1.5}.
#' @param legend Logical.  If \code{TRUE} (default), the site colour legend is
#'   displayed.
#' @param legend_text_size Numeric.  Font size of legend text in points.
#'   Default is \code{12}.
#' @param x_angle Numeric.  Rotation angle of x-axis tick labels in degrees.
#'   Default is \code{45}.
#' @param grayscale Logical.  If \code{TRUE}, overrides \code{site_colors} and
#'   renders all trajectories in greyscale.  Default is \code{FALSE}.
#' @param plot_title Character.  Plot title string.  Default is
#'   \code{"Shannon Entropy Trajectories"}.
#'
#' @return A named list with five elements:
#' \item{Data_Frame}{Long-format data frame with columns \code{sites}
#'   (factor), \code{entropies} (numeric), \code{class} (factor),
#'   \code{max_class} (integer), \code{period} (integer), and
#'   \code{coverage} (character, partition label).  Suitable for direct
#'   input to \code{\link{plot_site_class_trajectory}}.}
#' \item{Plot}{A \code{ggplot} object.  Augment with additional layers
#'   (e.g. \code{\link[ggplot2]{geom_vline}} for VOC emergence events) before
#'   printing or saving with \code{\link[ggplot2]{ggsave}}.}
#' \item{Colors}{Named character vector mapping each plotted site index
#'   (character) to its assigned colour string.  Pass as \code{site_colors}
#'   to subsequent calls to \code{plot_entropy_trajectories} for a consistent
#'   colour scheme across figures.}
#' \item{XBreaks}{Integer vector of partition period indices.  Pass as
#'   \code{xbreaks} to \code{\link{plot_site_class_trajectory}}.}
#' \item{XLabels}{Character vector of partition labels aligned with
#'   \code{XBreaks}.  Pass as \code{xlabels} to
#'   \code{\link{plot_site_class_trajectory}}.}
#'
#' @seealso
#' \code{\link{partition_time_windows}},
#' \code{\link{relabel_entropy_classes}},
#' \code{\link{plot_site_class_trajectory}}
#'
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous
#'   scale_y_continuous scale_color_manual scale_color_grey
#'   scale_linetype_manual scale_linewidth_manual labs theme_bw theme
#'   element_text unit guides guide_legend
#' @importFrom grDevices hcl.colors
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # Synthetic dataset: site 2 accumulates variability across partitions.
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
#' # Example 1: no relabeling — class labels are as assigned by Mclust.
#' result <- plot_entropy_trajectories(
#'   part_data  = part_data
#' )
#' print(result$Plot)

plot_entropy_trajectories <- function(part_data,
                                      sites            = NULL,
                                      labels           = NULL,
                                      site_colors      = NULL,
                                      by_group         = FALSE,
                                      groups_list      = NULL,
                                      line_type_groups = NULL,
                                      line_size_groups = NULL,
                                      transformation   = NULL,
                                      line_size        = 1.5,
                                      legend           = TRUE,
                                      legend_text_size = 12,
                                      x_angle          = 45,
                                      grayscale        = FALSE,
                                      plot_title       = "Shannon Entropy Trajectories") {

  # ---------------------------------------------------------------------------
  # 1. Validate part_data
  # ---------------------------------------------------------------------------
  required_fields <- c("Clusters", "Max_Entropy", "Dates_Labels", "N_partitions")
  missing_fields  <- setdiff(required_fields, names(part_data))
  if (length(missing_fields) > 0L)
    stop("`part_data` is missing required fields: ",
         paste(missing_fields, collapse = ", "),
         ". Ensure `part_data` is the direct output of partition_time_windows().", 
         call. = FALSE)

  n_part <- part_data$N_partitions

  if (n_part < 1L)
    stop("`part_data$N_partitions` is 0. There are no partitions to plot.", 
         call. = FALSE)

  # ---------------------------------------------------------------------------
  # 2. Partition labels
  # ---------------------------------------------------------------------------
  if (is.null(labels)) {
    labels <- part_data$Dates_Labels
  } else {
    labels <- as.character(labels)
    if (length(labels) != n_part)
      stop("`labels` must have length equal to `N_partitions` (", n_part, ").", 
           call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 3. Validate transformation
  # ---------------------------------------------------------------------------
  if (!is.null(transformation)) {
    if (!inherits(transformation, c("transform", "trans")))
      stop("`transformation` must be an object of class 'transform' produced by ",
           "scales::trans_new(), or NULL for the identity (no transformation).", 
           call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 4. Validate by_group arguments before the main loop
  # ---------------------------------------------------------------------------
  if (isTRUE(by_group)) {
    if (is.null(groups_list) || !is.list(groups_list) || length(groups_list) == 0L)
      stop("`groups_list` must be a non-empty list of integer vectors ",
           "when `by_group = TRUE`.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # 5. Assemble long-format data frame across all partitions
  # ---------------------------------------------------------------------------
  frames <- vector("list", n_part)

  for (i in seq_len(n_part)) {

    cl_i <- part_data$Clusters[[i]]
    df_i <- cl_i$DataFrame

    # Empty partition — skip; its absence creates a gap in the trajectory.
    if (is.null(df_i) || nrow(df_i) == 0L) next

    # max_class is taken directly from the pre-computed Max_Entropy field.
    # If the user pre-relabeled the partition DataFrames, they are responsible
    # for also updating Max_Entropy[i] to 1L before calling this function.
    me_val      <- part_data$Max_Entropy[i]
    max_class_i <- if (is.na(me_val)) NA_integer_ else as.integer(me_val)

    df_i$max_class <- max_class_i
    df_i$period    <- i
    df_i$coverage  <- labels[i]

    frames[[i]] <- df_i
  }

  frames <- Filter(Negate(is.null), frames)

  if (length(frames) == 0L)
    stop("No data are available for plotting: all partitions are empty.", 
         call. = FALSE)

  df_all           <- do.call(rbind, frames)
  rownames(df_all) <- NULL

  # ---------------------------------------------------------------------------
  # 6. Subset to requested sites
  # ---------------------------------------------------------------------------
  all_present <- sort(unique(as.integer(df_all$sites)))

  if (is.null(sites)) {
    sites <- all_present
  } else {
    sites  <- sort(as.integer(sites))
    absent <- setdiff(sites, all_present)
    if (length(absent) > 0L)
      warning(length(absent), " requested site(s) not found in any partition ",
              "(zero entropy, singleton removal, or empty window): ",
              paste(absent, collapse = ", "), ". These sites are omitted from the plot.")
    sites <- intersect(sites, all_present)
  }

  if (length(sites) == 0L)
    stop("No requested sites are present in any partition.", call. = FALSE)

  df_all <- df_all[df_all$sites %in% sites, , drop = FALSE]

  # Coerce for ggplot2: sites as ordered factor (preserves numeric sort order
  # in the legend), class and period as appropriate types.
  df_all$sites  <- factor(df_all$sites, levels = as.character(sort(sites)))
  df_all$class  <- as.factor(df_all$class)
  df_all$period <- as.integer(df_all$period)

  existing_sites <- levels(df_all$sites)     # character, sorted numerically
  n_sites_exist  <- length(existing_sites)

  # ---------------------------------------------------------------------------
  # 7. Assign site colours
  # ---------------------------------------------------------------------------
  # Auto-generate a qualitative palette; override any user-specified entries.
  auto_palette        <- grDevices::hcl.colors(n_sites_exist, palette = "Dark 2")
  names(auto_palette) <- existing_sites

  if (!is.null(site_colors)) {
    if (!is.character(site_colors) || is.null(names(site_colors)))
      stop("`site_colors` must be a named character vector with site indices as ",
           "names (e.g. c(\"681\" = \"#FB8072\", \"501\" = \"#80B1D3\")).", 
           call. = FALSE)
    sc_keys <- as.character(names(site_colors))
    matches <- intersect(sc_keys, existing_sites)
    if (length(matches) > 0L)
      auto_palette[matches] <- site_colors[sc_keys %in% matches]
  }

  final_palette <- auto_palette    # named character: site_char -> colour

  # ---------------------------------------------------------------------------
  # 8. x-axis breaks and labels
  # ---------------------------------------------------------------------------
  m_min    <- min(df_all$period)
  m_max    <- max(df_all$period)
  x_breaks <- seq.int(m_min, m_max)
  x_labels <- labels[x_breaks]

  # ---------------------------------------------------------------------------
  # 9. by_group: assign group membership and resolve line aesthetics
  # ---------------------------------------------------------------------------
  if (isTRUE(by_group)) {

    # Coerce group site indices to character for factor-level matching.
    groups_chr <- lapply(groups_list, function(g) as.character(as.integer(g)))

    # Auto-create remainder group for sites not covered by any explicit group.
    assigned  <- unlist(groups_chr, use.names = FALSE)
    remainder <- setdiff(existing_sites, assigned)
    if (length(remainder) > 0L)
      groups_chr[[length(groups_chr) + 1L]] <- remainder

    n_groups <- length(groups_chr)

    if (n_groups > 6L)
      stop("Total number of groups (explicit + automatic remainder) is ",
           n_groups, ". Maximum allowed is 6.", call. = FALSE)

    # Default line types: solid for group 1, dashed for all others.
    if (is.null(line_type_groups)) {
      line_type_groups <- c("solid", rep("dashed", n_groups - 1L))
    } else {
      if (length(line_type_groups) != n_groups)
        stop("`line_type_groups` must have length ", n_groups,
             " (number of explicit groups plus the automatic remainder group).", 
             call. = FALSE)
    }

    # Default line widths: 2 for group 1, 1 for all others.
    if (is.null(line_size_groups)) {
      line_size_groups <- c(2, rep(1, n_groups - 1L))
    } else {
      if (length(line_size_groups) != n_groups)
        stop("`line_size_groups` must have length ", n_groups,
             " (number of explicit groups plus the automatic remainder group).", 
             call. = FALSE)
    }

    # Assign group integer to each row, matched on the sites factor level.
    df_all$group <- NA_integer_
    for (g in seq_len(n_groups)) {
      mask               <- as.character(df_all$sites) %in% groups_chr[[g]]
      df_all$group[mask] <- g
    }
    df_all$group <- as.factor(df_all$group)

  }

  # ---------------------------------------------------------------------------
  # 10. y-axis label
  # ---------------------------------------------------------------------------
  y_label <- if (is.null(transformation)) {
    "Shannon Entropy"
  } else {
    paste0(transformation$name, "(Shannon Entropy)")
  }

  # ---------------------------------------------------------------------------
  # 11. Build ggplot
  # ---------------------------------------------------------------------------

  # --- 11a. Base plot and geom_line layer ------------------------------------
  # The .data pronoun from rlang is used throughout aes() to make column
  # references unambiguous to R CMD check's static analyser.
  if (isTRUE(by_group)) {
    p <- ggplot2::ggplot(
      df_all,
      ggplot2::aes(x         = .data$period,
                   y         = .data$entropies,
                   colour    = .data$sites,
                   linetype  = .data$group,
                   linewidth = .data$group)
    ) +
      ggplot2::geom_line() +
      ggplot2::scale_linetype_manual(values = line_type_groups,
                                     name   = "Group") +
      ggplot2::scale_linewidth_manual(values = line_size_groups,
                                      name   = "Group")
  } else {
    p <- ggplot2::ggplot(
      df_all,
      ggplot2::aes(x      = .data$period,
                   y      = .data$entropies,
                   colour = .data$sites)
    ) +
      ggplot2::geom_line(linewidth = line_size)
  }

  # --- 11b. Colour scale: grayscale overrides site_colors -------------------
  if (isTRUE(grayscale)) {
    p <- p + ggplot2::scale_color_grey(name = "Site")
  } else {
    p <- p + ggplot2::scale_color_manual(values = final_palette, name = "Site")
  }

  # --- 11c. Axes and y transformation ----------------------------------------
  p <- p + ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels)

  if (!is.null(transformation)) {
    p <- p + ggplot2::scale_y_continuous(trans = transformation,
                                         name  = y_label)
  } else {
    p <- p + ggplot2::scale_y_continuous(name = y_label)
  }

  # --- 11d. Labels, theme, guides --------------------------------------------
  p <- p +
    ggplot2::labs(title = plot_title,
                  x     = "Partition") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5),
      axis.text.x      = ggplot2::element_text(angle = x_angle, hjust = 1),
      legend.text      = ggplot2::element_text(size  = legend_text_size),
      legend.key.width = ggplot2::unit(2, "cm")
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(override.aes = list(linewidth = 3))
    )

  # Legend suppression applied last to override any theme defaults.
  if (!isTRUE(legend))
    p <- p + ggplot2::theme(legend.position = "none")

  # ---------------------------------------------------------------------------
  # 12. Return
  # ---------------------------------------------------------------------------
  list(
    Data_Frame = df_all,
    Plot       = p,
    Colors     = final_palette,
    XBreaks    = x_breaks,
    XLabels    = x_labels
  )
}
