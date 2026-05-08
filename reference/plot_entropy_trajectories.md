# Plot Shannon Entropy Trajectories

Plots per-site Shannon entropy as continuous trajectories across time
partitions for a selected set of sequence sites, using the output of
[`partition_time_windows`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md).

## Usage

``` r
plot_entropy_trajectories(
  part_data,
  sites = NULL,
  labels = NULL,
  site_colors = NULL,
  by_group = FALSE,
  groups_list = NULL,
  line_type_groups = NULL,
  line_size_groups = NULL,
  transformation = NULL,
  line_size = 1.5,
  legend = TRUE,
  legend_text_size = 12,
  x_angle = 45,
  grayscale = FALSE,
  plot_title = "Shannon Entropy Trajectories"
)
```

## Arguments

- part_data:

  Named list. Output of
  [`partition_time_windows`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md),
  optionally with per-partition `Clusters[[i]]$DataFrame` already
  relabeled by the user via
  [`relabel_entropy_classes`](https://vadimtyuryaev.github.io/ViralEntropR/reference/relabel_entropy_classes.md).
  Must contain elements `Clusters`, `Max_Entropy`, `Dates_Labels`, and
  `N_partitions`.

- sites:

  Integer vector. Site indices to include. Defaults to the union of all
  sites observed across all partitions (i.e. every site that has
  non-zero, non-singleton entropy in at least one partition window).

- labels:

  Character vector of length `N_partitions`. Partition labels used on
  the x-axis. Defaults to `part_data$Dates_Labels`.

- site_colors:

  Named character vector. Names are site indices as character strings
  (e.g. `"681"`); values are colour strings (e.g. `"#FB8072"`). Sites
  absent from `site_colors` receive automatically assigned colours.
  Default is `NULL` (all colours auto-assigned).

- by_group:

  Logical. If `TRUE`, maps line type and line width to site groups
  defined by `groups_list`. Default is `FALSE`.

- groups_list:

  List of integer vectors. Each element specifies the site indices
  belonging to one explicit group. Sites in `sites` not covered by any
  explicit group are automatically assigned to a remainder group
  appended as the final element. Total group count (explicit plus
  remainder) must not exceed 6. Required when `by_group = TRUE`.

- line_type_groups:

  Character vector. One line-type string per group (in order, including
  the automatic remainder group). Must have length equal to the total
  number of groups. Defaults to `"solid"` for the first group and
  `"dashed"` for all remaining groups.

- line_size_groups:

  Numeric vector. One line-width value per group (in order, including
  the automatic remainder group). Must have length equal to the total
  number of groups. Defaults to `2` for the first group and `1` for all
  remaining groups.

- transformation:

  Object of class `"transform"` or `"trans"` as returned by
  [`trans_new`](https://scales.r-lib.org/reference/new_transform.html),
  or `NULL` (identity, no transformation). Applied to the y-axis via
  [`scale_y_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html).
  Default is `NULL`.

- line_size:

  Numeric. Line width used when `by_group = FALSE`. Default is `1.5`.

- legend:

  Logical. If `TRUE` (default), the site colour legend is displayed.

- legend_text_size:

  Numeric. Font size of legend text in points. Default is `12`.

- x_angle:

  Numeric. Rotation angle of x-axis tick labels in degrees. Default is
  `45`.

- grayscale:

  Logical. If `TRUE`, overrides `site_colors` and renders all
  trajectories in greyscale. Default is `FALSE`.

- plot_title:

  Character. Plot title string. Default is
  `"Shannon Entropy Trajectories"`.

## Value

A named list with five elements:

- Data_Frame:

  Long-format data frame with columns `sites` (factor), `entropies`
  (numeric), `class` (factor), `max_class` (integer), `period`
  (integer), and `coverage` (character, partition label). Suitable for
  direct input to
  [`plot_site_class_trajectory`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_site_class_trajectory.md).

- Plot:

  A `ggplot` object. Augment with additional layers (e.g.
  [`geom_vline`](https://ggplot2.tidyverse.org/reference/geom_abline.html)
  for VOC emergence events) before printing or saving with
  [`ggsave`](https://ggplot2.tidyverse.org/reference/ggsave.html).

- Colors:

  Named character vector mapping each plotted site index (character) to
  its assigned colour string. Pass as `site_colors` to subsequent calls
  to `plot_entropy_trajectories` for a consistent colour scheme across
  figures.

- XBreaks:

  Integer vector of partition period indices. Pass as `xbreaks` to
  [`plot_site_class_trajectory`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_site_class_trajectory.md).

- XLabels:

  Character vector of partition labels aligned with `XBreaks`. Pass as
  `xlabels` to
  [`plot_site_class_trajectory`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_site_class_trajectory.md).

## Details

For each partition the function extracts the GMM clustering result from
`part_data$Clusters` and assembles a long-format data frame spanning all
selected sites across all partitions. Sites absent from a given
partition (removed by zero-entropy or singleton filtering, or because
the partition window was empty) are silently omitted from that
partition's trajectory and do not interrupt adjacent observations.

**Class relabeling.** This function does not perform any relabeling of
GMM class labels. If class 1 must denote the highest-entropy group
throughout the returned `$Data_Frame` (e.g. before passing it to
[`plot_site_class_trajectory`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_site_class_trajectory.md)),
the user should call
[`relabel_entropy_classes`](https://vadimtyuryaev.github.io/ViralEntropR/reference/relabel_entropy_classes.md)
on each partition's `Clusters[[i]]$DataFrame` and update
`Max_Entropy[i]` to `1L` prior to calling this function.

**Colour scheme.** Site colours are specified through `site_colors`, a
named character vector whose names are site indices (as character
strings) and whose values are valid R colour strings. Any site not
listed in `site_colors` receives an automatically assigned colour from
the HCL `"Dark 2"` qualitative palette. The final colour mapping is
returned as `$Colors` so the same scheme can be passed to subsequent
calls for cross-plot consistency.

**Group-stratified trajectories (`by_group = TRUE`).** When biological
groupings must be distinguished visually (e.g. defining SNP sites vs.
other mutation sites), `groups_list` partitions `sites` into explicitly
named groups. Any site not assigned to an explicit group is
automatically collected into a remainder group appended as the final
element of `groups_list`. Line type and line width are mapped to group
membership via `line_type_groups` and `line_size_groups`, both of which
must have length equal to the total number of groups (explicit plus the
automatic remainder). At most six groups are supported.

**`max_class` column.** The returned `$Data_Frame` carries a `max_class`
column recording the label of the highest-entropy GMM component for each
partition, taken directly from `part_data$Max_Entropy[i]`. This column
is consumed by
[`plot_site_class_trajectory`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_site_class_trajectory.md)
(red labels) and by downstream class-assignment tables.

## See also

[`partition_time_windows`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md),
[`relabel_entropy_classes`](https://vadimtyuryaev.github.io/ViralEntropR/reference/relabel_entropy_classes.md),
[`plot_site_class_trajectory`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_site_class_trajectory.md)

## Examples

``` r
# Synthetic dataset: site 2 accumulates variability across partitions.
# Partition entropies (bits):
#   Jan  s1: 0.000  s2: 0.000  — both removed (zero entropy)
#   Feb  s1: 0.000  s2: 0.722  — s1 removed, s2 retained
#   Mar  s1: 0.971  s2: 1.000  — both retained
#   Apr  s1: 0.881  s2: 0.722  — both retained
df <- data.frame(
  s1 = c(rep(1L, 10L),
         rep(1L, 10L),
         c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
         c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)),
  s2 = c(rep(1L, 10L),
         c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L),
         c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L),
         c(1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L)),
  Date = rep(
    seq(as.Date("2020-01-01"), by = "month", length.out = 4L),
    each = 10L
  )
)

part_data <- partition_time_windows(
  data          = df,
  n_sites       = 2L,
  window_length = 1L,
  window_type   = 3L,
  start_date    = "2020-01-01",
  end_date      = "2020-04-01"
)

# Example 1: no relabeling — class labels are as assigned by Mclust.
result <- plot_entropy_trajectories(
  part_data  = part_data
)
print(result$Plot)
```
