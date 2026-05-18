# Changelog

## ViralEntropR 0.6.1

- Additional fixes discovered during local verification of the
  resubmission:
  - Fix spurious “windows will be empty” warnings in
    [`partition_time_windows()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md).
    The warning predicate previously compared `end_date` against
    `max(data$Date)` directly, which fires whenever `end_date` is at all
    past `max(data$Date)` — even when the last window correctly contains
    the trailing data. The fix compares against the left boundary of the
    last window instead, so the warning fires only when one or more
    trailing windows would actually be empty. Mirror logic applies to
    `start_date`. The function’s partition output was always correct;
    only the warning fires more accurately now.
  - Update
    [`detect_changepoints_ecp()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_ecp.md)
    example: last example referenced `$cpLoc[[1]]`, which is the
    zero-change-point case (always empty under `ks.cp3o`’s indexing).
    The example line has been removed to avoid user confusion.
  - Improve example clarity for
    [`tabulate_site_evolution()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/tabulate_site_evolution.md):
    display the underlying data frame (`$table`) rather than the
    `kableExtra` HTML object (`$styled`), which prints as raw HTML text
    outside R Markdown contexts. The styled HTML still renders correctly
    in vignettes and interactive use.
  - Redesign the examples for
    [`plot_entropy_trajectories()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_entropy_trajectories.md)
    and
    [`plot_site_class_trajectory()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_site_class_trajectory.md)
    to use three-period synthetic datasets with distinct trajectory
    shapes and visible class changes across partitions, producing more
    pedagogically useful reference output. Added `expansion(add = 0.6)`
    to the
    [`plot_entropy_trajectories()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_entropy_trajectories.md)
    x-axis scale so rotated tick labels at the leftmost partition are no
    longer clipped.

## ViralEntropR 0.6.1

- Address CRAN initial-review feedback:
  - Expand “NCBI” on first use in DESCRIPTION.
  - Remove `@examples` blocks from internal functions
    (`get_site_counts`, `build_distance_matrix`, `get_variants`), which
    eliminates the `:::` usage previously required to call them.
  - Replace `\dontrun{}` with `\donttest{}` in the examples for
    [`extract_fasta_dates()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/extract_fasta_dates.md),
    [`extract_fasta_countries()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/extract_fasta_countries.md),
    [`fasta_to_char_matrix()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/fasta_to_char_matrix.md),
    and the `sarscov2_sample` dataset.
  - Remove [`getwd()`](https://rdrr.io/r/base/getwd.html) as a default
    for `save_path` in
    [`plot_site_class_trajectory()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_site_class_trajectory.md)
    and
    [`tabulate_site_evolution()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/tabulate_site_evolution.md);
    `save_path` is now `NULL` by default and must be supplied explicitly
    when `save = TRUE`.

## ViralEntropR 0.6.0

- Initial CRAN submission.
