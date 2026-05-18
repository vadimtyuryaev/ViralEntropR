# ViralEntropR 0.6.1

* Additional fixes discovered during local verification of the
  resubmission:
  - Fix spurious "windows will be empty" warnings in
    `partition_time_windows()`. The warning predicate previously
    compared `end_date` against `max(data$Date)` directly, which
    fires whenever `end_date` is at all past `max(data$Date)` — even
    when the last window correctly contains the trailing data. The
    fix compares against the left boundary of the last window
    instead, so the warning fires only when one or more trailing
    windows would actually be empty. Mirror logic applies to
    `start_date`. The function's partition output was always correct;
    only the warning fires more accurately now.
  - Update `detect_changepoints_ecp()` example: last example referenced 
    `$cpLoc[[1]]`, which is the zero-change-point case (always empty under
    `ks.cp3o`'s indexing). The example line has been removed to 
    avoid user confusion.
  - Improve example clarity for `tabulate_site_evolution()`: display
    the underlying data frame (`$table`) rather than the
    `kableExtra` HTML object (`$styled`), which prints as raw HTML
    text outside R Markdown contexts. The styled HTML still renders
    correctly in vignettes and interactive use.
  - Redesign the examples for `plot_entropy_trajectories()` and
    `plot_site_class_trajectory()` to use three-period synthetic
    datasets with distinct trajectory shapes and visible class
    changes across partitions, producing more pedagogically useful
    reference output. Added `expansion(add = 0.6)` to the
    `plot_entropy_trajectories()` x-axis scale so rotated tick
    labels at the leftmost partition are no longer clipped.

# ViralEntropR 0.6.1

* Address CRAN initial-review feedback:
  - Expand "NCBI" on first use in DESCRIPTION.
  - Remove `@examples` blocks from internal functions
    (`get_site_counts`, `build_distance_matrix`, `get_variants`), which
    eliminates the `:::` usage previously required to call them.
  - Replace `\dontrun{}` with `\donttest{}` in the examples for
    `extract_fasta_dates()`, `extract_fasta_countries()`,
    `fasta_to_char_matrix()`, and the `sarscov2_sample` dataset.
  - Remove `getwd()` as a default for `save_path` in
    `plot_site_class_trajectory()` and `tabulate_site_evolution()`;
    `save_path` is now `NULL` by default and must be supplied
    explicitly when `save = TRUE`.

# ViralEntropR 0.6.0

* Initial CRAN submission.
