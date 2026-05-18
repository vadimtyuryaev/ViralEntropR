## Resubmission

This is a patch release (0.6.1) of ViralEntropR addressing the items
raised in the CRAN initial-review feedback from Konstanze Lauseker
(2026-05-XX):

* Expanded "NCBI" on first use in the DESCRIPTION text.

* Removed `@examples` blocks from three internal functions
  (`get_site_counts`, `build_distance_matrix`, `get_variants`),
  eliminating the `:::` usage previously required to call them and
  removing examples documenting non-exported functions.

* Replaced `\dontrun{}` with `\donttest{}` in examples for
  `extract_fasta_dates()`, `extract_fasta_countries()`,
  `fasta_to_char_matrix()`, and the `sarscov2_sample` dataset. These
  examples depend on `Biostrings` (Suggests); `\donttest{}` is used
  because CRAN runs Suggests-dependent examples while end-user
  invocations may not have `Biostrings` installed.

* Removed `getwd()` as a default for `save_path` in
  `plot_site_class_trajectory()` and `tabulate_site_evolution()`.
  `save_path` is now `NULL` by default and must be supplied
  explicitly when `save = TRUE`; no function writes to disk by
  default or in any example, vignette, or test.

## Additional fixes discovered during local verification

While running `devtools::run_examples()` to verify the items above,
four small issues were found and addressed in this resubmission:

* `partition_time_windows()` previously emitted spurious "windows
  will be empty" warnings whenever `end_date` exceeded
  `max(data$Date)` by any amount, including cases where the last
  window correctly contained the trailing data. The predicate has
  been corrected to compare against the last window's left boundary,
  so the warning fires only when at least one window is genuinely
  empty. The partition output was always correct.

* `detect_changepoints_ecp()` example referenced `$cpLoc[[1]]`,
  which is the zero-change-point case (always empty under
  `ks.cp3o`'s indexing). The example line has been removed to
  avoid user confusion.

* `tabulate_site_evolution()` example previously inspected the
  `kableExtra` HTML object directly, which prints as raw HTML
  text outside R Markdown contexts; changed to display the
  underlying data frame.

* `plot_entropy_trajectories()` and `plot_site_class_trajectory()`
  examples have been redesigned with three-period synthetic datasets
  with distinct trajectory shapes, producing more informative
  reference output. The `plot_entropy_trajectories()` x-axis scale
  now uses `expansion(add = 0.6)` so rotated tick labels at the
  leftmost partition are not clipped.

No user-visible API changes beyond the `save_path` default; no new
dependencies; no policy-relevant changes.

## Test environments

* local: Windows 11, R 4.5.2
* GitHub Actions:
  - ubuntu-latest (R-devel, R-release, R-oldrel-1)
  - macos-latest (R-release)
  - windows-latest (R-release)
* win-builder: R-devel, R-release

## R CMD check results

0 errors | 0 warnings | 0 notes