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
  `fasta_to_char_matrix()`, and the `sarscov2_sample` dataset.
  These examples depend on `Biostrings` (Suggests); `\donttest{}`
  is used because CRAN runs Suggests-dependent examples while
  end-user invocations may not have `Biostrings` installed.
* Removed `getwd()` as a default for `save_path` in
  `plot_site_class_trajectory()` and `tabulate_site_evolution()`.
  `save_path` is now `NULL` by default and must be supplied
  explicitly when `save = TRUE`; no function writes to disk by
  default or in any example, vignette, or test.

No user-visible API changes beyond the `save_path` default; no new
dependencies; no policy-relevant changes.

## Test environments

* local: Windows 11, R 4.x.x
* GitHub Actions:
  - ubuntu-latest (R-devel, R-release, R-oldrel-1)
  - macos-latest (R-release)
  - windows-latest (R-release)
* win-builder: R-devel, R-release

## R CMD check results

0 errors | 0 warnings | 0 notes