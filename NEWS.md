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
