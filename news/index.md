# Changelog

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
