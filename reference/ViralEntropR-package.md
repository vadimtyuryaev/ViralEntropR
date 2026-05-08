# ViralEntropR: A Computational Pipeline for Entropy-Informed Detection of Emerging Viral Variants

A computational pipeline for detecting emerging variants in viral amino
acid sequence data, combining per-site Shannon entropy, Gaussian mixture
model site selection, Gower-distance Partitioning Around Medoids
clustering, Hellinger-distance quantification of distributional shifts,
and multivariate non-parametric change-point detection.

## Pipeline overview

The package supports a four-stage workflow:

- **Preprocessing.** Parse FASTA headers, filter ambiguous residues, and
  convert between integer and character representations of amino acid
  sequences under a 25-symbol alphabet. See
  [`extract_fasta_dates`](https://vadimtyuryaev.github.io/ViralEntropR/reference/extract_fasta_dates.md),
  [`extract_fasta_countries`](https://vadimtyuryaev.github.io/ViralEntropR/reference/extract_fasta_countries.md),
  [`fasta_to_char_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/fasta_to_char_matrix.md),
  [`filter_ambiguous_sequences`](https://vadimtyuryaev.github.io/ViralEntropR/reference/filter_ambiguous_sequences.md),
  [`encode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md),
  and
  [`decode_aa_sequence`](https://vadimtyuryaev.github.io/ViralEntropR/reference/decode_aa_sequence.md).

- **Site selection.** Compute per-site Shannon entropy across temporal
  partitions and cluster sites by entropy via Gaussian mixture models.
  See
  [`calculate_entropy`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_entropy.md),
  [`partition_time_windows`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md),
  [`cluster_sites_by_entropy`](https://vadimtyuryaev.github.io/ViralEntropR/reference/cluster_sites_by_entropy.md),
  and
  [`relabel_entropy_classes`](https://vadimtyuryaev.github.io/ViralEntropR/reference/relabel_entropy_classes.md).

- **Distributional analysis.** Quantify residue-composition shifts
  between time windows using the Hellinger distance. See
  [`calculate_hellinger_matrix`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_hellinger_matrix.md).

- **Change-point detection.** Identify temporal change points
  non-parametrically using energy statistics or wild binary
  segmentation. See
  [`detect_changepoints_ecp`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_ecp.md)
  and
  [`detect_changepoints_hdcp`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_hdcp.md).

## Visualisation and tabulation

- [`plot_entropy_trajectories`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_entropy_trajectories.md)
  — customizable multi-site Shannon entropy trajectories, summarizing
  evolutionary dynamics across time.

- [`plot_site_class_trajectory`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_site_class_trajectory.md)
  — single-site entropy trajectory with class-change markers for
  inspecting individual residues of interest.

- [`tabulate_site_evolution`](https://vadimtyuryaev.github.io/ViralEntropR/reference/tabulate_site_evolution.md)
  — per-site amino acid count and proportion tables across partitions,
  optionally rendered as styled HTML.

## Simulation

[`simulate_variant_evolution`](https://vadimtyuryaev.github.io/ViralEntropR/reference/simulate_variant_evolution.md)
provides a configurable multi-variant simulation engine with
user-specified emergence schedules, growth rates, mutation-rate
variability, and deleterious-mutation injection. Generates synthetic
time-series data with known ground truth for benchmarking detection
pipelines.

## Bundled data

- [`sarscov2_variants`](https://vadimtyuryaev.github.io/ViralEntropR/reference/sarscov2_variants.md)
  — curated metadata for twelve SARS-CoV-2 Variants of Concern and
  Variants of Interest, including WHO labels, Pango lineages, GISAID and
  Nextstrain clades, dates and countries of first detection, defining
  Spike-protein mutations and SNP sites, and 21 peer-reviewed references
  with DOIs.

- [`sarscov2_sample`](https://vadimtyuryaev.github.io/ViralEntropR/reference/sarscov2_sample.md)
  — a random sample of 100 NCBI Spike protein sequences for end-to-end
  testing without external downloads.

## External data

The complete preprocessed NCBI Spike protein dataset (137,132 sequences,
~181.5 MB uncompressed FASTA) underlying the package's real-data
pre-processing vignette is archived on Zenodo:
[doi:10.5281/zenodo.19040165](https://doi.org/10.5281/zenodo.19040165) .
The dataset can be read directly with
[`readAAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-io.html)
and processed end-to-end using the preprocessing toolkit; see
`vignette("preprocessing_pipeline", "ViralEntropR")` for the full
workflow.

## Vignettes

Three pre-rendered vignettes walk through the full workflow on real and
simulated data: `vignette("preprocessing_pipeline", "ViralEntropR")`,
`vignette("detecting_variants_simulation", "ViralEntropR")`, and
`vignette("clustering_accuracy", "ViralEntropR")`.

## See also

Useful links:

- <https://github.com/vadimtyuryaev/ViralEntropR>

- [doi:10.5281/zenodo.19040165](https://doi.org/10.5281/zenodo.19040165)

- Report bugs at <https://github.com/vadimtyuryaev/ViralEntropR/issues>

## Author

**Maintainer**: Vadim Tyuryaev <vadim.tyuryaev@gmail.com>
([ORCID](https://orcid.org/0009-0008-1361-6265))

Authors:

- Jane Heffernan <jmheffer@yorku.ca>

- Hanna Jankowski <hkj@yorku.ca>
