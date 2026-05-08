# Package index

## Package overview

- [`ViralEntropR`](https://vadimtyuryaev.github.io/ViralEntropR/reference/ViralEntropR-package.md)
  [`ViralEntropR-package`](https://vadimtyuryaev.github.io/ViralEntropR/reference/ViralEntropR-package.md)
  : ViralEntropR: A Computational Pipeline for Entropy-Informed
  Detection of Emerging Viral Variants

## Preprocessing

Fully vectorised FASTA-to-feature-matrix preprocessing layer.

- [`extract_fasta_dates()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/extract_fasta_dates.md)
  : Extract Dates from FASTA Sequence Names
- [`extract_fasta_countries()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/extract_fasta_countries.md)
  : Extract Countries from FASTA Sequence Names
- [`fasta_to_char_matrix()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/fasta_to_char_matrix.md)
  : Convert FASTA Object to Character Matrix
- [`filter_ambiguous_sequences()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/filter_ambiguous_sequences.md)
  : Remove Sequences Containing Ambiguous Residues
- [`encode_aa_sequence()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/encode_aa_sequence.md)
  : Encode Amino Acid Sequences
- [`decode_aa_sequence()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/decode_aa_sequence.md)
  : Decode Amino Acid Sequences

## Site selection

Per-site Shannon entropy, Gaussian mixture model classification, and
temporal partitioning into cumulative, sliding, or disjoint windows.

- [`calculate_entropy()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_entropy.md)
  : Calculate Shannon Entropy
- [`cluster_sites_by_entropy()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/cluster_sites_by_entropy.md)
  : Cluster a Univariate Numeric Vector by Gaussian Mixture Model
- [`relabel_entropy_classes()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/relabel_entropy_classes.md)
  : Relabel Entropy Classes
- [`partition_time_windows()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/partition_time_windows.md)
  : Partition Data into Time Windows

## Distance and change-point detection

Pairwise Hellinger distances and non-parametric change-point detection
(energy statistics; wild binary segmentation).

- [`calculate_hellinger_matrix()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/calculate_hellinger_matrix.md)
  : Calculate Hellinger Distance Matrix
- [`detect_changepoints_ecp()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_ecp.md)
  : Detect Temporal Change Points (ECP)
- [`detect_changepoints_hdcp()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/detect_changepoints_hdcp.md)
  : Detect Temporal Change Points (HDcpDetect)

## Visualisation and tabulation

Entropy trajectories, per-site class trajectories, and amino-acid
frequency tables.

- [`plot_entropy_trajectories()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_entropy_trajectories.md)
  : Plot Shannon Entropy Trajectories
- [`plot_site_class_trajectory()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/plot_site_class_trajectory.md)
  : Plot GMM Entropy Class Trajectory for a Single Site
- [`tabulate_site_evolution()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/tabulate_site_evolution.md)
  : Tabulate Site Frequency Evolution

## Simulation

Multi-variant evolutionary simulator for benchmarking detection.

- [`simulate_variant_evolution()`](https://vadimtyuryaev.github.io/ViralEntropR/reference/simulate_variant_evolution.md)
  : Simulate Viral Variant Evolution

## Data

Bundled datasets and example sequences for runnable examples and
reproducible analyses.

- [`sarscov2_variants`](https://vadimtyuryaev.github.io/ViralEntropR/reference/sarscov2_variants.md)
  : SARS-CoV-2 VOC/VOI Curated Variant Metadata
- [`sarscov2_sample`](https://vadimtyuryaev.github.io/ViralEntropR/reference/sarscov2_sample.md)
  : SARS-CoV-2 Surface Glycoprotein Sequences – NCBI Demo Sample
