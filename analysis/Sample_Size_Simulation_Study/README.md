Sample-Size Requirements for Entropy-Based Detection of Emerging Viral
Variants: A Simulation Study with the ViralEntropR Pipeline
================
Vadim Tyuryaev, Jane Heffernan, Hanna Jankowski
May 19, 2026

- [Overview](#overview)
- [Terminology](#terminology)
- [Research questions](#research-questions)
- [Pipeline overview](#pipeline-overview)
- [Empirical foundations](#empirical-foundations)
  - [Mutation counts per variant](#mutation-counts-per-variant)
  - [Position pools](#position-pools)
  - [Cross-variant collision
    properties](#cross-variant-collision-properties)
- [Deleterious-mutation noise](#deleterious-mutation-noise)
  - [Mechanism](#mechanism)
  - [Empirical justification of
    $p_{\text{del}} = 10^{-3}$](#empirical-justification-of-p_textdel--10-3)
- [Scenarios](#scenarios)
  - [Scenario 1: detect $V_2$ on a D614G-dominant
    background](#scenario-1-detect-v_2-on-a-d614g-dominant-background)
  - [Scenario 2: detect $V_3$ on a $V_1 + V_2$
    background](#scenario-2-detect-v_3-on-a-v_1--v_2-background)
  - [Scenario 3: detect $V_4$ (non-Omicron-like) on a $V_1 + V_2 + V_3$
    background](#scenario-3-detect-v_4-non-omicron-like-on-a-v_1--v_2--v_3-background)
  - [Scenario 4: detect $V_4$ (Omicron-like) on a $V_1 + V_2 + V_3$
    background](#scenario-4-detect-v_4-omicron-like-on-a-v_1--v_2--v_3-background)
- [Sampling structure](#sampling-structure)
  - [Axis 1: band](#axis-1-band)
  - [Axis 2: $N_{\text{ref}}$](#axis-2-n_textref)
  - [Axes 3, 4, and 5: ratios](#axes-3-4-and-5-ratios)
    - [Bands describe $N_{\text{ref}}$ scale
      only](#bands-describe-n_textref-scale-only)
  - [Cells and simulations per
    scenario](#cells-and-simulations-per-scenario)
  - [What is fixed and what is random within a
    cell](#what-is-fixed-and-what-is-random-within-a-cell)
  - [Why 30 replicates per cell](#why-30-replicates-per-cell)
  - [Why full factorial rather than stratified random
    sampling](#why-full-factorial-rather-than-stratified-random-sampling)
- [Sample-size sweep procedure](#sample-size-sweep-procedure)
  - [Ceiling](#ceiling)
  - [Grid construction](#grid-construction)
    - [Representative grid sizes](#representative-grid-sizes)
    - [Explicit grid values for representative
      ceilings](#explicit-grid-values-for-representative-ceilings)
  - [Grid-precision answers and within-cell
    aggregation](#grid-precision-answers-and-within-cell-aggregation)
  - [Per-replicate GMM-fit budget](#per-replicate-gmm-fit-budget)
- [Detection rule](#detection-rule)
- [Reproducibility](#reproducibility)
- [Quality control and post-hoc
  reconstruction](#quality-control-and-post-hoc-reconstruction)
- [Output structure](#output-structure)
- [Memory budget](#memory-budget)
- [Runtime budget](#runtime-budget)
- [Invocation](#invocation)
- [Methodology framing (ADEMP)](#methodology-framing-ademp)
- [Limitations](#limitations)
- [Session information](#session-information)
- [References](#references)

# Overview

This document is the analysis plan for the Sample-Size Simulation Study,
one of several end-to-end analyses accompanying the `ViralEntropR`
package and located in the `analysis/` folder.

The study addresses a single applied question. Genomic-surveillance
programmes continuously collect viral sequence samples and process them
through analytical pipelines that flag positions of unusual variability.
When a new viral variant begins circulating, how many sequences of that
variant must a surveillance pipeline observe before it can reliably flag
the variant’s mutation sites? This number is the *required sample size
for variant detection*, and its value depends on properties of both the
background population (how many wild-type reference sequences are
present, how many co-circulating variants are already established) and
the emerging variant itself (how many mutations it carries, and which
positions it mutates).

The pipeline evaluated here is the per-site Shannon-entropy plus
Gaussian Mixture Model site-selection step that forms the core of the
`ViralEntropR` analytical workflow. The study constructs synthetic
sequence populations under four scenarios of increasing background
complexity, sweeps the count of emerging-variant sequences upward, and
records the smallest count at which all of the variant’s mutation sites
land in the highest-entropy cluster of the model. All biological
parameters of the simulator are derived from empirical SARS-CoV-2 data.
The number of mutations per variant, the mutated Spike positions, and
the substitute amino acids drawn at each position come from
`sarscov2_variants` (see `?sarscov2_variants`), the curated catalogue of
WHO-labelled SARS-CoV-2 variants assembled from 21 peer-reviewed sources
and bundled with the package ([Tyuryaev, Heffernan, and Jankowski
2026](#ref-Tyuryaev2026)). The surveillance-level deleterious-mutation
rate $p_{\text{del}}$ is estimated directly from the post-processed NCBI
and GISAID US Spike feature matrices consumed downstream by the rest of
the `ViralEntropR` pipeline (Section 6.2).

The design is a full discrete factorial ([Montgomery
2017](#ref-Montgomery2017); [Baker et al. 2017](#ref-Baker2017)) across
structural axes (defined in the Terminology section below) with 30
replicates per cell. Each replicate draws biological parameters
independently from empirical distributions and applies independent
low-frequency deleterious-mutation noise to the simulated population.
The study does not include temporal partitioning, Hellinger-distance
computation, or change-point detection. Those components of the broader
`ViralEntropR` pipeline are evaluated separately in
`Entropy_Partitioning_Study/` and `CP_Detection_Study/` and discussed in
the
[`detecting_variants_simulation`](https://vadimtyuryaev.github.io/ViralEntropR/articles/detecting_variants_simulation.html)
vignette.

# Terminology

The following terms appear throughout the document. Biological terms are
defined for readers unfamiliar with viral genomics; statistical and
study-specific terms are defined for readers from other fields.

| Term | Definition |
|----|----|
| Spike protein | The 1,273-amino-acid surface glycoprotein of SARS-CoV-2. |
| Wild type | The unmutated reference Spike sequence [GenBank accession YP_009724390](https://www.ncbi.nlm.nih.gov/protein/1796318598). The background against which variants are compared. |
| Variant | A sequence carrying one or more amino-acid substitutions relative to the reference. |
| Mutation | An amino-acid substitution at a specific position, written as (reference residue)(position)(new residue). For example, D614G means aspartate (D) at position 614 replaced by glycine (G). |
| D614G | SARS-CoV-2 variant that effectively became the global background sequence after mid-2020 ([Korber et al. 2020](#ref-korber2020); [Plante et al. 2021](#ref-plante2020)). The simulator treats D614G as a deterministic baseline mutation in every variant of every scenario in this revised design. |
| Variant of Concern (VOC), Variant of Interest (VOI) | [World Health Organization classifications](https://www.who.int/publications/m/item/historical-working-definitions-and-primary-actions-for-sars-cov-2-variants). The `sarscov2_variants` catalogue contains 12 such variants. |
| Emerging variant | The variant whose detection the study evaluates. In each replicate, the emerging variant’s sample size is swept upward until detection succeeds. |
| Established (dominant) variant | A variant already circulating in the population when the emerging variant appears. Denoted $V_1$ (always D614G in this study), $V_2$, and $V_3$ in the multi-variant scenarios. |
| Reference background | The wild-type sequences in a simulated population, of size $N_{\text{ref}}$. |
| Per-site Shannon entropy | An information-theoretic measure of compositional uncertainty at a single sequence position, computed across all rows of the sequence matrix as $H_s = -\sum_{r} p_{s,r} \log_2 p_{s,r}$ where $p_{s,r}$ is the observed proportion of residue $r$ at position $s$ ([Shannon 1948](#ref-Shannon1948)). Entropy is zero at the endpoints where a single residue dominates (one variant fully replaced by another), grows with the number of co-occurring residues, and reaches its maximum $\log_2 K$ when $K$ residues are present in equal proportions — the configuration of maximum compositional uncertainty. The study computes one entropy value per Spike position (1,273 in total). |
| Gaussian Mixture Model (GMM) | A clustering method that fits a weighted sum of Gaussian components to a univariate sample (here, the 1,273 entropy values). The package’s `cluster_sites_by_entropy` function wraps `mclust::Mclust` ([Scrucca et al. 2016](#ref-scrucca2016)) to fit a one-dimensional GMM and classify each site into a component. |
| Class 1 (highest-entropy class) | After fitting the GMM and applying `relabel_entropy_classes`, the class containing the highest-entropy sites is labelled class 1. These are the sites most variable across the population. |
| Detection | The binary outcome of the study. A replicate is detected when every mutation site of the emerging variant appears in class 1 of the relabeled GMM output classes. |
| Replicate | One simulation run, indexed by `(scenario, cell, rep)`. Each cell contains 30 replicates with independent biological draws and independent noise patterns. |
| Cell | A specific combination of structural parameters. The set of cells per scenario forms a discrete grid (the factorial design). |
| Band | A category for the order of magnitude of the reference background size, $N_{\text{ref}}$: small (10 to 100), medium (100 to 1,000), large (1,000 to 10,000). The band describes the surveillance regime of the wild-type background only; the dominant variants $V_1$, $V_2$, and $V_3$ are not constrained to lie within the band. |
| $N_{\text{ref}}$ | Number of wild-type reference sequences in the simulated population. |
| $n_{V_1}$, $n_{V_2}$, $n_{V_3}$ | Number of $V_1$, $V_2$, $V_3$ sequences in the simulated population (in the scenarios where each is present). |
| $n_{\text{emerge}}$ | Number of emerging-variant sequences. The variable swept in each replicate. |
| $\text{ratio}_{V_i}$ | Multiplier used to derive $n_{V_i}$ from $N_{\text{ref}}$: $n_{V_i} = \text{round}(N_{\text{ref}} \times \text{ratio}_{V_i})$. The ratios are fixed values greater than 1, so the multipliers ensure $n_{V_i} > N_{\text{ref}}$ by construction (each dominant variant exceeds the reference background it competes against). |
| Sweep | The procedure of incrementing $n_{\text{emerge}}$ stepwise and testing detection at each step. This revised design uses a 50-point log-spaced grid across all three bands, with no refinement (grid-precision threshold reporting). |
| GISAID | [Global Initiative on Sharing All Influenza Data](https://gisaid.org). A public-access platform hosting SARS-CoV-2 (and other respiratory virus) genome submissions with curated metadata. The GISAID US Spike feature matrix used in this study is multiple-sequence-aligned to a 1,273-position reference and contains gap characters at deletion sites ([Shu and McCauley 2017](#ref-shu2017)). |
| NCBI | [National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov). The NIH institute hosting GenBank, the primary US-mirrored repository of SARS-CoV-2 sequence submissions. The NCBI US Spike feature matrix used in this study is unaligned and contains only standard amino-acid codes at each of the 1,273 Spike positions ([Sayers et al. 2022](#ref-sayers2022)). |
| FASTA | A plain-text format for biological sequence data in which each record begins with a single-line header (prefixed `>`) carrying the sequence identifier and metadata, followed by one or more lines containing the sequence itself in single-letter IUPAC amino-acid or nucleotide codes ([Pearson and Lipman 1988](#ref-PearsonLipman1988)). FASTA files are the input format for the `ViralEntropR` preprocessing pipeline, which parses, quality-filters, and encodes them into the integer feature matrices consumed downstream. |
| $p_{\text{del}}$ | Per-(sequence, position) Bernoulli probability that a residue is replaced with a uniformly chosen non-self standard amino acid. Models surveillance-level prevalence of transient deleterious substitutions. Default $10^{-3}$, calculated empirically based on 0.1%-cutoff rate from post-processed NCBI and GISADI feature matrices. See the `preprocessing_pipeline` vignette and `analysis/GISAID_data_preprocessing/` for the full preprocessing pipeline applied to the raw FASTA submissions before the feature matrix is built. |
| `MUT_COUNT_SET_11` | Empirical set of unique Spike-protein mutation counts across the 11 non-Omicron variants in `sarscov2_variants`. |
| `MUT_COUNT_OMICRON` | Omicron’s Spike-protein mutation count. |
| `POOL_11`, `POOL_12` | Empirical sets of Spike positions that are mutated in at least one variant. `POOL_11` is the union across the 11 non-Omicron variants; `POOL_12` is the union across all 12 variants. |
| `POOL_11_DRAW`, `POOL_12_DRAW` | `POOL_11` and `POOL_12` with position 614 (D614G) removed, used as random-draw pools when D614G is forced as a baseline mutation. |
| ADEMP framework | A standard reporting structure for simulation studies in statistical methodology: Aims, Data-generating mechanism, Estimands, Methods, Performance measures ([Morris, White, and Crowther 2019](#ref-morris2019)). |

# Research questions

**Q1.** How many de-novo emerging-variant sequences are required for the
GMM-based site-selection step to detect a new variant against a
non-trivial D614G-dominant background?

**Q2.** How does the required de-novo count change as the background
gains additional co-circulating variants of comparable or higher
abundance?

**Q3.** How does the required de-novo count depend on the emerging
variant’s mutational load, ranging from a small-to-moderate
non-Omicron-like profile (drawn from `MUT_COUNT_SET_11`) to an
Omicron-like profile of 33 mutations?

**Q4.** How robust are the inferred thresholds to the prevalence of
transient deleterious-mutation noise, varied across $p_{\text{del}}
\in \{0, 10^{-4}, 10^{-3}, 10^{-2}\}$ in a stratified-random subset of
cells?

# Pipeline overview

The simulator and the detection procedure operate together for each
replicate. The simulator builds a synthetic sequence population from the
structural parameters defined by the cell, then injects independent
low-frequency deleterious-mutation noise. The detection procedure sweeps
the number of emerging-variant sequences upward, testing at each step
whether the pipeline can identify the variant’s mutation sites.

The procedure for one replicate is as follows.

1.  Draw the per-replicate biological parameters. These are the mutation
    counts for $V_2$ and $V_3$ (in scenarios where each is present) and
    the emerging variant, the specific Spike positions each variant
    mutates, and the substitute amino acid placed at each mutated
    position. These draws are independent across the 30 replicates of a
    cell.

2.  Construct an integer-encoded sequence matrix at the cell’s ceiling
    row count. The matrix has one row per sequence and one column per
    Spike position (1,273 columns). Reference sequences contribute
    $N_{\text{ref}}$ rows; $V_1$ contributes $n_{V_1}$ rows; $V_2$ and
    $V_3$ each contribute their respective $n_{V_i}$ rows where present;
    the emerging variant contributes a number of rows equal to the
    cell’s ceiling.

3.  Apply deleterious-mutation noise once to the full matrix at rate
    $p_{\text{del}}$. Each cell of the matrix is independently
    overwritten with probability $p_{\text{del}}$ by a uniformly drawn
    non-self standard amino acid. The noisy matrix is the property of
    the simulated dataset; it is not regenerated across sweep points.

4.  For each candidate $n_{\text{emerge}}$ in the 50-point log-spaced
    sweep grid (in ascending order):

    - Treat the leading $N_{\text{ref}} + n_{V_1} + n_{V_2} + n_{V_3} +
      n_{\text{emerge}}$ rows of the matrix as the population at this
      sweep point.
    - Compute the per-site Shannon entropy across those rows.
    - Fit a one-dimensional Gaussian mixture (mclust defaults: E (equal
      variance) and V (variable variance) models evaluated together,
      $G \in 1{:}15$, BIC selection).
    - Relabel the fitted classes so that class 1 contains the
      highest-entropy sites.
    - Test detection: does every mutated position of the emerging
      variant belong to class 1? With strictly-equal degeneracy handlers
      for the $G = 1$ and all-identical (sentinel 999) cases.

5.  Record the smallest $n_{\text{emerge}}$ at which detection succeeded
    and stop the sweep. If the sweep reaches its ceiling without
    detection, record `NA`.

The sweep uses a 50-point log-spaced grid evaluated in ascending order
with no refinement; the smallest detecting grid point is the reported
threshold (grid-precision answer). See the Sample-size sweep procedure
section for explicit grid values per band.

# Empirical foundations

All biological parameters of the simulator derive from
`sarscov2_variants` and post-processed NCBI/GISAID feature matrices (see
Terminology section above).

## Mutation counts per variant

``` r
labels      <- unlist(sarscov2_variants$WHO_Label)
mut_sites   <- sarscov2_variants$Mutation_Sites
mut_strings <- sarscov2_variants$Spike_Mutations

mut_counts <- vapply(mut_sites, length, integer(1L))
names(mut_counts) <- labels

mut_counts_df <- data.frame(
  Variant                 = names(sort(mut_counts)),
  `Spike mutation count`  = unname(sort(mut_counts)),
  check.names             = FALSE
)

kbl(mut_counts_df,
    caption = "Spike-protein mutation counts per WHO-labelled variant.",
    align   = "lc") |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width        = FALSE)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Spike-protein mutation counts per WHO-labelled variant.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Variant
</th>

<th style="text-align:center;">

Spike mutation count
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Zeta
</td>

<td style="text-align:center;">

2
</td>

</tr>

<tr>

<td style="text-align:left;">

Lambda
</td>

<td style="text-align:center;">

3
</td>

</tr>

<tr>

<td style="text-align:left;">

Epsilon
</td>

<td style="text-align:center;">

4
</td>

</tr>

<tr>

<td style="text-align:left;">

Theta
</td>

<td style="text-align:center;">

4
</td>

</tr>

<tr>

<td style="text-align:left;">

Eta
</td>

<td style="text-align:center;">

8
</td>

</tr>

<tr>

<td style="text-align:left;">

Kappa
</td>

<td style="text-align:center;">

8
</td>

</tr>

<tr>

<td style="text-align:left;">

Beta
</td>

<td style="text-align:center;">

12
</td>

</tr>

<tr>

<td style="text-align:left;">

Gamma
</td>

<td style="text-align:center;">

12
</td>

</tr>

<tr>

<td style="text-align:left;">

Alpha
</td>

<td style="text-align:center;">

14
</td>

</tr>

<tr>

<td style="text-align:left;">

Iota
</td>

<td style="text-align:center;">

14
</td>

</tr>

<tr>

<td style="text-align:left;">

Delta
</td>

<td style="text-align:center;">

15
</td>

</tr>

<tr>

<td style="text-align:left;">

Omicron
</td>

<td style="text-align:center;">

33
</td>

</tr>

</tbody>

</table>

``` r
MUT_COUNT_SET_11  <- sort(unique(mut_counts[labels != "Omicron"]))
MUT_COUNT_OMICRON <- unname(mut_counts["Omicron"])
```

Two empirical sets of mutation counts follow from this table.
`MUT_COUNT_SET_11` is the set of unique counts across the 11 non-Omicron
variants and equals {2, 3, 4, 8, 12, 14, 15}. Per replicate, the
simulator draws independently from this set for the mutation count of
$V_2$ (in scenarios 2, 3, 4), $V_3$ (in scenarios 3, 4), and the
non-Omicron emerging variant (in scenarios 1, 2, 3); $V_1$ is excluded
from this draw because it is fixed at a single mutation, the D614G
baseline at position 614, in every replicate of every scenario.
`MUT_COUNT_OMICRON` equals 33 and is the deterministic mutation count of
the Omicron-like emerging variant in scenario 4.

## Position pools

``` r
POOL_11      <- sort(unique(unlist(mut_sites[labels != "Omicron"])))
POOL_12      <- sort(unique(unlist(mut_sites)))
POOL_OMICRON <- setdiff(POOL_12, POOL_11)

stopifnot(614L %in% POOL_11, 614L %in% POOL_12)

POOL_11_DRAW <- setdiff(POOL_11, 614L)
POOL_12_DRAW <- setdiff(POOL_12, 614L)

pools_df <- data.frame(
  Pool   = c("|POOL_11|",
             "|POOL_12|",
             "|POOL_12 \\\\ POOL_11|",
             "|POOL_11 \\\\ {614}|  (random-draw pool, non-Omicron variants)",
             "|POOL_12 \\\\ {614}|  (random-draw pool, Sc 4 emerging)"),
  Size   = c(length(POOL_11),
             length(POOL_12),
             length(POOL_OMICRON),
             length(POOL_11_DRAW),
             length(POOL_12_DRAW)),
  Notes  = c("union of mutation sites across the 11 non-Omicron variants",
             "union across all 12 variants",
             "Omicron-exclusive positions",
             "POOL_11 with position 614 (D614G) removed",
             "POOL_12 with position 614 removed"),
  check.names = FALSE
)

kbl(pools_df,
    caption   = "Empirical position pools, computed from sarscov2_variants$Mutation_Sites.",
    align     = "llr",
    col.names = c("Pool", "Size", "Notes")) |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width        = FALSE)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Empirical position pools, computed from
sarscov2_variants\$Mutation_Sites.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Pool
</th>

<th style="text-align:left;">

Size
</th>

<th style="text-align:right;">

Notes
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

\|POOL_11\|
</td>

<td style="text-align:left;">

53
</td>

<td style="text-align:right;">

union of mutation sites across the 11 non-Omicron variants
</td>

</tr>

<tr>

<td style="text-align:left;">

\|POOL_12\|
</td>

<td style="text-align:left;">

75
</td>

<td style="text-align:right;">

union across all 12 variants
</td>

</tr>

<tr>

<td style="text-align:left;">

\|POOL_12 \\ POOL_11\|
</td>

<td style="text-align:left;">

22
</td>

<td style="text-align:right;">

Omicron-exclusive positions
</td>

</tr>

<tr>

<td style="text-align:left;">

\|POOL_11 \\ {614}\| (random-draw pool, non-Omicron variants)
</td>

<td style="text-align:left;">

52
</td>

<td style="text-align:right;">

POOL_11 with position 614 (D614G) removed
</td>

</tr>

<tr>

<td style="text-align:left;">

\|POOL_12 \\ {614}\| (random-draw pool, Sc 4 emerging)
</td>

<td style="text-align:left;">

74
</td>

<td style="text-align:right;">

POOL_12 with position 614 removed
</td>

</tr>

</tbody>

</table>

Position 614 (D614G) is present in both pools because every catalogued
variant carries this mutation. The simulator treats D614G as a
deterministic baseline mutation in every variant of every scenario; each
variant’s remaining $n_{\text{muts}} - 1$ positions are drawn from the
appropriate `POOL_*_DRAW` pool.

## Cross-variant collision properties

``` r
all_sites_flat     <- unlist(mut_sites)
site_variant_count <- table(all_sites_flat)
max_v_per_site     <- max(site_variant_count)
site_with_max      <- names(site_variant_count)[site_variant_count == max_v_per_site]

all_muts_flat       <- unlist(mut_strings)
mut_variant_count   <- table(all_muts_flat)
non_d614g_counts    <- mut_variant_count[names(mut_variant_count) != "D614G"]
second_max_exact    <- max(non_d614g_counts)
muts_second_max     <- names(non_d614g_counts)[non_d614g_counts == second_max_exact]

mut_per_site        <- tapply(all_muts_flat,
                              gsub("[^0-9]", "", all_muts_flat),
                              function(m) length(unique(m)))
max_muts_per_site   <- max(mut_per_site)
sites_max_diversity <- names(mut_per_site)[mut_per_site == max_muts_per_site]

# Distinct mutations observed at the position(s) of maximum diversity.
distinct_muts_at_max <- sort(unique(all_muts_flat[
  gsub("[^0-9]", "", all_muts_flat) %in% sites_max_diversity
]))

collisions_df <- data.frame(
  Property = c(
    "Max distinct mutations at one position",
    "Position(s) with that maximum",
    "Distinct mutations at that position",
    "Max variants sharing one exact mutation, including D614G",
    "Mutation(s) at that maximum",
    "Max variants sharing one exact mutation, excluding D614G",
    "Mutation(s) at that maximum"),
  Value    = c(
    as.character(max_muts_per_site),
    paste(sites_max_diversity, collapse = ", "),
    paste(distinct_muts_at_max, collapse = ", "),
    as.character(max_v_per_site),
    paste(names(mut_variant_count)[mut_variant_count == max_v_per_site],
          collapse = ", "),
    as.character(second_max_exact),
    paste(muts_second_max, collapse = ", "))
)

kbl(collisions_df,
    caption = "Cross-variant collision statistics, computed from sarscov2_variants.",
    align   = "ll") |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width        = FALSE)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Cross-variant collision statistics, computed from sarscov2_variants.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Property
</th>

<th style="text-align:left;">

Value
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Max distinct mutations at one position
</td>

<td style="text-align:left;">

3
</td>

</tr>

<tr>

<td style="text-align:left;">

Position(s) with that maximum
</td>

<td style="text-align:left;">

484
</td>

</tr>

<tr>

<td style="text-align:left;">

Distinct mutations at that position
</td>

<td style="text-align:left;">

E484A, E484K, E484Q
</td>

</tr>

<tr>

<td style="text-align:left;">

Max variants sharing one exact mutation, including D614G
</td>

<td style="text-align:left;">

12
</td>

</tr>

<tr>

<td style="text-align:left;">

Mutation(s) at that maximum
</td>

<td style="text-align:left;">

D614G
</td>

</tr>

<tr>

<td style="text-align:left;">

Max variants sharing one exact mutation, excluding D614G
</td>

<td style="text-align:left;">

7
</td>

</tr>

<tr>

<td style="text-align:left;">

Mutation(s) at that maximum
</td>

<td style="text-align:left;">

E484K
</td>

</tr>

</tbody>

</table>

D614G’s 12-variant prevalence reflects inheritance rather than
independent convergent evolution. The next-most-shared exact mutation
appears in 7 variants and represents the biologically meaningful upper
bound on independent convergent evolution within the catalogue. With
three established variants plus an emerging variant all drawing from
`POOL_11_DRAW` (size 52), position collisions are expected and
biologically realistic.

The substitutions E484K mutation, E484Q mutation, and E484A mutation
occur at position 484 within the receptor-binding domain (RBD) of the
SARS-CoV-2 spike protein, a region critically involved in antibody
recognition and ACE2 interaction. Among these, E484K has received
particular attention because it is strongly associated with immune
escape and reduced neutralization by monoclonal antibodies, convalescent
sera, and vaccine-elicited antibodies ([Yang et al.
2022](#ref-Yang2022)). The mutation emerged independently in several
major variants including Beta (B.1.351) and Gamma (P.1), suggesting
convergent evolution under immune selective pressure. Multiple studies
have identified E484K as one of the most important antibody-escape
substitutions in SARS-CoV-2 evolution, substantially altering
antigenicity while maintaining viral fitness. E484Q, observed for
example in the Kappa lineage, appears to confer more moderate immune
evasion effects, whereas E484A became prominent in Omicron and
contributes to the broader immune escape profile of that lineage.

# Deleterious-mutation noise

## Mechanism

The simulator injects independent low-frequency substitution noise into
the population matrix after variant substitution and before the entropy
computation. Each cell of the matrix (one sequence at one position) is
independently overwritten with probability $p_{\text{del}}$ by a
substitute residue drawn uniformly from the 19 non-self standard amino
acids. Noise is applied to every population in the matrix (reference,
$V_1$, $V_2$, $V_3$, and emerging) and at every position (1 through
1,273). Noise can occasionally overwrite a variant-defining
substitution; this is the natural reading of “deleterious mutations may
occur at any position” and is biologically realistic at the chosen rate.

Mathematically, after deterministic variant substitution the matrix $M$
of size $n_{\text{total}} \times 1273$ is transformed by $$
M'_{ij} = \begin{cases}
\text{Uniform}\big(\{1, \ldots, 20\} \setminus \{M_{ij}\}\big) & \text{w.p. } p_{\text{del}} \\
M_{ij} & \text{w.p. } 1 - p_{\text{del}}
\end{cases}
$$ independently across all $(i, j)$ pairs. The expected number of noise
events per matrix is
$n_{\text{total}} \times 1273 \times p_{\text{del}}$.

## Empirical justification of $p_{\text{del}} = 10^{-3}$

Mutations in SARS-CoV-2 proteins are predominantly neutral or
deleterious, with strongly deleterious substitutions removed rapidly
from circulating populations by purifying selection ([Bloom and Neher
2023](#ref-bloomneher2023)).

A useful theoretical anchor comes from CirSeq-based direct measurement
of the SARS-CoV-2 mutation rate. Symons et al. ([2025](#ref-Symons2025))
report a genome-wide substitution rate of approximately
$1.5 \times 10^{-6}$ per nucleotide per viral replication cycle,
dominated by C→U transitions. For the 1,273-amino-acid Spike protein
(3,819 nucleotides) this rate implies:

- Per codon per replication cycle:
  $3 \times 1.5 \times 10^{-6} = 4.5 \times 10^{-6}$ (an upper bound on
  amino-acid-level mutations since some codon substitutions are
  synonymous).
- Per Spike sequence per replication cycle:
  $3{,}819 \times 1.5 \times 10^{-6}
  \approx 5.7 \times 10^{-3}$ mutations.
- Per whole genome per replication cycle:
  $29{,}903 \times 1.5 \times 10^{-6}
  \approx 0.045$ mutations.

Low-frequency mutations are commonly interpreted as transient variants
under purifying selection, since strongly deleterious substitutions
rarely persist at appreciable prevalence in population-level
surveillance datasets. A commonly used threshold is $0.5\%$
([Tonkin-Hill et al. 2021](#ref-TonkinHill2021)). Lower thresholds may
be appropriate in high-depth datasets when the objective is early
detection of emerging mutational structure rather than definitive
clinical variant calling ([Van Poelvoorde et al.
2021](#ref-VanPoelvoorde2021)). In the present work, a $0.1\%$ threshold
was selected to increase the sensitivity of the entropy-based framework
to weak transient variation that may contribute cumulatively to
genome-wide information-theoretic signals.

To anchor $p_{\text{del}}$ empirically, we compute the per-site noise
rate using `empirical_pdel()` (see `empirical_pdel.R` in this folder) on
the full post-processed surveillance feature matrices used throughout
the rest of the pipeline. The NCBI matrix contains 109,536 sequences
spanning January 2020 to September 2021; the GISAID matrix contains
129,371 sequences spanning January 2020 to March 2024. The GISAID matrix
is the aligned variant of the same US Spike collection, so its 1,273
columns also carry gap characters wherever the aligned consensus
contains a deletion; gap codes are tabulated alongside ordinary
amino-acid codes, so positions at which gaps reach or exceed the 0.1%
threshold report “Gap” as one of the residues above threshold and the
remaining sub-threshold residues, including any low-frequency gaps,
contribute to the noise count. The NCBI matrix is unaligned and carries
only standard amino-acid codes at each of the 1,273 positions.

At each evaluated position we classify residues by their cumulative
observed proportion: residues at or above $0.1\%$ are treated as
biologically established (variant-defining substitutions, deletions, or
majority consensus), and residues below $0.1\%$ are treated as transient
deleterious or sequencing-artifact substitutions that did not fix in any
circulating lineage. The threshold cleanly separates the bimodal
distribution of residue frequencies in both matrices: variant residues,
even those carried by a single WHO-labelled variant, appear in thousands
of sequences across the catalogue, whereas truly transient residues stay
well below their corresponding $0.1\%$ count on each dataset.

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

GISAID dataset: empirical noise rates at the 75 variant-defining Spike
positions (POOL_12), computed by empirical_pdel() on the post-processed
GISAID US Spike feature matrix. Residues listed above threshold are
those whose observed cumulative frequency at the position is at or above
0.1%; the noise count is the sum of all sub-threshold residues. The
per-site rate is the noise count divided by total observations. Gap
characters from the alignment are labelled ‘Gap’ rather than ‘-’ for
clarity.
</caption>

<thead>

<tr>

<th style="text-align:right;">

Site
</th>

<th style="text-align:right;">

Total
</th>

<th style="text-align:left;">

Residues above threshold
</th>

<th style="text-align:right;">

Noise count
</th>

<th style="text-align:right;">

Per-site noise rate
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

5
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

L+F
</td>

<td style="text-align:right;">

48
</td>

<td style="text-align:right;">

3.71e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

13
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

I+S
</td>

<td style="text-align:right;">

86
</td>

<td style="text-align:right;">

6.65e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

18
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+L+F
</td>

<td style="text-align:right;">

98
</td>

<td style="text-align:right;">

7.58e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

19
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+I+L+T
</td>

<td style="text-align:right;">

119
</td>

<td style="text-align:right;">

9.20e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

20
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+I+T
</td>

<td style="text-align:right;">

107
</td>

<td style="text-align:right;">

8.27e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

26
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

L+P+S+T+Gap
</td>

<td style="text-align:right;">

198
</td>

<td style="text-align:right;">

1.53e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

67
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

A+S+V
</td>

<td style="text-align:right;">

90
</td>

<td style="text-align:right;">

6.96e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

69
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

H+S+Y+Gap
</td>

<td style="text-align:right;">

227
</td>

<td style="text-align:right;">

1.75e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

70
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

G+I+F+V+Gap
</td>

<td style="text-align:right;">

189
</td>

<td style="text-align:right;">

1.46e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

80
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

A+D+G+P+Y
</td>

<td style="text-align:right;">

226
</td>

<td style="text-align:right;">

1.75e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

95
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

I+K+S+T
</td>

<td style="text-align:right;">

138
</td>

<td style="text-align:right;">

1.07e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

138
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

D+C+H+F+Y
</td>

<td style="text-align:right;">

228
</td>

<td style="text-align:right;">

1.76e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

142
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

D+G+F+Y+Gap
</td>

<td style="text-align:right;">

219
</td>

<td style="text-align:right;">

1.69e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

144
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

H+K+Y+Gap
</td>

<td style="text-align:right;">

362
</td>

<td style="text-align:right;">

2.80e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

145
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+H+K+Y+V+Gap
</td>

<td style="text-align:right;">

103
</td>

<td style="text-align:right;">

7.96e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

152
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+C+E+L+K+S+W
</td>

<td style="text-align:right;">

32
</td>

<td style="text-align:right;">

2.47e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

154
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

E+F+W
</td>

<td style="text-align:right;">

251
</td>

<td style="text-align:right;">

1.94e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

156
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+E+G+V+Gap
</td>

<td style="text-align:right;">

146
</td>

<td style="text-align:right;">

1.13e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

157
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

L+F+S+Y+V+Gap
</td>

<td style="text-align:right;">

45
</td>

<td style="text-align:right;">

3.48e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

158
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+G+S+Y+Gap
</td>

<td style="text-align:right;">

43
</td>

<td style="text-align:right;">

3.32e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

190
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+F+S+V
</td>

<td style="text-align:right;">

239
</td>

<td style="text-align:right;">

1.85e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

211
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+N+Gap
</td>

<td style="text-align:right;">

232
</td>

<td style="text-align:right;">

1.79e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

212
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

E+I+L+S
</td>

<td style="text-align:right;">

273
</td>

<td style="text-align:right;">

2.11e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

215
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

D+G+H+Y
</td>

<td style="text-align:right;">

258
</td>

<td style="text-align:right;">

1.99e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

222
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

A+S+V
</td>

<td style="text-align:right;">

154
</td>

<td style="text-align:right;">

1.19e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

241
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

L+Gap
</td>

<td style="text-align:right;">

23
</td>

<td style="text-align:right;">

1.78e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

242
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

L+Gap
</td>

<td style="text-align:right;">

214
</td>

<td style="text-align:right;">

1.65e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

243
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

A+S+Gap
</td>

<td style="text-align:right;">

177
</td>

<td style="text-align:right;">

1.37e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

253
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

D+G
</td>

<td style="text-align:right;">

248
</td>

<td style="text-align:right;">

1.92e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

258
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

L+W
</td>

<td style="text-align:right;">

157
</td>

<td style="text-align:right;">

1.21e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

339
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

D+G+H
</td>

<td style="text-align:right;">

256
</td>

<td style="text-align:right;">

1.98e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

346
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+I+K+S+T
</td>

<td style="text-align:right;">

53
</td>

<td style="text-align:right;">

4.10e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

371
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

L+F+S
</td>

<td style="text-align:right;">

38
</td>

<td style="text-align:right;">

2.94e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

373
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

P+S
</td>

<td style="text-align:right;">

98
</td>

<td style="text-align:right;">

7.58e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

375
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

F+S
</td>

<td style="text-align:right;">

151
</td>

<td style="text-align:right;">

1.17e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

384
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

L+P
</td>

<td style="text-align:right;">

136
</td>

<td style="text-align:right;">

1.05e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

417
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+K+T
</td>

<td style="text-align:right;">

95
</td>

<td style="text-align:right;">

7.34e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

440
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+K
</td>

<td style="text-align:right;">

97
</td>

<td style="text-align:right;">

7.50e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

446
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

G+S+V
</td>

<td style="text-align:right;">

202
</td>

<td style="text-align:right;">

1.56e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

452
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+Q+L+M
</td>

<td style="text-align:right;">

143
</td>

<td style="text-align:right;">

1.11e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

477
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+I+S
</td>

<td style="text-align:right;">

172
</td>

<td style="text-align:right;">

1.33e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

478
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+I+K+T
</td>

<td style="text-align:right;">

210
</td>

<td style="text-align:right;">

1.62e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

484
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

A+R+Q+E+K+V
</td>

<td style="text-align:right;">

183
</td>

<td style="text-align:right;">

1.41e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

490
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

F+P+S
</td>

<td style="text-align:right;">

194
</td>

<td style="text-align:right;">

1.50e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

493
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+Q
</td>

<td style="text-align:right;">

197
</td>

<td style="text-align:right;">

1.52e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

494
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

P+S
</td>

<td style="text-align:right;">

153
</td>

<td style="text-align:right;">

1.18e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

496
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

G+S
</td>

<td style="text-align:right;">

43
</td>

<td style="text-align:right;">

3.32e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

498
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+Q
</td>

<td style="text-align:right;">

26
</td>

<td style="text-align:right;">

2.01e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

501
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+T+Y
</td>

<td style="text-align:right;">

74
</td>

<td style="text-align:right;">

5.72e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

505
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

H+Y
</td>

<td style="text-align:right;">

13
</td>

<td style="text-align:right;">

1.00e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

516
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

E
</td>

<td style="text-align:right;">

125
</td>

<td style="text-align:right;">

9.66e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

547
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

I+K+T
</td>

<td style="text-align:right;">

40
</td>

<td style="text-align:right;">

3.09e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

570
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

A+D+V
</td>

<td style="text-align:right;">

183
</td>

<td style="text-align:right;">

1.41e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

614
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

D+G
</td>

<td style="text-align:right;">

71
</td>

<td style="text-align:right;">

5.49e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

655
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

H+Y
</td>

<td style="text-align:right;">

52
</td>

<td style="text-align:right;">

4.02e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

677
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

Q+H+P
</td>

<td style="text-align:right;">

170
</td>

<td style="text-align:right;">

1.31e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

679
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+K
</td>

<td style="text-align:right;">

192
</td>

<td style="text-align:right;">

1.48e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

681
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+H+P
</td>

<td style="text-align:right;">

182
</td>

<td style="text-align:right;">

1.41e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

701
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

A+S+V
</td>

<td style="text-align:right;">

75
</td>

<td style="text-align:right;">

5.80e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

716
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

I+T
</td>

<td style="text-align:right;">

30
</td>

<td style="text-align:right;">

2.32e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

764
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+K
</td>

<td style="text-align:right;">

30
</td>

<td style="text-align:right;">

2.32e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

796
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

D+H+Y
</td>

<td style="text-align:right;">

59
</td>

<td style="text-align:right;">

4.56e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

856
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+K
</td>

<td style="text-align:right;">

143
</td>

<td style="text-align:right;">

1.11e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

859
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+I+T
</td>

<td style="text-align:right;">

93
</td>

<td style="text-align:right;">

7.19e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

888
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

L+F
</td>

<td style="text-align:right;">

13
</td>

<td style="text-align:right;">

1.00e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

950
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+D+H
</td>

<td style="text-align:right;">

54
</td>

<td style="text-align:right;">

4.17e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

954
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

Q+H
</td>

<td style="text-align:right;">

116
</td>

<td style="text-align:right;">

8.97e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

957
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

R+Q
</td>

<td style="text-align:right;">

47
</td>

<td style="text-align:right;">

3.63e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

969
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+K
</td>

<td style="text-align:right;">

27
</td>

<td style="text-align:right;">

2.09e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

981
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

L+F
</td>

<td style="text-align:right;">

34
</td>

<td style="text-align:right;">

2.63e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

982
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

A+S
</td>

<td style="text-align:right;">

27
</td>

<td style="text-align:right;">

2.09e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

1027
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

I+T
</td>

<td style="text-align:right;">

43
</td>

<td style="text-align:right;">

3.32e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

1071
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

Q
</td>

<td style="text-align:right;">

194
</td>

<td style="text-align:right;">

1.50e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

1118
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

D+H+Y
</td>

<td style="text-align:right;">

48
</td>

<td style="text-align:right;">

3.71e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

1191
</td>

<td style="text-align:right;">

129371
</td>

<td style="text-align:left;">

N+K
</td>

<td style="text-align:right;">

56
</td>

<td style="text-align:right;">

4.33e-04
</td>

</tr>

</tbody>

</table>

Across the evaluated sites the pooled rate (9.86e-04), median per-site
rate (9.20e-04), and maximum per-site rate (2.80e-03 at site 144)
jointly establish the empirical range that supports
$p_{\text{del}} = 10^{-3}$ as a conservative central default.

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

GISAID dataset: summary statistics of per-site noise rates across the
evaluated variant-defining positions (POOL_12).
</caption>

<thead>

<tr>

<th style="text-align:left;">

Statistic
</th>

<th style="text-align:right;">

Value
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Pooled rate (total noise / total observations)
</td>

<td style="text-align:right;">

9.86e-04
</td>

</tr>

<tr>

<td style="text-align:left;">

Mean per-site rate
</td>

<td style="text-align:right;">

9.86e-04
</td>

</tr>

<tr>

<td style="text-align:left;">

Median per-site rate
</td>

<td style="text-align:right;">

9.20e-04
</td>

</tr>

<tr>

<td style="text-align:left;">

Minimum per-site rate (site 505)
</td>

<td style="text-align:right;">

1.00e-04
</td>

</tr>

<tr>

<td style="text-align:left;">

Maximum per-site rate (site 144)
</td>

<td style="text-align:right;">

2.80e-03
</td>

</tr>

</tbody>

</table>

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

NCBI dataset: empirical noise rates at the 53 non-Omicron
variant-defining Spike positions (POOL_11), computed by empirical_pdel()
on the post-processed NCBI US Spike feature matrix. Residues listed
above threshold are those whose observed cumulative frequency at the
position is at or above 0.1%; the noise count is the sum of all
sub-threshold residues.
</caption>

<thead>

<tr>

<th style="text-align:right;">

Site
</th>

<th style="text-align:right;">

Total
</th>

<th style="text-align:left;">

Residues above threshold
</th>

<th style="text-align:right;">

Noise count
</th>

<th style="text-align:right;">

Per-site noise rate
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:right;">

5
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

L+F
</td>

<td style="text-align:right;">

5
</td>

<td style="text-align:right;">

4.56e-05
</td>

</tr>

<tr>

<td style="text-align:right;">

13
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

I+S
</td>

<td style="text-align:right;">

14
</td>

<td style="text-align:right;">

1.28e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

18
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

L+F
</td>

<td style="text-align:right;">

23
</td>

<td style="text-align:right;">

2.10e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

19
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

R+T
</td>

<td style="text-align:right;">

72
</td>

<td style="text-align:right;">

6.57e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

20
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

N+I+T
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

3.65e-05
</td>

</tr>

<tr>

<td style="text-align:right;">

26
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

L+P+S
</td>

<td style="text-align:right;">

63
</td>

<td style="text-align:right;">

5.75e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

67
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

A+V
</td>

<td style="text-align:right;">

69
</td>

<td style="text-align:right;">

6.30e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

69
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

H
</td>

<td style="text-align:right;">

110
</td>

<td style="text-align:right;">

1.00e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

70
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

V
</td>

<td style="text-align:right;">

95
</td>

<td style="text-align:right;">

8.67e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

80
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

D+Y
</td>

<td style="text-align:right;">

72
</td>

<td style="text-align:right;">

6.57e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

95
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

I+T
</td>

<td style="text-align:right;">

68
</td>

<td style="text-align:right;">

6.21e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

138
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

D+H+Y
</td>

<td style="text-align:right;">

19
</td>

<td style="text-align:right;">

1.73e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

142
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

D+G+S+Y
</td>

<td style="text-align:right;">

17
</td>

<td style="text-align:right;">

1.55e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

144
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

K+Y
</td>

<td style="text-align:right;">

53
</td>

<td style="text-align:right;">

4.84e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

152
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

C+L+S+W
</td>

<td style="text-align:right;">

102
</td>

<td style="text-align:right;">

9.31e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

154
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

A+E+F
</td>

<td style="text-align:right;">

129
</td>

<td style="text-align:right;">

1.18e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

156
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

E+V
</td>

<td style="text-align:right;">

27
</td>

<td style="text-align:right;">

2.46e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

157
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

C+F+Y
</td>

<td style="text-align:right;">

162
</td>

<td style="text-align:right;">

1.48e-03
</td>

</tr>

<tr>

<td style="text-align:right;">

158
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

R+G+S
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

3.65e-05
</td>

</tr>

<tr>

<td style="text-align:right;">

190
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

R+S+V
</td>

<td style="text-align:right;">

40
</td>

<td style="text-align:right;">

3.65e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

215
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

D+Y
</td>

<td style="text-align:right;">

72
</td>

<td style="text-align:right;">

6.57e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

222
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

A+V
</td>

<td style="text-align:right;">

13
</td>

<td style="text-align:right;">

1.19e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

241
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

L
</td>

<td style="text-align:right;">

3
</td>

<td style="text-align:right;">

2.74e-05
</td>

</tr>

<tr>

<td style="text-align:right;">

242
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

L
</td>

<td style="text-align:right;">

18
</td>

<td style="text-align:right;">

1.64e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

243
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

A
</td>

<td style="text-align:right;">

70
</td>

<td style="text-align:right;">

6.39e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

253
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

D+G
</td>

<td style="text-align:right;">

24
</td>

<td style="text-align:right;">

2.19e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

258
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

L+W
</td>

<td style="text-align:right;">

19
</td>

<td style="text-align:right;">

1.73e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

384
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

P
</td>

<td style="text-align:right;">

101
</td>

<td style="text-align:right;">

9.22e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

417
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

K+T
</td>

<td style="text-align:right;">

9
</td>

<td style="text-align:right;">

8.22e-05
</td>

</tr>

<tr>

<td style="text-align:right;">

452
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

R+L
</td>

<td style="text-align:right;">

21
</td>

<td style="text-align:right;">

1.92e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

477
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

N+S
</td>

<td style="text-align:right;">

104
</td>

<td style="text-align:right;">

9.49e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

478
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

K+T
</td>

<td style="text-align:right;">

64
</td>

<td style="text-align:right;">

5.84e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

484
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

Q+E+K
</td>

<td style="text-align:right;">

9
</td>

<td style="text-align:right;">

8.22e-05
</td>

</tr>

<tr>

<td style="text-align:right;">

490
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

F
</td>

<td style="text-align:right;">

56
</td>

<td style="text-align:right;">

5.11e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

494
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

P+S
</td>

<td style="text-align:right;">

6
</td>

<td style="text-align:right;">

5.48e-05
</td>

</tr>

<tr>

<td style="text-align:right;">

501
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

N+T+Y
</td>

<td style="text-align:right;">

9
</td>

<td style="text-align:right;">

8.22e-05
</td>

</tr>

<tr>

<td style="text-align:right;">

516
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

E
</td>

<td style="text-align:right;">

17
</td>

<td style="text-align:right;">

1.55e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

570
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

A+D
</td>

<td style="text-align:right;">

93
</td>

<td style="text-align:right;">

8.49e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

614
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

D+G
</td>

<td style="text-align:right;">

14
</td>

<td style="text-align:right;">

1.28e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

655
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

H+Y
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

0.00e+00
</td>

</tr>

<tr>

<td style="text-align:right;">

677
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

Q+H+P
</td>

<td style="text-align:right;">

67
</td>

<td style="text-align:right;">

6.12e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

681
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

R+H+P
</td>

<td style="text-align:right;">

100
</td>

<td style="text-align:right;">

9.13e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

701
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

A+V
</td>

<td style="text-align:right;">

25
</td>

<td style="text-align:right;">

2.28e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

716
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

I+T
</td>

<td style="text-align:right;">

0
</td>

<td style="text-align:right;">

0.00e+00
</td>

</tr>

<tr>

<td style="text-align:right;">

859
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

I+T
</td>

<td style="text-align:right;">

32
</td>

<td style="text-align:right;">

2.92e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

888
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

F
</td>

<td style="text-align:right;">

7
</td>

<td style="text-align:right;">

6.39e-05
</td>

</tr>

<tr>

<td style="text-align:right;">

950
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

N+D
</td>

<td style="text-align:right;">

26
</td>

<td style="text-align:right;">

2.37e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

957
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

R+Q
</td>

<td style="text-align:right;">

18
</td>

<td style="text-align:right;">

1.64e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

982
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

A+S
</td>

<td style="text-align:right;">

21
</td>

<td style="text-align:right;">

1.92e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

1027
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

I+T
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

9.13e-06
</td>

</tr>

<tr>

<td style="text-align:right;">

1071
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

Q+L
</td>

<td style="text-align:right;">

103
</td>

<td style="text-align:right;">

9.40e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

1118
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

D+H
</td>

<td style="text-align:right;">

90
</td>

<td style="text-align:right;">

8.22e-04
</td>

</tr>

<tr>

<td style="text-align:right;">

1191
</td>

<td style="text-align:right;">

109536
</td>

<td style="text-align:left;">

N+K
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

3.65e-05
</td>

</tr>

</tbody>

</table>

For the NCBI dataset at these non-Omicron-defining sites the pooled rate
is 4.07e-04, the median per-site rate is 2.28e-04, and the maximum is
1.48e-03 at site 157.

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

NCBI dataset: summary statistics of per-site noise rates across the
evaluated non-Omicron variant-defining positions (POOL_11).
</caption>

<thead>

<tr>

<th style="text-align:left;">

Statistic
</th>

<th style="text-align:right;">

Value
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Pooled rate (total noise / total observations)
</td>

<td style="text-align:right;">

4.07e-04
</td>

</tr>

<tr>

<td style="text-align:left;">

Mean per-site rate
</td>

<td style="text-align:right;">

4.07e-04
</td>

</tr>

<tr>

<td style="text-align:left;">

Median per-site rate
</td>

<td style="text-align:right;">

2.28e-04
</td>

</tr>

<tr>

<td style="text-align:left;">

Minimum per-site rate (site 655)
</td>

<td style="text-align:right;">

0.00e+00
</td>

</tr>

<tr>

<td style="text-align:left;">

Maximum per-site rate (site 157)
</td>

<td style="text-align:right;">

1.48e-03
</td>

</tr>

</tbody>

</table>

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Genome-wide (all 1,273 Spike positions) empirical noise-rate summary on
GISAID and NCBI feature matrices at threshold = 0.1%. The pooled rate
equals the mean per-site rate by construction whenever every position
has identical total observations; the mean row is omitted to avoid
visual redundancy.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Statistic
</th>

<th style="text-align:right;">

GISAID
</th>

<th style="text-align:right;">

NCBI
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Pooled rate (= mean per-site rate)
</td>

<td style="text-align:right;">

6.85e-04
</td>

<td style="text-align:right;">

1.48e-04
</td>

</tr>

<tr>

<td style="text-align:left;">

Median per-site rate
</td>

<td style="text-align:right;">

5.18e-04
</td>

<td style="text-align:right;">

4.56e-05
</td>

</tr>

<tr>

<td style="text-align:left;">

Minimum per-site rate
</td>

<td style="text-align:right;">

0.00e+00 (site 1256)
</td>

<td style="text-align:right;">

0.00e+00 (site 1)
</td>

</tr>

<tr>

<td style="text-align:left;">

Maximum per-site rate
</td>

<td style="text-align:right;">

3.00e-03 (site 218)
</td>

<td style="text-align:right;">

1.55e-03 (site 1111)
</td>

</tr>

</tbody>

</table>

Interestingly, overall maximum occurs at spike residue 218 located
within the N-terminal domain (NTD) of the SARS-CoV-2 spike protein, a
region known to undergo substantial antigenic evolution and recurrent
mutation across multiple variants. Although position 218 itself is not
recognized as a major canonical immune-escape site, elevated variability
in this region may reflect broader selective pressures acting on
NTD-associated antigenic structure ([Tao et al. 2021](#ref-Tao2021)).

The genome-wide rate is consistently lower than the variant-defining-
site rate on the same dataset, by roughly 30%. This is the expected
sign: positions under positive selection accumulate slightly more
low-frequency variation per replication-cycle window than the genome
average ([Desai, Walczak, and Fisher 2013](#ref-Desai2013)). Both rates
lie below the chosen default $p_{\text{del}} =
10^{-3}$, and the maximum per-site rate stays at or below
$3 \times 10^{-3}$ across both datasets. The default is therefore above
the typical site (median, mean, pooled) and below the worst single-site
value on each dataset, providing a conservative central choice
consistent with all four empirical anchors evaluated above (GISAID-VOC,
NCBI-VOC, GISAID-overall, NCBI-overall).

The simulator uses $p_{\text{del}} = 10^{-3}$ as its default
deleterious-mutation rate, rounding the empirical estimate upward to a
conservative round value. The chosen value is consistent with each of
the four empirical anchors: no anchor site exceeds $2 \times 10^{-3}$ at
the VOC scale, and on the genome-wide scale the median per-site rate
lies below $10^{-3}$ on both datasets.

The deleterious-noise rate is treated as a surveillance-level prevalence
model for transient low-frequency variation rather than a direct
biological mutation-rate estimate. The rate captures the joint effect of
de-novo deleterious mutations and sequencing-artifact substitutions that
did not fix in any circulating lineage. A robustness sweep over
$p_{\text{del}} \in \{0, 10^{-4}, 10^{-3}, 10^{-2}\}$ on a
stratified-random subset of cells quantifies sensitivity of the inferred
thresholds to the noise rate.

# Scenarios

The reference sequence used by the simulator is the canonical pre-D614G
1,273-aa Spike (GenBank accession YP_009724390). Every variant in every
scenario carries D614G as its baseline mutation, modelling the
post-mid-2020 surveillance regime in which every detected variant
inherits D614G ([Korber et al. 2020](#ref-korber2020)).

Three labels distinguish how parameters vary across the study:

- **Fixed**: identical in every replicate of every cell.
- **Per-cell**: identical in all 30 replicates of one cell, varies
  across cells.
- **Per-replicate**: redrawn for each of the 30 replicates of a cell.

Note: the wild-type-only scenario (reference + emerging only) is omitted
from this design. Under that scenario the emerging variant is the sole
source of variability in the population, every mutation site shares an
identical entropy structure, and the entropy-GMM step is trivially
successful at $n_{\text{emerge}} = 2$. The omitted scenario contributes
no methodologically interesting information about the pipeline’s
sample-size requirements.

## Scenario 1: detect $V_2$ on a D614G-dominant background

| Population | Composition | Sample size |
|----|----|----|
| Reference | Canonical pre-D614G 1,273-aa Spike (fixed) | $N_{\text{ref}}$ per-cell |
| $V_1$ | D614G only (fixed) | $n_{V_1}$ per-cell |
| Emerging | D614G plus $(n_{\text{muts}} - 1)$ positions from `POOL_11 \ {614}`, | swept |
|  | $n_{\text{muts}} \in$ `MUT_COUNT_SET_11` (per-replicate) |  |

## Scenario 2: detect $V_3$ on a $V_1 + V_2$ background

| Population | Composition | Sample size |
|----|----|----|
| Reference | Canonical pre-D614G 1,273-aa Spike (fixed) | $N_{\text{ref}}$ per-cell |
| $V_1$ | D614G only (fixed) | $n_{V_1}$ per-cell |
| $V_2$ | D614G plus $(n_{V_2,\text{muts}} - 1)$ from `POOL_11 \ {614}`, | $n_{V_2}$ per-cell |
|  | $n_{V_2,\text{muts}} \in$ `MUT_COUNT_SET_11` (per-replicate) |  |
| Emerging | D614G plus $(n_{\text{muts}} - 1)$ from `POOL_11 \ {614}`, | swept |
|  | $n_{\text{muts}} \in$ `MUT_COUNT_SET_11` (per-replicate) |  |

## Scenario 3: detect $V_4$ (non-Omicron-like) on a $V_1 + V_2 + V_3$ background

| Population | Composition | Sample size |
|----|----|----|
| Reference | Canonical pre-D614G 1,273-aa Spike (fixed) | $N_{\text{ref}}$ per-cell |
| $V_1$ | D614G only (fixed) | $n_{V_1}$ per-cell |
| $V_2$ | D614G plus $(n_{V_2,\text{muts}} - 1)$ from `POOL_11 \ {614}`, | $n_{V_2}$ per-cell |
|  | $n_{V_2,\text{muts}} \in$ `MUT_COUNT_SET_11` (per-replicate) |  |
| $V_3$ | D614G plus $(n_{V_3,\text{muts}} - 1)$ from `POOL_11 \ {614}`, | $n_{V_3}$ per-cell |
|  | $n_{V_3,\text{muts}} \in$ `MUT_COUNT_SET_11` (per-replicate) |  |
| Emerging | D614G plus $(n_{\text{muts}} - 1)$ from `POOL_11 \ {614}`, | swept |
|  | $n_{\text{muts}} \in$ `MUT_COUNT_SET_11` (per-replicate) |  |

## Scenario 4: detect $V_4$ (Omicron-like) on a $V_1 + V_2 + V_3$ background

| Population | Composition | Sample size |
|----|----|----|
| Reference | Canonical pre-D614G 1,273-aa Spike (fixed) | $N_{\text{ref}}$ per-cell |
| $V_1$ | D614G only (fixed) | $n_{V_1}$ per-cell |
| $V_2$ | D614G plus $(n_{V_2,\text{muts}} - 1)$ from `POOL_11 \ {614}`, | $n_{V_2}$ per-cell |
|  | $n_{V_2,\text{muts}} \in$ `MUT_COUNT_SET_11` (per-replicate) |  |
| $V_3$ | D614G plus $(n_{V_3,\text{muts}} - 1)$ from `POOL_11 \ {614}`, | $n_{V_3}$ per-cell |
|  | $n_{V_3,\text{muts}} \in$ `MUT_COUNT_SET_11` (per-replicate) |  |
| Emerging | D614G plus 32 positions from `POOL_12 \ {614}` (fixed total of 33) | swept |

# Sampling structure

The design is a full discrete factorial across structural axes with 30
replicates per cell.

``` r
N_REF_GRID_SIZE   <- 10L
RATIOS            <- c(1.25, 1.5, 1.75, 2.5, 5)
N_REPS_PER_CELL   <- 30L
BASE_SEED         <- 2025L

BAND_RANGES <- list(
  small  = c(10L,   100L),
  medium = c(100L,  1000L),
  large  = c(1000L, 10000L)
)

n_ref_grid <- function(band) {
  rng <- BAND_RANGES[[band]]
  unique(round(exp(seq(log(rng[1]), log(rng[2]),
                       length.out = N_REF_GRID_SIZE))))
}

n_ref_grids_df <- data.frame(
  Band  = names(BAND_RANGES),
  Range = sapply(BAND_RANGES,
                 function(r) sprintf("[%d, %d]", r[1], r[2])),
  Grid  = sapply(names(BAND_RANGES),
                 function(b) paste(n_ref_grid(b), collapse = ", "))
)

kbl(n_ref_grids_df,
    caption = "Log-spaced N_ref grid per band.",
    align   = "lll") |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width        = FALSE)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Log-spaced N_ref grid per band.
</caption>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:left;">

Band
</th>

<th style="text-align:left;">

Range
</th>

<th style="text-align:left;">

Grid
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

small
</td>

<td style="text-align:left;">

small
</td>

<td style="text-align:left;">

\[10, 100\]
</td>

<td style="text-align:left;">

10, 13, 17, 22, 28, 36, 46, 60, 77, 100
</td>

</tr>

<tr>

<td style="text-align:left;">

medium
</td>

<td style="text-align:left;">

medium
</td>

<td style="text-align:left;">

\[100, 1000\]
</td>

<td style="text-align:left;">

100, 129, 167, 215, 278, 359, 464, 599, 774, 1000
</td>

</tr>

<tr>

<td style="text-align:left;">

large
</td>

<td style="text-align:left;">

large
</td>

<td style="text-align:left;">

\[1000, 10000\]
</td>

<td style="text-align:left;">

1000, 1292, 1668, 2154, 2783, 3594, 4642, 5995, 7743, 10000
</td>

</tr>

</tbody>

</table>

## Axis 1: band

3 categorical values describe the surveillance regime for the wild-type
reference background: small (early-pandemic or sparse sampling), medium
(established regional surveillance), and large (mature national or
archival surveillance). The four-order-of- magnitude span from 10 to
10,000 sequences matches empirical United States surveillance volumes
documented in Sayers et al. ([2022](#ref-sayers2022)) (NCBI) and Shu and
McCauley ([2017](#ref-shu2017)) (GISAID).

## Axis 2: $N_{\text{ref}}$

10 log-spaced integer points within each band, as shown in the table
above. Log-spaced points weight the lower (harder, smaller-sample) end
of each band as heavily as the upper end.

## Axes 3, 4, and 5: ratios

5 fixed multipliers over $N_{\text{ref}}$: {1.25, 1.5, 1.75, 2.5, 5}.

$$n_{V_i} = \text{round}(N_{\text{ref}} \times \text{ratio}_{V_i})$$

with no upper cap. The ratio reparameterisation enforces
$n_{V_i} > N_{\text{ref}}$ by construction because all multipliers are
greater than 1. Scenario 1 uses one ratio axis ($V_1$), scenario 2 uses
two ($V_1$, $V_2$), and scenarios 3 and 4 each use three ($V_1$, $V_2$,
$V_3$).

### Bands describe $N_{\text{ref}}$ scale only

A band is not a hard interval constraining the entire population. It
describes the surveillance volume of the wild-type reference background
alone. Co-circulating dominant variants typically achieve and maintain
sequence abundance well in excess of the prior reference-baseline
surveillance volume during their period of dominance ([Harvey et al.
2021](#ref-harvey2021)). Capping $n_{V_i}$ at the band’s upper bound
would systematically suppress the dominance signal.

## Cells and simulations per scenario

``` r
n_ratios   <- length(RATIOS)
N_REF_PTS  <- sum(sapply(names(BAND_RANGES),
                         function(b) length(n_ref_grid(b))))

cells_per_scenario <- c(
  Sc1 = N_REF_PTS * n_ratios,
  Sc2 = N_REF_PTS * n_ratios^2,
  Sc3 = N_REF_PTS * n_ratios^3,
  Sc4 = N_REF_PTS * n_ratios^3
)

sims_per_scenario <- cells_per_scenario * N_REPS_PER_CELL
total_cells       <- sum(cells_per_scenario)
total_sims        <- sum(sims_per_scenario)

cell_df <- data.frame(
  Scenario     = names(cells_per_scenario),
  `Axes used`  = c("band x N_ref x ratio_V1",
                   "band x N_ref x ratio_V1 x ratio_V2",
                   "band x N_ref x ratio_V1 x ratio_V2 x ratio_V3",
                   "band x N_ref x ratio_V1 x ratio_V2 x ratio_V3"),
  Cells        = cells_per_scenario,
  Replicates   = rep(N_REPS_PER_CELL, length(cells_per_scenario)),
  Simulations  = sims_per_scenario,
  check.names  = FALSE
)

kbl(cell_df,
    caption     = sprintf(
      "Cell counts per scenario. Total: %s cells, %s simulations.",
      format(total_cells, big.mark = ","),
      format(total_sims,  big.mark = ",")),
    align       = "llrrr",
    format.args = list(big.mark = ",")) |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width        = FALSE)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Cell counts per scenario. Total: 8,400 cells, 252,000 simulations.
</caption>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:left;">

Scenario
</th>

<th style="text-align:left;">

Axes used
</th>

<th style="text-align:right;">

Cells
</th>

<th style="text-align:right;">

Replicates
</th>

<th style="text-align:right;">

Simulations
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Sc1
</td>

<td style="text-align:left;">

Sc1
</td>

<td style="text-align:left;">

band x N_ref x ratio_V1
</td>

<td style="text-align:right;">

150
</td>

<td style="text-align:right;">

30
</td>

<td style="text-align:right;">

4,500
</td>

</tr>

<tr>

<td style="text-align:left;">

Sc2
</td>

<td style="text-align:left;">

Sc2
</td>

<td style="text-align:left;">

band x N_ref x ratio_V1 x ratio_V2
</td>

<td style="text-align:right;">

750
</td>

<td style="text-align:right;">

30
</td>

<td style="text-align:right;">

22,500
</td>

</tr>

<tr>

<td style="text-align:left;">

Sc3
</td>

<td style="text-align:left;">

Sc3
</td>

<td style="text-align:left;">

band x N_ref x ratio_V1 x ratio_V2 x ratio_V3
</td>

<td style="text-align:right;">

3,750
</td>

<td style="text-align:right;">

30
</td>

<td style="text-align:right;">

112,500
</td>

</tr>

<tr>

<td style="text-align:left;">

Sc4
</td>

<td style="text-align:left;">

Sc4
</td>

<td style="text-align:left;">

band x N_ref x ratio_V1 x ratio_V2 x ratio_V3
</td>

<td style="text-align:right;">

3,750
</td>

<td style="text-align:right;">

30
</td>

<td style="text-align:right;">

112,500
</td>

</tr>

</tbody>

</table>

## What is fixed and what is random within a cell

Fixed for all 30 replicates of a cell (bit-identical across reruns):

- `band`, $N_{\text{ref}}$, and the $n_{V_i}$ values derived from the
  cell’s ratios.
- Reference sequence: canonical pre-D614G 1,273-aa Spike.
- $V_1$ structure: D614G only, one mutation at position 614 with
  substitute G. Identical across all 30 replicates of all cells.

Per-replicate random within a cell (independent across the 30
replicates):

- $V_2$’s mutation count, drawn uniformly from `MUT_COUNT_SET_11`
  (scenarios 2 through 4).
- $V_2$’s non-614 positions: drawn without replacement from
  `POOL_11 \ {614}`. $V_2$’s substitute amino acids at non-614
  positions, uniform over the 19 non-reference standard residues,
  independent per position. (Scenarios 2 through 4.)
- $V_3$’s mutation count, positions, and substitutes: same scheme as
  $V_2$. (Scenarios 3 and 4.)
- Emerging variant’s mutation count, drawn uniformly from
  `MUT_COUNT_SET_11` (scenarios 1 through 3) or fixed at 33 (scenario
  4).
- Emerging variant’s positions: D614G plus draws from `POOL_11 \ {614}`
  (scenarios 1 through 3) or `POOL_12 \ {614}` (scenario 4). Substitute
  amino acids drawn uniformly per position.
- Deleterious-noise event pattern: independent per (sequence, position)
  cell of the matrix, Bernoulli at rate $p_{\text{del}}$, applied once
  at the cell’s ceiling matrix.
- The replicate’s RNG seed:
  `BASE_SEED + scenario * 10^8 + cell_id *   10^3 + rep_id`. The noise
  RNG uses `rep_seed + 1` as offset so that the same replicate at
  different $p_{\text{del}}$ values shares its variant placement but
  gets independent noise patterns. The seed is set inside the subprocess
  (or sequential iteration) before any per-replicate random draw.

## Why 30 replicates per cell

The 30-replicate design samples the within-cell distribution and
produces cell-mean confidence intervals of width
$\pm 1.96 \cdot \sigma / \sqrt{30}$, where $\sigma$ is the within-cell
standard deviation of $n_{\text{emerge}}^{\text{needed}}$ ([Cohen
1988](#ref-cohen1988)). Beyond approximately 30 replicates the precision
gains diminish as $1/\sqrt{n_{\text{reps}}}$, motivating 30 as a
conventional cost-precision tradeoff.

## Why full factorial rather than stratified random sampling

Both designs are defensible ([Cochran 1977](#ref-cochran1977); [Morris,
White, and Crowther 2019](#ref-morris2019)). Full factorial offers four
practical advantages: balanced evaluation across every grid cell (with
no random cell-selection bias to account for); direct marginal estimates
without weighted-average corrections; direct interaction estimates from
standard linear models on the balanced design; and clean visual
structure in heatmaps and faceted plots.

# Sample-size sweep procedure

For each replicate, the procedure sweeps $n_{\text{emerge}}$ from 2 to a
ceiling defined below using a 50-point log-spaced integer grid with no
refinement. Detection at grid point $g_k$ reports
$n_{\text{emerge}}^{\text{needed}} = g_k$ (grid-precision answer).

## Ceiling

The ceiling defines the maximum sample size of emerging-variant
sequences considered. For scenarios with co-circulating dominant
variants, the emerging variant is permitted to grow up to the size of
the largest dominant population. Beyond that point, the question changes
from detecting emergence to detecting dominance, which lies outside this
study’s scope.

| Scenario          | Ceiling                           |
|-------------------|-----------------------------------|
| Scenario 1        | $n_{V_1}$                         |
| Scenario 2        | $\max(n_{V_1}, n_{V_2})$          |
| Scenarios 3 and 4 | $\max(n_{V_1}, n_{V_2}, n_{V_3})$ |

## Grid construction

``` r
grid <- unique(round(exp(seq(log(2), log(ceiling), length.out = 50L))))
```

The same 50-point log spacing is applied in every band. Below
$n_{\text{emerge}} \approx 12$ the log-spaced points round to
consecutive integers and the grid is integer-precision; above that
threshold the relative gap between consecutive grid points grows
geometrically. Sweep evaluation halts at the first grid point that
detects; no refinement is performed.

### Representative grid sizes

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

50-point log-spaced grid: unique integer count per representative
ceiling.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Example
</th>

<th style="text-align:right;">

Ceiling
</th>

<th style="text-align:right;">

Unique points
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Sc1, N_ref=10, ratio_V1=1.25 (small band low)
</td>

<td style="text-align:right;">

12
</td>

<td style="text-align:right;">

11
</td>

</tr>

<tr>

<td style="text-align:left;">

All Sc, N_ref=100, ratio_V_i=5 (small band high)
</td>

<td style="text-align:right;">

500
</td>

<td style="text-align:right;">

44
</td>

</tr>

<tr>

<td style="text-align:left;">

All Sc, N_ref=1000, ratio_V_i=5 (medium band high)
</td>

<td style="text-align:right;">

5000
</td>

<td style="text-align:right;">

47
</td>

</tr>

<tr>

<td style="text-align:left;">

All Sc, N_ref=10000, ratio_V_i=5 (large band high)
</td>

<td style="text-align:right;">

50000
</td>

<td style="text-align:right;">

49
</td>

</tr>

</tbody>

</table>

The duplicate-collapse effect at very small ceilings (e.g., 12 unique
points reducing to 11 distinct integers after rounding) gives full
integer coverage of the small-end range.

### Explicit grid values for representative ceilings

The complete integer grids at the upper end of each band:

- **Ceiling 500** (small-band maximum, 44 unique points): 2, 3, 4, 5, 6,
  7, 8, 9, 10, 11, 12, 14, 15, 17, 19, 21, 24, 27, 30, 33, 37, 42, 47,
  53, 59, 66, 74, 82, 92, 103, 116, 129, 145, 162, 181, 203, 227, 254,
  285, 319, 357, 399, 447, 500.

- **Ceiling 5000** (medium-band maximum, 47 unique points): 2, 3, 4, 5,
  6, 7, 8, 10, 12, 14, 16, 19, 22, 26, 30, 35, 42, 49, 57, 67, 79, 92,
  108, 127, 149, 175, 205, 241, 282, 331, 389, 456, 535, 627, 736, 863,
  1013, 1188, 1394, 1635, 1918, 2250, 2640, 3097, 3633, 4262,

  5000. 

- **Ceiling 50000** (large-band maximum, 49 unique points): 2, 3, 4, 5,
  6, 7, 8, 10, 13, 16, 19, 24, 29, 36, 44, 55, 67, 83, 101, 125, 153,
  189, 232, 285, 351, 431, 530, 652, 801, 985, 1212, 1490, 1832, 2252,
  2770, 3405, 4187, 5148, 6330, 7784, 9571, 11768, 14469, 17791, 21875,
  26897, 33072, 40665, 50000.

## Grid-precision answers and within-cell aggregation

The 50-point log grid produces grid-precision threshold reports. The
worst-case relative grid gap per band, equal to $(r - 1)$ where
$r = (\text{ceiling}/2)^{1/49}$, is approximately 12% in the small band,
17% in the medium band, and 23% in the large band. Within-cell mean
estimates over 30 grid-aligned replicates are therefore biased upward by
at most half the local grid spacing, with relative bias bounded by 6%
(small), 8.5% (medium), and 11.5% (large) of the threshold value. The
within-cell variance estimate and the resulting confidence intervals are
unaffected by grid quantization, as the quantization is uniform across
replicates of one cell.

The Limitations section records this trade-off explicitly. Cell-mean
analysis remains valid; cell-median analysis (using the 15th or 16th
order statistic of 30 grid-aligned values) provides a quantization-
robust complement.

## Per-replicate GMM-fit budget

Maximum 50 GMM fits per replicate in the worst case (sweep traverses the
full grid without detection). In practice the sweep halts at the first
detecting grid point, so the per-replicate fit count typically ranges
from 1 to 30.

# Detection rule

A replicate is detected when all of the emerging variant’s mutation
sites are simultaneously in class 1 of the relabelled GMM, the
highest-entropy class after `relabel_entropy_classes` is applied. Two
degenerate paths return by a stricter equality test to prevent
false-positive detections:

- **All-identical entropy sentinel (class 999):** detection iff the
  surviving sites (those that pass the entropy filters `removez = TRUE`
  and `removesngl = TRUE`) are *exactly* the emerging-variant mutation
  positions.

- **Single-component fit ($G = 1$):** detection iff the surviving sites
  are *exactly* the emerging-variant mutation positions.

The strict-equality test prevents the degenerate paths from declaring
detection when established-variant sites also survive into the cluster.

The GMM parameters used by the detection rule are:

- `cluster_sites_by_entropy(entropies, nr = total_rows, G = 1:15)`, with
  `modelNames` omitted (the package default behaviour).
- Omitting `modelNames` invokes mclust’s default model search across
  both equal-variance (“E”) and variable-variance (“V”) univariate
  Gaussian mixtures, with BIC selecting the best (model, G) pair
  ([Scrucca et al. 2016](#ref-scrucca2016)). Forcing `modelNames = "V"`
  was tested but produced numerical degeneracies on inputs with exact
  entropy ties: the V model’s per-component variance estimate is
  undefined at zero-within-component variance, leaving $G = 1$ as the
  only well-defined BIC value. The simulator’s discrete population
  structure produces exact ties at all mutation sites of any single
  variant. The default search across E and V is therefore essential for
  correct cluster recovery in tied-data regimes; the E model is selected
  in those regimes and the V model is selected when the entropy
  distribution is sufficiently spread.
- $G = 1{:}15$ is a wider component-search range than the package’s
  pipeline default of $1{:}9$. The wider range is justified by this
  study’s much lighter per-fit cost (a single fit per detection test,
  with no temporal partitioning and no Hellinger matrix). BIC-based
  selection within the wider range provides finer entropy-class
  discrimination at negligible runtime cost.
- The chosen `(modelName, G)` is recorded per fit in each replicate’s
  RDS, supporting an audit of BIC behaviour in supplementary materials.
- Class 1 is the highest-entropy class after relabel, per the package
  convention. `relabel_entropy_classes` is the consumer’s
  responsibility; `cluster_sites_by_entropy` does not call it
  internally.

# Reproducibility

Three deterministic seeds suffice to reproduce the full study:

``` r
seed_cells          <- BASE_SEED
seed_replicate      <- BASE_SEED + scenario * 10^8 + cell_id * 10^3 + rep_id
seed_noise          <- seed_replicate + 1
```

The cell-table is precomputed once with `seed_cells` and serialised to
`outputs/cells_sc<N>.rds`. Each replicate’s variant-draw seed
(`seed_replicate`) is set inside the subprocess before mutation counts,
mutation positions, and substitute amino acids are drawn. The noise RNG
(`seed_noise`) is set immediately before the deleterious- noise
application; using a distinct offset means the same replicate at
different $p_{\text{del}}$ values gets the same variant placement but
independent noise patterns. The full study is reproducible from
`BASE_SEED` alone, modulo R version, package versions, and the FASTA
file checksum. All three quantities are recorded in the summary RDS via
`sessionInfo()`.

# Quality control and post-hoc reconstruction

A separate file, `reconstruct_replicate.R`, provides six post-hoc
reconstruction helpers used for quality control, failure-mode diagnosis,
and supplementary-figure generation after the production run completes.
They are not invoked by `simulation_study.R` itself; they exist to make
the on-disk results inspectable and the simulator deterministically
reproducible without storing the matrices themselves.

The motivation is storage cost. Saving every replicate’s
$(n_{\text{ref}} + \sum n_{V_i} + n_{\text{emerge}}) \times 1273$
integer matrix at production scale would consume approximately 20 TB of
raw disk (~2 TB compressed); see the Output structure section. We do not
save matrices. Every per-replicate RDS instead stores the seed, the
deleterious-noise rate, the variant placement parameters (via the cells
table), and the recorded sweep outcome. From these quantities, the full
simulator and detection pipeline are deterministically reproducible.

The six helpers are:

- `reconstruct_replicate_matrix(scenario, cell_id, rep_id, ...)`
  rebuilds the integer feature matrix at any chosen $n_{\text{emerge}}
  \in [0, \text{ceiling}]$. The output is bit-identical to the matrix
  the production run held in memory at the detecting sweep step. Useful
  for inspecting variant placement, the deleterious-noise pattern, or
  the row-by-row composition of any specific replicate.

- `reconstruct_replicate_pipeline(scenario, cell_id, rep_id, ...)`
  rebuilds the matrix and then replays the full detection pipeline at
  one sweep step: per-site entropy, `cluster_sites_by_entropy()` fit
  (including the `mclust::Mclust` object), `relabel_entropy_classes()`,
  and the per-emerging-position class assignments. Returns every
  intermediate object the production run discards.

- `list_na_replicates(config, scenarios, bands)` enumerates the
  replicates whose sweep returned NA. This is the natural starting point
  for any failure-mode diagnosis: the function returns a tidy data frame
  of (scenario, cell_id, rep_id) triples with the cell’s structural
  parameters attached, suitable for filtering and passing directly to
  the batch helpers.

- `reconstruct_na_replicates(na_table, what, workers)` is the batch
  driver. Given the output of `list_na_replicates()`, it loops the
  per-replicate reconstruction over every row, returning either the full
  pipeline objects (`what = "pipeline"`), the matrices only
  (`what = "matrix"`), or a tidy data frame of per-replicate diagnostics
  (`what = "summary"`) that records the detection-path label, the chosen
  $(\text{modelName}, G)$, the entropy values at the emerging positions,
  and the per-emerging-position class assignments.

- `reconstruct_noise_pattern(scenario, cell_id, rep_id, ...)` returns
  the exact set of $(\text{row}, \text{column})$ cells overwritten by
  the deleterious-noise step in one replicate. The output includes a
  boolean mask, a tidy data frame of noise events (with the pre-noise
  residue, the post-noise residue, the affected position, and the
  variant block each affected row originated from), per-row and
  per-position event counts, and the clean and noisy matrices for direct
  comparison. Computed by reconstructing the replicate twice — once
  without noise, once with — and differencing the two. Memory-heavy in
  the large band, where each call materialises two matrices of up to ~1
  GB each.

- `reconstruct_noise_patterns_batch(replicates_table, what, workers)` is
  the batch driver for the noise-pattern helper. It consumes the output
  of `list_na_replicates()` (or any equivalent subset table) and returns
  either a single tidy events data frame across all replicates
  (`what = "events"`), a per-replicate aggregate-counts data frame
  (`what = "summary"`), or the full noise-pattern objects keyed by
  replicate identifier (`what = "full"`). Reserved for surgical
  inspection rather than study-wide aggregation because of the
  per-replicate matrix cost; the typical use case is verifying noise
  uniformity across the variant blocks for a representative sample of
  cells, not enumerating every event in the production run.

Quality-control workflows the helpers support include:

1.  **Failure-mode diagnosis.**
    `reconstruct_na_replicates(..., what =    "summary")` returns the
    detection-path label (`multi_class`, `G1_collapse`, `sentinel_999`,
    `all_zero_entropy`) for every NA replicate, with the per-emerging-
    position class assignment counts. A
    `table(detection_path,    scenario)` of the result reveals whether
    NA replicates concentrate in a single failure mode (which would
    suggest a tightenable detection rule) or distribute across all four
    (which would indicate genuine geometric difficulty in the small-band
    hard corner of the design).

2.  **Per-site entropy auditing.** `reconstruct_replicate_pipeline()`
    exposes the 1,273-entry entropy vector at any sweep step. Plotting
    `sort(entropies)` and overlaying the emerging-variant positions
    shows directly where the variant’s sites sit in the entropy ranking,
    and at which $n_{\text{emerge}}$ they cross the threshold for class
    1 inclusion.

3.  **GMM fit auditing.** `reconstruct_replicate_pipeline()$cluster_fit`
    returns the full `mclust::Mclust` object, supporting inspection of
    per-component means, variances, mixing proportions, BIC values
    across the searched $(G, \text{modelName})$ pairs, and posterior
    probabilities. This is the basis for the supplementary materials
    that audit BIC behaviour under tied entropies (Section 9, Detection
    rule).

4.  **Noise-pattern auditing.** `reconstruct_noise_pattern()` confirms
    on a single replicate that the deleterious-noise mechanism is
    uniformly distributed across variant blocks and across the 1,273
    Spike positions, and that no spurious clustering or correlation has
    crept in. The `$by_position` and `$by_row_origin` summaries should
    match the binomial expectations
    $E[\text{events at position } s] = n_{\text{rows}} \cdot
    p_{\text{del}}$ and $E[\text{events in block } V_i] = n_{V_i}
    \cdot 1273 \cdot p_{\text{del}}$ up to Poisson noise. For
    replicate-level robustness checks across a representative subset of
    cells,
    `reconstruct_noise_patterns_batch(replicates_table,    what = "summary")`
    returns per-replicate event counts broken down by variant block;
    aggregating these counts (e.g. with
    `aggregate(n_events_total ~ band + scenario, ...)`) verifies that
    the empirical event rate matches the configured $p_{\text{del}}$
    uniformly across the design.

Because every reconstruction is deterministic from the recorded seed,
two researchers running the helpers on the same on-disk results will
obtain identical reconstructed matrices, entropies, and GMM fits. This
combined with the seed-recording convention in the per-replicate RDS
makes the simulation study fully reproducible without the ~20 TB matrix
archive.

# Output structure

    analysis/Sample_Size_Simulation_Study/
    |-- README.Rmd                  (this file)
    |-- README.html                 (knitted output)
    |-- references_sample_size.bib  (BibTeX source)
    |-- setup.R                     (libraries, config, helpers)
    |-- simulator.R                 (simulate_population_snapshot, apply_deleterious_noise)
    |-- detect_in_snapshot.R        (per-fit detection wrapper)
    |-- precompute_cells.R          (deterministic cell-table builder)
    |-- run_one_replicate.R         (per-replicate driver)
    |-- simulation_study.R          (orchestrator with end-of-run summary block)
    |-- empirical_pdel.R            (empirical p_del estimator on feature matrices)
    |-- reconstruct_replicate.R     (post-hoc QC: matrix + pipeline + noise reconstruction)
    |-- plot_results.R              (figures, results tables, Kaplan-Meier analysis)
    `-- outputs/
        |-- cells_sc<N>.rds                                 (4 files; deterministic cell tables)
        |-- replicates_sc<N>/sc<N>_cell<CCCC>_rep<RR>.rds   (per-replicate RDS, ~5 KB each)
        |-- summary_sc<N>.rds                               (tidy summary, one row per replicate)
        |-- error_log.txt                                   (orchestrator log; empty on a clean run)
        |-- empirical_pdel_cached.rds                       (GISAID VOC-site pdel results)
        |-- empirical_pdel_gisaid_meta.rds                  (GISAID dataset metadata: n_seqs, date range)
        |-- empirical_pdel_ncbi_cached.rds                  (NCBI VOC-site pdel results)
        |-- empirical_pdel_ncbi_meta.rds                    (NCBI dataset metadata: n_seqs, date range)
        |-- empirical_pdel_gisaid_overall_cached.rds        (GISAID genome-wide pdel results)
        |-- empirical_pdel_ncbi_overall_cached.rds          (NCBI genome-wide pdel results)
        |-- robustness_pdel/                                (optional; gated by RUN_ROBUSTNESS_SWEEP)
        |   `-- p_del_<value>/
        |       |-- replicates_sc<N>/sc<N>_cell<CCCC>_rep<RR>.rds
        |       `-- summary_sc<N>.rds
        |-- plots/                                          (built post hoc by generate_all_plots())
        |   |-- histograms.png                              (4 x 3 grid of per-cell distributions)
        |   |-- boxplots.png                                (notched boxplots per scenario, by band)
        |   |-- mean_ci.png                                 (mean +/- 95% BCa bootstrap CI)
        |   |-- detection_curves.png                        (empirical CDFs per (scenario, band))
        |   |-- km_curves.png                               (Kaplan-Meier curves with 95% CIs)
        |   `-- heatmap_sc<N>.png                           (per-scenario detection-rate heatmaps)
        `-- tables/                                         (built post hoc by generate_all_plots())
            |-- results_table.csv                           (per (scenario, band) summary)
            |-- results_table.rds                           (same, R-native)
            |-- km_summary.csv                              (Kaplan-Meier medians + 95% CIs)
            `-- km_summary.rds                              (same, R-native; logrank attribute attached)

Each per-replicate RDS (approximately 5 KB) holds: scenario, cell, and
replicate indices; structural parameters; the cell’s $n_{V_i}$ values;
the per-replicate deleterious-noise rate ($p_{\text{del}}$); the
emerging variant’s mutation count and positions; the executed sweep
grid; per-grid-point detection outcomes; per-grid-point chosen
`(modelName, G)` and class-1 site count; the detected flag and
$n_{\text{emerge}}^{\text{needed}}$; the seed and walltime; and R and
`ViralEntropR` version strings. The summary RDS written at the end of
each scenario is a tidy data frame with one row per replicate, ready for
`ggplot` and tidyverse aggregation.

# Memory budget

``` r
INT_BYTES <- 4L
COLS      <- 1273L

worst_per_band <- function(band) {
  n_ref_max <- max(n_ref_grid(band))
  n_v_max   <- round(n_ref_max * max(RATIOS))
  # Sc3/4 worst case: ref + 3 dominants + emerging at ceiling, all at max
  rows      <- n_ref_max + 3 * n_v_max + n_v_max
  matrix_gb <- rows * COLS * INT_BYTES / 1024^3
  list(n_ref     = n_ref_max,
       n_v_max   = n_v_max,
       ceiling   = n_v_max,
       rows      = rows,
       matrix_gb = matrix_gb)
}

mem_df <- do.call(rbind, lapply(names(BAND_RANGES), function(b) {
  w <- worst_per_band(b)
  data.frame(
    Band               = b,
    `N_ref max`        = format(w$n_ref,   big.mark = ","),
    `n_V_i max`        = format(w$n_v_max, big.mark = ","),
    `n_emerge ceiling` = format(w$ceiling, big.mark = ","),
    `Matrix rows`      = format(w$rows,    big.mark = ","),
    `Matrix size`      = sprintf("%.2f GB", w$matrix_gb),
    check.names        = FALSE
  )
}))

kbl(mem_df,
    caption = "Worst-case integer-matrix dimensions per band (scenarios 3 and 4, all ratios at max = 5).",
    align   = "lrrrrr") |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width        = FALSE)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Worst-case integer-matrix dimensions per band (scenarios 3 and 4, all
ratios at max = 5).
</caption>

<thead>

<tr>

<th style="text-align:left;">

Band
</th>

<th style="text-align:right;">

N_ref max
</th>

<th style="text-align:right;">

n_V_i max
</th>

<th style="text-align:right;">

n_emerge ceiling
</th>

<th style="text-align:right;">

Matrix rows
</th>

<th style="text-align:right;">

Matrix size
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

small
</td>

<td style="text-align:right;">

100
</td>

<td style="text-align:right;">

500
</td>

<td style="text-align:right;">

500
</td>

<td style="text-align:right;">

2,100
</td>

<td style="text-align:right;">

0.01 GB
</td>

</tr>

<tr>

<td style="text-align:left;">

medium
</td>

<td style="text-align:right;">

1,000
</td>

<td style="text-align:right;">

5,000
</td>

<td style="text-align:right;">

5,000
</td>

<td style="text-align:right;">

21,000
</td>

<td style="text-align:right;">

0.10 GB
</td>

</tr>

<tr>

<td style="text-align:left;">

large
</td>

<td style="text-align:right;">

10,000
</td>

<td style="text-align:right;">

50,000
</td>

<td style="text-align:right;">

50,000
</td>

<td style="text-align:right;">

210,000
</td>

<td style="text-align:right;">

1.00 GB
</td>

</tr>

</tbody>

</table>

Per-replicate peak memory in the large band is dominated by the sequence
matrix. The deleterious-noise step operates in-place on the matrix and
does not increase peak memory. Adding the per-site entropy vector, the
GMM internal state at `G = 1:15`, and the R session overhead, the
per-subprocess peak is approximately 1.5 GB in the large band, several
hundred MB in the medium band, and well under 1 GB in the small band.

Per-orchestrator peak memory with `N_WORKERS` concurrent subprocesses:

``` r
large_peak_gb <- worst_per_band("large")$matrix_gb + 0.5

hw_df <- data.frame(
  Hardware       = c("128 GB server, dedicated",
                     "16-core laptop, 32 GB",
                     "8-core laptop, 16 GB",
                     "4-core laptop, 8 GB",
                     "single-thread fallback"),
  `N_WORKERS`    = c(16L, 8L, 4L, 2L, 1L),
  `Total peak`   = sprintf("%.1f GB",
                           large_peak_gb * c(16, 8, 4, 2, 1)),
  check.names = FALSE
)

kbl(hw_df,
    caption = "Per-orchestrator peak memory by hardware tier (large-band worst case).",
    align   = "lrr") |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width        = FALSE)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Per-orchestrator peak memory by hardware tier (large-band worst case).
</caption>

<thead>

<tr>

<th style="text-align:left;">

Hardware
</th>

<th style="text-align:right;">

N_WORKERS
</th>

<th style="text-align:right;">

Total peak
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

128 GB server, dedicated
</td>

<td style="text-align:right;">

16
</td>

<td style="text-align:right;">

23.9 GB
</td>

</tr>

<tr>

<td style="text-align:left;">

16-core laptop, 32 GB
</td>

<td style="text-align:right;">

8
</td>

<td style="text-align:right;">

12.0 GB
</td>

</tr>

<tr>

<td style="text-align:left;">

8-core laptop, 16 GB
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

6.0 GB
</td>

</tr>

<tr>

<td style="text-align:left;">

4-core laptop, 8 GB
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

3.0 GB
</td>

</tr>

<tr>

<td style="text-align:left;">

single-thread fallback
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

1.5 GB
</td>

</tr>

</tbody>

</table>

These are peak values, encountered only at the upper-right corner of the
large-band grid (Sc3 or Sc4 with all ratios at the maximum). Most cells
consume far less memory.

# Runtime budget

``` r
PER_REP_TYP_S <- 2.5  # estimated typical per-replicate seconds (50-point grid)

hw_runtime_df <- data.frame(
  Hardware         = c("128 GB server, dedicated",
                       "16-core laptop, 32 GB",
                       "8-core laptop, 16 GB",
                       "4-core laptop, 8 GB",
                       "single-thread fallback"),
  `N_WORKERS`      = c(16L, 8L, 4L, 2L, 1L),
  `Wall-clock typ.`= sprintf("%.1f h",
                             total_sims * PER_REP_TYP_S /
                             c(16, 8, 4, 2, 1) / 3600),
  check.names = FALSE
)

kbl(hw_runtime_df,
    caption = sprintf(
      "Estimated wall-clock for %s total simulations, assuming approximately %.1f seconds typical per replicate. To be calibrated after a pilot run.",
      format(total_sims, big.mark = ","), PER_REP_TYP_S),
    align   = "lrr") |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width        = FALSE)
```

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Estimated wall-clock for 252,000 total simulations, assuming
approximately 2.5 seconds typical per replicate. To be calibrated after
a pilot run.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Hardware
</th>

<th style="text-align:right;">

N_WORKERS
</th>

<th style="text-align:right;">

Wall-clock typ.
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

128 GB server, dedicated
</td>

<td style="text-align:right;">

16
</td>

<td style="text-align:right;">

10.9 h
</td>

</tr>

<tr>

<td style="text-align:left;">

16-core laptop, 32 GB
</td>

<td style="text-align:right;">

8
</td>

<td style="text-align:right;">

21.9 h
</td>

</tr>

<tr>

<td style="text-align:left;">

8-core laptop, 16 GB
</td>

<td style="text-align:right;">

4
</td>

<td style="text-align:right;">

43.8 h
</td>

</tr>

<tr>

<td style="text-align:left;">

4-core laptop, 8 GB
</td>

<td style="text-align:right;">

2
</td>

<td style="text-align:right;">

87.5 h
</td>

</tr>

<tr>

<td style="text-align:left;">

single-thread fallback
</td>

<td style="text-align:right;">

1
</td>

<td style="text-align:right;">

175.0 h
</td>

</tr>

</tbody>

</table>

Per-replicate disk footprint is approximately 5 KB. Total study disk is
approximately 1.3 GB across the 252,000 per-replicate RDS files, plus
the four summary RDS files and the plots directory.

# Invocation

The orchestrator runs all four scenarios sequentially. Each scenario
writes to its own subdirectory; scenarios do not interfere. Resume
behaviour is implicit: at startup, the orchestrator skips any replicate
whose RDS already exists in the output directory.

``` bash
# Default settings (auto-detected N_WORKERS):
nohup Rscript analysis/Sample_Size_Simulation_Study/simulation_study.R \
    > analysis/Sample_Size_Simulation_Study/outputs/run_log.txt 2>&1 &

# Explicit parallel run with 16 workers:
N_WORKERS=16 nohup Rscript simulation_study.R > run.log 2>&1 &

# Sequential debug run (laptop):
N_WORKERS=1 Rscript simulation_study.R

# Sequential debug + robustness sweep:
N_WORKERS=1 RUN_ROBUSTNESS_SWEEP=1 Rscript simulation_study.R
```

Single-replicate debug:

``` r
source("analysis/Sample_Size_Simulation_Study/setup.R")
source("analysis/Sample_Size_Simulation_Study/simulator.R")
source("analysis/Sample_Size_Simulation_Study/detect_in_snapshot.R")
source("analysis/Sample_Size_Simulation_Study/run_one_replicate.R")

config <- build_config()
config$STUDY_DIR <- "analysis/Sample_Size_Simulation_Study"
ref_seq_int <- load_reference_sequence("data-raw/Ref_Seq.fasta")
empirical   <- build_empirical_constants()

result <- run_one_replicate(
  scenario    = 1L, cell_id = 1L, rep_id = 1L,
  cells       = readRDS("analysis/Sample_Size_Simulation_Study/outputs/cells_sc1.rds"),
  ref_seq_int = ref_seq_int, empirical = empirical,
  config      = config
)
```

# Methodology framing (ADEMP)

Following Morris, White, and Crowther ([2019](#ref-morris2019)), the
study is structured according to the ADEMP framework for simulation
studies in statistical methodology evaluation:

| Component | Specification |
|----|----|
| **A**ims | Quantify sample-size requirements for entropy-based GMM site selection under realistic background compositions, mutational loads, and surveillance-level deleterious-mutation noise. |
| **D**ata-generating mechanism | A snapshot simulator parameterised by `sarscov2_variants` (for variant structure) and by empirical GISAID per-position cumulative frequencies (for the deleterious-noise rate). |
| **E**stimands | Mean and median of $n_{\text{emerge}}^{\text{needed}}$ per cell, per band, per scenario, marginalised over band and mutation count; per-cell detection rate. |
| **M**ethods | `cluster_sites_by_entropy(entropies, nr, G = 1:15)` with mclust’s default E and V model search; `relabel_entropy_classes`; strict (all-sites-in-class-1) detection with strict-equality handlers for the $G = 1$ and 999-sentinel paths. |
| **P**erformance | Detection rate as a function of $n_{\text{emerge}}$; median $n_{\text{emerge}}^{\text{needed}}$ per cell; 30-replicate confidence intervals per cell; per-cell NA rate; robustness of inferred thresholds across $p_{\text{del}} \in \{0, 10^{-4}, 10^{-3}, 10^{-2}\}$ on a stratified-random subset of cells. |

The empirical grounding of every biological parameter
(`sarscov2_variants`-derived mutation counts, position pools, and the
D614G baseline; GISAID- and NCBI-derived noise rate) keeps the
simulation parameters in contact with surveillance data and avoids the
recurrent concern that simulation studies inflate effect sizes by using
parameters disconnected from real data ([Boulesteix, Lauer, and Eugster
2013](#ref-boulesteix2013)).

# Limitations

**(a) Snapshot model with no temporal dynamics.** The sample sizes
reported here apply to single-window detection. Temporal aggregation
effects (cumulative versus sliding versus disjoint partitioning) are
explored in `analysis/Entropy_Partitioning_Study/`,
`analysis/CP_Detection_Study/`and the `detecting_variants_simulation`
vignette.

**(b) Uniform substitute amino-acid draw.** Substitute residues are
drawn uniformly over the 19 non-reference standard amino acids. Real
mutations favour biochemically conservative changes and immune-escape
residues ([Harvey et al. 2021](#ref-harvey2021)); this bias is
unmodelled. The simplification is justified for the per-site entropy
signal, which depends on the presence of variation at a position rather
than the identity of the substitute, but is acknowledged as a limitation
for extensions to per-residue analyses.

**(c) Strict detection rule.** Partial-fraction matching, for example
detecting at least 75% of mutation sites in class 1, is not reported as
a primary output. Per-replicate metadata in the RDS files supports
post-hoc computation of fractional-match thresholds.

**(d) Fixed reference sequence.** All replicates share the canonical
pre-D614G 1,273-aa Spike [GenBank
YP_009724390](https://www.ncbi.nlm.nih.gov/protein/1796318598).
Sensitivity to alternative reference choices is not tested.

**(e) Position 614 trivially detected.** Every variant carries D614G, so
position 614 has very high cross-population entropy and is reliably in
class 1. The reported $n_{\text{emerge}}^{\text{needed}}$ therefore
reflects the harder challenge of detecting the variant’s
$n_{\text{muts}} - 1$ other mutations on top of the D614G baseline.

**(f) Independent random position draws across variants.** With three
established variants plus an emerging variant all drawing from
`POOL_11 \ {614}` (size 52), position collisions are expected. This is
biologically real (the next-most-shared empirical mutation appears in 7
variants) and is intentionally allowed by the simulator.

**(g) Grid-precision threshold reporting.** The 50-point log-spaced
sweep returns $n_{\text{emerge}}^{\text{needed}}$ at grid precision. The
worst-case relative grid gap is approximately 12% in the small band, 17%
in the medium band, and 23% in the large band, giving an upper bound of
approximately 6%, 8.5%, and 11.5% on the bias of the within-cell mean
estimate. Within-cell variance estimates and the resulting confidence
intervals are unaffected by the quantization. Cell-median analysis
(using the 15th or 16th order statistic of 30 grid-aligned values)
provides a quantization-robust complement.

**(h) Wild-type-only scenario omitted.** The pre-D614G scenario
(reference plus a single emerging variant on a wild-type background) is
not included; under that scenario the entropy-GMM step is trivially
successful at $n_{\text{emerge}} = 2$ because the emerging variant is
the sole source of variability in the population and every mutation site
has identical population structure.

# Session information

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26200)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8    LC_MONETARY=English_Canada.utf8
    ## [4] LC_NUMERIC=C                    LC_TIME=English_Canada.utf8    
    ## 
    ## time zone: America/Toronto
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] kableExtra_1.4.0   knitr_1.51         ViralEntropR_0.6.1 devtools_2.4.6     usethis_3.2.1     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1    viridisLite_0.4.3   dplyr_1.2.0         farver_2.1.2        Biostrings_2.78.0  
    ##  [6] S7_0.2.1            fastmap_1.2.0       xopen_1.0.1         digest_0.6.39       timechange_0.3.0   
    ## [11] lifecycle_1.0.5     ellipsis_0.3.2      processx_3.8.6      magrittr_2.0.4      compiler_4.5.2     
    ## [16] rlang_1.1.7         sass_0.4.10         tools_4.5.2         yaml_2.3.12         prettyunits_1.2.0  
    ## [21] labeling_0.4.3      pkgbuild_1.4.8      mclust_6.1.2        curl_7.0.0          xml2_1.5.1         
    ## [26] RColorBrewer_1.1-3  pkgload_1.5.0       withr_3.0.2         purrr_1.2.1         BiocGenerics_0.56.0
    ## [31] desc_1.4.3          grid_4.5.2          stats4_4.5.2        roxygen2_7.3.3      ggplot2_4.0.2      
    ## [36] scales_1.4.0        cli_3.6.5           rmarkdown_2.30      crayon_1.5.3        generics_0.1.4     
    ## [41] remotes_2.5.0       otel_0.2.0          rstudioapi_0.18.0   httr_1.4.8          commonmark_2.0.0   
    ## [46] sessioninfo_1.2.3   ecp_3.1.6           cachem_1.1.0        stringr_1.6.0       XVector_0.50.0     
    ## [51] vctrs_0.7.1         HDcpDetect_0.1.0    jsonlite_2.0.0      callr_3.7.6         IRanges_2.44.0     
    ## [56] rcmdcheck_1.4.0     S4Vectors_0.48.0    systemfonts_1.3.1   testthat_3.3.2      jquerylib_0.1.4    
    ## [61] glue_1.8.0          codetools_0.2-20    ps_1.9.1            lubridate_1.9.4     stringi_1.8.7      
    ## [66] gtable_0.3.6        tibble_3.3.1        pillar_1.11.1       htmltools_0.5.9     Seqinfo_1.0.0      
    ## [71] brio_1.1.5          R6_2.6.1            textshaping_1.0.4   rprojroot_2.1.1     evaluate_1.0.5     
    ## [76] lattice_0.22-9      memoise_2.0.1       bslib_0.10.0        Rcpp_1.1.1          svglite_2.2.2      
    ## [81] xfun_0.56           fs_1.6.6            zoo_1.8-15          pkgconfig_2.0.3

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Baker2017" class="csl-entry">

Baker, Timothy B., Linda M. Collins, Robin Mermelstein, et al. 2017.
“Implementing Clinical Research Using Factorial Designs.” *Nicotine &
Tobacco Research* 19 (2): 132–41. <https://doi.org/10.1093/ntr/ntw075>.

</div>

<div id="ref-bloomneher2023" class="csl-entry">

Bloom, J. D., and R. A. Neher. 2023. “Fitness Effects of Mutations to
SARS-CoV-2 Proteins.” *Virus Evolution* 9 (2): vead055.
<https://doi.org/10.1093/ve/vead055>.

</div>

<div id="ref-boulesteix2013" class="csl-entry">

Boulesteix, A.-L., S. Lauer, and M. J. A. Eugster. 2013. “A Plea for
Neutral Comparison Studies in Computational Sciences.” *PLoS ONE* 8 (4):
e61562. <https://doi.org/10.1371/journal.pone.0061562>.

</div>

<div id="ref-cochran1977" class="csl-entry">

Cochran, W. G. 1977. *Sampling Techniques*. 3rd ed. New York: Wiley.

</div>

<div id="ref-cohen1988" class="csl-entry">

Cohen, J. 1988. *Statistical Power Analysis for the Behavioral
Sciences*. 2nd ed. Hillsdale, NJ: Lawrence Erlbaum Associates.

</div>

<div id="ref-Desai2013" class="csl-entry">

Desai, Michael M., Aleksandra M. Walczak, and Daniel S. Fisher. 2013.
“Genetic Diversity and the Structure of Genealogies in Rapidly Adapting
Populations.” *Genetics* 193 (2): 565–85.
<https://doi.org/10.1534/genetics.112.147678>.

</div>

<div id="ref-harvey2021" class="csl-entry">

Harvey, W. T., A. M. Carabelli, B. Jackson, R. K. Gupta, E. C. Thomson,
E. M. Harrison, C. Ludden, et al. 2021. “SARS-CoV-2 Variants, Spike
Mutations and Immune Escape.” *Nature Reviews Microbiology* 19: 409–24.
<https://doi.org/10.1038/s41579-021-00573-0>.

</div>

<div id="ref-korber2020" class="csl-entry">

Korber, B., W. M. Fischer, S. Gnanakaran, H. Yoon, J. Theiler, W.
Abfalterer, N. Hengartner, et al. 2020. “Tracking Changes in SARS-CoV-2
Spike: Evidence That D614G Increases Infectivity of the COVID-19 Virus.”
*Cell* 182 (4): 812–827.e19.
<https://doi.org/10.1016/j.cell.2020.06.043>.

</div>

<div id="ref-Montgomery2017" class="csl-entry">

Montgomery, Douglas C. 2017. *Design and Analysis of Experiments*. 9th
ed. Hoboken, NJ: Wiley.

</div>

<div id="ref-morris2019" class="csl-entry">

Morris, T. P., I. R. White, and M. J. Crowther. 2019. “Using Simulation
Studies to Evaluate Statistical Methods.” *Statistics in Medicine* 38
(11): 2074–2102. <https://doi.org/10.1002/sim.8086>.

</div>

<div id="ref-neher2022" class="csl-entry">

Neher, R. A. 2022. “Contributions of Adaptation and Purifying Selection
to SARS-CoV-2 Evolution.” *Virus Evolution* 8 (2): veac113.
<https://doi.org/10.1093/ve/veac113>.

</div>

<div id="ref-PearsonLipman1988" class="csl-entry">

Pearson, William R., and David J. Lipman. 1988. “Improved Tools for
Biological Sequence Comparison.” *Proceedings of the National Academy of
Sciences* 85 (8): 2444–48. <https://doi.org/10.1073/pnas.85.8.2444>.

</div>

<div id="ref-plante2020" class="csl-entry">

Plante, J. A., Y. Liu, J. Liu, H. Xia, B. A. Johnson, K. G. Lokugamage,
X. Zhang, et al. 2021. “Spike Mutation D614G Alters SARS-CoV-2 Fitness.”
*Nature* 592: 116–21. <https://doi.org/10.1038/s41586-020-2895-3>.

</div>

<div id="ref-sayers2022" class="csl-entry">

Sayers, E. W., E. E. Bolton, J. R. Brister, K. Canese, J. Chan, D. C.
Comeau, R. Connor, et al. 2022. “Database Resources of the National
Center for Biotechnology Information.” *Nucleic Acids Research* 50 (D1):
D20–26. <https://doi.org/10.1093/nar/gkab1112>.

</div>

<div id="ref-scrucca2016" class="csl-entry">

Scrucca, L., M. Fop, T. B. Murphy, and A. E. Raftery. 2016. “Mclust 5:
Clustering, Classification and Density Estimation Using Gaussian Finite
Mixture Models.” *The R Journal* 8 (1): 289–317.
<https://doi.org/10.32614/RJ-2016-021>.

</div>

<div id="ref-Shannon1948" class="csl-entry">

Shannon, Claude E. 1948. “A Mathematical Theory of Communication.” *Bell
System Technical Journal* 27 (3): 379–423.
<https://doi.org/10.1002/j.1538-7305.1948.tb01338.x>.

</div>

<div id="ref-shu2017" class="csl-entry">

Shu, Y., and J. McCauley. 2017. “GISAID: Global Initiative on Sharing
All Influenza Data — from Vision to Reality.” *Eurosurveillance* 22
(13): 30494. <https://doi.org/10.2807/1560-7917.ES.2017.22.13.30494>.

</div>

<div id="ref-Symons2025" class="csl-entry">

Symons, J., C. Chung, B. M. Verheijen, S. J. Shemtov, D. de Jong, G.
Amatngalim, M. Nijhuis, M. Vermulst, and J. F. Gout. 2025. “The
Mutational Landscape of SARS-CoV-2 Provides New Insight into Viral
Evolution and Fitness.” *Nature Communications* 16 (1): 6425.
<https://doi.org/10.1038/s41467-025-61555-x>.

</div>

<div id="ref-Tao2021" class="csl-entry">

Tao, Kai, Phillip L. Tzou, Jiratchaya Nouhin, et al. 2021. “The
Biological and Clinical Significance of Emerging SARS-CoV-2 Variants.”
*Nature Reviews Genetics* 22 (12): 757–73.
<https://doi.org/10.1038/s41576-021-00408-x>.

</div>

<div id="ref-TonkinHill2021" class="csl-entry">

Tonkin-Hill, Gerry, Inigo Martincorena, Roberto Amato, Andrew R. J.
Lawson, Moritz Gerstung, Ian Johnston, David K. Jackson, et al. 2021.
“Patterns of Within-Host Genetic Diversity in SARS-CoV-2.” *eLife* 10:
e66857. <https://doi.org/10.7554/eLife.66857>.

</div>

<div id="ref-Tyuryaev2026" class="csl-entry">

Tyuryaev, Vadim, Jane Heffernan, and Hanna Jankowski. 2026.
*ViralEntropR: A Computational Pipeline for Entropy-Informed Detection
of Emerging Viral Variants*.
<https://github.com/vadimtyuryaev/ViralEntropR>.

</div>

<div id="ref-VanPoelvoorde2021" class="csl-entry">

Van Poelvoorde, Linde A. E., X Saelens, K Roose, et al. 2021. “Strategy
and Performance Evaluation of Low-Frequency Variant Calling for
SARS-CoV-2 Using Targeted Deep Illumina Sequencing.” *Frontiers in
Microbiology* 12: 747458. <https://doi.org/10.3389/fmicb.2021.747458>.

</div>

<div id="ref-Yang2022" class="csl-entry">

Yang, Wen Ting et al. 2022. “SARS-CoV-2 E484K Mutation Narrative Review:
Epidemiology, Immune Escape, Clinical Implications, and Future
Considerations.” *Infection and Drug Resistance* 15: 373–85.
<https://doi.org/10.2147/IDR.S344099>.

</div>

</div>
