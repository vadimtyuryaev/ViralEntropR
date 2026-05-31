ViralEntropR: Functional Clustering of Per-Site Shannon Entropy Curves
and Pairwise Hellinger Trajectories Around WHO Variant Emergence Dates
================
Vadim Tyuryaev, Jane Heffernan, Hanna Jankowski
May 28, 2026

- [Overview](#overview)
- [Terminology](#terminology)
- [Research questions](#research-questions)
- [Datasets](#datasets)
- [Ground truth](#ground-truth)
- [Methodology](#methodology)
  - [Per-cell window construction](#per-cell-window-construction)
  - [Three windowing strategies](#three-windowing-strategies)
  - [Per-window GMM site reduction](#per-window-gmm-site-reduction)
  - [Per-site entropy curves and Hellinger
    trajectories](#per-site-entropy-curves-and-hellinger-trajectories)
  - [Functional clustering via
    fdahclust](#functional-clustering-via-fdahclust)
  - [Per-frame animation pipeline](#per-frame-animation-pipeline)
    - [Three clustering-quality
      guards](#three-clustering-quality-guards)
    - [Oversized-cluster visualisation
      gate](#oversized-cluster-visualisation-gate)
    - [Dual-panel rendering](#dual-panel-rendering)
  - [Evaluation metrics](#evaluation-metrics)
- [Computational budget and runtime](#computational-budget-and-runtime)
- [Reproducibility](#reproducibility)
- [Implementation summary](#implementation-summary)
- [Per-cell output schema](#per-cell-output-schema)
- [Supplementary tables (data
  dictionary)](#supplementary-tables-data-dictionary)
  - [tab_01_internal_quality.csv](#tab_01_internal_qualitycsv)
  - [tab_02_cluster_membership.csv](#tab_02_cluster_membershipcsv)
  - [tab_03_stability.csv](#tab_03_stabilitycsv)
  - [tab_04_cross_strategy_ari.csv](#tab_04_cross_strategy_aricsv)
  - [tab_05_snp_enrichment.csv](#tab_05_snp_enrichmentcsv)
  - [tab_06_per_frame.csv](#tab_06_per_framecsv)
- [Limitations](#limitations)
- [ADEMP summary](#ademp-summary)
- [Invocation](#invocation)
- [Project structure](#project-structure)
- [Session information](#session-information)
- [References](#references)

# Overview

This study treats per-site Shannon entropy and pairwise Hellinger
distance, computed across temporal partitions of the SARS-CoV-2 Spike
protein, as **functional data** (Ramsay and Silverman 2005) and clusters
the resulting per-site curves and trajectories around the US first-
detection dates of the WHO-labelled SARS-CoV-2 variants. The clustering
is performed in a Functional Data Analysis (FDA) framework using
agglomerative hierarchical clustering (Sangalli et al. 2010) as
implemented in the `fdacluster` package, with the number of clusters
chosen by silhouette width (Rousseeuw 1987) over the grid
$k \in \{2, 3, 4, 5, 6\}$.

The design is factorial. **Twelve** WHO-labelled variants — five
Variants of Concern (Alpha, Beta, Gamma, Delta, Omicron) and seven
Variants of Interest (Epsilon, Zeta, Eta, Iota, Theta, Kappa, Lambda) —
are crossed with two independent surveillance datasets (NCBI US and
GISAID US) and three temporal-partitioning strategies (2-month sliding,
2-month disjoint, 1-month cumulative). A cell is one
`(variant, dataset, strategy)` combination. The cross is
$12 \times 2 \times 3 = 72$ cells in design.

Three terminal cell statuses occur in production:

- **`ok`** (52 cells) — full clustering pipeline ran end-to-end.
- **`underpowered_bin`** (17 cells) — at least one bin in the cell’s
  window contains fewer than `MIN_SEQUENCES_PER_BIN = 30L` sequences, so
  the per-site entropy and Hellinger estimates would be unreliable for
  that bin and the cell is excluded from clustering. Concentrated on
  low-prevalence VOIs (Epsilon, Eta, Iota, Zeta), where regionally
  sparse US deposition cannot support the bin-density threshold.
- **`out_of_coverage`** (3 cells) — the variant’s US detection date lies
  outside the dataset’s coverage. Affects only Omicron × NCBI US:
  Omicron was first detected in the US in December 2021, but the NCBI US
  feature matrix ends in September 2021.

Each `ok` cell yields **three** functional clustering solutions — one on
Shannon entropy curves, one on Hellinger trajectories anchored at the
window-start partition $T_1$, and one on Hellinger trajectories anchored
at the partition containing the *predecessor* variant’s US detection
date (denoted $T_{\text{pred}}$). The study therefore reports
$52 \times 3 = 156$ cell-level clustering solutions, together with the
corresponding static publication figures and animated visualisations.

The animation pipeline is **per-frame**: rather than fitting a single
cell-level clustering and rendering the same cluster assignments across
every bin, the pipeline fits a fresh GMM and fresh `fdahclust` on every
five-partition fitting sub-window (centred at each interior partition of
the cell). Each frame thus shows the clustering as it would have been
concluded *if a researcher had stopped sampling at that moment in time*.
The resulting GIFs reveal when the variant’s defining mutations converge
into a single cluster and how stable that convergence is across time.
See §[Methodology](#methodology) for details.

The headline question is descriptive rather than predictive: **do the
Spike positions retained by Gaussian Mixture Model (GMM) entropy site
selection organise into temporally coherent groups around each variant’s
emergence, and is that organisation reproducible across the three
windowing strategies and two independent surveillance datasets?** The
evaluation comprises four metric families covering internal clustering
quality (silhouette width), robustness to site resampling (bootstrap
Jaccard stability), methodological concordance (cross-strategy Adjusted
Rand Index), and biological plausibility (Fisher exact enrichment of
defining-SNP sites within clusters).

# Terminology

The following terms appear throughout the document. Definitions are
written for readers approaching the work from either the statistical,
virological, or functional-data-analysis side; cross-references to
standard sources are given where applicable.

| Term | Definition |
|----|----|
| Spike protein | The 1,273-amino-acid surface glycoprotein of SARS-CoV-2. The primary target of immune surveillance and of this study. Reference: GenBank accession YP_009724390. |
| Variant | A SARS-CoV-2 lineage carrying defining amino-acid substitutions, classified by the WHO. This study analyses all twelve WHO-labelled variants with non-missing `Date_First_Detected_US` in the `sarscov2_variants` catalogue: the five Variants of Concern (VOCs: Alpha, Beta, Gamma, Delta, Omicron) and seven Variants of Interest (VOIs: Epsilon, Zeta, Eta, Iota, Theta, Kappa, Lambda). |
| Surveillance dataset | A real-world feature matrix of dated Spike-protein sequences. This study uses NCBI US (Sayers et al. 2022) (137,132 raw sequences, 109,536 post-processed) and GISAID US (Shu and McCauley 2017) (129,371 post-processed). |
| Feature matrix | A `data.frame` with one row per dated sequence, one column for the collection date, and 1,273 integer-coded columns for the 1,273 Spike positions, encoded under a 25-symbol amino-acid alphabet (20 standard residues, the three IUPAC ambiguous codes B, Z, X, the stop codon `*`, and the alignment gap `-`). |
| Cell | One `(variant, dataset, strategy)` combination. The study has 72 cells in design; 69 are attempted (3 are `out_of_coverage`), of which 52 reach `ok` status and 17 are marked `underpowered_bin`. |
| Detection date $d_v$ | The first US detection date of variant $v$, taken from `sarscov2_variants$Date_First_Detected_US` and parsed to a first-of-month `Date` by `parse_month_year()`. |
| Per-cell window | The contiguous date range $[d_v - K \cdot b_{\max}, \; d_v + K \cdot b_{\max}]$ snapped to first-of-month boundaries, where $K = 5$ is the number of bins on each side and $b_{\max} = 2$ months is the largest bin width used. The window is clipped to the dataset’s coverage at runtime. When clipping at one edge reduces the realised $K$ on that side, the cell records both `K_pre_realised` and `K_post_realised`. |
| Partition / bin | A contiguous date sub-range within the per-cell window over which residue frequencies and the per-site entropy are computed. The three strategies differ in how the window is partitioned into bins (see below). Bins are denoted $T_1, T_2, \ldots, T_n$ in temporal order. |
| Strategy 1: 2-month sliding | A 2-month bin advances one month at a time. A 20-month window yields $n = 19$ overlapping bins. Adjacent bins share 1 month of data and differ by 1 month. `partition_time_windows(window_type = 2L, window_length = 2L)`. |
| Strategy 2: 2-month disjoint | Non-overlapping 2-month bins. A 20-month window yields $n = 10$ bins, each containing 2 months of distinct data. `partition_time_windows(window_type = 3L, window_length = 2L)`. |
| Strategy 3: 1-month cumulative | An expanding bin that begins at the window’s start and grows by one month at each step. Bin $T_k$ contains all sequences from window-start through month $k$. A 20-month window yields $n = 20$ nested bins. `partition_time_windows(window_type = 1L, window_length = 1L)`. |
| Reference bin $T_1$ | The first bin of the window. For Hellinger trajectories its residue distribution serves as one of two anchor distributions (see *Dual anchors* below). |
| Predecessor anchor $T_{\text{pred}}$ | The partition containing the US detection date of the focal variant’s *predecessor* in the canonical Spike-evolution lineage (Alpha → D614G, Beta/Gamma/Theta/Kappa/Lambda → Alpha, Delta → Beta, Omicron → Delta, Epsilon/Zeta/Eta/Iota → D614G). When the predecessor’s detection date falls outside the cell window, the anchor is snapped to the nearest in-window partition and the cell records `anchor_Tpred$snapped = TRUE`. |
| Shannon entropy $H(s, T_k)$ | The Shannon entropy of the residue distribution at Spike position $s$ in bin $T_k$, in bits: $H(s, T_k) = -\sum_{a} \hat p_a \log_2 \hat p_a$, where $\hat p_a$ is the empirical proportion of residue $a$ at position $s$ among the bin’s sequences. Computed by `ViralEntropR::calculate_entropy()`. |
| Per-site entropy curve | The vector $\mathbf{e}_s = (H(s, T_1), H(s, T_2), \ldots, H(s, T_n))$. Functional clustering treats $\mathbf{e}_s$ as the discrete sampling of a continuous function of time and clusters sites with similar functions. |
| Hellinger distance $D_H(P, Q)$ | A bounded metric on probability distributions: $D_H(P, Q) = \frac{1}{\sqrt 2}\sqrt{\sum_a (\sqrt{p_a} - \sqrt{q_a})^2}$. Bounded in $[0, 1]$ when normalised; this study uses the unnormalised form, bounded in $[0, \sqrt 2]$. Reference: (Vaart 1998, ch. 14). |
| Per-site Hellinger trajectory | The vector $\mathbf{h}_s = (D_H(P_{s, T_2}, P_{s, T_a}), \ldots, D_H(P_{s, T_n}, P_{s, T_a}))$, of length $n - 1$, where $T_a \in \{T_1, T_{\text{pred}}\}$ is the chosen anchor. Computed by `ViralEntropR::calculate_hellinger_matrix()` with the anchor passed as `reference_idx`. |
| GMM (Gaussian Mixture Model) | A probabilistic clustering model that fits a weighted sum of $G$ Gaussian components to a univariate sample. This study uses `mclust::Mclust` (Scrucca et al. 2016) via `ViralEntropR::cluster_sites_by_entropy()` over $G \in \{1, \ldots, 15\}$ components and selects the best fit by Bayesian Information Criterion. |
| Class 1 (highest-entropy sites) | After `relabel_entropy_classes()` standardises labels so that class 1 = highest mean entropy, this is the set of Spike positions in the highest-entropy GMM component for the cell’s window. Functional clustering is applied only to these class-1 sites; sites in classes 2 to $G$ are not analysed. |
| FDA (Functional Data Analysis) | A statistical framework that treats observations indexed by an ordered variable (here, time) as samples of an underlying smooth function. References: (Ramsay and Silverman 2005; Kokoszka and Reimherr 2017). |
| `fdahclust` | The hierarchical agglomerative clustering routine in the `fdacluster` package (Sangalli et al. 2010), which builds a dendrogram of functional curves under a user-specified linkage and metric, with optional curve alignment. This study uses `warping_class = "none"` (no curve alignment) and the default $L^2$ metric. |
| log1p transformation | The transformation $x \mapsto \log(1 + x)$, applied to both entropy curves and Hellinger trajectories before functional clustering. Reduces the leverage of high-value tails while preserving order, maps zero to zero, and avoids $\log 0$ divergence. Reference: (Box and Cox 1964). |
| Silhouette width $s_i$ | A standard internal clustering quality measure. For each clustered curve $i$, $s_i = (b_i - a_i) / \max(a_i, b_i)$ where $a_i$ is its mean within-cluster distance and $b_i$ is the smallest mean between-cluster distance. Average $\bar s$ ranges in $[-1, 1]$; $> 0.5$ decisive, $0.25$–$0.5$ reasonable, $< 0.25$ weak. Reference: (Rousseeuw 1987). |
| Bootstrap Jaccard stability | The stability of a clustering under site resampling (Hennig 2007). Sites are resampled with replacement $B = 100$ times, the clustering is recomputed on each resample, and for each original cluster the maximum Jaccard similarity against any resample cluster is recorded. Per-cluster mean Jaccard $\geq 0.75$ indicates a stable cluster; $0.60$–$0.75$ a patterned cluster; $< 0.60$ instability. |
| ARI (Adjusted Rand Index) | A measure of agreement between two clusterings of the same items, corrected for chance. ARI $= 1$ indicates perfect agreement; $0$ indicates agreement at the chance level; $< 0$ indicates worse-than-chance agreement (Hubert and Arabie 1985). This study reports cross-strategy ARI to quantify whether the choice of partitioning strategy materially changes the clustering. |
| SNP enrichment (Fisher) | A one-tailed Fisher exact test of independence between cluster assignment and defining-SNP-site membership, computed on the cell’s class-1 sites only (Agresti 1992). Tests whether the variant’s defining mutations concentrate disproportionately in one cluster. Multiple-comparison correction uses Bonferroni across `N_TOTAL_ENRICHMENT_TESTS = 220` (12 variants × 2 datasets × 3 strategies × 3 metrics minus the 9 metric-level slots in the 3 `out_of_coverage` cells, plus a small safety margin). |
| `underpowered_bin` | Cell status assigned when at least one bin in the per-cell window contains fewer than `MIN_SEQUENCES_PER_BIN = 30L` sequences. The cell records the worst offending bin in `status_reason` and is excluded from cell-level FDA, all downstream metric tables, and animations. |
| `out_of_coverage` | Cell status assigned when the variant’s US detection date $d_v$ lies outside the dataset’s coverage range. Affects only Omicron × NCBI US. |

# Research questions

- **Q1.** Around each variant’s US first-detection date, do the GMM-
  selected class-1 Spike positions partition into temporally coherent
  groups of entropy curves and Hellinger trajectories under
  agglomerative functional clustering?
- **Q2.** How sensitive is the resulting clustering to the choice of
  temporal-partitioning strategy (2-month sliding vs 2-month disjoint vs
  1-month cumulative)? The cross-strategy ARI quantifies the answer per
  cell.
- **Q3.** Is the clustering reproducible across two independent
  surveillance datasets (NCBI US and GISAID US) for the variants whose
  detection dates fall within both datasets’ coverage?
- **Q4.** Does the clustering organise the variant’s defining-SNP sites
  into a single cluster, providing biological plausibility for the
  cluster structure?
- **Q5.** How robust are the clusters to perturbation of the input
  curves by site resampling? Bootstrap Jaccard stability addresses the
  question per cluster within each cell.
- **Q6.** *When* do the defining sites converge into a single cluster?
  The per-frame design records, for each `(cell, metric)`, the first
  frame at which the SNP / Mutation / class-1 convergence predicates
  become TRUE; reported in
  `tab_01_internal_quality.csv$convergence_first_*`.

# Datasets

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Surveillance datasets analysed.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Dataset
</th>

<th style="text-align:left;">

Date range
</th>

<th style="text-align:right;">

Months
</th>

<th style="text-align:left;">

Sequences (post-processed)
</th>

<th style="text-align:left;">

Source
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

NCBI US
</td>

<td style="text-align:left;">

January 2020 – September 2021
</td>

<td style="text-align:right;">

21
</td>

<td style="text-align:left;">

109,536
</td>

<td style="text-align:left;">

Zenodo DOI 10.5281/zenodo.19040165
</td>

</tr>

<tr>

<td style="text-align:left;">

GISAID US
</td>

<td style="text-align:left;">

January 2020 – March 2024
</td>

<td style="text-align:right;">

51
</td>

<td style="text-align:left;">

129,371
</td>

<td style="text-align:left;">

GISAID EpiCov, processed via `analysis/GISAID_data_preprocessing/`
</td>

</tr>

</tbody>

</table>

Both feature matrices are integer-encoded under the 25-symbol alphabet,
sorted ascending by collection date, and date-clipped to the first of
each collection month at preprocessing time (Sayers et al. 2022; Shu and
McCauley 2017). Site columns occupy positions 1 through 1,273; the
`Date` column follows. The clipping is automatic via
`AL_df$Date <- as.Date(format(corrected_dates, "%Y-%m-01"))` in both
preprocessing pipelines, which means every sequence’s `Date` is
unambiguously assignable to a bin in any window whose start and end are
themselves month boundaries.

# Ground truth

The ground-truth catalogue is the `sarscov2_variants` dataset bundled
with `ViralEntropR`, restricted to the twelve WHO-labelled variants with
a non-missing `Date_First_Detected_US`. For each variant, three pieces
of catalogue metadata enter the analysis:

- **Detection date $d_v$** = `Date_First_Detected_US` after
  `parse_month_year()` conversion to a first-of-month `Date`. Used as
  the centre of the per-cell window.
- **Defining-SNP sites** = `Defining_SNP_Sites[[v]]`, the integer vector
  of Spike positions whose amino-acid substitutions define the variant
  in standard nomenclature. Used as the biological reference set in the
  SNP enrichment test (Q4).
- **Mutation sites** = `Mutation_Sites[[v]]`, the broader set of
  positions reported as mutated in the variant’s defining lineage. Used
  as a secondary convergence target in the per-frame analysis.

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Variant detection dates, predecessor map, and dataset-coverage status.
Omicron’s December 2021 first-US-detection date lies outside NCBI’s
September 2021 coverage, dropping three cells (3 strategies × Omicron ×
NCBI) from the 72-cell design.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Variant
</th>

<th style="text-align:left;">

Class
</th>

<th style="text-align:left;">

US first detection
</th>

<th style="text-align:left;">

Predecessor
</th>

<th style="text-align:left;">

In NCBI US window
</th>

<th style="text-align:left;">

In GISAID US window
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Alpha
</td>

<td style="text-align:left;">

VOC
</td>

<td style="text-align:left;">

Dec 2020
</td>

<td style="text-align:left;">

D614G
</td>

<td style="text-align:left;">

yes
</td>

<td style="text-align:left;">

yes
</td>

</tr>

<tr>

<td style="text-align:left;">

Beta
</td>

<td style="text-align:left;">

VOC
</td>

<td style="text-align:left;">

Jan 2021
</td>

<td style="text-align:left;">

Alpha
</td>

<td style="text-align:left;">

yes
</td>

<td style="text-align:left;">

yes
</td>

</tr>

<tr>

<td style="text-align:left;">

Gamma
</td>

<td style="text-align:left;">

VOC
</td>

<td style="text-align:left;">

Jan 2021
</td>

<td style="text-align:left;">

Alpha
</td>

<td style="text-align:left;">

yes
</td>

<td style="text-align:left;">

yes
</td>

</tr>

<tr>

<td style="text-align:left;">

Delta
</td>

<td style="text-align:left;">

VOC
</td>

<td style="text-align:left;">

Mar 2021
</td>

<td style="text-align:left;">

Beta
</td>

<td style="text-align:left;">

yes
</td>

<td style="text-align:left;">

yes
</td>

</tr>

<tr>

<td style="text-align:left;">

Omicron
</td>

<td style="text-align:left;">

VOC
</td>

<td style="text-align:left;">

Dec 2021
</td>

<td style="text-align:left;">

Delta
</td>

<td style="text-align:left;">

no (out-of-coverage)
</td>

<td style="text-align:left;">

yes
</td>

</tr>

<tr>

<td style="text-align:left;">

Epsilon
</td>

<td style="text-align:left;">

VOI
</td>

<td style="text-align:left;">

Mar 2020
</td>

<td style="text-align:left;">

D614G
</td>

<td style="text-align:left;">

yes
</td>

<td style="text-align:left;">

yes
</td>

</tr>

<tr>

<td style="text-align:left;">

Zeta
</td>

<td style="text-align:left;">

VOI
</td>

<td style="text-align:left;">

Mar 2020
</td>

<td style="text-align:left;">

D614G
</td>

<td style="text-align:left;">

yes
</td>

<td style="text-align:left;">

yes
</td>

</tr>

<tr>

<td style="text-align:left;">

Eta
</td>

<td style="text-align:left;">

VOI
</td>

<td style="text-align:left;">

Nov 2020
</td>

<td style="text-align:left;">

D614G
</td>

<td style="text-align:left;">

yes
</td>

<td style="text-align:left;">

yes
</td>

</tr>

<tr>

<td style="text-align:left;">

Iota
</td>

<td style="text-align:left;">

VOI
</td>

<td style="text-align:left;">

Nov 2020
</td>

<td style="text-align:left;">

D614G
</td>

<td style="text-align:left;">

yes
</td>

<td style="text-align:left;">

yes
</td>

</tr>

<tr>

<td style="text-align:left;">

Theta
</td>

<td style="text-align:left;">

VOI
</td>

<td style="text-align:left;">

Mar 2021
</td>

<td style="text-align:left;">

Alpha
</td>

<td style="text-align:left;">

yes
</td>

<td style="text-align:left;">

yes
</td>

</tr>

<tr>

<td style="text-align:left;">

Kappa
</td>

<td style="text-align:left;">

VOI
</td>

<td style="text-align:left;">

Apr 2021
</td>

<td style="text-align:left;">

Alpha
</td>

<td style="text-align:left;">

yes
</td>

<td style="text-align:left;">

yes
</td>

</tr>

<tr>

<td style="text-align:left;">

Lambda
</td>

<td style="text-align:left;">

VOI
</td>

<td style="text-align:left;">

Aug 2021
</td>

<td style="text-align:left;">

Alpha
</td>

<td style="text-align:left;">

yes
</td>

<td style="text-align:left;">

yes
</td>

</tr>

</tbody>

</table>

# Methodology

## Per-cell window construction

For each `(variant, dataset)` pair, a candidate window is constructed as
$[d_v - K \cdot b_{\max}, \; d_v + K \cdot b_{\max}]$, where $K = 5$ and
$b_{\max} = 2$ months, giving a 20-month nominal window. The candidate
window is then snapped to first-of-month boundaries and clipped to the
dataset’s available date range. A cell is marked `out_of_coverage` and
excluded from production when the detection date itself lies outside the
dataset’s coverage; this affects only Omicron × NCBI US.

A further bin-density check is applied at the strategy level. Once the
strategy has partitioned the (clipped) window into bins, any cell with
at least one bin containing fewer than `MIN_SEQUENCES_PER_BIN = 30L`
sequences is marked `underpowered_bin` and excluded from cell-level FDA
and all downstream tables. The `status_reason` field records the worst
offending bin (e.g., *“2 / 16 bins have \< 30 seqs (worst = 2)”*), which
lets downstream readers distinguish a near-miss from a structural
failure.

## Three windowing strategies

Within the per-cell window, three partitioning strategies are applied in
parallel. Each strategy generates its own bin grid, its own per-site
entropy curves, and its own Hellinger trajectories. The intent is to
make the methodological knob “how do we discretise time?” an *explicit
experimental factor* rather than a hidden assumption.

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Three windowing strategies applied to each (variant, dataset) cell.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Strategy
</th>

<th style="text-align:left;">

partition_time_windows args
</th>

<th style="text-align:left;">

Bins per 20-mo window
</th>

<th style="text-align:left;">

Properties
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

2-month sliding
</td>

<td style="text-align:left;">

window_length = 2L, window_type = 2L
</td>

<td style="text-align:left;">

19
</td>

<td style="text-align:left;">

Overlapping; adjacent bins share 1 month; high temporal resolution at
the cost of dependence between adjacent bins.
</td>

</tr>

<tr>

<td style="text-align:left;">

2-month disjoint
</td>

<td style="text-align:left;">

window_length = 2L, window_type = 3L
</td>

<td style="text-align:left;">

10
</td>

<td style="text-align:left;">

Independent; non-overlapping 2-month chunks; the cleanest statistical
interpretation.
</td>

</tr>

<tr>

<td style="text-align:left;">

1-month cumulative
</td>

<td style="text-align:left;">

window_length = 1L, window_type = 1L
</td>

<td style="text-align:left;">

20
</td>

<td style="text-align:left;">

Nested; each bin contains all previous data; gives a smoothed,
monotone-coverage trajectory.
</td>

</tr>

</tbody>

</table>

## Per-window GMM site reduction

For each cell, the GMM site-reduction step of the ViralEntropR pipeline
is applied *within the cell’s window only* — no oracle leakage from the
rest of the dataset is permitted. The procedure mirrors
`analysis/Changepoint_Detection_Study/helpers_hellinger.R::fit_window_gmm()`:

1.  Per-site Shannon entropy is computed across all sequences in the
    window (one scalar per Spike position).
2.  Invariant sites (entropy $\approx 0$) and singleton sites (entropy
    $= H(\frac{1}{n_w}, \frac{n_w - 1}{n_w})$, where $n_w$ is the number
    of sequences in the window) are filtered.
3.  `mclust::Mclust` is fitted over $G \in \{1, \ldots, 15\}$ with model
    selection by BIC.
4.  `relabel_entropy_classes()` standardises labels so that class 1 =
    the GMM component with the highest mean entropy.
5.  The class-1 site set is extracted. The cell is marked
    `reduced_skipped` (and excluded from the metric tables) when the GMM
    returns the all-identical 999 sentinel, when class 1 is empty, or
    when class 1 contains fewer than `MIN_CLASS1_SITES_FOR_FDA = 4`
    sites. The minimum-site threshold is set so that silhouette over
    $k \in \{2, \ldots, 6\}$ is well- defined at the upper end of the
    grid.

The GMM is fit *once* per cell (not once per strategy), because the GMM
operates on across-window entropy and is therefore strategy- agnostic.
The class-1 site set is then propagated to all three strategies of that
cell.

## Per-site entropy curves and Hellinger trajectories

For each strategy of each cell, with the class-1 site set fixed:

- **Per-site entropy curve** $\mathbf{e}_s$ is the vector of per-bin
  Shannon entropies at position $s$. Computed by `apply()`-ing
  `ViralEntropR::calculate_entropy()` over the bin’s sequence matrix.
  Shape: $(n_{\text{class-1 sites}}, n_{\text{bins}})$ matrix of entropy
  values.
- **Per-site Hellinger trajectory** $\mathbf{h}_s$ is the vector of
  pairwise Hellinger distances from a designated *anchor bin* to every
  other bin in temporal order, at position $s$. Computed by
  `ViralEntropR::calculate_hellinger_matrix()` with the anchor bin
  passed as `reference_idx` and the class-1 sites passed as the `sites`
  argument. Shape: $(n_{\text{class-1 sites}}, n_{\text{bins}} - 1)$.

**Dual anchors.** Each cell produces *two* Hellinger trajectories rather
than one:

1.  **$T_1$ anchor.** The anchor is the window-start partition. This is
    the conventional choice and answers “how has the residue
    distribution at site $s$ evolved relative to the start of the
    observation window?”.
2.  **$T_{\text{pred}}$ anchor.** The anchor is the partition containing
    the predecessor variant’s US detection date. The predecessor map
    encodes the canonical Spike-evolution lineage: Alpha and the early
    VOIs (Epsilon, Zeta, Eta, Iota) anchor on D614G (April 2020, (Korber
    et al. 2020)); Beta, Gamma, Theta, Kappa, and Lambda anchor on
    Alpha; Delta anchors on Beta; Omicron anchors on Delta. When the
    predecessor’s detection date lies outside the cell’s window, the
    anchor is snapped to the nearest in-window partition and the cell
    records `anchor_Tpred$snapped = TRUE`.

Both matrices are then log1p-transformed elementwise prior to
clustering. The transformation $x \mapsto \log(1 + x)$ leaves zeros
fixed, is order-preserving, and substantially reduces the leverage of
high-entropy / high-Hellinger sites on the $L^2$ distance between
curves.

## Functional clustering via fdahclust

Each $(n_{\text{class-1 sites}} \times n_{\text{time points}})$ matrix
is then passed to `fdacluster::fdahclust()` with the following
parameters:

- `warping_class = "none"`: no curve alignment is performed. The time
  axis is bin index, not biological time within a bin, and we do not
  want clusters to be defined by within-bin temporal warping.
- Default $L^2$ metric between curves.
- `centroid_type = "mean"`, `linkage = "complete"`.
- Number of clusters $k$ chosen by average silhouette width over
  $k \in \{2, 3, 4, 5, 6\}$.

`fdahclust` returns a hierarchical agglomerative clustering object with
cluster assignments per site. The cell’s clustering result is serialised
together with the entropy and Hellinger matrices, the chosen $k$, the
silhouette value, and the GMM and window metadata.

## Per-frame animation pipeline

In addition to the cell-level (full-window) clustering described above,
the pipeline computes a *per-frame* clustering used by the animation
system. A **frame** corresponds to a five-partition fitting sub-window
centred at an interior partition of the cell’s window. With $h =$
`FITTING_WINDOW_HALFWIDTH = 2`, the set of valid frame centres is
$\{h+1, \ldots, n_{\text{bins}} - h\}$, giving 15 frames for sliding-2m
cells (19 bins), 6 frames for disjoint-2m cells (10 bins), and 16 frames
for cumulative-1m cells (20 bins).

For each frame, three independent re-fits run for the three metrics
(entropy, Hellinger $T_1$, Hellinger $T_{\text{pred}}$):

1.  **Per-frame GMM.** All sequences in the five partitions centred at
    the frame’s central bin are pooled, per-site Shannon entropy is
    recomputed, and a fresh `Mclust` fit selects the per-frame class-1
    site set. The cell-level class-1 set is *not* propagated.
2.  **Per-frame curve matrix.** Per-site entropy curves or Hellinger
    trajectories are computed for the frame’s class-1 sites *over the
    full cell window* (not just the five partitions). The fitting window
    selects which sites enter clustering; the rendered curves span the
    entire window so the reader can see what happens outside the fitting
    slice.
3.  **Per-frame `fdahclust` + silhouette.** $k$ is selected by mean
    silhouette over the grid $\{2, \ldots, 6\}$, exactly as in the
    cell-level fit. Per-site silhouette widths are recorded so that
    marginal sites can be rendered separately in the visualisation
    (state-1 sites with $s_i < 0$).
4.  **Convergence record.** Three boolean predicates are evaluated:
    `converged_snp` (all variant-defining SNP sites that are present in
    this frame’s class-1 share a single cluster), `converged_mut` (same
    for all variant-defining mutation sites), and `converged_class1`
    (all well-clustered class-1 sites of this frame share a single
    cluster). The “first frame at which the predicate becomes TRUE” is
    the headline number reported per cell.

### Three clustering-quality guards

A frame’s clustering is treated as informative only when all three of
the following guards pass:

- **Guard A.** The per-frame GMM must produce $G > 1$. When $G = 1$ the
  frame is *null*: there is no class-1 site set and the bottom panel
  renders all variant-defining sites red.
- **Guard B.** The frame’s mean silhouette width must be $\geq 0.25$,
  the Kaufman–Rousseeuw “reasonable structure” threshold (Kaufman and
  Rousseeuw 1990). Frames with $\bar s < 0.25$ are marked
  `structurally_null` and convergence predicates return FALSE regardless
  of cluster assignment patterns.
- **Guard C.** Per-site silhouette $s_i \geq 0$. Class-1 sites with
  $s_i < 0$ are marginal — they straddle clusters and would be drawn to
  a different cluster centroid on a small input perturbation. Such sites
  are rendered red and excluded from the convergence predicates.

### Oversized-cluster visualisation gate

A separate visualisation-only gate handles the case where the clustering
produces a single very-large cluster that would otherwise dominate the
plot. A cluster is **oversized** when its size exceeds $M =$
`MAX_CLUSTER_SIZE = 33` (Omicron’s defining mutation count, the largest
among the 12 WHO variants and therefore a biologically motivated upper
bound on what a real variant signature could look like). Sites assigned
to an oversized cluster are rendered in light grey in both panels. The
convergence predicates also require all qualifying sites to share a
cluster that is *not* oversized — i.e., an “all sites in one giant
cluster” event does not count as convergence.

The gate is purely visual. The underlying RDS payload and the
`tab_02_cluster_membership.csv` table record the unfiltered cluster
assignments. Readers and downstream analyses see the full picture; the
plot is just rendered without the visually overwhelming cluster.

### Dual-panel rendering

Each frame is rendered as a two-panel composite. The **top panel** shows
the entropy or Hellinger curves over the full cell window: all 1,273
Spike positions appear as faint light-grey background lines, the frame’s
class-1 sites are drawn over them in cluster colours (red if marginal,
light grey if oversized), a vertical dashed line marks the variant’s
detection bin, and two open triangles ($\triangleright$ /
$\triangleleft$) mark the fitting sub-window boundaries. The **bottom
panel** is a two-row track of the variant’s defining sites along
positions 1–1,273: row “SNP” is `Defining_SNP_Sites`, row “Mutation” is
`Mutation_Sites`, and each site is coloured by its frame state (cluster
colour if well-clustered class-1, light grey if oversized cluster, red
otherwise). The frames are composed into an animated GIF per
`(cell, metric)` at `ANIMATION_FPS = 1L`, and a static “summary” PNG
using the full- window cell-level fit is also written per
`(cell, metric)`.

## Evaluation metrics

Four metric families are computed per cell, per metric type (entropy,
Hellinger $T_1$, Hellinger $T_{\text{pred}}$).

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Four metric families and their cross-cell granularity.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Family
</th>

<th style="text-align:left;">

Metric
</th>

<th style="text-align:left;">

Reference
</th>

<th style="text-align:left;">

Granularity
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Internal quality
</td>

<td style="text-align:left;">

Silhouette width $\bar s$
</td>

<td style="text-align:left;">

(Rousseeuw 1987)
</td>

<td style="text-align:left;">

One value per `(cell, metric)`. Reported in
`tab_01_internal_quality.csv`.
</td>

</tr>

<tr>

<td style="text-align:left;">

Stability under resampling
</td>

<td style="text-align:left;">

Bootstrap Jaccard stability
</td>

<td style="text-align:left;">

(Hennig 2007)
</td>

<td style="text-align:left;">

One value per `(cell, metric, cluster)`. Reported in
`tab_03_stability.csv`.
</td>

</tr>

<tr>

<td style="text-align:left;">

Cross-strategy concordance
</td>

<td style="text-align:left;">

Adjusted Rand Index (ARI)
</td>

<td style="text-align:left;">

(Hubert and Arabie 1985)
</td>

<td style="text-align:left;">

One value per `(variant, dataset, metric, strategy_pair)`. Reported in
`tab_04_cross_strategy_ari.csv`.
</td>

</tr>

<tr>

<td style="text-align:left;">

Biological plausibility
</td>

<td style="text-align:left;">

Fisher exact SNP enrichment
</td>

<td style="text-align:left;">

(Agresti 1992)
</td>

<td style="text-align:left;">

One Fisher $p$-value and Bonferroni-adjusted $p$-value per
`(cell, metric)`. Reported in `tab_05_snp_enrichment.csv`.
</td>

</tr>

</tbody>

</table>

A fifth table, `tab_02_cluster_membership.csv`, is descriptive (not a
metric) and lists each cluster’s site composition and SNP overlap; it is
the bridge between the clustering and its biological interpretation. A
sixth table, `tab_06_per_frame.csv`, surfaces the per-frame diagnostics
used by the animation pipeline. All six CSVs are written by
`aggregate_tables.R` from the per-cell RDS files in `outputs/cells/`.
Their column-level data dictionaries are given in §[Supplementary tables
(data dictionary)](#supplementary-tables-data-dictionary).

The choice of metrics targets four independent failure modes. A
clustering can be:

- *Decisively separated but unstable* (high silhouette, low Jaccard) —
  the cluster boundaries are well-resolved on the observed data but
  would re-organise under modest perturbation.
- *Stable but uninformative* (high Jaccard, low silhouette) — the
  algorithm consistently produces the same partition, but the partition
  is poorly motivated.
- *Internally good but strategy-dependent* (high silhouette and Jaccard,
  low cross-strategy ARI) — the chosen partitioning strategy materially
  determines the conclusion.
- *Internally good but biologically incongruent* (high silhouette and
  Jaccard, non-significant SNP enrichment) — the clusters are real
  statistically but do not align with the variant’s known mutation
  profile.

Reporting all four allows reviewers to triangulate cluster meaning
rather than relying on any single quality index.

# Computational budget and runtime

The numbers below are **measured** from the production run on the shared
server `tcell` (125 GiB RAM, 60-core class; Ubuntu Linux), not pre-run
estimates.

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Measured wall-clock budget per analysis phase. Phase 1 totals 17 h 48
min, dominated by Omicron × `cumulative_1m` and Delta × `cumulative_1m`
on GISAID US.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Phase
</th>

<th style="text-align:left;">

Wall-clock time (measured)
</th>

<th style="text-align:left;">

N workers
</th>

<th style="text-align:left;">

Outcome
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Phase 1a: cell production (first pass)
</td>

<td style="text-align:left;">

11 h 49 min
</td>

<td style="text-align:left;">

12 (callr::r_bg)
</td>

<td style="text-align:left;">

67 / 72 cells written; 5 SIGKILL (OOM) on Omicron / Alpha / Beta / Kappa
/ Lambda × heavy strategies
</td>

</tr>

<tr>

<td style="text-align:left;">

Phase 1b: cell production (OOM-recovery pass)
</td>

<td style="text-align:left;">

5 h 59 min
</td>

<td style="text-align:left;">

5 (callr::r_bg)
</td>

<td style="text-align:left;">

5 / 5 OOM-affected cells re-run cleanly
</td>

</tr>

<tr>

<td style="text-align:left;">

Phase 2: aggregate_tables.R
</td>

<td style="text-align:left;">

6 s
</td>

<td style="text-align:left;">

1
</td>

<td style="text-align:left;">

6 supplementary CSVs written (216 / 5058 / 469 / 141 / 216 / 1764 rows)
</td>

</tr>

<tr>

<td style="text-align:left;">

Phase 3: plot_results.R
</td>

<td style="text-align:left;">

≈ 40–50 min
</td>

<td style="text-align:left;">

1
</td>

<td style="text-align:left;">

1,920 PNGs written (52 ok cells × 3 metrics × ≈ 12 frames + summary
panels)
</td>

</tr>

<tr>

<td style="text-align:left;">

Phase 4: render_existing.R (post-hoc cluster_membership CSVs)
</td>

<td style="text-align:left;">

40 min 32 s
</td>

<td style="text-align:left;">

8 (xargs -P 8)
</td>

<td style="text-align:left;">

156 `cluster_membership_&lt;metric&gt;.csv` files written (52 ok cells ×
3 metrics)
</td>

</tr>

<tr>

<td style="text-align:left;">

Phase 5: animate_results.R
</td>

<td style="text-align:left;">

single-pass; runtime recorded in the run log
</td>

<td style="text-align:left;">

1
</td>

<td style="text-align:left;">

Up to 156 `anim_&lt;metric&gt;.gif` files written
</td>

</tr>

</tbody>

</table>

**Memory footprint per worker.** The feature matrix is loaded once at
subprocess startup (≈ 140–165 MB), and the per-cell working set is
bounded by the window sub-matrix (≈ 10–40 MB depending on coverage) plus
the GMM fit, the three `fdacluster` cell-level objects, the per-frame
`fdacluster` / `Mclust` intermediates, and the full 1273-site entropy
and Hellinger matrices (≈ 200 KB each). Sustained worker RSS during the
first-pass run ranged 3–6.5 GB, with transient spikes substantially
above that on Omicron × cumulative_1m × GISAID US and a handful of other
heavy cells — five of these spikes exceeded the 10 GB/worker headroom
implied by `N_WORKERS = 12` on the 125 GB server and were terminated by
the kernel OOM-killer. These five cells were re-run in the Phase 1b
recovery pass at `N_WORKERS = 5` (≈ 25 GB/worker headroom), which
succeeded for all five. **Recommended production setting on hosts with ≥
64 GB RAM: `N_WORKERS = 6` or `N_WORKERS = 8`**, which both gives every
cell enough headroom for the worst case and avoids the over-subscription
regime that caused the OOM kills.

**Storage footprint per cell (RDS):** ≈ 600 KB for the three full-
1273-site matrices, ≈ 150 KB for the three class-1 matrices, ≈ 100 KB
for the per-frame records (~15 frames × 3 metrics × ~25 fields each), ≈
50 KB for the cell-level fits / stability / enrichment. Total ≈ 1 MB per
cell × 72 cells ≈ 70 MB of cell RDS files. The OOM-recovered cells
observed in the production run ranged 270 KB (NCBI × Alpha ×
cumulative_1m, terminal status `underpowered_bin` after the recovery) to
1.06 MB (GISAID × Lambda × sliding_2m).

**Figure / animation storage:** ≈ 1,920 frame PNGs × 500 KB + 156
summary PNGs × 500 KB + up to 156 GIFs × 1 MB ≈ 1.2 GB total in
`outputs/plots/` and `outputs/gifs/`.

# Reproducibility

- **Deterministic seeding.** Each cell receives a deterministic seed:
  `BASE_SEED + dataset_id*1e8 + variant_id*1e6`. With
  `BASE_SEED = 2025L`, `dataset_id ∈ {1, 2}`, and
  `variant_id ∈ {1, …, 12}`, the maximum seed is well within
  `.Machine$integer.max` and is used to seed both the cell-level GMM EM
  initialisation, the per-frame GMM fits (each with an additional
  `+100 + centre_bin` offset for the frame), and the bootstrap Jaccard
  resampling.
- **Atomic writes.** All RDS and CSV files are written via
  `save_rds_atomic()` / `write_csv_atomic()`, which rename the file from
  a per-PID temporary, ensuring no partially written file can be loaded
  by a downstream step.
- **Schema-aware resume.** On restart, `fda_analysis.R` scans
  `outputs/cells/cell_*.rds`, loads each, and skips any cell whose
  `schema_version` matches the current run. Older v1 RDS files are
  flagged stale and re-run. Downstream scripts (`plot_results.R`,
  `animate_results.R`, `aggregate_tables.R`, `render_existing.R`)
  regenerate their outputs from the cell RDS files and are therefore
  cheap to re-run.
- **No oracle leakage.** Per-cell GMM is fit on the cell’s window only,
  not on the full dataset. Per-frame GMM is fit on the frame’s five
  partitions only. Class-1 site selection is therefore locally
  determined at every level; differences across cells reflect biology,
  not training protocol.
- **Version pinning.** R and `ViralEntropR` version strings are recorded
  in every per-cell RDS and reproduced in the session- information block
  at the end of this README.

# Implementation summary

The implementation follows the same module conventions established in
`Changepoint_Detection_Study/`:

- **Concurrency.** `N_WORKERS = 1L` runs sequentially in-process with
  `tryCatch` per cell; `N_WORKERS > 1L` dispatches via `callr::r_bg`
  with a worker pool; `N_WORKERS = NULL` auto-detects with a
  RAM-budgeted cap of 64.
- **I/O.** Atomic RDS and CSV writes; the per-cell schema is documented
  under §[Per-cell output schema](#per-cell-output-schema) below.
- **Resume.** On restart, the orchestrator scans `outputs/cells/` for
  existing files and re-runs only the missing ones.
- **Logging.** Structured log via `log_msg()` to
  `outputs/logs/orchestrator_<timestamp>.log`; per-cell errors caught
  and recorded in `outputs/logs/error_log.txt` without interrupting the
  orchestrator. All `stop()` and `warning()` calls use `call. = FALSE`
  per package convention.

Environment overrides at orchestrator startup:

``` bash
N_WORKERS=8 VIRAL_FDA_STUDY_DIR=/path/to/FDA_Analysis Rscript fda_analysis.R
```

The production runs used `N_WORKERS = 12` for the first pass and
`N_WORKERS = 5` for the OOM-recovery pass. On hosts with ≥ 64 GB RAM,
`N_WORKERS = 6`–`8` is recommended as the steady-state setting (see
§[Computational budget and runtime](#computational-budget-and-runtime)
for the empirical justification).

# Per-cell output schema

Each `outputs/cells/cell_<dataset>__<variant>__<strategy>.rds` contains
a named list with the following 61 fields:

``` r
list(
  schema_version,                       # 2L

  # Cell identifiers
  dataset, variant, strategy,
  dataset_id, variant_id, strategy_id,

  # Window construction
  detection_date,                       # parse_month_year(...)
  window_start, window_end,             # Date, snapped to month boundary
  K_pre_realised, K_post_realised,      # may differ from K = 5 at edges
  n_seqs_window,                        # total sequences in window

  # Predecessor metadata
  predecessor_name,                     # "D614G", another variant, or NA
  predecessor_us_date,                  # Date or NA
  predecessor_status,                   # "ok" / "no_predecessor" / "snapped" / ...

  # Bin grid (strategy-dependent)
  bin_starts, bin_ends, bin_labels, bin_counts,
  n_bins, detection_bin,

  # Hellinger anchor indices
  anchor_T1    = list(bin_idx, snapped, distance_months),
  anchor_Tpred = list(bin_idx, snapped, distance_months),

  # GMM site reduction (full-window, shared across metrics)
  gmm_meta = list(G, modelName, n_class1_sites,
                  class1_sites, reduced_skipped, reason),
  snp_sites_truth,                      # integer vector
  mutation_sites_truth,                 # integer vector

  # FULL 1273-site matrices for background panels
  entropy_full_matrix,                  # 1273 x n_bins
  hellinger_T1_full_matrix,             # 1273 x (n_bins - 1)
  hellinger_Tpred_full_matrix,          # 1273 x (n_bins - 1) or NULL

  # Class-1 matrices, raw and log1p-transformed (three metrics)
  entropy_class1_matrix,
  hellinger_T1_class1_matrix,
  hellinger_Tpred_class1_matrix,
  entropy_class1_matrix_transformed,
  hellinger_T1_class1_matrix_transformed,
  hellinger_Tpred_class1_matrix_transformed,
  transformation,                       # "log1p" or "identity"

  # Cell-level FDA fits (three metrics)
  fda_entropy             = list(fit, n_clusters, silhouette_mean,
                                  silhouette_per_k, cluster_assignments,
                                  transformation, status, status_reason),
  fda_hellinger_T1        = list(...),
  fda_hellinger_Tpred     = list(...),

  # Bootstrap Jaccard stability (three metrics)
  stability_entropy        = list(n_bootstrap, per_cluster_jaccard,
                                   n_stable_clusters, n_failed_bootstrap),
  stability_hellinger_T1   = list(...),
  stability_hellinger_Tpred = list(...),

  # SNP Fisher enrichment (three metrics)
  snp_enrichment_entropy        = list(fisher_p, fisher_p_bonferroni,
                                        n_snp_in_class1, n_class1_total,
                                        snp_overlap_per_cluster,
                                        max_overlap_cluster, max_overlap_OR,
                                        contingency_table, status),
  snp_enrichment_hellinger_T1   = list(...),
  snp_enrichment_hellinger_Tpred = list(...),

  # Per-frame records
  frame_centres,                        # integer vector of centre bin indices
  frames_entropy          = list(...),  # list of per-frame records
  frames_hellinger_T1     = list(...),
  frames_hellinger_Tpred  = list(...),

  # Convergence summary: first centre_bin at which each predicate holds
  convergence_entropy            = list(first_snp, first_mut, first_class1),
  convergence_hellinger_T1       = list(first_snp, first_mut, first_class1),
  convergence_hellinger_Tpred    = list(first_snp, first_mut, first_class1),

  # Status and run metadata
  status,                               # "ok", "out_of_coverage",
                                        # "underpowered_bin", "reduced_skipped",
                                        # "too_few_sites", "error"
  status_reason,                        # character or NA
  walltime_s, seed_cell, seed_strategy,
  r_version, package_version, run_timestamp
)
```

Each element of `frames_<metric>` is itself a named list with shape
defined in `helpers_frames.R::run_one_frame()` and includes the
per-frame class-1 set, curve matrix (raw and transformed), cluster
assignments, per-site silhouette widths, oversized-cluster flags,
convergence flags, and the GMM status of the frame.

The cross-strategy ARI is computed by `aggregate_tables.R` from the
three same-`(dataset, variant, metric)` cell files; it is not stored per
cell.

# Supplementary tables (data dictionary)

`aggregate_tables.R` rolls the 72 per-cell RDS files into six
supplementary CSVs under `outputs/tables/`. The six tables differ in
granularity: some have one row per cell (`tab_01`, `tab_05`), some one
row per cell × cluster (`tab_03`), one per cell × site (`tab_02`), one
per `(dataset, variant, metric, strategy_pair)` (`tab_04`), and one per
cell × frame (`tab_06`). Row counts reported below are from the
production run on `tcell`.

The per-cell `cluster_membership_<metric>.csv` files written by
`render_existing.R` under
`outputs/plots/<dataset>/<variant>/<strategy>/` are a different artefact
— they expose the per-site cluster assignment as it is *rendered* in the
GIF, including the per-frame attribution. The aggregate
`tab_02_cluster_membership.csv` is the cell-level analogue.

## tab_01_internal_quality.csv

**Granularity:** one row per `(cell, metric)`. **Row count:** 216 = 72
cells × 3 metrics. **Purpose:** the headline cell-level quality table —
silhouette width, GMM model summary, chosen $k$, per-metric status, and
the three convergence-first-frame indices. This is the table to consult
first when interpreting a cell.

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Column dictionary for `tab_01_internal_quality.csv`.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Column
</th>

<th style="text-align:left;">

Type
</th>

<th style="text-align:left;">

Description
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

dataset
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Dataset name: `NCBI_US` or `GISAID_US`.
</td>

</tr>

<tr>

<td style="text-align:left;">

variant
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

WHO variant label, one of the twelve variants in `sarscov2_variants`.
</td>

</tr>

<tr>

<td style="text-align:left;">

strategy
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Windowing strategy: `sliding_2m` (2-month bins, 1-month step),
`disjoint_2m` (non-overlapping 2-month bins), or `cumulative_1m`
(expanding bins, 1-month step).
</td>

</tr>

<tr>

<td style="text-align:left;">

metric
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Curve type clustered. `entropy` = per-site Shannon-entropy curve (one
entropy value per partition, in bits). `hellinger_T1` = per-site
Hellinger trajectory anchored at the window-start partition `T1`
(distance from each subsequent partition’s residue distribution to
`T1`’s). `hellinger_Tpred` = per-site Hellinger trajectory anchored at
`T_pred`, the partition containing the predecessor variant’s US
detection date.
</td>

</tr>

<tr>

<td style="text-align:left;">

cell_status
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Cell-level terminal status. `ok` = clustering ran end-to-end.
`underpowered_bin` = at least one bin had fewer than
`MIN_SEQUENCES_PER_BIN = 30L` sequences. `out_of_coverage` = the
variant’s US detection date is outside the dataset’s coverage.
`reduced_skipped` = the GMM fit was degenerate. `too_few_sites` =
class-1 had fewer than `MIN_CLASS1_SITES_FOR_FDA = 4L` sites. `error` =
an uncaught exception.
</td>

</tr>

<tr>

<td style="text-align:left;">

cell_status_reason
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Free-text reason matching the status (e.g.,
`"2 / 16 bins have &lt; 30 seqs (worst = 2)."` for `underpowered_bin`);
`NA` for `ok` cells.
</td>

</tr>

<tr>

<td style="text-align:left;">

predecessor_name
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Variant whose US detection date anchors `hellinger_Tpred`. `D614G` is
the canonical April-2020 reference used for Alpha and the early VOIs
(Epsilon, Zeta, Eta, Iota); other variants chain through the lineage
(Beta/Gamma/Theta/Kappa/Lambda → Alpha, Delta → Beta, Omicron → Delta).
</td>

</tr>

<tr>

<td style="text-align:left;">

n_class1_sites
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of class-1 (highest-entropy) Spike positions selected by the
cell’s window-level GMM; `NA` when the cell was rejected before the GMM
step (`underpowered_bin`, `out_of_coverage`).
</td>

</tr>

<tr>

<td style="text-align:left;">

gmm_G
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of GMM components selected by BIC over `G` in 1, …, 15; `NA` when
the GMM did not run.
</td>

</tr>

<tr>

<td style="text-align:left;">

gmm_modelName
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

`mclust` covariance model name selected by BIC (e.g., `V` =
variable-variance univariate, `E` = equal-variance, `X` =
single-component degenerate).
</td>

</tr>

<tr>

<td style="text-align:left;">

n_bins
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of partitions the strategy produced for this cell (after window
clipping). For an unclipped 20-month window: 19 for `sliding_2m`, 10 for
`disjoint_2m`, 20 for `cumulative_1m`.
</td>

</tr>

<tr>

<td style="text-align:left;">

n_frames
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of fitting frames computed for this metric. Each frame is a
five-partition fitting sub-window centred on an interior partition.
</td>

</tr>

<tr>

<td style="text-align:left;">

k_selected
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of clusters selected for this metric by mean silhouette over `k`
in 2, …, 6; `NA` when the metric was skipped (cell not `ok`, or class-1
too small).
</td>

</tr>

<tr>

<td style="text-align:left;">

silhouette_mean
</td>

<td style="text-align:left;">

num
</td>

<td style="text-align:left;">

Mean silhouette width across all clustered sites at the chosen `k`.
Ranges in -1 to 1; the Kaufman & Rousseeuw thresholds are 0.25
(reasonable), 0.50 (decisive), 0.70 (strong).
</td>

</tr>

<tr>

<td style="text-align:left;">

metric_status
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Per-metric terminal status (`ok`, `skipped_too_few_sites`, `error`, …).
Distinct from `cell_status`: a cell can be `ok` while a particular
metric (e.g., `hellinger_Tpred`) is skipped because its anchor falls
outside the window.
</td>

</tr>

<tr>

<td style="text-align:left;">

convergence_first_snp
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Index of the first frame centre at which all defining-SNP sites present
in the frame’s class-1 share a single, non-oversized cluster (Guard
predicate `converged_snp`). `NA` if convergence never occurred within
this metric’s frames.
</td>

</tr>

<tr>

<td style="text-align:left;">

convergence_first_mut
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Same as `convergence_first_snp` but for the variant’s broader
defining-mutation site set (`converged_mut`).
</td>

</tr>

<tr>

<td style="text-align:left;">

convergence_first_class1
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Same as `convergence_first_snp` but for all well-clustered class-1 sites
of the frame (`converged_class1`).
</td>

</tr>

<tr>

<td style="text-align:left;">

walltime_s
</td>

<td style="text-align:left;">

num
</td>

<td style="text-align:left;">

Per-cell wall-clock seconds (the three metric rows for one cell repeat
this value).
</td>

</tr>

</tbody>

</table>

## tab_02_cluster_membership.csv

**Granularity:** one row per `(cell, metric, site)`, restricted to cells
with `status == "ok"` and metrics whose `fda_*$status == "ok"`. **Row
count:** 5,058. **Purpose:** the bridge between the clustering and the
biology — for each clustered site, which cluster it lives in and whether
it is a defining-SNP or broader-mutation site for the variant.

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Column dictionary for `tab_02_cluster_membership.csv`.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Column
</th>

<th style="text-align:left;">

Type
</th>

<th style="text-align:left;">

Description
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

dataset
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Dataset name: `NCBI_US` or `GISAID_US`.
</td>

</tr>

<tr>

<td style="text-align:left;">

variant
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

WHO variant label.
</td>

</tr>

<tr>

<td style="text-align:left;">

strategy
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Windowing strategy (`sliding_2m`, `disjoint_2m`, or `cumulative_1m`).
</td>

</tr>

<tr>

<td style="text-align:left;">

metric
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Curve type clustered. `entropy` = per-site Shannon-entropy curve over
time (bits). `hellinger_T1` = per-site Hellinger trajectory from each
subsequent partition’s residue distribution to the window-start
partition `T1`’s. `hellinger_Tpred` = same, but anchored at `T_pred`
(the partition containing the predecessor variant’s US detection date).
</td>

</tr>

<tr>

<td style="text-align:left;">

cluster_id
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Cluster ID assigned by `fdacluster::fdahclust` (1-based). The chosen `k`
(number of clusters) is in `tab_01_internal_quality.csv$k_selected`;
valid IDs run from 1 through that `k`.
</td>

</tr>

<tr>

<td style="text-align:left;">

site
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Spike position (1 to 1,273) belonging to the cluster.
</td>

</tr>

<tr>

<td style="text-align:left;">

is_snp
</td>

<td style="text-align:left;">

lgl
</td>

<td style="text-align:left;">

`TRUE` if `site` is in the variant’s defining-SNP set
(`sarscov2_variants$Defining_SNP_Sites[[variant]]`).
</td>

</tr>

<tr>

<td style="text-align:left;">

is_mutation_site
</td>

<td style="text-align:left;">

lgl
</td>

<td style="text-align:left;">

`TRUE` if `site` is in the variant’s broader defining-mutation set
(`sarscov2_variants$Mutation_Sites[[variant]]`).
</td>

</tr>

</tbody>

</table>

## tab_03_stability.csv

**Granularity:** one row per `(cell, metric, cluster_id)`, restricted to
cells with `status == "ok"` and metrics whose bootstrap procedure ran.
**Row count:** 469. **Purpose:** per-cluster Hennig (2007) bootstrap
Jaccard stability. The `is_stable` flag uses
`STABLE_JACCARD_THRESHOLD = 0.75`; rerunning aggregation with a
different threshold (e.g., 0.60) would re-classify clusters.

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Column dictionary for `tab_03_stability.csv`.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Column
</th>

<th style="text-align:left;">

Type
</th>

<th style="text-align:left;">

Description
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

dataset
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Dataset name: `NCBI_US` or `GISAID_US`.
</td>

</tr>

<tr>

<td style="text-align:left;">

variant
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

WHO variant label.
</td>

</tr>

<tr>

<td style="text-align:left;">

strategy
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Windowing strategy (`sliding_2m`, `disjoint_2m`, or `cumulative_1m`).
</td>

</tr>

<tr>

<td style="text-align:left;">

metric
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Curve type clustered (see `tab_01$metric` for definitions).
</td>

</tr>

<tr>

<td style="text-align:left;">

cluster_id
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Cluster ID (1-based, identifies the same cluster as in
`tab_02_cluster_membership.csv`).
</td>

</tr>

<tr>

<td style="text-align:left;">

mean_jaccard
</td>

<td style="text-align:left;">

num
</td>

<td style="text-align:left;">

Hennig (2007) bootstrap Jaccard stability: the mean over `B = 100`
site-resampling bootstrap rounds of the maximum Jaccard similarity
between this original cluster’s site set and any cluster of the
resampled fit. Ranges in 0 to 1; closer to 1 = more stable.
</td>

</tr>

<tr>

<td style="text-align:left;">

is_stable
</td>

<td style="text-align:left;">

lgl
</td>

<td style="text-align:left;">

`TRUE` iff `mean_jaccard &gt;= STABLE_JACCARD_THRESHOLD = 0.75`
(Hennig’s “stable” tier). The 0.60-0.75 range is “patterned”; below 0.60
is unstable.
</td>

</tr>

<tr>

<td style="text-align:left;">

n_bootstrap
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of bootstrap rounds attempted (`N_BOOTSTRAP = 100L`).
</td>

</tr>

<tr>

<td style="text-align:left;">

n_failed_bootstrap
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of bootstrap rounds in which the recompute failed (e.g., the
resampled site set was too small for `fdahclust` or yielded a degenerate
fit).
</td>

</tr>

</tbody>

</table>

## tab_04_cross_strategy_ari.csv

**Granularity:** one row per
`(dataset, variant, metric, strategy_pair)`, restricted to
`(dataset, variant)` groups where both strategies of the pair returned
`status == "ok"` and the common class-1 site set has at least 3
elements. **Row count:** 141. **Purpose:** Adjusted Rand Index between
the two cluster assignments restricted to their common class-1 sites,
quantifying methodological concordance (Q2). The three pairs are
`(sliding_2m, disjoint_2m)`, `(sliding_2m, cumulative_1m)`, and
`(disjoint_2m, cumulative_1m)`.

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Column dictionary for `tab_04_cross_strategy_ari.csv`.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Column
</th>

<th style="text-align:left;">

Type
</th>

<th style="text-align:left;">

Description
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

dataset
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Dataset name: `NCBI_US` or `GISAID_US`.
</td>

</tr>

<tr>

<td style="text-align:left;">

variant
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

WHO variant label.
</td>

</tr>

<tr>

<td style="text-align:left;">

metric
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Curve type clustered (see `tab_01$metric` for definitions).
</td>

</tr>

<tr>

<td style="text-align:left;">

strategy_a
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Name of the first strategy in the pair (`sliding_2m`, `disjoint_2m`, or
`cumulative_1m`).
</td>

</tr>

<tr>

<td style="text-align:left;">

strategy_b
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Name of the second strategy in the pair.
</td>

</tr>

<tr>

<td style="text-align:left;">

n_sites_common
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of class-1 sites present in both strategies (size of the
intersection of the two class-1 site sets). Rows with fewer than 3
common sites are dropped because ARI on so few items is uninformative.
</td>

</tr>

<tr>

<td style="text-align:left;">

k_a
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of clusters selected for `strategy_a` (matches
`tab_01_internal_quality.csv$k_selected` for that
`(dataset, variant, strategy_a, metric)` row).
</td>

</tr>

<tr>

<td style="text-align:left;">

k_b
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of clusters selected for `strategy_b`.
</td>

</tr>

<tr>

<td style="text-align:left;">

ARI
</td>

<td style="text-align:left;">

num
</td>

<td style="text-align:left;">

Adjusted Rand Index \[`mclust::adjustedRandIndex`\] between the two
cluster assignments restricted to the `n_sites_common` common sites. `1`
= perfect agreement; `0` = chance agreement; negative = worse than
chance.
</td>

</tr>

</tbody>

</table>

## tab_05_snp_enrichment.csv

**Granularity:** one row per `(cell, metric)`. **Row count:** 216.
**Purpose:** Fisher exact test of whether the variant’s defining- SNP
sites cluster together more than would be expected by chance (Q4). The
Bonferroni correction uses `N_TOTAL_ENRICHMENT_TESTS = 220` (the design
ceiling).

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Column dictionary for `tab_05_snp_enrichment.csv`.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Column
</th>

<th style="text-align:left;">

Type
</th>

<th style="text-align:left;">

Description
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

dataset
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Dataset name: `NCBI_US` or `GISAID_US`.
</td>

</tr>

<tr>

<td style="text-align:left;">

variant
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

WHO variant label.
</td>

</tr>

<tr>

<td style="text-align:left;">

strategy
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Windowing strategy (`sliding_2m`, `disjoint_2m`, or `cumulative_1m`).
</td>

</tr>

<tr>

<td style="text-align:left;">

metric
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Curve type clustered (see `tab_01$metric` for definitions).
</td>

</tr>

<tr>

<td style="text-align:left;">

cell_status
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Cell-level terminal status (same domain as `tab_01$cell_status`).
</td>

</tr>

<tr>

<td style="text-align:left;">

enrichment_status
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Per-metric enrichment status: `ok` (test ran), `no_snp_in_class1` (no
defining-SNP site landed in the class-1 set), `fisher_failed`
(degenerate contingency table), `skipped` (cell not `ok` or metric not
`ok`), …
</td>

</tr>

<tr>

<td style="text-align:left;">

n_snp_in_class1
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Count of defining-SNP sites that appear in this cell’s class-1 site set
(subset of the variant’s `Defining_SNP_Sites`).
</td>

</tr>

<tr>

<td style="text-align:left;">

n_class1_total
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Total number of sites in the cell’s class-1 site set (matches
`tab_01$n_class1_sites`).
</td>

</tr>

<tr>

<td style="text-align:left;">

n_snp_total_truth
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Total number of defining-SNP sites for the variant in the ground-truth
catalogue (`length(sarscov2_variants$Defining_SNP_Sites[[variant]])`).
</td>

</tr>

<tr>

<td style="text-align:left;">

fisher_p
</td>

<td style="text-align:left;">

num
</td>

<td style="text-align:left;">

One-tailed Fisher exact `p`-value for the most-enriched cluster, testing
the null that defining-SNP sites are distributed across clusters in
proportion to cluster size.
</td>

</tr>

<tr>

<td style="text-align:left;">

fisher_p_bonferroni
</td>

<td style="text-align:left;">

num
</td>

<td style="text-align:left;">

Bonferroni-adjusted `p`: `min(1, fisher_p * N_TOTAL_ENRICHMENT_TESTS)`
with `N_TOTAL_ENRICHMENT_TESTS = 220` (the design ceiling on the number
of enrichment tests across the study).
</td>

</tr>

<tr>

<td style="text-align:left;">

max_overlap_cluster
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Cluster ID with the largest defining-SNP overlap; `NA` if no overlap.
</td>

</tr>

<tr>

<td style="text-align:left;">

max_overlap_OR
</td>

<td style="text-align:left;">

num
</td>

<td style="text-align:left;">

Odds ratio for the most-enriched cluster’s 2x2 contingency table (`Inf`
when one off-diagonal cell is zero, i.e., perfect or near-perfect
overlap).
</td>

</tr>

</tbody>

</table>

## tab_06_per_frame.csv

**Granularity:** one row per `(cell, metric, frame)`, restricted to
cells with `status == "ok"`. **Row count:** 1,764. **Purpose:** the
per-frame diagnostics that drive the animations — the table to consult
when interpreting whether a specific frame in a GIF crossed the
convergence threshold or was rendered null by one of the three quality
guards.

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Column dictionary for `tab_06_per_frame.csv`.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Column
</th>

<th style="text-align:left;">

Type
</th>

<th style="text-align:left;">

Description
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

dataset
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Dataset name: `NCBI_US` or `GISAID_US`.
</td>

</tr>

<tr>

<td style="text-align:left;">

variant
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

WHO variant label.
</td>

</tr>

<tr>

<td style="text-align:left;">

strategy
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Windowing strategy (`sliding_2m`, `disjoint_2m`, or `cumulative_1m`).
</td>

</tr>

<tr>

<td style="text-align:left;">

metric
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Curve type clustered (see `tab_01$metric` for definitions).
</td>

</tr>

<tr>

<td style="text-align:left;">

frame_idx
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Sequential frame index (1-based) within this `(cell, metric)`. Frames
are emitted for every interior partition with at least
`FITTING_WINDOW_HALFWIDTH = 2` partitions on either side.
</td>

</tr>

<tr>

<td style="text-align:left;">

centre_bin
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Index of the partition (bin) at the centre of this frame’s fitting
sub-window. Used as the frame’s identifier in the GIF and as the
time-axis position of the convergence indices in `tab_01`.
</td>

</tr>

<tr>

<td style="text-align:left;">

fitting_start_bin
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Index of the first partition in the frame’s 5-partition fitting
sub-window (`centre_bin - 2`).
</td>

</tr>

<tr>

<td style="text-align:left;">

fitting_end_bin
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Index of the last partition in the frame’s 5-partition fitting
sub-window (`centre_bin + 2`).
</td>

</tr>

<tr>

<td style="text-align:left;">

gmm_status
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Status of the per-frame GMM refit: `ok`, `G_eq_1` (Guard A fails: GMM
selects a single component, so no class-1 set exists), `failed`
(numerical failure in `mclust::Mclust`), or `skipped`.
</td>

</tr>

<tr>

<td style="text-align:left;">

gmm_G
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of GMM components selected by BIC over `G` in 1, …, 15 for the
frame’s pooled-sequence entropy.
</td>

</tr>

<tr>

<td style="text-align:left;">

gmm_n_sequences
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of sequences pooled across the frame’s 5 partitions (the basis
for the per-frame entropy and GMM).
</td>

</tr>

<tr>

<td style="text-align:left;">

fit_status
</td>

<td style="text-align:left;">

chr
</td>

<td style="text-align:left;">

Status of the per-frame `fdacluster::fdahclust` refit (`ok`,
`skipped_too_few_sites`, `error`, …).
</td>

</tr>

<tr>

<td style="text-align:left;">

n_class1_sites
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of class-1 sites in this frame’s per-frame GMM (not the
cell-level class-1 set; per-frame GMMs are independent of the cell-level
fit).
</td>

</tr>

<tr>

<td style="text-align:left;">

n_clusters
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of clusters selected for this frame’s `fdahclust` fit by
silhouette over `k` in 2, …, 6.
</td>

</tr>

<tr>

<td style="text-align:left;">

silhouette_mean
</td>

<td style="text-align:left;">

num
</td>

<td style="text-align:left;">

Mean silhouette width at the frame’s chosen `k`.
</td>

</tr>

<tr>

<td style="text-align:left;">

structurally_null
</td>

<td style="text-align:left;">

lgl
</td>

<td style="text-align:left;">

Guard B flag. `TRUE` iff
`silhouette_mean &lt; MIN_MEAN_SILHOUETTE = 0.25` (Kaufman & Rousseeuw
1990 “reasonable structure” floor). Frames with
`structurally_null = TRUE` cannot satisfy any convergence predicate
regardless of cluster pattern.
</td>

</tr>

<tr>

<td style="text-align:left;">

n_snp_in_class1
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of defining-SNP sites that are in this frame’s class-1 site set.
</td>

</tr>

<tr>

<td style="text-align:left;">

n_mut_in_class1
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of defining-mutation sites that are in this frame’s class-1 site
set.
</td>

</tr>

<tr>

<td style="text-align:left;">

converged_snp
</td>

<td style="text-align:left;">

lgl
</td>

<td style="text-align:left;">

`TRUE` iff all SNP sites in this frame’s class-1 share a single cluster
whose size does not exceed `MAX_CLUSTER_SIZE = 33` (i.e., not an
“oversized” cluster).
</td>

</tr>

<tr>

<td style="text-align:left;">

converged_mut
</td>

<td style="text-align:left;">

lgl
</td>

<td style="text-align:left;">

`TRUE` iff the same holds for the variant’s broader defining-mutation
sites.
</td>

</tr>

<tr>

<td style="text-align:left;">

converged_class1
</td>

<td style="text-align:left;">

lgl
</td>

<td style="text-align:left;">

`TRUE` iff all class-1 sites with positive per-site silhouette (Guard C:
well-clustered sites) share a single non-oversized cluster.
</td>

</tr>

<tr>

<td style="text-align:left;">

qualifying_cluster_snp
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Cluster ID into which the SNP sites converged for this frame; `NA` if
`converged_snp = FALSE`.
</td>

</tr>

<tr>

<td style="text-align:left;">

qualifying_cluster_mut
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Cluster ID into which the mutation sites converged; `NA` if
`converged_mut = FALSE`.
</td>

</tr>

<tr>

<td style="text-align:left;">

qualifying_cluster_class1
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Cluster ID into which the class-1 sites converged; `NA` if
`converged_class1 = FALSE`.
</td>

</tr>

<tr>

<td style="text-align:left;">

anchor_idx
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Hellinger anchor partition index used to compute this frame’s curves
(`bin_idx` of `T1` or `T_pred`). `NA` for the `entropy` metric, which
has no anchor.
</td>

</tr>

<tr>

<td style="text-align:left;">

n_oversized_clusters
</td>

<td style="text-align:left;">

int
</td>

<td style="text-align:left;">

Number of clusters in this frame whose size exceeds
`MAX_CLUSTER_SIZE = 33` (Omicron’s defining-mutation count, the largest
among the 12 WHO variants). Oversized clusters are rendered light grey
in the GIF panels and are excluded from convergence predicates.
</td>

</tr>

</tbody>

</table>

# Limitations

1.  **17 underpowered cells.** Cells in which any bin contained fewer
    than `MIN_SEQUENCES_PER_BIN = 30L` sequences are excluded from
    clustering. This concentrates on the four low-prevalence VOIs
    (Epsilon, Eta, Iota, Zeta), where a single early-pandemic bin
    regularly contained 2–8 sequences despite the variant’s overall
    window count exceeding 20,000. One Alpha × cumulative_1m × GISAID US
    cell also crossed the threshold (one of 20 bins had 6 sequences).
    Lowering `MIN_SEQUENCES_PER_BIN` would recover these cells at the
    cost of less reliable per-bin frequency estimates; a sensitivity
    sweep across thresholds is a natural follow-up.
2.  **3 out-of-coverage cells.** Omicron × NCBI US is structurally
    excluded across all three strategies because NCBI US coverage ends
    in September 2021 and Omicron was first detected in the US in
    December 2021. Mitigated by the parallel run on GISAID US, which
    covers Omicron in all three strategies.
3.  **US-only, Spike-only.** Other geographies and other SARS-CoV-2
    proteins are not analysed. Generalisation outside Spike is not
    tested.
4.  **Window-length fixed at $\pm 10$ months.** The choice $K = 5$ bins
    at $b_{\max} = 2$ months is informed by the earlier “5 points
    before, 5 after” convention but is not itself a parameter sweep. A
    sensitivity sweep on window length is a natural follow-up.
5.  **Three-strategy panel non-exhaustive.** Other discretisations exist
    (3-month disjoint, 1-month sliding) and were considered but excluded
    to keep the cell count tractable.
6.  **No explicit “true” cluster partition.** Functional clustering has
    no ground truth in the strict sense of change-point detection (where
    the truth set is the WHO catalogue’s detection dates). The SNP
    enrichment test (Q4) is a biological proxy — a *plausibility* check,
    not a *validation* of cluster correctness. Captions and discussion
    phrase this distinction explicitly.
7.  **Hellinger reference is dual-anchored, not swept.** The study
    anchors the Hellinger trajectories at the window-start partition
    $T_1$ *and* at the partition containing the predecessor variant’s US
    detection ($T_{\text{pred}}$). The map of predecessors is curated
    (Alpha → D614G; Beta / Gamma / Theta / Kappa / Lambda → Alpha; Delta
    → Beta; Omicron → Delta; Epsilon / Zeta / Eta / Iota → D614G).
    Alternative ancestry mappings (e.g., by Pango parent rather than by
    WHO label) are not swept and would be a natural follow-up.
8.  **Class-1 site reduction is a sequential bottleneck.** The full-
    window cell-level fits depend on the GMM having identified the
    correct high-entropy component (class 1). A small high-entropy class
    1 may exclude positions whose entropy is just below the class
    boundary, and a misidentified class 1 propagates to all downstream
    metrics. The per-frame design partly mitigates this by re-fitting
    the GMM in every fitting sub-window, so a class-1 misidentification
    at a single moment in time does not corrupt the neighbouring frames.
    A formal sensitivity analysis that takes the union of class-1 and
    class-2 sites is a supplementary check reserved for follow-up work.
9.  **`fdahclust` warping disabled.** Curve alignment (DTW-style
    warping) is not used. The justification is methodological: the
    x-axis is bin index, *which at the population scale used here is
    calendar time*. Warping the calendar would correspond to claiming
    that the biological clock of the virus ran faster or slower in some
    bins than in others — a claim with no clear evidentiary basis at the
    level of national surveillance pools used in this study.
10. **Fisher exact assumes site independence under linkage
    disequilibrium.** The Fisher exact SNP enrichment test treats
    cluster assignment of each site as independent of every other site,
    which is violated by linkage disequilibrium across Spike positions
    co-segregating in real lineages. The reported $p$-values should
    therefore be read as descriptive concentration statistics rather
    than as strict null-hypothesis tests; a permutation-style
    sensitivity analysis that resamples cluster labels while preserving
    the SNP-site count is planned as a supplementary check.
11. **Bin-width robustness untested.** Bin widths are fixed by strategy.
    A sensitivity sweep on bin width within a strategy is a natural
    follow-up (mirroring the planned bin-width sweep on
    `Changepoint_Detection_Study`).
12. **GIF frame rate fixed at 1 fps.** Visual presentation of animations
    is a deliverable, not a tuning target; the frame rate was chosen to
    be slow enough that the eye can track cluster re-organisation across
    frames.
13. **N_WORKERS sensitivity to memory pressure.** The first-pass run at
    `N_WORKERS = 12` on the 125 GB shared server lost five cells to the
    kernel OOM-killer (Omicron × cumulative_1m, Alpha × cumulative_1m,
    Beta / Kappa / Lambda × sliding_2m). These were recovered cleanly at
    `N_WORKERS = 5`. The Hellinger pairwise computation for the heaviest
    cells (longest cumulative strategies on GISAID US) is the dominant
    allocation; users on smaller hosts should use a conservative
    `N_WORKERS` setting.

# ADEMP summary

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

ADEMP framing (Morris, White, and Crowther 2019). Aim, Data-generating
mechanism, Estimand, Methods, Performance measures.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Component
</th>

<th style="text-align:left;">

Description
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Aim
</td>

<td style="text-align:left;">

Cluster per-site Shannon entropy curves and pairwise Hellinger
trajectories (dual-anchored at the window start T1 and at the partition
containing the predecessor variant’s US detection Tpred) from real
SARS-CoV-2 Spike surveillance data, around the US first-detection dates
of twelve WHO-labelled variants, under three temporal-partitioning
strategies, and quantify the internal quality, resampling stability,
cross-strategy concordance, and biological plausibility of the resulting
clusters at both the cell level (full-window fit) and per frame
(five-partition sliding fits).
</td>

</tr>

<tr>

<td style="text-align:left;">

Data-generating mechanism
</td>

<td style="text-align:left;">

Per-variant 20-month windows of two real surveillance datasets (NCBI US
and GISAID US), partitioned by three strategies, with GMM-based class-1
site selection per window (cell-level) and per five-partition fitting
slice (per-frame), and no oracle leakage across cells.
</td>

</tr>

<tr>

<td style="text-align:left;">

Estimand
</td>

<td style="text-align:left;">

Per-(variant, dataset, strategy, metric) cluster partition of class-1
Spike sites at the cell level and at every interior fitting frame,
characterised by the chosen number of clusters, silhouette width,
per-cluster bootstrap Jaccard stability, cross-strategy Adjusted Rand
Index, Fisher exact SNP enrichment p-value, and the first-frame
convergence indices for the three convergence predicates (SNP / Mutation
/ class-1).
</td>

</tr>

<tr>

<td style="text-align:left;">

Methods
</td>

<td style="text-align:left;">

fdacluster::fdahclust under warping_class = ‘none’, L2 metric,
centroid_type = ‘mean’, linkage = ‘complete’, k chosen by silhouette
over k in {2, …, 6}, log1p input transformation; mclust::Mclust over G
in {1, …, 15} for class-1 site selection at both cell and frame levels.
</td>

</tr>

<tr>

<td style="text-align:left;">

Performance measures
</td>

<td style="text-align:left;">

Average silhouette width (internal quality), bootstrap Jaccard stability
with B = 100 (robustness), cross-strategy ARI (methodological
concordance), Fisher exact SNP enrichment p-value with Bonferroni
correction across N_TOTAL_ENRICHMENT_TESTS = 220 (biological
plausibility), and the first-frame convergence indices per metric
(longitudinal interpretability).
</td>

</tr>

</tbody>

</table>

# Invocation

``` bash
# Working directory
cd analysis/FDA_Analysis

# Optional: pre-flight parse check (no model evaluation).
Rscript -e 'for (f in c("setup.R","helpers_fda.R","helpers_frames.R",
                         "helpers_panels.R","run_one_cell.R","fda_analysis.R",
                         "plot_results.R","animate_results.R",
                         "aggregate_tables.R","render_existing.R"))
              parse(file=f); cat("parse OK\n")'

# Phase 1: cell production (orchestrator). Recommended on a 64+ GB host:
TS=$(date +%Y%m%d_%H%M%S)
LOG="outputs/logs/orchestrator_${TS}.log"
N_WORKERS=8 nohup Rscript fda_analysis.R > "$LOG" 2>&1 &
echo $! > outputs/logs/orchestrator.pid
disown

# Phase 2-5 (after Phase 1 completes):
Rscript aggregate_tables.R   2>&1 | tee outputs/logs/aggregate_${TS}.log
Rscript plot_results.R       2>&1 | tee outputs/logs/plot_${TS}.log

# render_existing.R writes per-cell cluster_membership_*.csv files
# (parallelisable with xargs -P; see analysis/FDA_Analysis/render_existing.R).
Rscript -e '
  cells <- list.files("outputs/cells", pattern = "\\.rds$", full.names = TRUE)
  rows <- lapply(cells, function(p) {
    r <- readRDS(p)
    if (!is.null(r$status) && r$status == "ok")
      sprintf("%s %s %s", r$dataset, r$variant, r$strategy)
    else NULL
  })
  writeLines(unlist(rows), "ok_cells.txt")'
cat ok_cells.txt | xargs -P 8 -L 1 bash -c '
  SMOKE_DATASET=$0 SMOKE_VARIANT=$1 SMOKE_STRATEGY=$2 \
    Rscript render_existing.R \
    > "outputs/logs/render_${0}__${1}__${2}.log" 2>&1
  echo "[$(date +%H:%M:%S)] $0 $1 $2 exit=$?"'

nohup Rscript animate_results.R > outputs/logs/animate_${TS}.log 2>&1 &
```

Single-cell debug:

``` r
setwd("analysis/FDA_Analysis")
source("setup.R")
source("helpers_fda.R")
source("helpers_frames.R")
source("run_one_cell.R")

config <- build_config()
truth  <- load_truth_catalogue()
fm     <- load_feature_matrix(config$DATASETS$GISAID_US$feature_rds,
                                expected_n_sites = config$N_SITES)
result <- run_one_cell(
  dataset_name  = "GISAID_US",
  variant_name  = "Delta",
  strategy_name = "disjoint_2m",
  feature_matrix = fm,
  truth_df       = truth,
  config         = config
)
str(result, max.level = 2L)
```

# Project structure

    analysis/FDA_Analysis/
    |-- README.Rmd                          (this file)
    |-- README.md                           (knitted GitHub output)
    |-- README.html                         (knitted HTML output)
    |-- references_fda.bib                  (BibTeX source)
    |-- captions_and_discussion.txt         (figure captions + main-text discussion)
    |
    |-- setup.R                             (config, helpers, paths, seeds, predecessor map)
    |-- helpers_fda.R                       (per-cell window slicing, GMM, entropy, Hellinger, fdahclust)
    |-- helpers_frames.R                    (per-frame GMM, fdahclust, silhouette, convergence)
    |-- helpers_panels.R                    (dual-panel rendering for static plots and GIF frames)
    |-- run_one_cell.R                      (per-cell driver, including per-frame loop)
    |-- fda_analysis.R                      (orchestrator)
    |
    |-- plot_results.R                      (per-frame PNGs + per-cell summary PNGs)
    |-- animate_results.R                   (per-cell × per-metric GIF animations)
    |-- aggregate_tables.R                  (6 supplementary CSV tables)
    |-- render_existing.R                   (per-cell cluster_membership CSV side-output)
    |-- smoke_test.R                        (single-cell smoke test)
    |
    `-- outputs/
        |-- cells/cell_<dataset>__<variant>__<strategy>.rds   (72 per-cell result files)
        |-- logs/
        |   |-- orchestrator_<timestamp>.log
        |   |-- aggregate_<timestamp>.log
        |   |-- plot_<timestamp>.log
        |   |-- animate_<timestamp>.log
        |   |-- render_*.log                  (per-cell render_existing.R logs)
        |   `-- error_log.txt                 (per-cell errors, normally empty)
        |
        |-- plots/<dataset>/<variant>/<strategy>/
        |   |-- frame_entropy_<centre:03d>.png            (one per frame)
        |   |-- frame_hellinger_T1_<centre:03d>.png       (one per frame)
        |   |-- frame_hellinger_Tpred_<centre:03d>.png    (one per frame)
        |   |-- summary_entropy.png                       (full-window cell fit)
        |   |-- summary_hellinger_T1.png
        |   |-- summary_hellinger_Tpred.png
        |   |-- frames.manifest.csv                       (frame inventory per cell)
        |   |-- cluster_membership_entropy.csv            (rendered cluster attribution per site)
        |   |-- cluster_membership_hellinger_T1.csv
        |   `-- cluster_membership_hellinger_Tpred.csv
        |
        |-- gifs/<dataset>/<variant>/<strategy>/
        |   |-- anim_entropy.gif
        |   |-- anim_hellinger_T1.gif
        |   `-- anim_hellinger_Tpred.gif
        |
        `-- tables/
            |-- tab_01_internal_quality.csv     (216 rows; per cell × metric)
            |-- tab_02_cluster_membership.csv   (5,058 rows; per cell × metric × site)
            |-- tab_03_stability.csv            (469 rows; per cell × metric × cluster)
            |-- tab_04_cross_strategy_ari.csv   (141 rows; per dataset × variant × metric × strategy pair)
            |-- tab_05_snp_enrichment.csv       (216 rows; per cell × metric)
            `-- tab_06_per_frame.csv            (1,764 rows; per cell × metric × frame)

Each per-cell RDS holds the cell identifiers, window metadata, GMM fit
summary, the entropy and Hellinger matrices on the class-1 site set
(both raw and log1p-transformed, both anchors), the three
`fdacluster::fdahclust` fits with their chosen number of clusters and
silhouette values, the three per-cluster bootstrap Jaccard stability
vectors, the three per-cluster SNP overlap counts and Fisher test
statistics, the per-frame records for all three metrics, the three
convergence summaries, run metadata (status, wall time, seeds), and R /
`ViralEntropR` version strings. The six supplementary CSVs are computed
by `aggregate_tables.R` from these per-cell files.

# Session information

    #> R version 4.5.2 (2025-10-31 ucrt)
    #> Platform: x86_64-w64-mingw32/x64
    #> Running under: Windows 11 x64 (build 26200)
    #> 
    #> Matrix products: default
    #>   LAPACK version 3.12.1
    #> 
    #> locale:
    #> [1] LC_COLLATE=English_Canada.utf8  LC_CTYPE=English_Canada.utf8   
    #> [3] LC_MONETARY=English_Canada.utf8 LC_NUMERIC=C                   
    #> [5] LC_TIME=English_Canada.utf8    
    #> 
    #> time zone: America/Toronto
    #> tzcode source: internal
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] kableExtra_1.4.0
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] vctrs_0.7.1        svglite_2.2.2      cli_3.6.5          knitr_1.51        
    #>  [5] rlang_1.1.7        xfun_0.56          stringi_1.8.7      otel_0.2.0        
    #>  [9] textshaping_1.0.4  glue_1.8.0         htmltools_0.5.9    scales_1.4.0      
    #> [13] rmarkdown_2.30     evaluate_1.0.5     fastmap_1.2.0      yaml_2.3.12       
    #> [17] lifecycle_1.0.5    stringr_1.6.0      compiler_4.5.2     RColorBrewer_1.1-3
    #> [21] rstudioapi_0.18.0  systemfonts_1.3.1  farver_2.1.2       digest_0.6.39     
    #> [25] viridisLite_0.4.3  R6_2.6.1           magrittr_2.0.4     tools_4.5.2       
    #> [29] xml2_1.5.1

# References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Agresti1992" class="csl-entry">

Agresti, Alan. 1992. “A Survey of Exact Inference for Contingency
Tables.” *Statistical Science* 7 (1): 131–53.
<https://doi.org/10.1214/ss/1177011454>.

</div>

<div id="ref-BoxCox1964" class="csl-entry">

Box, G. E. P., and D. R. Cox. 1964. “An Analysis of Transformations.”
*Journal of the Royal Statistical Society, Series B (Methodological)* 26
(2): 211–43. <https://doi.org/10.1111/j.2517-6161.1964.tb00553.x>.

</div>

<div id="ref-Hennig2007" class="csl-entry">

Hennig, Christian. 2007. “Cluster-Wise Assessment of Cluster Stability.”
*Computational Statistics & Data Analysis* 52 (1): 258–71.
<https://doi.org/10.1016/j.csda.2006.11.025>.

</div>

<div id="ref-HubertArabie1985" class="csl-entry">

Hubert, Lawrence, and Phipps Arabie. 1985. “Comparing Partitions.”
*Journal of Classification* 2 (1): 193–218.
<https://doi.org/10.1007/BF01908075>.

</div>

<div id="ref-KaufmanRousseeuw1990" class="csl-entry">

Kaufman, Leonard, and Peter J. Rousseeuw. 1990. *Finding Groups in Data:
An Introduction to Cluster Analysis*. New York: Wiley.
<https://doi.org/10.1002/9780470316801>.

</div>

<div id="ref-KokoszkaReimherr2017" class="csl-entry">

Kokoszka, Piotr, and Matthew Reimherr. 2017. *Introduction to Functional
Data Analysis*. Texts in Statistical Science. Boca Raton, FL: Chapman &
Hall/CRC. <https://doi.org/10.1201/9781315117416>.

</div>

<div id="ref-Korber2020" class="csl-entry">

Korber, Bette, Will M. Fischer, Sandrasegaram Gnanakaran, Hyejin Yoon,
James Theiler, Werner Abfalterer, Nick Hengartner, et al. 2020.
“Tracking Changes in SARS-CoV-2 Spike: Evidence That D614G Increases
Infectivity of the COVID-19 Virus.” *Cell* 182 (4): 812–827.e19.
<https://doi.org/10.1016/j.cell.2020.06.043>.

</div>

<div id="ref-Morris2019" class="csl-entry">

Morris, Tim P., Ian R. White, and Michael J. Crowther. 2019. “Using
Simulation Studies to Evaluate Statistical Methods.” *Statistics in
Medicine* 38 (11): 2074–2102. <https://doi.org/10.1002/sim.8086>.

</div>

<div id="ref-RamsaySilverman2005" class="csl-entry">

Ramsay, J. O., and B. W. Silverman. 2005. *Functional Data Analysis*.
2nd ed. Springer Series in Statistics. New York: Springer.

</div>

<div id="ref-Rousseeuw1987" class="csl-entry">

Rousseeuw, Peter J. 1987. “Silhouettes: A Graphical Aid to the
Interpretation and Validation of Cluster Analysis.” *Journal of
Computational and Applied Mathematics* 20: 53–65.
<https://doi.org/10.1016/0377-0427(87)90125-7>.

</div>

<div id="ref-Sangalli2010" class="csl-entry">

Sangalli, Laura M., Piercesare Secchi, Simone Vantini, and Valeria
Vitelli. 2010. “K-Mean Alignment for Curve Clustering.” *Computational
Statistics & Data Analysis* 54 (5): 1219–33.
<https://doi.org/10.1016/j.csda.2009.12.008>.

</div>

<div id="ref-Sayers2022" class="csl-entry">

Sayers, Eric W., Evan E. Bolton, J. Rodney Brister, Kathi Canese,
Jessica Chan, Donald C. Comeau, Ryan Connor, et al. 2022. “Database
Resources of the National Center for Biotechnology Information.”
*Nucleic Acids Research* 50 (D1): D20–26.
<https://doi.org/10.1093/nar/gkab1112>.

</div>

<div id="ref-Scrucca2016" class="csl-entry">

Scrucca, Luca, Michael Fop, T. Brendan Murphy, and Adrian E. Raftery.
2016. “Mclust 5: Clustering, Classification and Density Estimation Using
Gaussian Finite Mixture Models.” *The R Journal* 8 (1): 289–317.
<https://doi.org/10.32614/RJ-2016-021>.

</div>

<div id="ref-ShuMcCauley2017" class="csl-entry">

Shu, Yuelong, and John McCauley. 2017. “GISAID: Global Initiative on
Sharing All Influenza Data — from Vision to Reality.” *Eurosurveillance*
22 (13): 30494. <https://doi.org/10.2807/1560-7917.ES.2017.22.13.30494>.

</div>

<div id="ref-vanDerVaart1998" class="csl-entry">

Vaart, Aad W. van der. 1998. *Asymptotic Statistics*. Cambridge Series
in Statistical and Probabilistic Mathematics. Cambridge: Cambridge
University Press. <https://doi.org/10.1017/CBO9780511802256>.

</div>

</div>
