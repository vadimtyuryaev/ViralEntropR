ViralEntropR: A Stratified Bootstrap Benchmark of Change-Point Detection
Methods on SARS-CoV-2 Spike Hellinger Trajectories
================
Vadim Tyuryaev, Jane Heffernan, Hanna Jankowski
May 20, 2026

- [Overview](#overview)
- [Terminology](#terminology)
- [Research questions](#research-questions)
- [Datasets](#datasets)
- [Ground truth](#ground-truth)
- [Methodology](#methodology)
  - [Stratified bootstrap design](#stratified-bootstrap-design)
  - [Bin alignment](#bin-alignment)
  - [Hellinger trajectory](#hellinger-trajectory)
  - [Per-window GMM dimensionality reduction
    (secondary)](#per-window-gmm-dimensionality-reduction-secondary)
  - [The four change-point detection
    methods](#the-four-change-point-detection-methods)
    - [Output harmonisation](#output-harmonisation)
  - [Truth-to-detection matching](#truth-to-detection-matching)
  - [Metrics](#metrics)
  - [Sensitivity analyses](#sensitivity-analyses)
- [Computational budget and runtime](#computational-budget-and-runtime)
- [Reproducibility](#reproducibility)
- [Implementation summary](#implementation-summary)
- [Per-replicate output schema](#per-replicate-output-schema)
- [Limitations](#limitations)
- [ADEMP summary](#ademp-summary)
- [Invocation](#invocation)
- [Project structure](#project-structure)
- [Session information](#session-information)
- [References](#references)

# Overview

This study benchmarks four high-dimensional change-point detection
methods against known SARS-CoV-2 variant emergence dates on two
independent surveillance datasets. The four methods span two algorithmic
families: energy-statistic methods from the **`ecp`** package (Matteson
and James 2014; James and Matteson 2014) (`ks.cp3o` and `e.agglo`) and
high-dimensional binary segmentation methods from the **`HDcpDetect`**
package (Fryzlewicz 2014; Olshen et al. 2004) (`binary.segmentation` and
`wild.binary.segmentation`). All four operate on multivariate
Hellinger-distance trajectories derived from per-site amino acid
distributions in the SARS-CoV-2 Spike protein.

Detection is evaluated against the WHO-labelled variant-of-concern (VOC)
and variant-of-interest (VOI) emergence dates in
`sarscov2_variants$Date_First_Detected_US`. The primary endpoint is the
per-window precision, recall, F1 score, and temporal localisation error
(TLE) (Truong, Oudre, and Vayatis 2020; Aminikhanghahi and Cook 2017),
with the headline question being which method, if any, dominates across
the operating regime of real-world genomic surveillance.

A secondary endpoint addresses whether GMM-based entropy dimensionality
reduction — the central feature-selection step of the ViralEntropR
pipeline — preserves change-point detectability when applied per-window
(i.e., without oracle knowledge of the full surveillance timeline).

# Terminology

The following terms appear throughout the document. Definitions are
written for readers approaching the work from either the statistical or
virological side; cross-references to standard sources are given where
applicable.

| Term | Definition |
|----|----|
| Spike protein | The 1,273-amino-acid surface glycoprotein of SARS-CoV-2. The primary target of immune surveillance and of this study. Reference: GenBank accession YP_009724390. |
| Variant (VOC / VOI) | A SARS-CoV-2 lineage carrying defining amino-acid substitutions and classified by the WHO as a Variant of Concern or Variant of Interest. The `sarscov2_variants` catalogue contains 12 such variants with US first-detection dates. |
| Surveillance dataset | A real-world feature matrix of dated Spike-protein sequences. This study uses NCBI US (Sayers et al. 2022) and GISAID US (Shu and McCauley 2017). |
| Feature matrix | A data.frame with one row per dated sequence, one column for the collection date, and 1,273 integer-coded columns for the 1,273 Spike positions (25-residue alphabet: 20 standard amino acids, B, Z, X, \*, gap). |
| Window (sub-window) | A contiguous date range within the dataset, with start and end on first-of-month boundaries. The bootstrap design draws 5,001 random windows per dataset. |
| Bin | A 2-month non-overlapping partition of a window. A window of $L$ months contains $L/2$ bins, denoted $T_1, T_2, \ldots, T_{L/2}$. |
| Reference bin ($T_1$) | The first bin of a window. Its per-site amino-acid distribution serves as the reference against which subsequent bins are compared. Truth events falling in $T_1$ are dropped because no Hellinger value exists for the reference bin itself. |
| Hellinger trajectory | A matrix of pairwise Hellinger distances from $T_1$ to each subsequent bin $T_2, T_3, \ldots, T_{L/2}$, computed per Spike position. Shape: $n_\text{sites} \times (L/2 - 1)$. Methods see the transposed shape: $(L/2 - 1) \times n_\text{sites}$. |
| Truth event | A first-US-detection date from the `sarscov2_variants` catalogue intersected with the window’s date range, mapped to a Hellinger row by the bin-shift convention (variant in bin $T_k \rightarrow$ Hellinger row $k - 1$). |
| Detection | An interior change-point index returned by one of the four methods, in the harmonised convention “first row of the new segment”. |
| Tier | Stratification of windows by length: **Short** (6, 8, or 10 months), **Medium** (12, 14, or 16 months), **Long** (18, 20, 22, 24, 26, or 28 months — dataset-capped). Each tier receives 1,667 accepted windows per dataset. |
| Stratified bootstrap | The sampling scheme: per dataset, per tier, draw 1,667 candidate windows uniformly at random over feasible (length, start) combinations, accept those passing the bin-count and truth-bin thresholds. |
| GMM (Gaussian Mixture Model) | A clustering method that fits a weighted sum of Gaussian components to a univariate sample. This study uses `mclust::Mclust` (Scrucca et al. 2016) via `ViralEntropR::cluster_sites_by_entropy` to cluster per-site Shannon entropies into $G \in \{1, \ldots, 15\}$ components and select the best fit by BIC. |
| Class 1 (GMM-reduced sites) | After `relabel_entropy_classes`, the highest-entropy GMM component. The secondary analysis recomputes Hellinger trajectories using only class-1 sites; the primary analysis uses all 1,273. |
| `ecp` | The CRAN package implementing energy-distance-based change-point detection (James and Matteson 2014). |
| `ks.cp3o` (ECP-1) | An `ecp` method that finds at most $K$ change points by dynamic programming over the energy-distance objective. $K$ is a user-supplied upper bound; this study uses both a data-driven dynamic-$K$ (= number of Hellinger rows minus 2) and a sweep over fixed $K \in \{1, 2, 3, 5, 10\}$. |
| `e.agglo` (ECP-2) | A $K$-free `ecp` method that begins with each time point in its own segment and agglomeratively merges adjacent segments to maximise the within-segment energy distance. Configured with `alpha = 1` and `penalty = function(cps) 0` (no segment-count penalty). |
| `HDcpDetect` | The CRAN package implementing high-dimensional binary segmentation and wild binary segmentation with permutation-based or asymptotic significance testing (Fryzlewicz 2014). |
| `binary.segmentation` (HDcp-BS) | Recursive binary partitioning: at each step, find the single split point maximising a high-dimensional CUSUM statistic; test for significance; recurse on each daughter segment if significant. Genealogy: Olshen et al. (2004) in the genomics context; extended to high-dimensional theory in Fryzlewicz (2014). |
| `wild.binary.segmentation` (HDcp-WBS) | A randomised variant of binary segmentation that draws $M$ random sub-intervals at each recursive step and selects the candidate with the maximum CUSUM statistic across them. This study uses $M = 100$. |
| $K$ | The change-point-count budget for `ks.cp3o`. With `minsize = 2`, the algorithm requires the input series to have at least $(K + 1) \times 2$ rows. |
| $M$ | The number of random sub-intervals drawn at each recursive step of WBS (Fryzlewicz 2014). |
| `minsize` | The minimum segment length permitted by `ks.cp3o`. Fixed at 2 (one bin) throughout this study. The methods in the `HDcpDetect` family inherit an analogous internal minimum that produces NA-propagation failures on very short trajectories. |
| Bin-shift convention | A variant emerging in partition $T_k$ appears in the Hellinger trajectory at row index $k - 1$, because $T_1$ is the reference and is dropped. The truth-to-bin mapping applies this $-1$ shift before matching. |
| Truth-detection matching | Truth-first greedy with consumption: for each unique truth bin (in temporal order), claim the nearest available detection within $\pm 1$ bin tolerance; mark consumed detections; remaining unmatched detections are FPs, remaining unclaimed truths are FNs. |
| TP / FP / FN | Counts produced by the matching: true positive (a truth claimed a detection), false positive (a detection went unclaimed), false negative (a truth went unclaimed). |
| Precision / Recall / F1 | Standard derivations: $P = \text{TP} / (\text{TP} + \text{FP})$, $R = \text{TP} / (\text{TP} + \text{FN})$, $F_1 = 2 P R / (P + R)$. Degenerate-case conventions follow Truong, Oudre, and Vayatis (2020). |
| TLE (Temporal Localisation Error) | Mean absolute bin-distance between matched (truth, detection) pairs. One bin = 2 months. NA when no matches. |
| K-sweep | The sensitivity sweep for `ks.cp3o`: per replicate, the algorithm runs at the dynamic $K$ and at fixed $K \in \{1, 2, 3, 5, 10\}$. Saturation patterns inform the operationally minimum $K$. |
| Truth-shift sweep | A separate sensitivity sweep: the truth set is shifted by $\{-2, -1, 0, +1, +2\}$ bins (= $\{-4, -2, 0, +2, +4\}$ months) and matching is re-run. Detection sets are unchanged across shifts; only the bin-mapped truth indices move. |
| Replicate | One per-window detection result, indexed by `(dataset, tier, run_id)`. Each per-replicate RDS file contains the Hellinger trajectory metadata, the GMM fit, the four-method detection sets, and the metric vectors at the primary site set and the GMM-reduced site set. |
| Full-dataset detection | A separate analysis applying the four methods to the entire dataset’s Hellinger trajectory (not sub-windowed) and matching against the full catalogue. Produces the timeline figures (figs 7–10). |

# Research questions

- **Q1.** Which of the four change-point detection methods achieves the
  highest precision, recall, F1, and lowest TLE on Hellinger
  trajectories computed from real surveillance data?
- **Q2.** Does the ranking of methods depend on surveillance window
  length (short 6–10 months, medium 12–16 months, long 18+ months)?
- **Q3.** Does GMM-based entropy reduction (per-window class-1 site
  selection) preserve change-point detection accuracy relative to using
  the full 1,273-site Spike protein?
- **Q4.** How sensitive are method rankings to ground-truth date
  perturbations of $\pm 4$ months?
- **Q5.** For the $K$-requiring method `ks.cp3o`, how does fixed
  $K \in \{1, 2, 3, 5, 10\}$ compare against the data-driven
  `dynamic_k = nrow(window) - 2`?

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

Sequences
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

Zenodo DOI 10.5281/zenodo.19040165 (post-processed)
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

Both feature matrices are integer-encoded (25-residue alphabet: 20
standard amino acids + B, Z, X, \*, gap), sorted by collection date, and
date-clipped to the first of each collection month at preprocessing time
(Sayers et al. 2022; Shu and McCauley 2017). The clipping is automatic
by virtue of
`AL_df$Date <- as.Date(format(corrected_dates, "%Y-%m-01"))` in both
preprocessing pipelines, which means every sequence’s `Date` is
unambiguously assignable to a 2-month bin in any window whose start and
end dates are themselves month boundaries.

The GISAID date range substantially exceeds NCBI’s; this is reflected in
the per-tier window-length distribution (the long tier on NCBI admits at
most lengths in $\{18, 20\}$ months whereas GISAID admits the full
$\{18, 20, 22, 24, 26, 28\}$).

# Ground truth

The ground-truth catalogue is the `sarscov2_variants` dataset bundled
with `ViralEntropR`, restricted to WHO-labelled VOCs and VOIs with a
non-missing `Date_First_Detected_US`. Variants under monitoring (VUMs),
de-escalated variants, and entries with NA US-detection dates are
excluded. The restriction is conservative: it limits attention to
variants whose surveillance significance is uncontested and whose
emergence dates are peer-curated. Cross-references via Nextstrain
(Hadfield et al. 2018) and outbreak.info (Gangavarapu et al. 2023)
inform the truth-shift sensitivity analysis to bracket the temporal
uncertainty in first-detection dates.

The number of evaluable truth events per dataset is bounded above by the
size of this catalogue intersected with each dataset’s date range: all
12 VOC/VOI entries fall within both NCBI’s and GISAID’s date ranges (the
catalogue is sparse post-Omicron; no additional named variants emerged
in `sarscov2_variants` after December 2021).

# Methodology

## Stratified bootstrap design

For each dataset, 5,001 random sub-windows are drawn, stratified by
window length:

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Stratified sampling tiers with target acceptance counts.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Tier
</th>

<th style="text-align:left;">

Window length (months)
</th>

<th style="text-align:left;">

Bins
</th>

<th style="text-align:left;">

Runs per tier per dataset
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Short
</td>

<td style="text-align:left;">

6, 8, 10
</td>

<td style="text-align:left;">

3, 4, 5
</td>

<td style="text-align:left;">

1,667
</td>

</tr>

<tr>

<td style="text-align:left;">

Medium
</td>

<td style="text-align:left;">

12, 14, 16
</td>

<td style="text-align:left;">

6, 7, 8
</td>

<td style="text-align:left;">

1,667
</td>

</tr>

<tr>

<td style="text-align:left;">

Long
</td>

<td style="text-align:left;">

18, 20, 22, 24, 26, 28 (dataset-capped)
</td>

<td style="text-align:left;">

9, 10, 11, 12, 13, 14
</td>

<td style="text-align:left;">

1,667
</td>

</tr>

</tbody>

</table>

Lengths are restricted to even integers so the bin count = length / 2 is
integer. NCBI’s 21-month range caps long-tier lengths at $\{18, 20\}$;
GISAID’s 51-month range supports the full set.

Within each tier, candidate windows are drawn by:

1.  Sampling a uniform random window length from the tier’s feasible
    set.
2.  Sampling a uniform random start offset (in months) from
    $\{0, 1, \ldots, n\}$ where $n$ is the number of feasible start
    offsets (`data_end - length` in months).
3.  Counting truth events whose `Date_First_Detected_US` falls in
    $[\texttt{start}, \texttt{end})$, dropping any in the $T_1$ bin,
    deduplicating to **unique bins** (so two variants emerging in the
    same 2-month bin count as one truth bin), and verifying at least two
    distinct truth bins remain.
4.  Counting sequences per bin; verifying every bin contains at least
    **30 sequences**.

Candidates failing any check are rejected and re-drawn until each tier
accumulates 1,667 accepted windows. Total per dataset: **5,001 accepted
windows** (1,667 $\times$ 3). Window draws are deterministic in the
master seed.

## Bin alignment

Every sequence’s `Date` is the first day of its collection month. Window
start and end dates are also month boundaries by construction, and
`(end - start)` is constrained to be even (in months). Bins are formed
left-to-right from `start` in 2-month steps; the last bin ends exactly
on `end`. No partial bins, no rounding choices, no off-by-one ambiguity.

## Hellinger trajectory

For each accepted window:

- Sequences with `Date $\in [\texttt{start}, \texttt{end})$` are
  retained.
- The retained sequences are partitioned into the bin grid via
  `bin_row_ranges()` (a `findInterval`-based partitioner that handles
  empty bins correctly).
- Per-bin amino-acid frequency tables are computed; pairwise Hellinger
  distances from $T_1$ to $T_2, T_3, \ldots, T_n$ are computed via
  `ViralEntropR::calculate_hellinger_matrix()` (`aa_levels = 25L`,
  `normalized = FALSE`).
- The returned `sites $\times$ time_steps` matrix is transposed to
  `time_steps $\times$ sites` for change-point input.

## Per-window GMM dimensionality reduction (secondary)

For the secondary analysis, **GMM-based site selection is refitted per
window** to avoid oracle leakage:

1.  Per-site Shannon entropy is computed on the window’s sequence set
    (not the full dataset).
2.  `cluster_sites_by_entropy()` fits a Gaussian mixture (default model
    search across `E` and `V`, $G \in \{1, \ldots, 15\}$) (Scrucca et
    al. 2016).
3.  `relabel_entropy_classes()` ensures class 1 denotes the
    highest-entropy group.
4.  Class-1 sites are passed to a second `calculate_hellinger_matrix()`
    call, producing the reduced-site Hellinger trajectory.

This per-window approach mirrors what a surveillance analyst would do
prospectively (entropy class assignment based only on data seen to date)
and is the honest comparison for assessing whether the entropy reduction
preserves change-point detectability.

**Edge cases handled silently:**

- If the per-window GMM returns the all-identical sentinel (class 999),
  the reduced-site analysis is skipped for that replicate and
  `reduced_skipped = TRUE` is recorded.
- If the GMM returns $G = 1$, class-1 sites are the union of all sites
  that passed the `removez = TRUE` and `removesngl = TRUE` filters (a
  sparser-than-full but non-trivial subset).

## The four change-point detection methods

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Methods compared.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Code
</th>

<th style="text-align:left;">

Function
</th>

<th style="text-align:left;">

Family
</th>

<th style="text-align:left;">

Requires K?
</th>

<th style="text-align:left;">

Key params
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

ECP-1
</td>

<td style="text-align:left;">

`ecp::ks.cp3o`
</td>

<td style="text-align:left;">

Dynamic programming (energy distance)
</td>

<td style="text-align:left;">

Yes (data-driven or fixed)
</td>

<td style="text-align:left;">

`K = nrow - 2`, `minsize = 2`, K-sweep in {1,2,3,5,10}
</td>

</tr>

<tr>

<td style="text-align:left;">

ECP-2
</td>

<td style="text-align:left;">

`ecp::e.agglo`
</td>

<td style="text-align:left;">

Agglomerative hierarchical (energy distance)
</td>

<td style="text-align:left;">

No
</td>

<td style="text-align:left;">

`alpha = 1`, `penalty = function(cps) 0`
</td>

</tr>

<tr>

<td style="text-align:left;">

HDcp-BS
</td>

<td style="text-align:left;">

`HDcpDetect::binary.segmentation`
</td>

<td style="text-align:left;">

Binary segmentation (high-dim CUSUM)
</td>

<td style="text-align:left;">

No
</td>

<td style="text-align:left;">

defaults
</td>

</tr>

<tr>

<td style="text-align:left;">

HDcp-WBS
</td>

<td style="text-align:left;">

`HDcpDetect::wild.binary.segmentation`
</td>

<td style="text-align:left;">

Wild binary segmentation (high-dim CUSUM)
</td>

<td style="text-align:left;">

No
</td>

<td style="text-align:left;">

`M = 100`
</td>

</tr>

</tbody>

</table>

`ks.cp3o` is the only $K$-requiring method. The primary analysis uses
`dynamic_k = TRUE` (sets `K = nrow(window) - 2` per call). A
**sensitivity sweep over fixed $K \in \{1, 2, 3, 5, 10\}$** runs inside
the same replicate at no additional Hellinger cost.

### Output harmonisation

The four methods return change points in three different shapes:

- `ks.cp3o`: integer vector in `$estimates` (first row of new segment)
- `e.agglo`: integer vector in `$estimates` (first row of new segment)
- `binary.segmentation`: two-column matrix with `FoundList` (last row of
  old segment) and `pvalues`
- `wild.binary.segmentation`: integer vector or the character string
  `"No Change Points Found"`

`cp_methods.R` provides a per-method extractor that returns a sorted
integer vector of interior CP indices in the harmonised convention
“first row of new segment”, with boundaries stripped to
$[2, n_\text{rows}]$. **The convention shift for the two HDcp methods is
a $+1$ index adjustment** verified empirically against the canonical
examples in each package’s help page. P-values from
`binary.segmentation` are preserved in the per-replicate RDS for
post-hoc filtering but are not used in the primary metrics.

## Truth-to-detection matching

For a window with bin-mapped truth set $T = \{t_1, \ldots, t_m\}$ and a
method’s detection set $D = \{d_1, \ldots, d_n\}$, matching proceeds via
**truth-first greedy with consumption**:

    For each t in T (in temporal order):
      Find d* in D minimising |d* - t|, subject to |d* - t| <= 1 (+/- 1 bin tolerance)
      If d* exists:
        Mark (t, d*) as a true positive
        Remove d* from D (consumption)
      Else:
        Mark t as a false negative
    Remaining unmatched detections in D are false positives.

The $\pm 1$ bin tolerance permits a one-bin (two-month) offset between
detection and truth, which is the standard convention in the
change-point evaluation literature (Aminikhanghahi and Cook 2017;
Truong, Oudre, and Vayatis 2020). The consumption rule prevents one
detection from being credited against multiple truths.

## Metrics

Per replicate, per method, per site set:

- **Precision** = TP / (TP + FP)
- **Recall** = TP / (TP + FN)
- **F1** = 2 $\cdot$ Precision $\cdot$ Recall / (Precision + Recall)
- **Temporal localisation error (TLE)** = mean absolute bin distance
  between matched (truth, detection) pairs, in bin units (one bin = 2
  months). NA when no matches.

Convention for degenerate cases (Truong, Oudre, and Vayatis 2020):

- Zero detections, $\geq 1$ truth $\rightarrow$ P = 1, R = 0, F1 = 0,
  TLE = NA
- Zero truths after $T_1$ filtering $\rightarrow$ window rejected at
  sampling time

## Sensitivity analyses

Two sensitivity analyses run alongside the primary detection at
near-zero additional cost:

**K-sensitivity sweep** (for `ks.cp3o` only): replicates re-evaluate
detection with $K \in \{1, 2, 3, 5, 10\}$ alongside the dynamic-$K$
primary. Each fixed-$K$ result is saved with its own
$(P, R, F_1, \text{TLE})$.

**Truth-shift sensitivity:** the truth set is shifted by
$\{-2, -1, 0, +1, +2\}$ bins (corresponding to $\{-4, -2, 0, +2, +4\}$
months) and matching is re-run. The detection set is unchanged across
shifts; only the bin-mapped truth indices change. Shifted-truth metrics
are saved alongside the primary.

# Computational budget and runtime

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Computational budget.
</caption>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:left;">

Value
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Datasets
</td>

<td style="text-align:left;">

2 (NCBI US, GISAID US)
</td>

</tr>

<tr>

<td style="text-align:left;">

Windows per dataset
</td>

<td style="text-align:left;">

5,001 (1,667 per tier x 3 tiers)
</td>

</tr>

<tr>

<td style="text-align:left;">

Site sets
</td>

<td style="text-align:left;">

2 (full 1,273 sites, GMM-reduced class 1)
</td>

</tr>

<tr>

<td style="text-align:left;">

Methods
</td>

<td style="text-align:left;">

4 (ECP-1, ECP-2, HDcp-BS, HDcp-WBS)
</td>

</tr>

<tr>

<td style="text-align:left;">

K-sweep variants
</td>

<td style="text-align:left;">

5 (K in {1, 2, 3, 5, 10}, ks.cp3o only)
</td>

</tr>

<tr>

<td style="text-align:left;">

Truth shifts
</td>

<td style="text-align:left;">

5 ({-4, -2, 0, +2, +4} months)
</td>

</tr>

<tr>

<td style="text-align:left;">

Method invocations
</td>

<td style="text-align:left;">

approx. 5,001 x 2 x 9 = 90,018 primary + 5x truth-shift per dataset
</td>

</tr>

</tbody>

</table>

**Measured runtime on a 24-core Linux workstation:**

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

Measured wall-clock times at N_WORKERS = 24.
</caption>

<thead>

<tr>

<th style="text-align:left;">

Stage
</th>

<th style="text-align:left;">

Wall-clock
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

`precompute_windows.R` (per dataset, manifest)
</td>

<td style="text-align:left;">

\< 30 s
</td>

</tr>

<tr>

<td style="text-align:left;">

`simulation_study.R` benchmark, NCBI US
</td>

<td style="text-align:left;">

approx. 2.5 hours (24 workers via callr::r_bg)
</td>

</tr>

<tr>

<td style="text-align:left;">

`simulation_study.R` benchmark, GISAID US
</td>

<td style="text-align:left;">

approx. 2 hours (24 workers via callr::r_bg)
</td>

</tr>

<tr>

<td style="text-align:left;">

`run_full_dataset_detection.R` (both datasets)
</td>

<td style="text-align:left;">

approx. 5 minutes
</td>

</tr>

<tr>

<td style="text-align:left;">

`plot_results.R` (10 figures)
</td>

<td style="text-align:left;">

approx. 30 seconds
</td>

</tr>

<tr>

<td style="text-align:left;">

`aggregate_tables.R` (4 CSV tables)
</td>

<td style="text-align:left;">

approx. 20 seconds
</td>

</tr>

<tr>

<td style="text-align:left;">

Total end-to-end
</td>

<td style="text-align:left;">

approx. 4.5 hours
</td>

</tr>

</tbody>

</table>

Per-subprocess peak memory is approximately 700 MB (feature matrix
load + window slice + Hellinger working set + four method states + GMM
fit). With 24 concurrent workers, peak resident set sizes to
approximately 17 GB. The per-method walltime is methodologically
irrelevant: medians are 4 ms (`e_agglo`), 1 ms (`hdcp_*`,
`ks_cp3o_dynamic`), and sub-millisecond for each fixed-$K$ `ks_cp3o`
variant; the Hellinger computation and the GMM fit dominate
per-replicate time.

**Total disk footprint** is approximately 500 MB:

- 2 $\times$ window manifests: approx. 20 MB
- 10,002 per-replicate RDS files (approx. 25–50 KB each): approx. 350 MB
- 2 summary RDS files: approx. 10 MB
- 2 full-dataset detection RDS files: approx. 5 MB
- 10 figure PNGs: approx. 2 MB
- 4 CSV tables: \< 1 MB
- Logs and metadata: \< 5 MB

# Reproducibility

- `BASE_SEED = 2025L`.
- `seed_for_replicate(dataset_id, tier_id, run_id) = BASE_SEED +   dataset_id \times 10^8 + tier_id \times 10^6 + run_id \times 10^3`.
- Window manifests are pre-computed deterministically by
  `precompute_windows.R` and saved as `outputs/windows_<dataset>.rds`
  with full metadata: tier, run_id, window_start, window_end, n_bins,
  bin_counts, truth_dates, truth_labels, truth_hellinger_idx.
- Per-replicate RDS files are written atomically (`save_rds_atomic`) and
  follow the path convention
  `outputs/replicates_<dataset>/run_<NNNN>.rds`.
- The orchestrator’s resume-scan reads on-disk replicate files at
  startup and skips already-completed work, making `simulation_study.R`
  safe to re-run after interruption.
- R and `ViralEntropR` version stamps are recorded per replicate.

# Implementation summary

The implementation follows the same patterns established in
`Sample_Size_Simulation_Study/`:

- **Concurrency:** `N_WORKERS = 1L` runs sequentially in-process with
  `tryCatch`; `N_WORKERS > 1L` dispatches via `callr::r_bg` with a
  worker pool; `N_WORKERS = NULL` auto-detects with a RAM-budgeted cap
  of 24.
- **Memory:** feature matrices loaded once per subprocess (approx.
  140–165 MB resident); per-replicate working set bounded to approx. 700
  MB (window sub-matrix + Hellinger + four method states + GMM fit).
- **I/O:** atomic RDS writes; per-replicate metadata schema documented
  in `run_one_replicate.R`.
- **Resume:** on restart, the orchestrator scans
  `outputs/replicates_<dataset>/` for existing files and re-runs only
  the missing ones.
- **Logging:** structured log via `log_msg()` to `outputs/run_log.txt`;
  per-replicate errors caught and recorded in `outputs/error_log.txt`
  without interrupting the run.

Environment overrides at orchestrator startup:

``` bash
N_WORKERS=24 RUN_TRUTH_SHIFT_SWEEP=1 Rscript simulation_study.R
```

# Per-replicate output schema

Each `outputs/replicates_<dataset>/run_<NNNN>.rds` contains a named list
with at minimum:

``` r
list(
  dataset, tier, tier_id, run_in_tier, cell_id,
  window_start, window_end, n_bins, n_seqs,
  truths_raw, truths_dropped_T1, truths_effective,
  truth_dates, truth_labels, truth_hellinger_idx,
  gmm_meta = list(G, modelName, n_class1_sites,
                  class1_sites, reduced_skipped),
  reduced_status,                     # "ok" or "skipped"
  metrics_full = list(                # 1 entry per method config
    ks_cp3o_dynamic = list(detected_cps, status, error_msg,
                           walltime_s, extra,
                           metrics_primary = list(P, R, F1, TLE,
                                                  TP, FP, FN),
                           metrics_truth_shift),
    ks_cp3o_K1      = list(...),
    ks_cp3o_K2      = list(...),
    ks_cp3o_K3      = list(...),
    ks_cp3o_K5      = list(...),
    ks_cp3o_K10     = list(...),
    e_agglo         = list(...),
    hdcp_binseg     = list(..., pvalues = ...),
    hdcp_wbs        = list(...)
  ),
  metrics_reduced = list(...),        # same shape, NULL if reduced_status = "skipped"
  walltime_s, seed_replicate,
  r_version, package_version, run_timestamp
)
```

# Limitations

1.  **Short-tier detection ceiling.** With `minsize = 2` (the `ks.cp3o`
    minimum and the implicit `HDcpDetect` minimum), windows of 3–4 bins
    cannot detect any change point; windows of 5–6 bins can detect at
    most one. `e.agglo` is unaffected. This is a structural property of
    the methods, not the data; it appears in the results as a short-tier
    recall ceiling for three of the four methods.
2.  **US-only, Spike-only.** Generalisation to other geographies and to
    other SARS-CoV-2 proteins is not tested.
3.  **Ground-truth precision.** `Date_First_Detected_US` is the first US
    detection date; whether emergence dates differ from sequential
    detectability dates is addressed only via the truth-shift
    sensitivity sweep ($\pm 4$ months), not via independent ground-truth
    sources (Hadfield et al. 2018; Gangavarapu et al. 2023).
4.  **Truth catalogue sparsity post-Omicron.** The `sarscov2_variants`
    catalogue ends at Omicron (December 2021). GISAID’s 2022–2024 window
    contains no additional named variants; sub-lineage transitions
    within Omicron BA.\* are not separately catalogued and therefore not
    evaluable as distinct truth events. Long-tier GISAID windows landing
    in the post-Omicron tail are typically rejected at sampling time for
    insufficient unique truth bins.
5.  **Method scope.** Four methods are benchmarked; many others exist
    (PELT, BCP, MOSUM, NOT, breakfast). The four chosen represent
    distinct algorithmic families with established support for
    high-dimensional input.
6.  **Per-window GMM stability.** Empirically (90,018 fits across both
    datasets) the convergence rate was 100 %, but very-short-window fits
    may be sensitive to bin-count variability. Edge cases are handled by
    the `reduced_skipped` flag in `gmm_meta`.
7.  **Bin width fixed at 2 months.** Robustness to bin-width choice is
    not addressed by this study.
8.  **K-sweep range.** Fixed $K \in \{1, 2, 3, 5, 10\}$ explores
    low-to-moderate $K$. On our data, $K = 10$ saturates `ks.cp3o`’s
    performance on NCBI’s long tier and continues to climb on GISAID’s,
    suggesting the saturation point depends on dataset temporal span.
9.  **Bootstrap independence.** Window draws are stratified random but
    not strictly independent (overlapping windows are not excluded).
    Statistical inference treats the per-method distribution of metrics
    across windows as the unit of evidence.

# ADEMP summary

<table class="table table-striped table-hover table-condensed" style="width: auto !important; margin-left: auto; margin-right: auto;">

<caption>

ADEMP framing (Morris, White, and Crowther 2019).
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

Quantify and rank the change-point detection accuracy of four
high-dimensional methods on real SARS-CoV-2 surveillance Hellinger
trajectories, with secondary attention to GMM site reduction.
</td>

</tr>

<tr>

<td style="text-align:left;">

Data-generating mechanism
</td>

<td style="text-align:left;">

Stratified random sub-windows of two real surveillance datasets (NCBI US
and GISAID US); per-window Hellinger trajectories on the 1,273-site
Spike representation and on per-window GMM-reduced class-1 sites.
</td>

</tr>

<tr>

<td style="text-align:left;">

Estimand
</td>

<td style="text-align:left;">

Per-method, per-tier, per-site-set distribution of precision, recall,
F1, and TLE against bin-mapped WHO VOC/VOI emergence dates.
</td>

</tr>

<tr>

<td style="text-align:left;">

Methods
</td>

<td style="text-align:left;">

ecp::ks.cp3o (dynamic-K + fixed K in {1, 2, 3, 5, 10}), ecp::e.agglo,
HDcpDetect::binary.segmentation, HDcpDetect::wild.binary.segmentation (M
= 100).
</td>

</tr>

<tr>

<td style="text-align:left;">

Performance measures
</td>

<td style="text-align:left;">

Truth-first greedy matching with +/- 1 bin tolerance; P, R, F1, TLE per
replicate; truth-shift sensitivity at {-4, -2, 0, +2, +4} months.
</td>

</tr>

</tbody>

</table>

# Invocation

``` bash
# Full benchmark + post-processing pipeline, 24 parallel workers,
# truth-shift sweep enabled, logs to file:
cd analysis/Changepoint_Detection_Study

N_WORKERS=24 RUN_TRUTH_SHIFT_SWEEP=1 nohup Rscript simulation_study.R \
    > outputs/run_log.txt 2>&1 &

# After the benchmark completes:
nohup Rscript run_full_dataset_detection.R \
    > outputs/full_dataset_log.txt 2>&1 &

# Then the plotting and aggregation layers:
nohup Rscript plot_results.R     > outputs/plot_log.txt   2>&1 &
nohup Rscript aggregate_tables.R > outputs/tables_log.txt 2>&1 &
```

Single-replicate debug:

``` r
setwd("analysis/Changepoint_Detection_Study")
source("setup.R")
source("helpers_windows.R")
source("helpers_hellinger.R")
source("cp_methods.R")
source("metrics.R")
source("run_one_replicate.R")

config   <- build_config()
manifest <- readRDS("outputs/windows_NCBI_US.rds")
result   <- run_one_replicate_for_row(
  manifest = manifest, row_idx = 1L, config = config
)
```

# Project structure

    analysis/Changepoint_Detection_Study/
    |-- README.Rmd                          (this file)
    |-- README.md                           (knitted GitHub output)
    |-- README.html                         (knitted HTML output)
    |-- references_changepoint.bib          (BibTeX source)
    |-- captions_and_discussion.txt         (figure captions + main-text discussion)
    |
    |-- setup.R                             (config, helpers, paths, seeds)
    |-- helpers_windows.R                   (stratified window sampler, bin alignment)
    |-- helpers_hellinger.R                 (Hellinger + per-window GMM wrappers)
    |-- cp_methods.R                        (four-method harmonisation layer + K-sweep)
    |-- metrics.R                           (truth-matching, P/R/F1/TLE, truth-shift sweep)
    |
    |-- precompute_windows.R                (deterministic window manifest builder)
    |-- run_one_replicate.R                 (per-window driver)
    |-- simulation_study.R                  (5,001-window benchmark orchestrator)
    |
    |-- run_full_dataset_detection.R        (single full-trajectory detection per dataset)
    |-- plot_results.R                      (10 publication-quality figures)
    |-- aggregate_tables.R                  (4 supplementary CSV tables)
    |
    `-- outputs/
        |-- windows_NCBI_US.rds                  (window manifest)
        |-- windows_GISAID_US.rds                (window manifest)
        |-- replicates_NCBI_US/run_NNNN.rds      (5,001 per-window result files)
        |-- replicates_GISAID_US/run_NNNN.rds    (5,001 per-window result files)
        |-- summary_NCBI_US.rds                  (long-format aggregated metrics)
        |-- summary_GISAID_US.rds                (long-format aggregated metrics)
        |-- full_dataset_detection_NCBI_US.rds   (full-dataset trajectory + 9 method results)
        |-- full_dataset_detection_GISAID_US.rds (full-dataset trajectory + 9 method results)
        |-- run_log.txt                          (orchestrator log)
        |-- full_dataset_log.txt                 (full-dataset detection log)
        |-- plot_log.txt                         (plotting log)
        |-- tables_log.txt                       (aggregation log)
        |-- error_log.txt                        (per-replicate errors, normally empty)
        |
        |-- plots/
        |   |-- fig_01_f1_by_method_tier_dataset.png
        |   |-- fig_02_operating_curve.png
        |   |-- fig_03_k_sweep_saturation.png
        |   |-- fig_04_full_vs_reduced_scatter.png
        |   |-- fig_05_failure_heatmap.png
        |   |-- fig_06_cross_dataset_rank_scatter.png
        |   |-- fig_07_timeline_NCBI_US.png            (full sites)
        |   |-- fig_08_timeline_GISAID_US.png          (full sites)
        |   |-- fig_09_timeline_NCBI_US_reduced.png    (GMM class-1 reduced sites)
        |   `-- fig_10_timeline_GISAID_US_reduced.png  (GMM class-1 reduced sites)
        |
        `-- tables/
            |-- tab_01_method_tier_dataset_summary.csv
            |-- tab_02_operating_curve.csv
            |-- tab_03_k_sweep_saturation.csv
            `-- tab_04_full_vs_reduced_delta.csv

Each per-replicate RDS (approximately 25–50 KB) holds: window metadata;
truth bin indices; per-window GMM fit summary; the four-method detection
sets and metrics on the full 1,273-site Hellinger trajectory; the same
on the GMM class-1 reduced Hellinger trajectory; the truth-shift
sensitivity sweep; and R and `ViralEntropR` version strings. The summary
RDS written by the orchestrator at the end of each dataset is a
long-format data frame with one row per `(replicate, method, site_set)`
combination, ready for `ggplot` and tidyverse aggregation.

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

<div id="ref-AminikhanghahiCook2017" class="csl-entry">

Aminikhanghahi, Samaneh, and Diane J. Cook. 2017. “A Survey of Methods
for Time Series Change Point Detection.” *Knowledge and Information
Systems* 51 (2): 339–67. <https://doi.org/10.1007/s10115-016-0987-z>.

</div>

<div id="ref-Fryzlewicz2014" class="csl-entry">

Fryzlewicz, Piotr. 2014. “Wild Binary Segmentation for Multiple
Change-Point Detection.” *The Annals of Statistics* 42 (6): 2243–81.
<https://doi.org/10.1214/14-AOS1245>.

</div>

<div id="ref-Gangavarapu2023" class="csl-entry">

Gangavarapu, Karthik, Alaa Abdel Latif, Julia L. Mullen, Manar
Alkuzweny, Emory Hufbauer, Ginger Tsueng, Emily Haag, et al. 2023.
“<span class="nocase">outbreak.info</span> Genomic Reports: Scalable and
Dynamic Surveillance of SARS-CoV-2 Variants and Mutations.” *Nature
Methods* 20 (4): 512–22. <https://doi.org/10.1038/s41592-023-01769-3>.

</div>

<div id="ref-Hadfield2018" class="csl-entry">

Hadfield, James, Colin Megill, Sidney M. Bell, John Huddleston, Barney
Potter, Charlton Callender, Pavel Sagulenko, Trevor Bedford, and Richard
A. Neher. 2018. “Nextstrain: Real-Time Tracking of Pathogen Evolution.”
*Bioinformatics* 34 (23): 4121–23.
<https://doi.org/10.1093/bioinformatics/bty407>.

</div>

<div id="ref-Matteson2014" class="csl-entry">

James, Nicholas A., and David S. Matteson. 2014. “Ecp: An R Package for
Nonparametric Multiple Change Point Analysis of Multivariate Data.”
*Journal of Statistical Software* 62 (7): 1–25.
<https://doi.org/10.18637/jss.v062.i07>.

</div>

<div id="ref-MattesonJames2014" class="csl-entry">

Matteson, David S., and Nicholas A. James. 2014. “A Nonparametric
Approach for Multiple Change Point Analysis of Multivariate Data.”
*Journal of the American Statistical Association* 109 (505): 334–45.
<https://doi.org/10.1080/01621459.2013.849605>.

</div>

<div id="ref-Morris2019" class="csl-entry">

Morris, Tim P., Ian R. White, and Michael J. Crowther. 2019. “Using
Simulation Studies to Evaluate Statistical Methods.” *Statistics in
Medicine* 38 (11): 2074–2102. <https://doi.org/10.1002/sim.8086>.

</div>

<div id="ref-Olshen2004" class="csl-entry">

Olshen, Adam B., E. S. Venkatraman, Robert Lucito, and Michael Wigler.
2004. “Circular Binary Segmentation for the Analysis of Array-Based DNA
Copy Number Data.” *Biostatistics* 5 (4): 557–72.
<https://doi.org/10.1093/biostatistics/kxh008>.

</div>

<div id="ref-Sayers2022" class="csl-entry">

Sayers, Eric W., Evan E. Bolton, J. Rodney Brister, Kathi Canese,
Jessica Chan, Donald C. Comeau, Ryan Connor, et al. 2022. “Database
Resources of the National Center for Biotechnology Information.”
*Nucleic Acids Research* 50 (D1): D20–26.
<https://doi.org/10.1093/nar/gkab1112>.

</div>

<div id="ref-Scrucca2016" class="csl-entry">

Scrucca, Luca, Michael Fop, Thomas Brendan Murphy, and Adrian E.
Raftery. 2016. “Mclust 5: Clustering, Classification and Density
Estimation Using Gaussian Finite Mixture Models.” *The R Journal* 8 (1):
289–317. <https://doi.org/10.32614/RJ-2016-021>.

</div>

<div id="ref-ShuMcCauley2017" class="csl-entry">

Shu, Yuelong, and John McCauley. 2017. “GISAID: Global Initiative on
Sharing All Influenza Data – from Vision to Reality.” *Eurosurveillance*
22 (13): 30494. <https://doi.org/10.2807/1560-7917.ES.2017.22.13.30494>.

</div>

<div id="ref-TruongOudreVayatis2020" class="csl-entry">

Truong, Charles, Laurent Oudre, and Nicolas Vayatis. 2020. “Selective
Review of Offline Change Point Detection Methods.” *Signal Processing*
167: 107299. <https://doi.org/10.1016/j.sigpro.2019.107299>.

</div>

<div id="ref-Tyuryaev2026" class="csl-entry">

Tyuryaev, Vadim, Jane Heffernan, and Hanna Jankowski. 2026.
*ViralEntropR: A Computational Pipeline for Entropy-Informed Detection
of Emerging Viral Variants*.
<https://github.com/vadimtyuryaev/ViralEntropR>.

</div>

</div>
