
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ViralEntropR

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/ViralEntropR)](https://CRAN.R-project.org/package=ViralEntropR)
[![R-CMD-check](https://github.com/vadimtyuryaev/ViralEntropR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vadimtyuryaev/ViralEntropR/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**ViralEntropR** is an R package for the efficient preprocessing of
large FASTA archives and computational surveillance of emerging viral
variants. It provides a fully vectorized preprocessing layer, including
header parsing, country and date extraction, ambiguity filtering, and
amino-acid encoding, together with the core computational components for
downstream analysis: per-site Shannon entropy (Shannon 1948),
time-window partitioning, Gaussian-mixture-based entropy clustering
(Everitt et al. 2011; Scrucca et al. 2016), pairwise Hellinger distances
(van der Vaart 1998) between empirical amino-acid distributions, and
non-parametric change-point detection (Fryzlewicz 2014; Matteson and
James 2014). A controlled variant-simulation engine is included for
benchmarking detection performance against known ground truth.

The package was developed alongside a proposed analysis pipeline for
emerging-variant detection, demonstrated in the bundled vignettes:
GMM-driven selection of high-entropy sites, partitioning around medoids
on the Gower distance (Gower 1971; Kaufman and Rousseeuw 1990; Rousseeuw
1987) over the selected sites, t-SNE visualization of cluster geometry
(van der Maaten and Hinton 2008), non-parametric change-point detection
validation, and hierarchical agglomerative clustering (HAC) (Everitt et
al. 2011; Murtagh and Legendre 2014; Sangalli et al. 2010; Tucker et al.
2013) of per-site entropy curves and time-indexed Hellinger-distance
trajectories, treated within a functional-data-analysis framework
(Ramsay and Silverman 2005) as discrete realizations of continuous
functions of time. The pipeline was evaluated on filtered,
post-processed SARS-CoV-2 Spike-protein sequences: 137,132 raw NCBI
records (Sayers et al. 2022) reduced to **109,536** US-filtered
sequences (archived on Zenodo at DOI
[10.5281/zenodo.19040165](https://doi.org/10.5281/zenodo.19040165)), and
16.7 million raw GISAID records (Shu and McCauley 2017) reduced to
**129,371** unique US-filtered sequences. The methods are general enough
to apply to any aligned amino-acid or nucleotide time series.

## Pipeline

                           ┌─────────────────────────┐
                           │       FASTA file        │
                           └─────────────────────────┘
                                        │
                                        ▼
                    ┌───────────────────────────────────────┐
                    │     All sequences of equal width?     │
                    └───────────────────────────────────────┘
                       │                                 │
                      yes                               no
                       │                                 │
                       │                                 ▼
                       │                 ┌───────────────────────────────┐
                       │                 │       Multiple-sequence       │
                       │                 │        alignment, e.g.        │
                       │                 │         msa::msa() /          │
                       │                 │     DECIPHER::AlignSeqs()     │
                       │                 └───────────────────────────────┘
                       │                                 │
                       └────────────────┬────────────────┘
                                        ▼
                        ┌───────────────────────────────┐
                        │     Aligned sequence set      │
                        └───────────────────────────────┘
                                        │
                                        ▼
              ┌───────────────────────────────────────────────────┐
              │                Preprocessing layer                │
              │  ───────────────────────────────────────────────  │
              │           extract_fasta_dates()                   │
              │           extract_fasta_countries()               │
              │           fasta_to_char_matrix()                  │
              │           filter_ambiguous_sequences()            │
              │           encode_aa_sequence()                    │
              └───────────────────────────────────────────────────┘
                                        │
                                        ▼
                    ┌───────────────────────────────────────┐
                    │            Feature matrix:            │
                    │        m sequences  ×  n sites  +     │
                    │               Date  +  Country        │
                    └───────────────────────────────────────┘
                                        │
                                        ▼
              ┌───────────────────────────────────────────────────┐
              │             partition_time_windows()              │
              │                                                   │
              │     cumulative / sliding / disjoint               │
              │     per-window entropies + GMM,                   │
              │     via cluster_sites_by_entropy() internally     │
              └───────────────────────────────────────────────────┘
                                        │
                       ┌────────────────┴────────────────┐
                       ▼                                 ▼
       ┌───────────────────────────────┐ ┌───────────────────────────────┐
       │      Entropy / GMM-class      │ │ calculate_hellinger_matrix()  │
       │         trajectories          │ │                               │
       │                               │ │                               │
       │ relabel_entropy_classes()     │ │     pairwise per-window       │
       │ plot_entropy_trajectories()   │ │     residue-distribution      │
       │ plot_site_class_trajectory()  │ │     distances from the        │
       │ tabulate_site_evolution()     │ │     reference window          │
       └───────────────────────────────┘ └───────────────────────────────┘
                       │                                 │
                       └────────────────┬────────────────┘
                                        ▼
              ┌───────────────────────────────────────────────────┐
              │       Non-parametric change-point detection       │
              │      applied to entropy AND Hellinger series      │
              │  ───────────────────────────────────────────────  │
              │            detect_changepoints_ecp()              │
              │            detect_changepoints_hdcp()             │
              └───────────────────────────────────────────────────┘
                                        │
                       ┌────────────────┴────────────────┐
                       ▼                                 ▼
       ┌───────────────────────────────┐ ┌───────────────────────────────┐
       │   Gower + Silhouette + PAM    │ │   HAC of entropy curves and   │
       │     + t-SNE visualization     │ │    Hellinger trajectories     │
       │                               │ │                               │
       │       cluster::daisy(),       │ │    fdacluster::fdahclust()    │
       │       cluster::pam(),         │ │                               │
       │       Rtsne::Rtsne()          │ │                               │
       │                               │ │                               │
       │  orchestrated in vignettes;   │ │  orchestrated in vignettes;   │
       │  not exported by the package  │ │  not exported by the package  │
       └───────────────────────────────┘ └───────────────────────────────┘
                       │                                 │
                       └────────────────┬────────────────┘
                                        ▼
    ┌───────────────────────────────────────────────────────────────────────┐
    │                           Pipeline outputs                            │
    │  ───────────────────────────────────────────────────────────────────  │
    │    • Post-processed feature matrix                                    │
    │    • Partitioned time windows with GMM-ranked sites per partition     │
    │    • Hellinger distance matrices                                      │
    │    • Clustered sites (candidate variants)                             │
    │    • Change points (variant emergence times)                          │
    │    • Frequency tables + entropy class changes                         │
    │      (sites of high variability / under selection)                    │
    └───────────────────────────────────────────────────────────────────────┘

## Installation

You can install the released version from CRAN:

``` r
install.packages("ViralEntropR")
```

The development version is on GitHub:

``` r
# install.packages("remotes")
remotes::install_github("vadimtyuryaev/ViralEntropR")
```

To build the vignettes locally, install with:

``` r
remotes::install_github("vadimtyuryaev/ViralEntropR", build_vignettes = TRUE)
```

## Quick start

The package ships with a 100-sequence sample of NCBI SARS-CoV-2 Spike
sequences in `inst/extdata/` for runnable examples:

``` r
library(ViralEntropR)

# 1. Load the bundled sample (100 NCBI Spike sequences).
sample_path <- system.file("extdata", "sarscov2_sample.fasta.gz",
                           package = "ViralEntropR")
fasta <- Biostrings::readAAStringSet(sample_path)

# 2. Extract dates and countries from FASTA headers.
#    NCBI Virus format: option = 4 (date at end), position = 2 (between pipes).
dates_result     <- extract_fasta_dates(fasta, option = 4)
countries_result <- extract_fasta_countries(fasta, position = 2)
# diagnostic: were any headers unparseable?
dates_result$message       
countries_result$message

# 3. Drop sequences whose date or country could not be parsed.
#    `missing_id` is NA when nothing failed and an integer vector of indices
#    otherwise — strip NAs from the union before use.
ids      <- c(dates_result$missing_id, countries_result$missing_id)
drop_ids <- unique(ids[!is.na(ids)])
keep     <- setdiff(seq_len(length(fasta)), drop_ids)

fasta           <- fasta[keep]
corrected_dates <- dates_result$corrected_dates[keep]
countries       <- countries_result$countries[keep]

# 4. Convert to a character matrix, drop sequences with ambiguous residues
#    (B, J, X, Z), and align metadata to the surviving rows.
char_mat        <- fasta_to_char_matrix(fasta)
filtered        <- filter_ambiguous_sequences(char_mat, option = 2)
keep_idx        <- setdiff(seq_len(nrow(char_mat)), filtered$DeletedSeqId)
corrected_dates <- corrected_dates[keep_idx]
countries       <- countries[keep_idx]

# 5. Integer-encode under the 25-symbol ViralEntropR alphabet.
int_mat <- encode_aa_sequence(filtered$FilteredMatrix)

# 6. Assemble the analysis-ready data frame: sites 1..1273 + Date + Country.
AL_df           <- as.data.frame(int_mat)
colnames(AL_df) <- as.character(seq_len(ncol(int_mat)))
AL_df$Date      <- as.Date(format(corrected_dates, "%Y-%m-01"))
AL_df$Country   <- countries
AL_df           <- AL_df[order(AL_df$Date), ]

# 7. Per-site Shannon entropy + GMM-based site classification.
#    Sites of zero entropy and singletons are excluded. 
ent         <- apply(int_mat, 2, calculate_entropy)
cls         <- cluster_sites_by_entropy(ent, nr = nrow(int_mat))
cls_labeled <- relabel_entropy_classes(cls$DataFrame)
head(cls_labeled, 10)
```

For the full pipeline — including time-window partitioning, Hellinger
distances, change-point detection, and visualization — see the bundled
vignettes.

## Vignettes

Three pre-rendered vignettes walk through the full workflow:

``` r
browseVignettes("ViralEntropR")
```

| Vignette | Topic |
|----|----|
| `preprocessing_pipeline` | Complete preprocessing of raw NCBI SARS-CoV-2 Spike-protein FASTA: two-pass date extraction (`yyyy-mm-dd` then `yyyy-mm`), country extraction from pipe-delimited headers, ambiguity filtering (B/J/X/Z), 25-symbol integer encoding, and assembly of an analysis-ready data frame keyed on `Date` and `Country`. Applied to the full 137,132-sequence NCBI Spike archive on Zenodo. |
| `detecting_variants_simulation` | End-to-end variant detection demonstration on a controlled synthetic benchmark: four variants emerging over 24 months under pairwise competition with stochastic growth; per-partition entropies and Hellinger distances; three complementary change-point methods — `ks.cp3o` (dynamic programming, exact globally optimal), `detect_changepoints_ecp()` (expanding-window for online surveillance), and `e.agglo` (agglomerative hierarchical, K-free) — compared against the known emergence schedule. |
| `clustering_accuracy` | Empirical evaluation of Entropy-GMM-Gower-Silhouette-PAM clustering against a labelled ground truth: wild-type period (May–June 2020) versus Delta-dominant period (July–August 2021); GMM-driven site selection, PAM on the Gower distance with silhouette-based k (optimal number of clusters) selection over k = 2 … 40, 2D t-SNE embedding, precision / recall / F1 at selected k levels, and medoid analysis cross-referenced against the curated mutation catalog of SARS-CoV-2 mutations in `sarscov2_variants`. |

## Analysis

### Completed

Three standalone scripts in `analysis/` (excluded from the installed
package via `.Rbuildignore`; available by cloning the GitHub
repository):

| Script (in `analysis/`) | Topic |
|----|----|
| `GISAID_data_preprocessing/` | Reproducible six-stage preprocessing of the GISAID Spike-protein archive: load, date / country extraction, US country and date-window filter, ambiguity filtering, ClustalOmega alignment of short sequences (width \< 1273 aa) merged with full-length sequences, encoding, and final data-frame assembly. |
| `GMM_transformation_study/` | Empirical robustness study of GMM-based site selection under six monotone transformations of per-site entropy against the untransformed baseline; selected sites cross-referenced against VOC SNP catalogues (Alpha / Beta / Gamma / Delta for NCBI; with Omicron added for GISAID data). Conclusion: entropy-based GMM site selection is invariant to the choice of transformation. |
| `Entropy_Partitioning_Study/` | Site-specific Shannon entropy dynamics at variant-defining mutation sites under three temporal partitioning strategies — cumulative (1-month expanding from a fixed origin), sliding 2-month (one-month overlap), and disjoint 2-month — for six variants (Alpha (B.1.1.7), Beta (B.1.351), Gamma (P.1), Delta (B.1.617.2), Omicron (B.1.1.529), D614G). Four output types are produced and each is interpretable: entropy trajectories identify sites whose variability rises or falls together over time, providing the basis for the functional-data-analysis treatment; class-trajectory plots track each site’s GMM class assignment across partitions, with sites under positive selection expected to enter the highest-entropy class earlier and to escalate through classes faster than neutral sites; amino-acid frequency tables record which substitutions appear at each defining SNP position in each partition; and GMM class-assignment tables summarise per-site class labels per partition, with `Nseq` and `Nclust` summary rows giving the per-partition sample size and the number of fitted GMM components. |

### Forthcoming analyses

Three additional analyses are in progress.

| Script (in `analysis/`) | Topic |
|----|----|
| `Sample_Size_Simulation_Study/` | Simulation study quantifying the number of sequences required to detect an emerging variant and to recover its associated change point under the proposed pipeline. Synthetic populations are generated by `simulate_variant_evolution()` across scenarios that vary the number of variants, their emergence schedule, per-variant growth multipliers, deleterious-mutation rate, and random-substitution event probability. For each scenario the proposed pipeline is applied end-to-end, and two summary statistics are recorded per iteration: the sample size at which **all** distinct variant-defining mutation-site sets are simultaneously detected, and the per-scenario change-point accuracy. |
| `FDA_Analysis/` | Hierarchical agglomerative clustering of per-site Shannon entropy curves and pairwise Hellinger distance trajectories within a functional-data-analysis framework, for the five WHO-designated SARS-CoV-2 variants of concern (Alpha (B.1.1.7), Beta (B.1.351), Gamma (P.1), Delta (B.1.617.2), Omicron (B.1.1.529)). Trajectories are treated as discrete realizations of continuous functions of time and clustered using HAC framework. Outputs include animated `.gif` visualizations of the clustered entropy and Hellinger trajectories themselves, showing how cluster membership and curve shape evolve before, at, and after official detection dates for the variants in question. |
| `CP_Detection_Study/` | Systematic evaluation of the two non-parametric change-point detection methods exported by the package — `detect_changepoints_ecp()` and `detect_changepoints_hdcp()` — applied to the post-processed NCBI and GISAID Spike-protein time series. Detected change points are compared against the documented first-detection dates of WHO-designated variants of concern, taken from `sarscov2_variants$Date_First_Detected_US`. Detection accuracy is reported across the three partitioning strategies (cumulative, sliding, disjoint). |

## Data

- **`sarscov2_variants`** — bundled metadata for twelve WHO-designated
  SARS-CoV-2 variants of concern and variants of interest, plus 21
  peer-reviewed literature references with DOIs. Loaded automatically
  with the package; see `?sarscov2_variants`.
- **`inst/extdata/sarscov2_sample.fasta.gz`** — 100 randomly sampled
  NCBI Spike-protein sequences for runnable examples (see the *Quick
  start* above).
- **`data-raw/NCBI_US_unaligned_feature_matrix_1273aa.rds`** — the
  post-processed NCBI feature matrix (US-filtered Spike protein, 109,536
  sequences × 1,273 sites with `Date` and `Country` columns) produced by
  the `preprocessing_pipeline` vignette. Tracked in the repository for
  direct downstream use.
- **Source NCBI Spike-protein archive** — the input FASTA underlying the
  preprocessing pipeline (137,132 sequences, ~181.5 MB uncompressed) is
  archived on Zenodo: <https://doi.org/10.5281/zenodo.19040165>.
- **GISAID data are not redistributed.** The raw GISAID Spike-protein
  archive (`spikeprot0410.fasta`, ~22.6 GB) and the derived
  post-processed feature matrix
  (`GISAID_US_aligned_feature_matrix_1273aa.rds`, 129,371 US-filtered
  sequences) are excluded from the repository and from Zenodo in
  accordance with the [GISAID Database Access
  Agreement](https://gisaid.org/terms-of-use/). Users may reproduce
  these objects independently by registering at
  [GISAID](https://gisaid.org/), downloading the equivalent release, and
  re-running `analysis/GISAID_data_preprocessing/`.

## Installing on R server and computational notes

The vignettes and analysis scripts depend on `kableExtra` for HTML table
rendering. On a fresh Linux R server `install.packages("kableExtra")`
may fail because its dependencies `systemfonts` and `svglite` need
system-level `cairo`, `freetype`, and `fontconfig` development headers.
If installation fails, install those headers first
(`apt install libcairo2-dev libfontconfig1-dev libfreetype-dev libxt-dev`
on Debian / Ubuntu, or the equivalent `-devel` packages on RHEL /
CentOS) and retry; if the server administrator has already installed
`kableExtra` system-wide, no action is needed. `Biostrings` is similarly
required and must be installed from Bioconductor
(`BiocManager::install("Biostrings")`).

Three scripts currently in `analysis/`, as well as scripts in progress
that will appear there, have substantial memory or runtime footprints
and should be run on a server rather than a laptop. Estimated peak RAM
varies between approximately 15 GB and 112 GB depending on the script
and the dataset to which it is applied.

The three pre-rendered vignettes and the `sarscov2_sample.fasta.gz`
quick-start are designed to run on a laptop without any of the above
caveats.

## Citation

If you use **ViralEntropR** in your research, please cite both the
software and the underlying methods:

``` r
citation("ViralEntropR")
```

BibTeX entry:

    @Manual{,
      title = {ViralEntropR: A Computational Pipeline for Entropy-Informed Detection of Emerging Viral Variants},
      author = {Vadim Tyuryaev and Jane Heffernan and Hanna Jankowski},
      note = {R package version 0.6.0},
      url = {https://github.com/vadimtyuryaev/ViralEntropR},
      }

## Acknowledgments

We gratefully acknowledge the authors from the originating laboratories
responsible for obtaining the specimens, and the submitting laboratories
that generated and shared the genetic sequence data via GISAID and NCBI
Virus, on which this research is based.

## License

MIT © Vadim Tyuryaev, Jane Heffernan, Hanna Jankowski. See
[`LICENSE.md`](LICENSE.md) for the full text.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Everitt2011" class="csl-entry">

Everitt, B. S., Landau, S., Leese, M., and Stahl, D. (2011), *Cluster
analysis*, Wiley. <https://doi.org/10.1002/9780470977811>.

</div>

<div id="ref-FerratyVieu2006" class="csl-entry">

Ferraty, F., and Vieu, P. (2006), *Nonparametric functional data
analysis: Theory and practice*, Springer.
<https://doi.org/10.1007/0-387-36620-2>.

</div>

<div id="ref-Fryzlewicz2014" class="csl-entry">

Fryzlewicz, P. (2014), “Wild binary segmentation for multiple
change-point detection,” *The Annals of Statistics*, 42, 2243–2281.
<https://doi.org/10.1214/14-AOS1245>.

</div>

<div id="ref-Gower1971" class="csl-entry">

Gower, J. C. (1971), “A general coefficient of similarity and some of
its properties,” *Biometrics*, 27, 857–871.
<https://doi.org/10.2307/2528823>.

</div>

<div id="ref-Kaufman1990" class="csl-entry">

Kaufman, L., and Rousseeuw, P. J. (1990), *Finding groups in data: An
introduction to cluster analysis*, Wiley.
<https://doi.org/10.1002/9780470316801>.

</div>

<div id="ref-Li2015" class="csl-entry">

Li, X., Jankowski, H., Wang, X., and Heffernan, J. M. (2015), “A method
for clustering hemagglutinin influenza protein sequences,” in *BIOMAT
2014*, World Scientific, pp. 272–285.
<https://doi.org/10.1142/9789814667944_0018>.

</div>

<div id="ref-MattesonJames2014" class="csl-entry">

Matteson, D. S., and James, N. A. (2014), “A nonparametric approach for
multiple change point analysis of multivariate data,” *Journal of the
American Statistical Association*, 109, 334–345.
<https://doi.org/10.1080/01621459.2013.849605>.

</div>

<div id="ref-MurtaghLegendre2014" class="csl-entry">

Murtagh, F., and Legendre, P. (2014), “Ward’s hierarchical agglomerative
clustering method: Which algorithms implement ward’s criterion?”
*Journal of Classification*, 31, 274–295.
<https://doi.org/10.1007/s00357-014-9161-z>.

</div>

<div id="ref-RamsaySilverman2005" class="csl-entry">

Ramsay, J. O., and Silverman, B. W. (2005), *Functional data analysis*,
Springer. <https://doi.org/10.1007/b98888>.

</div>

<div id="ref-Rousseeuw1987" class="csl-entry">

Rousseeuw, P. J. (1987), “Silhouettes: A graphical aid to the
interpretation and validation of cluster analysis,” *Journal of
Computational and Applied Mathematics*, 20, 53–65.
<https://doi.org/10.1016/0377-0427(87)90125-7>.

</div>

<div id="ref-Sangalli2010" class="csl-entry">

Sangalli, L. M., Secchi, P., Vantini, S., and Vitelli, V. (2010),
“K-mean alignment for curve clustering,” *Computational Statistics &
Data Analysis*, 54, 1219–1233.
<https://doi.org/10.1016/j.csda.2009.12.008>.

</div>

<div id="ref-Sayers2022" class="csl-entry">

Sayers, E. W., Bolton, E. E., Brister, J. R., Canese, K., Chan, J.,
Comeau, D. C., Connor, R., Funk, K., Kelly, C., Kim, S., Madej, T.,
Marchler-Bauer, A., Lanczycki, C., Lathrop, S., Lu, Z., Thibaud-Nissen,
F., Murphy, T., Phan, L., Skripchenko, Y., Tse, T., Wang, J., Williams,
R., Trawick, B. W., Pruitt, K. D., and Sherry, S. T. (2022), “Database
resources of the national center for biotechnology information,”
*Nucleic Acids Research*, 50, D20–D26.
<https://doi.org/10.1093/nar/gkab1112>.

</div>

<div id="ref-Scrucca2016" class="csl-entry">

Scrucca, L., Fop, M., Murphy, T. B., and Raftery, A. E. (2016),
“<span class="nocase">mclust</span> 5: Clustering, classification and
density estimation using Gaussian finite mixture models,” *The R
Journal*, 8, 289–317. <https://doi.org/10.32614/RJ-2016-021>.

</div>

<div id="ref-Shannon1948" class="csl-entry">

Shannon, C. E. (1948), “A mathematical theory of communication,” *Bell
System Technical Journal*, 27, 379–423.
<https://doi.org/10.1002/j.1538-7305.1948.tb01338.x>.

</div>

<div id="ref-ShuMcCauley2017" class="csl-entry">

Shu, Y., and McCauley, J. (2017), “GISAID: Global initiative on sharing
all influenza data – from vision to reality,” *Eurosurveillance*, 22,
30494. <https://doi.org/10.2807/1560-7917.ES.2017.22.13.30494>.

</div>

<div id="ref-Sievers2011" class="csl-entry">

Sievers, F., Wilm, A., Dineen, D., Gibson, T. J., Karplus, K., Li, W.,
Lopez, R., McWilliam, H., Remmert, M., Söding, J., Thompson, J. D., and
Higgins, D. G. (2011), “Fast, scalable generation of high-quality
protein multiple sequence alignments using Clustal Omega,” *Molecular
Systems Biology*, 7, 539. <https://doi.org/10.1038/msb.2011.75>.

</div>

<div id="ref-Tucker2013" class="csl-entry">

Tucker, J. D., Wu, W., and Srivastava, A. (2013), “Generative models for
functional data using phase and amplitude separation,” *Computational
Statistics & Data Analysis*, 61, 50–66.
<https://doi.org/10.1016/j.csda.2012.12.001>.

</div>

<div id="ref-Tyuryaev2026" class="csl-entry">

Tyuryaev, V., Jankowski, H., and Heffernan, J. M. (2026), “SARS-CoV-2
surface glycoprotein sequences, NCBI data hub, October 2021
(ViralEntropR archive),” Zenodo.
<https://doi.org/10.5281/zenodo.19040165>.

</div>

<div id="ref-vanderMaaten2008" class="csl-entry">

van der Maaten, L., and Hinton, G. (2008), “[Visualizing data using
t-SNE](http://jmlr.org/papers/v9/vandermaaten08a.html),” *Journal of
Machine Learning Research*, 9, 2579–2605.

</div>

<div id="ref-vanDerVaart1998" class="csl-entry">

van der Vaart, A. W. (1998), *Asymptotic statistics*, Cambridge series
in statistical and probabilistic mathematics, Cambridge University
Press. <https://doi.org/10.1017/CBO9780511802256>.

</div>

</div>
