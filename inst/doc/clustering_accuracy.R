knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  warning   = FALSE,
  message   = FALSE,
  fig.width = 9,
  fig.height = 6
)

# Set working directory to package root for local knitting
knitr::opts_knit$set(root.dir = "C:/YORK_PhD/RESEARCH/PAPERS/GitHub/ViralEntropR")

library(ViralEntropR)
library(dplyr)
library(ggplot2)
library(scales)       # pretty_breaks for silhouette plot
library(cluster)      # pam(), daisy()
library(Rtsne)        # t-SNE visualisation
library(future.apply) # parallel PAM silhouette loop
library(kableExtra)   # compact HTML tables

# ── User-modifiable parameters ───────────────────────────────────────────────

# Data
COUNTRY         <- "USA"          # filter to this country (NULL = all countries)
UNIQUE_SEQS     <- TRUE           # TRUE = deduplicate sequences before PAM
                                  # FALSE = keep all sequences (including duplicates)
N_SITES         <- 1273L          # number of Spike protein sites

# Time windows
PERIOD1_START   <- "2020-05-01"   # Period 1: wild-type phase
PERIOD1_END     <- "2020-07-01"
PERIOD2_START   <- "2021-07-01"   # Period 2: Delta-dominant phase
PERIOD2_END     <- "2021-09-01"

# PAM silhouette sweep
K_MAX           <- 40L                # maximum number of PAM clusters to evaluate
K_HIGHLIGHT     <- c(2L, 3L, 4L, 9L)  # PAM k values to plot in detail

# t-SNE parameters
TSNE_PERPLEXITY <- 10L
TSNE_SEED       <- 42L

# Parallelism (future.apply)
N_WORKERS       <- max(1L, parallel::detectCores() - 1L)

# -- Replace this path with your actual preprocessed data object -------------
# AL_df <- readRDS("path/to/AL_Cat_preprocessed.rds")
# For demonstration we load from the bundled sample -- for full results
# use the complete preprocessed NCBI dataset.
# ---------------------------------------------------------------------------

fasta_path <- file.path("data-raw", "sequences.fasta")
fasta      <- Biostrings::readAAStringSet(fasta_path)

# --- Dates ------------------------------------------------------------------
# Pass 1: identify sequences with year-only dates
dates_ymd <- extract_fasta_dates(
  fasta,
  custom_pattern = "(?<=\\|)[0-9]{4}-(0?[1-9]|1[0-2])-(0?[1-9]|[12][0-9]|3[01]|00)"
)

# Pass 2: extract yyyy-mm, remove year-only
dates_result <- extract_fasta_dates(
  fasta,
  custom_pattern = "(?<=\\|)[0-9]{4}-(0[1-9]|1[0-2])",
  date_format    = "%Y-%m"
)

if (!is.na(dates_result$missing_id[1]) && length(dates_result$missing_id) > 0) {
  fasta        <- fasta[-dates_result$missing_id]
  dates_result <- extract_fasta_dates(
    fasta,
    custom_pattern = "(?<=\\|)[0-9]{4}-(0[1-9]|1[0-2])",
    date_format    = "%Y-%m"
  )
}

# --- Countries --------------------------------------------------------------
countries_result <- extract_fasta_countries(fasta, position = 2)

if (!is.na(countries_result$missing_id[1]) && length(countries_result$missing_id) > 0) {
  fasta            <- fasta[-countries_result$missing_id]
  countries_result <- extract_fasta_countries(fasta, position = 2)
  dates_result     <- extract_fasta_dates(
    fasta,
    custom_pattern = "(?<=\\|)[0-9]{4}-(0[1-9]|1[0-2])",
    date_format    = "%Y-%m"
  )
}

# --- Convert, filter, encode ------------------------------------------------
char_mat       <- fasta_to_char_matrix(fasta)
filtered       <- filter_ambiguous_sequences(char_mat, option = 2)
clean_char_mat <- filtered$FilteredMatrix
deleted_rows   <- filtered$DeletedSeqId

corrected_dates <- dates_result$corrected_dates[-deleted_rows]
countries       <- countries_result$countries[-deleted_rows]

int_mat <- encode_aa_sequence(clean_char_mat)

# --- Assemble data frame ----------------------------------------------------
n_sites         <- ncol(int_mat)
AL_df           <- as.data.frame(int_mat)
colnames(AL_df) <- as.character(seq_len(n_sites))
AL_df[]         <- lapply(AL_df, as.integer)
AL_df$Date      <- as.Date(format(corrected_dates, "%Y-%m-01"))
AL_df$Country   <- countries
AL_df           <- AL_df[order(AL_df$Date), ]
rownames(AL_df) <- NULL

sprintf("Final data frame: %d sequences x %d sites", nrow(AL_df), n_sites)


# Filter by country if specified
if (!is.null(COUNTRY)) {
  AL_df <- AL_df[!is.na(AL_df$Country) & AL_df$Country == COUNTRY, ]
  sprintf("Sequences after country filter (%s): %d", COUNTRY, nrow(AL_df))
}

part1 <- partition_time_windows(
  data          = AL_df,
  n_sites       = N_SITES,
  window_length = 2L,
  window_type   = 2L,          # sliding (2-month window)
  start_date    = PERIOD1_START,
  end_date      = PERIOD1_END,
  verbose       = FALSE
)

part2 <- partition_time_windows(
  data          = AL_df,
  n_sites       = N_SITES,
  window_length = 2L,
  window_type   = 2L,
  start_date    = PERIOD2_START,
  end_date      = PERIOD2_END,
  verbose       = FALSE
)

df_p1 <- part1$Partitions[[1]]
df_p2 <- part2$Partitions[[1]]

sprintf("Period 1 (%s to %s): %d sequences",
        PERIOD1_START, PERIOD1_END, nrow(df_p1))
sprintf("Period 2 (%s to %s): %d sequences",
        PERIOD2_START, PERIOD2_END, nrow(df_p2))

# Helper: extract highest-entropy sites from a cluster DataFrame
top_sites <- function(clust_df) {
  sort(clust_df[clust_df$class == max(clust_df$class), ]$sites)
}

# Period 1 top sites
sites_p1 <- top_sites(part1$Clusters[[1]]$DataFrame)

cat("Period 1 highest-entropy sites:\n")
print(sites_p1)

# Period 2 top sites
sites_p2 <- top_sites(part2$Clusters[[1]]$DataFrame)

cat("\nPeriod 2 highest-entropy sites:\n")
print(sites_p2)

# Combine periods -- rows are sequences, columns are sites
df_combined           <- rbind(df_p1, df_p2)
rownames(df_combined) <- seq_len(nrow(df_combined))

# Period 1 size -- used throughout for labelling
n_p1 <- nrow(df_p1)
n_p2 <- nrow(df_p2)
n_total <- nrow(df_combined)

cat(sprintf("Combined: %d sequences (%d Period 1 + %d Period 2)\n",
            n_total, n_p1, n_p2))

# Entropy on site columns only
entrp_combined <- apply(df_combined[, seq_len(N_SITES), drop = FALSE],
                        2, calculate_entropy)

clust_combined <- cluster_sites_by_entropy(entrp_combined, nr = n_total)
clust_combined_df <- clust_combined$DataFrame

selected_sites <- sort(
  clust_combined_df[clust_combined_df$class == max(clust_combined_df$class), ]$sites
)

cat(sprintf("\nGMM selected %d sites from combined data:\n",
            length(selected_sites)))
print(selected_sites)

# ── Helper: pull SNP catalogue for a named variant ──────────────────────────
get_variant_snps <- function(label) {
  idx   <- which(unlist(sarscov2_variants$WHO_Label) == label)
  sites <- sarscov2_variants$Defining_SNP_Sites[[idx]]
  snps  <- sarscov2_variants$Defining_SNPs[[idx]]
  list(label = label, sites = sites, snps = snps)
}

# Retrieve catalogues for the four VOCs of interest
var_delta <- get_variant_snps("Delta")
var_alpha <- get_variant_snps("Alpha")
var_beta  <- get_variant_snps("Beta")
var_gamma <- get_variant_snps("Gamma")

# Convenience aliases (used throughout downstream chunks)
delta_sites <- var_delta$sites;  delta_snps <- var_delta$snps
alpha_sites <- var_alpha$sites;  alpha_snps <- var_alpha$snps
beta_sites  <- var_beta$sites;   beta_snps  <- var_beta$snps
gamma_sites <- var_gamma$sites;  gamma_snps <- var_gamma$snps

# ── Print SNP tables ─────────────────────────────────────────────────────────
for (v in list(var_delta, var_alpha, var_beta, var_gamma)) {
  cat(sprintf("\n%s defining SNP sites:\n", v$label))
  print(data.frame(SNP = v$snps, Site = v$sites))
}

# ── Overlap with GMM-selected sites ─────────────────────────────────────────
overlap_list <- lapply(list(var_delta, var_alpha, var_beta, var_gamma),
                       function(v) intersect(selected_sites, v$sites))
names(overlap_list) <- c("Delta", "Alpha", "Beta", "Gamma")

# Keep Delta overlap as scalar for backward compatibility
overlap <- overlap_list[["Delta"]]

# voc_list defined here so it is available to all downstream chunks
# (delta-flag, pam-and-plots, mode-profile, contrast-gmm)
voc_list <- list(Delta = var_delta, Alpha = var_alpha,
                 Beta  = var_beta,  Gamma = var_gamma)

cat("\nGMM-selected sites overlapping with VOC defining SNP sites:\n")
for (nm in names(overlap_list)) {
  ov <- overlap_list[[nm]]
  cat(sprintf("  %-6s: %s\n",
              nm,
              if (length(ov) > 0) paste(ov, collapse = ", ") else "none"))
}

# ── Note on shared mutations across VOCs ─────────────────────────────────────
# Identify sites present in the SNP catalogues of more than one VOC
all_voc_sites <- lapply(voc_list, function(v) v$sites)
site_freq     <- table(unlist(all_voc_sites))
shared_sites  <- as.integer(names(site_freq[site_freq > 1L]))

cat("\nSites shared by more than one VOC catalogue:\n")
for (s in sort(shared_sites)) {
  vocs_with_site <- names(voc_list)[
    vapply(voc_list, function(v) s %in% v$sites, logical(1L))
  ]
  # Collect the actual SNP names at this site for each VOC
  snp_labels <- vapply(vocs_with_site, function(vn) {
    v   <- voc_list[[vn]]
    idx <- which(v$sites == s)
    paste(v$snps[idx], collapse = "/")
  }, character(1L))
  cat(sprintf("  Site %-5d: %s\n", s,
              paste(sprintf("%s (%s)", vocs_with_site, snp_labels),
                    collapse = ", ")))
}

# Subset to selected sites and convert to factors (required for Gower distance)
AL_Cat_sub <- df_combined[, as.character(selected_sites), drop = FALSE]
AL_Cat_sub[] <- lapply(AL_Cat_sub, as.integer)

if (UNIQUE_SEQS) {
  # Deduplicate -- record original row indices for labelling
  AL_Cat_fac    <- AL_Cat_sub %>% mutate_if(is.integer, as.factor)
  unique_idx    <- which(!duplicated(AL_Cat_fac))
  AL_Cat_fac    <- AL_Cat_fac[unique_idx, , drop = FALSE]
  orig_period   <- ifelse(unique_idx <= n_p1, 1L, 2L)
  cat(sprintf("Unique sequences: %d (from %d total)\n",
              nrow(AL_Cat_fac), n_total))
} else {
  AL_Cat_fac  <- AL_Cat_sub %>% mutate_if(is.integer, as.factor)
  orig_period <- ifelse(seq_len(n_total) <= n_p1, 1L, 2L)
  cat(sprintf("All sequences retained: %d\n", nrow(AL_Cat_fac)))
}

n_p1_eff <- sum(orig_period == 1L)
n_p2_eff <- sum(orig_period == 2L)
cat(sprintf("Period 1 sequences in PAM input: %d\n", n_p1_eff))
cat(sprintf("Period 2 sequences in PAM input: %d\n", n_p2_eff))

gower_dist <- daisy(AL_Cat_fac, metric = "gower")

future::plan(future::multisession, workers = N_WORKERS)

sil_width <- c(NA, future_sapply(2:K_MAX, function(k) {
  pam(gower_dist, diss = TRUE, k = k)$silinfo$avg.width
}))

future::plan(future::sequential)

# Silhouette plot -- grid lines, x-ticks every 5
# Drop k=1 (undefined, NA) before plotting to avoid missing-value warnings.
# na.rm = TRUE on geom calls provides an additional safety net.
sil_df <- data.frame(k = seq_len(K_MAX), width = sil_width)
sil_df <- sil_df[!is.na(sil_df$width), ]   # removes k=1 row

ggplot(sil_df, aes(x = k, y = width)) +
  geom_line(colour = "black", na.rm = TRUE) +
  geom_point(colour = "darkorange", shape = 1, size = 2, na.rm = TRUE) +
  geom_vline(xintercept = K_HIGHLIGHT, linetype = "dashed",
             colour = "steelblue", alpha = 0.6) +
  scale_x_continuous(limits       = c(2L, K_MAX),
                     breaks        = seq(5, K_MAX, by = 5),
                     minor_breaks  = seq(2, K_MAX, by = 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(title    = "Silhouette Analysis: Optimal Number of Clusters",
       subtitle  = sprintf("Blue dashed lines mark k = %s (highlighted solutions)",
                           paste(K_HIGHLIGHT, collapse = ", ")),
       x = "Number of clusters",
       y = "Average Silhouette Width") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey85"),
        panel.grid.minor = element_line(colour = "grey93"),
        plot.title       = element_text(hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, size = 9),
        axis.title.x     = element_text(colour = "steelblue"),
        axis.title.y     = element_text(colour = "darkorange"))

optimal_k <- which.max(sil_width)
cat(sprintf("Optimal k by silhouette: %d  (width = %.4f)\n",
            optimal_k, max(sil_width, na.rm = TRUE)))

set.seed(TSNE_SEED)
tsne_obj <- Rtsne(gower_dist, is_distance = TRUE,
                  perplexity = TSNE_PERPLEXITY)

tsne_base <- tsne_obj$Y %>%
  as.data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(period = orig_period)

# ── Helper: count how many of a variant's SNPs are present in a sequence ────
# seq_int: named integer vector (site name -> encoded amino acid value)
# voc:     list with $sites and $snps (from get_variant_snps())
# Returns: integer count of matching SNP alleles
count_snp_matches <- function(seq_int, voc) {
  # Only consider SNP sites that are present in the selected_sites subset
  overlap_sites <- intersect(as.character(voc$sites), names(seq_int))
  if (length(overlap_sites) == 0L) return(0L)

  expected_codes <- vapply(as.integer(overlap_sites), function(s) {
    snp_idx <- which(voc$sites == s)
    aa_char <- substring(voc$snps[snp_idx], nchar(voc$snps[snp_idx]))
    as.integer(encode_aa_sequence(matrix(aa_char))[1L, 1L])
  }, integer(1L))
  names(expected_codes) <- overlap_sites

  sum(vapply(overlap_sites, function(s) {
    isTRUE(seq_int[[s]] == expected_codes[[s]])
  }, logical(1L)))
}

# ── Apply to all Period 1 sequences ─────────────────────────────────────────
# p1_rows_int: row indices into AL_Cat_int that correspond to Period 1
p1_rows_int <- which(orig_period == 1L)

# Build integer version of the deduplicated feature matrix for display/mode calc
AL_Cat_int <- AL_Cat_fac %>%
  mutate_if(is.factor, function(x) as.integer(as.character(x)))

# Convert to named integer matrix for fast row-wise access
AL_Cat_int_mat <- as.matrix(AL_Cat_int)
colnames(AL_Cat_int_mat) <- colnames(AL_Cat_int)

# voc_list <- list(Delta = var_delta, Alpha = var_alpha,
#                  Beta  = var_beta,  Gamma = var_gamma)

# match_counts: matrix [n_p1_int x 4] -- SNP match count per Period 1 sequence
match_counts <- do.call(cbind, lapply(voc_list, function(voc) {
  vapply(p1_rows_int, function(r) {
    count_snp_matches(AL_Cat_int_mat[r, ], voc)
  }, integer(1L))
}))
rownames(match_counts) <- p1_rows_int

# ── Helper: return names of matched SNPs (not just count) ───────────────────
get_matched_snps <- function(seq_int, voc) {
  overlap_sites <- intersect(as.character(voc$sites), names(seq_int))
  if (length(overlap_sites) == 0L) return(character(0L))
  vapply(overlap_sites, function(s) {
    snp_idx   <- which(voc$sites == as.integer(s))
    aa_char   <- substring(voc$snps[snp_idx], nchar(voc$snps[snp_idx]))
    exp_code  <- as.integer(encode_aa_sequence(matrix(aa_char))[1L, 1L])
    if (isTRUE(seq_int[[s]] == exp_code)) voc$snps[snp_idx] else NA_character_
  }, character(1L)) |> (\(x) x[!is.na(x)])()
}

# ── Summary table with SNP detail ───────────────────────────────────────────
cat("Period 1 sequences by number of VOC SNP matches\n")
cat("(only SNPs at GMM-selected sites are considered)\n\n")

for (voc_nm in names(voc_list)) {
  voc         <- voc_list[[voc_nm]]
  n_overlap   <- length(intersect(as.character(voc$sites),
                                  as.character(selected_sites)))
  counts      <- match_counts[, voc_nm]
  cat(sprintf("%s  (%d defining SNPs, %d overlap with selected sites):\n",
              voc_nm, length(voc$sites), n_overlap))
  if (n_overlap == 0L) {
    cat("  No overlap with GMM-selected sites -- skip.\n\n")
    next
  }
  for (m in seq(n_overlap, 0L)) {
    seq_idx <- which(counts == m)
    if (length(seq_idx) == 0L && m < n_overlap) next
    cat(sprintf("  %d/%d SNP match(es): %d sequence(s)",
                m, n_overlap, length(seq_idx)))
    if (m > 0L && length(seq_idx) > 0L) {
      # Collect matched SNP names for each sequence at this tier
      snp_sets <- lapply(p1_rows_int[seq_idx], function(r) {
        get_matched_snps(AL_Cat_int_mat[r, ], voc)
      })
      # Unique SNP patterns
      snp_patterns <- sort(unique(vapply(snp_sets,
                                         paste, character(1L), collapse = "+")))
      cat(sprintf("  [SNP pattern(s): %s]", paste(snp_patterns, collapse = " | ")))
    }
    cat("\n")
  }
  cat("\n")
}


# Reusable plot function: takes PAM clustering vector and medoid row indices,
# produces a t-SNE plot with three visual layers:
#   1. Coloured shapes by cluster
#   2. Black cross overlay for Period 1 (wild-type) sequences
#   3. Filled red diamond for each cluster medoid
plot_tsne <- function(pam_clust, k, medoid_rows,
                      tsne_df    = tsne_base,
                      period_vec = orig_period,
                      xlim_fixed = NULL,
                      ylim_fixed = NULL) {

  df <- tsne_df %>%
    mutate(cluster   = factor(pam_clust),
           is_p1     = (period_vec == 1L),
           is_medoid = row_number() %in% medoid_rows)

  p <- ggplot(df, aes(x = X, y = Y, colour = cluster, shape = cluster)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_shape_manual(values = seq(0, k)) +
    # Black cross overlay for all Period 1 sequences
    geom_point(data         = subset(df, is_p1),
               aes(x = X, y = Y),
               colour       = "black", shape = 4, size = 2,
               inherit.aes  = FALSE) +
    # Semi-transparent red diamond: alpha reveals the medoid is a real data point.
    geom_point(data         = subset(df, is_medoid),
               aes(x = X, y = Y),
               colour       = "red", shape = 18, size = 5,
               alpha        = 0.45,
               inherit.aes  = FALSE) +
    labs(title    = sprintf("PAM k = %d  |  t-SNE visualisation", k),
         subtitle = "Crosses = Period 1 (wild-type);  Red diamonds = cluster medoids",
         x = "t-SNE dimension 1",
         y = "t-SNE dimension 2") +
    theme_bw() +
    theme(plot.title    = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 9))

  # Apply fixed axis ranges if supplied -- ensures identical plot area across
  # all k values regardless of legend size differences
  if (!is.null(xlim_fixed) && !is.null(ylim_fixed))
    p <- p + coord_cartesian(xlim = xlim_fixed, ylim = ylim_fixed)

  print(p)
  invisible(p)
}

# Helper: render kableExtra table in results='asis' chunks inside loops
kt <- function(tbl) cat(as.character(tbl), "\n")
# Fit PAM for all highlighted k values
pam_fits <- lapply(K_HIGHLIGHT, function(k) {
  pam(gower_dist, diss = TRUE, k = k)
})
names(pam_fits) <- paste0("k", K_HIGHLIGHT)

# Precompute fixed axis ranges from the common t-SNE coordinates so that all
# plots have an identical data area regardless of legend size differences.
x_pad    <- diff(range(tsne_base$X)) * 0.05
y_pad    <- diff(range(tsne_base$Y)) * 0.05
xlim_fix <- range(tsne_base$X) + c(-x_pad, x_pad)
ylim_fix <- range(tsne_base$Y) + c(-y_pad, y_pad)

for (k in K_HIGHLIGHT) {
  fit         <- pam_fits[[paste0("k", k)]]
  medoid_rows <- fit$id.med

  # ── Subsection heading ──────────────────────────────────────────────────
  cat(sprintf("\n## PAM k = %d\n\n### t-SNE Plot\n\n", k))

  # ── Plot ─────────────────────────────────────────────────────────────────
  plot_tsne(fit$clustering, k, medoid_rows,
            xlim_fixed = xlim_fix, ylim_fixed = ylim_fix)

  # ── Medoid sequences kable ───────────────────────────────────────────────
  cat(sprintf("\n### Medoid Sequences (k = %d)\n\n", k))
  medoid_aa <- decode_aa_sequence(
    as.matrix(AL_Cat_int[medoid_rows, , drop = FALSE])
  )
  medoid_df           <- as.data.frame(medoid_aa)
  colnames(medoid_df) <- as.character(selected_sites)
  rownames(medoid_df) <- paste0("cl_", seq_len(k))

  # Red cell_spec wherever cl_i (i > 1) differs from cl_1 at that site.
  # Direct [row, col_name] access -- avoids any as.character(df[row,]) pitfalls.
  if (k > 1L) {
    for (ri in 2L:k) {
      for (nm in colnames(medoid_df)) {
        v1 <- medoid_df[1L, nm]; vi <- medoid_df[ri, nm]
        if (vi != v1)
          medoid_df[ri, nm] <- cell_spec(vi, format = "html",
                                         color = "red", bold = TRUE)
      }
    }
  }

  kt(
    kbl(medoid_df, format = "html", escape = FALSE,
        caption = sprintf("Medoid sequences (k = %d) | red = differs from cl_1", k)) |>
      kable_styling(bootstrap_options = c("condensed", "hover"),
                    font_size = 11, full_width = FALSE) |>
      row_spec(1L, background = "#d6eaf8")
  )

  # ── VOC SNP content kable ────────────────────────────────────────────────
  cat(sprintf("\n### VOC SNP Content of Medoids (k = %d)\n\n", k))

  voc_tbl <- do.call(rbind, lapply(seq_len(k), function(cl_i) {
    row_i <- medoid_rows[cl_i]
    vapply(names(voc_list), function(voc_nm) {
      voc   <- voc_list[[voc_nm]]
      n_ov  <- length(intersect(as.character(voc$sites),
                                as.character(selected_sites)))
      n_m   <- count_snp_matches(AL_Cat_int_mat[row_i, ], voc)
      sprintf("%d / %d", n_m, n_ov)
    }, character(1L))
  }))
  voc_tbl           <- as.data.frame(voc_tbl)
  rownames(voc_tbl) <- paste0("cl_", seq_len(k))

  kt(
    kbl(voc_tbl, format = "html",
        caption = sprintf("VOC SNP matches per medoid (matched / overlap sites), k = %d", k)) |>
      kable_styling(bootstrap_options = c("condensed", "hover"),
                    font_size = 11, full_width = FALSE)
  )
  cat("\n")
}

get_mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

# Store mode matrices for Step 7 contrast analysis.
mode_results <- list()

for (k in K_HIGHLIGHT) {
  fit         <- pam_fits[[paste0("k", k)]]
  clust_vec   <- fit$clustering
  medoid_rows <- fit$id.med
  k_levels    <- sort(unique(clust_vec))

  mode_mat <- do.call(rbind, lapply(k_levels, function(cl) {
    rows <- which(clust_vec == cl)
    vapply(colnames(AL_Cat_int), function(s) get_mode(AL_Cat_int[rows, s]),
           integer(1L))
  }))
  rownames(mode_mat) <- paste0("cl_", k_levels)
  colnames(mode_mat) <- colnames(AL_Cat_int)

  mode_aa           <- decode_aa_sequence(as.matrix(mode_mat))
  rownames(mode_aa) <- paste0("cl_", k_levels)
  colnames(mode_aa) <- as.character(selected_sites)

  medoid_int          <- as.matrix(AL_Cat_int[medoid_rows, , drop = FALSE])
  medoid_aa           <- decode_aa_sequence(medoid_int)
  rownames(medoid_aa) <- paste0("med_", k_levels)
  colnames(medoid_aa) <- as.character(selected_sites)

  mode_results[[paste0("k", k)]] <- list(
    mode_int  = mode_mat, k_levels = k_levels, clust_vec = clust_vec)

  same_flag <- all(mode_mat == medoid_int)
  obs_note  <- if (same_flag)
    sprintf("**Mode = Medoid for all clusters (k = %d): tight cohesion confirmed.**", k)
  else {
    diff_s <- which(apply(mode_mat != medoid_int, 2, any))
    sprintf("**Sites where mode differs from medoid (k = %d): %s**",
            k, paste(selected_sites[diff_s], collapse = ", "))
  }

  cat(sprintf("\n### k = %d\n\n%s\n\n", k, obs_note))

  # Build combined df with cell_spec red highlighting where mode != medoid
  combined_df <- as.data.frame(rbind(mode_aa, medoid_aa), stringsAsFactors = FALSE)
  n_cl <- length(k_levels)
  for (col_i in seq_along(selected_sites)) {
    col_nm <- as.character(selected_sites[col_i])
    for (ci in seq_along(k_levels)) {
      mv <- mode_aa[ci, col_i]; dv <- medoid_aa[ci, col_i]
      if (mv != dv) {
        combined_df[ci,        col_nm] <- cell_spec(mv, color="red", bold=TRUE)
        combined_df[n_cl + ci, col_nm] <- cell_spec(dv, color="red", bold=TRUE)
      }
    }
  }

  kt(
    kbl(combined_df, format="html", escape=FALSE,
        caption=sprintf("Mode (cl_* = blue) / Medoid (med_* = yellow) | red = differs  [k=%d]", k)) |>
      kable_styling(bootstrap_options=c("condensed","hover"), font_size=11, full_width=FALSE) |>
      row_spec(seq_len(n_cl),           background="#d6eaf8") |>
      row_spec(seq(n_cl+1L, n_cl*2L),   background="#fef9e7")
  )
  cat("\n")
}

compute_metrics <- function(clust_vec, period_vec, n_p1_eff) {
  p1_idx <- which(period_vec == 1L)
  cl_purity <- tapply(clust_vec[p1_idx], clust_vec[p1_idx], length)
  positive_cluster <- as.integer(names(which.max(cl_purity)))
  pred_positive <- clust_vec == positive_cluster
  true_positive <- period_vec == 1L
  TP <- sum( pred_positive &  true_positive)
  FP <- sum( pred_positive & !true_positive)
  FN <- sum(!pred_positive &  true_positive)
  TN <- sum(!pred_positive & !true_positive)
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  recall    <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  f1        <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
               2 * precision * recall / (precision + recall) else NA_real_
  list(k = length(unique(clust_vec)), positive_cluster = positive_cluster,
       TP = TP, FP = FP, FN = FN, TN = TN,
       precision = round(precision, 4), recall = round(recall, 4),
       f1 = round(f1, 4))
}

metrics_list <- lapply(K_HIGHLIGHT, function(k) {
  compute_metrics(pam_fits[[paste0("k", k)]]$clustering, orig_period, n_p1_eff)
})
metrics_df <- do.call(rbind, lapply(metrics_list, as.data.frame))
print(metrics_df)

# annotate_site_full(s): compact VOC annotation for Spike site s.
# Uses ALL 12 variants in sarscov2_variants (not just the 4 in voc_list).
# Distinguishes:
#   "VOC-SNP"  = site is a defining SNP for that VOC (Defining_SNP_Sites)
#   "VOC"      = site is mutated in that VOC but is NOT a defining SNP
#   ""         = not found in any VOC catalogue
# Examples: 452 -> "Alpha/Delta/Epsilon/Iota/Kappa-SNP"
#            95 -> "Delta/Iota/Kappa/Omicron"
#           190 -> "Gamma-SNP; Iota"
annotate_site_full <- function(s) {
  labels    <- unlist(sarscov2_variants$WHO_Label)
  snp_sites <- sarscov2_variants$Defining_SNP_Sites
  mut_sites <- sarscov2_variants$Mutation_Sites
  snp_vocs  <- labels[vapply(snp_sites, function(ss) s %in% ss, logical(1L))]
  mut_only  <- setdiff(labels[vapply(mut_sites, function(ms) s %in% ms, logical(1L))],
                       snp_vocs)
  parts <- c(
    if (length(snp_vocs) > 0L) paste0(paste(snp_vocs, collapse = "/"), "-SNP"),
    if (length(mut_only)  > 0L) paste(mut_only, collapse = "/"))
  if (length(parts) == 0L) "" else paste(parts, collapse = "; ")
}

max_possible <- sum(vapply(K_HIGHLIGHT, function(k) {
  length(unique(pam_fits[[paste0("k", k)]]$clustering)) - 1L
}, integer(1L)))

site_diff_counts <- integer(length(selected_sites))
names(site_diff_counts) <- as.character(selected_sites)
contrast_rows <- list()

for (k in K_HIGHLIGHT) {
  mr        <- mode_results[[paste0("k", k)]]
  clust_vec <- mr$clust_vec
  k_levels  <- mr$k_levels
  mode_int  <- mr$mode_int

  p1_idx <- which(orig_period == 1L)
  purity <- tapply(clust_vec[p1_idx], clust_vec[p1_idx], length)
  wt_cl  <- as.integer(names(which.max(purity)))
  wt_row <- which(k_levels == wt_cl)

  for (cl in k_levels[k_levels != wt_cl]) {
    cl_row        <- which(k_levels == cl)
    diff_vec      <- mode_int[cl_row, ] != mode_int[wt_row, ]
    diff_sites_cl <- selected_sites[diff_vec]

    site_diff_counts[as.character(diff_sites_cl)] <-
      site_diff_counts[as.character(diff_sites_cl)] + 1L

    voc_str <- if (length(diff_sites_cl) > 0L)
      vapply(diff_sites_cl, annotate_site_full, character(1L))
    else character(0L)

    contrast_rows[[length(contrast_rows)+1L]] <- data.frame(
      k=k, wt=wt_cl, vs=cl, n=length(diff_sites_cl),
      sites = if (length(diff_sites_cl)>0L) paste(diff_sites_cl, collapse=", ") else "-",
      VOC   = if (length(diff_sites_cl)>0L)
        paste(sprintf("%s[%s]", diff_sites_cl, voc_str), collapse=", ") else "-",
      stringsAsFactors=FALSE)
  }
}

# Ranked frequency table
freq_df <- data.frame(
  Site = as.integer(names(site_diff_counts)),
  Freq = as.integer(site_diff_counts),
  Pct  = round(100 * site_diff_counts / max_possible, 1),
  stringsAsFactors = FALSE)
freq_df <- freq_df[freq_df$Freq > 0L, ]
if (nrow(freq_df) > 0L) {
  freq_df <- freq_df[order(-freq_df$Freq), ]
  freq_df$VOC <- vapply(freq_df$Site, annotate_site_full, character(1L))
  cat(sprintf("\n**Site differentiation across all %d (k, cluster) comparisons:**\n\n",
              max_possible))
  kt(
    kbl(freq_df, format="html",
        col.names=c("Site","Freq",sprintf("%% of %d comparisons", max_possible),"VOC"),
        caption="Sites ranked by differentiation frequency (wild-type vs other clusters)") |>
      kable_styling(bootstrap_options=c("condensed","hover","striped"),
                    font_size=11, full_width=FALSE) |>
      column_spec(3,
        color = ifelse(freq_df$Pct >= 50, "darkred", "black"),
        bold  = freq_df$Pct >= 50)
  )
  cat("\n")
} else {
  cat("\n*No differentiating sites found (all clusters share the same modal residues).*\n\n")
}

# Per-k breakdown
if (length(contrast_rows) > 0L) {
  contrast_df <- do.call(rbind, contrast_rows)
  cat("\n**Per-(k, cluster) breakdown:**\n\n")
  kt(
    kbl(contrast_df, format="html",
        caption="Differentiating sites: wild-type cluster vs each other cluster") |>
      kable_styling(bootstrap_options=c("condensed","hover"), font_size=11, full_width=FALSE)
  )
}

# Borderline analysis for k=4 only. Boundary regions from visual inspection:
#   Region A: X > 15  AND  Y <= 0          (right border, near cluster 2)
#   Region B: X > -10 AND  X < 10  AND  Y < -10  (center bottom border, clusters 3/4)
# Applies to ALL cluster-1 sequences regardless of time period.

k         <- 4L
fit       <- pam_fits[[paste0("k", k)]]
clust_vec <- fit$clustering
med_rows  <- fit$id.med

wt_cl       <- 1L
foreign_cls <- setdiff(sort(unique(clust_vec)), wt_cl)
wt_rows     <- which(clust_vec == wt_cl)

cat(sprintf("\n### k = %d  (wild-type = cluster %d)\n\n", k, wt_cl))

# Apply coordinate filter
region_A <- tsne_base$X[wt_rows] > 15  & tsne_base$Y[wt_rows] <= 0
region_B <- tsne_base$X[wt_rows] > -10 & tsne_base$X[wt_rows] < 10 &
            tsne_base$Y[wt_rows] < -10
is_bl   <- region_A | region_B
bl_rows <- wt_rows[is_bl]

cat(sprintf("%d wt-cluster sequences fall in the boundary t-SNE region(s).\n\n",
            length(bl_rows)))

if (length(bl_rows) == 0L) {
  cat("*No sequences in the defined boundary regions.*\n\n")
} else {

  # t-SNE distances to own and nearest foreign medoid
  wt_med_x <- tsne_base$X[med_rows[wt_cl]]
  wt_med_y <- tsne_base$Y[med_rows[wt_cl]]
  d_own    <- sqrt((tsne_base$X[bl_rows] - wt_med_x)^2 +
                   (tsne_base$Y[bl_rows] - wt_med_y)^2)

  fgn_dist_mat <- vapply(foreign_cls, function(cl) {
    fx <- tsne_base$X[med_rows[cl]]; fy <- tsne_base$Y[med_rows[cl]]
    sqrt((tsne_base$X[bl_rows] - fx)^2 + (tsne_base$Y[bl_rows] - fy)^2)
  }, numeric(length(bl_rows)))

  d_fgn  <- if (length(foreign_cls) == 1L) as.numeric(fgn_dist_mat)
             else apply(fgn_dist_mat, 1L, min)
  margin <- d_own - d_fgn

  period_lbl <- ifelse(orig_period[bl_rows] == 1L, "P1", "P2")

  dist_tbl <- data.frame(
    global_idx = bl_rows,
    period     = period_lbl,
    tsne_X     = round(tsne_base$X[bl_rows], 1),
    tsne_Y     = round(tsne_base$Y[bl_rows], 1),
    d_own      = round(d_own,  2),
    d_fgn      = round(d_fgn, 2),
    margin     = round(margin, 2)
  )
  dist_tbl <- dist_tbl[order(dist_tbl$margin), ]

  kt(
    kbl(dist_tbl, format = "html",
        col.names = c("Row idx", "Period", "t-SNE X", "t-SNE Y",
                      "d(own med)", "d(nearest fgn med)", "Margin"),
        caption = sprintf(
          "Boundary-region wt-cluster sequences, sorted by margin [k=%d, wt=cl%d]",
          k, wt_cl)) |>
      kable_styling(bootstrap_options = c("condensed","hover"),
                    font_size = 11, full_width = FALSE) |>
      column_spec(7L,
                  color = ifelse(dist_tbl$margin < 0, "red", "black"),
                  bold  = dist_tbl$margin < 0)
  )
  cat("\n")

  # VOC SNP content table
  n_ref         <- 1L + length(foreign_cls)
  all_rows_bl   <- c(med_rows[wt_cl], med_rows[foreign_cls], dist_tbl$global_idx)
  row_labels_bl <- c(
    sprintf("wt_med (cl%d)", wt_cl),
    sprintf("fgn_med (cl%d)", foreign_cls),
    sprintf("border[row%d] %s", dist_tbl$global_idx, dist_tbl$period)
  )
  voc_bl <- do.call(rbind, lapply(all_rows_bl, function(r) {
    vapply(names(voc_list), function(voc_nm) {
      voc  <- voc_list[[voc_nm]]
      n_ov <- length(intersect(as.character(voc$sites),
                               as.character(selected_sites)))
      n_m  <- count_snp_matches(AL_Cat_int_mat[r, ], voc)
      sprintf("%d / %d", n_m, n_ov)
    }, character(1L))
  }))
  voc_bl           <- as.data.frame(voc_bl)
  rownames(voc_bl) <- row_labels_bl
  n_bl             <- nrow(voc_bl) - n_ref

  tbl_voc <- kbl(voc_bl, format = "html",
      caption = sprintf(
        "VOC SNP content: wt medoid | foreign medoids | borderline sequences [k=%d]",
        k)) |>
    kable_styling(bootstrap_options = c("condensed","hover"),
                  font_size = 11, full_width = FALSE) |>
    row_spec(1L,             background = "#d6eaf8") |>
    row_spec(seq(2L, n_ref), background = "#fef9e7")
  if (n_bl > 0L)
    tbl_voc <- row_spec(tbl_voc, seq(n_ref + 1L, nrow(voc_bl)),
                        background = "#fce4ec")
  kt(tbl_voc)
  cat("\n")

  # Sequence comparison table
  own_med_aa <- as.data.frame(decode_aa_sequence(
    as.matrix(AL_Cat_int[med_rows[wt_cl], , drop = FALSE])))
  foreign_aa <- do.call(rbind, lapply(foreign_cls, function(cl)
    as.data.frame(decode_aa_sequence(
      as.matrix(AL_Cat_int[med_rows[cl], , drop = FALSE])))))
  border_aa  <- as.data.frame(decode_aa_sequence(
    as.matrix(AL_Cat_int[bl_rows, , drop = FALSE])))

  compare_df <- rbind(own_med_aa, foreign_aa, border_aa)
  colnames(compare_df) <- as.character(selected_sites)
  rownames(compare_df) <- c(
    sprintf("wt_med (cl%d)", wt_cl),
    sprintf("fgn_med (cl%d)", foreign_cls),
    sprintf("border[row%d] %s margin=%.2f",
            dist_tbl$global_idx, dist_tbl$period, dist_tbl$margin)
  )

  out_df <- compare_df
  for (ri in 2L:nrow(compare_df)) {
    for (nm in colnames(compare_df)) {
      v1 <- compare_df[1L, nm]; vi <- compare_df[ri, nm]
      if (!is.na(vi) && vi != v1)
        out_df[ri, nm] <- cell_spec(vi, format = "html", color = "red", bold = TRUE)
    }
  }
  kt(
    kbl(out_df, format = "html", escape = FALSE,
        caption = sprintf(
          "wt medoid (blue) | foreign medoids (yellow) | borderline (pink) | red = differs from wt medoid  [k=%d]",
          k)) |>
      kable_styling(bootstrap_options = c("condensed","hover"),
                    font_size = 11, full_width = FALSE) |>
      row_spec(1L,             background = "#d6eaf8") |>
      row_spec(seq(2L, n_ref), background = "#fef9e7") |>
      row_spec(seq(n_ref + 1L, nrow(out_df)), background = "#fce4ec")
  )
  cat("\n")
}

# Entropy on all N_SITES sites across all sequences in AL_df
entrp_full <- apply(AL_df[, seq_len(N_SITES), drop = FALSE], 2L, calculate_entropy)

clust_full    <- cluster_sites_by_entropy(entrp_full, nr = nrow(AL_df))
clust_full_rl <- relabel_entropy_classes(clust_full$DataFrame)

selected_sites_full <- sort(
  clust_full_rl[clust_full_rl$class == 1L, ]$sites
)

cat(sprintf("Full-dataset GMM selected %d sites from %d sequences:\n",
            length(selected_sites_full), nrow(AL_df)))
print(selected_sites_full)

# Overlap with VOC defining SNP sites -- uses same get_variant_snps() as Step 2.
# The difference vs. the two-period selected_sites is expected: entropy was
# computed on 109k sequences with different frequency distributions, so GMM
# selects a different (potentially overlapping) site subset.
# Delta two-period: 5 overlap sites (19,452,478,681,950).
# Delta full-dataset overlap depends on selected_sites_full computed above.
overlap_full <- lapply(voc_list, function(v)
  intersect(selected_sites_full, v$sites))
names(overlap_full) <- names(voc_list)

cat("\nGMM-selected sites overlapping with VOC defining SNP sites (full dataset):\n")
for (nm in names(overlap_full)) {
  ov <- overlap_full[[nm]]
  cat(sprintf("  %-6s: %s\n", nm,
              if (length(ov) > 0L) paste(ov, collapse = ", ") else "none"))
}

# Subset to selected sites and convert to factors
AL_full_fac <- AL_df[, as.character(selected_sites_full), drop = FALSE]
AL_full_fac[] <- lapply(AL_full_fac, function(x) as.factor(as.integer(x)))

# Deduplicate before Gower: 109k sequences produce an infeasible distance
# matrix (~48 GB). Unique sequences are sufficient for PAM/t-SNE and
# dramatically reduce computation while preserving the full sequence diversity.
full_unique_idx  <- which(!duplicated(AL_full_fac))
AL_full_fac_uniq <- AL_full_fac[full_unique_idx, , drop = FALSE]

cat(sprintf("Unique sequences for full-dataset clustering: %d (from %d total)\n",
            nrow(AL_full_fac_uniq), nrow(AL_full_fac)))

gower_dist_full <- daisy(AL_full_fac_uniq, metric = "gower")

future::plan(future::multisession, workers = N_WORKERS)

sil_width_full <- c(NA, future_sapply(2:K_MAX, function(k) {
  pam(gower_dist_full, diss = TRUE, k = k)$silinfo$avg.width
}))

future::plan(future::sequential)

sil_df_full <- data.frame(k = seq_len(K_MAX), width = sil_width_full)
sil_df_full <- sil_df_full[!is.na(sil_df_full$width), ]

ggplot(sil_df_full, aes(x = k, y = width)) +
  geom_line(colour = "black") +
  geom_point(colour = "darkorange", shape = 1, size = 2) +
  geom_vline(xintercept = c(6L, 7L, 9L), linetype = "dashed", colour = "steelblue", alpha = 0.6) +
  scale_x_continuous(limits = c(2L, K_MAX),
                     breaks = seq(5, K_MAX, by = 5),
                     minor_breaks = seq(2, K_MAX, by = 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  labs(title    = "Silhouette Analysis: Full Dataset (All Sequences)",
       subtitle = "Blue dashed lines mark k = 6, 7, 9",
       x = "Number of clusters", y = "Average Silhouette Width") +
  theme_bw() +
  theme(panel.grid.major  = element_line(colour = "grey85"),
        panel.grid.minor  = element_line(colour = "grey93"),
        plot.title        = element_text(hjust = 0.5),
        plot.subtitle     = element_text(hjust = 0.5, size = 9),
        axis.title.x      = element_text(colour = "steelblue"),
        axis.title.y      = element_text(colour = "darkorange"))

optimal_k_full <- which.max(sil_width_full)
cat(sprintf("Optimal k by silhouette (full dataset): %d  (width = %.4f)\n",
            optimal_k_full, max(sil_width_full, na.rm = TRUE)))

K_HIGHLIGHT_FULL <- c(6L, 7L, 9L)

pam_fits_full <- lapply(K_HIGHLIGHT_FULL, function(k)
  pam(gower_dist_full, diss = TRUE, k = k))
names(pam_fits_full) <- paste0("k", K_HIGHLIGHT_FULL)

for (k in K_HIGHLIGHT_FULL) {
  cat(sprintf("k=%d cluster sizes: ", k))
  cat(paste(table(pam_fits_full[[paste0("k",k)]]$clustering), collapse=" | "), "\n")
}

set.seed(TSNE_SEED)
tsne_full_obj <- Rtsne(gower_dist_full, is_distance = TRUE,
                       perplexity = TSNE_PERPLEXITY)

tsne_full_base <- as.data.frame(tsne_full_obj$Y) |>
  setNames(c("X", "Y"))

# Integer matrix for medoid decoding and VOC matching (built once, reused)
AL_full_int     <- AL_full_fac_uniq |>
  mutate_if(is.factor, function(x) as.integer(as.character(x)))
AL_full_int_mat <- as.matrix(AL_full_int)

# VOC list for full dataset (Alpha, Beta, Delta, Gamma only)
voc_list_full <- list(
  Delta = get_variant_snps("Delta"),
  Alpha = get_variant_snps("Alpha"),
  Beta  = get_variant_snps("Beta"),
  Gamma = get_variant_snps("Gamma")
)

# Fixed axis ranges for identical plot area across all k values
x_pad_f    <- diff(range(tsne_full_base$X)) * 0.05
y_pad_f    <- diff(range(tsne_full_base$Y)) * 0.05
xlim_fix_f <- range(tsne_full_base$X) + c(-x_pad_f,  x_pad_f)
ylim_fix_f <- range(tsne_full_base$Y) + c(-y_pad_f,  y_pad_f)

for (k in K_HIGHLIGHT_FULL) {
  fit_f      <- pam_fits_full[[paste0("k", k)]]
  med_f      <- fit_f$id.med
  cl_sizes_f <- as.integer(table(fit_f$clustering))

  # ── t-SNE plot ─────────────────────────────────────────────────────────────
  cat(sprintf("\n## Full Dataset PAM k = %d\n\n### t-SNE Plot\n\n", k))

  df_f <- tsne_full_base |>
    mutate(cluster   = factor(fit_f$clustering),
           is_medoid = row_number() %in% med_f)

  p_f <- ggplot(df_f, aes(x = X, y = Y, colour = cluster, shape = cluster)) +
    geom_point(size = 1.5, alpha = 0.5) +
    scale_shape_manual(values = seq(0, k)) +
    geom_point(data        = subset(df_f, is_medoid),
               aes(x = X, y = Y),
               colour      = "red", shape = 18, size = 5, alpha = 0.45,
               inherit.aes = FALSE) +
    coord_cartesian(xlim = xlim_fix_f, ylim = ylim_fix_f) +
    labs(title    = sprintf("PAM k = %d | Full Dataset t-SNE", k),
         subtitle = sprintf("Unique sequences: %d | Red diamonds = cluster medoids",
                            nrow(AL_full_fac_uniq)),
         x = "t-SNE dimension 1", y = "t-SNE dimension 2") +
    theme_bw() +
    theme(plot.title    = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
  print(p_f)

  # ── Medoid sequences table ──────────────────────────────────────────────────
  cat(sprintf("\n### Medoid Sequences (k = %d)\n\n", k))

  med_aa <- decode_aa_sequence(as.matrix(AL_full_int[med_f, , drop = FALSE]))
  med_df <- as.data.frame(med_aa)
  colnames(med_df) <- as.character(selected_sites_full)
  rownames(med_df) <- paste0("cl_", seq_len(k))
  med_df <- cbind(N = cl_sizes_f, med_df)

  if (k > 1L) {
    for (ri in 2L:k) {
      for (nm in as.character(selected_sites_full)) {
        v1 <- med_df[1L, nm]; vi <- med_df[ri, nm]
        if (!is.na(vi) && vi != v1)
          med_df[ri, nm] <- cell_spec(vi, format = "html", color = "red", bold = TRUE)
      }
    }
  }
  kt(
    kbl(med_df, format = "html", escape = FALSE,
        caption = sprintf(
          "Medoid sequences (k=%d, full dataset) | N = cluster size | red = differs from cl_1",
          k)) |>
      kable_styling(bootstrap_options = c("condensed","hover"),
                    font_size = 11, full_width = FALSE) |>
      row_spec(1L, background = "#d6eaf8")
  )

  # ── VOC SNP content table ───────────────────────────────────────────────────
  cat(sprintf("\n### VOC SNP Content (k = %d)\n\n", k))

  voc_tbl_f <- do.call(rbind, lapply(seq_len(k), function(cl_i) {
    row_i <- med_f[cl_i]
    vapply(names(voc_list_full), function(voc_nm) {
      voc  <- voc_list_full[[voc_nm]]
      n_ov <- length(intersect(as.character(voc$sites),
                               as.character(selected_sites_full)))
      n_m  <- count_snp_matches(AL_full_int_mat[row_i, ], voc)
      sprintf("%d / %d", n_m, n_ov)
    }, character(1L))
  }))
  voc_tbl_f           <- as.data.frame(voc_tbl_f)
  rownames(voc_tbl_f) <- paste0("cl_", seq_len(k))

  kt(
    kbl(voc_tbl_f, format = "html",
        caption = sprintf(
          "VOC SNP matches per medoid (matched / overlap sites), full dataset k=%d",
          k)) |>
      kable_styling(bootstrap_options = c("condensed","hover"),
                    font_size = 11, full_width = FALSE)
  )
  cat("\n")
}

cat("== Analysis Summary ==\n\n")
cat(sprintf("Country filter:         %s\n",
            ifelse(is.null(COUNTRY), "none (all countries)", COUNTRY)))
cat(sprintf("Unique sequences:       %s\n", UNIQUE_SEQS))
cat(sprintf("Period 1 sequences:     %d  (%s to %s)\n",
            n_p1_eff, PERIOD1_START, PERIOD1_END))
cat(sprintf("Period 2 sequences:     %d  (%s to %s)\n",
            n_p2_eff, PERIOD2_START, PERIOD2_END))
cat(sprintf("GMM-selected sites:     %d  [%s]\n",
            length(selected_sites),
            paste(selected_sites, collapse = ", ")))
cat("\nVOC SNP overlap with GMM-selected sites:\n")
for (nm in names(overlap_list)) {
  ov <- overlap_list[[nm]]
  cat(sprintf("  %-6s: %s\n",
              nm,
              if (length(ov) > 0) paste(ov, collapse = ", ") else "none"))
}
cat(sprintf("\nOptimal k by silhouette:   %d\n", optimal_k))
cat(sprintf("Optimal k by F1 score:     %d\n\n",
            K_HIGHLIGHT[which.max(metrics_df$f1)]))
cat("Classification metrics (k = 2, 3, 4, 9):\n")
print(metrics_df[, c("k", "positive_cluster", "TP", "FP", "FN",
                      "precision", "recall", "f1")])
