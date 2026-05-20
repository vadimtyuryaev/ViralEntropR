################################################################################
## aggregate_tables.R — Changepoint_Detection_Study
##
## Produces four CSV tables for the manuscript and supplementary materials.
##
##   tab_01_method_tier_dataset_summary.csv
##       One row per (dataset, tier, method, site_set). Columns: F1 mean,
##       median, IQR, TLE mean, n_detected mean, failure_rate, n_reps.
##
##   tab_02_operating_curve.csv
##       One row per (dataset, method, truths_effective). Columns:
##       F1 mean, F1 95% CI, n_reps. Long format suitable for re-plotting.
##
##   tab_03_k_sweep_saturation.csv
##       One row per (dataset, tier, K). Columns: F1 mean, F1 95% CI, n_reps.
##
##   tab_04_full_vs_reduced_delta.csv
##       One row per (dataset, tier, method). Columns: F1 full, F1 reduced,
##       paired delta (full - reduced), Wilcoxon signed-rank p-value, n_pairs.
##
## All numeric columns are rounded to 4 decimal places.
## All CSVs are UTF-8 with comma separators and headers.
################################################################################

suppressPackageStartupMessages({
  source("setup.R")
})

config    <- build_config()
table_dir <- file.path(output_dir(config), "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("Loading benchmark summary tables...")
sN <- readRDS(summary_path(config, "NCBI_US"))
sG <- readRDS(summary_path(config, "GISAID_US"))
sN$dataset_label <- "NCBI_US"
sG$dataset_label <- "GISAID_US"
sAll <- rbind(sN, sG)

# Convenience: 95% CI half-width via 1.96 × SE.
ci_half <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2L) return(NA_real_)
  1.96 * sd(x) / sqrt(length(x))
}

write_csv_atomic <- function(df, path) {
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  utils::write.csv(df, tmp, row.names = FALSE, quote = TRUE,
                   na = "NA", fileEncoding = "UTF-8")
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    file.remove(tmp)
  }
  log_msg(sprintf("Saved: %s (%d rows)", path, nrow(df)))
}

# ── 1. tab_01_method_tier_dataset_summary ────────────────────────────────────

log_msg("Building tab_01: method × tier × dataset summary...")
agg_by <- function(df, FUN) {
  aggregate(
    list(value = df$F1),
    list(dataset = df$dataset_label, tier = df$tier,
         method = df$method, site_set = df$site_set),
    FUN, na.rm = TRUE
  )
}

agg_mean   <- agg_by(sAll, mean);   names(agg_mean)[5L]   <- "F1_mean"
agg_median <- agg_by(sAll, median); names(agg_median)[5L] <- "F1_median"
agg_q25    <- agg_by(sAll, function(x, na.rm) quantile(x, 0.25, na.rm = na.rm))
names(agg_q25)[5L]    <- "F1_Q25"
agg_q75    <- agg_by(sAll, function(x, na.rm) quantile(x, 0.75, na.rm = na.rm))
names(agg_q75)[5L]    <- "F1_Q75"

agg_tle <- aggregate(
  list(TLE_mean = sAll$TLE),
  list(dataset = sAll$dataset_label, tier = sAll$tier,
       method = sAll$method, site_set = sAll$site_set),
  mean, na.rm = TRUE
)

agg_n <- aggregate(
  list(n_detected_mean = sAll$n_detected),
  list(dataset = sAll$dataset_label, tier = sAll$tier,
       method = sAll$method, site_set = sAll$site_set),
  mean
)

agg_fail <- aggregate(
  list(failure_rate = sAll$status == "failed"),
  list(dataset = sAll$dataset_label, tier = sAll$tier,
       method = sAll$method, site_set = sAll$site_set),
  mean
)

agg_count <- aggregate(
  list(n_reps = sAll$F1),
  list(dataset = sAll$dataset_label, tier = sAll$tier,
       method = sAll$method, site_set = sAll$site_set),
  length
)

keys <- c("dataset", "tier", "method", "site_set")
tab_01 <- Reduce(function(a, b) merge(a, b, by = keys, all = TRUE),
                 list(agg_mean, agg_median, agg_q25, agg_q75,
                      agg_tle,  agg_n,      agg_fail, agg_count))

# Round to 4 decimals where applicable.
num_cols <- c("F1_mean", "F1_median", "F1_Q25", "F1_Q75",
              "TLE_mean", "n_detected_mean", "failure_rate")
for (col in num_cols)
  if (col %in% colnames(tab_01))
    tab_01[[col]] <- round(tab_01[[col]], 4)

tab_01 <- tab_01[order(tab_01$dataset, tab_01$tier,
                       tab_01$site_set, tab_01$method), , drop = FALSE]
write_csv_atomic(tab_01,
                 file.path(table_dir, "tab_01_method_tier_dataset_summary.csv"))

# ── 2. tab_02_operating_curve ────────────────────────────────────────────────

log_msg("Building tab_02: operating curve F1 vs truths_effective...")
op_df <- subset(sAll, site_set == "full" &
                       method %in% c("e_agglo", "ks_cp3o_dynamic"))

tab_02_mean <- aggregate(
  list(F1_mean = op_df$F1),
  list(dataset = op_df$dataset_label, method = op_df$method,
       truths_effective = op_df$truths_effective),
  mean, na.rm = TRUE
)
tab_02_ci <- aggregate(
  list(F1_CI_half = op_df$F1),
  list(dataset = op_df$dataset_label, method = op_df$method,
       truths_effective = op_df$truths_effective),
  ci_half
)
tab_02_n <- aggregate(
  list(n_reps = op_df$F1),
  list(dataset = op_df$dataset_label, method = op_df$method,
       truths_effective = op_df$truths_effective),
  function(x) sum(!is.na(x))
)
tab_02 <- Reduce(function(a, b) merge(a, b,
                                      by = c("dataset", "method",
                                             "truths_effective"),
                                      all = TRUE),
                 list(tab_02_mean, tab_02_ci, tab_02_n))
tab_02$F1_mean    <- round(tab_02$F1_mean, 4)
tab_02$F1_CI_half <- round(tab_02$F1_CI_half, 4)
tab_02 <- tab_02[order(tab_02$dataset, tab_02$method,
                       tab_02$truths_effective), , drop = FALSE]
write_csv_atomic(tab_02, file.path(table_dir, "tab_02_operating_curve.csv"))

# ── 3. tab_03_k_sweep_saturation ─────────────────────────────────────────────

log_msg("Building tab_03: K-sweep saturation...")
ks_df <- subset(sAll,
                site_set == "full" &
                grepl("^ks_cp3o_K[0-9]+$", as.character(method)))
ks_df$K <- as.integer(sub("ks_cp3o_K", "", as.character(ks_df$method)))

tab_03_mean <- aggregate(
  list(F1_mean = ks_df$F1),
  list(dataset = ks_df$dataset_label, tier = ks_df$tier, K = ks_df$K),
  mean, na.rm = TRUE
)
tab_03_ci <- aggregate(
  list(F1_CI_half = ks_df$F1),
  list(dataset = ks_df$dataset_label, tier = ks_df$tier, K = ks_df$K),
  ci_half
)
tab_03_n <- aggregate(
  list(n_reps = ks_df$F1),
  list(dataset = ks_df$dataset_label, tier = ks_df$tier, K = ks_df$K),
  function(x) sum(!is.na(x))
)
tab_03 <- Reduce(function(a, b) merge(a, b,
                                      by = c("dataset", "tier", "K"),
                                      all = TRUE),
                 list(tab_03_mean, tab_03_ci, tab_03_n))
tab_03$F1_mean    <- round(tab_03$F1_mean,    4)
tab_03$F1_CI_half <- round(tab_03$F1_CI_half, 4)
tab_03 <- tab_03[order(tab_03$dataset, tab_03$tier, tab_03$K),
                 , drop = FALSE]
write_csv_atomic(tab_03, file.path(table_dir, "tab_03_k_sweep_saturation.csv"))

# ── 4. tab_04_full_vs_reduced_delta ──────────────────────────────────────────
#
# Per (dataset, tier, method): paired test of full vs reduced F1, where the
# pairing is per cell_id (same window, same method, two site sets).

log_msg("Building tab_04: full vs reduced paired Δ-F1...")

paired_block <- function(df_block) {
  full <- df_block[df_block$site_set == "full",
                   c("cell_id", "F1")]
  red  <- df_block[df_block$site_set == "reduced",
                   c("cell_id", "F1")]
  m <- merge(full, red, by = "cell_id", suffixes = c("_full", "_reduced"))
  m <- m[!is.na(m$F1_full) & !is.na(m$F1_reduced), , drop = FALSE]
  if (nrow(m) < 5L)
    return(c(F1_full_mean = NA_real_, F1_reduced_mean = NA_real_,
             delta_mean = NA_real_, wilcoxon_p = NA_real_,
             n_pairs = nrow(m)))
  delta <- m$F1_full - m$F1_reduced
  wp <- if (length(unique(delta)) > 1L) {
    suppressWarnings(
      wilcox.test(m$F1_full, m$F1_reduced, paired = TRUE)$p.value
    )
  } else NA_real_
  c(F1_full_mean    = round(mean(m$F1_full),    4),
    F1_reduced_mean = round(mean(m$F1_reduced), 4),
    delta_mean      = round(mean(delta),        4),
    wilcoxon_p      = if (is.na(wp)) NA_real_ else signif(wp, 4),
    n_pairs         = nrow(m))
}

groups <- unique(sAll[, c("dataset_label", "tier", "method")])
tab_04 <- do.call(rbind, lapply(seq_len(nrow(groups)), function(i) {
  g       <- groups[i, , drop = FALSE]
  block   <- sAll[sAll$dataset_label == g$dataset_label &
                  sAll$tier          == g$tier &
                  sAll$method        == g$method, , drop = FALSE]
  stats   <- paired_block(block)
  cbind(g, as.data.frame(t(stats)))
}))
tab_04 <- tab_04[order(tab_04$dataset_label, tab_04$tier, tab_04$method),
                 , drop = FALSE]
names(tab_04)[1] <- "dataset"
write_csv_atomic(tab_04,
                 file.path(table_dir, "tab_04_full_vs_reduced_delta.csv"))

log_msg("==============================================================")
log_msg(sprintf("All tables saved to: %s", table_dir))
log_msg("==============================================================")
