################################################################################
## helpers_windows.R — Changepoint_Detection_Study
##
## Stratified random sub-window sampling on first-of-month dates. All windows
## start and end on month boundaries; durations are even integers (months) so
## bin count = duration / bin_months is integer.
################################################################################

# ── 1. Date arithmetic (timezone-independent, first-of-month semantics) ─────
#
# We decompose Date inputs via format() and recompose via sprintf(). This
# avoids POSIXlt's timezone-dependent month rollover and is bit-identical
# across systems regardless of LC_TIME or TZ environment variables.

add_months <- function(date, n_months) {
  if (length(date) == 0L) return(date)
  y     <- as.integer(format(date, "%Y"))
  m0    <- as.integer(format(date, "%m")) - 1L     # 0-indexed month
  total <- (y * 12L + m0) + as.integer(n_months)
  new_y <- total %/% 12L
  new_m <- (total %% 12L) + 1L                     # back to 1-indexed
  as.Date(sprintf("%04d-%02d-01", new_y, new_m))
}

# Whole-month difference between two first-of-month dates: date_to − date_from.
diff_months <- function(date_to, date_from) {
  y_to   <- as.integer(format(date_to,   "%Y"))
  m_to   <- as.integer(format(date_to,   "%m"))
  y_from <- as.integer(format(date_from, "%Y"))
  m_from <- as.integer(format(date_from, "%m"))
  12L * (y_to - y_from) + (m_to - m_from)
}

# ── 2. Bin endpoints for a window ────────────────────────────────────────────

bin_endpoints <- function(window_start, window_end, bin_months = 2L) {
  if (window_end <= window_start)
    stop("`window_end` must be strictly after `window_start`.", call. = FALSE)
  total_months <- diff_months(window_end, window_start)
  if (total_months %% bin_months != 0L)
    stop("Window length (", total_months,
         " months) is not divisible by bin_months (",
         bin_months, ").", call. = FALSE)
  n_bins <- total_months %/% bin_months
  add_months(window_start, seq.int(0L, n_bins) * bin_months)
}

# ── 3. Per-bin sequence counts ───────────────────────────────────────────────
#
# Assumes `dates` is sorted ascending (preprocessing guarantees this).
# findInterval(d, endpoints, left.open = FALSE) returns k iff
# endpoints[k] ≤ d < endpoints[k+1]. Dates before window_start map to 0 and
# dates ≥ window_end map to n_bins+1; tabulate(idx, nbins = n_bins) drops both.

count_sequences_per_bin <- function(dates, window_start, window_end,
                                    bin_months = 2L) {
  endpoints <- bin_endpoints(window_start, window_end, bin_months)
  n_bins    <- length(endpoints) - 1L
  if (length(dates) == 0L) return(integer(n_bins))
  idx    <- findInterval(dates, endpoints, left.open = FALSE)
  tabulate(idx, nbins = n_bins)
}

# Row-index ranges for each bin (used to slice the feature matrix).
# Returns a list of length n_bins; element k is the integer vector of row
# indices whose Date falls in [endpoints[k], endpoints[k+1]).

bin_row_ranges <- function(dates, window_start, window_end, bin_months = 2L) {
  endpoints <- bin_endpoints(window_start, window_end, bin_months)
  n_bins    <- length(endpoints) - 1L
  bin_idx   <- findInterval(dates, endpoints, left.open = FALSE)
  # bin_idx values:
  #   0           → date < endpoints[1]
  #   1..n_bins   → date in [endpoints[k], endpoints[k+1])
  #   n_bins + 1  → date >= endpoints[n_bins+1]
  lapply(seq_len(n_bins), function(k) which(bin_idx == k))
}

# ── 4. Ground-truth filtering and bin mapping ────────────────────────────────

truths_in_window <- function(truth_df, window_start, window_end) {
  keep <- truth_df$Date_First_Detected_US >= window_start &
          truth_df$Date_First_Detected_US <  window_end
  truth_df[keep, , drop = FALSE]
}

# Map truth dates → unique Hellinger row indices (k - 1 for bin T_k).
# Truths in T1 (bin 1) are dropped — they have no Hellinger column.
# Returns a list with:
#   hellinger_idx : sorted unique integer vector of Hellinger row indices
#   dropped_T1    : count of truth EVENTS falling in T1
#   events_raw    : total truth events in the window (before T1 drop)
#   bins_unique   : length of hellinger_idx (= n_truths for matching)

map_truth_to_hellinger_idx <- function(truth_dates, window_start,
                                       bin_months = 2L) {
  if (length(truth_dates) == 0L)
    return(list(hellinger_idx = integer(0L),
                dropped_T1    = 0L,
                events_raw    = 0L,
                bins_unique   = 0L))

  months_from_start <- diff_months(truth_dates, window_start)
  bin_idx <- (months_from_start %/% bin_months) + 1L  # 1-indexed
  drop_mask <- bin_idx == 1L
  kept_hellinger <- sort(unique(bin_idx[!drop_mask] - 1L))

  list(
    hellinger_idx = kept_hellinger,
    dropped_T1    = sum(drop_mask),
    events_raw    = length(truth_dates),
    bins_unique   = length(kept_hellinger)
  )
}

# ── 5. Uniform candidate draw ────────────────────────────────────────────────

draw_one_window <- function(tier_lengths, data_start, max_window_end) {

  feasible <- tier_lengths[
    tier_lengths <= diff_months(max_window_end, data_start)
  ]
  if (length(feasible) == 0L) return(NULL)

  length_months <- sample(feasible, 1L)
  max_start     <- add_months(max_window_end, -length_months)
  n_options     <- diff_months(max_start, data_start) + 1L
  if (n_options <= 0L) return(NULL)

  start_offset <- sample(seq_len(n_options), 1L) - 1L
  window_start <- add_months(data_start, start_offset)
  window_end   <- add_months(window_start, length_months)

  list(
    window_start  = window_start,
    window_end    = window_end,
    length_months = length_months,
    n_bins        = length_months %/% 2L
  )
}

# ── 6. Acceptance criteria ───────────────────────────────────────────────────
#
# Accept iff:
#   (a) every bin contains ≥ MIN_SEQUENCES_PER_BIN sequences;
#   (b) the window contains ≥ MIN_TRUTHS_PER_WINDOW UNIQUE non-T1 truth bins.

accept_window <- function(candidate, dates, truth_df, config) {

  # Bin-level sequence counts.
  bin_counts <- count_sequences_per_bin(
    dates, candidate$window_start, candidate$window_end,
    bin_months = config$BIN_MONTHS
  )
  min_count <- min(bin_counts)
  if (min_count < config$MIN_SEQUENCES_PER_BIN)
    return(list(accepted = FALSE,
                reason   = sprintf("min bin count %d < %d",
                                   min_count, config$MIN_SEQUENCES_PER_BIN)))

  # Truth bins after T1 drop and deduplication.
  truths <- truths_in_window(truth_df, candidate$window_start,
                             candidate$window_end)
  mapping <- map_truth_to_hellinger_idx(
    truths$Date_First_Detected_US,
    candidate$window_start,
    bin_months = config$BIN_MONTHS
  )
  if (mapping$bins_unique < config$MIN_TRUTHS_PER_WINDOW)
    return(list(accepted = FALSE,
                reason   = sprintf("unique truth bins %d < %d (T1 drop = %d)",
                                   mapping$bins_unique,
                                   config$MIN_TRUTHS_PER_WINDOW,
                                   mapping$dropped_T1)))

  # Subset truths to those NOT in T1 (kept events) for diagnostic columns.
  months_to_truth <- diff_months(truths$Date_First_Detected_US,
                                 candidate$window_start)
  kept_event_mask <- (months_to_truth %/% config$BIN_MONTHS) + 1L != 1L
  kept_truths     <- truths[kept_event_mask, , drop = FALSE]

  list(
    accepted             = TRUE,
    reason               = NA_character_,
    window_start         = candidate$window_start,
    window_end           = candidate$window_end,
    length_months        = candidate$length_months,
    n_bins               = candidate$n_bins,
    bin_counts           = bin_counts,
    truths_raw           = mapping$events_raw,
    truths_dropped_T1    = mapping$dropped_T1,
    truths_effective     = mapping$bins_unique,   # UNIQUE non-T1 bins
    truth_dates          = kept_truths$Date_First_Detected_US,
    truth_labels         = kept_truths$WHO_Label,
    truth_hellinger_idx  = mapping$hellinger_idx
  )
}

# ── 7. Stratified draw loop (one tier) ───────────────────────────────────────

draw_tier <- function(tier_name, tier_lengths, data_start, max_window_end,
                      dates, truth_df, n_target, seed, config,
                      max_attempts = 100000L) {

  set.seed(seed)
  accepted <- vector("list", n_target)
  n_kept   <- 0L
  attempts <- 0L
  reasons  <- character(0L)

  while (n_kept < n_target && attempts < max_attempts) {
    attempts <- attempts + 1L
    cand <- draw_one_window(tier_lengths, data_start, max_window_end)
    if (is.null(cand)) {
      reasons <- c(reasons, "no feasible length")
      next
    }
    res <- accept_window(cand, dates, truth_df, config)
    if (!isTRUE(res$accepted)) {
      reasons <- c(reasons, res$reason)
      next
    }
    n_kept <- n_kept + 1L
    res$tier        <- tier_name
    res$run_in_tier <- n_kept
    accepted[[n_kept]] <- res
  }

  if (n_kept < n_target)
    warning(sprintf(
      "Tier %s: only %d/%d windows accepted after %d attempts.",
      tier_name, n_kept, n_target, attempts), call. = FALSE)

  list(
    accepted  = accepted[seq_len(n_kept)],
    n_target  = n_target,
    n_kept    = n_kept,
    attempts  = attempts,
    reasons   = reasons
  )
}
