################################################################################
## metrics.R — Changepoint_Detection_Study
##
## Truth-first greedy matching with consumption, and the four metrics:
## precision, recall, F1, temporal localisation error (TLE).
################################################################################

# ── 1. Truth-first greedy matching with consumption ──────────────────────────
#
# truth_idx and detected_idx are integer vectors of Hellinger row indices.
# truth_idx is deduplicated and sorted at function entry (in case the caller
# passes duplicates from co-located events). detected_idx is left as-is.
#
# For each unique truth t (in temporal order):
#   - Find the nearest available detection d with |d − t| ≤ tolerance.
#   - Tie-break in favour of the smaller d (earlier in time).
#   - Mark (t, d) as a TP and remove d from the available pool.
#   - If no such d exists, t is a FN.
# Remaining detections are FPs.
#
# Returns list(tp_pairs, fn, fp).

match_truth_to_detections <- function(truth_idx, detected_idx, tolerance = 1L) {

  truth_idx    <- sort(unique(as.integer(truth_idx)))
  detected_idx <- as.integer(detected_idx)

  if (length(truth_idx) == 0L)
    return(list(
      tp_pairs = data.frame(truth = integer(0L), detected = integer(0L),
                            abs_dist = integer(0L)),
      fn       = integer(0L),
      fp       = detected_idx
    ))

  available <- rep(TRUE, length(detected_idx))
  tp_truth  <- integer(0L)
  tp_detect <- integer(0L)
  tp_dist   <- integer(0L)
  fn_idx    <- integer(0L)

  for (t in truth_idx) {
    cand_pos <- which(available)
    if (length(cand_pos) == 0L) { fn_idx <- c(fn_idx, t); next }

    cand_vals <- detected_idx[cand_pos]
    dists     <- abs(cand_vals - t)
    min_d     <- min(dists)
    if (min_d > tolerance) { fn_idx <- c(fn_idx, t); next }

    # Tie-break: smallest detection index among those at min distance.
    tied <- which(dists == min_d)
    best_local  <- tied[which.min(cand_vals[tied])]
    best_global <- cand_pos[best_local]

    tp_truth  <- c(tp_truth,  t)
    tp_detect <- c(tp_detect, detected_idx[best_global])
    tp_dist   <- c(tp_dist,   as.integer(min_d))
    available[best_global] <- FALSE
  }

  list(
    tp_pairs = data.frame(truth = tp_truth, detected = tp_detect,
                          abs_dist = tp_dist),
    fn       = fn_idx,
    fp       = detected_idx[available]
  )
}

# ── 2. P / R / F1 / TLE from a match result ──────────────────────────────────
#
# Degenerate-case conventions (Truong, Oudre & Vayatis 2020):
#   zero detections, ≥ 1 truth      → P = 1, R = 0, F1 = 0, TLE = NA
#   ≥ 1 detection, zero truths      → P = 0, R = NA, F1 = NA, TLE = NA
#   zero detections, zero truths    → all NA (vacuous; should not occur)

compute_metrics <- function(match_result, n_truths, n_detected) {

  TP <- nrow(match_result$tp_pairs)
  FP <- length(match_result$fp)
  FN <- length(match_result$fn)

  if (TP + FP != n_detected)
    stop("Internal metric mismatch: TP + FP != n_detected.", call. = FALSE)
  if (TP + FN != n_truths)
    stop("Internal metric mismatch: TP + FN != n_truths.", call. = FALSE)

  if (n_truths == 0L && n_detected == 0L)
    return(list(P = NA_real_, R = NA_real_, F1 = NA_real_,
                TLE = NA_real_, TP = 0L, FP = 0L, FN = 0L))
  if (n_detected == 0L)
    return(list(P = 1.0, R = 0.0, F1 = 0.0,
                TLE = NA_real_, TP = TP, FP = FP, FN = FN))
  if (n_truths == 0L)
    return(list(P = 0.0, R = NA_real_, F1 = NA_real_,
                TLE = NA_real_, TP = TP, FP = FP, FN = FN))

  P  <- TP / (TP + FP)
  R  <- TP / (TP + FN)
  F1 <- if (P + R == 0) 0 else 2 * P * R / (P + R)
  TLE <- if (TP == 0L) NA_real_ else mean(match_result$tp_pairs$abs_dist)

  list(P = P, R = R, F1 = F1, TLE = TLE,
       TP = TP, FP = FP, FN = FN)
}

# ── 3. Apply matching + metrics to one detection set ─────────────────────────

evaluate_one <- function(truth_idx, detected_idx, tolerance = 1L) {
  m <- match_truth_to_detections(truth_idx, detected_idx, tolerance)
  # n_truths is the count of UNIQUE truths after dedup (same as length of
  # match_result inputs effectively); n_detected is detected_idx length.
  n_truths_unique <- length(unique(as.integer(truth_idx)))
  compute_metrics(m, n_truths = n_truths_unique,
                  n_detected = length(detected_idx))
}

# ── 4. Truth-shift sensitivity sweep ─────────────────────────────────────────
#
# For each shift s in `shifts` (in bins), shift truth_idx by s and re-match.
# Shifted truths outside [1, n_hellinger_rows] are dropped from that shift's
# evaluation. Detection set is unchanged across shifts.

compute_truth_shift_metrics <- function(truth_idx, detected_idx,
                                        n_hellinger_rows,
                                        shifts = c(-2L, -1L, 0L, 1L, 2L),
                                        tolerance = 1L) {
  out <- vector("list", length(shifts))
  names(out) <- sprintf("shift_%+d", shifts)
  for (k in seq_along(shifts)) {
    s <- shifts[k]
    shifted <- as.integer(truth_idx) + as.integer(s)
    shifted <- shifted[shifted >= 1L & shifted <= n_hellinger_rows]
    out[[k]] <- list(
      shift_bins      = as.integer(s),
      truth_idx_kept  = shifted,
      n_truths_kept   = length(unique(shifted)),
      metrics         = evaluate_one(shifted, detected_idx, tolerance)
    )
  }
  out
}

# ── 5. Evaluate every method's detection set ─────────────────────────────────

evaluate_methods <- function(method_results, truth_idx, n_hellinger_rows,
                             tolerance = 1L, truth_shifts = NULL) {
  out <- method_results
  for (key in names(method_results)) {
    detected <- method_results[[key]]$detected_cps
    out[[key]]$metrics_primary <- evaluate_one(truth_idx, detected, tolerance)
    if (!is.null(truth_shifts) && length(truth_shifts) > 0L) {
      out[[key]]$metrics_truth_shift <- compute_truth_shift_metrics(
        truth_idx          = truth_idx,
        detected_idx       = detected,
        n_hellinger_rows   = n_hellinger_rows,
        shifts             = truth_shifts,
        tolerance          = tolerance
      )
    }
  }
  out
}
