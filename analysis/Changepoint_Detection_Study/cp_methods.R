################################################################################
## cp_methods.R — Changepoint_Detection_Study
##
## Four change-point detection methods, harmonised to return a sorted integer
## vector of interior CP indices in the Hellinger time series.
##
## Methods:
##   ECP-1     : ecp::ks.cp3o          (K-required, dynamic + K-sweep ∈ {2,5,10})
##   ECP-2     : ecp::e.agglo          (K-free hierarchical agglomerative)
##   HDcp-BS   : HDcpDetect::binary.segmentation         (via wrapper)
##   HDcp-WBS  : HDcpDetect::wild.binary.segmentation    (via wrapper, M = 100)
##
## tryCatch policy: catch ERRORS only. Warnings from underlying methods are
## suppressed (via suppressWarnings) but do not abort the call — some methods
## emit benign warnings (e.g., K capped) while still returning valid output.
##
## Input convention: `dat_t` is a T × d matrix (rows = time points, cols =
## sites). For a Hellinger matrix `H` of shape sites × time_steps, callers
## pass `t(H)`.
################################################################################

# ── 1. Interior-CP filter ────────────────────────────────────────────────────
#
# A CP at row k means "row k starts a new segment relative to row k − 1".
# Valid interior indices: 2..n_rows. Index 1 is the trivial start (every
# series begins a "segment" at row 1) and is stripped. Anything outside
# [2, n_rows] is also stripped (out-of-range boundary). Duplicates removed,
# result sorted.

.strip_boundaries <- function(idx, n_rows) {
  if (length(idx) == 0L) return(integer(0L))
  idx <- as.integer(idx)
  idx <- idx[!is.na(idx)]
  idx <- idx[idx >= 2L & idx <= n_rows]
  sort(unique(idx))
}

# ── 2. ks.cp3o (energy distance, dynamic programming) ────────────────────────

run_ks_cp3o <- function(dat_t, K, minsize) {

  n  <- nrow(dat_t)
  t0 <- Sys.time()

  result <- list(
    method        = "ks_cp3o",
    K             = as.integer(K),
    detected_cps  = integer(0L),
    status        = "ok",
    error_msg     = NA_character_,
    walltime_s    = NA_real_,
    extra         = list()
  )

  fit <- tryCatch(
    suppressWarnings(
      ecp::ks.cp3o(Z = dat_t, K = as.integer(K),
                   minsize = as.integer(minsize), verbose = FALSE)
    ),
    error = function(e) e
  )

  result$walltime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(fit, "error")) {
    result$status    <- "failed"
    result$error_msg <- conditionMessage(fit)
    return(result)
  }

  result$detected_cps <- .strip_boundaries(fit$estimates, n)
  result$extra        <- list(n_cps_internal = length(fit$estimates))
  result
}

# ── 3. e.agglo (energy distance, hierarchical agglomerative, K-free) ─────────

run_e_agglo <- function(dat_t, alpha = 1, penalty = function(cps) 0) {

  n  <- nrow(dat_t)
  t0 <- Sys.time()

  result <- list(
    method        = "e_agglo",
    detected_cps  = integer(0L),
    status        = "ok",
    error_msg     = NA_character_,
    walltime_s    = NA_real_,
    extra         = list()
  )

  fit <- tryCatch(
    suppressWarnings(
      ecp::e.agglo(X = dat_t, member = seq_len(n),
                   alpha = alpha, penalty = penalty)
    ),
    error = function(e) e
  )

  result$walltime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(fit, "error")) {
    result$status    <- "failed"
    result$error_msg <- conditionMessage(fit)
    return(result)
  }

  result$detected_cps <- .strip_boundaries(fit$estimates, n)
  result$extra        <- list(n_cps_internal = length(fit$estimates))
  result
}

# ── 4. HDcpDetect: binary segmentation (preserves p-values) ──────────────────

run_hdcp_binseg <- function(dat_t) {

  n  <- nrow(dat_t)
  t0 <- Sys.time()

  result <- list(
    method        = "hdcp_binseg",
    detected_cps  = integer(0L),
    status        = "ok",
    error_msg     = NA_character_,
    walltime_s    = NA_real_,
    pvalues       = numeric(0L),
    extra         = list()
  )

  fit <- tryCatch(
    suppressWarnings(detect_changepoints_hdcp(
      data_matrix    = dat_t,
      min_window     = 1L,
      max_window     = n,
      n_timesteps    = 0L,
      rolling_window = FALSE,
      wild           = FALSE
    )),
    error = function(e) e
  )

  result$walltime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(fit, "error")) {
    result$status    <- "failed"
    result$error_msg <- conditionMessage(fit)
    return(result)
  }

  hdcp_out <- fit$HDcp_list[[1L]]
  if (is.null(hdcp_out)) {
    result$status    <- "failed"
    result$error_msg <- "HDcp_list[[1]] is NULL"
    return(result)
  }

  # binary.segmentation returns a matrix with cols FoundList and pvalues.
  # Empty result is a 0-row matrix.
  if (is.matrix(hdcp_out) && nrow(hdcp_out) > 0L) {
    idx <- as.integer(hdcp_out[, "FoundList"])
    pv  <- as.numeric(hdcp_out[, "pvalues"])
  } else {
    idx <- integer(0L)
    pv  <- numeric(0L)
  }

  # CONVENTION ALIGNMENT: HDcpDetect returns the LAST row of the OLD segment;
  # ECP (ks.cp3o, e.agglo) and our truth-to-Hellinger mapping both use the
  # FIRST row of the NEW segment. Add +1 to align. Verified empirically on
  # the canonical 20-row baseline/variant example from each package's help
  # page: ECP returns 11, HDcp binary.segmentation returns 10, truth=11.
  idx_aligned <- idx + 1L
  kept        <- .strip_boundaries(idx_aligned, n)
  result$detected_cps <- kept
  if (length(kept) > 0L)
    result$pvalues <- pv[match(kept, idx_aligned)]
  result$extra <- list(n_cps_internal = length(idx))
  result
}

# ── 5. HDcpDetect: wild binary segmentation ──────────────────────────────────

run_hdcp_wbs <- function(dat_t, M = 100L) {

  n  <- nrow(dat_t)
  t0 <- Sys.time()

  result <- list(
    method        = "hdcp_wbs",
    M             = as.integer(M),
    detected_cps  = integer(0L),
    status        = "ok",
    error_msg     = NA_character_,
    walltime_s    = NA_real_,
    extra         = list()
  )

  fit <- tryCatch(
    suppressWarnings(detect_changepoints_hdcp(
      data_matrix    = dat_t,
      min_window     = 1L,
      max_window     = n,
      n_timesteps    = 0L,
      rolling_window = FALSE,
      wild           = TRUE,
      M              = as.integer(M)
    )),
    error = function(e) e
  )

  result$walltime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(fit, "error")) {
    result$status    <- "failed"
    result$error_msg <- conditionMessage(fit)
    return(result)
  }

  hdcp_out <- fit$HDcp_list[[1L]]
  if (is.null(hdcp_out)) {
    result$status    <- "failed"
    result$error_msg <- "HDcp_list[[1]] is NULL"
    return(result)
  }

  # wild.binary.segmentation returns a plain numeric vector of CP indices
  # (no p-values). Defensive: handle matrix form just in case.
  idx <- if (is.numeric(hdcp_out) || is.integer(hdcp_out)) {
    as.integer(hdcp_out)
  } else if (is.matrix(hdcp_out)) {
    as.integer(hdcp_out[, 1L])
  } else {
    integer(0L)
  }

  # CONVENTION ALIGNMENT: see run_hdcp_binseg for full explanation.
  # HDcpDetect → ECP convention via +1 shift before boundary stripping.
  result$detected_cps <- .strip_boundaries(idx + 1L, n)
  result$extra        <- list(n_cps_internal = length(idx))
  result
}

# ── 6. Dispatch all four methods (plus K-sweep for ks.cp3o) ──────────────────

run_all_methods <- function(dat_t, config) {

  n <- nrow(dat_t)
  K_dynamic <- max(1L, as.integer(n - 2L))

  results <- list()

  results$ks_cp3o_dynamic <- run_ks_cp3o(
    dat_t, K = K_dynamic, minsize = config$KS_CP3O_MINSIZE
  )

  for (K in config$KS_CP3O_K_SWEEP) {
    key <- sprintf("ks_cp3o_K%d", K)
    results[[key]] <- run_ks_cp3o(
      dat_t, K = K, minsize = config$KS_CP3O_MINSIZE
    )
  }

  results$e_agglo     <- run_e_agglo(dat_t, alpha = config$E_AGGLO_ALPHA,
                                     penalty = function(cps) 0)
  results$hdcp_binseg <- run_hdcp_binseg(dat_t)
  results$hdcp_wbs    <- run_hdcp_wbs(dat_t, M = config$HDCP_M)

  results
}
