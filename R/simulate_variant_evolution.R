#' @title Simulate Viral Variant Evolution
#' @description Generates a synthetic, ground-truth-traceable time series of
#'   viral amino-acid sequences in which a reference strain is progressively
#'   challenged by emerging variants under stochastic growth and pairwise
#'   competition. Designed as a controllable signal generator for benchmarking
#'   the entropy, Hellinger, and change-point components of the
#'   \pkg{ViralEntropR} pipeline; the dynamics are deliberately simplified
#'   rather than a faithful biological model of evolutionary fitness.
#' @export
#'
#' @details
#' Simulates the trajectory of a viral population over discrete monthly time
#' steps. All stochastic components are reproducible under a fixed
#' \code{seed} and are designed to expose a known emergence pattern that
#' downstream detection methods can be evaluated against.
#'
#' \strong{Simulation Phases:}
#' \enumerate{
#'   \item \strong{Reference Phase:} A stable period in which only the
#'     reference strain exists. Reference months are resampled with
#'     replacement to assemble \code{n_ref_months} monthly batches.
#'   \item \strong{Variant Emergence:} New variants are introduced at the
#'     monthly intervals specified by \code{variant_intervals}.
#'   \item \strong{Pairwise Growth:} Each month the newest variant grows
#'     alongside its immediate predecessor. Both grow stochastically: each
#'     variant's next-month count is its previous-month count times
#'     \code{Uniform(mult*(1-variability), mult*(1+variability))}. The two
#'     are not coupled by a shared frequency constraint — only by the hard
#'     monthly quota
#'     \code{ceiling(n_sequences_total / (n_periods - n_ref_months))} that
#'     trims any overshoot. Older, non-paired variants enter the fill pool
#'     with sampling weights \code{2^(i+1)}.
#'   \item \strong{Capped "Deleterious" Sites:} A Bernoulli-random subset of
#'     each variant's mutated positions is flagged with probability
#'     \code{deleterious_rate}. At each post-emergence month the count of
#'     variant rows still carrying the flagged allele is capped at
#'     \code{n_deleterious_limit} by reverting overflow sequences to the
#'     reference allele \emph{at the flagged site only}. This is a
#'     per-site cap, not a genotype-level fitness penalty.
#'   \item \strong{Random Substitution Events:} Each month independently
#'     fires a Bernoulli coin flip with probability
#'     \code{prob_deletion_event}. If it fires, the last
#'     \code{n_rows_to_delete} rows of that month's batch receive a single
#'     non-reference amino-acid substitution at a randomly chosen site
#'     (sites 1 and 2 excluded), and those rows are flagged with
#'     \code{Delet = "Yes"}. Parameter names retain "deletion" / "delete"
#'     for backward compatibility; the operation itself is a substitution.
#' }
#'
#' @param ref_sequences Either a single character string giving the reference
#'   amino-acid sequence, or a data.frame with a \code{Date} column and site
#'   columns named with integer position labels (\code{"1"}, \code{"2"},
#'   ..., \code{as.character(L)}; one column per site). The sequence length
#'   \code{L} must be at least 3, since variant mutations and random
#'   substitution events are sampled from sites \eqn{>2} (sites 1 and 2 are
#'   reserved for reference variability when \code{ref_variability = TRUE}).
#' @param n_ref_months Integer. Duration of the initial reference phase
#'   (months). Default \code{3}.
#' @param start_date Character or Date. Simulation start.
#' @param end_date Character or Date. Simulation end (inclusive month).
#' @param variants_config Integer vector. Number of mutations per variant,
#'   e.g., \code{c(2, 6)}.
#' @param variant_intervals Integer vector. Months between consecutive
#'   variant emergences (length \code{length(variants_config) - 1}).
#' @param n_new_mutations Integer. Number of de-novo sequence \emph{copies}
#'   introduced at a variant's first appearance. Despite the name, this
#'   controls the count of seeded sequences, not the number of mutations
#'   carried by the variant (which is set per-variant by
#'   \code{variants_config}). Default \code{1}.
#' @param mutation_rate Numeric scalar or vector. Monthly growth multiplier
#'   per variant. Default \code{2}.
#' @param mutation_rate_variability Numeric in \code{[0, 1)}. Fractional
#'   spread around \code{mutation_rate}: growth is drawn from
#'   \code{Uniform(mult*(1-v), mult*(1+v))}. Set to \code{0} for
#'   deterministic growth. Default \code{0.25}.
#' @param deleterious_rate Numeric in \code{[0, 1]}. Per-mutation probability
#'   of being flagged as "deleterious". Default \code{0}.
#' @param n_deleterious_limit Integer. Per-site cap on the number of variant
#'   rows allowed to retain a flagged ("deleterious") allele in a given
#'   month; overflow is reverted to the reference allele at the flagged
#'   site only. See \strong{Capped "Deleterious" Sites} in Details.
#'   Default \code{1}.
#' @param n_sequences_total Integer. Total sequences generated over the full
#'   mutation phase (spread evenly per month). Default \code{140}.
#' @param ref_variability Logical. If \code{TRUE}, introduces low-level
#'   variability at sites 1 and 2 of the reference pool. Requires
#'   \code{n_ref_sequences >= 6 * n_ref_months} for the variability rows
#'   to be correctly captured by the multi-variant fill pool; smaller
#'   values trigger a warning. Default \code{FALSE}.
#' @param n_ref_sequences Integer. When \code{ref_sequences} is a single
#'   string, the number of reference rows to create. Default \code{100}.
#' @param prob_deletion_event Numeric in \code{[0, 1]}. Per-month probability
#'   of triggering a substitution event affecting \code{n_rows_to_delete}
#'   rows. (Parameter name retained for backward compatibility; the
#'   operation is substitution, not literal deletion. See Details.)
#'   Default \code{0}.
#' @param n_rows_to_delete Integer. Number of rows affected when a
#'   substitution event fires (see \code{prob_deletion_event}). Default
#'   \code{0}.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#'   Default \code{NULL}.
#'
#' @return An object of class \code{"viralSim"}, a named list:
#' \item{Simulation_Output}{Data frame of all sequences, one column per site
#'   (named \code{"1"}, ..., \code{as.character(L)}) plus \code{Variant}
#'   (per-row strain label), \code{Phase} (\code{"R"} or \code{"RV"} for
#'   reference rows depending on whether the row matches the wild-type
#'   sequence; \code{"M1"}, \code{"M2"}, ... in mutation-phase months —
#'   note that \code{"Mn"} is a \emph{month label} and applies to every
#'   row in that month regardless of strain), \code{Date}, \code{Period}
#'   (1-indexed month), and \code{Delet} (\code{"Yes"} for rows affected
#'   by substitution events, \code{"No"} otherwise).}
#' \item{Variant_Details}{List of per-variant metadata. Each element has
#'   \code{em} (emergence-month index), \code{pos} (mutated site indices),
#'   \code{flags} (logical vector of "deleterious" flags per mutated site),
#'   \code{vseq} (full mutated reference sequence as character vector),
#'   \code{mult} (growth multiplier), \code{last} (count in the most recent
#'   month), \code{cum} (cumulative count), and \code{repl} (logical flag
#'   set when the variant has been displaced by its successor).}
#' \item{Simulation_Dates}{Date vector of all simulated months.}
#' \item{Baseline_Ref_Sequence}{Character. The wild-type reference string.}
#' \item{Delet_Records}{Named list of substitution events keyed by date
#'   (\code{"YYYY-MM-DD"}). Each entry holds \code{period}, \code{date},
#'   \code{site}, \code{old_aa}, \code{new_aa}, and \code{rows} (absolute
#'   row indices in \code{Simulation_Output} that were modified).}
#' \item{Pool}{Fill-pool data frame from the final multi-variant competition
#'   period (includes \code{._weight} column); \code{NULL} if at most one
#'   variant was ever active.}
#'
#' @importFrom stats rbinom runif
#' @importFrom utils tail
#'
#' @seealso \code{\link{partition_time_windows}},
#'   \code{\link{calculate_hellinger_matrix}},
#'   \code{\link{detect_changepoints_ecp}},
#'   \code{\link{detect_changepoints_hdcp}}
#'
#' @examples
#' ref_seq <- "MKTIIALSYIFCLVFADYKDDDDK"
#'
#' sim <- simulate_variant_evolution(
#'   ref_sequences   = ref_seq,
#'   n_ref_months    = 3,
#'   start_date      = "2021-01-01",
#'   end_date        = "2021-12-01",
#'   variants_config = c(3, 5),
#'   variant_intervals = c(4),
#'   n_sequences_total = 50,
#'   mutation_rate   = 1.5,
#'   seed            = 123
#' )
#' head(sim$Simulation_Output)
simulate_variant_evolution <- function(
    ref_sequences,
    n_ref_months             = 3L,
    start_date,
    end_date,
    variants_config,
    variant_intervals,
    n_new_mutations          = 1L,
    mutation_rate            = 2,
    mutation_rate_variability = 0.25,
    deleterious_rate         = 0,
    n_deleterious_limit      = 1L,
    n_sequences_total        = 140L,
    ref_variability          = FALSE,
    n_ref_sequences          = 100L,
    prob_deletion_event      = 0,
    n_rows_to_delete         = 0L,
    seed                     = NULL
) {
  
  # ── 0. Input validation ─────────────────────────────────────────────────
  
  # ref_sequences: single character string OR data.frame with Date column
  if (!((is.character(ref_sequences) && length(ref_sequences) == 1L) ||
        (is.data.frame(ref_sequences) && "Date" %in% names(ref_sequences))))
    stop("`ref_sequences` must be either a single character string or a ",
         "data.frame with a `Date` column.", call. = FALSE)
  
  if (is.character(ref_sequences) && nchar(ref_sequences) < 3L)
    stop("`ref_sequences` must encode a sequence of length >= 3 ",
         "(sites 1 and 2 are reserved for reference variability).",
         call. = FALSE)
  if (is.data.frame(ref_sequences) &&
      sum(grepl("^[0-9]+$", names(ref_sequences))) < 3L)
    stop("`ref_sequences` data.frame must contain at least 3 site columns ",
         "named with integer position labels (\"1\", \"2\", ...).",
         call. = FALSE)
  
  # Date arguments
  start_date <- tryCatch(as.Date(start_date),
                         error = function(e)
                           stop("`start_date` is not a valid date.",
                                call. = FALSE))
  end_date   <- tryCatch(as.Date(end_date),
                         error = function(e)
                           stop("`end_date` is not a valid date.",
                                call. = FALSE))
  if (anyNA(c(start_date, end_date)))
    stop("`start_date` and `end_date` must both be valid dates.",
         call. = FALSE)
  if (start_date > end_date)
    stop("`start_date` must be on or before `end_date`.", call. = FALSE)
  
  # Integer-count arguments (non-negative scalars)
  is_count <- function(x) is.numeric(x) && length(x) == 1L &&
    !is.na(x) && x >= 0 && x == as.integer(x)
  if (!is_count(n_ref_months) || n_ref_months < 1L)
    stop("`n_ref_months` must be a positive integer.", call. = FALSE)
  for (nm in c("n_new_mutations", "n_deleterious_limit",
               "n_sequences_total", "n_ref_sequences", "n_rows_to_delete")) {
    if (!is_count(get(nm)))
      stop(sprintf("`%s` must be a non-negative integer.", nm),
           call. = FALSE)
  }
  
  # Timeline must extend beyond the reference phase
  n_months_total <- length(seq(start_date, end_date, by = "month"))
  if (n_months_total <= n_ref_months)
    stop("The simulation spans only ", n_months_total, " month(s) but ",
         "`n_ref_months` is ", n_ref_months, "; the timeline must include ",
         "at least one period beyond the reference phase.", call. = FALSE)
  
  # variants_config: positive integer vector, length >= 1
  if (!is.numeric(variants_config) || length(variants_config) < 1L ||
      anyNA(variants_config) || any(variants_config < 1) ||
      any(variants_config != as.integer(variants_config)))
    stop("`variants_config` must be a positive-integer vector of length >= 1.",
         call. = FALSE)
  
  # variant_intervals: positive integers when more than one variant requested
  nv_arg <- length(variants_config)
  if (nv_arg > 1L) {
    if (!is.numeric(variant_intervals) || length(variant_intervals) < 1L ||
        anyNA(variant_intervals) || any(variant_intervals < 1) ||
        any(variant_intervals != as.integer(variant_intervals)))
      stop("`variant_intervals` must be a positive-integer vector of length ",
           ">= 1 when `variants_config` has more than one variant.",
           call. = FALSE)
  }
  
  # mutation_rate: positive numeric, scalar or matched to variants_config
  if (!is.numeric(mutation_rate) || anyNA(mutation_rate) ||
      any(mutation_rate <= 0) ||
      !(length(mutation_rate) == 1L || length(mutation_rate) == nv_arg))
    stop("`mutation_rate` must be a positive numeric scalar or a vector of ",
         "length `length(variants_config)`.", call. = FALSE)
  
  # Probability-like arguments
  is_prob <- function(x, strict_upper = FALSE)
    is.numeric(x) && length(x) == 1L && !is.na(x) && x >= 0 &&
    (if (strict_upper) x < 1 else x <= 1)
  if (!is_prob(mutation_rate_variability, strict_upper = TRUE))
    stop("`mutation_rate_variability` must be a numeric scalar in [0, 1).",
         call. = FALSE)
  if (!is_prob(deleterious_rate))
    stop("`deleterious_rate` must be a numeric scalar in [0, 1].",
         call. = FALSE)
  if (!is_prob(prob_deletion_event))
    stop("`prob_deletion_event` must be a numeric scalar in [0, 1].",
         call. = FALSE)
  
  # ref_variability: single logical
  if (!is.logical(ref_variability) || length(ref_variability) != 1L ||
      is.na(ref_variability))
    stop("`ref_variability` must be a single logical value.", call. = FALSE)
  
  # Multi-variant fill-pool integrity: the hardcoded row indices c(1, 2, 5)
  # used in the fill-pool construction (Section 7c) assume each ref-phase
  # month carries at least 6 rows when `ref_variability = TRUE`, so that
  # idx 1 = wild-type, idx 2-4 = site-1 variability, idx 5-6 = site-2
  # variability. Below this density the variability rows can be pulled from
  # the wrong month chunk and the fill pool is mis-composed.
  if (ref_variability && n_ref_sequences < 6L * n_ref_months)
    warning("`ref_variability = TRUE` assumes >= 6 reference rows per month ",
            "(have ", n_ref_sequences, " across ", n_ref_months,
            " months). Site-1/site-2 variability rows may not be correctly ",
            "captured in the multi-variant fill pool.", call. = FALSE)
  
  # seed: NULL or single numeric value
  if (!is.null(seed) &&
      !(is.numeric(seed) && length(seed) == 1L && !is.na(seed)))
    stop("`seed` must be NULL or a single numeric value.", call. = FALSE)
  
  if (!is.null(seed)) set.seed(seed)
  
  aa <- c("A","R","N","D","C","E","Q","G","H","I",
          "L","K","M","F","P","S","T","W","Y","V")
  
  # ── 1. Internal helpers ──────────────────────────────────────────────────
  
  # del_cap: revert excess deleterious-mutation sequences back to reference
  del_cap <- function(mat, vidx, vdet, limit, ref_chars) {
    if (!is.matrix(mat)) {
      v <- as.vector(mat)
      mat <- matrix(v, nrow = 1L, ncol = length(v))
    }
    vd <- vdet[[vidx]]
    if (!any(vd$flags)) return(mat)
    for (p in vd$pos[vd$flags]) {
      allele <- vd$vseq[p]
      cnt    <- sum(mat[, p] == allele)
      if (cnt > limit) {
        cands <- which(mat[, p] == allele)
        mat[sample(cands, cnt - limit), p] <- ref_chars[p]
      }
    }
    mat
  }
  
  # ── 2. Build reference data frame ────────────────────────────────────────
  
  if (is.character(ref_sequences) && length(ref_sequences) == 1L) {
    base_seq <- ref_sequences
    L        <- nchar(base_seq)
    ref_df   <- as.data.frame(
      matrix(strsplit(base_seq, "")[[1]], nrow = n_ref_sequences, ncol = L,
             byrow = TRUE),
      stringsAsFactors = FALSE
    )
    colnames(ref_df) <- as.character(seq_len(L))
    ref_df$Date <- rep(
      seq(start_date, by = "month", length.out = n_ref_months),
      length.out = n_ref_sequences
    )
  } else {
    ref_df   <- ref_sequences
    seq_num  <- grep("^[0-9]+$", names(ref_df))
    base_seq <- paste(ref_df[1L, seq_num], collapse = "")
    L        <- length(seq_num)
  }
  
  seq_cols <- as.character(seq_len(L))
  
  # Add Month for reference-phase construction (kept until ref_phase is built)
  ref_df$Month <- format(as.Date(ref_df$Date), "%Y-%m")
  
  # ── 3. Reference variability (2 RNG draws if TRUE) ───────────────────────
  
  if (ref_variability) {
    sub1 <- sample(setdiff(aa, substr(base_seq, 1L, 1L)), 1L)
    sub2 <- sample(setdiff(aa, substr(base_seq, 2L, 2L)), 1L)
    for (m in unique(ref_df$Month)) {
      idx <- sort(which(ref_df$Month == m))
      if (length(idx) >= 4L) ref_df[idx[2L:4L], seq_cols[1L]] <- sub1
      if (length(idx) >= 6L) ref_df[idx[5L:6L], seq_cols[2L]] <- sub2
    }
  }
  
  # ── 4. Simulation timeline ───────────────────────────────────────────────
  
  # seq.Date inclusive of both endpoints gives the full set of monthly time
  # points. Note that lubridate::interval() %/% period() returns the
  # *distance* between two dates (e.g. 23 for Jan-2020 -> Dec-2021), which
  # is one short of the 24 monthly time points actually traversed.
  sim_dates <- seq(start_date, end_date, by = "month")
  n_periods <- length(sim_dates)
  
  template_cols <- c(seq_cols, "Variant", "Phase", "Date",
                     "Period", "Delet", "._weight")
  
  ensure_cols <- function(d) {
    miss <- setdiff(template_cols, names(d))
    d[miss] <- NA
    d[, template_cols, drop = FALSE]
  }
  
  # ── 5. Build reference phase via random month sampling (1 RNG draw) ──────
  
  choice    <- sample(unique(ref_df$Month), n_ref_months, replace = TRUE)
  ref_phase <- do.call(rbind, lapply(seq_len(n_ref_months), function(i) {
    chunk          <- ref_df[ref_df$Month == choice[i], ]
    chunk$Month    <- NULL
    chunk$Date     <- sim_dates[i]
    chunk$Period   <- i
    chunk$Variant  <- "Reference"
    chunk$Phase    <- ifelse(
      apply(chunk[, seq_cols, drop = FALSE], 1L,
            function(r) paste(r, collapse = "")) == base_seq,
      "R", "RV"
    )
    chunk$Delet    <- "No"
    ensure_cols(chunk)
  }))
  ref_phase[["._weight"]] <- 1
  
  output <- ref_phase
  
  # ── 6. Define variants (RNG: pos + AA per variant + rbinom flags) ─────────
  
  nv    <- length(variants_config)
  gaps  <- rep(variant_intervals, length.out = max(0L, nv - 1L))
  em    <- numeric(nv)
  em[1] <- n_ref_months + 1L
  if (nv > 1L) for (j in 2L:nv) em[j] <- em[j - 1L] + gaps[j - 1L]
  
  valid           <- which(em <= n_periods)
  em              <- em[valid]
  variants_config <- variants_config[valid]
  nv_valid        <- length(em)
  rates <- if (length(mutation_rate) == 1L) rep(mutation_rate, nv) else mutation_rate
  
  ref_seq_chars <- as.character(ref_phase[1L, seq_cols])
  
  variant_details <- vector("list", nv_valid)
  vrow_pool       <- vector("list", nv_valid)  # pre-built rows for fill pool
  
  for (i in seq_len(nv_valid)) {
    nm      <- variants_config[[i]]
    pos     <- sample(setdiff(seq_len(L), 1L:2L), nm)
    new_aas <- sapply(pos, function(p) sample(setdiff(aa, ref_seq_chars[p]), 1L))
    vseq    <- ref_seq_chars; vseq[pos] <- new_aas
    flags   <- stats::rbinom(nm, 1L, deleterious_rate) == 1L
    variant_details[[i]] <- list(
      em   = em[i], pos = pos, flags = flags, vseq = vseq,
      mult = rates[i], last = 0L, cum = 0L, repl = FALSE
    )
    # Pre-built pool row (used when this variant is neither cur nor cmp)
    vr                   <- as.data.frame(
      matrix(vseq, 1L, dimnames = list(NULL, seq_cols)),
      stringsAsFactors = FALSE
    )
    vr$Variant           <- paste0("Variant_", i)
    vr$Phase             <- paste0("M", i)
    vr$Date              <- as.Date(NA)
    vr$Period            <- NA_integer_
    vr$Delet             <- "No"
    vr[["._weight"]]     <- 2^(i + 1L)
    vrow_pool[[i]]       <- ensure_cols(vr)
  }
  
  n_fixed_per_month <- ceiling(n_sequences_total /
                                 (n_periods - n_ref_months))
  delet_records     <- list()
  pool              <- NULL   # updated each competition period; returned for inspection
  
  # ── 7. Mutation phase loop ───────────────────────────────────────────────
  
  for (t in (n_ref_months + 1L):n_periods) {
    
    today     <- sim_dates[t]
    row_start <- nrow(output) + 1L
    
    # Identify most-recently-emerged non-replaced variant
    active_idx <- 0L
    for (v in seq_len(nv_valid)) {
      if (t >= variant_details[[v]]$em && !variant_details[[v]]$repl)
        active_idx <- v
    }
    
    # ── 7a. Reference-only (pre-variant) ──────────────────────────────────
    if (active_idx == 0L) {
      # Dead code under standard parameters (loop starts at em[1] == n_ref_months+1)
      # but kept for correctness if parameters ever create a gap
      samp <- sample(nrow(ref_phase), n_fixed_per_month, TRUE,
                     ref_phase[["._weight"]])
      df              <- ref_phase[samp, template_cols]
      df$Date         <- today
      df$Period       <- t
      
      # ── 7b. Single variant ────────────────────────────────────────────────
    } else if (active_idx == 1L) {
      
      vd    <- variant_details[[1L]]
      n_cur <- if (t == vd$em) n_new_mutations else
        max(1L, ceiling(vd$last *
                          stats::runif(1L,
                                       vd$mult * (1 - mutation_rate_variability),
                                       vd$mult * (1 + mutation_rate_variability))))
      variant_details[[1L]]$last <- n_cur
      variant_details[[1L]]$cum  <- variant_details[[1L]]$cum + n_cur
      
      m_cur <- matrix(vd$vseq, nrow = n_cur, ncol = L, byrow = TRUE,
                      dimnames = list(NULL, seq_cols))
      left  <- n_fixed_per_month - n_cur
      
      if (left > 0L) {
        samp <- sample(nrow(ref_phase), left, TRUE, ref_phase[["._weight"]])
        mat  <- rbind(m_cur, as.matrix(ref_phase[samp, seq_cols]))
        lbl  <- c(rep("Variant_1", n_cur), ref_phase$Variant[samp])
      } else {
        mat <- m_cur
        lbl <- rep("Variant_1", n_cur)
      }
      
      idx_v      <- seq_len(min(n_cur, nrow(mat)))
      mat[idx_v, ] <- del_cap(mat[idx_v, , drop = FALSE], 1L,
                              variant_details, n_deleterious_limit,
                              ref_seq_chars)
      
      df         <- as.data.frame(mat, stringsAsFactors = FALSE)
      df$Variant <- lbl; df$Phase <- "M1"
      df$Date    <- today; df$Period <- t; df$Delet <- "No"
      
      # ── 7c. Multi-variant competition ─────────────────────────────────────
    } else {
      
      cur    <- active_idx
      cmp    <- cur - 1L
      vd_cur <- variant_details[[cur]]
      vd_cmp <- variant_details[[cmp]]
      
      # Stochastic growth for both competing variants
      n_cur <- if (t == vd_cur$em) n_new_mutations else
        max(1L, ceiling(vd_cur$last *
                          stats::runif(1L,
                                       vd_cur$mult * (1 - mutation_rate_variability),
                                       vd_cur$mult * (1 + mutation_rate_variability))))
      variant_details[[cur]]$last <- n_cur
      variant_details[[cur]]$cum  <- variant_details[[cur]]$cum + n_cur
      
      n_cmp <- max(1L, ceiling(vd_cmp$last *
                                 stats::runif(1L,
                                              vd_cmp$mult * (1 - mutation_rate_variability),
                                              vd_cmp$mult * (1 + mutation_rate_variability))))
      variant_details[[cmp]]$last <- n_cmp
      variant_details[[cmp]]$cum  <- variant_details[[cmp]]$cum + n_cmp
      
      m_cur <- matrix(vd_cur$vseq, nrow = n_cur, ncol = L, byrow = TRUE,
                      dimnames = list(NULL, seq_cols))
      m_cmp <- matrix(vd_cmp$vseq, nrow = n_cmp, ncol = L, byrow = TRUE,
                      dimnames = list(NULL, seq_cols))
      mat   <- rbind(m_cur, m_cmp)
      lbl   <- c(rep(paste0("Variant_", cur), n_cur),
                 rep(paste0("Variant_", cmp), n_cmp))
      
      # Fill remainder from reference + older (non-competing) variants
      left <- n_fixed_per_month - nrow(mat)
      if (left > 0L) {
        pool <- if (ref_variability) ref_phase[c(1L, 2L, 5L), ]
        else                 ref_phase[1L, , drop = FALSE]
        if (cur >= 3L) {
          for (j in seq_len(cur - 2L)) {
            vr          <- vrow_pool[[j]]
            vr$Date     <- today
            vr$Period   <- t
            vr$Phase    <- paste0("M", cur)
            pool        <- rbind(pool, vr)
          }
        }
        samp <- sample(nrow(pool), left, TRUE, pool[["._weight"]])
        mat  <- rbind(mat, as.matrix(pool[samp, seq_cols]))
        lbl  <- c(lbl, pool$Variant[samp])
      }
      
      # Trim to quota if cur+cmp overshoot
      if (nrow(mat) > n_fixed_per_month) {
        mat <- mat[seq_len(n_fixed_per_month), ]
        lbl <- lbl[seq_len(n_fixed_per_month)]
      }
      
      # Apply deleterious caps to cur and cmp rows
      idx_c <- seq_len(min(n_cur, nrow(mat)))
      mat[idx_c, ] <- del_cap(mat[idx_c, , drop = FALSE], cur,
                              variant_details, n_deleterious_limit,
                              ref_seq_chars)
      cs <- min(n_cur + 1L, nrow(mat))
      ce <- min(n_cur + n_cmp, nrow(mat))
      if (cs <= ce)
        mat[cs:ce, ] <- del_cap(mat[cs:ce, , drop = FALSE], cmp,
                                variant_details, n_deleterious_limit,
                                ref_seq_chars)
      
      df         <- as.data.frame(mat, stringsAsFactors = FALSE)
      df$Variant <- lbl; df$Phase <- paste0("M", cur)
      df$Date    <- today; df$Period <- t; df$Delet <- "No"
      
      # Mark predecessor as replaced once it drops out completely
      if (!any(lbl == paste0("Variant_", cmp)))
        variant_details[[cmp]]$repl <- TRUE
    }
    
    # ── 7d. Coin-flip substitution event ──────────────────────────────────
    if (stats::runif(1L) < prob_deletion_event) {
      site_idx          <- sample(setdiff(seq_len(L), 1L:2L), 1L)
      old_aa            <- df[1L, as.character(site_idx)]
      new_aa            <- sample(setdiff(aa, old_aa), 1L)
      n_rows            <- min(nrow(df), n_rows_to_delete)
      if (n_rows > 0L) {
        tr                              <- tail(seq_len(nrow(df)), n_rows)
        df[tr, as.character(site_idx)]  <- new_aa
        df$Delet[tr]                    <- "Yes"
        delet_records[[as.character(today)]] <- list(
          period = t, date = today, site = site_idx,
          old_aa = old_aa, new_aa = new_aa,
          rows   = row_start + tr - 1L
        )
      }
    }
    
    output <- rbind(output, ensure_cols(df))
  }
  
  # ── 8. Finalise output ───────────────────────────────────────────────────
  
  output[["._weight"]] <- NULL
  rownames(output)     <- NULL
  names(output)[seq_len(L)] <- as.character(seq_len(L))
  
  structure(
    list(
      Simulation_Output     = output,
      Variant_Details       = variant_details,
      Simulation_Dates      = sim_dates,
      Baseline_Ref_Sequence = base_seq,
      Delet_Records         = delet_records,
      Pool                  = pool
    ),
    class = "viralSim"
  )
}