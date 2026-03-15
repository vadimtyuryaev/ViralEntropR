#' @title Simulate Viral Variant Evolution
#' @description Generates synthetic time-series of viral sequences evolving from
#'   a reference strain, with realistic stochastic growth and multi-variant
#'   competition dynamics.
#' @export
#'
#' @details
#' Simulates the evolutionary dynamics of a viral population over discrete
#' monthly time steps.
#'
#' \strong{Simulation Phases:}
#' \enumerate{
#'   \item \strong{Reference Phase:} A stable period where only the reference
#'     strain exists. Reference months are sampled randomly from the reference
#'     pool to introduce baseline variability.
#'   \item \strong{Variant Emergence:} New variants are introduced at specified
#'     monthly intervals.
#'   \item \strong{Competition and Growth:} Each month the newest variant
#'     competes directly against its immediate predecessor ("pairwise"
#'     competition). Growth is stochastic — each variant's count is drawn from
#'     \code{Uniform(mult*(1-variability), mult*(1+variability))} — so setting
#'     \code{mutation_rate_variability > 0} introduces realistic fluctuation.
#'     Older (non-competing) variants contribute to the fill pool with
#'     exponentially increasing sampling weights.
#'   \item \strong{Deleterious Mutations:} Variants can carry fitness-reducing
#'     mutations capped at \code{n_deleterious_limit} sequences.
#'   \item \strong{Random Deletion Events:} Each month independently draws a
#'     Bernoulli coin flip; if it fires, \code{n_rows_to_delete} sequences at a
#'     randomly chosen site receive a random amino acid substitution.
#' }
#'
#' @param ref_sequences Character string (single reference sequence) or
#'   data.frame with a \code{Date} column and site columns.
#' @param n_ref_months Integer. Duration of the initial reference phase
#'   (months). Default \code{3}.
#' @param start_date Character or Date. Simulation start.
#' @param end_date Character or Date. Simulation end (inclusive month).
#' @param variants_config Integer vector. Number of mutations per variant,
#'   e.g., \code{c(2, 6)}.
#' @param variant_intervals Integer vector. Months between consecutive variant
#'   emergences (length \code{length(variants_config) - 1}).
#' @param n_new_mutations Integer. De-novo sequences produced at a variant's
#'   first appearance. Default \code{1}.
#' @param mutation_rate Numeric scalar or vector. Monthly growth multiplier per
#'   variant. Default \code{2}.
#' @param mutation_rate_variability Numeric in \code{[0, 1)}. Fractional
#'   spread around \code{mutation_rate}: growth is drawn from
#'   \code{Uniform(mult*(1-v), mult*(1+v))}. Set to \code{0} for deterministic
#'   growth. Default \code{0.25}.
#' @param deleterious_rate Numeric in \code{[0, 1]}. Per-mutation probability
#'   of being deleterious. Default \code{0}.
#' @param n_deleterious_limit Integer. Maximum sequences allowed to carry a
#'   deleterious mutation. Default \code{1}.
#' @param n_sequences_total Integer. Total sequences generated over the full
#'   mutation phase (spread evenly per month). Default \code{140}.
#' @param ref_variability Logical. If \code{TRUE}, introduces low-level
#'   variability at sites 1 and 2 of the reference pool. Default \code{FALSE}.
#' @param n_ref_sequences Integer. When \code{ref_sequences} is a single
#'   string, the number of reference rows to create. Default \code{100}.
#' @param prob_deletion_event Numeric in \code{[0, 1]}. Per-month probability
#'   of injecting a block of sequences with a random deletion at one site.
#'   Default \code{0}.
#' @param n_rows_to_delete Integer. Sequences affected when a deletion event
#'   fires. Default \code{0}.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#'   Default \code{NULL}.
#'
#' @return An object of class \code{"viralSim"}, a named list:
#' \item{Simulation_Output}{Data frame of all sequences with columns for each
#'   site plus \code{Variant}, \code{Phase}, \code{Date}, \code{Period},
#'   \code{Delet}.}
#' \item{Variant_Details}{List of per-variant metadata (emergence month,
#'   mutated positions, flags, sequence, growth rate).}
#' \item{Simulation_Dates}{Date vector of all simulated months.}
#' \item{Baseline_Ref_Sequence}{Character. The wild-type reference string.}
#' \item{Delet_Records}{Named list of deletion events, one entry per affected
#'   month.}
#' \item{Pool}{The fill pool data frame from the final multi-variant
#'   competition period (includes \code{._weight} column); \code{NULL} if
#'   only one variant was ever active.}
#'
#' @importFrom stats rbinom runif
#' @importFrom utils tail
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

  if (!is.null(seed)) set.seed(seed)

  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)

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

  # Use seq.Date with explicit end so both endpoints are included.
  # lubridate::interval() %/% period() gives the *distance* between two dates
  # (23 for Jan-2020 → Dec-2021), which is one short of the 24 monthly
  # time-points we need.  seq() inclusive matches simulate_variants_new_11.
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

    # ── 7d. Coin-flip deletion event ──────────────────────────────────────
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
