#' @title Simulate Viral Variant Evolution
#' @description Generates synthetic time-series of viral sequences evolving
#'   from a reference strain, with configurable multi-variant competition
#'   dynamics.
#' @export
#'
#' @details
#' Simulates the evolutionary dynamics of a viral population over discrete
#' monthly time steps. Designed to benchmark variant detection algorithms by
#' producing ground-truth datasets with known emergence events.
#'
#' \strong{Simulation Phases:}
#' \enumerate{
#'   \item \strong{Reference Phase}: A stable period where only the reference
#'     strain circulates.
#'   \item \strong{Single-Variant Emergence}: First variant grows stochastically
#'     against a reference background.
#'   \item \strong{Multi-Variant Competition}: When a second variant emerges,
#'     both the new and previous variants grow simultaneously; earlier variants
#'     fill the remaining background. A variant is marked as \emph{replaced}
#'     (\code{repl = TRUE}) when it no longer appears in a time step.
#'   \item \strong{Deleterious Mutations}: Variants can acquire fitness costs
#'     that cap their spread via \code{del_cap}.
#'   \item \strong{Random Deletion Events}: Simulates sequencing errors or
#'     dropouts at a given probability per period.
#' }
#'
#' @param ref_sequences Character string (single sequence) or data frame with
#'   numeric site columns and a \code{Date} column.
#' @param n_ref_months Integer. Duration of the initial reference phase in
#'   months. Default is \code{3}.
#' @param start_date Date or character. Simulation start date.
#' @param end_date Date or character. Simulation end date.
#' @param variants_config List of scalars (or integer vector). Each element
#'   gives the number of mutations introduced by the corresponding variant
#'   (e.g. \code{list(2, 3)} — variant 1 has 2 mutations, variant 2 has 3).
#' @param variant_intervals Integer vector. Gap in months between consecutive
#'   variant emergences. Recycled if shorter than \code{length(variants_config) - 1}.
#' @param n_new_mutations Integer. Number of sequences seeded at a variant's
#'   first appearance. Default is \code{1}.
#' @param mutation_rate Numeric scalar or vector. Growth multiplier per variant.
#'   If scalar, applied to all variants. Default is \code{2}.
#' @param mutation_rate_variability Numeric. Relative spread around the
#'   multiplier for stochastic growth:
#'   \code{runif(1, mult*(1-v), mult*(1+v))}. Default is \code{0.25}.
#' @param deleterious_rate Numeric \eqn{[0,1]}. Probability each mutation
#'   site is deleterious. Default is \code{0}.
#' @param n_deleterious_limit Integer. Maximum sequences allowed at a
#'   deleterious site before capping. Default is \code{1}.
#' @param n_sequences_total Integer. Target total sequences per month.
#'   Default is \code{140}.
#' @param ref_variability Logical. If \code{TRUE}, injects low-level noise into
#'   the reference pool at positions 1 and 2. Default is \code{FALSE}.
#' @param n_ref_sequences Integer. Number of reference sequences in the pool.
#'   Default is \code{100}.
#' @param prob_deletion_event Numeric \eqn{[0,1]}. Probability of a random
#'   deletion event per period. Default is \code{0}.
#' @param n_rows_to_delete Integer. Sequences affected per deletion event.
#'   Default is \code{0}.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#'
#' @return An object of class \code{"viralSim"} — a named list with:
#' \item{Simulation_Output}{Data frame of all sequences across all time periods,
#'   with columns for each site, plus \code{Variant}, \code{Phase},
#'   \code{Date}, \code{Period}, and \code{Delet}.}
#' \item{Variant_Details}{List of per-variant tracking objects (\code{em},
#'   \code{pos}, \code{flags}, \code{vseq}, \code{mult}, \code{last},
#'   \code{cum}, \code{repl}).}
#' \item{Simulation_Dates}{Date vector of all simulated months.}
#' \item{Baseline_Ref_Sequence}{Character string of the reference sequence.}
#' \item{Delet_Records}{Named list of deletion events, keyed by date.}
#' \item{Pool}{The background sampling pool used in the last multi-variant
#'   period, or \code{NULL} if only one variant was ever active.}
#'
#' @importFrom lubridate interval period
#' @importFrom stats rbinom runif
#' @importFrom utils tail
#'
#' @examples
#' ref_seq = "MKTIIALSYIFCLVFA"
#' sim = simulate_variant_evolution(
#'   ref_sequences     = ref_seq,
#'   n_ref_months      = 3,
#'   start_date        = "2021-01-01",
#'   end_date          = "2021-12-01",
#'   variants_config   = list(2, 3),
#'   variant_intervals = 4,
#'   n_sequences_total = 50,
#'   mutation_rate     = 1.5,
#'   seed              = 123
#' )
#' head(sim$Simulation_Output)
#' sim$Variant_Details[[1]]$cum   # cumulative count for variant 1
simulate_variant_evolution = function(
    ref_sequences,
    n_ref_months              = 3,
    start_date,
    end_date,
    variants_config,
    variant_intervals,
    n_new_mutations           = 1,
    mutation_rate             = 2,
    mutation_rate_variability = 0.25,
    deleterious_rate          = 0,
    n_deleterious_limit       = 1,
    n_sequences_total         = 140,
    ref_variability           = FALSE,
    n_ref_sequences           = 100,
    prob_deletion_event       = 0,
    n_rows_to_delete          = 0,
    seed                      = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  start_date = as.Date(start_date)
  end_date   = as.Date(end_date)
  
  aa = c("A","R","N","D","C","E","Q","G","H","I",
          "L","K","M","F","P","S","T","W","Y","V")
  
  # --- Internal helpers ------------------------------------------------------
  
  # Coerce to matrix safely, preserving dimnames
  ensure_matrix = function(mat) {
    if (!is.matrix(mat)) {
      vals = as.vector(mat)
      nms  = if (!is.null(names(vals))) names(vals) else seq_along(vals)
      mat  = matrix(vals, nrow = 1L, ncol = length(vals),
                     dimnames = list(NULL, nms))
    }
    mat
  }
  
  # Cap deleterious variant alleles to at most `limit` sequences per flagged site
  del_cap = function(mat, vidx, vdet, limit, ref_seq_chars) {
    mat = ensure_matrix(mat)
    vd  = vdet[[vidx]]
    if (!any(vd$flags)) return(mat)
    for (p in vd$pos[vd$flags]) {
      allele = vd$vseq[p]
      cands  = which(mat[, p] == allele)
      cnt    = length(cands)
      if (cnt > limit && cnt > (cnt - limit)) {
        mat[sample(cands, cnt - limit), p] = ref_seq_chars[p]
      }
    }
    mat
  }
  
  # Build a single-row data frame for a variant sequence (used in pool)
  mk_row = function(seq_vec, label, seq_cols, wt, template_cols) {
    df = as.data.frame(
      matrix(seq_vec, 1L, dimnames = list(NULL, seq_cols)),
      stringsAsFactors = FALSE
    )
    df$Variant   = label
    df$Phase     = paste0("M", as.numeric(sub("Variant_", "", label)))
    df$Date      = as.Date(NA)
    df$Period    = NA_integer_
    df$Delet     = "No"
    df$._weight  = wt
    df[, template_cols, drop = FALSE]
  }
  
  # Ensure all template_cols exist in df; closes over template_cols
  ensure_cols = function(df) {
    miss = setdiff(template_cols, names(df))
    df[miss] = NA
    df[, template_cols, drop = FALSE]
  }
  
  # Stochastic growth count
  grow = function(last, mult, variability) {
    max(1L, ceiling(last * stats::runif(
      1L,
      mult * (1 - variability),
      mult * (1 + variability)
    )))
  }
  
  # --- Reference setup -------------------------------------------------------
  if (is.character(ref_sequences) && length(ref_sequences) == 1L) {
    base_seq = ref_sequences
    L        = nchar(base_seq)
    ref_df   = as.data.frame(
      matrix(strsplit(base_seq, "")[[1]], nrow = n_ref_sequences, ncol = L,
             byrow = TRUE),
      stringsAsFactors = FALSE
    )
    colnames(ref_df) = as.character(seq_len(L))
    ref_df$Date = rep(
      seq.Date(start_date, by = "month", length.out = n_ref_months),
      length.out = n_ref_sequences
    )
  } else {
    ref_df   = ref_sequences
    seq_num  = grep("^[0-9]+$", names(ref_df))
    base_seq = paste(ref_df[1L, seq_num], collapse = "")
    L        = nchar(base_seq)
  }
  seq_cols       = setdiff(names(ref_df), "Date")
  K              = length(seq_cols)
  ref_df$Date    = as.Date(ref_df$Date)
  ref_df$Month   = format(ref_df$Date, "%Y-%m")
  
  # --- Reference variability -------------------------------------------------
  if (ref_variability) {
    sub1 = sample(setdiff(aa, substr(base_seq, 1L, 1L)), 1L)
    sub2 = sample(setdiff(aa, substr(base_seq, 2L, 2L)), 1L)
    for (m in unique(ref_df$Month)) {
      idx = sort(which(ref_df$Month == m))
      # Use string column names to avoid corrupting metadata columns
      if (length(idx) >= 4L) ref_df[idx[2:4], "1"] = sub1
      if (length(idx) >= 6L) ref_df[idx[5:6], "2"] = sub2
    }
  }
  ref_df$Month = NULL
  
  # --- Simulation timeline ---------------------------------------------------
  # Use period() throughout — avoids lubridate::months() namespace export issue
  n_months_sim = lubridate::interval(start_date, end_date) %/%
    lubridate::period(1L, "months")
  sim_dates  = seq.Date(start_date, by = "month", length.out = n_months_sim)
  n_periods  = length(sim_dates)
  
  template_cols = c(seq_cols, "Variant", "Phase", "Date",
                     "Period", "Delet", "._weight")
  
  # --- Build reference phase -------------------------------------------------
  # Randomly sample reference months (reproduces original stochastic behaviour)
  choice = sample(unique(format(ref_df$Date, "%Y-%m")),
                   n_ref_months, replace = TRUE)
  
  ref_phase = do.call(rbind, lapply(seq_len(n_ref_months), function(i) {
    df         = ref_df[format(ref_df$Date, "%Y-%m") == choice[i], ]
    df$Date    = sim_dates[i]
    df$Period  = i
    df$Variant = "Reference"
    df
  }))
  # Label as "R" (exact reference) or "RV" (reference with variability noise)
  ref_phase$Phase = ifelse(
    apply(ref_phase[, seq_cols], 1L, paste, collapse = "") == base_seq,
    "R", "RV"
  )
  ref_phase$Delet    = "No"
  ref_phase          = ensure_cols(ref_phase)
  ref_phase$._weight[ref_phase$Variant == "Reference"] = 1
  
  # --- Define variants -------------------------------------------------------
  # variants_config: list of scalars (number of mutations per variant)
  variants_config = as.list(variants_config)   # accept vector or list
  nv   = length(variants_config)
  gaps = rep(variant_intervals, length.out = max(0L, nv - 1L))
  emerg        = numeric(nv)
  emerg[1L]    = n_ref_months + 1L
  if (nv > 1L) for (i in 2L:nv) emerg[i] = emerg[i - 1L] + gaps[i - 1L]
  valid        = emerg <= n_periods
  emerg        = emerg[valid]
  variants_config = variants_config[valid]
  rates        = if (length(mutation_rate) == 1L) rep(mutation_rate, nv) else mutation_rate
  
  ref_seq_chars = as.character(ref_phase[1L, seq_cols])
  
  # Build per-variant tracking objects (vdet) and single-row pool rows (vrow)
  n_valid = length(emerg)
  vdet    = vector("list", n_valid)
  vrow    = vector("list", n_valid)
  
  for (i in seq_len(n_valid)) {
    nm   = as.integer(variants_config[[i]])      # number of mutations
    pos  = sample(setdiff(seq_len(K), 1:2), nm)
    new  = sapply(pos, function(p)
      sample(setdiff(aa, ref_seq_chars[p]), 1L))
    vseq = ref_seq_chars; vseq[pos] = new
    flags = stats::rbinom(nm, 1L, deleterious_rate) == 1L
    vdet[[i]] = list(em   = emerg[i], pos  = pos, flags = flags,
                      vseq = vseq,     mult = rates[i],
                      last = 0L,       cum  = 0L, repl  = FALSE)
    vrow[[i]] = mk_row(vseq, paste0("Variant_", i),
                        seq_cols, wt = 2^(i + 1L), template_cols)
  }
  
  n_fixed = ceiling(n_sequences_total / (n_periods - n_ref_months))
  
  # --- Simulate mutation phase -----------------------------------------------
  output        = ref_phase
  delet_records = list()
  pool          = NULL     # retained for return value
  
  for (t in (n_ref_months + 1L):n_periods) {
    today     = sim_dates[t]
    row_start = nrow(output) + 1L
    
    # Highest-indexed variant that has emerged and has not been replaced
    active = max(which(sapply(vdet, function(v) t >= v$em && !v$repl)),
                  0L)
    
    # ── Case 0: pre-variant period ──────────────────────────────────────── #
    if (active == 0L) {
      df          = ref_phase[sample(nrow(ref_phase), n_fixed, TRUE,
                                      ref_phase$._weight), template_cols]
      df$Date     = today
      df$Period   = t
      
      # ── Case 1: single active variant ───────────────────────────────────── #
    } else if (active == 1L) {
      n1 = if (t == vdet[[1L]]$em) n_new_mutations else
        grow(vdet[[1L]]$last, vdet[[1L]]$mult, mutation_rate_variability)
      vdet[[1L]]$last = n1
      vdet[[1L]]$cum  = vdet[[1L]]$cum + n1
      
      m1   = matrix(vdet[[1L]]$vseq, nrow = n1, ncol = K, byrow = TRUE,
                     dimnames = list(NULL, seq_cols))
      left = n_fixed - n1
      if (left > 0L) {
        samp = sample(nrow(ref_phase), left, TRUE, ref_phase$._weight)
        mat  = rbind(m1, as.matrix(ref_phase[samp, seq_cols]))
        lbl  = c(rep("Variant_1", n1), ref_phase$Variant[samp])
      } else {
        mat = m1; lbl = rep("Variant_1", n1)
      }
      mat[seq_len(min(n1, nrow(mat))), ] =
        del_cap(mat[seq_len(min(n1, nrow(mat))), , drop = FALSE],
                1L, vdet, n_deleterious_limit, ref_seq_chars)
      
      df          = as.data.frame(mat, stringsAsFactors = FALSE)
      df$Variant  = lbl
      df$Phase    = "M1"
      df$Date     = today; df$Period = t; df$Delet = "No"
      
      # ── Case 2+: multi-variant competition ──────────────────────────────── #
    } else {
      cur = active; cmp = cur - 1L
      
      nc = if (t == vdet[[cur]]$em) n_new_mutations else
        grow(vdet[[cur]]$last, vdet[[cur]]$mult, mutation_rate_variability)
      vdet[[cur]]$last = nc; vdet[[cur]]$cum = vdet[[cur]]$cum + nc
      
      np = if (t == vdet[[cmp]]$em) n_new_mutations else
        grow(vdet[[cmp]]$last, vdet[[cmp]]$mult, mutation_rate_variability)
      vdet[[cmp]]$last = np; vdet[[cmp]]$cum = vdet[[cmp]]$cum + np
      
      mc  = matrix(vdet[[cur]]$vseq, nrow = nc, ncol = K, byrow = TRUE,
                    dimnames = list(NULL, seq_cols))
      mp  = matrix(vdet[[cmp]]$vseq, nrow = np, ncol = K, byrow = TRUE,
                    dimnames = list(NULL, seq_cols))
      mat = rbind(mc, mp)
      lbl = c(rep(paste0("Variant_", cur), nc),
               rep(paste0("Variant_", cmp), np))
      
      left = n_fixed - nrow(mat)
      if (left > 0L) {
        # Background pool: ref rows + any superseded older variants
        pool = if (ref_variability) ref_phase[c(1L,2L,5L), ] else
          ref_phase[1L, , drop = FALSE]
        if (cur >= 3L) {
          for (j in seq_len(cur - 2L)) {
            vr         = vrow[[j]]
            vr$Date    = today; vr$Period = t
            vr$Phase   = paste0("M", cur)
            pool       = rbind(pool, vr)
          }
        }
        samp = sample(nrow(pool), left, TRUE, pool$._weight)
        mat  = rbind(mat, as.matrix(pool[samp, seq_cols]))
        lbl  = c(lbl, pool$Variant[samp])
      }
      if (nrow(mat) > n_fixed) {
        mat = mat[seq_len(n_fixed), , drop = FALSE]
        lbl = lbl[seq_len(n_fixed)]
      }
      
      # Apply del_cap to current, then to previous variant rows
      idxc = seq_len(min(nc, nrow(mat)))
      mat[idxc, ] = del_cap(mat[idxc, , drop = FALSE],
                             cur, vdet, n_deleterious_limit, ref_seq_chars)
      cs = min(nc + 1L, nrow(mat)); ce = min(nc + np, nrow(mat))
      if (cs <= ce)
        mat[cs:ce, ] = del_cap(mat[cs:ce, , drop = FALSE],
                                cmp, vdet, n_deleterious_limit, ref_seq_chars)
      
      df          = as.data.frame(mat, stringsAsFactors = FALSE)
      df$Variant  = lbl
      df$Phase    = paste0("M", cur)
      df$Date     = today; df$Period = t; df$Delet = "No"
      
      # Mark previous variant as replaced if it no longer appears
      if (!any(lbl == paste0("Variant_", cmp))) vdet[[cmp]]$repl = TRUE
    }
    
    # ── Random deletion event ───────────────────────────────────────────── #
    if (stats::runif(1L) < prob_deletion_event) {
      site_idx = sample(setdiff(seq_along(seq_cols), 1:2), 1L)
      old_aa   = df[1L, site_idx]
      new_aa   = sample(setdiff(aa, old_aa), 1L)
      to_mut   = utils::tail(seq_len(nrow(df)), n_rows_to_delete)
      df[to_mut, site_idx]  = new_aa
      df$Delet[to_mut]      = "Yes"
      delet_records[[as.character(today)]] = list(
        period = t, date = today,
        site   = site_idx, old_aa = old_aa, new_aa = new_aa,
        rows   = row_start + to_mut - 1L
      )
    }
    
    output = rbind(output, ensure_cols(df))
  }
  
  # --- Final cleanup ---------------------------------------------------------
  output             = output[, setdiff(names(output), "._weight")]
  names(output)[seq_len(K)] = as.character(seq_len(K))
  rownames(output)   = NULL
  
  structure(
    list(
      Simulation_Output     = output,
      Variant_Details       = vdet,
      Simulation_Dates      = sim_dates,
      Baseline_Ref_Sequence = base_seq,
      Delet_Records         = delet_records,
      Pool                  = pool
    ),
    class = "viralSim"
  )
}