#' @title Run Variant Detection Sample Size Study
#' @description Orchestrates a comprehensive simulation benchmark to evaluate
#'   the sensitivity and accuracy of the variant detection pipeline across
#'   biological scenarios and sample sizes.
#'
#' @details
#' For each scenario and each value of \code{n_new_mut_seq_vec} the function:
#' \enumerate{
#'   \item \strong{Simulates} a viral population via
#'     \code{\link{simulate_variant_evolution}}.
#'   \item \strong{Encodes and partitions} the sequences into time windows via
#'     \code{\link{encode_aa_sequence}} and
#'     \code{\link{partition_time_windows}}.
#'   \item \strong{Computes Hellinger distances} across sites and partitions via
#'     \code{\link{calculate_hellinger_matrix}}.
#'   \item \strong{Detects change points} globally and per-variant via
#'     \code{\link[ecp]{e.agglo}}.
#'   \item \strong{Evaluates detection} by checking whether each variant's
#'     mutation sites appear in the highest-entropy cluster at the expected
#'     partition (or later).
#' }
#'
#' Two complementary detection signals are recorded per variant:
#' \itemize{
#'   \item \strong{Site-based}: mutation sites found in top entropy cluster.
#'   \item \strong{Change-point-based}: \code{e.agglo} applied globally
#'     (full timeline) and locally (up to \code{part_em - 1}).
#' }
#'
#' When \code{mc.cores > 1} and the platform is not Windows, the inner sweep
#' over \code{n_new_mut_seq_vec} is parallelised using
#' \code{parallel::mclapply} (Linux fork, copy-on-write, no export overhead).
#' Each forked worker handles one \code{n_new_seq} value independently,
#' reducing peak memory from O(length(n_new_mut_seq_vec)) to O(1) per call.
#' On Windows or when \code{mc.cores = 1}, a sequential \code{lapply} is
#' used automatically.
#'
#' @param ref_seq Character. The reference amino acid sequence string.
#' @param variants_list List. Each element is an integer vector of mutation
#'   counts per variant for that scenario (e.g. \code{list(c(2), c(2,3))}).
#' @param mutation_rate_list List. Each element is a numeric vector of growth
#'   multipliers matching \code{variants_list}.
#' @param n_new_mut_seq_vec Integer vector. Range of initial mutant sequence
#'   counts to sweep over.
#' @param num_months_ref_seq Integer. Duration of the reference phase in months.
#' @param variant_interval Integer vector. Months between consecutive variant
#'   emergences; recycled as needed.
#' @param start_date Character or Date. Simulation start. Default
#'   \code{"2020-01-01"}.
#' @param end_date Character or Date. Simulation end. Default
#'   \code{"2021-12-01"}.
#' @param mutation_rate_variability Numeric. Stochastic spread around the
#'   growth multiplier. Default \code{0.0}.
#' @param deleterious_rate Numeric \eqn{[0,1]}. Probability a mutation is
#'   deleterious. Default \code{0.0}.
#' @param ref_variability Logical. Inject low-level noise into reference pool.
#'   Default \code{TRUE}.
#' @param n_seq_per_month Integer. Target sequences per month. Default
#'   \code{1000}.
#' @param prob_delet Numeric. Probability of a random deletion event per
#'   period. Default \code{0}.
#' @param n_delet Integer. Sequences affected per deletion event (also used as
#'   the deleterious cap limit). Default \code{0}.
#' @param sliding_window_length Integer. Partition window size in months.
#'   Default \code{2}.
#' @param window_option Integer. Window type passed to
#'   \code{\link{partition_time_windows}} (1=Cumulative, 2=Sliding,
#'   3=Non-overlapping). Default \code{3}.
#' @param save_html Logical. If \code{TRUE}, saves an HTML table to disk.
#'   Default \code{TRUE}.
#' @param output_file Character. Filename for the HTML report. Default
#'   \code{"variants_detection_study.html"}.
#' @param mc.cores Integer. Number of cores for the inner
#'   \code{n_new_mut_seq_vec} sweep via \code{parallel::mclapply}.
#'   Defaults to \code{1L} (sequential). Values greater than 1 are silently
#'   clamped to 1 on Windows. On Linux HPC set to the number of available
#'   cores (e.g. \code{12L}).
#' @param ... Additional arguments passed to
#'   \code{\link{partition_time_windows}} (and on to
#'   \code{\link{cluster_sites_by_entropy}}).
#'
#' @return A named list:
#' \item{Sim_List}{Nested list \code{[[scenario]][[n_new_seq]]} of
#'   \code{viralSim} objects.}
#' \item{Part_Data}{Nested list \code{[[scenario]][[n_new_seq]]} of
#'   \code{partition_time_windows} outputs.}
#' \item{Results}{Data frame with one row per (scenario × n_new_seq × variant),
#'   containing site-level and change-point detection results with HTML
#'   highlighting.}
#' \item{First_Detect_Sample}{Integer vector (length = n_scenarios). First
#'   \code{n_new_seq} at which all variants in the scenario were simultaneously
#'   detected with distinct site sets.}
#' \item{CP_Accuracy}{Numeric vector (length = n_scenarios). Mean proportion
#'   of variants whose actual CP was correctly identified.}
#' \item{Table}{The \code{kableExtra} HTML table object.}
#'
#' @examples
#' \dontrun{
#' res <- run_detection_study(
#'   ref_seq               = "MKTIIALSYI",
#'   variants_list         = list(c(2)),
#'   mutation_rate_list    = list(c(2.0)),
#'   n_new_mut_seq_vec     = c(2L, 5L),
#'   num_months_ref_seq    = 2L,
#'   variant_interval      = 3L,
#'   start_date            = "2021-01-01",
#'   end_date              = "2021-10-01",
#'   n_seq_per_month       = 50L,
#'   sliding_window_length = 2L,
#'   window_option         = 3L,
#'   save_html             = FALSE,
#'   mc.cores              = 1L
#' )
#'
#' print(res$Results[, c("variant", "num_denovo_sequences",
#'                        "actual_partition", "detection_partition",
#'                        "time_to_detection")])
#' cat("CP Accuracy:", res$CP_Accuracy, "\n")
#' cat("First distinct detection at n_new =",
#'     res$First_Detect_Sample, "sequences\n")
#' }
#'
#' @importFrom kableExtra kable_styling save_kable
#' @importFrom ecp e.agglo
#' @importFrom knitr kable
#' @importFrom utils tail
#' @importFrom magrittr %>%
#'
#' @export
run_detection_study <- function(
    ref_seq,
    variants_list,
    mutation_rate_list,
    n_new_mut_seq_vec,
    num_months_ref_seq,
    variant_interval,
    start_date                = "2020-01-01",
    end_date                  = "2021-12-01",
    mutation_rate_variability = 0.0,
    deleterious_rate          = 0.0,
    ref_variability           = TRUE,
    n_seq_per_month           = 1000,
    prob_delet                = 0,
    n_delet                   = 0,
    sliding_window_length     = 2,
    window_option             = 3,
    save_html                 = TRUE,
    output_file               = "variants_detection_study.html",
    mc.cores                  = 1L,
    ...
) {
  
  # --- 0. Setup --------------------------------------------------------------
  n_scen                <- length(variants_list)
  first_distinct_detect <- rep(NA_real_, n_scen)
  cp_accuracy           <- rep(NA_real_, n_scen)
  sim_res_list          <- vector("list", n_scen)
  part_data_list        <- vector("list", n_scen)
  all_rows              <- vector("list", n_scen)   # collected before rbind
  
  n_col      <- nchar(ref_seq)
  sim_dates  <- seq.Date(as.Date(start_date), as.Date(end_date), by = "month")
  total_seqs <- (length(sim_dates) - num_months_ref_seq) * n_seq_per_month
  n_ref      <- num_months_ref_seq * n_seq_per_month
  
  # Clamp mc.cores to 1 on Windows — forking is not available
  mc.cores <- if (.Platform$OS.type == "windows") 1L
  else max(1L, as.integer(mc.cores))
  
  # --- 1. Shared helper functions --------------------------------------------
  
  sort_string <- function(x) {
    if (is.na(x) || x == "" || x == "-") return(x)
    paste(sort(as.numeric(unlist(strsplit(x, ",")))), collapse = ",")
  }
  
  highlight_sites_advanced <- function(actual, detected, deleterious) {
    a  <- unlist(strsplit(actual,   ","))
    d  <- unlist(strsplit(detected, ","))
    ah <- vapply(a, function(s)
      if (s %in% d) sprintf('<span style="background-color:yellow;">%s</span>', s) else s,
      character(1L))
    dh <- vapply(d, function(s)
      if (s %in% a) sprintf('<span style="background-color:yellow;">%s</span>', s) else s,
      character(1L))
    list(sites    = paste(ah, collapse = ","),
         detected = paste(dh, collapse = ","),
         deleterious = deleterious)
  }
  
  highlight_actual_cp <- function(actual, detected_var_raw) {
    if (is.na(actual) || actual == "") return(as.character(actual))
    det_vals <- suppressWarnings(
      as.numeric(unlist(strsplit(detected_var_raw, ",")))
    )
    if (actual %in% det_vals)
      sprintf('<span style="background-color:yellow;">%s</span>', actual)
    else as.character(actual)
  }
  
  highlight_detected_cp_var <- function(detected_var_raw, actual) {
    parts <- unlist(strsplit(detected_var_raw, ","))
    out   <- vapply(parts, function(p) {
      if (nzchar(p) && !is.na(actual) &&
          suppressWarnings(as.numeric(p)) == actual)
        sprintf('<span style="background-color:yellow;">%s</span>', p)
      else p
    }, character(1L))
    paste(out, collapse = ",")
  }
  
  # --- 2. HTML caption -------------------------------------------------------
  caption_text <- paste0(
    "Variants Detection Study<br/>",
    paste(vapply(seq_along(variants_list), function(i)
      paste0("Scenario ", i, ": [",
             paste(variants_list[[i]],      collapse = ","), "] ; [",
             paste(mutation_rate_list[[i]], collapse = ","), "]"),
      character(1L)), collapse = "<br/>")
  )
  
  # --- 3. Setup -------------------------------------------------------------
  extra_args <- list(...)
  
  # --- 4. Scenario loop (sequential — scenarios are typically few) -----------
  for (scenario_idx in seq_along(variants_list)) {
    
    number_of_variants <- variants_list[[scenario_idx]]
    mutation_rates     <- mutation_rate_list[[scenario_idx]]
    nv                 <- length(number_of_variants)
    vi                 <- variant_interval
    if (length(vi) < nv - 1L)
      vi <- c(vi, rep(utils::tail(vi, 1L), nv - 1L - length(vi)))
    
    # ── Per-n_new_seq worker (closes over scenario variables) ──────────────
    run_one_n_new <- function(n_new_seq) {
      
      # 1. Simulate -----------------------------------------------------------
      sim_res <- simulate_variant_evolution(
        ref_sequences             = ref_seq,
        n_ref_months              = num_months_ref_seq,
        start_date                = start_date,
        end_date                  = end_date,
        variants_config           = number_of_variants,
        variant_intervals         = vi,
        n_new_mutations           = n_new_seq,
        mutation_rate             = mutation_rates,
        mutation_rate_variability = mutation_rate_variability,
        deleterious_rate          = deleterious_rate,
        n_deleterious_limit       = max(0L, n_delet),
        n_sequences_total         = total_seqs,
        ref_variability           = ref_variability,
        n_ref_sequences           = n_ref,
        prob_deletion_event       = prob_delet,
        n_rows_to_delete          = n_delet
      )
      
      # 2. Encode + partition -------------------------------------------------
      mat_sim  <- sim_res$Simulation_Output
      enc_mat  <- encode_aa_sequence(as.matrix(mat_sim[, seq_len(n_col)]))
      AL_df    <- as.data.frame(enc_mat)
      AL_df[]  <- lapply(AL_df, as.integer)
      AL_df$Date <- as.Date(format(as.Date(mat_sim$Date), "%Y-%m-01"))
      
      part_data <- do.call(partition_time_windows,
                           c(list(data          = AL_df,
                                  n_sites       = n_col,
                                  window_length = sliding_window_length,
                                  window_type   = window_option,
                                  start_date    = start_date,
                                  end_date      = end_date),
                             extra_args))
      
      # 3. Relabel all cluster DataFrames once upfront -------------------------
      relabeled_clusters <- lapply(part_data$Clusters, function(cl) {
        cl$DataFrame <- relabel_entropy_classes(cl$DataFrame)
        cl
      })
      
      # 4. Hellinger matrix + transpose for e.agglo --------------------------
      hell_mat <- calculate_hellinger_matrix(
        partitions = part_data$Partitions,
        sites      = seq_len(n_col),
        aa_levels  = 25L
      )
      dat_t  <- t(hell_mat)
      
      # 5. Global change point detection (full timeline) ----------------------
      max_cp <- length(part_data$Clusters) - 1L
      if (max_cp >= 2L) {
        cp_all_raw <- ecp::e.agglo(
          X      = dat_t,
          member = seq_len(nrow(dat_t)),
          alpha  = 1,
          penalty = function(cps) 0
        )$estimates
        cp_all <- cp_all_raw[cp_all_raw > 1L & cp_all_raw <= max_cp]
      } else {
        cp_all <- integer(0L)
      }
      detected_cp_all <- paste(sort(unique(cp_all)), collapse = ",")
      
      # 6. Per-variant evaluation ---------------------------------------------
      details            <- sim_res$Variant_Details
      delet              <- sim_res$Delet_Records
      n_variants         <- length(details)
      detected_flags_all <- logical(n_variants)
      sites_list_all     <- vector("list", n_variants)
      rows_this          <- vector("list", n_variants)
      
      for (v_idx in seq_along(details)) {
        info  <- details[[v_idx]]
        sites <- sort(info$pos)
        em    <- info$em
        
        part_em <- ceiling(em / sliding_window_length)
        
        if (part_em > length(relabeled_clusters)) {
          warning(sprintf(
            "Variant %d: part_em (%d) exceeds number of partitions (%d). ",
            v_idx, part_em, length(relabeled_clusters),
            "Check that end_date is set one month beyond the last simulation ",
            "month. Clamping part_em to the last available partition."
          ))
          part_em <- length(relabeled_clusters)
        }
        
        cl_em     <- relabeled_clusters[[part_em]]$DataFrame
        top_sites <- if (nrow(cl_em) > 0L)
          cl_em[as.numeric(cl_em$class) == 1L, ]$sites
        else integer(0L)
        
        date_em          <- as.character(sim_res$Simulation_Dates[em])
        deleterious_site <- if (date_em %in% names(delet))
          delet[[date_em]]$site else NA
        
        dp <- part_em
        if (!all(setdiff(sites, deleterious_site) %in% top_sites)) {
          dp <- NA
          n_clusters <- length(relabeled_clusters)
          if (part_em < n_clusters) {
            for (i in (part_em + 1L):n_clusters) {
              cl_i  <- relabeled_clusters[[i]]$DataFrame
              top_i <- if (nrow(cl_i) > 0L)
                cl_i[as.numeric(cl_i$class) == 1L, ]$sites
              else integer(0L)
              if (all(setdiff(sites, deleterious_site) %in% top_i)) {
                dp <- i; break
              }
            }
          }
          if (is.na(dp)) dp <- "-"
        }
        
        detected_flags_all[v_idx] <- isTRUE(dp == part_em)
        detected_sites_str        <- paste(sort(top_sites), collapse = ",")
        sites_list_all[[v_idx]]   <- detected_sites_str
        
        if (part_em > 2L && max_cp >= 2L) {
          cp_var_raw <- ecp::e.agglo(
            X       = dat_t[seq_len(part_em - 1L), , drop = FALSE],
            member  = seq_len(part_em - 1L),
            alpha   = 1,
            penalty = function(cps) 0
          )$estimates
          cp_var <- cp_var_raw[
            cp_var_raw > 1L &
              cp_var_raw <= max_cp &
              cp_var_raw <= (part_em - 1L)
          ]
          detected_cp_var <- paste(sort(unique(cp_var)), collapse = ",")
        } else {
          detected_cp_var <- ""
        }
        
        del_str           <- if (!is.na(deleterious_site))
          as.character(deleterious_site) else ""
        detection_month   <- if (is.numeric(dp)) dp * sliding_window_length else NA
        time_to_detection <- if (is.numeric(dp)) detection_month - em       else NA
        
        rows_this[[v_idx]] <- data.frame(
          scenario             = scenario_idx,
          variant              = paste0("Variant_", v_idx),
          num_denovo_sequences = n_new_seq,
          sites_raw            = paste(sites, collapse = ","),
          detected_sites_raw   = detected_sites_str,
          deleterious_sites    = del_str,
          num_months_ref_seq   = num_months_ref_seq,
          actual_month         = em,
          detection_month      = detection_month,
          time_to_detection    = time_to_detection,
          actual_partition     = part_em,
          detection_partition  = dp,
          actual_cp_raw        = part_em - 1L,
          detected_cp_var_raw  = detected_cp_var,
          detected_cp_all      = detected_cp_all,
          stringsAsFactors     = FALSE
        )
      }
      
      list(
        rows           = do.call(rbind, rows_this),
        fully_detected = all(detected_flags_all),
        sites_list     = sites_list_all,
        sim_res        = sim_res,
        part_data      = part_data
      )
    } # end run_one_n_new
    
    # ── Dispatch: parallel on Linux, sequential on Windows -----------------
    # mc.cores > 1 only ever reached on Linux (clamped above on Windows).
    # mclapply forks the current process — all variables in scope are
    # available in each child via copy-on-write with zero export overhead.
    # Each worker handles exactly one n_new_seq value, so peak memory per
    # worker is O(1 simulation) rather than O(length(n_new_mut_seq_vec)).
    iter_results <- if (mc.cores > 1L) {
      parallel::mclapply(
        n_new_mut_seq_vec,
        run_one_n_new,
        mc.cores    = mc.cores,
        mc.set.seed = TRUE
      )
    } else {
      lapply(n_new_mut_seq_vec, run_one_n_new)
    }
    
    # Normalise any crashed workers (try-error from mclapply) to NULL rows
    # so downstream rbind does not encounter unexpected types.
    iter_results <- lapply(iter_results, function(r) {
      if (inherits(r, "try-error"))
        list(rows = NULL, fully_detected = FALSE,
             sites_list = list(), sim_res = NULL, part_data = NULL)
      else r
    })
    
    # ── Accumulate results ────────────────────────────────────────────────
    sim_res_list[[scenario_idx]]   <- setNames(
      lapply(iter_results, `[[`, "sim_res"),
      as.character(n_new_mut_seq_vec)
    )
    part_data_list[[scenario_idx]] <- setNames(
      lapply(iter_results, `[[`, "part_data"),
      as.character(n_new_mut_seq_vec)
    )
    all_rows[[scenario_idx]] <- do.call(rbind, lapply(iter_results, `[[`, "rows"))
    
    # ── First distinct detection threshold ────────────────────────────────
    for (k in seq_along(n_new_mut_seq_vec)) {
      if (iter_results[[k]]$fully_detected) {
        vec <- unlist(strsplit(
          unlist(iter_results[[k]]$sites_list), ","
        ))
        if (length(unique(vec)) == length(vec) &&
            is.na(first_distinct_detect[scenario_idx])) {
          first_distinct_detect[scenario_idx] <- n_new_mut_seq_vec[k]
        }
      }
    }
    
  } # end scenario loop
  
  # --- 5. Combine all rows (one rbind, not incremental) ---------------------
  results_df <- do.call(rbind, all_rows)
  rownames(results_df) <- NULL
  
  # --- 6. CP accuracy per scenario ------------------------------------------
  for (i in seq_len(n_scen)) {
    sub     <- results_df[results_df$scenario == i, ]
    correct <- mapply(function(acp, dcp_raw) {
      vals <- suppressWarnings(
        as.numeric(unlist(strsplit(dcp_raw, ",")))
      )
      isTRUE(acp %in% vals)
    }, sub$actual_cp_raw, sub$detected_cp_var_raw)
    cp_accuracy[i] <- mean(correct, na.rm = TRUE)
  }
  
  # --- 7. Post-processing: sort and HTML highlighting -----------------------
  results_df$sites          <- vapply(results_df$sites_raw,          sort_string, character(1L))
  results_df$detected_sites <- vapply(results_df$detected_sites_raw, sort_string, character(1L))
  
  results_df$actual_cp <- mapply(
    highlight_actual_cp,
    results_df$actual_cp_raw,
    results_df$detected_cp_var_raw
  )
  results_df$detected_cp_var <- mapply(
    highlight_detected_cp_var,
    results_df$detected_cp_var_raw,
    results_df$actual_cp_raw
  )
  
  hl <- mapply(
    highlight_sites_advanced,
    results_df$sites,
    results_df$detected_sites,
    results_df$deleterious_sites,
    SIMPLIFY = FALSE
  )
  results_df$sites          <- vapply(hl, `[[`, character(1L), "sites")
  results_df$detected_sites <- vapply(hl, `[[`, character(1L), "detected")
  
  # --- 8. HTML table --------------------------------------------------------
  display_cols <- c(
    "scenario", "variant", "num_denovo_sequences",
    "sites", "detected_sites", "deleterious_sites",
    "num_months_ref_seq", "actual_month", "detection_month", "time_to_detection",
    "actual_partition", "detection_partition",
    "actual_cp", "detected_cp_var", "detected_cp_all"
  )
  
  html_table <- knitr::kable(
    results_df[, display_cols],
    format       = "html",
    escape       = FALSE,
    table.attr   = "class='table table-bordered'",
    caption      = caption_text
  ) %>%
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive")
    )
  
  if (isTRUE(save_html) && !is.null(output_file))
    kableExtra::save_kable(html_table, file = output_file)
  
  # --- 9. Return ------------------------------------------------------------
  list(
    Sim_List            = sim_res_list,
    Part_Data           = part_data_list,
    Results             = results_df,
    First_Detect_Sample = first_distinct_detect,
    CP_Accuracy         = cp_accuracy,
    Table               = html_table
  )
}