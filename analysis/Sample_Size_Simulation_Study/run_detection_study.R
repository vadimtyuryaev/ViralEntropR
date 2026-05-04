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
#' \strong{Memory note:} When \code{return_sims = FALSE} (the default),
#' simulation outputs and partition data are not retained after each
#' \code{n_new_seq} iteration, keeping peak memory proportional to the size
#' of a single simulation rather than \code{length(n_new_mut_seq_vec)}
#' simulations. Set \code{return_sims = TRUE} only when inspecting
#' intermediate results interactively (e.g. in a vignette), and only with
#' short \code{n_new_mut_seq_vec} vectors.
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
#' @param return_sims Logical. If \code{TRUE}, the returned list includes
#'   \code{Sim_List} and \code{Part_Data} (full simulation and partition
#'   objects for every \code{n_new_seq} value). If \code{FALSE} (default),
#'   these are \code{NULL}, keeping memory usage proportional to a single
#'   simulation. Use \code{TRUE} only for interactive inspection with short
#'   \code{n_new_mut_seq_vec} vectors.
#' @param ... Additional arguments passed to
#'   \code{\link{partition_time_windows}} (and on to
#'   \code{\link{cluster_sites_by_entropy}}).
#'
#' @return A named list:
#' \item{Sim_List}{Nested list \code{[[scenario]][[n_new_seq]]} of
#'   \code{viralSim} objects, or \code{NULL} when \code{return_sims = FALSE}.}
#' \item{Part_Data}{Nested list \code{[[scenario]][[n_new_seq]]} of
#'   \code{partition_time_windows} outputs, or \code{NULL} when
#'   \code{return_sims = FALSE}.}
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
#'   return_sims           = FALSE
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
    return_sims               = FALSE,
    ...
) {
  
  # --- 0. Setup --------------------------------------------------------------
  n_scen                <- length(variants_list)
  first_distinct_detect <- rep(NA_real_, n_scen)
  cp_accuracy           <- rep(NA_real_, n_scen)
  all_rows              <- vector("list", n_scen)
  
  # Sim_List and Part_Data are only populated when return_sims = TRUE.
  # When FALSE (default / simulation study), they remain NULL, keeping
  # peak memory proportional to one simulation at a time rather than
  # length(n_new_mut_seq_vec) simulations simultaneously.
  sim_res_list   <- if (isTRUE(return_sims)) vector("list", n_scen) else NULL
  part_data_list <- if (isTRUE(return_sims)) vector("list", n_scen) else NULL
  
  n_col      <- nchar(ref_seq)
  sim_dates  <- seq.Date(as.Date(start_date), as.Date(end_date), by = "month")
  total_seqs <- (length(sim_dates) - num_months_ref_seq) * n_seq_per_month
  n_ref      <- num_months_ref_seq * n_seq_per_month
  
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
    list(sites       = paste(ah, collapse = ","),
         detected    = paste(dh, collapse = ","),
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
  
  # --- 3. Extra args ---------------------------------------------------------
  extra_args <- list(...)
  
  # --- 4. Scenario loop ------------------------------------------------------
  for (scenario_idx in seq_along(variants_list)) {
    
    number_of_variants <- variants_list[[scenario_idx]]
    mutation_rates     <- mutation_rate_list[[scenario_idx]]
    nv                 <- length(number_of_variants)
    vi                 <- variant_interval
    if (length(vi) < nv - 1L)
      vi <- c(vi, rep(utils::tail(vi, 1L), nv - 1L - length(vi)))
    
    # Pre-allocate accumulators for this scenario
    rows_this_scenario  <- vector("list", length(n_new_mut_seq_vec))
    fully_detected_vec  <- logical(length(n_new_mut_seq_vec))
    sites_list_per_n    <- vector("list", length(n_new_mut_seq_vec))
    
    if (isTRUE(return_sims)) {
      sim_res_list[[scenario_idx]]   <- vector("list", length(n_new_mut_seq_vec))
      part_data_list[[scenario_idx]] <- vector("list", length(n_new_mut_seq_vec))
    }
    
    # --- 4a. n_new_seq loop (sequential) -------------------------------------
    for (k in seq_along(n_new_mut_seq_vec)) {
      
      n_new_seq <- n_new_mut_seq_vec[k]
      
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
      mat_sim    <- sim_res$Simulation_Output
      enc_mat    <- encode_aa_sequence(as.matrix(mat_sim[, seq_len(n_col)]))
      AL_df      <- as.data.frame(enc_mat)
      AL_df[]    <- lapply(AL_df, as.integer)
      AL_df$Date <- as.Date(format(as.Date(mat_sim$Date), "%Y-%m-01"))
      
      part_data <- do.call(partition_time_windows,
                           c(list(data          = AL_df,
                                  n_sites       = n_col,
                                  window_length = sliding_window_length,
                                  window_type   = window_option,
                                  start_date    = start_date,
                                  end_date      = end_date),
                             extra_args))
      
      # 3. Relabel all cluster DataFrames once upfront ------------------------
      # Class 1 is the highest-entropy group after relabeling. Applied here,
      # once per part_data, before any class access in the variant loop.
      relabeled_clusters <- lapply(part_data$Clusters, function(cl) {
        cl$DataFrame <- relabel_entropy_classes(cl$DataFrame)
        cl
      })
      
      # 4. Hellinger matrix ---------------------------------------------------
      hell_mat <- calculate_hellinger_matrix(
        partitions = part_data$Partitions,
        sites      = seq_len(n_col),
        aa_levels  = 25L
      )
      dat_t <- t(hell_mat)
      
      # 5. Global change point detection (full timeline) ----------------------
      max_cp <- length(part_data$Clusters) - 1L
      if (max_cp >= 2L) {
        cp_all_raw <- ecp::e.agglo(
          X       = dat_t,
          member  = seq_len(nrow(dat_t)),
          alpha   = 1,
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
        
        # Guard: part_em can exceed available partitions when end_date equals
        # (rather than one month beyond) the last simulation month.
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
        
        # Forward search: find the first partition where all non-deleterious
        # mutation sites appear in the top-entropy cluster.
        # An explicit range check guards against the seq(n+1, n) pitfall in R
        # which returns a descending sequence rather than an empty one.
        dp <- part_em
        if (!all(setdiff(sites, deleterious_site) %in% top_sites)) {
          dp         <- NA
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
        
        # Per-variant CP: e.agglo on Hellinger sequence up to part_em - 1
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
      } # end variant loop
      
      # Accumulate lightweight results rows
      rows_this_scenario[[k]] <- do.call(rbind, rows_this)
      fully_detected_vec[k]   <- all(detected_flags_all)
      sites_list_per_n[[k]]   <- sites_list_all
      
      # Store heavy objects only when explicitly requested
      if (isTRUE(return_sims)) {
        sim_res_list[[scenario_idx]][[k]]   <- sim_res
        part_data_list[[scenario_idx]][[k]] <- part_data
      }
      
      # Free heavy objects immediately when not needed — this keeps peak
      # memory proportional to one simulation rather than n_new_mut_seq_vec
      # simulations, which is the root cause of the 30 GB accumulation.
      rm(sim_res, mat_sim, enc_mat, AL_df, part_data,
         relabeled_clusters, hell_mat, dat_t, rows_this,
         details, delet)
      
    } # end n_new_seq loop
    
    # Combine rows for this scenario
    all_rows[[scenario_idx]] <- do.call(rbind, rows_this_scenario)
    
    # First distinct detection threshold:
    # First n_new_seq where all variants detected AND site sets are mutually
    # distinct (no site attributed to more than one variant simultaneously).
    for (k in seq_along(n_new_mut_seq_vec)) {
      if (fully_detected_vec[k]) {
        vec <- unlist(strsplit(unlist(sites_list_per_n[[k]]), ","))
        if (length(unique(vec)) == length(vec) &&
            is.na(first_distinct_detect[scenario_idx])) {
          first_distinct_detect[scenario_idx] <- n_new_mut_seq_vec[k]
        }
      }
    }
    
    # Name the stored lists by n_new_seq value for easy access
    if (isTRUE(return_sims)) {
      names(sim_res_list[[scenario_idx]])   <- as.character(n_new_mut_seq_vec)
      names(part_data_list[[scenario_idx]]) <- as.character(n_new_mut_seq_vec)
    }
    
  } # end scenario loop
  
  # --- 5. Combine all rows ---------------------------------------------------
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
  
  # --- 7. Post-processing: sort site strings and apply HTML highlighting ----
  results_df$sites          <- vapply(results_df$sites_raw,
                                      sort_string, character(1L))
  results_df$detected_sites <- vapply(results_df$detected_sites_raw,
                                      sort_string, character(1L))
  
  results_df$actual_cp <- mapply(highlight_actual_cp,
                                 results_df$actual_cp_raw,
                                 results_df$detected_cp_var_raw)
  results_df$detected_cp_var <- mapply(highlight_detected_cp_var,
                                       results_df$detected_cp_var_raw,
                                       results_df$actual_cp_raw)
  
  hl <- mapply(highlight_sites_advanced,
               results_df$sites,
               results_df$detected_sites,
               results_df$deleterious_sites,
               SIMPLIFY = FALSE)
  results_df$sites          <- vapply(hl, `[[`, character(1L), "sites")
  results_df$detected_sites <- vapply(hl, `[[`, character(1L), "detected")
  
  # --- 8. HTML table ---------------------------------------------------------
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
  
  # --- 9. Return -------------------------------------------------------------
  list(
    Sim_List            = sim_res_list,
    Part_Data           = part_data_list,
    Results             = results_df,
    First_Detect_Sample = first_distinct_detect,
    CP_Accuracy         = cp_accuracy,
    Table               = html_table
  )
}