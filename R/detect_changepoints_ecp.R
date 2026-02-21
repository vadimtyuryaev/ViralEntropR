#' @title Detect Temporal Change Points (ECP)
#' @description Runs Energy Change Point detection on time-series matrices (typically distance matrices).
#' @export
#'
#' @details
#' This function iterates through time steps (expanding or rolling window) and applies
#' the \code{\link[ecp]{ks.cp3o}} algorithm to detect distributional changes in the data matrix.
#'
#' @param data_matrix A numeric matrix (Time x Features). Usually a transposed Hellinger distance matrix.
#' @param min_window Integer. The starting window size (lower bound index).
#' @param max_window Integer. The ending window size offset.
#' @param n_timesteps Integer. Number of steps to iterate forward.
#' @param rolling_window Logical. If \code{TRUE}, both the start and end of the window move forward (sliding).
#' If \code{FALSE} (default), the start remains fixed and the window expands.
#' @param dynamic_k Logical. If \code{TRUE}, sets the maximum number of change points \code{K}
#' to \code{window_length - 1} dynamically for each iteration.
#' @param n_cores Integer. Number of cores for parallel execution. Default is
#'   \code{1L} (sequential). Uses \code{mclapply} on Unix/macOS and
#'   \code{parLapply} on Windows.
#' @param ... Additional arguments passed to \code{\link[ecp]{ks.cp3o}}.
#'
#' @return A named list matching the original \code{generate_ECP_list_modified}
#'   structure:
#' \item{Points_List}{List of integer vectors \code{c(start, end)} giving the
#'   row index window used at each step.}
#' \item{ECP_list}{List of full \code{ks.cp3o} result objects, one per step.}
#' \item{ECP_est_list}{List of change point estimate vectors, one per step.}
#'
#' @importFrom ecp ks.cp3o
#' @importFrom parallel mclapply parLapply makeCluster stopCluster clusterExport
#'
#' @examples
#' set.seed(123)
#' baseline = matrix(rnorm(50, mean = 0, sd = 0.1), nrow = 10, ncol = 5)
#' variant  = matrix(rnorm(50, mean = 3, sd = 0.1), nrow = 10, ncol = 5)
#' data_mat = rbind(baseline, variant)
#'
#' # Single step over the full window; expect a change point near row 11
#' res = detect_changepoints_ecp(
#'   data_matrix = data_mat,
#'   min_window  = 1,
#'   max_window  = 20,
#'   n_timesteps = 0,
#'   minsize     = 2
#' )
#' print(res$ECP_est_list[[1]])
detect_changepoints_ecp = function(data_matrix,
                                    min_window,
                                    max_window,
                                    n_timesteps,
                                    rolling_window = FALSE,
                                    dynamic_k      = FALSE,
                                    n_cores        = 1L,
                                    ...) {
  
  # --- 1. Validation ---------------------------------------------------------
  if (!is.matrix(data_matrix)) stop("`data_matrix` must be a matrix.")
  if (nrow(data_matrix) < (max_window + n_timesteps)) {
    warning("`data_matrix` has fewer rows than the requested window extent.")
  }
  
  # --- 2. Per-step worker ----------------------------------------------------
  # Rolling: both start and end advance by 1 each step (j tracks start offset,
  # matching the original loop where j increments after each iteration).
  # Expanding (non-rolling): start is fixed, end grows by 1 each step.
  #
  # dynamic_k mirrors the original adj_K formula:
  #   K = (max_window - min_window - 1) + i  =  nrow(subset) - 2
  # This is deliberately one less than nrow-1 to stay conservative.
  run_step = function(i) {
    start_idx = if (rolling_window) (min_window + i) else min_window
    end_idx   = max_window + i
    
    if (end_idx > nrow(data_matrix)) {
      warning("Step ", i, ": end_idx (", end_idx,
              ") exceeds matrix rows; skipping.")
      return(list(window = c(start_idx, end_idx),
                  estimate = NULL,
                  full     = NULL))
    }
    
    subset_mat = data_matrix[start_idx:end_idx, , drop = FALSE]
    
    # Build ks.cp3o argument list; only inject K when dynamic_k is TRUE,
    # preserving the ability to pass K manually via ... when it is FALSE.
    args = list(Z = subset_mat, ...)
    if (dynamic_k) {
      # Original formula: upper - lower - 1 + i = nrow(subset) - 2
      args$K = nrow(subset_mat) - 2L
    }
    
    ecp_res = do.call(ecp::ks.cp3o, args)
    
    list(
      window   = c(start_idx, end_idx),
      estimate = ecp_res$estimates,
      full     = ecp_res
    )
  }
  
  # --- 3. Execute: parallel or sequential ------------------------------------
  indices      = 0:n_timesteps
  use_parallel = n_cores > 1L && length(indices) > 1L
  
  if (use_parallel) {
    if (.Platform$OS.type == "windows") {
      cl = parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(
        cl,
        varlist = c("data_matrix", "min_window", "max_window",
                    "rolling_window", "dynamic_k"),
        envir = environment()
      )
      results_list = parallel::parLapply(cl, indices, run_step)
    } else {
      results_list = parallel::mclapply(indices, run_step,
                                         mc.cores = n_cores)
    }
  } else {
    results_list = lapply(indices, run_step)
  }
  
  # --- 4. Return — names match original generate_ECP_list_modified -----------
  list(
    Points_List  = lapply(results_list, `[[`, "window"),
    ECP_list     = lapply(results_list, `[[`, "full"),
    ECP_est_list = lapply(results_list, `[[`, "estimate")
  )
}