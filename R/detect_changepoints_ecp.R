#' @title Detect Temporal Change Points (ECP)
#' @description Runs Energy Change Point detection on time-series matrices
#'   (typically Hellinger distance matrices).
#' @export
#'
#' @details
#' Iterates through time steps using an expanding or rolling window and applies
#' the \code{\link[ecp]{ks.cp3o}} algorithm to detect distributional changes.
#'
#' @param data_matrix A numeric matrix (Time x Features). Usually a transposed
#'   Hellinger distance matrix.
#' @param min_window Integer. Starting row index of the initial window.
#' @param max_window Integer. Ending row index of the initial window.
#' @param n_timesteps Integer. Number of additional steps to iterate forward.
#' @param rolling_window Logical. If \code{TRUE}, both start and end advance by
#'   1 each step (sliding). If \code{FALSE} (default), start is fixed and end
#'   expands.
#' @param dynamic_k Logical. If \code{TRUE}, sets \code{K} to
#'   \code{nrow(subset) - 2} dynamically for each step.
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
#'
#' @examples
#' set.seed(123)
#' baseline = matrix(rnorm(50, mean = 0, sd = 0.1), nrow = 10, ncol = 5)
#' variant  = matrix(rnorm(50, mean = 3, sd = 0.1), nrow = 10, ncol = 5)
#' data_mat = rbind(baseline, variant)
#'
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
                                   ...) {
  
  # --- 1. Validation ---------------------------------------------------------
  if (!is.matrix(data_matrix)) stop("`data_matrix` must be a matrix.")
  if (nrow(data_matrix) < (max_window + n_timesteps))
    warning("`data_matrix` has fewer rows than the requested window extent.")
  
  # --- 2. Per-step worker ----------------------------------------------------
  # Capture ... once so it does not need to be re-evaluated on every iteration.
  # dynamic_k mirrors the original adj_K formula: nrow(subset) - 2.
  extra_args = list(...)
  
  run_step = function(i) {
    start_idx = if (rolling_window) (min_window + i) else min_window
    end_idx   = max_window + i
    
    if (end_idx > nrow(data_matrix)) {
      warning("Step ", i, ": end_idx (", end_idx,
              ") exceeds matrix rows; skipping.")
      return(list(window   = c(start_idx, end_idx),
                  estimate = NULL,
                  full     = NULL))
    }
    
    subset_mat = data_matrix[start_idx:end_idx, , drop = FALSE]
    
    args = c(list(Z = subset_mat), extra_args)
    if (dynamic_k) args$K = nrow(subset_mat) - 2L
    
    ecp_res = do.call(ecp::ks.cp3o, args)
    
    list(window   = c(start_idx, end_idx),
         estimate = ecp_res$estimates,
         full     = ecp_res)
  }
  
  # --- 3. Execute ------------------------------------------------------------
  results_list = lapply(0:n_timesteps, run_step)
  
  # --- 4. Return — names match original generate_ECP_list_modified -----------
  list(
    Points_List  = lapply(results_list, `[[`, "window"),
    ECP_list     = lapply(results_list, `[[`, "full"),
    ECP_est_list = lapply(results_list, `[[`, "estimate")
  )
}