#' @title Detect Temporal Change Points (HDcpDetect)
#' @description Runs high-dimensional change point detection on time-series
#'   matrices using either Binary Segmentation or Wild Binary Segmentation.
#' @export
#'
#' @details
#' Iterates through time steps using an expanding or rolling window and applies
#' either \code{\link[HDcpDetect]{binary.segmentation}} or
#' \code{\link[HDcpDetect]{wild.binary.segmentation}}.
#'
#' Window behaviour mirrors \code{\link{detect_changepoints_ecp}}:
#' \itemize{
#'   \item \strong{Expanding} (\code{rolling_window = FALSE}): start index is
#'     fixed at \code{min_window}; end index grows by 1 each step.
#'   \item \strong{Rolling} (\code{rolling_window = TRUE}): both start and end
#'     advance by 1 each step, keeping the window length constant.
#' }
#'
#' @param data_matrix A numeric matrix (Time x Features). Rows are time points,
#'   columns are features (e.g. Hellinger distances per site).
#' @param min_window Integer. Starting row index of the initial window.
#' @param max_window Integer. Ending row index of the initial window.
#' @param n_timesteps Integer. Number of additional steps to iterate forward.
#' @param rolling_window Logical. If \code{TRUE}, the window slides. If
#'   \code{FALSE} (default), the window expands.
#' @param wild Logical. If \code{TRUE}, uses
#'   \code{\link[HDcpDetect]{wild.binary.segmentation}}. If \code{FALSE}
#'   (default), uses \code{\link[HDcpDetect]{binary.segmentation}}.
#' @param ... Additional arguments passed to the chosen segmentation function.
#'
#' @return A named list matching the original \code{generate_HDcp_list}
#'   structure:
#' \item{Points_List}{List of integer vectors \code{c(start, end)} giving the
#'   row index window used at each step.}
#' \item{HDcp_list}{List of full segmentation result objects, one per step.
#'   Steps where \code{end_idx} exceeds the matrix rows return \code{NULL}
#'   with a warning.}
#'
#' @importFrom HDcpDetect binary.segmentation wild.binary.segmentation
#'
#' @examples
#' set.seed(42)
#' baseline = matrix(rnorm(50, mean = 0, sd = 0.1), nrow = 10, ncol = 5)
#' variant  = matrix(rnorm(50, mean = 3, sd = 0.1), nrow = 10, ncol = 5)
#' data_mat = rbind(baseline, variant)
#'
#' res = detect_changepoints_hdcp(
#'   data_matrix = data_mat,
#'   min_window  = 1,
#'   max_window  = 20,
#'   n_timesteps = 0
#' )
#' print(res$HDcp_list[[1]])
detect_changepoints_hdcp = function(data_matrix,
                                    min_window,
                                    max_window,
                                    n_timesteps,
                                    rolling_window = FALSE,
                                    wild           = FALSE,
                                    ...) {
  
  # --- 1. Validation ---------------------------------------------------------
  if (!is.matrix(data_matrix)) stop("`data_matrix` must be a matrix.")
  if (nrow(data_matrix) < (max_window + n_timesteps))
    warning("`data_matrix` has fewer rows than the requested window extent.")
  
  # --- 2. Per-step worker ----------------------------------------------------
  extra_args = list(...)
  
  run_step = function(i) {
    start_idx = if (rolling_window) (min_window + i) else min_window
    end_idx   = max_window + i
    
    if (end_idx > nrow(data_matrix)) {
      warning("Step ", i, ": end_idx (", end_idx,
              ") exceeds matrix rows; skipping.")
      return(list(window = c(start_idx, end_idx),
                  hdcp   = NULL))
    }
    
    subset_mat = data_matrix[start_idx:end_idx, , drop = FALSE]
    
    hdcp_res = if (wild) {
      do.call(HDcpDetect::wild.binary.segmentation,
              c(list(subset_mat), extra_args))
    } else {
      do.call(HDcpDetect::binary.segmentation,
              c(list(subset_mat), extra_args))
    }
    
    list(window = c(start_idx, end_idx),
         hdcp   = hdcp_res)
  }
  
  # --- 3. Execute ------------------------------------------------------------
  results_list = lapply(0:n_timesteps, run_step)
  
  # --- 4. Return — names match original generate_HDcp_list ------------------
  list(
    Points_List = lapply(results_list, `[[`, "window"),
    HDcp_list   = lapply(results_list, `[[`, "hdcp")
  )
}