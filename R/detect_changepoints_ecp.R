#' @title Detect Temporal Change Points (ECP)
#' @description Runs Energy Change Point detection on time-series matrices
#'   (typically Hellinger distance matrices) over one or more time windows.
#' @export
#'
#' @details
#' The function applies \code{\link[ecp]{ks.cp3o}} repeatedly to slices of
#' \code{data_matrix}, advancing the slice forward by one row at each step.
#' A total of \code{n_timesteps + 1} detections are performed: one on the
#' initial window \code{[min_window, max_window]}, then \code{n_timesteps}
#' additional detections on subsequent windows.
#'
#' Two window-advancement modes are supported:
#' \itemize{
#'   \item \strong{Expanding} (\code{rolling_window = FALSE}, default): the
#'     start index is fixed at \code{min_window}; the end index grows by 1
#'     each step. Each iteration uses one more row than the previous.
#'     Natural for online surveillance, where change-point detection is
#'     re-run as new data accumulates.
#'   \item \strong{Rolling} (\code{rolling_window = TRUE}): both start and
#'     end indices advance by 1 each step, keeping the window length
#'     constant. Natural for retrospective windowed analysis, where local
#'     change-point structure is examined in fixed-size segments.
#' }
#'
#' @param data_matrix A numeric matrix (Time x Features). Rows are time
#'   points, columns are features (e.g. Hellinger distances per site).
#' @param min_window Integer. Starting row index of the initial window.
#' @param max_window Integer. Ending row index of the initial window. Must
#'   not exceed \code{nrow(data_matrix)}.
#' @param n_timesteps Integer >= 0. Number of \emph{additional} detections
#'   to run after the initial window. The function performs
#'   \code{n_timesteps + 1} detections in total. Set to \code{0} for a
#'   single detection on the initial window.
#' @param rolling_window Logical. If \code{TRUE}, both window endpoints
#'   advance per step (sliding window of constant length). If \code{FALSE}
#'   (default), only the end advances (expanding window).
#' @param dynamic_k Logical. If \code{TRUE}, sets \code{K} to
#'   \code{nrow(subset) - 2L} at each step, so the algorithm considers
#'   the maximum number of change points permitted by the window length
#'   (subject to \code{minsize}). Useful for exploratory analysis where
#'   the true number of change points is unknown a priori. Most idiomatic
#'   when combined with \code{n_timesteps > 0}, where the maximum
#'   admissible K legitimately scales with window size, but also valid
#'   for single-window exploratory use. If \code{FALSE} (default), the
#'   value of \code{K} supplied via \code{...} is honoured at every step;
#'   if no \code{K} is supplied, the \code{ks.cp3o} default of 1 is used.
#' @param ... Additional arguments passed to \code{\link[ecp]{ks.cp3o}}.
#'   Common choices include \code{K} (maximum number of change points,
#'   default 1) and \code{minsize} (minimum segment size, default 30).
#'   Note that \code{K} is overridden when \code{dynamic_k = TRUE}; a
#'   warning is emitted in that case.
#'
#' @return A named list:
#' \item{Points_List}{List of integer vectors \code{c(start, end)} giving
#'   the row index window used at each step.}
#' \item{ECP_list}{List of full \code{ks.cp3o} result objects, one per
#'   step. Inspect \code{$cpLoc} for the optimal change-point locations
#'   at each candidate count, and \code{$gofM} for the goodness-of-fit
#'   curve, to apply post-hoc filtering.}
#' \item{ECP_est_list}{List of change-point estimate vectors (the
#'   algorithm's own selection, equivalent to \code{$estimates}).}
#'
#' @section Warnings and errors:
#' \describe{
#'   \item{\code{stop()}}{Triggered if \code{max_window} exceeds
#'     \code{nrow(data_matrix)}: the initial window cannot be
#'     constructed and no detection can run.}
#'   \item{\code{warning()}}{Triggered (and execution continues) if
#'     (a) \code{max_window + n_timesteps} exceeds
#'     \code{nrow(data_matrix)}, in which case the offending later
#'     iterations are skipped and \code{NULL} entries appear in the
#'     result lists; (b) \code{dynamic_k = TRUE} and \code{K} is also
#'     supplied via \code{...}, in which case the user-supplied
#'     \code{K} is silently overridden by \code{nrow(subset) - 2L}.}
#' }
#'
#' @importFrom ecp ks.cp3o
#'
#' @examples
#' set.seed(123)
#' baseline = matrix(rnorm(50, mean = 0, sd = 0.1), nrow = 10, ncol = 5)
#' variant  = matrix(rnorm(50, mean = 3, sd = 0.1), nrow = 10, ncol = 5)
#' data_mat = rbind(baseline, variant)
#'
#' # Single-window detection: one true change point at row 11.
#' res = detect_changepoints_ecp(
#'   data_matrix = data_mat,
#'   min_window  = 1,
#'   max_window  = 20,
#'   n_timesteps = 0,
#'   minsize     = 5
#' )
#' print(res$ECP_est_list[[1]])
#'
#' # Strongest single change point: inspect $cpLoc[[1]] in the full result.
#' print(res$ECP_list[[1]]$cpLoc[[1]])
detect_changepoints_ecp = function(data_matrix,
                                   min_window,
                                   max_window,
                                   n_timesteps,
                                   rolling_window = FALSE,
                                   dynamic_k      = FALSE,
                                   ...) {
  
  # --- 1. Validation ---------------------------------------------------------
  if (!is.matrix(data_matrix))
    stop("`data_matrix` must be a matrix.", call. = FALSE)
  
  if (!is.numeric(n_timesteps) || length(n_timesteps) != 1L || n_timesteps < 0)
    stop("`n_timesteps` must be a single non-negative integer.",
         call. = FALSE)
  
  n_rows = nrow(data_matrix)
  
  # Hard guard: if the initial window itself is impossible, no detection
  # can run. Returning a list of NULLs would be silently degraded output;
  # an explicit error forces the caller to fix the call.
  if (max_window > n_rows)
    stop("`max_window` (", max_window, ") exceeds nrow(data_matrix) (",
         n_rows, "). The initial window [", min_window, ", ",
         max_window, "] cannot be constructed.",
         call. = FALSE)
  
  # Soft guard: if later iterations would go past the matrix, warn with
  # concrete numbers so the caller can decide whether to fix or accept
  # the degraded output. Skipped iterations return NULL entries.
  needed = max_window + n_timesteps
  if (n_rows < needed) {
    n_skipped = needed - n_rows
    warning(sprintf(
      "`data_matrix` has %d rows but max_window (%d) + n_timesteps (%d) requires %d. %d step(s) will be skipped and return NULL.",
      n_rows, max_window, n_timesteps, needed, n_skipped),
      call. = FALSE)
  }
  
  # Capture ... once so it does not need to be re-evaluated on every iteration.
  extra_args = list(...)
  
  # Guard against silent K override: if the caller supplied K via ... and
  # also set dynamic_k = TRUE, the dynamic-K branch overrides the
  # user-supplied K. Warn once at function entry rather than once per window.
  if (isTRUE(dynamic_k) && "K" %in% names(extra_args)) {
    warning(
      "Argument `K` supplied via `...` is ignored because ",
      "`dynamic_k = TRUE`; K is set to nrow(subset) - 2L for each window. ",
      "Set `dynamic_k = FALSE` to honour the user-supplied K.",
      call. = FALSE
    )
  }
  
  # --- 2. Per-step worker ----------------------------------------------------
  # dynamic_k mirrors the original adj_K formula: nrow(subset) - 2.
  
  run_step = function(i) {
    start_idx = if (rolling_window) (min_window + i) else min_window
    end_idx   = max_window + i
    
    if (end_idx > n_rows) {
      warning("Step ", i, ": end_idx (", end_idx,
              ") exceeds matrix rows (", n_rows, "); skipping.",
              call. = FALSE)
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