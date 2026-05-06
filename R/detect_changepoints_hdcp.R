#' @title Detect Temporal Change Points (HDcpDetect)
#' @description Runs high-dimensional change point detection on time-series
#'   matrices (typically Hellinger distance matrices) over one or more time
#'   windows, using either binary segmentation or wild binary segmentation.
#' @export
#'
#' @details
#' The function applies \code{\link[HDcpDetect]{binary.segmentation}} or
#' \code{\link[HDcpDetect]{wild.binary.segmentation}} repeatedly to slices of
#' \code{data_matrix}, advancing the slice forward by one row at each step.
#' A total of \code{n_timesteps + 1} detections are performed: one on the
#' initial window \code{[min_window, max_window]}, then \code{n_timesteps}
#' additional detections on subsequent windows.
#'
#' Two window-advancement modes are supported, mirroring
#' \code{\link{detect_changepoints_ecp}}:
#' \itemize{
#'   \item \strong{Expanding} (\code{rolling_window = FALSE}, default): the
#'     start index is fixed at \code{min_window}; the end index grows by 1
#'     each step. Each iteration uses one more row than the previous.
#'     Natural for online surveillance, where change-point detection is
#'     re-run as new data accumulates.
#'   \item \strong{Rolling} (\code{rolling_window = TRUE}): both start and
#'     end indices advance by 1 each step, keeping the window length
#'     constant. Natural for retrospective windowed analysis.
#' }
#'
#' \strong{Choice of segmentation method.} Binary segmentation
#' (\code{wild = FALSE}, default) is the classical recursive method: it
#' finds the most likely change point, splits the series, and recurses.
#' Wild binary segmentation (\code{wild = TRUE}) is the variant of
#' Fryzlewicz (2014) that draws random subintervals to improve detection
#' of multiple closely-spaced change points; it accepts an additional
#' \code{M} argument controlling the number of random intervals.
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
#' @param wild Logical. If \code{TRUE}, uses
#'   \code{\link[HDcpDetect]{wild.binary.segmentation}}; if \code{FALSE}
#'   (default), uses \code{\link[HDcpDetect]{binary.segmentation}}. See
#'   Details for the methodological distinction.
#' @param ... Additional arguments passed to the chosen segmentation
#'   function. Note that \code{M} (number of random intervals) is only
#'   accepted by \code{wild.binary.segmentation}; passing it with
#'   \code{wild = FALSE} produces an "unused argument" error.
#'
#' @return A named list:
#' \item{Points_List}{List of integer vectors \code{c(start, end)} giving
#'   the row index window used at each step.}
#' \item{HDcp_list}{List of full segmentation result objects, one per
#'   step. Steps where \code{end_idx} exceeds the matrix rows return
#'   \code{NULL} with a warning.}
#'
#' @section Warnings and errors:
#' \describe{
#'   \item{\code{stop()}}{Triggered if \code{max_window} exceeds
#'     \code{nrow(data_matrix)}: the initial window cannot be
#'     constructed and no detection can run.}
#'   \item{\code{warning()}}{Triggered (and execution continues) if
#'     \code{max_window + n_timesteps} exceeds \code{nrow(data_matrix)},
#'     in which case the offending later iterations are skipped and
#'     \code{NULL} entries appear in \code{HDcp_list}.}
#' }
#'
#' @seealso \code{\link{detect_changepoints_ecp}} for the energy-statistic
#'   alternative; \code{\link{calculate_hellinger_matrix}} for the typical
#'   upstream input.
#'
#' @importFrom HDcpDetect binary.segmentation wild.binary.segmentation
#'
#' @examples
#' set.seed(42)
#' baseline = matrix(rnorm(50, mean = 0, sd = 0.1), nrow = 10, ncol = 5)
#' variant  = matrix(rnorm(50, mean = 3, sd = 0.1), nrow = 10, ncol = 5)
#' data_mat = rbind(baseline, variant)
#'
#' # Single-window detection with the default binary segmentation.
#' res = detect_changepoints_hdcp(
#'   data_matrix = data_mat,
#'   min_window  = 1,
#'   max_window  = 20,
#'   n_timesteps = 0
#' )
#' print(res$HDcp_list[[1]])
#'
#' # Wild binary segmentation with a custom number of random intervals.
#' res_wild = detect_changepoints_hdcp(
#'   data_matrix = data_mat,
#'   min_window  = 1,
#'   max_window  = 20,
#'   n_timesteps = 0,
#'   wild        = TRUE,
#'   M           = 100
#' )
#' print(res_wild$HDcp_list[[1]])
detect_changepoints_hdcp = function(data_matrix,
                                    min_window,
                                    max_window,
                                    n_timesteps,
                                    rolling_window = FALSE,
                                    wild           = FALSE,
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
  
  # --- 2. Per-step worker ----------------------------------------------------
  
  run_step = function(i) {
    start_idx = if (rolling_window) (min_window + i) else min_window
    end_idx   = max_window + i
    
    if (end_idx > n_rows) {
      warning("Step ", i, ": end_idx (", end_idx,
              ") exceeds matrix rows (", n_rows, "); skipping.",
              call. = FALSE)
      return(list(window = c(start_idx, end_idx),
                  hdcp   = NULL))
    }
    
    subset_mat = data_matrix[start_idx:end_idx, , drop = FALSE]
    
    hdcp_res = tryCatch({
      utils::capture.output({
        res_inner = if (wild) {
          do.call(HDcpDetect::wild.binary.segmentation,
                  c(list(subset_mat), extra_args))
        } else {
          do.call(HDcpDetect::binary.segmentation,
                  c(list(subset_mat), extra_args))
        }
      })
      res_inner
    }, error = function(e) {
      warning("Step ", i, ": HDcpDetect call failed: ",
              conditionMessage(e), "; returning NULL.",
              call. = FALSE)
      NULL
    })
    
    list(window = c(start_idx, end_idx),
         hdcp   = hdcp_res)
  }
  
  # --- 3. Execute ------------------------------------------------------------
  results_list = lapply(0:n_timesteps, run_step)
  
  # --- 4. Return -------------------------------------------------------------
  list(
    Points_List = lapply(results_list, `[[`, "window"),
    HDcp_list   = lapply(results_list, `[[`, "hdcp")
  )
}