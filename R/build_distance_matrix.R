#' @title Build Distance/Entropy Matrix
#' @description Converts the output of \code{\link{calculate_hellinger_matrix}}
#'   or the legacy \code{get_hellinger_dist} / entropy list into a numeric
#'   \code{sites × time_steps} matrix suitable for change point detection.
#'
#' @details
#' In the current pipeline \code{\link{calculate_hellinger_matrix}} already
#' returns a \code{sites × time_steps} matrix, so calling this function is
#' usually unnecessary — transposing that matrix directly gives \code{dat_t}
#' for \code{\link[ecp]{e.agglo}} or \code{\link[ecp]{ks.cp3o}}:
#'
#' \preformatted{
#'   hell_mat = calculate_hellinger_matrix(partitions, sites = seq_len(n_sites))
#'   dat_t    = t(hell_mat)   # time_steps × sites — ready for ECP
#' }
#'
#' This function is retained for \strong{backward compatibility} with code that
#' was written against the old \code{get_hellinger_dist} list format, and as a
#' convenience wrapper that also accepts entropy lists.
#'
#' \strong{Input formats accepted:}
#' \enumerate{
#'   \item A numeric matrix (e.g. from \code{calculate_hellinger_matrix}) —
#'     returned as-is with site rownames normalised.
#'   \item The legacy named list from \code{get_hellinger_dist} containing
#'     \code{$Sites} and \code{$Hellinger_Distances} (list of vectors).
#'   \item An entropy list containing \code{$Sites} and \code{$Entropies}
#'     (list of numeric vectors), used when no Hellinger distances are present.
#' }
#'
#' @param data_input Either:
#'   \itemize{
#'     \item A numeric matrix (\code{sites × time_steps}) as returned by
#'       \code{\link{calculate_hellinger_matrix}}.
#'     \item A named list with elements \code{Sites} and one of
#'       \code{Hellinger_Distances} or \code{Entropies}, as returned by the
#'       legacy \code{get_hellinger_dist} function.
#'   }
#'
#' @return A numeric matrix with:
#'   \item{Rows}{Sites (1 to \code{max(sites)}, sparse if sites are
#'     non-contiguous).}
#'   \item{Columns}{Time steps (partitions T2, T3, … or entropy periods).}
#'   Row names are character site indices. Column names are preserved from
#'   the input where available.
#'
#' @seealso \code{\link{calculate_hellinger_matrix}},
#'   \code{\link{detect_changepoints_ecp}}
#'
#' @keywords internal
build_distance_matrix = function(data_input) {
  
  # --- 1. Matrix input: already in the right format -------------------------
  # Returned by calculate_hellinger_matrix — just normalise rownames and return
  if (is.matrix(data_input)) {
    if (is.null(rownames(data_input))) {
      rownames(data_input) = as.character(seq_len(nrow(data_input)))
    }
    return(data_input)
  }
  
  # --- 2. Legacy list input -------------------------------------------------
  if (!is.list(data_input)) {
    stop("`data_input` must be a matrix or a named list (see ?build_distance_matrix).")
  }
  
  # Determine data source: Hellinger distances take priority over Entropies
  has_hellinger = length(data_input$Hellinger_Distances) > 0
  dat_temp      = if (has_hellinger) data_input$Hellinger_Distances
  else               data_input$Entropies
  
  if (is.null(dat_temp) || length(dat_temp) == 0) {
    stop("Input list must contain non-empty `Hellinger_Distances` or `Entropies`.")
  }
  
  sites    = sort(unique(as.numeric(as.character(data_input$Sites))))
  n_cols   = length(dat_temp[[1L]])            # time steps per site
  data_mat = matrix(0, nrow = max(sites), ncol = n_cols)
  rownames(data_mat) = as.character(seq_len(max(sites)))
  
  # Preserve column names from the first element if available
  if (!is.null(names(dat_temp[[1L]]))) {
    colnames(data_mat) = names(dat_temp[[1L]])
  }
  
  # Fill rows by site index (sparse-safe: zeros remain for absent sites)
  for (i in seq_along(sites)) {
    loc              = sites[i]
    match_idx        = which(data_input$Sites == loc)
    data_mat[loc, ]  = dat_temp[[match_idx]]
  }
  
  data_mat
}