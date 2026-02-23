#' @title Calculate Hellinger Distance Matrix
#' @description Computes the Hellinger distance between the amino acid distribution
#'   of a reference time point (first partition) and all subsequent time points,
#'   for each requested site.
#'
#' @details
#' The Hellinger distance between two discrete distributions \eqn{P} and \eqn{Q} is:
#' \deqn{H(P, Q) = \sqrt{\sum_{i=1}^{k} (\sqrt{p_i} - \sqrt{q_i})^2}}
#' When \code{normalized = TRUE} the result is scaled by \eqn{1/\sqrt{2}}, bounding
#' the distance to \eqn{[0, 1]}. Otherwise the range is \eqn{[0, \sqrt{2}]}.
#'
#' Internally, amino acid counts per site per partition are tabulated using
#' \code{\link{get_site_counts}}, which uses \code{\link[base]{tabulate}} for fast
#' integer counting. Column-wise proportions are then computed with
#' \code{\link[base]{sweep}} and distances with \code{\link[base]{colSums}} —
#' fully vectorized over partitions for each site.
#'
#' @param partitions A list of data frames (one per time window). Each data frame
#'   must have numeric-encoded amino acid sequences as columns (integers 1 to
#'   \code{aa_levels}).
#' @param sites Integer vector. Indices of the sites to analyse. Defaults to all
#'   sites (\code{seq_len(ncol(partitions[[1]]))}).
#' @param aa_levels Integer. Alphabet size (e.g. 25 for standard amino acid
#'   encoding). Default is 25.
#' @param normalized Logical. If \code{TRUE}, scales distances by
#'   \eqn{1/\sqrt{2}} to bound the result in \eqn{[0, 1]}. Default is
#'   \code{FALSE}.
#' @param labels Character vector of partition labels. Length must equal
#'   \code{length(partitions)}. Defaults to \code{"T1"}, \code{"T2"}, …
#' @param save_freq_tables Logical. If \code{TRUE}, the return value also
#'   includes raw count tables and proportion tables for each site. Default is
#'   \code{FALSE}.
#' @param n_cores Integer. Number of cores for parallel processing of sites via
#'   \code{future.apply::future_lapply}. Default is \code{1L} (sequential).
#'   Sets \code{future::multisession} internally when \code{n_cores > 1} and
#'   restores \code{future::sequential} on exit.
#'
#' @return When \code{save_freq_tables = FALSE} (default): a numeric matrix with
#'   rows = sites and columns = time steps T2, T3, … (distances relative to T1).
#'   Row names are the site indices; column names are taken from
#'   \code{labels[-1]}.
#'
#'   When \code{save_freq_tables = TRUE}: a named list with elements
#'   \code{Sites}, \code{Hellinger_Distances}, \code{Frequency_Tables}, and
#'   \code{Proportions_Tables}.
#'
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#'
#' @examples
#' p1 = data.frame(s1 = c(1, 1, 1, 1, 1), s2 = c(20, 20, 20, 20, 20))
#' p2 = data.frame(s1 = c(20, 20, 20, 20, 20), s2 = c(20, 20, 20, 20, 20))
#' p3 = data.frame(s1 = c(1, 1, 20, 20, 20), s2 = c(20, 20, 20, 20, 20))
#' parts = list(T1 = p1, T2 = p2, T3 = p3)
#'
#' # All sites, unnormalized
#' result = calculate_hellinger_matrix(parts, sites = 1:2)
#' print(result)
#'
#' # With freq tables saved
#' result2 = calculate_hellinger_matrix(parts, sites = 1:2,
#'                                       save_freq_tables = TRUE)
#' print(result2$Hellinger_Distances)
#' @export
calculate_hellinger_matrix = function(partitions,
                                      sites       = seq_len(ncol(partitions[[1]])),
                                      aa_levels   = 25L,
                                      normalized  = FALSE,
                                      labels      = paste0("T", seq_along(partitions)),
                                      save_freq_tables = FALSE,
                                      n_cores     = 1L) {
  
  # --- 1. Input validation ---------------------------------------------------
  n_partitions = length(partitions)
  if (n_partitions < 2L) stop("At least 2 partitions are required.")
  if (length(labels) != n_partitions) {
    stop("`labels` must have the same length as `partitions`.")
  }
  
  n_sites     = length(sites)
  norm_factor = if (normalized) (1 / sqrt(2)) else 1
  
  # --- 2. Per-site worker ----------------------------------------------------
  compute_site = function(site_idx) {
    counts = get_site_counts(partitions, site_idx, aa_levels)  # aa_levels x n_partitions
    
    # Column-wise proportions; guard against all-zero columns (empty partitions)
    cs           = colSums(counts)
    cs[cs == 0L] = 1L
    props        = sweep(counts, 2L, cs, "/")
    
    p0     = props[, 1L]
    others = props[, -1L, drop = FALSE]
    
    dists = norm_factor * sqrt(colSums((sqrt(others) - sqrt(p0))^2))
    
    list(dists  = dists,
         counts = counts,
         props  = props)
  }
  
  # --- 3. Execute: parallel or sequential ------------------------------------
  use_parallel = n_cores > 1L && n_sites > 1L
  
  if (use_parallel) {
    future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(future::sequential), add = TRUE)
    
    site_results = future.apply::future_lapply(
      sites,
      compute_site,
      future.seed = TRUE
    )
  } else {
    site_results = lapply(sites, compute_site)
  }
  
  # --- 4. Assemble output matrix ---------------------------------------------
  dist_matrix           <- do.call(rbind, lapply(site_results, `[[`, "dists"))
  rownames(dist_matrix) <- as.character(sites)
  colnames(dist_matrix) <- labels[-1]          # T2, T3, ...
  storage.mode(dist_matrix) <- "double"        # ensure pure numeric for ecp/ks.cp3o
  
  # --- 5. Return -------------------------------------------------------------
  # Default: return the matrix directly so callers can do t(hell_mat) without
  # unpacking a list. When freq tables are requested, return the full named
  # list for backward compatibility with build_distance_matrix / legacy code.
  if (save_freq_tables) {
    list(
      Sites               = sites,
      Hellinger_Distances = dist_matrix,
      Frequency_Tables    = lapply(site_results, `[[`, "counts"),
      Proportions_Tables  = lapply(site_results, `[[`, "props")
    )
  } else {
    dist_matrix
  }
}