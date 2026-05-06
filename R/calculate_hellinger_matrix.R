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
#' \code{\link{get_site_counts}} (built on \code{\link[base]{tabulate}}).
#' Per-partition proportions and Hellinger distances are then computed by
#' fully vectorised matrix operations — no inner loop over partitions.
#'
#' @param partitions A list of data frames, one per time window. Each data frame
#'   must have numeric-encoded amino acid sequences as columns (integers 1 to
#'   \code{aa_levels}). This is typically the \code{$Partitions} element 
#'   returned by \code{\link{partition_time_windows}}.
#' @param sites Integer vector. Indices of the sites to analyse. Defaults to all
#'   sites (\code{seq_len(ncol(partitions[[1]]))}).
#' @param aa_levels Integer. Alphabet size. Must match the encoding used
#'   when the partitions were created: \code{\link{encode_aa_sequence}}
#'   produces values 1–25 by default (20 standard residues, three
#'   ambiguous codes, \code{*}, \code{-}). Default is \code{25L}.
#' @param normalized Logical. If \code{TRUE}, scales distances by
#'   \eqn{1/\sqrt{2}} to bound the result in \eqn{[0, 1]}. Default is
#'   \code{FALSE}.
#' @param labels Character vector of partition labels. Length must equal
#'   \code{length(partitions)}. The first label names the reference
#'   partition; subsequent labels become column names of the returned
#'   matrix. Defaults to \code{"T1"}, \code{"T2"}, …
#' @param include_freq_tables Logical. If \code{TRUE}, the return value also
#'   includes raw count tables and proportion tables for each site. Default is
#'   \code{FALSE}.
#'
#' @return When \code{include_freq_tables = FALSE} (default): a numeric matrix
#'   with rows corresponding to \code{sites} and columns corresponding to
#'   \code{labels[-1]}, each entry being the Hellinger distance from the
#'   reference partition (\code{labels[1]}) at that site. Row names are
#'   the site indices; column names are taken from \code{labels[-1]}.
#'   With default labels, the reference partition is \code{"T1"} and the
#'   matrix has columns \code{"T2"}, \code{"T3"}, ….
#'
#'   When \code{include_freq_tables = TRUE}: a named list with elements
#'   \code{Sites}, \code{Hellinger_Distances}, \code{Frequency_Tables}, and
#'   \code{Proportions_Tables}.
#'   
#' @seealso
#' \code{\link{partition_time_windows}} for producing temporally
#'   partitioned data, \code{\link{encode_aa_sequence}} for integer
#'   encoding consistent with \code{aa_levels}, and
#'   \code{\link{detect_changepoints_ecp}} or
#'   \code{\link{detect_changepoints_hdcp}} for downstream change-point
#'   detection on the returned matrix.
#'
#' @examples
#' # Toy 3-partition data: site 1 starts homogeneous (all Alanine, 1),
#' # acquires Valine (20) across two later partitions; site 2 stays
#' # constant (all Valine throughout).
#' p1 = data.frame(s1 = c(1, 1, 1, 1, 1),  s2 = c(20, 20, 20, 20, 20))
#' p2 = data.frame(s1 = c(20, 20, 20, 20, 20), s2 = c(20, 20, 20, 20, 20))
#' p3 = data.frame(s1 = c(1, 1, 20, 20, 20), s2 = c(20, 20, 20, 20, 20))
#' parts = list(T1 = p1, T2 = p2, T3 = p3)
#'
#' result = calculate_hellinger_matrix(parts, sites = 1:2)
#' print(result)
#' # Site 1 has nonzero distance in T2 and T3 (composition shift).
#' # Site 2 has zero distance throughout (constant).
#'
#' # With raw frequency tables for inspection.
#' result2 = calculate_hellinger_matrix(parts, sites = 1:2,
#'                                       include_freq_tables = TRUE)
#' print(result2$Frequency_Tables[[1]])
#' @export
calculate_hellinger_matrix = function(partitions,
                                      sites            = seq_len(ncol(partitions[[1]])),
                                      aa_levels        = 25L,
                                      normalized       = FALSE,
                                      labels           = paste0("T", seq_along(partitions)),
                                      include_freq_tables = FALSE) {
  
  # --- 1. Input validation ---------------------------------------------------
  n_partitions = length(partitions)
  if (n_partitions < 2L) stop("At least 2 partitions are required.")
  if (length(labels) != n_partitions)
    stop("`labels` must have the same length as `partitions`.")
  
  norm_factor = if (normalized) (1 / sqrt(2)) else 1
  
  # --- 2. Per-site computation -----------------------------------------------
  # compute_site is fully vectorized over partitions: sweep/colSums operate on
  # the full aa_levels x n_partitions matrix in one call per site.
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
  
  # --- 3. Execute over all sites ---------------------------------------------
  site_results = lapply(sites, compute_site)
  
  # --- 4. Assemble output matrix ---------------------------------------------
  dist_matrix              <- do.call(rbind, lapply(site_results, `[[`, "dists"))
  rownames(dist_matrix)    <- as.character(sites)
  colnames(dist_matrix)    <- labels[-1]          # T2, T3, ...
  storage.mode(dist_matrix) <- "double"           # ensure pure numeric for ecp/ks.cp3o
  
  # --- 5. Return -------------------------------------------------------------
  if (include_freq_tables) {
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