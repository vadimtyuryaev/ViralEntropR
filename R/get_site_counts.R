#' @title Get Site Counts Matrix (Internal)
#' @description Extract and tabulate amino acid counts for a specific site
#'   across all partitions.
#'
#' @details Built on \code{\link[base]{tabulate}} for fast integer counting.
#'   Consumed internally by \code{\link{calculate_hellinger_matrix}}.
#'
#' @param partitions List of data frames, one per time window. Each must
#'   have integer-encoded amino acid columns.
#' @param site_index Integer. Index of the site (column) to count.
#' @param alphabet_size Integer. Size of the amino acid alphabet
#'   (typically \code{25} for the package's default encoding).
#' @return Integer matrix with \code{alphabet_size} rows (one per amino
#'   acid code, 1 through \code{alphabet_size}) and one column per
#'   partition. Each cell is the count of that amino acid at
#'   \code{site_index} in that partition.
#' @keywords internal
get_site_counts = function(partitions, site_index, alphabet_size) {
  site_vecs = lapply(partitions, `[[`, site_index)
  vapply(site_vecs, tabulate, integer(alphabet_size), nbins = alphabet_size)
}