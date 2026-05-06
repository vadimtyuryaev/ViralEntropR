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
#'
#' @examples
#' # 1. Create dummy partitions.
#' # Partition 1: site 1 has three 'A's (code 1).
#' p1 = data.frame(s1 = c(1, 1, 1))
#' # Partition 2: site 1 has one 'A' (1) and two 'R's (2).
#' p2 = data.frame(s1 = c(1, 2, 2))
#' parts = list(p1, p2)
#'
#' # 2. Get counts for site 1.
#' # Internal function — accessed via the triple-colon operator.
#' counts = ViralEntropR:::get_site_counts(parts, site_index = 1,
#'                                          alphabet_size = 25)
#'
#' # Row 1 (A) is [3, 1]; row 2 (R) is [0, 2].
#' print(counts[1:5, ])
get_site_counts = function(partitions, site_index, alphabet_size) {
  site_vecs = lapply(partitions, `[[`, site_index)
  vapply(site_vecs, tabulate, integer(alphabet_size), nbins = alphabet_size)
}