#' @title Get Site Counts Matrix (Internal)
#' @description Extract and tabulate amino acid counts for a specific site across all partitions.
#' @param partitions List of data frames.
#' @param site_index Integer.
#' @param alphabet_size Integer.
#' @return Numeric matrix (rows=AA, cols=Time).
#' @keywords internal
#' 
#' @examples
#' # 1. Create dummy partitions
#' # Partition 1: Site 1 has three 'A's (code 1)
#' p1 = data.frame(s1 = c(1, 1, 1))
#' 
#' # Partition 2: Site 1 has one 'A' (1) and two 'R's (2)
#' p2 = data.frame(s1 = c(1, 2, 2))
#' 
#' parts = list(p1, p2)
#' 
#' # 2. Get counts for Site 1
#' # Note: Since this is internal, we call it directly if loaded, 
#' # or use ViralEntroR:::get_site_counts if installed.
#' counts = ViralEntropR:::get_site_counts(parts, site_index = 1, alphabet_size = 25)
#' 
#' # Expect Row 1 (A) to be [3, 1]
#' # Expect Row 2 (R) to be [0, 2]
#' print(counts[1:5, ])
get_site_counts = function(partitions, site_index, alphabet_size) {
  site_vecs = lapply(partitions, `[[`, site_index)
  vapply(site_vecs, tabulate, integer(alphabet_size), nbins = alphabet_size)
}