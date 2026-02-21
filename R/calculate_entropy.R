#' @title Calculate Shannon Entropy
#' @description Computes the Shannon entropy of a categorical vector.
#'
#' @details
#' Entropy is calculated as \eqn{H(X) = -\sum p(x) \log_b p(x)}, where \eqn{p(x)}
#' is the proportion of observations belonging to category \eqn{x}.
#'
#' @param vctr A vector (character, factor, or integer) representing categorical data.
#' @param base A numeric scalar. The base of the logarithm. Default is 2.
#' @param precision Integer. The number of decimal places to round the result to. Default is 6.
#'
#' @return A numeric scalar representing the entropy. Returns 0 if the vector contains
#' only one unique value or has length 0.
#'
#' @export
#' 
#' @examples
#' seq_vec = c("A", "A", "T", "G", "C", "A")
#' calculate_entropy(seq_vec)
#'
#' # Pure homogeneity
#' calculate_entropy(rep("A", 10))
#' 
calculate_entropy = function(vctr, base = 2, precision = 6) {
  if (length(vctr) < 2) return(0)
  
  # Calculate proportions
  counts = table(vctr)
  probs  = as.vector(counts) / sum(counts)
  
  # Calculate Entropy (handle 0 probabilities implicitly)
  e = -sum(probs * log(probs, base = base))
  
  round(e, precision)
}