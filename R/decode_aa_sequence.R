#' @title Decode Amino Acid Sequences
#' @description Converts an integer-encoded matrix of amino acids back to its
#'   character representation under the package's 25-symbol alphabet. Inverse
#'   of \code{\link{encode_aa_sequence}}.
#'
#' @details
#' Decoding is a single vectorised lookup against the fixed alphabet
#' (\code{A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V} for
#' the twenty standard residues, then \code{B, Z, X, *, -} for ambiguous
#' codes, stop codons, and gaps). The input matrix is flattened, indexed
#' against the alphabet vector in one operation, and reshaped — there is
#' no per-row or per-element loop.
#'
#' Values outside the valid range \code{1:25} (including \code{0}, which
#' \code{\link{encode_aa_sequence}} produces for unrecognised characters)
#' are returned as the sentinel string \code{"0"}. This preserves
#' round-trip consistency with the encoder: encoding then decoding any
#' character originally outside the alphabet yields \code{"0"} rather
#' than throwing an error. \code{NA} values are also mapped to
#' \code{"0"}.
#'
#' Row and column names of the input matrix are preserved on the output.
#'
#' @param matrix_input Numeric matrix of integer-encoded amino acids,
#'   typically the output of \code{\link{encode_aa_sequence}}. Values
#'   outside \code{1:25} (including \code{0} and \code{NA}) decode to
#'   the sentinel string \code{"0"}. A non-matrix input is coerced via
#'   \code{\link[base]{as.matrix}}.
#'
#' @return A character matrix of the same dimensions as \code{matrix_input},
#'   with the same \code{dimnames}. Each cell is either a one-character
#'   amino acid code from the 25-symbol alphabet, or the sentinel
#'   \code{"0"} for out-of-range and missing values.
#'
#' @seealso \code{\link{encode_aa_sequence}} for the inverse operation;
#'   \code{\link{fasta_to_char_matrix}} for the FASTA-to-character-matrix
#'   step that typically precedes encoding.
#'
#' @export
#'
#' @examples
#' # 1. Decode a numeric matrix.
#' num_mat = matrix(c(1, 2, 25, 10), nrow = 2, byrow = TRUE)
#' decoded = decode_aa_sequence(num_mat)
#' print(decoded)
#' # 1 -> "A", 2 -> "R", 25 -> "-", 10 -> "I"
#'
#' # 2. Round-trip consistency check (excluding unknowns).
#' orig = matrix(c("A", "C", "W", "G"), nrow = 2)
#' enc  = encode_aa_sequence(orig)
#' dec  = decode_aa_sequence(enc)
#' all.equal(orig, dec)
#'
#' # 3. Out-of-range and NA values both decode to the sentinel "0".
#' decode_aa_sequence(matrix(c(0, NA, 30, 5), nrow = 2))
decode_aa_sequence = function(matrix_input) {
  aa_chars = c("A","R","N","D","C","Q","E","G","H","I",
               "L","K","M","F","P","S","T","W","Y","V",
               "B","Z","X","*","-")
  
  if (!is.matrix(matrix_input)) matrix_input = as.matrix(matrix_input)
  
  # Mask out-of-range and NA values so indexing preserves length:
  # 0 would collapse a position; NA preserves it (and we map NA -> "0" below).
  input_vec  = as.vector(matrix_input)
  safe_idx   = input_vec
  safe_idx[is.na(safe_idx) | safe_idx < 1L | safe_idx > 25L] = NA_integer_
  
  decoded_vec = aa_chars[safe_idx]
  decoded_vec[is.na(decoded_vec)] = "0"
  
  matrix(decoded_vec, nrow = nrow(matrix_input), ncol = ncol(matrix_input),
         dimnames = dimnames(matrix_input))
}