#' @title Decode Amino Acid Sequences
#' @description Converts numeric matrices back to amino acid characters.
#' @param matrix_input Numeric matrix.
#' @return Character matrix. Values outside 1-25 (including 0) are returned as
#'   \code{"0"}, matching the behaviour of the original \code{numbers_to_letters_vectorized}.
#' @export
#'
#' @examples
#' # 1. Decode a numeric matrix
#' num_mat = matrix(c(1, 2, 25, 10), nrow = 2, byrow = TRUE)
#' decoded = decode_aa_sequence(num_mat)
#' print(decoded)
#'
#' # 2. Round-trip consistency check (excluding unknowns)
#' orig = matrix(c("A", "C", "W", "G"), nrow = 2)
#' enc = encode_aa_sequence(orig)
#' dec = decode_aa_sequence(enc)
#' all.equal(orig, dec)
decode_aa_sequence = function(matrix_input) {
  aa_chars = c("A","R","N","D","C","Q","E","G","H","I",
                "L","K","M","F","P","S","T","W","Y","V",
                "B","Z","X","*","-")
  
  if (!is.matrix(matrix_input)) matrix_input = as.matrix(matrix_input)
  input_vec = as.vector(matrix_input)
  decoded_vec = aa_chars[input_vec]
  decoded_vec[is.na(decoded_vec) | input_vec < 1 | input_vec > 25] = "0"
  
  matrix(decoded_vec, nrow = nrow(matrix_input), ncol = ncol(matrix_input),
         dimnames = dimnames(matrix_input))
}