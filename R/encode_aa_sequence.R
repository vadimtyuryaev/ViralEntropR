#' @title Encode Amino Acid Sequences
#' @description Converts character matrices of amino acids to numeric representations.
#' @export
#'
#' @param matrix_input A character matrix of amino acid codes.
#' @return A numeric matrix of the same dimensions as input, where mapped values
#'   are 1-25 and any character not in the standard alphabet maps to 0.
#'
#' @examples
#' # 1. Encode a simple matrix of sequences
#' seq_mat = matrix(c("A", "R", "N", "D"), nrow = 2, byrow = TRUE)
#' encoded = encode_aa_sequence(seq_mat)
#'
#' # 2. Handle gaps and unknown characters
#' # '-' maps to 25, '?' (not in alphabet) maps to 0
#' gapped_mat = matrix(c("A", "-", "G", "?"), nrow = 1)
#' encode_aa_sequence(gapped_mat)
encode_aa_sequence = function(matrix_input) {
  # Internal mapping - 25 standard/ambiguous codes
  aa_map = c(
    "A" = 1,   "R" = 2,   "N" = 3,   "D" = 4,   "C" = 5,
    "Q" = 6,   "E" = 7,   "G" = 8,   "H" = 9,   "I" = 10,
    "L" = 11,  "K" = 12,  "M" = 13,  "F" = 14,  "P" = 15,
    "S" = 16,  "T" = 17,  "W" = 18,  "Y" = 19,  "V" = 20,
    "B" = 21,  "Z" = 22,  "X" = 23,  "*" = 24,  "-" = 25
  )
  
  if (!is.matrix(matrix_input)) matrix_input = as.matrix(matrix_input)
  
  # Normalize input to uppercase
  input_vec = toupper(as.vector(matrix_input))
  
  # Map characters to numbers
  # Values not in the map (including your specified unknown_char) become NA, then 0
  encoded_vec = aa_map[input_vec]
  encoded_vec[is.na(encoded_vec)] = 0
  
  matrix(encoded_vec, 
         nrow = nrow(matrix_input), 
         ncol = ncol(matrix_input),
         dimnames = dimnames(matrix_input))
}