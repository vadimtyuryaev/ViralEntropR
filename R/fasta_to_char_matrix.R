#' @title Convert FASTA Object to Character Matrix
#' @description Creates a character matrix from an \code{AAStringSet} object
#'   loaded via \code{\link[Biostrings]{readAAStringSet}}. Rows represent
#'   sequences; columns represent residue positions (sites).
#'
#' @details
#' The function expects an aligned \code{AAStringSet} (all sequences of equal
#' width). It can also handle unaligned objects: sequences shorter than the
#' maximum width are padded with \code{"-"} on the right.
#'
#' The conversion is fully vectorized: all sequences are coerced to a single
#' character string, split simultaneously, and reshaped into a matrix in one
#' operation — replacing the original row-by-row \code{strsplit} loop for a
#' substantial speed improvement on large inputs.
#'
#' @param fsta An \code{AAStringSet} object.
#'
#' @return A character matrix with \code{length(fsta)} rows and
#'   \code{max(nchar(as.character(fsta)))} columns.
#'
#' @examples
#' \dontrun{
#' fasta  <- Biostrings::readAAStringSet("sequences.fasta")
#' mat    <- fasta_to_char_matrix(fasta)
#' dim(mat)
#' mat[1:3, 1:10]
#' }
#'
#' @export
fasta_to_char_matrix <- function(fsta) {
  
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required. Install it with:\n",
         "  install.packages('BiocManager')\n",
         "  BiocManager::install('Biostrings')",
         call. = FALSE)
  }

  n        <- length(fsta)
  seqs     <- as.character(fsta)          # named character vector, one string per seq
  widths   <- nchar(seqs)
  m        <- max(widths)

  # For aligned inputs all widths are equal and the reshape is exact.
  # For unaligned inputs, pad shorter sequences with "-" so every string
  # reaches width m before splitting — this preserves the original behaviour
  # of filling the pre-allocated matrix (which used 0; here "-" is more
  # meaningful for amino acid data).
  if (any(widths < m)) {
    seqs <- formatC(seqs, width = -m, flag = "-")   # left-align, pad right with spaces
    # replace padding spaces with gap character
    seqs <- gsub(" ", "-", seqs, fixed = TRUE)
  }

  # Split all sequences at once: strsplit returns a list of length n, each
  # element a character vector of length m. unlist + matrix reshapes in one
  # step — no per-row loop, no intermediate list kept alive.
  seq_matrix <- matrix(
    unlist(strsplit(seqs, "", fixed = TRUE), use.names = FALSE),
    nrow  = n,
    ncol  = m,
    byrow = TRUE
  )

  seq_matrix
}
