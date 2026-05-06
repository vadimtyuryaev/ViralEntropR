#' @title Convert FASTA Object to Character Matrix
#' @description Converts an \code{AAStringSet} object loaded via
#'   \code{\link[Biostrings]{readAAStringSet}} into a character matrix where
#'   rows are sequences and columns are residue positions (sites). Inverse
#'   structural transformation of \code{\link{encode_aa_sequence}}'s expected
#'   input shape.
#'
#' @details
#' \strong{Alignment.} The function expects an aligned \code{AAStringSet}
#' — all sequences of equal width. Unaligned input is accepted and shorter
#' sequences are right-padded with the gap character \code{"-"} to match
#' the longest sequence, but downstream entropy-based analysis assumes
#' positional homology across rows; if sequences in the input are not
#' biologically aligned, results from per-site computations will not be
#' meaningful. For unaligned input, run a multiple-sequence alignment
#' (e.g. \code{msa::msa()} or \code{DECIPHER::AlignSeqs()}) before calling
#' this function.
#'
#' \strong{Performance.} Conversion is fully vectorised: all sequences are
#' coerced to a single character string vector, split simultaneously, and
#' reshaped into a matrix in one operation. No per-row loop, no
#' intermediate list of split sequences kept alive — substantially faster
#' than per-row \code{strsplit} on large inputs (100k+ sequences).
#'
#' @param fsta An \code{AAStringSet} object, typically the output of
#'   \code{\link[Biostrings]{readAAStringSet}}. May be aligned or
#'   unaligned (see Details).
#'
#' @return A character matrix with \code{length(fsta)} rows and
#'   \code{max(nchar(as.character(fsta)))} columns. Each cell contains a
#'   single-character amino acid code from the input sequences (or the
#'   gap character \code{"-"} for padded positions in unaligned input).
#'   The matrix has no row or column names; sequence names from the
#'   \code{AAStringSet} are not carried over. An empty input
#'   (\code{length(fsta) == 0}) returns an empty 0-by-0 character matrix.
#'
#' @seealso \code{\link{encode_aa_sequence}} for converting the resulting
#'   character matrix to an integer-encoded matrix;
#'   \code{\link{filter_ambiguous_sequences}} for removing rows containing
#'   ambiguous residues; \code{\link[Biostrings]{readAAStringSet}} for
#'   loading FASTA files into the input format.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Convert the bundled sample to a character matrix.
#' path  <- system.file("extdata", "sarscov2_sample.fasta.gz",
#'                       package = "ViralEntropR")
#' fasta <- Biostrings::readAAStringSet(path)
#' mat   <- fasta_to_char_matrix(fasta)
#' dim(mat)
#' mat[1:3, 1:10]
#' }
fasta_to_char_matrix <- function(fsta) {
  
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Package 'Biostrings' is required. Install it with:\n",
         "  install.packages('BiocManager')\n",
         "  BiocManager::install('Biostrings')",
         call. = FALSE)
  }
  
  n <- length(fsta)
  
  # Empty input: return empty matrix, schema-stable for downstream consumers.
  if (n == 0L) {
    return(matrix(character(0), nrow = 0L, ncol = 0L))
  }
  
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
