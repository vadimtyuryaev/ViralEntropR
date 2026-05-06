#' @title Remove Sequences Containing Ambiguous Residues
#' @description Removes rows (sequences) that contain at least one ambiguous
#'   amino acid residue (B, J, X, or Z) — and, under integer-encoded input,
#'   any unrecognised character — from a sequence matrix. Accepts both
#'   integer-encoded matrices and character matrices.
#'
#' @details
#' \strong{What is removed.} Sequences are flagged for removal if any of
#' their residue positions contain one of the four IUPAC ambiguous codes:
#' \itemize{
#'   \item \code{B} — Aspartate / Asparagine.
#'   \item \code{J} — Leucine / Isoleucine.
#'   \item \code{X} — any residue.
#'   \item \code{Z} — Glutamate / Glutamine.
#' }
#'
#' \strong{What is NOT removed.} Standard alignment gaps (\code{-},
#' integer code \code{25}) are retained — gaps represent known absences
#' rather than uncertain identities and are typically positionally
#' meaningful in aligned data. Sequences containing only canonical 20
#' amino acids and gaps are kept.
#'
#' \strong{How input mode is handled.}
#' \itemize{
#'   \item \code{option = 1} (integer-encoded input from
#'     \code{\link{encode_aa_sequence}}): rows are removed if any cell
#'     equals \code{0} (unrecognised — including J, NA, empty, lowercase
#'     mismatches, byte-order marks, and other characters that fell
#'     outside the encoding alphabet), \code{21} (B), \code{22} (Z), or
#'     \code{23} (X). The \code{0} sentinel acts as a catch-all for
#'     anything not in the 25-symbol alphabet.
#'   \item \code{option = 2} (character input from
#'     \code{\link{fasta_to_char_matrix}}): rows are removed if any cell
#'     is exactly \code{"B"}, \code{"J"}, \code{"X"}, or \code{"Z"}.
#'     Unrecognised characters in character input (e.g. lowercase letters,
#'     \code{NA}, empty strings) are NOT caught at this stage; encode
#'     first if you need that catch-all behaviour.
#' }
#'
#' \strong{Performance.} Detection is fully vectorised: a single logical
#' matrix comparison followed by \code{\link[base]{rowSums}} counts
#' ambiguous residues per sequence in one C-level call, replacing the
#' original row-by-row loop for a substantial speed improvement on
#' large matrices (100k+ rows).
#'
#' @param NumMatrix A matrix. Rows are sequences, columns are sites. Either
#'   integer-encoded (\code{option = 1}) or character (\code{option = 2}).
#'   Despite the name, character matrices are also accepted under
#'   \code{option = 2}.
#' @param option Integer. \code{1} (default) for integer-encoded matrices
#'   produced by \code{\link{encode_aa_sequence}}; \code{2} for character
#'   matrices produced by \code{\link{fasta_to_char_matrix}}.
#'
#' @return A named list:
#' \item{OriginalDim}{Character string reporting the number of input
#'   sequences.}
#' \item{NewDim}{Character string reporting the number of sequences
#'   remaining after filtering.}
#' \item{NumberAmbiguous}{Character string reporting the number of
#'   sequences that contained at least one ambiguous residue.}
#' \item{RangeAmbiguous}{Character string reporting the min and max count
#'   of ambiguous residues per removed sequence, or
#'   \code{"No ambiguous sequences found"} when none were removed.}
#' \item{DeletedSeqId}{Integer vector of row indices that were removed.
#'   Empty integer vector if nothing was removed.}
#' \item{FilteredMatrix}{The filtered matrix with ambiguous rows removed,
#'   preserving the original column structure and storage mode.}
#'
#' @seealso \code{\link{encode_aa_sequence}} and
#'   \code{\link{fasta_to_char_matrix}} for producing the typical input;
#'   \code{\link{decode_aa_sequence}} for inspecting the surviving
#'   sequences in character form.
#'
#' @export
#'
#' @examples
#' # Synthetic example: 50 sequences, 10 sites, drawn from canonical residues.
#' set.seed(1)
#' m <- matrix(sample(1:20, 500, replace = TRUE), nrow = 50, ncol = 10)
#' # Inject ambiguous codes into 3 specific rows: 21 (B), 23 (X), 0 (unrecognised).
#' m[c(3, 17, 42), sample(1:10, 3)] <- c(21, 23, 0)
#' result <- filter_ambiguous_sequences(m, option = 1)
#' cat(result$NumberAmbiguous, "\n")
#' cat(result$RangeAmbiguous, "\n")
#' dim(result$FilteredMatrix)
#'
#' # Character-mode example.
#' chr <- matrix(c("M", "K", "T", "I", "I", "X", "K", "T", "I", "I"),
#'                nrow = 2, byrow = TRUE)
#' filter_ambiguous_sequences(chr, option = 2)$DeletedSeqId
filter_ambiguous_sequences <- function(NumMatrix, option = 1) {

  if (!option %in% 1:2) stop("`option` must be 1 (numeric) or 2 (character).", 
                             call. = FALSE)
  
  if (!is.matrix(NumMatrix)) stop("`NumMatrix` must be a matrix.", 
                                  call. = FALSE)

  # Vectorized detection: build a logical matrix of ambiguous positions,
  # then rowSums counts ambiguous residues per sequence in one C-level call.
  # This replaces the original row-by-row for loop — O(n*m) but fully in C.
  if (option == 1) {
    ambig_mat <- NumMatrix == 21L | NumMatrix == 22L |
                 NumMatrix == 23L | NumMatrix == 0L
  } else {
    ambig_mat <- NumMatrix == "B" | NumMatrix == "J" |
                 NumMatrix == "X" | NumMatrix == "Z"
  }

  ambig_per_row <- rowSums(ambig_mat)   # integer vector, length = nrow
  id_delete     <- which(ambig_per_row > 0L)
  l_w           <- ambig_per_row[id_delete]

  filtered <- if (length(id_delete) > 0L) NumMatrix[-id_delete, , drop = FALSE]
              else NumMatrix

  range_str <- if (length(l_w) > 0L)
    sprintf(
      "Number of ambiguous protein characters per sequence varies between %d and %d",
      min(l_w), max(l_w)
    )
  else
    "No ambiguous sequences found"

  list(
    OriginalDim     = sprintf("Number of original sequences is %d",  nrow(NumMatrix)),
    NewDim          = sprintf("Number of filtered sequences is %d",  nrow(filtered)),
    NumberAmbiguous = sprintf(
      "Number of sequences containing at least one of B, X, Z or J characters is %d",
      length(id_delete)
    ),
    RangeAmbiguous  = range_str,
    DeletedSeqId    = id_delete,
    FilteredMatrix  = filtered
  )
}
