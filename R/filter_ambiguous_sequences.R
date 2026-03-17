#' @title Remove Sequences Containing Ambiguous Amino Acids
#' @description Removes rows (sequences) that contain at least one ambiguous
#'   amino acid residue (B, J, X, or Z) from a sequence matrix. Accepts both
#'   numeric-encoded and character matrices.
#'
#' @details
#' Ambiguous residues map to the following numeric codes under the package
#' encoding: B = 21, Z = 22, X = 23, J = 0 (unrecognised/gap).
#'
#' The detection step is fully vectorized: a single logical matrix comparison
#' followed by \code{\link[base]{rowSums}} replaces the original row-by-row
#' loop, giving a substantial speed improvement on large matrices (100k+ rows).
#'
#' @param NumMatrix A matrix. Rows are sequences, columns are sites. Either
#'   integer-encoded (\code{option = 1}) or character (\code{option = 2}).
#' @param option Integer. \code{1} (default) for numeric-encoded matrices;
#'   \code{2} for character matrices.
#'
#' @return A named list:
#' \item{OriginalDim}{Character string reporting the number of input sequences.}
#' \item{NewDim}{Character string reporting the number of sequences after
#'   filtering.}
#' \item{NumberAmbiguous}{Character string reporting the number of sequences
#'   that contained at least one ambiguous residue.}
#' \item{RangeAmbiguous}{Character string reporting the min and max count of
#'   ambiguous residues per removed sequence, or \code{"No ambiguous sequences
#'   found"} when none were removed.}
#' \item{DeletedSeqId}{Integer vector of row indices that were removed.}
#' \item{FilteredMatrix}{The filtered matrix with ambiguous rows removed.}
#'
#' @examples
#' set.seed(1)
#' m <- matrix(sample(1:25, 500, replace = TRUE), nrow = 50, ncol = 10)
#' # Inject some ambiguous codes
#' m[c(3, 17, 42), sample(1:10, 3)] <- c(21, 23, 0)
#' result <- filter_ambiguous_sequences(m, option = 1)
#' cat(result$NumberAmbiguous, "\n")
#' cat(result$RangeAmbiguous, "\n")
#' dim(result$FilteredMatrix)
#'
#' @export
filter_ambiguous_sequences <- function(NumMatrix, option = 1) {

  if (!option %in% 1:2) stop("`option` must be 1 (numeric) or 2 (character).")
  if (!is.matrix(NumMatrix)) stop("`NumMatrix` must be a matrix.")

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
