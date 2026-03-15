#' @title Download the Full SARS-CoV-2 Surface Glycoprotein Dataset
#' @description Downloads the complete SARS-CoV-2 surface glycoprotein amino
#'   acid sequence dataset (n = 137,132 sequences) from its Zenodo archive.
#'   The file is cached locally after the first download so subsequent calls
#'   return immediately.
#'
#' @details
#' The dataset was originally downloaded from the NCBI SARS-CoV-2 Data Hub on
#' October 12, 2021 (taxid: 2697049; complete sequences; surface glycoprotein
#' only), yielding 137,132 sequences. It is archived on Zenodo under
#' DOI \doi{10.5281/zenodo.XXXXXXX} and released under CC0 1.0 Universal.
#'
#' The file is cached in a persistent, user-level directory that survives
#' across R sessions (see \code{\link[tools]{R_user_dir}}). To find the cache
#' location:
#' \preformatted{
#' tools::R_user_dir("ViralEntropR", which = "cache")
#' }
#'
#' @param destdir Character or \code{NULL}. Directory to save the file.
#'   Defaults to a persistent package-specific cache via
#'   \code{\link[tools]{R_user_dir}}. Supply a path to override.
#' @param force Logical. If \code{TRUE}, re-downloads even if a cached copy
#'   already exists. Default \code{FALSE}.
#'
#' @return Invisibly returns the full file path to the downloaded
#'   \code{.fasta.gz} file. The path can be passed directly to
#'   \code{\link[Biostrings]{readAAStringSet}}.
#'
#' @references
#' Tyuryaev V, et al. (2025). SARS-CoV-2 surface glycoprotein sequences for
#' ViralEntropR. \emph{Zenodo}. \doi{10.5281/zenodo.XXXXXXX}
#'
#' Hatcher EL, Zhdanov SA, Bao Y, et al. (2017).
#' Virus Variation Resource — improved response to emergent viral outbreaks.
#' \emph{Nucleic Acids Research}, 45(D1), D482–D490.
#' \doi{10.1093/nar/gkw1065}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{fasta_to_char_matrix}} — convert to character matrix
#'   \item \code{\link{extract_fasta_dates}} — extract collection dates
#'   \item \code{\link{extract_fasta_countries}} — extract country metadata
#'   \item \code{\link{filter_ambiguous_sequences}} — remove ambiguous residues
#'   \item \code{\link{encode_aa_sequence}} — encode for pipeline entry
#' }
#'
#' @importFrom tools R_user_dir
#' @export
#'
#' @examples
#' \dontrun{
#' # First call: downloads ~35 MB and caches locally
#' path <- download_sarscov2_data()
#'
#' # Subsequent calls: returns cached path instantly
#' path <- download_sarscov2_data()
#'
#' # Load into R
#' fasta <- Biostrings::readAAStringSet(path)
#' length(fasta)   # 137132
#'
#' # Force a fresh download (e.g. after a Zenodo version update)
#' path <- download_sarscov2_data(force = TRUE)
#'
#' # Save to a custom directory instead of the default cache
#' path <- download_sarscov2_data(destdir = "~/my_data")
#' }
download_sarscov2_data <- function(destdir = NULL, force = FALSE) {
  
  # --- 1. Resolve cache directory --------------------------------------------
  if (is.null(destdir)) {
    destdir <- tools::R_user_dir("ViralEntropR", which = "cache")
  }
  if (!dir.exists(destdir)) {
    dir.create(destdir, recursive = TRUE)
  }
  
  dest <- file.path(destdir, "sarscov2_spike_ncbi_20211012.fasta")
  
  # Zenodo direct download URL — update XXXXXXX after deposit
  url <- paste0(
    "https://zenodo.org/record/XXXXXXX/files/",
    "sarscov2_spike_ncbi_20211012.fasta?download=1"
  )
  
  # --- 2. Return cached file if available ------------------------------------
  if (file.exists(dest) && !isTRUE(force)) {
    message("Using cached file: ", dest)
    message("Use force = TRUE to re-download.")
    return(invisible(dest))
  }
  
  # --- 3. Download -----------------------------------------------------------
  message(
    "Downloading SARS-CoV-2 surface glycoprotein dataset from Zenodo.\n",
    "  n = 137,132 sequences | ~173 MB uncompressed\n",
    "  DOI: 10.5281/zenodo.XXXXXXX\n",
    "  This is a one-time download. File will be cached at:\n  ", dest
  )
  
  tryCatch(
    utils::download.file(url, dest, mode = "wb", quiet = FALSE),
    error = function(e) {
      # Clean up partial file if download failed mid-way
      if (file.exists(dest)) unlink(dest)
      stop(
        "Download failed. Please check your internet connection.\n",
        "You can also download the file manually from:\n",
        "  https://zenodo.org/record/XXXXXXX\n",
        "and save it to: ", dest, "\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  # --- 4. Basic integrity check ----------------------------------------------
  fsize <- file.size(dest)
  if (fsize < 100e6) {   # less than 100 MB almost certainly means a failed download
    unlink(dest)
    stop(
      "Downloaded file appears incomplete (", round(fsize / 1e6, 1), " MB).\n",
      "Please try again or download manually from:\n",
      "  https://zenodo.org/record/XXXXXXX",
      call. = FALSE
    )
  }
  
  message(sprintf("Download complete (%.1f MB).", fsize / 1e6))
  invisible(dest)
}
