#' @title Cluster Sequence Sites by Entropy
#' @description Identifies groups of sites with similar variability using
#'   Gaussian Mixture Models (GMM).
#'
#' @details
#' Groups sites based on their Shannon entropy values via
#' \code{\link[mclust]{Mclust}}. Preprocessing steps remove invariant sites
#' (entropy = 0) and/or singleton sites (entropy corresponding to exactly one
#' differing sequence across \code{nr} rows).
#'
#' Class assignment rules (applied in priority order):
#' \itemize{
#'   \item \strong{No rows remaining after filtering}: empty DataFrame returned
#'     with a zero-length \code{class} column (consistent schema).
#'   \item \strong{Single row remaining}: class \code{1} assigned directly;
#'     Mclust is not called (undefined on 1 observation).
#'   \item \strong{All entropies identical}: class \code{999} for all sites
#'     (sentinel — one undifferentiated group).
#'   \item \strong{Normal Mclust result}: raw class labels \code{1, 2, ..., G}.
#'     Note these are Mclust's own arbitrary integer labels — call
#'     \code{\link{relabel_entropy_classes}} on the returned DataFrame to
#'     standardize so that class 1 = highest entropy group.
#'   \item \strong{Genuine Mclust failure}: empty DataFrame returned (same
#'     schema as the no-rows case). Treats the partition as uninformative,
#'     which is scientifically correct — a failed clustering produces no
#'     classifiable information.
#' }
#'
#' \strong{Note on \code{mclustBIC}:} \code{mclust::Mclust} internally calls
#' \code{mclustBIC} via character-string dispatch in a context that searches
#' the calling package's namespace rather than mclust's own. The
#' \code{@importFrom mclust Mclust mclustBIC} directive pulls both symbols
#' into this package's namespace, making \code{mclustBIC} findable without
#' any \code{library(mclust)} call by the user.
#'
#' @param entropies Numeric vector of entropy values, one per site.
#' @param nr Integer. Total number of sequences used to compute the entropies.
#'   Required when \code{removesngl = TRUE}.
#' @param nsites Integer. Expected number of sites. If it mismatches
#'   \code{length(entropies)}, the actual length is used with a warning.
#'   Default is \code{length(entropies)}.
#' @param precision Integer. Decimal places for rounding during singleton
#'   threshold comparison and the all-identical uniqueness check. Default
#'   is \code{6}.
#' @param removez Logical. If \code{TRUE}, removes sites with entropy = 0
#'   (invariant sites), using a small tolerance \code{1e-9} to absorb
#'   floating-point near-zeros. Default is \code{TRUE}.
#' @param removesngl Logical. If \code{TRUE}, removes sites whose entropy
#'   equals the singleton value (one differing sequence out of \code{nr}).
#'   Uses tolerance-based comparison. Default is \code{TRUE}.
#' @param transfr A function, or an object of class \code{transform} with a
#'   \code{$transform()} method, applied to entropies before clustering.
#'   Default is \code{NULL} (no transformation).
#' @param verbose Logical. If \code{TRUE}, emits diagnostic warnings for
#'   non-fatal events (empty partitions, Mclust failures, etc.). Default is
#'   \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[mclust]{Mclust}}.
#'
#' @return A named list with two elements:
#' \item{FitObject}{The raw \code{Mclust} result, or a minimal
#'   \code{list(classification = integer(0L))} when clustering was bypassed
#'   or failed.}
#' \item{DataFrame}{A data frame with columns \code{sites} (original site
#'   indices), \code{entropies} (values after any transformation), and
#'   \code{class} (GMM cluster label). The \code{class} column is
#'   \strong{always} present in every return path, including zero-row
#'   DataFrames. Downstream consumers need only guard on \code{nrow(df) > 0}
#'   before accessing class values. Raw Mclust labels are returned as-is;
#'   call \code{\link{relabel_entropy_classes}} to standardize label ordering.}
#'
#' @importFrom mclust Mclust mclustBIC
#' @export
#'
#' @examples
#' # Clear bimodal structure: 5 low-entropy + 5 high-entropy sites.
#' set.seed(42)
#' entropies <- c(rnorm(5, mean = 0.1, sd = 0.01),
#'                rnorm(5, mean = 1.5, sd = 0.1))
#' result <- cluster_sites_by_entropy(entropies, nr = 100,
#'                                    removez = FALSE, removesngl = FALSE)
#' print(result$DataFrame)
#'
#' # Single-row edge case: class = 1 assigned directly.
#' res1 <- cluster_sites_by_entropy(0.35, nr = 50, removesngl = FALSE)
#' print(res1$DataFrame)
#'
#' # All-identical edge case: class = 999 (sentinel, one undifferentiated group).
#' res2 <- cluster_sites_by_entropy(c(0.35, 0.35, 0.35), nr = 50,
#'                                  removesngl = FALSE)
#' print(res2$DataFrame)
cluster_sites_by_entropy <- function(entropies,
                                     nr,
                                     nsites     = length(entropies),
                                     precision  = 6L,
                                     removez    = TRUE,
                                     removesngl = TRUE,
                                     transfr    = NULL,
                                     verbose    = FALSE,
                                     ...) {

  # --- Helper: empty result with consistent schema ---------------------------
  # Used for the no-rows-remaining case and genuine Mclust failure.
  # The class column is always a zero-length integer vector so downstream
  # consumers can rely on a consistent DataFrame schema and only need to
  # guard on nrow(df) > 0, not also on "class" %in% names(df).
  empty_result <- function(df) {
    df$class <- integer(0L)
    list(FitObject = list(classification = integer(0L)),
         DataFrame  = df)
  }

  # --- 1. Validate nsites ----------------------------------------------------
  if (length(entropies) != nsites) {
    if (verbose)
      warning("Length of entropies (", length(entropies),
              ") does not match nsites (", nsites,
              "). Using length(entropies) as nsites.")
    nsites <- length(entropies)
  }

  # --- 2. Build and sort data frame ------------------------------------------
  df <- data.frame(sites = seq_len(nsites), entropies = entropies)
  df <- df[order(df$entropies), , drop = FALSE]

  # --- 3. Remove invariant sites (entropy ~ 0) --------------------------------
  # Tolerance-based: absorbs floating-point near-zeros that behave as invariant
  # sites (e.g., entropy = 1.2e-10 from rounding in calculate_entropy).
  if (removez) {
    df <- df[df$entropies > 1e-9, , drop = FALSE]
  }

  # --- 4. Remove singleton sites ---------------------------------------------
  # Singleton entropy = H when exactly 1 sequence out of nr differs.
  # Tolerance-based comparison avoids floating-point false negatives.
  if (removesngl) {
    if (missing(nr))
      stop("`nr` must be provided when `removesngl = TRUE`.")
    p1   <- 1 / nr
    p2   <- (nr - 1) / nr
    sngl <- -(p1 * log2(p1) + p2 * log2(p2))
    tol  <- 10^(-precision)
    df   <- df[abs(df$entropies - sngl) > tol, , drop = FALSE]
  }

  # --- 5. Apply optional transformation --------------------------------------
  if (!is.null(transfr)) {
    if (is.function(transfr)) {
      df$entropies <- transfr(df$entropies)
    } else if (inherits(transfr, "transform")) {
      df$entropies <- transfr$transform(df$entropies)
    } else {
      stop("`transfr` must be a function or an object of class `transform`.")
    }
  }

  # --- 6. Strip non-finite values --------------------------------------------
  # Inf / NaN / NA cause Mclust to fail with an uninformative error.
  # Arises when a transformation amplifies near-zero entropies.
  bad <- !is.finite(df$entropies)
  if (any(bad)) {
    if (verbose)
      warning(sum(bad), " non-finite entropy value(s) removed before clustering.")
    df <- df[!bad, , drop = FALSE]
  }

  # --- 7. Edge case: no rows remaining ---------------------------------------
  if (nrow(df) == 0L) {
    if (verbose)
      warning("No sites remaining after filtering; returning empty DataFrame.")
    return(empty_result(df))
  }

  # --- 8. Edge case: single row ----------------------------------------------
  # Mclust is undefined on 1 observation. class = 1 is the only valid answer.
  if (nrow(df) == 1L) {
    df$class <- 1L
    return(list(FitObject = list(classification = 1L),
                DataFrame  = df))
  }

  # --- 9. Edge case: all remaining entropies identical -----------------------
  # Round before uniqueness check to absorb floating-point noise.
  # Sentinel class 999: one undifferentiated group — structurally valid but
  # uninformative for variant detection. Preserved as-is by
  # relabel_entropy_classes so the sentinel meaning is never lost.
  if (length(unique(round(df$entropies, precision))) == 1L) {
    df$class <- 999L
    return(list(FitObject = list(classification = rep(999L, nrow(df))),
                DataFrame  = df))
  }

  # --- 10. Cluster with Mclust -----------------------------------------------
  # mclustBIC is imported explicitly via @importFrom mclust Mclust mclustBIC.
  # This resolves the "could not find function 'mclustBIC'" error that occurs
  # when mclust is in Imports but not attached: in certain mclust versions,
  # Mclust calls mclustBIC via character-string dispatch in a context that
  # searches the calling package's namespace. Importing mclustBIC makes it
  # findable there without any library(mclust) call by the user.
  #
  # suppressWarnings: Mclust emits EM convergence warnings on difficult data;
  # these flood the console during simulation studies with many calls.
  #
  # On genuine failure after the namespace fix (truly degenerate or pathological
  # input), tryCatch returns NULL and step 11 returns an empty DataFrame —
  # treating the partition as uninformative, which is scientifically correct.
  
  # fit <- tryCatch({
  #   suppressWarnings(mclust::Mclust(df$entropies, ...))
  # }, error = function(e) {
  #   if (verbose) warning("Mclust error: ", conditionMessage(e))
  #   NULL
  # })
  
  fit <- tryCatch({
    suppressWarnings({
      utils::capture.output(
        fit_inner <- mclust::Mclust(df$entropies, ...)
      )
      fit_inner
    })
  }, error = function(e) {
    if (verbose) warning("Mclust error: ", conditionMessage(e))
    NULL
  })

  # --- 11. Handle genuine Mclust failure -------------------------------------
  # NULL fit = genuine failure. After the @importFrom fix the namespace error
  # no longer triggers, so reaching here means truly degenerate input data.
  # Return empty DataFrame — treat partition as uninformative.
  if (is.null(fit)) {
    if (verbose)
      warning("Mclust returned NULL; treating partition as uninformative.")
    return(empty_result(df))
  }

  # --- 12. Safety check on classification length ----------------------------
  # fit$classification is a named integer vector of length nrow(df) in all
  # modern mclust versions. This guard is a safety net only.
  if (is.null(fit$classification) ||
      length(fit$classification) != nrow(df)) {
    if (verbose)
      warning("Unexpected classification vector; treating partition as uninformative.")
    return(empty_result(df))
  }

  # --- 13. Attach classification — guaranteed present -----------------------
  # Raw Mclust labels (arbitrary integers 1..G). Call relabel_entropy_classes
  # on the returned DataFrame to standardize so that class 1 = highest entropy.
  df$class <- fit$classification

  list(FitObject = fit, DataFrame = df)
}
