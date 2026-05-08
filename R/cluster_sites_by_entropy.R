#' @title Cluster a Univariate Numeric Vector by Gaussian Mixture Model
#' @description Wraps \code{\link[mclust]{Mclust}} for unsupervised clustering
#'   of a univariate numeric vector, with preprocessing rules and edge-case
#'   handling tailored to per-site Shannon entropy values from viral sequence
#'   data, which is the package's primary use case, but applicable to any
#'   univariate data the user wishes to cluster by GMM.
#'
#' @details
#' In the package's typical use, sites are clustered by their Shannon entropy
#' to identify groups of residue positions with similar variability across
#' a sequence collection. Two preprocessing rules apply when clustering
#' entropies: \code{removez = TRUE} drops invariant sites (entropy = 0), and
#' \code{removesngl = TRUE} drops singleton sites whose entropy corresponds
#' to exactly one differing sequence across \code{nr} rows.
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
#'     These are Mclust's own integer labels, ordered by increasing component
#'     mean (univariate Mclust orders components by mean) — call
#'     \code{relabel_entropy_classes()} on the returned data frame to obtain
#'     application-friendly class labels (highest-entropy class \code{= 1},
#'     lowest-entropy class \code{= G}).
#'     \code{\link{relabel_entropy_classes}} on the returned DataFrame to
#'     standardise so that class 1 = highest-entropy group.
#'   \item \strong{Mclust failure}: empty DataFrame returned (same schema as
#'     the no-rows case), treating the partition as uninformative.
#' }
#'
#' @param entropies Numeric vector to cluster. In the package's primary use
#'   case these are per-site Shannon entropy values, but any univariate
#'   numeric vector is accepted.
#' @param nr Integer. Total number of sequences from which the entropies
#'   were computed. Required only when \code{removesngl = TRUE}.
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
#'   call \code{\link{relabel_entropy_classes}} to standardise label ordering.}
#'
#' @seealso
#' \code{\link{calculate_entropy}} for computing per-site entropy values,
#'   \code{\link{relabel_entropy_classes}} for standardising the returned
#'   class labels, and \code{\link{partition_time_windows}}, which calls
#'   this function on each temporal partition.
#'
#' @importFrom mclust Mclust mclustBIC
#' @export
#'
#' @examples
#' # Clear bimodal structure: 5 low-entropy + 5 high-entropy sites.
#' set.seed(42)
#' entropies <- c(rnorm(5, mean = 0.1, sd = 0.01),
#'                rnorm(5, mean = 1.5, sd = 0.1))
#' result <- cluster_sites_by_entropy(entropies, removez = FALSE,
#'                                    removesngl = FALSE)
#' print(result$DataFrame)
#'
#' # Single-row edge case: class = 1 assigned directly.
#' res1 <- cluster_sites_by_entropy(0.35, removesngl = FALSE)
#' print(res1$DataFrame)
#'
#' # All-identical edge case: class = 999 (sentinel, one undifferentiated group).
#' res2 <- cluster_sites_by_entropy(c(0.35, 0.35, 0.35), removesngl = FALSE)
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

  # --- 3. Remove invariant sites (entropy ~ 0) -------------------------------
  if (removez) {
    df <- df[df$entropies > 1e-9, , drop = FALSE]
  }

  # --- 4. Remove singleton sites ---------------------------------------------
  # Singleton entropy = H when exactly 1 sequence out of nr differs.
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
  # Inf / NaN / NA would otherwise cause Mclust to fail with an opaque error.
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
  # Class 999 is the sentinel for an undifferentiated group; preserved as-is
  # by relabel_entropy_classes so the sentinel meaning is never lost.
  if (length(unique(round(df$entropies, precision))) == 1L) {
    df$class <- 999L
    return(list(FitObject = list(classification = rep(999L, nrow(df))),
                DataFrame  = df))
  }

  # --- 10. Cluster with Mclust -----------------------------------------------
  # suppressWarnings: silence Mclust EM convergence warnings during heavy use.
  # capture.output: silence Mclust's internal cat() during BIC search.
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
  if (is.null(fit)) {
    if (verbose)
      warning("Mclust returned NULL; treating partition as uninformative.")
    return(empty_result(df))
  }

  # --- 12. Safety check on classification length ----------------------------
  if (is.null(fit$classification) ||
      length(fit$classification) != nrow(df)) {
    if (verbose)
      warning("Unexpected classification vector; treating partition as uninformative.")
    return(empty_result(df))
  }

  # --- 13. Attach classification — guaranteed present -----------------------
  # Raw Mclust labels (1..G); call relabel_entropy_classes() to standardise.
  df$class <- fit$classification

  list(FitObject = fit, DataFrame = df)
}
