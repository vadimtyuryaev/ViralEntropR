#' @title Cluster Sequence Sites by Entropy
#' @description Identifies groups of sites with similar variability using
#'   Gaussian Mixture Models (GMM).
#'
#' @details
#' Groups sites based on their Shannon entropy values via
#' \code{\link[mclust]{Mclust}}. Preprocessing steps remove invariant sites
#' (entropy = 0) and/or singleton sites (entropy corresponding to exactly one
#' differing sequence across \code{nr} rows). The returned \code{DataFrame}
#' always contains a \code{class} column — this is guaranteed regardless of
#' Mclust's behaviour.
#'
#' Class assignment rules:
#' \itemize{
#'   \item Normal case: Mclust class labels (integers 1, 2, ..., G).
#'   \item All entropies identical after filtering: class \code{999} for all sites.
#'   \item Single site remaining: class \code{1}.
#'   \item Mclust returns NULL or G = 1 (trivial solution): class \code{1} for all sites.
#'   \item No sites remaining: empty DataFrame with no \code{class} column.
#' }
#'
#' @param entropies Numeric vector of entropy values, one per site.
#' @param nr Integer. Total number of sequences used to compute the entropies.
#'   Required when \code{removesngl = TRUE}.
#' @param nsites Integer. Expected number of sites. If it mismatches
#'   \code{length(entropies)}, the actual length is used with a warning.
#'   Default is \code{length(entropies)}.
#' @param precision Integer. Decimal places for rounding during singleton
#'   threshold comparison. Default is \code{6}.
#' @param removez Logical. If \code{TRUE}, removes sites with entropy = 0
#'   (invariant sites). Default is \code{TRUE}.
#' @param removesngl Logical. If \code{TRUE}, removes sites whose entropy equals
#'   the singleton value (one differing sequence out of \code{nr}). Default is
#'   \code{TRUE}.
#' @param transfr A function, or an object of class \code{transform} with a
#'   \code{$transform()} method, applied to entropies before clustering. Default
#'   is \code{NULL} (no transformation).
#' @param verbose Logical. If \code{TRUE}, emits diagnostic warnings for
#'   non-fatal events (empty partitions, Mclust fallbacks, etc.). Default is
#'   \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[mclust]{Mclust}}.
#'
#' @return A named list with two elements:
#' \item{FitObject}{The raw \code{Mclust} result, or a minimal
#'   \code{list(classification = <vector>)} when clustering was bypassed.}
#' \item{DataFrame}{A data frame with columns \code{sites} (original site
#'   indices), \code{entropies} (values after any transformation), and
#'   \code{class} (GMM cluster label). The \code{class} column is \strong{always}
#'   present except when no sites survive filtering.}
#'
#' @importFrom mclust Mclust
#'
#' @examples
#' # Clear bimodal structure: 5 low-entropy + 5 high-entropy sites.
#' # removesngl = FALSE avoids stochastic filtering of the random values.
#' set.seed(42)
#' entropies = c(rnorm(5, mean = 0.1, sd = 0.01),
#'                rnorm(5, mean = 1.5, sd = 0.1))
#' result = cluster_sites_by_entropy(entropies, nr = 100,
#'                                    removez = FALSE, removesngl = FALSE)
#' print(result$DataFrame)           # always has 'class' column
#' print(result$FitObject$classification)
#'
#' @export
cluster_sites_by_entropy = function(entropies,
                                     nr,
                                     nsites    = length(entropies),
                                     precision = 6L,
                                     removez   = TRUE,
                                     removesngl = TRUE,
                                     transfr   = NULL,
                                     verbose   = FALSE,
                                     ...) {
  
  # --- 1. Validate nsites ----------------------------------------------------
  if (length(entropies) != nsites) {
    warning("Length of entropies (", length(entropies),
            ") does not match nsites (", nsites,
            "). Using length(entropies) as nsites.")
    nsites = length(entropies)
  }
  
  # --- 2. Build and sort data frame ------------------------------------------
  df = data.frame(sites = seq_len(nsites), entropies = entropies)
  df = df[order(df$entropies), ]
  
  # --- 3. Remove zeros -------------------------------------------------------
  # Use a small tolerance rather than exact equality to catch floating-point
  # near-zeros that behave as invariant sites.
  if (removez) {
    df = df[df$entropies > 1e-9, , drop = FALSE]
  }
  
  # --- 4. Remove singletons --------------------------------------------------
  # Singleton entropy = H when exactly 1 sequence differs from the other nr-1.
  # Use tolerance-based comparison to avoid floating-point false negatives.
  if (removesngl) {
    if (missing(nr)) {
      stop("`nr` must be provided when `removesngl = TRUE`.")
    }
    p1   = 1 / nr
    p2   = (nr - 1) / nr
    sngl = -(p1 * log2(p1) + p2 * log2(p2))
    tol  = 10^(-precision)
    df   = df[abs(df$entropies - sngl) > tol, , drop = FALSE]
  }
  
  # --- 5. Apply optional transformation --------------------------------------
  if (!is.null(transfr)) {
    if (is.function(transfr)) {
      df$entropies = transfr(df$entropies)
    } else if (inherits(transfr, "transform")) {
      df$entropies = transfr$transform(df$entropies)
    } else {
      stop("Please provide a function or an object of class `transform` for `transfr`.")
    }
  }
  
  # --- 6. Strip non-finite values --------------------------------------------
  # Inf / NaN / NA crash Mclust or produce garbage. Remove before any further
  # checks. Arises legitimately when a site has only one observed amino acid
  # in a small partition (log(0) = -Inf in the entropy formula).
  bad = !is.finite(df$entropies)
  if (any(bad)) {
    if (verbose) warning(sum(bad), " non-finite entropy value(s) removed.")
    df = df[!bad, , drop = FALSE]
  }
  
  # --- 7. Edge case: no rows remaining ---------------------------------------
  # Return early with no class column — there is nothing to classify.
  if (nrow(df) == 0L) {
    if (verbose)
      warning("No sites remaining after filtering; returning empty DataFrame.")
    return(list(FitObject = list(classification = numeric(0L)),
                DataFrame  = df))
  }
  
  # --- 8. Edge case: single row ----------------------------------------------
  # Mclust is undefined on 1 observation. Assign class 1 directly.
  if (nrow(df) == 1L) {
    df$class = 1L
    return(list(FitObject = list(classification = 1L),
                DataFrame  = df))
  }
  
  # --- 9. Edge case: all remaining entropies identical -----------------------
  # Round before uniqueness check to absorb floating-point noise.
  # Use sentinel class 999 — convention inherited from original entropy_clust.
  if (length(unique(round(df$entropies, precision))) == 1L) {
    df$class = 999L
    return(list(FitObject = list(classification = rep(999L, nrow(df))),
                DataFrame  = df))
  }
  
  # --- 10. Cluster with Mclust -----------------------------------------------
  # suppressWarnings: Mclust emits EM convergence warnings on difficult data.
  # These flood the console at 1000s of calls in simulation studies.
  fit = tryCatch({
    suppressWarnings(mclust::Mclust(df$entropies, verbose = verbose, ...))
  }, error = function(e) {
    if (verbose) warning("Mclust error: ", conditionMessage(e))
    NULL
  })
  
  # --- 11. Normalise Mclust output to a guaranteed classification vector -----
  # This is the critical fix vs entropy_clust_modified:
  # we ALWAYS produce a valid classification vector of length nrow(df) before
  # the cbind, so the class column is unconditionally present — mirroring the
  # behaviour of the original entropy_clust.
  #
  # Cases handled:
  #   NULL fit           : Mclust could not fit any model at all.
  #   NULL / length-0    : G=1 trivial solution — Mclust considers the problem
  #                        solved with one component and omits the vector.
  #   Length mismatch    : safety net; should not occur after the guards above.
  #
  # In all fallback cases, class 1 is the neutral assignment (one group).
  if (is.null(fit)) {
    if (verbose) warning("Mclust returned NULL; assigning all sites to class 1.")
    classification = rep(1L, nrow(df))
    fit = list(classification = classification)
  } else if (is.null(fit$classification) ||
             length(fit$classification) == 0L) {
    if (verbose) warning("Mclust G=1 solution; assigning all sites to class 1.")
    classification = rep(1L, nrow(df))
    fit$classification = classification
  } else if (length(fit$classification) != nrow(df)) {
    if (verbose) warning("Classification length mismatch; assigning all sites to class 1.")
    classification = rep(1L, nrow(df))
    fit$classification = classification
  } else {
    classification = fit$classification
  }
  
  # --- 12. Attach classification — always present (mirrors entropy_clust) ----
  # cbind preserves row order; class labels are integers matching Mclust groups.
  df = cbind(df, class = classification)
  
  list(FitObject = fit, DataFrame = df)
}