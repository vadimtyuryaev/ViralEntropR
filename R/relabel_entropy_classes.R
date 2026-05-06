#' @title Relabel Entropy Classes
#' @description Relabels GMM cluster labels so that class 1 is always the
#'   highest-entropy group, class 2 the second highest, and so on.
#'
#' @details
#' For univariate input (the package's typical use case), \code{\link[mclust]{Mclust}}
#' orders cluster components by increasing mean. Applied to per-site Shannon
#' entropy values, this places the highest-entropy group at the highest
#' label number (label \code{G} for \code{G} fitted components). This
#' function flips that convention so that label \code{1} always denotes
#' the highest-entropy group, which is more natural for downstream
#' filtering ("class 1 = top") and visual labeling. Cluster identities
#' and means are unchanged; only the integer labels are remapped.
#'
#' \strong{Sentinel preservation.} The class label \code{999L} marks
#' undifferentiated groups (all entropies equal in
#' \code{\link{cluster_sites_by_entropy}}) and is never relabeled. If a
#' data frame contains a mixture of real GMM labels and one or more
#' \code{999} entries, only the real labels are ranked and overwritten;
#' the \code{999} rows pass through unchanged.
#'
#' \strong{No-op return paths.} The input is returned unchanged (with at
#' most a \code{num_classes} column added) in four situations:
#' \itemize{
#'   \item Missing required columns (\code{class} or \code{entropies}) —
#'     warns and returns input.
#'   \item Zero rows.
#'   \item All rows already carry the sentinel \code{999L}.
#'   \item Only one class label is present (relabeling is a no-op).
#' }
#'
#' @param df A data frame containing clustering results. Must have columns
#'   \code{class} (integer cluster labels from
#'   \code{\link{cluster_sites_by_entropy}}) and \code{entropies}
#'   (numeric Shannon entropy values per row).
#'
#' @return The data frame with a relabeled \code{class} column (integer)
#'   and an added \code{num_classes} column reporting the count of
#'   distinct labels present after relabeling.
#'
#' @seealso \code{\link{cluster_sites_by_entropy}} for the upstream
#'   clustering step; \code{\link{plot_entropy_trajectories}} and
#'   \code{\link{plot_site_class_trajectory}} for downstream consumers
#'   that rely on the relabeled class convention;
#'   \code{\link{partition_time_windows}} which calls
#'   \code{cluster_sites_by_entropy} per window.
#'
#' @importFrom stats aggregate setNames
#' @export
#'
#' @examples
#' # Standard case: three classes, ranked by mean entropy.
#' df <- data.frame(
#'   sites     = 1:6,
#'   entropies = c(0.1, 0.1, 0.5, 0.5, 1.2, 1.3),
#'   class     = c(1L, 1L, 2L, 2L, 3L, 3L)
#' )
#' relabel_entropy_classes(df)
#' # Class 3 (highest mean entropy 1.25) -> 1; class 2 (0.5) -> 2;
#' # class 1 (0.1) -> 3.
#'
#' # Sentinel preservation: 999 rows pass through unchanged even when
#' # mixed with real classes.
#' df_mixed <- data.frame(
#'   sites     = 1:4,
#'   entropies = c(0.5, 1.2, 0.3, 0.4),
#'   class     = c(1L, 2L, 999L, 1L)
#' )
#' relabel_entropy_classes(df_mixed)
#' # Class 2 (highest) -> 1; class 1 -> 2; class 999 stays 999.
relabel_entropy_classes <- function(df) {
  
  # --- Guards: cases where relabeling is undefined or a no-op ---------------
  if (!all(c("class", "entropies") %in% names(df))) {
    warning("relabel_entropy_classes: 'class' or 'entropies' column missing ",
            "- returning input unchanged.", call. = FALSE)
    return(df)
  }
  
  if (nrow(df) == 0L)
    return(df)
  
  # Sentinel: all entropies identical, one undifferentiated group.
  # Relabeling would overwrite 999 -> 1 and destroy the sentinel.
  if (all(df$class == 999L)) {
    df$num_classes <- 1L
    return(df)
  }
  
  # Single class present — relabeling one group is a no-op.
  if (length(unique(df$class)) == 1L) {
    df$num_classes <- 1L
    return(df)
  }
  
  # --- Relabel: rank classes by mean entropy descending ---------------------
  # Sentinel rows (999) are excluded from ranking and pass through
  # unchanged. Only the real GMM labels are reordered.
  is_sentinel <- df$class == 999L
  
  if (any(is_sentinel) && length(unique(df$class[!is_sentinel])) <= 1L) {
    # Mixed sentinel + at most one real class: no real ranking to perform.
    df$num_classes <- length(unique(df$class))
    return(df)
  }
  
  df_real          <- if (any(is_sentinel)) df[!is_sentinel, , drop = FALSE] else df
  class_stats      <- aggregate(entropies ~ class, data = df_real, FUN = mean)
  class_stats$rank <- rank(-class_stats$entropies, ties.method = "first")
  mapping          <- setNames(class_stats$rank, class_stats$class)
  
  # Coerce to integer (not factor) for a consistent class column type across
  # all return paths. as.numeric(factor) returns level codes which happen to
  # equal the label values here, but plain integer is unambiguous and safe
  # for all downstream comparisons (==, %in%, is.integer).
  if (any(is_sentinel)) {
    df$class[!is_sentinel] <- as.integer(mapping[as.character(df$class[!is_sentinel])])
    # Sentinel rows keep their 999 — no assignment needed.
  } else {
    df$class <- as.integer(mapping[as.character(df$class)])
  }
  
  df$num_classes <- length(unique(df$class))
  
  df
}