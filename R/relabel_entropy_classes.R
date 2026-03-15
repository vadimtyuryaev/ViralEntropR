#' @title Relabel Entropy Classes
#' @description Relabels GMM cluster labels so that class 1 is always the
#'   highest entropy group, class 2 the second highest, and so on.
#' @export
#'
#' @param df A data frame containing clustering results. Must have columns
#'   \code{class} and \code{entropies}.
#'
#' @return The data frame with a relabeled \code{class} column (integer) and a
#'   \code{num_classes} column. Returned unchanged in four no-op cases: missing
#'   columns, zero rows, all class \code{999} (sentinel), or a single class.
#'
#' @importFrom stats aggregate setNames
#'
#' @examples
#' df <- data.frame(
#'   sites     = 1:6,
#'   entropies = c(0.1, 0.1, 0.5, 0.5, 1.2, 1.3),
#'   class     = c(1, 1, 2, 2, 3, 3)
#' )
#' # Class 3 (highest entropy) -> 1, Class 1 (lowest entropy) -> 3
#' relabel_entropy_classes(df)
relabel_entropy_classes <- function(df) {

  # --- Guards: cases where relabeling is undefined or a no-op ---------------
  if (!all(c("class", "entropies") %in% names(df))) {
    warning("relabel_entropy_classes: 'class' or 'entropies' column missing ",
            "- returning input unchanged.")
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
  class_stats      <- aggregate(entropies ~ class, data = df, FUN = mean)
  class_stats$rank <- rank(-class_stats$entropies, ties.method = "first")
  mapping          <- setNames(class_stats$rank, class_stats$class)
  # Coerce to integer (not factor) for a consistent class column type across
  # all return paths. as.numeric(factor) returns level codes which happen to
  # equal the label values here, but plain integer is unambiguous and safe
  # for all downstream comparisons (==, %in%, is.integer).
  df$class         <- as.integer(mapping[as.character(df$class)])
  df$num_classes   <- length(unique(df$class))

  df
}
