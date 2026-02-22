#' @title Relabel Entropy Classes
#' @description Standardizes cluster labels so that Class 1 represents the highest entropy group.
#' @export
#' 
#' @details
#' Takes a data frame containing \code{class} and \code{entropy} (or similar) columns.
#' It calculates the mean entropy for each class, ranks them in descending order,
#' and re-assigns class labels: 1 = Highest Entropy, 2 = Second Highest, etc.
#'
#' @param df A data frame containing clustering results. Must have columns \code{class} and \code{entropy}.
#'
#' @return The data frame with a modified \code{class} column (factor).
#' 
#' @importFrom stats aggregate setNames
#' 
#' @examples
#' # 1. Create dummy clustering data
#' df <- data.frame(
#'   sites    = 1:6,
#'   entropies = c(0.1, 0.1, 0.5, 0.5, 1.2, 1.3),
#'   class    = c(1, 1, 2, 2, 3, 3)
#' )
#'
#' print("Original:")
#' print(df)
#'
#' # 2. Relabel: Class 3 (high entropy) -> Class 1, Class 1 (low entropy) -> Class 3
#' df_new <- relabel_entropy_classes(df)
#' print("Relabelled:")
#' print(df_new)
relabel_entropy_classes = function(df) {
  if (!all(c("class", "entropies") %in% names(df))) {
    warning("relabel_entropy_classes: 'class' or 'entropies' column missing - returning input unchanged.")
    return(df)
  }
  
  # Identify unique classes and their mean entropy
  # Using aggregate to be safe against base R behavior
  class_stats = aggregate(entropies ~ class, data = df, FUN = mean)
  
  # Rank classes: Descending order of entropy (1 is highest)
  class_stats$rank = rank(-class_stats$entropies, ties.method = "first")
  
  # Create a mapping vector
  # names = old class, values = new class
  mapping = setNames(class_stats$rank, class_stats$class)
  
  # Apply mapping
  df$class = factor(mapping[as.character(df$class)])
  
  # num_classes: how many distinct entropy groups Mclust identified.
  # Note: maxclass is intentionally NOT added here. After relabeling, class 1
  # is always the highest-entropy group by definition, so a maxclass column
  # would always equal 1 and carry no information. Downstream code should
  # use num_classes to know how many clusters were found.
  df$num_classes <- length(unique(df$class))
  
  return(df)
}