#' @title Partition Data into Time Windows
#' @description Splits a data frame of time-stamped sequences into discrete time
#'   windows, computes per-site entropy and Mclust clustering for each window.
#'
#' @details
#' Three windowing strategies are supported via \code{window_type}:
#' \itemize{
#'   \item \strong{1 — Cumulative}: Start is fixed; end expands by
#'     \code{window_length} months each step.
#'   \item \strong{2 — Sliding}: A window of \code{window_length} months slides
#'     one month at a time. Produces
#'     \code{total_months - window_length + 1} chunks.
#'   \item \strong{3 — Non-overlapping (Jumping)}: Consecutive non-overlapping
#'     windows of \code{window_length} months.
#' }
#'
#' @param data Data frame. Must contain a \code{Date} column coercible to
#'   \code{Date}, plus numeric site columns.
#' @param n_sites Integer. Number of site columns (assumed to be columns
#'   \code{1:n_sites} of \code{data}).
#' @param window_length Integer. Window duration in months. Default is \code{1}.
#' @param window_type Integer. Windowing strategy: \code{1} = Cumulative,
#'   \code{2} = Sliding, \code{3} = Non-overlapping. Default is \code{3}.
#' @param start_date Date or character. Window start. Defaults to earliest date
#'   in \code{data}.
#' @param end_date Date or character. Window end (cutoff). Defaults to latest
#'   date in \code{data}.
#' @param date_format Character. Format string for date labels. Default is
#'   \code{"\%b-\%Y"}.
#' @param verbose Logical. If \code{TRUE}, prints a progress bar and completion
#'   summary. Default is \code{FALSE}.
#' @param ... Additional arguments passed to
#'   \code{\link{cluster_sites_by_entropy}}.
#'
#' @return A named list:
#' \item{Partitions}{List of data frame chunks, one per window.}
#' \item{Entropies}{List of numeric entropy vectors, one per window.}
#' \item{Clusters}{List of clustering results from
#'   \code{\link{cluster_sites_by_entropy}}, one per window.}
#' \item{Max_Entropy}{Numeric vector. Maximum cluster class label per window
#'   (equals number of clusters found; \code{NA} if window was empty).}
#' \item{Dates_Labels}{Character vector of window label strings.}
#' \item{N_partitions}{Integer. Total number of windows.}
#'
#' @importFrom lubridate %m+% interval period
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' dates <- seq(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
#' df <- data.frame(
#'   Date = rep(dates, each = 5),
#'   s1 = 1L,
#'   s2 = c(rep(1L, 15), rep(2L, 15))
#' )
#' res <- partition_time_windows(df, n_sites = 2, window_length = 2,
#'                               window_type = 3, verbose = TRUE)
#' print(res$Dates_Labels)
#' print(res$Max_Entropy)
#'
#' @export
partition_time_windows <- function(data,
                                   n_sites,
                                   window_length = 1,
                                   window_type   = 3,
                                   start_date    = NULL,
                                   end_date      = NULL,
                                   date_format   = "%b-%Y",
                                   verbose       = FALSE,
                                   ...) {
  
  # --- 1. Validation ---------------------------------------------------------
  if (!"Date" %in% names(data)) {
    stop("The column with dates is missing or is not properly named.")
  }
  
  data <- data[order(data$Date), ]
  
  if (is.null(start_date)) start_date <- min(data$Date)
  if (is.null(end_date))   end_date   <- max(data$Date)
  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)
  
  start_date_dspl <- format(start_date, format = date_format)
  
  # --- 2. Compute n_chunks per window type -----------------------------------
  total_months <- lubridate::interval(start_date, end_date) %/%
    lubridate::period(1, "months")
  
  if (window_type == 1) {
    n_months <- total_months
    if (n_months %% window_length != 0)
      n_months <- n_months - (n_months %% window_length)
    n_chunks <- n_months / window_length
    
  } else if (window_type == 2) {
    n_chunks <- total_months - window_length + 1
    
  } else {
    n_chunks <- total_months %/% window_length
  }
  
  # --- 3. Per-window processing ----------------------------------------------
  extra_args <- list(...)
  
  if (verbose) {
    window_type_name <- switch(as.character(window_type),
                               "1" = "cumulative",
                               "2" = "sliding",
                               "3" = "non-overlapping",
                               "unknown")
    message(sprintf("Partitioning %d %s window%s ...",
                    n_chunks, window_type_name,
                    if (n_chunks == 1L) "" else "s"))
    pb <- utils::txtProgressBar(min = 0, max = n_chunks, style = 3, width = 50)
  }
  
  raw <- vector("list", n_chunks)
  
  for (i in seq_len(n_chunks)) {
    
    # Window bounds
    if (window_type == 1) {
      win_start <- start_date
      win_end   <- lubridate::`%m+%`(start_date,
                                     lubridate::period(window_length * i, "months"))
      
    } else if (window_type == 2) {
      win_start <- lubridate::`%m+%`(start_date,
                                     lubridate::period(i - 1, "months"))
      win_end   <- lubridate::`%m+%`(win_start,
                                     lubridate::period(window_length, "months"))
      
    } else {
      win_start <- lubridate::`%m+%`(start_date,
                                     lubridate::period(window_length * (i - 1), "months"))
      win_end   <- lubridate::`%m+%`(win_start,
                                     lubridate::period(window_length, "months"))
    }
    
    chunk <- data[data$Date >= win_start & data$Date < win_end &
                    data$Date <= end_date, , drop = FALSE]
    
    # Date label
    if (window_type == 1) {
      dates_disp <- if (nrow(chunk) > 0)
        format(as.Date(max(chunk$Date)), format = date_format)
      else
        format(win_end, format = date_format)
      label <- paste(start_date_dspl, dates_disp, sep = " - ")
    } else {
      label <- paste(
        format(win_start, format = date_format),
        format(lubridate::`%m+%`(win_end, lubridate::period(-1, "months")),
               format = date_format),
        sep = " - "
      )
    }
    
    # Entropy + clustering
    if (nrow(chunk) > 0) {
      entrp_all <- apply(chunk[, seq_len(n_sites), drop = FALSE], 2,
                         calculate_entropy)
      all_clust <- do.call(cluster_sites_by_entropy,
                           c(list(entrp_all,
                                  nr     = nrow(chunk),
                                  nsites = n_sites),
                             extra_args))
      
      # max_ent = max raw Mclust class label = number of clusters found.
      # Computed before any downstream relabeling.
      # relabel_entropy_classes belongs in the downstream consumer.
      max_ent <- if (nrow(all_clust$DataFrame) > 0 &&
                     "class" %in% names(all_clust$DataFrame))
        max(as.numeric(all_clust$DataFrame$class), na.rm = TRUE)
      else
        NA_integer_
      
      raw[[i]] <- list(chunk   = chunk,
                       entrp   = entrp_all,
                       clust   = all_clust,
                       max_ent = max_ent,
                       label   = label)
    } else {
      raw[[i]] <- list(chunk   = chunk,
                       entrp   = rep(0, n_sites),
                       clust   = list(FitObject = list(classification = numeric(0)),
                                      DataFrame = data.frame()),
                       max_ent = NA_integer_,
                       label   = label)
    }
    
    if (verbose) utils::setTxtProgressBar(pb, i)
  }
  
  if (verbose) {
    close(pb)
    message(sprintf("Partitioning complete: %d partition%s generated (%s to %s).",
                    n_chunks,
                    if (n_chunks == 1L) "" else "s",
                    format(start_date, date_format),
                    format(end_date, date_format)))
  }
  
  # --- 4. Unpack results -----------------------------------------------------
  list(
    Partitions   = lapply(raw, `[[`, "chunk"),
    Entropies    = lapply(raw, `[[`, "entrp"),
    Clusters     = lapply(raw, `[[`, "clust"),
    Max_Entropy  = vapply(raw, `[[`, numeric(1), "max_ent"),
    Dates_Labels = vapply(raw, `[[`, character(1), "label"),
    N_partitions = n_chunks
  )
}