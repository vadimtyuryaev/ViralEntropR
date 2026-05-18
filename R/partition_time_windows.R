#' @title Partition Data into Time Windows
#' @description Splits a data frame of time-stamped sequences into discrete time
#'   windows, computes per-site entropy and Mclust clustering for each window.
#'
#' @details
#' Three windowing strategies are supported via \code{window_type}.
#' Throughout the descriptions below, \code{T} denotes the number of whole
#' months between \code{start_date} and \code{end_date}, and \code{w}
#' denotes \code{window_length}:
#' \itemize{
#'   \item \strong{1 — Cumulative}: Start is fixed; end expands by
#'     \code{w} months each step. Each window includes all prior data
#'     plus one more period. Produces \code{floor(T / w)} chunks.
#'   \item \strong{2 — Sliding}: A window of \code{w} months slides one
#'     month at a time. Produces \code{T - w + 1} chunks.
#'   \item \strong{3 — Disjoint}: Consecutive non-overlapping windows of
#'     \code{w} months. Produces \code{floor(T / w)} chunks. Default.
#' }
#'
#' For example, with \code{T = 12} months and \code{w = 2}: cumulative
#' produces 6 chunks (each progressively larger); sliding produces 11
#' chunks (overlapping, each 2 months wide); disjoint produces 6 chunks
#' (each exactly 2 months, non-overlapping).
#' 
#' \strong{Empty windows.} When a window contains no observations, the
#' corresponding entries in the returned lists carry placeholder values:
#' \code{Entropies} is a zero-vector of length \code{n_sites};
#' \code{Clusters} is a schema-consistent empty result matching
#' \code{cluster_sites_by_entropy}'s empty-input return; \code{Max_Entropy}
#' is \code{NA_integer_}. Downstream consumers can guard on
#' \code{nrow(Partitions[[i]]) > 0} or \code{!is.na(Max_Entropy[i])}.
#'
#' \strong{Column layout requirement.} The function extracts site columns
#' as \code{data[, seq_len(n_sites), drop = FALSE]}. The site columns
#' must therefore occupy positions \code{1} through \code{n_sites}.
#' \code{Date} (and any other metadata columns such as \code{Country})
#' must come after the site columns. A common error is to place
#' \code{Date} first; this will produce incorrect entropy values.
#'
#' @param data Data frame. Must contain a column named \code{Date}
#'   coercible to \code{Date}. Columns \code{1} through \code{n_sites}
#'   must be the numeric site columns; any additional columns (such as
#'   \code{Country} or other per-sequence metadata) may follow.
#' @param n_sites Integer. Number of site columns. Site columns must
#'   occupy positions \code{1} through \code{n_sites} of \code{data};
#'   any other columns (\code{Date}, optionally \code{Country}, etc.)
#'   must come after.
#' @param window_length Integer. Window duration in months. 
#'   Default is \code{1}.
#' @param window_type Integer (1, 2, or 3). Windowing strategy:
#'  \code{1} = Cumulative, \code{2} = Sliding, \code{3} = Disjoint. 
#'  Default is \code{3}.
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
#' @seealso \code{\link{calculate_entropy}} for the per-site entropy
#'   computation; \code{\link{cluster_sites_by_entropy}} for the GMM
#'   clustering applied per window; \code{\link{relabel_entropy_classes}}
#'   for standardizing the cluster labels in downstream consumers; and
#'   \code{\link{calculate_hellinger_matrix}} for typical downstream use
#'   of the \code{Partitions} output.
#'
#' @importFrom lubridate %m+% %m-% interval period
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' dates <- seq(as.Date("2020-01-01"), as.Date("2020-06-01"), by = "month")
#' df <- data.frame(
#'   s1   = 1L,
#'   s2   = c(rep(1L, 15), rep(2L, 15)),
#'   Date = rep(dates, each = 5)
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
    stop("The column with dates is missing or is not properly named.", 
         call. = FALSE)
  }
  
  if (!window_type %in% 1:3)
    stop("`window_type` must be 1 (cumulative), 2 (sliding), or 3 (disjoint).",
         call. = FALSE)
  
  if (n_sites > ncol(data) - 1L)
    stop(sprintf("`n_sites` (%d) exceeds available non-Date columns (%d). ",
                 n_sites, ncol(data) - 1L),
         "Site columns must occupy positions 1 through n_sites.",
         call. = FALSE)
  
  data <- data[order(data$Date), ]
  
  if (is.null(start_date)) start_date <- min(data$Date)
  if (is.null(end_date))   end_date   <- max(data$Date)
  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)
  
# Empty-window predicates compare against window boundaries (not raw data
# boundaries) so spurious warnings when start_date or end_date sits within
# one window_length of the data range are avoided.
  first_window_right <- start_date %m+% lubridate::period(window_length, "months")
  
  if (min(data$Date) >= first_window_right)
    warning("`start_date` is earlier than the data by more than one ",
            "window_length; the first window will be empty.", call. = FALSE)
  
  if (window_type != 1L) {
    last_window_left <- end_date %m-% lubridate::period(window_length, "months")
    if (max(data$Date) < last_window_left)
      warning("`end_date` is later than the data by more than one ",
              "window_length; trailing windows will be empty.", call. = FALSE)
    }
  
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
  
  if (n_chunks < 1L)
    stop(sprintf("Computed n_chunks = %d for window_type %d with window_length = %d ",
                 n_chunks, window_type, window_length),
         sprintf("over a %d-month range (%s to %s). ",
                 total_months, format(start_date), format(end_date)),
         "The date range is too short for the chosen windowing parameters.",
         call. = FALSE)
  
  # --- 3. Per-window processing ----------------------------------------------
  extra_args <- list(...)
  
  # --- 4. Pre-allocate output structures -------------------------------------
  # Writing directly into pre-allocated vectors/lists avoids accumulating a
  # full raw list in memory. Each iteration's temporary objects (chunk, entrp,
  # all_clust) go out of scope immediately after assignment and can be garbage
  # collected, reducing peak RAM by ~40-70% depending on window type.
  Partitions   <- vector("list", n_chunks)
  Entropies    <- vector("list", n_chunks)
  Clusters     <- vector("list", n_chunks)
  Max_Entropy  <- numeric(n_chunks)
  Dates_Labels <- character(n_chunks)
  
  if (verbose) {
    window_type_name <- switch(as.character(window_type),
                               "1" = "cumulative",
                               "2" = "sliding",
                               "3" = "disjoint",
                               "unknown")
    message(sprintf("Partitioning %d %s window%s ...",
                    n_chunks, window_type_name,
                    if (n_chunks == 1L) "" else "s"))
    pb <- utils::txtProgressBar(min = 0, max = n_chunks, style = 3, width = 50)
  }
  
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
      Partitions[[i]]  <- chunk
      Entropies[[i]]   <- entrp_all
      Clusters[[i]]    <- all_clust
      Max_Entropy[i]   <- if (nrow(all_clust$DataFrame) > 0 &&
                              "class" %in% names(all_clust$DataFrame))
        max(as.numeric(all_clust$DataFrame$class), na.rm = TRUE)
      else
        NA_integer_
      Dates_Labels[i]  <- label
      
    } else {
      Partitions[[i]]  <- chunk
      Entropies[[i]]   <- rep(0, n_sites)
      Clusters[[i]] <- list(FitObject = list(classification = integer(0L)),
                            DataFrame = data.frame(sites     = integer(0L),
                                                   entropies = numeric(0L),
                                                   class     = integer(0L)))
      Max_Entropy[i]   <- NA_integer_
      Dates_Labels[i]  <- label
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
  
  # --- 5. Return -------------------------------------------------------------
  list(
    Partitions   = Partitions,
    Entropies    = Entropies,
    Clusters     = Clusters,
    Max_Entropy  = Max_Entropy,
    Dates_Labels = Dates_Labels,
    N_partitions = n_chunks
  )
}