knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  warning  = FALSE,
  message  = FALSE
)

library(ViralEntropR)
library(Biostrings)
library(dplyr)
library(ggplot2)

# # Run once -- downloads ~173 MB and caches locally.
# # Subsequent calls return the cached path instantly.
# fasta_path <- download_sarscov2_data()
# fasta_path   # path to the cached file

# -- TEMPORARY: read directly from data-raw/ until Zenodo DOI is finalized ---
# Once download_sarscov2_data() is operational, replace these two lines with:
#   fasta_path <- download_sarscov2_data()
#   fasta      <- Biostrings::readAAStringSet(fasta_path)
fasta_path <- file.path("data-raw", "sequences.fasta")
fasta      <- Biostrings::readAAStringSet(fasta_path)
# ---------------------------------------------------------------------------

sprintf("Total sequences loaded:  %d", length(fasta))
sprintf("Unique sequences:        %d", length(unique(fasta)))

# Print 5 randomly sampled headers to see the field layout
set.seed(42)
cat(paste(sample(names(fasta), 5), collapse = "\n\n"))


# The sarscov2_variants object is available immediately after library(ViralEntropR).

# --- Overview of all 12 variants -------------------------------------------

data.frame(
  WHO_Label      = unlist(sarscov2_variants$WHO_Label),
  Pango_Lineage  = unlist(sarscov2_variants$Pango_Lineage),
  GISAID_Clade   = sarscov2_variants$GISAID_Clade,
  Nextstrain     = sarscov2_variants$Nextstrain_Clade,
  First_Detected = sarscov2_variants$Date_First_Detected,
  First_US       = sarscov2_variants$Date_First_Detected_US
)

# --- Cross-reference mutation sites between variants -----------------------
# Which VOCs/VOIs share mutation site 501 (N501Y -- key RBD mutation)?

idx_501 <- which(sapply(sarscov2_variants$Mutation_Sites,
                        function(s) 501 %in% s))

cat("Variants carrying a mutation at site 501:\n")
cat(paste(unlist(sarscov2_variants$WHO_Label)[idx_501], collapse = ", "), "\n")

# # To learn more
# ?sarscov2_variants

# sarscov2_sample.fasta.gz contains 100 random sequences from the full NCBI dataset.
# Useful for testing pipeline code without downloading the full 173 MB file.

path_sample  <- system.file("extdata", "sarscov2_sample.fasta.gz",
                     package = "ViralEntropR")
fasta_sample <- Biostrings::readAAStringSet(path_sample)

sprintf("Sample sequences: %d", length(fasta_sample))

# check NAMES structure
sample(fasta_sample@ranges@NAMES,1)

# Countries represented in the sample
sort(table(extract_fasta_countries(fasta_sample, position = 2)$countries),
     decreasing = TRUE)

# # To learn more
# ?sarscov2_sample

sprintf("All sequences are length 1,273: %s",
        all(width(fasta) == 1273L))

# Pass 1: attempt full yyyy-mm-dd extraction
dates_ymd <- extract_fasta_dates(
  fasta,
  custom_pattern = "(?<=\\|)[0-9]{4}-(0?[1-9]|1[0-2])-(0?[1-9]|[12][0-9]|3[01]|00)"
)

sprintf("Sequences with full yyyy-mm-dd dates: %d",
        sum(!is.na(dates_ymd$raw_date_strings)))
sprintf("Sequences missing full date:          %d",
        sum(is.na(dates_ymd$raw_date_strings)))

# Inspect headers of sequences that lack a full yyyy-mm-dd date
if (length(dates_ymd$missing_id) > 0) {
  cat("Sample headers with partial/missing dates:\n")
  cat(paste(head(names(fasta)[dates_ymd$missing_id], 5), collapse = "\n"))
}

# Pass 2: extract yyyy-mm (covers both yyyy-mm-dd and yyyy-mm submissions)
dates_result <- extract_fasta_dates(
  fasta,
  custom_pattern = "(?<=\\|)[0-9]{4}-(0[1-9]|1[0-2])",
  date_format    = "%Y-%m"
)

sprintf("Sequences with parseable yyyy-mm dates: %d", 
        sum(!is.na(dates_result$raw_date_strings)))
sprintf("Sequences still missing a date:         %d",
        sum(is.na(dates_result$raw_date_strings)))

# Sequences with year-only dates -- these cannot be placed in a monthly window
if (length(dates_result$missing_id) > 0) {
  cat("Headers with year-only dates (to be removed):\n")
  cat(paste(head(names(fasta)[dates_result$missing_id], 5), collapse = "\n"))
}

# Remove sequences with year-only dates
if (length(dates_result$missing_id) > 0) {
  fasta        <- fasta[-dates_result$missing_id]
  dates_result <- extract_fasta_dates(
    fasta,
    custom_pattern = "(?<=\\|)[0-9]{4}-(0[1-9]|1[0-2])",
    date_format    = "%Y-%m"
  )
}

cat("Date extraction status after removal:", dates_result$message, "\n")
sprintf("Sequences remaining: %d", length(fasta))

dates_df_before <- data.frame(
  Date = as.Date(dates_result$corrected_dates)
) %>%
  mutate(YearMonth = format(Date, "%Y-%m"))

ggplot(dates_df_before, aes(x = YearMonth)) +
  geom_bar(stat = "count", fill = "steelblue") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    plot.title  = element_text(hjust = 0.5)
  ) +
  labs(
    x     = "Collection Month",
    y     = "Number of Sequences",
    title = "Sequence Count by Collection Month (before ambiguous residue filtering)"
  )

# Extract country field (position 2 = between first and second pipe)
countries_result <- extract_fasta_countries(fasta, position = 2)

sprintf("Sequences with parseable country: %d", 
        sum(!is.na(countries_result$countries)))
sprintf("Sequences with missing country:   %d",
        sum(is.na(countries_result$countries)))

# Inspect headers with missing countries before removing
if (length(countries_result$missing_id) > 0) {
  cat("Sample headers with missing country field:\n")
  cat(paste(head(names(fasta)[countries_result$missing_id], 5), collapse = "\n"))
  
  # Remove sequences with missing countries
  fasta            <- fasta[-countries_result$missing_id]
  countries_result <- extract_fasta_countries(fasta, position = 2)
}

cat("Country extraction status after removal:", countries_result$message, "\n")
sprintf("Sequences remaining: %d", length(fasta))

# Re-extract dates after removing sequences with missing countries
dates_result <- extract_fasta_dates(
  fasta,
  custom_pattern = "(?<=\\|)[0-9]{4}-(0[1-9]|1[0-2])",
  date_format    = "%Y-%m"
)

# Confirm no remaining missing values in either field
sprintf("Missing dates remaining:     %d", sum(is.na(dates_result$corrected_dates)))
sprintf("Missing countries remaining: %d", sum(is.na(countries_result$countries)))

# Country distribution -- sequences per country of collection
countries_result$countries %>%
  table() %>%
  sort(decreasing = TRUE) 

char_mat <- fasta_to_char_matrix(fasta)

sprintf("Character matrix dimensions: %d sequences x %d sites",
        nrow(char_mat), ncol(char_mat))

# Inspect a small region
char_mat[1:5, 1:10]

filtered <- filter_ambiguous_sequences(char_mat, option = 2)

cat(filtered$OriginalDim,     "\n")
cat(filtered$NewDim,          "\n")
cat(filtered$NumberAmbiguous, "\n")
cat(filtered$RangeAmbiguous,  "\n")

# Extract the clean matrix and align metadata to retained rows
clean_char_mat  <- filtered$FilteredMatrix
deleted_rows    <- filtered$DeletedSeqId

# Remove deleted rows from date and country vectors
corrected_dates <- dates_result$corrected_dates[-deleted_rows]
countries       <- countries_result$countries[-deleted_rows]

sprintf("Clean matrix: %d sequences x %d sites",
        nrow(clean_char_mat), ncol(clean_char_mat))

int_mat <- encode_aa_sequence(clean_char_mat)

# Verify only codes 1-20 remain
code_range <- range(int_mat)
sprintf("Encoded value range: %d to %d  (expected: 1 to 20 after filtering)",
        code_range[1], code_range[2])

n_sites <- ncol(int_mat)

# Build the data frame
AL_df           <- as.data.frame(int_mat)
colnames(AL_df) <- as.character(seq_len(n_sites))
AL_df[]         <- lapply(AL_df, as.integer)

# Standardise to first of each month for consistent monthly partitioning
AL_df$Date    <- as.Date(format(corrected_dates, "%Y-%m-01"))
AL_df$Country <- countries

# Sort by date -- required by partition_time_windows()
AL_df         <- AL_df[order(AL_df$Date), ]
rownames(AL_df) <- NULL

sprintf("Final data frame: %d sequences x %d site columns + Date + Country",
        nrow(AL_df), n_sites)

# Preview first few rows -- site columns and metadata
AL_df[1:10, c(1:5, n_sites + 1, n_sites + 2)]

# Country distribution in the clean final dataset
AL_df$Country %>% table() %>% sort(decreasing = TRUE)

ggplot(AL_df, aes(x = format(Date, "%Y-%m"))) +
  geom_bar(stat = "count", fill = "steelblue") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    plot.title  = element_text(hjust = 0.5)
  ) +
  labs(
    x     = "Collection Month",
    y     = "Number of Sequences",
    title = "Sequence Count by Collection Month (clean dataset, ready for analysis)"
  )

sprintf("Date range: %s to %s",
        format(min(AL_df$Date), "%Y-%m"),
        format(max(AL_df$Date), "%Y-%m"))

sprintf("Total sequences in final dataset: %d", nrow(AL_df))

# Example: partition into non-overlapping 2-month windows
part_data <- partition_time_windows(
  data          = AL_df,
  n_sites       = n_sites,
  window_length = 2,
  window_type   = 3,      # non-overlapping / jumping windows
  verbose       = FALSE
)

sprintf("Number of 2-month partitions created: %d", part_data$N_partitions)
cat("\nPartition date labels:\n")
cat(paste(part_data$Dates_Labels, collapse = "\n"))

# Entropy-based site clustering for the first partition
cat("\nSite clustering for partition 1 (reference period):\n")
part_data$Clusters[[1]]$DataFrame

# Entropy-based site clustering for the last partition
cat("\nSite clustering for partition 10 (reference period):\n")
part_data$Clusters[[10]]$DataFrame
