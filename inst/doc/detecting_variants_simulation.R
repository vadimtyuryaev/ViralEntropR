knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  warning  = FALSE,
  message  = FALSE
)

library(ViralEntropR)
library(dplyr)
library(lubridate)
library(ecp)
library(knitr)
library(kableExtra)

ref_sequence <- "MKTIIALSYIFCLVFADYKDDDDK"
sprintf("There are %d sites in the reference sequence", nchar(ref_sequence))

sim_result <- simulate_variant_evolution(
  ref_sequences             = ref_sequence,            # single sequence string
  n_ref_months              = 4,                       # 4-month reference phase
  start_date                = "2020-01-01",
  end_date                  = "2021-12-01",
  variants_config           = c(2, 3, 4, 2),           # 4 variants: 2,3,4,2 mutations
  variant_intervals         = c(4, 5, 5),              # V2 at month 9; V3 at 14; V4 at 19
  n_new_mutations           = 2,                       # 2 de-novo sequences at emergence
  mutation_rate             = c(1.25, 1.25, 1.25, 1.25), # per-variant growth multipliers
  mutation_rate_variability = 0,                       # deterministic growth (no spread)
  deleterious_rate          = 0,                       # no within-variant deleterious mutations
  n_deleterious_limit       = 0,
  n_sequences_total         = 1000,                    # 1000/(24-4) = 50 per period
  ref_variability           = TRUE,                    # variability at sites 1 and 2
  n_ref_sequences           = 200,                     # 200/4 = 50 reference rows per month
  prob_deletion_event       = 0.2,                     # 20% per-month deletion probability
  n_rows_to_delete          = 3,                       # 3 sequences affected if event fires
  seed                      = 2025
)

sim_result$Pool[, c(25, 26, 30)]

print("Variant emergence period:")
sim_result$Variant_Details[[1]]$em

print("Variant 1 sites:")
sim_result$Variant_Details[[1]]$pos %>% sort()

# Month 9 translates into partition 5.
print("Variant 2 emergence period:")
sim_result$Variant_Details[[2]]$em

print("Variant 2 sites:")
sim_result$Variant_Details[[2]]$pos %>% sort()

# Month 14 translates into partition 7.
print("Variant 3 emergence period:")
sim_result$Variant_Details[[3]]$em

print("Variant 3 sites:")
sim_result$Variant_Details[[3]]$pos %>% sort()

# Month 19 translates into partition 10.
print("Variant 4 emergence period:")
sim_result$Variant_Details[[4]]$em

print("Variant 4 sites:")
sim_result$Variant_Details[[4]]$pos %>% sort()

sprintf("The number of deleterious periods is %d",
        length(sim_result$Delet_Records))

for (i in seq_along(sim_result$Delet_Records)) {
  cat(sprintf("Deleterious Period %d: month %d, site %d\n",
              i,
              sim_result$Delet_Records[[i]]$period,
              sim_result$Delet_Records[[i]]$site))
}

sim_result$Simulation_Output[1:10, ]

sim_result$Simulation_Output[sim_result$Simulation_Output$Period == 4, ][1:10, ]

sim_result$Simulation_Output[
  sim_result$Simulation_Output$Period == 5,
  c(sim_result$Variant_Details[[1]]$pos %>% sort(), 25:29)
][1:10, ]

sim_result$Simulation_Output[
  sim_result$Simulation_Output$Period == 6,
  c(sim_result$Variant_Details[[1]]$pos %>% sort(), 25:29)
][1:10, ]

sim_result$Simulation_Output[
  sim_result$Simulation_Output$Period == 8,
  c(sim_result$Variant_Details[[1]]$pos %>% sort(), 25:29)
][1:10, ]

sim_result$Simulation_Output[sim_result$Simulation_Output$Period == 9, ]

sim_result$Simulation_Output[sim_result$Simulation_Output$Period == 13, ]

sim_result$Simulation_Output[sim_result$Simulation_Output$Period == 14, ]

toy_example <- sim_result$Simulation_Output
n_col       <- nchar(ref_sequence)

colnames(toy_example)[seq_len(n_col)] <- as.character(seq_len(n_col))

# Encode character sequences to integers (A=1 … V=20, ambiguous 21–25, unknown=0)
char_mat   <- as.matrix(toy_example[, seq_len(n_col)])
AL_Cat_mix <- as.data.frame(encode_aa_sequence(char_mat))
AL_Cat_mix[] <- lapply(AL_Cat_mix, as.integer)

# Standardise Date column to the first day of each month
AL_Cat_mix$Date <- as.Date(format(as.Date(toy_example$Date), "%Y-%m-01"))

AL_Cat_mix[1:10, ]

part_data <- partition_time_windows(
  data          = AL_Cat_mix,
  n_sites       = n_col,
  window_length = 2,
  window_type   = 3,       # non-overlapping / jumping
  start_date    = "2020-01-01",
  end_date      = "2022-01-01"
)

part_data$Partitions[[1]]$Date %>% unique() %>% sort()

part_data$Partitions[[part_data$N_partitions]]$Date %>% unique() %>% sort()

part_data$Clusters[[1]]$DataFrame

part_data$Clusters[[2]]$DataFrame

part_data$Clusters[[3]]$DataFrame

part_data$Clusters[[4]]$DataFrame

part_data$Clusters[[5]]$DataFrame

part_data$Clusters[[6]]$DataFrame

part_data$Clusters[[7]]$DataFrame

part_data$Clusters[[8]]$DataFrame

part_data$Clusters[[9]]$DataFrame

part_data$Clusters[[10]]$DataFrame

part_data$Clusters[[11]]$DataFrame

part_data$Clusters[[12]]$DataFrame

all_sites <- seq_len(n_col)

hellinger_mat <- calculate_hellinger_matrix(
  partitions = part_data$Partitions,
  sites      = all_sites,
  aa_levels  = 25L   # 25-character alphabet: 20 standard AAs + B, Z, X, *, -
)

dim(hellinger_mat)
hellinger_mat

dat_mat_t <- t(hellinger_mat)
dat_mat_t

cp_m1 <- ks.cp3o(
  Z       = dat_mat_t,
  K       = 11,        # upper bound: at most 11 CPs in 11 partition differences
  minsize = 2,         # each segment must span at least 2 partitions
  verbose = FALSE
)

print("Optimal change point locations (cp_m1$estimates):")
cp_m1$estimates

print("Full solution matrix across k = 1 to K change points (cp_m1$cpLoc):")
cp_m1$cpLoc

lower     <- 1
upper     <- 3
timesteps <- nrow(dat_mat_t) - upper

ECP_res <- detect_changepoints_ecp(
  data_matrix    = dat_mat_t,
  min_window     = lower,
  max_window     = upper,
  n_timesteps    = timesteps,
  rolling_window = FALSE,
  dynamic_k      = TRUE,
  minsize        = 2,
  verbose        = FALSE
)

# Union of all change points detected across all window sizes.
# Indicates which true change points the algorithm recovers
# at any stage of data accumulation (overall sensitivity).
ecp_indices <- ECP_res$ECP_est_list %>%
  unlist() %>%
  .[. > 1] %>%
  unique() %>%
  sort()

print("Change points detected at any window size (union):")
ecp_indices

# Early detection timeline: for each expanding window size, which true
# change points are detectable with the data available up to that point?
# Ground truth: variants emerged at partition indices 3, 5, 7, 10,
# corresponding to Hellinger change point indices 2, 4, 6, 9.
true_cps <- c(2, 4, 6, 9)

detection_timeline <- data.frame(
  window_size  = seq(upper, nrow(dat_mat_t)),
  detected_cps = sapply(ECP_res$ECP_est_list, function(cps) {
    detected <- intersect(true_cps, cps[cps > 1])
    if (length(detected) == 0L) "none" else paste(sort(detected), collapse = ", ")
  })
)

print("Early detection timeline — true change points recovered per window size:")
detection_timeline

cp_m3 <- e.agglo(
  X       = dat_mat_t,
  member  = seq_len(nrow(dat_mat_t)),  # naive: treat each time point as its own segment
  alpha   = 1,
  penalty = function(cps) { 0 }        # no penalty: let the data drive the segmentation
)

print("Change points estimated by method 3 (naive, no penalty):")
cp_m3$estimates

print("Segment membership per time point:")
cp_m3$cluster
