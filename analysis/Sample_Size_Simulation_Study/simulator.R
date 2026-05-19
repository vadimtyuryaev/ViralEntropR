# =============================================================================
# simulator.R
# =============================================================================
#
# Purpose-built simulator for the Sample-Size Simulation Study.
#
# Produces an integer-encoded amino-acid matrix representing one snapshot of
# a viral sequence population: reference rows plus zero or more variant
# blocks. Adds optional surveillance-level deleterious-mutation noise after
# variant substitution. No temporal axis, no growth dynamics, no
# change-point logic; the snapshot is consumed by detect_in_snapshot.R.
#
# All functions are vectorised at the variant level. Substitute amino-acid
# draws use a precomputed 20x19 lookup matrix (.SUBSTITUTE_MATRIX in setup.R),
# avoiding per-position setdiff() calls.
#
# Author : Vadim Tyuryaev
# =============================================================================

# -----------------------------------------------------------------------------
# simulate_population_snapshot
# -----------------------------------------------------------------------------
# Builds the population matrix.
#
# Inputs
#   ref_seq_int    : integer vector of length L. The encoded reference
#                    sequence (output of encode_aa_sequence on the canonical
#                    Spike protein).
#   n_ref          : integer scalar. Number of identical reference rows.
#   variants_spec  : list of variant specifications. Each element is a list
#                    with three components:
#                      $positions   : integer vector of mutated positions
#                      $substitutes : integer vector of substitute codes,
#                                     same length as $positions
#                      $n_sequences : integer scalar, number of identical
#                                     rows of this variant
#                    Variant order in the list determines row order in the
#                    output: reference rows first, then variants_spec[[1]],
#                    then variants_spec[[2]], etc.
#   seed           : integer or NULL. Optional reseed before any internal
#                    randomness. The simulator itself is deterministic given
#                    its inputs; the seed is forwarded for caller convenience.
#
# Returns a list with
#   $matrix     : integer matrix (total_rows x L) with values in 1..20
#                 (only standard amino acids; the simulator never produces
#                 ambiguous codes; deleterious noise is applied separately
#                 via apply_deleterious_noise()).
#   $row_labels : integer vector of length total_rows. Each row's source:
#                 0 = reference, k = index of variant in variants_spec.
#
# Implementation notes
#   - Pre-allocates the full matrix in one C-level call via
#     matrix(ref_seq_int, nrow = total_rows, ncol = L, byrow = TRUE).
#   - Variant substitutions overwrite full column slices in vectorised
#     assignments mat[rows, position] <- substitute.
#   - No per-row loops.
# -----------------------------------------------------------------------------
simulate_population_snapshot <- function(ref_seq_int,
                                         n_ref,
                                         variants_spec,
                                         seed = NULL) {

  # ----- Input validation --------------------------------------------------
  if (!is.integer(ref_seq_int))
    stop("`ref_seq_int` must be an integer vector.", call. = FALSE)
  if (length(ref_seq_int) < 1L)
    stop("`ref_seq_int` must have positive length.", call. = FALSE)
  n_ref <- as.integer(n_ref)
  if (length(n_ref) != 1L || is.na(n_ref) || n_ref < 0L)
    stop("`n_ref` must be a single non-negative integer.", call. = FALSE)
  if (!is.null(variants_spec) && !is.list(variants_spec))
    stop("`variants_spec` must be NULL or a list of variant specifications.",
         call. = FALSE)

  L <- length(ref_seq_int)
  n_variants <- length(variants_spec)

  # ----- Per-variant row counts --------------------------------------------
  n_per_var <- if (n_variants > 0L) {
    vapply(variants_spec, function(v) as.integer(v$n_sequences), integer(1L))
  } else {
    integer(0L)
  }
  total_rows <- n_ref + sum(n_per_var)

  if (total_rows == 0L) {
    return(list(
      matrix     = matrix(integer(0L), nrow = 0L, ncol = L),
      row_labels = integer(0L)
    ))
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  # ----- Pre-fill matrix with the reference sequence on every row ----------
  # matrix(byrow = TRUE) is fully vectorised at the C level: one allocation,
  # one fill pass. Storage mode is integer because ref_seq_int is integer.
  mat <- matrix(ref_seq_int, nrow = total_rows, ncol = L, byrow = TRUE)
  storage.mode(mat) <- "integer"

  row_labels <- integer(total_rows)
  if (n_ref > 0L) row_labels[seq_len(n_ref)] <- 0L

  # ----- Apply per-variant substitutions -----------------------------------
  row_offset <- n_ref
  for (i in seq_len(n_variants)) {
    v <- variants_spec[[i]]
    n_v <- n_per_var[i]
    if (n_v == 0L) next

    positions   <- as.integer(v$positions)
    substitutes <- as.integer(v$substitutes)
    if (length(positions) != length(substitutes))
      stop(sprintf(
        "Variant %d: length(positions) (%d) != length(substitutes) (%d).",
        i, length(positions), length(substitutes)), call. = FALSE)

    rows <- seq.int(row_offset + 1L, row_offset + n_v)
    row_labels[rows] <- i

    # Vectorised per-position assignment. The inner loop is over the number
    # of mutations (typically 2..15, max 33 for Omicron), not over rows.
    for (k in seq_along(positions)) {
      mat[rows, positions[k]] <- substitutes[k]
    }

    row_offset <- row_offset + n_v
  }

  list(matrix = mat, row_labels = row_labels)
}

# -----------------------------------------------------------------------------
# apply_deleterious_noise
# -----------------------------------------------------------------------------
# Adds independent low-frequency substitution noise to an already-built
# population matrix. Models the surveillance-level prevalence of transient
# deleterious mutations that did not fix in any circulating lineage.
#
# Inputs
#   mat   : integer matrix of size (n_rows x n_cols), values in 1..20 at
#           every cell (output of simulate_population_snapshot before any
#           noise has been applied). Other code values are permitted but
#           are passed through unchanged (the noise mechanism only acts
#           on cells whose current residue is in 1..20).
#   p_del : numeric scalar in [0, 1). Per-(sequence, position) Bernoulli
#           probability that a cell is overwritten with a uniformly chosen
#           non-self residue from 1..20. Default 1e-3, anchored empirically
#           by the 0.1%-frequency threshold analysis on 21 SARS-CoV-2
#           variant-defining Spike positions (pooled rate 9.6e-4).
#   seed  : integer or NULL. Optional reseed before any internal randomness.
#
# Returns: integer matrix of the same shape, with noise applied. Storage
# mode is integer.
#
# Implementation notes
#   - The expected number of noise events per matrix is
#     n_rows * n_cols * p_del. For a small-band cell with ~40 rows and
#     1273 columns, p_del = 1e-3 gives ~0.05 events per matrix. For a
#     large-band cell with ~50,000 rows, p_del = 1e-3 gives ~63,650 events
#     per matrix. The implementation scales linearly with the number of
#     events, not with matrix size.
#   - The per-cell Bernoulli is realised in one rbinom() call drawing
#     n_rows * n_cols values, then converted to row/column indices via
#     arithmetic on the affected linear positions. There is no per-cell
#     loop in R.
#   - Substitutes are looked up via .SUBSTITUTE_MATRIX. A cell currently
#     holding residue r is replaced with one of the 19 non-r residues
#     uniformly. Noise can therefore (rarely) flip a variant's defining
#     substitution to a different residue or back to the reference; this
#     is the natural reading of "deleterious mutations may occur at any
#     position including variant-defining positions" and is biologically
#     realistic at p_del = 1e-3.
#   - When p_del = 0 the function is a no-op and returns mat unchanged.
# -----------------------------------------------------------------------------
apply_deleterious_noise <- function(mat, p_del = 1e-3, seed = NULL) {

  if (!is.matrix(mat) || !is.integer(mat))
    stop("`mat` must be an integer matrix.", call. = FALSE)
  p_del <- as.numeric(p_del)
  if (length(p_del) != 1L || is.na(p_del) || p_del < 0 || p_del >= 1)
    stop("`p_del` must be a single numeric value in [0, 1).", call. = FALSE)

  if (p_del == 0) return(mat)

  if (!is.null(seed)) set.seed(as.integer(seed))

  n_cells <- length(mat)
  if (n_cells == 0L) return(mat)

  # Vectorised Bernoulli over every cell.
  noise_flags <- rbinom(n_cells, size = 1L, prob = p_del) == 1L
  if (!any(noise_flags)) return(mat)

  affected_idx <- which(noise_flags)
  current_resi <- mat[affected_idx]

  # Restrict to cells currently in 1..20 (the standard alphabet). Cells
  # outside that range are passed through unchanged. Under the simulator's
  # invariants this restriction is a no-op (the matrix is pure 1..20), but
  # it is defensive armour in case noise is ever applied to a matrix that
  # contains ambiguous codes.
  in_range <- current_resi >= 1L & current_resi <= 20L
  if (!any(in_range)) return(mat)

  affected_idx <- affected_idx[in_range]
  current_resi <- current_resi[in_range]

  # Draw a uniform substitute from .SUBSTITUTE_MATRIX[current_resi, ].
  col_idx        <- sample.int(19L, length(affected_idx), replace = TRUE)
  substitute_resi <- .SUBSTITUTE_MATRIX[cbind(current_resi, col_idx)]

  mat[affected_idx] <- substitute_resi
  storage.mode(mat) <- "integer"

  mat
}

# -----------------------------------------------------------------------------
# draw_variant_spec
# -----------------------------------------------------------------------------
# Builds one variant specification (positions, substitutes, n_sequences) by
# drawing biological parameters per the study's empirical rules.
#
# Inputs
#   n_mutations  : integer scalar. Total number of mutations the variant
#                  carries, including any forced D614G baseline.
#   force_d614g  : logical scalar. If TRUE, the variant's first mutation is
#                  fixed to position 614 with substitute G; the remaining
#                  n_mutations - 1 positions are drawn from draw_pool.
#                  If FALSE, all n_mutations positions are drawn from
#                  draw_pool (currently unused; reserved for future
#                  extensions).
#   draw_pool    : integer vector of candidate positions. Should be
#                  POOL_11_DRAW (Sc1 V1, V2; Sc2 V1, V2, V3; Sc3 V1, V2, V3,
#                  emerging) or POOL_12_DRAW (Sc4 emerging).
#   n_sequences  : integer scalar. Number of identical rows of this variant
#                  to include in the population (passed straight through).
#   ref_seq_int  : integer vector. Encoded reference sequence; used to
#                  ensure substitutes differ from the reference residue at
#                  each chosen position.
#
# Returns a list with $positions, $substitutes, $n_sequences ready for
# simulate_population_snapshot().
#
# Implementation notes
#   - Substitutes drawn via a single matrix lookup against the precomputed
#     .SUBSTITUTE_MATRIX (20x19) defined in setup.R.
#   - No setdiff() inside any loop; the substitute pool is precomputed.
# -----------------------------------------------------------------------------
draw_variant_spec <- function(n_mutations,
                              force_d614g,
                              draw_pool,
                              n_sequences,
                              ref_seq_int) {

  n_mutations <- as.integer(n_mutations)
  n_sequences <- as.integer(n_sequences)
  draw_pool   <- as.integer(draw_pool)

  if (length(n_mutations) != 1L || is.na(n_mutations) || n_mutations < 1L)
    stop("`n_mutations` must be a single positive integer.", call. = FALSE)
  if (!is.logical(force_d614g) || length(force_d614g) != 1L)
    stop("`force_d614g` must be a single logical.", call. = FALSE)
  if (length(draw_pool) < 1L)
    stop("`draw_pool` must be non-empty.", call. = FALSE)
  if (length(n_sequences) != 1L || is.na(n_sequences) || n_sequences < 0L)
    stop("`n_sequences` must be a single non-negative integer.", call. = FALSE)

  if (force_d614g) {
    positions_baseline   <- 614L
    substitutes_baseline <- .AA_TO_INT[["G"]]
    n_remaining          <- n_mutations - 1L
  } else {
    positions_baseline   <- integer(0L)
    substitutes_baseline <- integer(0L)
    n_remaining          <- n_mutations
  }

  if (n_remaining > length(draw_pool))
    stop(sprintf(
      "Cannot draw %d positions from a pool of size %d.",
      n_remaining, length(draw_pool)), call. = FALSE)

  if (n_remaining > 0L) {
    # Sample positions without replacement from the pool
    pool_idx          <- sample.int(length(draw_pool), n_remaining, replace = FALSE)
    positions_random  <- sort(draw_pool[pool_idx])

    # Substitute lookup: at each position with reference residue r in 1..20,
    # pick one of the 19 non-r residues. Vectorised matrix indexing.
    ref_at_pos        <- ref_seq_int[positions_random]
    col_idx           <- sample.int(19L, n_remaining, replace = TRUE)
    substitutes_random <- .SUBSTITUTE_MATRIX[cbind(ref_at_pos, col_idx)]
  } else {
    positions_random   <- integer(0L)
    substitutes_random <- integer(0L)
  }

  list(
    positions   = as.integer(c(positions_baseline, positions_random)),
    substitutes = as.integer(c(substitutes_baseline, substitutes_random)),
    n_sequences = n_sequences
  )
}

# -----------------------------------------------------------------------------
# load_reference_sequence
# -----------------------------------------------------------------------------
# Reads a FASTA file containing the canonical Spike protein, encodes it to
# the package's 25-symbol integer alphabet, and returns the integer vector.
#
# Inputs
#   fasta_path : path to a FASTA file containing one or more sequences.
#                The first sequence is used.
#   expected_L : integer. Expected sequence length. Default 1273.
#
# Returns: integer vector of length expected_L with codes in 1..25.
#
# Errors
#   - Missing file.
#   - Missing Biostrings package.
#   - Empty FASTA.
#   - Length mismatch.
#   - Any non-standard residue (code 0 or 21..25) in the reference, which
#     would invalidate the simulator's assumption that ref_seq_int[p] is in
#     1..20 for substitute lookup.
# -----------------------------------------------------------------------------
load_reference_sequence <- function(fasta_path, expected_L = 1273L) {
  if (!file.exists(fasta_path))
    stop("Reference FASTA not found: ", fasta_path, call. = FALSE)

  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("Package 'Biostrings' is required to read the reference FASTA. ",
         "Install via BiocManager::install('Biostrings').", call. = FALSE)

  fasta <- Biostrings::readAAStringSet(fasta_path)
  if (length(fasta) == 0L)
    stop("FASTA file contains no sequences: ", fasta_path, call. = FALSE)

  ref_chars <- strsplit(as.character(fasta[[1L]]), "", fixed = TRUE)[[1L]]
  if (length(ref_chars) != expected_L)
    stop(sprintf("Reference sequence length %d does not match expected %d.",
                 length(ref_chars), expected_L), call. = FALSE)

  ref_int_mat <- ViralEntropR::encode_aa_sequence(matrix(ref_chars, nrow = 1L))
  ref_int     <- as.integer(ref_int_mat[1L, ])

  bad <- which(ref_int < 1L | ref_int > 20L)
  if (length(bad) > 0L)
    stop(sprintf(
      "Reference sequence contains non-standard residues at %d position(s) ",
      length(bad)),
      "(codes outside 1..20). The simulator requires a fully standard ",
      "reference sequence; positions are: ",
      paste(utils::head(bad, 20L), collapse = ", "),
      if (length(bad) > 20L) sprintf(" ... and %d more.", length(bad) - 20L) else "",
      call. = FALSE)

  ref_int
}
