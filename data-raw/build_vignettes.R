## data-raw/build_vignettes.R
##
## Reproducibly rebuilds the three pre-rendered ViralEntropR vignettes for
## the R.rsp::asis vignette engine.
##
## For each Rmd in data-raw/, the script:
##   1. Renders the Rmd to a temporary HTML output.
##   2. Copies the rendered HTML to vignettes/<name>.html.
##   3. Writes vignettes/<name>.html.asis with a wrapper whose
##      %\VignetteIndexEntry{} is read programmatically from the YAML title.
##
## HOW TO RUN (from the package root):
##   source("data-raw/build_vignettes.R")
##
## Prerequisites:
##   - All packages loaded by the vignettes are installed locally
##     (Biostrings, Rtsne, future.apply, scales, cluster, kableExtra, ...).
##   - The full preprocessing data referenced inside the vignettes is
##     available at the paths the Rmds expect (data-raw/sequences.fasta,
##     etc.). This script does not subset or mock data — it produces the
##     production HTML.
##
## Output:
##   vignettes/<name>.html
##   vignettes/<name>.html.asis
##
## NOT included on the package install path (data-raw/ is in .Rbuildignore).
## ---------------------------------------------------------------------------

stopifnot(requireNamespace("rmarkdown", quietly = TRUE))
stopifnot(requireNamespace("here",      quietly = TRUE))

# Snapshot the user's workspace so the cleanup at the end of this script
# only removes objects this script itself created. Leading dot keeps the
# snapshot hidden from default ls(). Cleanup runs only on successful
# completion: top-level on.exit() does not fire on script exit (it
# attaches to the enclosing function, which does not exist at top level),
# so a render failure mid-loop will leak a few helper objects into
# .GlobalEnv — re-running the script after a fix clears them.
.objs_before_build_vignettes <- ls(envir = .GlobalEnv)

# -- 0. Resolve package root --------------------------------------------------

pkg_root <- tryCatch(here::here(),
                     error = function(e) getwd())

if (!file.exists(file.path(pkg_root, "DESCRIPTION")))
  stop("Cannot locate package root (no DESCRIPTION at: ", pkg_root, ").",
       call. = FALSE)

message("Package root: ", pkg_root)

# -- 1. Vignette inventory ----------------------------------------------------
# basename (no extension) of each vignette. The Rmd source must live at
# data-raw/<name>.Rmd, and outputs are written to vignettes/<name>.html and
# vignettes/<name>.html.asis.

vignette_names <- c(
  "preprocessing_pipeline",
  "detecting_variants_simulation",
  "clustering_accuracy"
)

# -- 2. Helper: read YAML title from an Rmd -----------------------------------
# Fallback: if YAML parsing fails, derive a human-readable title from the
# basename. The wrapper is still valid; only the browseVignettes() display
# string degrades.

extract_yaml_title <- function(rmd_path) {
  yaml <- tryCatch(
    rmarkdown::yaml_front_matter(rmd_path),
    error = function(e) NULL
  )
  if (is.null(yaml) || is.null(yaml$title) || !nzchar(yaml$title)) {
    warning("Could not parse YAML title from ", rmd_path,
            " — using basename fallback.", call. = FALSE)
    return(tools::toTitleCase(gsub("_", " ", tools::file_path_sans_ext(basename(rmd_path)))))
  }
  # Strip embedded backslashes / braces that would break the
  # %\VignetteIndexEntry{} parser. Keep ASCII-printable text only.
  ttl <- gsub("[{}\\\\]", "", yaml$title)
  trimws(ttl)
}

# -- 3. Helper: write one .html.asis wrapper ----------------------------------

write_asis_wrapper <- function(asis_path, title) {
  lines <- c(
    sprintf("%%\\VignetteIndexEntry{%s}", title),
    "%\\VignetteEngine{R.rsp::asis}",
    "%\\VignetteKeyword{HTML}",
    "%\\VignetteKeyword{vignette}"
  )
  writeLines(lines, asis_path)
  message("  wrote ", asis_path)
}

# -- 4. Per-vignette rebuild --------------------------------------------------

vign_dir <- file.path(pkg_root, "vignettes")
src_dir  <- file.path(pkg_root, "data-raw")
dir.create(vign_dir, showWarnings = FALSE, recursive = TRUE)

for (nm in vignette_names) {
  
  rmd_path  <- file.path(src_dir,  paste0(nm, ".Rmd"))
  html_out  <- file.path(vign_dir, paste0(nm, ".html"))
  asis_out  <- file.path(vign_dir, paste0(nm, ".html.asis"))
  
  if (!file.exists(rmd_path)) {
    warning("Skipping '", nm, "': source not found at ", rmd_path,
            call. = FALSE)
    next
  }
  
  message("\n== Rendering: ", nm, " ==")
  tmp_html <- tempfile(fileext = ".html")
  
  rmarkdown::render(
    input         = rmd_path,
    output_file   = tmp_html,
    output_format = "html_document",
    knit_root_dir = pkg_root,
    quiet         = FALSE
  )
  
  # Move rendered HTML into vignettes/, overwriting any previous build.
  ok <- file.copy(tmp_html, html_out, overwrite = TRUE)
  if (!isTRUE(ok))
    stop("Failed to copy rendered HTML to ", html_out, call. = FALSE)
  unlink(tmp_html)
  message("  wrote ", html_out)
  
  # Write/refresh the .html.asis wrapper with a YAML-derived title.
  title <- extract_yaml_title(rmd_path)
  write_asis_wrapper(asis_out, title)
}

message("\nAll vignettes rebuilt.\n")
message("Verify with:")
message("  devtools::install(build_vignettes = TRUE)")
message("  browseVignettes('ViralEntropR')")

# --- 5. Workspace cleanup --------------------------------------------------
# Remove only objects this script added; pre-existing user objects untouched.
# Does not run on render failure; see snapshot comment near the top.
.objs_added <- setdiff(ls(envir = .GlobalEnv), .objs_before_build_vignettes)
rm(list = c(.objs_added, ".objs_added"), envir = .GlobalEnv)
gc()
