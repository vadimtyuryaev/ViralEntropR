library(devtools)
library(here)

load_all(here())
setwd(here())

#----------------------- variant detection ------------------------------------#

rmarkdown::render(
  "data-raw/detecting_variants_simulation.Rmd",
  output_dir    = "inst/doc",
  knit_root_dir = getwd()
)
knitr::purl("data-raw/detecting_variants_simulation.Rmd",
            output        = "inst/doc/detecting_variants_simulation.R",
            documentation = 0)
file.copy("data-raw/detecting_variants_simulation.Rmd",
          "inst/doc/detecting_variants_simulation.Rmd",
          overwrite = TRUE)

#------------------------ pre-processing --------------------------------------#

rmarkdown::render(
  "data-raw/preprocessing_pipeline.Rmd",
  output_dir    = "inst/doc",
  knit_root_dir = getwd()
)

knitr::purl("data-raw/preprocessing_pipeline.Rmd",
            output = "inst/doc/preprocessing_pipeline.R",
            documentation = 0)

file.copy("data-raw/preprocessing_pipeline.Rmd",
          "inst/doc/preprocessing_pipeline.Rmd",
          overwrite = TRUE)

#----------------------- clustering accuracy ----------------------------------#

rmarkdown::render(
  "data-raw/clustering_accuracy.Rmd",
  output_dir    = "inst/doc",
  knit_root_dir = getwd()
)

knitr::purl("data-raw/clustering_accuracy.Rmd",
            output = "inst/doc/clustering_accuracy.R",
            documentation = 0)

file.copy("data-raw/clustering_accuracy.Rmd",
          "inst/doc/clustering_accuracy.Rmd",
          overwrite = TRUE)

#---------------------- copy Rmd sources to vignettes/ (required for CRAN) ----#

for (f in c("detecting_variants_simulation.Rmd",
            "preprocessing_pipeline.Rmd",
            "clustering_accuracy.Rmd")) {
  file.copy(
    file.path("data-raw", f),
    file.path("vignettes", f),
    overwrite = TRUE
  )
}

#---------------------- rebuild package documentation -------------------------#

devtools::document()

#---------------------- clean up session --------------------------------------#
rm(list = ls())