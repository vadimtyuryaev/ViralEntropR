setwd("C:/YORK_PhD/RESEARCH/PAPERS/GitHub/ViralEntropR")

#---------------- simulation --------------------------------------------------#



#---------------- pre-processing ----------------------------------------------#

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

#---------------- clustering accuracy -----------------------------------------#

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