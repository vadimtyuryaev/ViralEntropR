## Test environments

* local: Windows 11, R 4.x.x
* GitHub Actions:
  - ubuntu-latest (R-devel, R-release, R-oldrel-1)
- macos-latest (R-release)
- windows-latest (R-release)
* win-builder: R-devel, R-release
* R-hub v2: linux, macos-arm64, windows, donttest, nosuggests

## R CMD check results

0 errors | 0 warnings | 0 notes

## Submission notes

First submission of ViralEntropR.

Biostrings (Bioconductor) is declared in Suggests because it is referenced
only in @examples blocks. The package's exported functions operate on
AAStringSet objects but do not themselves call into Bioconductor packages,
so installation of ViralEntropR does not require BiocManager.

The bundled vignettes are pre-rendered via R.rsp::asis because they
reference large external datasets archived on Zenodo (DOI:
10.5281/zenodo.19040165) and on GISAID, which cannot be redistributed
inside the tarball.

## Downstream dependencies

There are currently no downstream dependencies for this package.