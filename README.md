
<!-- README.md is generated from README.Rmd. Please edit that file -->
# MarkerGenes

<!-- badges: start -->
[![Lifecycle: dormant](https://img.shields.io/badge/lifecycle-dormant-blue.svg)](https://www.tidyverse.org/lifecycle/#dormant) \[\[DOI\] <!-- badges: end -->

This repository contains gene cell-type/tissue specificity matrices, which can be used to determine the specificity of a gene to a particular cell type or tissue.

## Installation instructions

To access the datasets within this package, please clone the repository to your local directory.

To use functions within the package, install from github. This can be done using the following lines of code:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

devtools::install_github("RHReynolds/MarkerGenes")
```

Please note that there is no plan to ever submit this code to `CRAN` or `Bioconductor`. This code was developed for personal use.

## Usage

For details, please refer to the vignette.

# License

The code in this repository is released under an MIT license. This repository is distributed in the hope that it will be useful to the wider community, but without any warranty of any kind. Please see the [LICENSE](LICENSE) file for more details.

## Citation

If you use any specificity matrices from this repository, please make sure you cite the original publication the specificity values were derived from. For details, please refer to: [dataset\_metadata.html](docs/dataset_metadata.html)

Below is the citation output from using `citation('MarkerGenes')` in R. Please run this yourself to check for any updates on how to cite **MarkerGenes**.

``` r
print(citation("MarkerGenes"), bibtex = TRUE)
```
