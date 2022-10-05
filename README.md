
<!-- README.md is generated from README.Rmd. Please edit that file -->
# MarkerGenes

<!-- badges: start -->
[![Lifecycle: dormant](https://img.shields.io/badge/lifecycle-dormant-blue.svg)](https://www.tidyverse.org/lifecycle/#dormant) [![DOI](https://zenodo.org/badge/183661929.svg)](https://zenodo.org/badge/latestdoi/183661929) <!-- badges: end -->

This repository contains gene cell-type/tissue specificity matrices, which can be used to determine the specificity of a gene to a particular cell type or tissue.

## Installation instructions

To access the datasets within this package, either clone the repository to your local directory or use the following code chunk (and edit accordingly):

``` r

args <- 
  list(
    url = "https://github.com/RHReynolds/MarkerGenes/raw/master/specificity_matrices_new/ctd_aibsMultipleCrtxSmrtSeq.rda",
    file_name = "ctd_aibsMultipleCrtxSmrtSeq.rda",
    out_dir = here::here("tmp")
  )

# create temporary directory
dir.create(args$out_dir)

# check if file exists and download if it doesn't
if (!file.exists(file.path(args$out_dir, args$file_name))) {

  download.file(
    url = args$url,
    destfile = file.path(args$out_dir, args$file_name)
  )

}
```

To use functions within the package, install from github. This can be done using the following lines of code:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

devtools::install_github("RHReynolds/MarkerGenes")
```

Please note that there is no plan to ever submit this code to `CRAN` or `Bioconductor`. This code was developed for personal use.

## Usage

For details, please refer to the [vignette](https://rhreynolds.github.io/MarkerGenes/articles/MarkerGenes.html).

# License

The code in this repository is released under an MIT license. This repository is distributed in the hope that it will be useful to the wider community, but without any warranty of any kind. Please see the [LICENSE](https://github.com/RHReynolds/MarkerGenes/tree/master/LICENSE.md) for more details.

## Citation

If you use any specificity matrices from this repository, please make sure you cite the original publication the specificity values were derived from. For details, please refer to: [dataset\_metadata.html](https://rhreynolds.github.io/MarkerGenes/articles/articles/dataset_metadata.html).

Below is the citation output from using `citation('MarkerGenes')` in R. Please run this yourself to check for any updates on how to cite **MarkerGenes**.

``` r
utils::citation("MarkerGenes")
#> 
#> Reynolds RH (2022). _MarkerGenes_.
#> https://github.com/RHReynolds/MarkerGenes - R package version 0.99.0,
#> <URL: https://github.com/RHReynolds/MarkerGenes>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {MarkerGenes},
#>     author = {Regina H. Reynolds},
#>     year = {2022},
#>     url = {https://github.com/RHReynolds/MarkerGenes},
#>     note = {https://github.com/RHReynolds/MarkerGenes - R package version 0.99.0},
#>   }
```

# Code contents

Within this repository you will find:

<table>
<colgroup>
<col width="11%" />
<col width="88%" />
</colgroup>
<thead>
<tr class="header">
<th>Directory</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/data">data</a></td>
<td>Package data.</td>
</tr>
<tr class="even">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/data-raw">data-raw</a></td>
<td>Scripts related to generation of data available in repository</td>
</tr>
<tr class="odd">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/docs">docs</a></td>
<td>Source files for repository website</td>
</tr>
<tr class="even">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/flat_lists">flat_lists</a></td>
<td>Marker genes from various sources, as described in folder readme.txt</td>
</tr>
<tr class="odd">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/inst">inst</a></td>
<td>External data used in <a href="https://github.com/RHReynolds/MarkerGenes/tree/master/data-raw">data-raw</a></td>
</tr>
<tr class="even">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/man">man</a></td>
<td>Function documentation</td>
</tr>
<tr class="odd">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/metadata">metadata</a></td>
<td>Metadata related to specificity matrices/dataframes</td>
</tr>
<tr class="even">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/nohup_logs">nohup_logs</a></td>
<td>For any scripts that were run outside of an <code>.Rmd</code> (e.g. scripts from the <a href="https://github.com/RHReynolds/MarkerGenes/tree/master/data-raw">data-raw</a> directory), a log file was recorded and can be accessed here</td>
</tr>
<tr class="odd">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/R">R</a></td>
<td>Repository functions</td>
</tr>
<tr class="even">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/results">results</a></td>
<td>Results relating to .Rmds in <a href="https://github.com/RHReynolds/MarkerGenes/tree/master/workflows">workflows</a></td>
</tr>
<tr class="odd">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/specificity_df">specificity_df</a></td>
<td>Folder of specificity dataframes</td>
</tr>
<tr class="even">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/specificity_matrices">specificity_matrices</a></td>
<td>Folder of specificity matrices compatible with the defunct <a href="https://bioconductor.riken.jp/packages/3.5/bioc/html/EWCE.html">EWCE v1.3.0</a> available on Bioconductor v3.5</td>
</tr>
<tr class="odd">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/specificity_matrices">specificity_matrices_new</a></td>
<td>Folder of specificity matrices compatible with the new <a href="https://nathanskene.github.io/EWCE/index.html">EWCE</a> that is available in Bioconductor&gt;=3.14</td>
</tr>
<tr class="even">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/vignettes">vignettes</a></td>
<td>Repository vignette</td>
</tr>
<tr class="odd">
<td><a href="https://github.com/RHReynolds/MarkerGenes/tree/master/workflows">workflows</a></td>
<td>Miscellaneous .Rmds</td>
</tr>
</tbody>
</table>
