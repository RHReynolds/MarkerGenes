# MarkerGenes
Marker gene resource for cell-type labelling.

1. [Installing the package](#install)
2. [Brief description](#description)
3. [Specificity matrices](#matrices)

## Installing the package <a name="install"></a>
To use, install from github. This can be done using the following lines of code:

``` r
install.packages("devtools")
library(devtools)
install_github("RHReynolds/MarkerGenes")
```

## Brief description <a name="description"></a>
Resource that contains:
- A number of brain-relevant cell-type markers in the form of:
    - [Flat lists of genes.](flat_lists/)
    - [Specificity matrices.](specificity_matrices/)
- For metadata on the available lists/matrices, please refer to:  [MarkerLists.html](workflows/MarkerLists.html).
- Also contains a number of functions that:
    - Allow users to query specificity matrices ([`query_gene_ctd()`](R/query_gene_ctd.R) - for more details, see [query_gene_ctd_tutorial.html](workflows/query_gene_ctd_tutorial.html)) and summarise the output of this ([`summarise_specificity_plot()`](R/summarise_specificity_plot.R)).
    - Run EWCE with multiple gene lists as input and multiple specificity matrices derived from the same organism ([`run_ewce()`](R/run_ewce.R)).

## Specificity matrices <a name="matrices"></a>
- Specificity of a gene represents the proportion of the total expression of a gene found in one cell type as compared to that in all cell types (i.e., the mean expression in one cell type divided by the mean expression in all cell types). If the expression of a gene is shared between two or more cell types, it will get a lower specificity measure.
- For a description of how it's calculated, please refer to [Skene et al. 2016](https://www.frontiersin.org/articles/10.3389/fnins.2016.00016/full)
- **Note**: specificity matrices can be used together with the package [EWCE](https://github.com/NathanSkene/EWCE), which allows you to determine whether a gene list has higher than expected expression in a particular cell type compared to all others in the specificity matrix.

### Structure of a specificity matrix
-	Each .rda file represents one single-cell dataset and typically contains two lists (some datasets may only have level 1 cell types). These two lists represent:
    - Level 1 cell types – broad cell type categories e.g. astrocyte, glutamatergic neuron, etc.
    - Level 2 cell types – subtypes of the broader cell type categories
- With each of the two lists there are three additional lists, which typically have the same name.
     - annot: a character vector that provides the level 1 or level 2 cell type assigned to each single cell, which was used in the calculation
     - mean_exp: a matrix of mean expression for all genes (rows) across each cell type (columns)
     - specificity: a matrix of the specificity values for all genes (rows) across each cell type (columns)
