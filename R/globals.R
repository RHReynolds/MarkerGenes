# bypass R CMD Check notes, related to tidyverse non-standard evaluation
# https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
utils::globalVariables(c(
  ".",
  "CellType",
  "Gene",
  "MGI_symbol",
  "Species",
  "Specificity",
  "Study",
  "Study_species",
  "data",
  "desc",
  "ensembl_gene_id",
  "hgnc_symbol",
  "mmusculus_homolog_associated_gene_name",
  "mmusculus_homolog_ensembl_gene",
  "mmusculus_homolog_orthology_type",
  "mouse_to_human_orthologs"
))
