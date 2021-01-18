#' Convert a gene list for \code{query_gene_ctd} to correct species.
#'
#' Allows user to convert between species dependent on species in the query gene
#' vector and the ctd matrices of interest.
#'
#' @param genes Query genes as a vector of MGI or HGNC symbols.
#' @param genelistSpecies Either 'mouse' or 'human' depending on whether MGI or
#'   HGNC symbols are used.
#' @param ctdSpecies Either 'mouse', 'human', or 'both', depending on the ctd
#'   datasets being used. Are they derived purely from mouse, human or both
#'   species?
#'
#' @return A vector with the converted genes.
#' @export
#'

convert_between_species <-function(genes, genelistSpecies = c("mouse", "human"), ctdSpecies = c("mouse", "human")) {

  data("mouse_to_human_orthologs", package = "MarkerGenes", envir = environment())

  ## Filter for only one to one orthologs
  one2one <- mouse_to_human_orthologs %>%
    dplyr::filter(mmusculus_homolog_orthology_type == "ortholog_one2one")

  # If gene lists and ctd are from different species then convert the gene list species to match the species of the ctd
  if (genelistSpecies == "human" & ctdSpecies == "mouse") {
    genes <- one2one %>%
      dplyr::filter(hgnc_symbol %in% genes) %>%
      .[["mmusculus_homolog_associated_gene_name"]]

  }

  if (genelistSpecies == "mouse" & ctdSpecies == "human") {
    genes <- one2one %>%
      dplyr::filter(mmusculus_homolog_associated_gene_name %in% genes) %>%
      .[["hgnc_symbol"]]

  }

  # Check that sufficient genes are still present in the target list
  if (length(genes) < 1) {
    stop("No genes left in the list. Check that correct species set in genelistSpecies and ctdSpecies, and be sure to provide genes as MGI or HGNC symbols.")
  }

  return(genes)

}

#' Convert ensembl IDs for human/mouse to HGNC symbols/MGI symbols,
#' respectively.
#'
#' Allows user to convert ensembl IDs to symbols for human or mouse.
#'
#' @param genes Query genes as a vector of ensembl IDs.
#' @param genelistSpecies Either 'mouse' or 'human'.
#'
#' @return Returns a vector of gene symbols, either HGNC or MGI depending on
#'   selected species.
#' @export
#'

convert_between_ensembl_and_symbols <- function(genes, genelistSpecies = c("mouse", "human")){

  data("mouse_to_human_orthologs", package = "MarkerGenes", envir = environment())

  # # Test code
  # genes <- c("ENSMUSG00000065298", mouse_to_human_orthologs[1:10,] %>% .[["ensembl_gene_id"]])

  if (genelistSpecies == "mouse") {
    # Check that inputted ensembl IDs are from mouse
    if (any(stringr::str_detect(genes, "ENSMUSG") == FALSE)) {
      stop("Somes genes in the inputted vector do not not contain 'ENSMUSG' prefix. Are you sure all inputted gene IDs are from mouse? Or that correct species selected?")
    }

    genes <- mouse_to_human_orthologs %>%
      dplyr::filter(mmusculus_homolog_ensembl_gene %in% genes) %>%
      .[["mmusculus_homolog_associated_gene_name"]]

  }

  if (genelistSpecies == "human") {
    # Check that inputted ensembl IDs are from human
    if (any(stringr::str_detect(genes, "ENSG") == FALSE)) {
      stop("Somes genes in the inputted vector do not not contain 'ENSG' prefix. Are you sure all inputted gene IDs are from human? Or that correct species selected?")
    }

    genes <- mouse_to_human_orthologs %>%
      dplyr::filter(ensembl_gene_id %in% genes) %>%
      .[["hgnc_symbol"]]

  }

  return(genes)

}
