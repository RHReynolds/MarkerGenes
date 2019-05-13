#' Query genes in specificity matrices.
#'
#' Allows user to query a vector of genes in a celltype_data file with mean
#' expression and specificity matrices.
#'
#' @param genes Query genes as a vector.
#' @param ... Names of input ctds. An individual ctd, should have been generated
#'   using \code{generate.celltype.data} from the EWCE package, and contains
#'   mean expression and specificity matrices from a scRNAseq study.
#' @param celltypeLevel Either 1 or 2, depending on whether user wants to query
#'   level 1 or 2 cell types.
#' @param genelistSpecies Either 'mouse' or 'human' depending on whether MGI or
#'   HGNC symbols are used.
#' @param ctdSpecies Either 'mouse' or 'human' depending on the ctd datasets
#'   being used. Species must be the same across all input ctds.
#'
#' @return Outputs a dataframe, with mean expression and specificity per gene
#'   from each study.
#' @export
#'

query_gene_ctd <- function(genes, ... , celltypeLevel = c(1, 2),
                           genelistSpecies = c("mouse", "human"),
                           ctdSpecies = c("mouse", "human")) {

  # Extract names of ctd inputs to name elements of list
  # Need to remove first unnamed argument, which is the function name, and named arguments.
  argument_names <- as.list(match.call()) %>%
    within(., rm(genes, celltypeLevel, genelistSpecies, ctdSpecies))
  argument_names <- argument_names[-1] %>%
    as.character()

  # Create list
  ctd_list <- setNames(list(...), argument_names)

  # Comment out when exporting package
  source("/home/rreynolds/projects/MarkerGenes/R/conversion_functions.R")

  # Check genelist and ctd species are the same. If not convert genes to same species as ctd.
  if (genelistSpecies == "human" & ctdSpecies == "mouse" |
      genelistSpecies == "mouse" & ctdSpecies == "human") {
    genes <- convert_between_species(genes, genelistSpecies, ctdSpecies)
  }

  # Genes now converted, and can perform query.
  for (i in 1:length(ctd_list)) {
    ctd <- ctd_list[[i]]

    # Filter specificity by gene list
    filtered_specificity <- ctd[[celltypeLevel]]$specificity %>%
      as_tibble(., rownames = "Gene") %>%
      filter(Gene %in% genes) %>%
      gather(key = "CellType", value = "Specificity",-Gene)

    # Filter mean expression by gene list
    filtered_mean_exp <- ctd[[celltypeLevel]]$mean_exp %>%
      as_tibble(., rownames = "Gene") %>%
      filter(Gene %in% genes) %>%
      gather(key = "CellType", value = "Mean_Expression",-Gene)

    # Join both filtered tibbles
    joint_table <- filtered_specificity %>%
      inner_join(filtered_mean_exp, by = c("Gene", "CellType")) %>%
      dplyr::mutate(Study = names(ctd_list[i])) %>%
      dplyr::mutate(Study_species = ctdSpecies) %>%
      dplyr::select(Study, Study_species, everything()) %>%
      arrange(Gene, desc(Specificity))

    # Add tibble to list
    if (i == 1) {
      master_df <- joint_table
    } else {
      master_df <- master_df %>%
        bind_rows(joint_table)
    }

  }

  return(master_df)

}
