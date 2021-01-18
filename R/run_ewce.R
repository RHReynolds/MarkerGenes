#' Run EWCE with control for GC content and transcript length.
#'
#' Perform EWCE bootstrapping across several gene lists in more than one
#' specificity matrix (aka. ctd).
#'
#' @param list_of_genes Named list containing vectors with gene IDs
#'   corresponding to each gene list to be tested.
#' @param ... Names of input ctds. An individual ctd, should have been generated
#'   using \code{generate.celltype.data} from the EWCE package, and contains
#'   mean expression and specificity matrices from a scRNAseq study.
#' @param celltypeLevel Choose whether to run only with level 1 cell types,
#'   level 2 cell types or both. 1 for only level 1, or "both" for both.
#' @param reps Number of repeats that should be performed in bootstrapping.
#' @param genelistSpecies Has to be human, as geneSizeControl is set to TRUE,
#'   and this requires extracting GC content and transcript length from the
#'   human ENSEMBL.
#' @param sctSpecies Either 'mouse' or 'human' depending on the ctd datasets
#'   being used. Species must be the same across all input ctds.
#' @param mouse_to_human Default is to use the default dataframe supplied by
#'   within MarkerGenes packages, which contains all Human->Mouse orthologs for
#'   all human genes. If user wishes to supply their own dataframe, must contain
#'   the columns, 'hgnc_symbol' and 'MGI_symbol'.
#'
#' @return A dataframe of EWCE results for each of the gene lists in each study.
#' @export

run_ewce_controlled <- function(list_of_genes, ..., celltypeLevel = c(1, 2, "both"),
                                reps, genelistSpecies = c("human"),sctSpecies = c("mouse", "human"),
                                mouse_to_human = NULL) {

  # Extract names of ctd inputs to name elements of list
  # Need to remove first unnamed argument, which is the function name, and named arguments.
  argument_names <- as.list(match.call()) %>%
    within(., rm(list_of_genes, celltypeLevel, reps, genelistSpecies, sctSpecies, mouse_to_human))
  argument_names <- argument_names[-1] %>%
    as.character()

  # Create list
  ctd_list <- setNames(list(...), argument_names)

  # Setting analysis parameters
  # Set seed so if run again, results will be very similar
  set.seed(1234)
  reps = reps # <- For publishable analysis use >10000. In Skene et al. 2016, 100,000 used.

  # Loading mouse to human orthologs
  if(is.null(mouse_to_human)){
    data("mouse_to_human_orthologs", package = "MarkerGenes", envir = environment())
    ## Filter for only one to one orthologs
    m2h <- mouse_to_human_orthologs %>%
      dplyr::filter(mmusculus_homolog_orthology_type == "ortholog_one2one") %>%
      dplyr::select(hgnc_symbol, mmusculus_homolog_associated_gene_name) %>%
      dplyr::rename(MGI_symbol = mmusculus_homolog_associated_gene_name) %>%
      dplyr::filter(hgnc_symbol != "") %>%
      dplyr::distinct(hgnc_symbol, MGI_symbol)
  } else{
    m2h <- mouse_to_human
  }

  # Currently a loop, but need to parallelise at some point
  if(celltypeLevel == 1){

    for (i in seq_along(list_of_genes)) {
      list_name <- names(list_of_genes)[i]

      for (j in seq_along(ctd_list)) {
        ctd_name <- names(ctd_list)[j]

        print(str_c("Performing EWCE analysis for: ", list_name, " in ", ctd_name))

        hits = list_of_genes[[i]] %>% unique()

        # The next step is to determine the most suitable background set. The experimental methods used to find these gene are all genome wide,
        # so there is no restriction imposed as a result of that. Thus our initial background set is the set of all mouse/human genes.
        # Not all human genes have mouse orthologs however, so we need to drop all genes from the target and background set
        # which do not have mouse orthologs.
        bg = unique(c(hits, m2h$hgnc_symbol))

        # Level 1 cell types: Bootstrap significance testing controlling for transcript length and GC content
        cont_results =
          EWCE::bootstrap.enrichment.test(
          sct_data = ctd_list[[j]],
          hits = hits,
          bg = bg,
          reps = reps,
          annotLevel = 1,
          geneSizeControl = TRUE,
          genelistSpecies = genelistSpecies,
          sctSpecies = sctSpecies)

        # Assembling results dataframe
        if (i == 1 & j == 1) {
          Master_df <- cont_results$results %>%
            dplyr::mutate(GeneSet = list_name,
                          Study = ctd_name %>% str_replace(., "ctd_", ""))
        } else{
          Master_df <- Master_df %>%
            dplyr::bind_rows(
              cont_results$results %>%
                dplyr::mutate(
                  GeneSet = list_name,
                  Study = ctd_name %>% str_replace(., "ctd_", "")
                )
            )
        }

      }

    }

  }

  if(celltypeLevel == 2){

    for (i in seq_along(list_of_genes)) {
      list_name <- names(list_of_genes)[i]

      for (j in seq_along(ctd_list)) {
        ctd_name <- names(ctd_list)[j]

        print(str_c("Performing EWCE analysis for: ", list_name, " in ", ctd_name))

        hits = list_of_genes[[i]] %>% unique()

        # Background set
        bg = unique(c(hits, m2h$hgnc_symbol))

        # Level 1 cell types: Bootstrap significance testing controlling for transcript length and GC content
        cont_results =
          EWCE::bootstrap.enrichment.test(
          sct_data = ctd_list[[j]],
          hits = hits,
          bg = bg,
          reps = reps,
          annotLevel = 2,
          geneSizeControl = TRUE,
          genelistSpecies = genelistSpecies,
          sctSpecies = sctSpecies)

        # Assembling results dataframe
        if (i == 1 & j == 1) {
          Master_df <- cont_results$results %>%
            dplyr::mutate(GeneSet = list_name,
                          Study = ctd_name %>% str_replace(., "ctd_", ""))
        } else{
          Master_df <- Master_df %>%
            dplyr::bind_rows(
              cont_results$results %>%
                dplyr::mutate(
                  GeneSet = list_name,
                  Study = ctd_name %>% str_replace(., "ctd_", "")
                )
            )
        }

      }

    }

  }

  if(celltypeLevel == "both"){

    for (i in seq_along(list_of_genes)) {
      list_name <- names(list_of_genes)[i]

      for (j in seq_along(ctd_list)) {
        ctd_name <- names(ctd_list)[j]

        print(str_c("Performing EWCE analysis for: ", list_name, " in ", ctd_name))

        hits = list_of_genes[[i]] %>% unique()

        # Background set
        bg = unique(c(hits, m2h$hgnc_symbol))

        # Level 1 cell types: Bootstrap significance testing controlling for transcript length and GC content
        cont_results =
          EWCE::bootstrap.enrichment.test(
          sct_data = ctd_list[[j]],
          hits = hits,
          bg = bg,
          reps = reps,
          annotLevel = 1,
          geneSizeControl = TRUE,
          genelistSpecies = genelistSpecies,
          sctSpecies = sctSpecies
        )

        # Level 2 cell types: Bootstrap significance testing controlling for transcript length and GC content
        cont_results_level2 =
          EWCE::bootstrap.enrichment.test(
          sct_data = ctd_list[[j]],
          hits = hits,
          bg = bg,
          reps = reps,
          annotLevel = 2,
          geneSizeControl = TRUE,
          genelistSpecies = genelistSpecies,
          sctSpecies = sctSpecies
        )

        # Assembling results dataframe
        if (i == 1 & j == 1) {
          Master_df <- cont_results$results %>%
            dplyr::mutate(GeneSet = list_name,
                          Study = ctd_name %>% str_replace(., "ctd_", "")) %>%
            dplyr::bind_rows(
              cont_results_level2$results %>%
                dplyr::mutate(
                  GeneSet = list_name,
                  Study = ctd_name %>% str_replace(., "ctd_", "")
                )
            )
        } else{
          Master_df <- Master_df %>%
            dplyr::bind_rows(
              cont_results$results %>%
                dplyr::mutate(
                  GeneSet = list_name,
                  Study = ctd_name %>% str_replace(., "ctd_", "")
                )
            ) %>%
            dplyr::bind_rows(
              cont_results_level2$results %>%
                dplyr::mutate(
                  GeneSet = list_name,
                  Study = ctd_name %>% str_replace(., "ctd_", "")
                )
            )

        }

      }

    }

  }

  return(Master_df)

}
