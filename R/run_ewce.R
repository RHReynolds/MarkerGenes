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
#'
#' @return A dataframe of EWCE results for each of the gene lists in each study.
#'
#' @example run_ewce_controlled(list_of_genes, reps = 10000, celltypeLevel = 1,
#'   genelistSpecies = "human", sctSpecies = "human", ctd_DRONC_human,
#'   ctd_AIBS2019)
#'
#' @export

run_ewce_controlled <- function(list_of_genes, ..., celltypeLevel = c(1, 2, "both"), reps, genelistSpecies = c("human"),sctSpecies = c("mouse", "human")) {

  # EWCE package necessary to run. Check if installed, and if not, install.
  if(require(EWCE)){
    print("EWCE is loaded correctly")
  } else {
    print("trying to install EWCE")
    library(devtools)
    install_github("nathanskene/ewce")
    if(require(EWCE)){
      print("EWCE installed and loaded")
    } else {
      stop("could not install EWCE")
    }
  }

  # Extract names of ctd inputs to name elements of list
  # Need to remove first unnamed argument, which is the function name, and named arguments.
  argument_names <- as.list(match.call()) %>%
    within(., rm(list_of_genes, celltypeLevel, reps, genelistSpecies, sctSpecies))
  argument_names <- argument_names[-1] %>%
    as.character()

  # Create list
  ctd_list <- setNames(list(...), argument_names)

  # Setting analysis parameters
  # Set seed so if run again, results will be very similar
  set.seed(1234)
  reps = reps # <- For publishable analysis use >10000. In Skene et al. 2016, 100,000 used.

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
        bg = unique(c(hits, m2h$HGNC.symbol))

        # Level 1 cell types: Bootstrap significance testing controlling for transcript length and GC content
        cont_results = bootstrap.enrichment.test(
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
            bind_rows(
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
        bg = unique(c(hits, m2h$HGNC.symbol))

        # Level 1 cell types: Bootstrap significance testing controlling for transcript length and GC content
        cont_results = bootstrap.enrichment.test(
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
            bind_rows(
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
        bg = unique(c(hits, m2h$HGNC.symbol))

        # Level 1 cell types: Bootstrap significance testing controlling for transcript length and GC content
        cont_results = bootstrap.enrichment.test(
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
        cont_results_level2 = bootstrap.enrichment.test(
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
            bind_rows(
              cont_results_level2$results %>%
                dplyr::mutate(
                  GeneSet = list_name,
                  Study = ctd_name %>% str_replace(., "ctd_", "")
                )
            )
        } else{
          Master_df <- Master_df %>%
            bind_rows(
              cont_results$results %>%
                dplyr::mutate(
                  GeneSet = list_name,
                  Study = ctd_name %>% str_replace(., "ctd_", "")
                )
            ) %>%
            bind_rows(
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
