#---Load Libraries------------------------------------------------------------------------------------------------------------------------####
library(Matrix)
library(tidyverse)
library(stringr)
library(EWCE)
library(limma)
library(data.table)

# 2020/03/31 - Files used in this script can be downloaded from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140231
# Script sourced using:
# nohup Rscript /home/rreynolds/projects/MarkerGenes/data-raw/Agarwal2020_GeneratingSImatrix.R &>/home/rreynolds/projects/MarkerGenes/nohup_logs/EWCE_SNIG2020.log&

#---Preparing file for SI calculation-----------------------------------------------------------------------------------------------------####

#---Load in mapping of nuclei to cell type--------------------

cell_ids <- setNames(c(list(read_delim("/home/rreynolds/data/snRNAseq_SNIG/metadata_SNatlas.txt", delim = "\t") %>%
                              dplyr::rename(level_2 = level2class_new, level1 = level1_class)),
                       list(read_delim("/home/rreynolds/data/snRNAseq_SNIG/metadata_CORTEXatlas.txt", delim = "\t") %>%
                              dplyr::mutate(level1 = str_replace_all(level1, " ", "_")))),
                     c("SNIG", "CRTX"))

#---Prepare file dataframe------------------------------------
# Files named such that each set of barcodes, genes and matrix can be identified either by unique "GSM" tag or Sample_*
# Post-mortem individuals (total of 5) identified by C*/N*, with tissue denoted by C (cortex) and N (substantia nigra)
sample_df <- data.frame(file_path = list.files("/home/rreynolds/data/snRNAseq_SNIG/GSE140231_RAW/", full.names = T),
           file_name = list.files("/home/rreynolds/data/snRNAseq_SNIG/GSE140231_RAW/") %>%
             str_replace(., "\\..*", "")
           ) %>%
  tidyr::separate(file_name, into = c("gsm_id", "sample", "sample_id", "tissue_individual", "file_type")) %>%
  dplyr::mutate(tissue = case_when(str_detect(tissue_individual, "N") ~ "SNIG",
                                   str_detect(tissue_individual, "C") ~ "CRTX"),
                individual = str_extract(tissue_individual, "[[:digit:]]")) %>%
  dplyr::select(-sample)

# Let's check there is a total of 5 individuals
print("Check that there are a total of 5 individuals")
c(sample_df %>%
    dplyr::distinct(individual) %>% nrow()) == 5

# Total of 12 samples were sequenced, cortex + substantia nigra from each individual + 2 substantia nigra replicates
# Let's check there are 7 substantia nigra samples and 5 cortical samples
print("Check that there are a total of 7 SNIG and 5 CRTX samples")
sample_df %>%
  dplyr::distinct(tissue_individual, tissue) %>%
  dplyr::group_by(tissue) %>%
  dplyr::summarise(n = n())

#---Check genes match between mtx-----------------------------

# Isolate unique GSM tag
gsm <- sample_df %>% .[["gsm_id"]] %>% unique()

for(i in seq_along(gsm)){

  # Load genes
  genes <- read.table(sample_df %>%
                        dplyr::filter(gsm_id == gsm[i], file_type == "genes") %>%
                        .[["file_path"]] %>%
                        as.character(),
                      header = FALSE,
                      sep = "\t") %>%
    dplyr::select(V2)

  colnames(genes) <- gsm[i]

  if(i == 1){

    genes_df <- genes

  } else{

    genes_df <- genes_df %>%
      dplyr::bind_cols(genes)

  }

}


# Iteratively check that all columns are identical
# Need to -1 from ncol(genes_df), as only need to check that second to last column and last column are identical
print("Check that gene names match across all sample matrices")
sapply(1:c(ncol(genes_df)-1), function(x){

  # This will check that column x is the same as the next column (x + 1)
  identical(genes_df[, x+1],
            genes_df[, x])

}) %>%
  all()

#---Loop and load files---------------------------------------

# Initiate empty lists and vectors for loop
mtx_SNIG <- matrix(nrow = nrow(genes_df))
mtx_CRTX <- matrix(nrow = nrow(genes_df))

# Loop over GSM tags and load in data
for(i in seq_along(gsm)){

  # Load matrix
  mtx <- Matrix::readMM(sample_df %>%
                          dplyr::filter(gsm_id == gsm[i], file_type == "matrix") %>%
                          .[["file_path"]] %>%
                          as.character()) %>%
    as.matrix()

  # Load barcodes
  # Change to format matching barcode mapping
  barcodes <- str_c(sample_df %>%
                      dplyr::filter(gsm_id == gsm[i]) %>%
                      .[["tissue_individual"]] %>% unique(), "_",
                    c(scan(sample_df %>%
                             dplyr::filter(gsm_id == gsm[i], file_type == "barcodes") %>%
                             .[["file_path"]] %>%
                             as.character(),
                           what = '',
                           sep='\n') %>%
                        str_replace("-.*", "")))

  # Load genes
  genes <- read.table(sample_df %>%
                        dplyr::filter(gsm_id == gsm[i], file_type == "genes") %>%
                        .[["file_path"]] %>%
                        as.character(),
                      header = FALSE,
                      sep = "\t") %>%
    dplyr::rename(ens_id = V1, hgnc_symbol = V2)

  # Add rownames and colnames to matrix
  rownames(mtx) <- genes[["hgnc_symbol"]]
  colnames(mtx) <- barcodes

  # Save mtx dimensions
  dim_df <- data.frame(gsm_id = gsm[i],
                       n_genes = dim(mtx)[1],
                       n_barcodes = dim(mtx)[2])

  # Save tissue for conditional statements below
  tissue <- sample_df %>% dplyr::filter(gsm_id == gsm[i]) %>% .[["tissue"]] %>% unique()

  if(i == 1){

    master_dim_df <- dim_df

    if(tissue == "SNIG"){

      mtx_SNIG <- mtx

    } else{

      mtx_CRTX <- mtx

    }

  } else{

    master_dim_df <- master_dim_df %>%
      dplyr::bind_rows(dim_df)

    if(tissue == "SNIG"){

      mtx_SNIG <- mtx_SNIG %>%
        cbind(mtx)

    } else{

      mtx_CRTX <- mtx_CRTX %>%
        cbind(mtx)

    }

  }

}

#---EWCE------------------------------------------------------

# EWCE requires:
# 1. Matrix with data, with rows as genes and columns as cells.
# 2. Annotation dataframe with a minimum of cell_id, level1class and level2class.

# Generate annotation dataframe
cell_ids <- cell_ids %>%
  lapply(., function(x){

    x %>%
      dplyr::select(cell_id, level1class = level1, level2class = level_2)

  })

# Limit matrices to known cell ids from provided cell id dataframes
mtx_CRTX <- mtx_CRTX[,cell_ids$CRTX$cell_id]
mtx_SNIG <- mtx_SNIG[,cell_ids$SNIG$cell_id]

# Create list with metadat and matrix.
# To use with EWCE, must ensure matrix is named 'exp' and metadat is 'annot'.
CRTX <- setNames(list(mtx_CRTX, cell_ids$CRTX),
                 c("exp", "annot"))

SNIG <- setNames(list(mtx_SNIG, cell_ids$SNIG),
                 c("exp", "annot"))

saveRDS(CRTX, file="/home/rreynolds/data/snRNAseq_SNIG/EWCE_data/CRTX_DataForEWCE.Rds")
saveRDS(SNIG, file="/home/rreynolds/data/snRNAseq_SNIG/EWCE_data/SNIG_DataForEWCE.Rds")

#---EWCE: Calculating specificity matrices-----------------------------------------------------------------------------------------------####

Sys.time()
EWCE_data <- setNames(c(list(readRDS("/home/rreynolds/data/snRNAseq_SNIG/EWCE_data/CRTX_DataForEWCE.Rds")),
                        list(readRDS("/home/rreynolds/data/snRNAseq_SNIG/EWCE_data/SNIG_DataForEWCE.Rds"))),
                      c("CRTX", "SNIG"))
print("Data loaded.")

for(i in seq_along(EWCE_data)){

  print(str_c("Working on: ", names(EWCE_data)[i]))

  # Drop genes which do not show significant evidence of varying between level 2 celltypes (based on ANOVA)
  exp_dropped <- drop.uninformative.genes(exp = EWCE_data[[i]]$exp,
                                          level2annot = EWCE_data[[i]]$annot$level2class)

  print("Uninformative genes now removed.")

  annotLevels <- list(level1class = EWCE_data[[i]]$annot$level1class,
                      level2class = EWCE_data[[i]]$annot$level2class)

  print("Annotation levels assigned.")

  # Calculate cell type averages and specificity for each gene
  fNames_allcells <- generate.celltype.data(exp= exp_dropped,
                                            annotLevels = annotLevels,
                                            groupName = str_c("Agarwal2020_", names(EWCE_data)[i]))
  print("Cell type averages and specificity calculated.")


}
