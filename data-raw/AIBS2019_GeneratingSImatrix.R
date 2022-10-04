# Description: create AIBS ctd

# Dataset includes:
# Single-nucleus transcriptomes from 49,495 nuclei across multiple human cortical areas.
# Individual layers of cortex were dissected from tissues covering the middle temporal gyrus (MTG), anterior cingulate cortex
# (ACC; also known as the ventral division of medial prefrontal cortex, A24), primary visual cortex (V1C), primary motor cortex
# (M1C), primary somatosensory cortex (S1C) and primary auditory cortex (A1C) derived from human brain.
# Nuclei were dissociated and sorted using the neuronal marker NeuN. Nuclei were sampled from post-mortem and
# neurosurgical (MTG only) donor brains and expression was profiled with SMART-Seq v4 RNA-sequencing.

# Data was downloaded from: https://portal.brain-map.org/atlases-and-data/rnaseq/human-multiple-cortical-areas-smart-seq
# Used exon + intron counts

# Load packages -----------------------------------------------------------

library(dplyr)
library(EWCE)
library(ewceData)
library(forcats)
library(Matrix)
library(readxl)
library(stringr)
library(vroom)

# Set arguments -----------------------------------------------------------

args <-
  list(
    data_dir = file.path("/home/rreynolds/data/scRNAseq_AIBS/multiple_cortical_smartseq")
  )

args <-
  c(
    args,
    list(
      path_to_mtx = file.path(args$data_dir, "matrix.csv"),
      path_to_metadata = file.path(args$data_dir, "metadata.csv"),
      cores = 10
    )
  )

# Load data ---------------------------------------------------------------

print(stringr::str_c(Sys.time(), " - starting script..."))

mtx <-
  vroom::vroom(args$path_to_mtx)

metadata <-
  readr::read_csv(args$path_to_metadata)

# Main --------------------------------------------------------------------

sample_name <- mtx$sample_name
genes <- colnames(mtx[!names(mtx) %in% c("sample_name")])

print(stringr::str_c(Sys.time(), " - converting csv to sparse matrix..."))

# Sparse matrix
mtx <-
  mtx[!names(mtx) %in% c("sample_name")] %>%
  as.matrix() %>%
  Matrix::Matrix(sparse = T)

rownames(mtx) <- sample_name

# Remove outlier calls (outlier based on metadata)
mtx <-
  mtx[!rownames(mtx) %in%
        c(
          metadata %>%
            dplyr::filter(outlier_call == TRUE) %>%
            dplyr::pull(sample_name)
        ),]

# Create annotation labels
# Small number of endothelial (70)/pericyte (32)/VLMC(11)
# Therefore, chose to combine these 3 to "vascular cells" for level 1 annotation
annotation <-
  metadata %>%
  dplyr::filter(outlier_call == FALSE) %>%
  dplyr::mutate(
    sample_name =
      forcats::fct_relevel(
        sample_name,
        rownames(mtx)
      ),
    level1class =
      dplyr::case_when(
        class_label == "Non-neuronal" &
          subclass_label %in% c("Oligodendrocyte", "OPC", "Astrocyte", "Microglia") ~ subclass_label,
        class_label == "Non-neuronal" &
          subclass_label %in% c("Endothelial", "Pericyte", "VLMC") ~ "Vascular_cells",
        TRUE ~ stringr::str_c(
          class_label,
          "_",
          EWCE::fix_celltype_names(subclass_label)
        )
      ),
    level2class =
      cluster_label %>%
      EWCE::fix_celltype_names()
  ) %>%
  dplyr::select(
    sample_name, level1class, level2class
  ) %>%
  dplyr::arrange(sample_name)

# Double-check row order is equivalent
print("Do rownames of matrix have same order as metadata sample names?")
all(annotation$sample_name == rownames(mtx))

# Fix bad gene symbols and drop uninformative genes
exp <-
  mtx %>%
  t() %>%
  EWCE::fix_bad_hgnc_symbols(dropNonHGNC = TRUE) %>%
  EWCE::drop_uninformative_genes(
    level2annot = annotation$level2class,
    mtc_method = "BH",
    adj_pval_thresh = 1e-05,
    input_species = "human",
    output_species = "human",
    no_cores = args$cores
  )

rm(mtx)

EWCE::generate_celltype_data(
  exp = exp,
  annotLevels =
    list(
      level1class = annotation$level1class,
      level2class = annotation$level2class
    ),
  groupName = "aibsMultipleCrtxSmrtSeq",
  savePath = args$data_dir,
  no_cores = args$cores
)

print(stringr::str_c(Sys.time(), " - script done!"))
