#---Load Libraries------------------------------------------------------------------------------------------------------------------------####
library(tidyverse)
library(EWCE)
library(limma)
library(data.table)

# 2019/01/07 - Files used in this script can be downloaded from https://portal.brain-map.org/atlases-and-data/rnaseq as part of their October 2018 release
# Note their MTG release is dated slightly differently i.e. June 2018, as opposed to the other two brain regions which are dated October 2018

#---Preparing file for SI calculation-----------------------------------------------------------------------------------------------------####
# gene_rows <- read_delim(file = "/home/rreynolds//human_MTG_2018-06-14_genes-rows.txt", delim = ",")
# # sample_columns <- read_csv(file = "/home/rreynolds/data/scRNAseq_AIBS/MTG/human_MTG_2018-06-14_samples-columns.csv")
# sample_columns <- read_csv(file = "/home/rreynolds/data/scRNAseq_AIBS/MTG/human_MTG_2018-06-14_samples-columns_level1added.csv")
#
# # Only use exons, on account of introns likely reflecting levels of pre-mRNA.
# exons <- fread(file = "/home/rreynolds/data/scRNAseq_AIBS/MTG/human_MTG_2018-06-14_exon-matrix.csv")
# # introns <- fread(file = "/home/rreynolds/data/scRNAseq_AIBS/MTG/human_MTG_2018-06-14_intron-matrix.csv")
#
# # EWCE requires:
# # 1. Matrix with data, with rows as genes and columns as cells.
# # 2. Annotation dataframe with a minimum of cell_id, level1class and level2class.
#
# # # Currently 4 classes: Glutamatergic, GABAergic, non-neuronal, no class
# # # Remove no class from metadata.
# # # Split non-neuronal into common cell types e.g. oligos, astrocytes, etc.
# # sample_columns <- sample_columns %>%
# #   unite(level1, c("class", "cluster"), sep = ":", remove = FALSE) %>%
# #   mutate(level1 = str_replace(level1, "Glutamatergic:.*", "Glutamatergic"),
# #          level1 = str_replace(level1, "GABAergic:.*", "GABAergic"),
# #          level1 = str_replace(level1, "Non-neuronal:Oligo.*", "Oligodendrocyte"),
# #          level1 = str_replace(level1, "Non-neuronal:OPC.*", "OPC"),
# #          level1 = str_replace(level1, "Non-neuronal:Astro.*", "Astrocyte"),
# #          level1 = str_replace(level1, "Non-neuronal:Micro.*", "Microglia"),
# #          level1 = str_replace(level1, "Non-neuronal:Endo.*", "Endothelial cell"))
# # write_csv(sample_columns, path = "/home/rreynolds/data/scRNAseq_AIBS/MTG/human_MTG_2018-06-14_samples-columns_level1added.csv")
#
# # Generate annotation dataframe
# metadat <- sample_columns %>%
#   dplyr::filter(!class == "no class") %>%
#   dplyr::select(sample_name, class, level1, cluster) %>%
#   mutate(cell_id = sample_name,
#          level1class = level1,
#          level2class = str_replace_all(cluster, " ", "_")) %>%
#   dplyr::select(-sample_name, -level1, -cluster)
#
# # Remove cells from data table with no class, convert gene identifiers from entrez to gene id, and convert to matrix.
# noclass <- sample_columns %>% dplyr::filter(class == "no class") %>% .[["sample_name"]] # amounts to 325 cells
# exons <- exons %>%
#   inner_join(gene_rows %>%
#                dplyr::select(gene, entrez_id),
#              by = c("V1" = "entrez_id")) %>%
#   dplyr::select(gene, everything(), -V1)
#
# matrix.please<-function(dataframe) {
#   # create matrix
#   m<-as.matrix(dataframe[,-1])
#
#   # add row names using first column of dataframe
#   rownames(m)<-dataframe[,1]
#
#   return(m)
# }
#
# exp <- exons %>%
#   dplyr::select(-one_of(noclass)) %>%
#   matrix.please()
#
# # Create list with metadat and matrix.
# # To use with EWCE, must ensure matrix is named 'exp' and metadat is 'annot'.
# AIBS2019 <- list(exp, metadat)
# names(AIBS2019) <- c("exp", "annot")
#
# save(AIBS2019,file="/home/rreynolds/data/scRNAseq_AIBS/MTG/AIBS2019_DataForEWCE.Rda")

#---EWCE: Calculating specificity matrices-----------------------------------------------------------------------------------------------####
# Sourced using:
# nohup Rscript /home/rreynolds/projects/MarkerGenes/R/20190501_AIBS2019_GeneratingSImatrix.R &>/home/rreynolds/projects/MarkerGenes/nohup_logs/EWCE_AIBS2019.log&
# Drop genes which do not show significant evidence of varying between level 2 celltypes (based on ANOVA)
Sys.time()
load("/home/rreynolds/data/scRNAseq_AIBS/MTG/AIBS2019_DataForEWCE.Rda")
print("Data loaded.")

exp_DROPPED = drop.uninformative.genes(exp=AIBS2019$exp,level2annot = AIBS2019$annot$level2class)
print("Uninformative genes now removed.")

annotLevels = list(level1class=AIBS2019$annot$level1class,level2class=AIBS2019$annot$level2class)
print("Annotation levels assigned.")

# Remove unneccessary files
rm(AIBS2019)
print("AIBS2019 file now removed.")

# Calculate cell type averages and specificity for each gene
fNames_AIBS2019 = generate.celltype.data(exp=exp_DROPPED,annotLevels=annotLevels,groupName="AIBS2019")
print("Cell type averages and specificity calculated.")
