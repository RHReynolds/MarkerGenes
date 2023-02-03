# Script to make new ctd expression data from Zeisel 2018 -------------------------------------------------------------------------------####
# 1. Need to make expression x cell type table for input into make CTD for new version of EWCE ---#
# 2. Convert expression data to ctd under new EWCE version functions

#---Load Libraries and data--------------------------------------------------------------------------------------------------------------####
library(loomR)# needed to open expression data from Zeisel 2018
library(EWCE)
library(tidyverse)
library(limma)
library(readxl)
library(tidyverse)
library(sctransform)
library(stringr)
# loom file from: http://mousebrain.org/adolescent/downloads.html
lfile <- connect(filename = "./ctd/l5_all.loom", mode = "r+")

#---Main---------------------------------------------------------------------------------------------------------------------------------####
# EWCE requires:
# (i). Matrix with data, with rows as genes and columns as cells.
# (ii). Annotation dataframe with a minimum of cell__id, level1class (taxonomy rank 4) and level2class (ClusterName).

# 1. Form expression x cell type data
# Pull necessary metadata from the column attributes
cell_id <- lfile$col.attrs$CellID[]
Class <- lfile$col.attrs$Class[]
level1class <- lfile$col.attrs$TaxonomyRank4[]
level2class <- lfile$col.attrs$ClusterName[]
region <- lfile$col.attrs$Region[] # All CNS regions


# Create data frame from attributes
metadat <- as.data.frame(cell_id) %>%
  bind_cols(as.data.frame(Class)) %>%
  bind_cols(as.data.frame(level1class)) %>%
  bind_cols(as.data.frame(level2class)) %>%
  bind_cols(as.data.frame(region))

# Filter according to desired taxonomy level 4 cells
metadat_filtered_level1 <- metadat %>%
  filter(!Class == "PeripheralGlia") %>%
  filter(!Class == "Vascular") %>%
  filter(!Class == "Ependymal") %>%
  filter(!level1class == "Oligodendrocyte precursor cells") %>%
  filter(!level1class == "Spinal cord inhibitory neurons") %>%
  filter(!level1class == "Spinal cord excitatory neurons") %>%
  filter(!level1class == "Hindbrain neurons") %>%
  filter(!level1class == "Di- and mesencephalon excitatory neurons") %>%
  filter(!level1class == "Di- and mesencephalon inhibitory neurons") %>%
  filter(!level1class == "Peripheral sensory peptidergic neurons") %>%
  filter(!level1class == "Peripheral sensory neurofilament neurons") %>%
  filter(!level1class == "Peripheral sensory non-peptidergic neurons") %>%
  filter(!level1class == "Enteric neurons") %>%
  filter(!level1class == "Sympathetic cholinergic neurons") %>%
  filter(!level1class == "Sympathetic noradrenergic neurons") %>%
  filter(!level1class == "Olfactory inhibitory neurons") %>% # No olfactory bulb in humans, so don't include
  filter(!level1class == "Olfactory ensheathing cells") %>%
  filter(!level1class == "Dentate gyrus radial glia-like cells") %>%
  filter(!level1class == "Subventricular zone radial glia-like cells") %>%
  filter(!level1class == "Enteric glia") %>%
  filter(!level2class == "MOL3") %>%
  filter(!level2class == "MOL2") %>%
  filter(!level2class == "NFOL2") %>%
  filter(!level2class == "NFOL1") %>%
  filter(!level2class == "COP2")

# Pull gene names for naming of matrix
Gene.names <-  lfile$row.attrs$Gene[]
CellID.names <- metadat_filtered_level1 %>% .[["cell_id"]]

# Pull gene expression for CNS clusters run through S-LDSC
# data.subset <- lfile[["matrix"]][lfile$col.attrs$TaxonomyRank4[] == "Oligodendrocytes" |lfile$col.attrs$TaxonomyRank4[] == "Cholinergic and monoaminergic neurons" | lfile$col.attrs$TaxonomyRank4[] == "Telencephalon projecting excitatory neurons" | lfile$col.attrs$TaxonomyRank4[] == "Telencephalon inhibitory interneurons" | lfile$col.attrs$TaxonomyRank4[] == "Peptidergic neurons" | lfile$col.attrs$TaxonomyRank4[] == "Di- and mesencephalon excitatory neurons" | lfile$col.attrs$TaxonomyRank4[] == "Glutamatergic neuroblasts" | lfile$col.attrs$TaxonomyRank4[] == "Hindbrain neurons" | lfile$col.attrs$TaxonomyRank4[] == "Spinal cord excitatory neurons" | lfile$col.attrs$TaxonomyRank4[] == "Telencephalon projecting inhibitory neurons" | lfile$col.attrs$TaxonomyRank4[] == "Non-glutamatergic neuroblasts" | lfile$col.attrs$TaxonomyRank4[] == "Dentate gyrus radial glia-like cells" | lfile$col.attrs$TaxonomyRank4[] == "Subventricular zone radial glia-like cells" | lfile$col.attrs$TaxonomyRank4[] == "Oligodendrocyte precursor cells" | lfile$col.attrs$TaxonomyRank4[] == "Ependymal cells" | lfile$col.attrs$TaxonomyRank4[] == "Subcommissural organ hypendymal cells" | lfile$col.attrs$TaxonomyRank4[] == "Dentate gyrus granule neurons" | lfile$col.attrs$TaxonomyRank4[] == "Cerebellum neurons" | lfile$col.attrs$TaxonomyRank4[] == "Di- and mesencephalon inhibitory neurons" | lfile$col.attrs$TaxonomyRank4[] == "Spinal cord inhibitory neurons" | lfile$col.attrs$TaxonomyRank4[] == "Vascular and leptomeningeal cells" | lfile$col.attrs$TaxonomyRank4[] == "Vascular smooth muscle cells" | lfile$col.attrs$TaxonomyRank4[] == "Pericytes" | lfile$col.attrs$TaxonomyRank4[] == "Vascular endothelial cells" | lfile$col.attrs$TaxonomyRank4[] == "Microglia" | lfile$col.attrs$TaxonomyRank4[] == "Perivascular macrophages" | lfile$col.attrs$TaxonomyRank4[] == "Astrocytes" | lfile$col.attrs$TaxonomyRank4[] == "Choroid epithelial cells", ]
data.subset <- lfile[["matrix"]][lfile$col.attrs$TaxonomyRank4[] == "Cholinergic and monoaminergic neurons" | lfile$col.attrs$TaxonomyRank4[] == "Telencephalon projecting excitatory neurons" | lfile$col.attrs$TaxonomyRank4[] == "Telencephalon inhibitory interneurons" | lfile$col.attrs$TaxonomyRank4[] == "Peptidergic neurons" | lfile$col.attrs$TaxonomyRank4[] == "Glutamatergic neuroblasts" | lfile$col.attrs$TaxonomyRank4[] == "Telencephalon projecting inhibitory neurons" | lfile$col.attrs$TaxonomyRank4[] == "Non-glutamatergic neuroblasts" | lfile$col.attrs$TaxonomyRank4[] == "Dentate gyrus granule neurons" | lfile$col.attrs$TaxonomyRank4[] == "Cerebellum neurons" | lfile$col.attrs$TaxonomyRank4[] == "Microglia" | lfile$col.attrs$TaxonomyRank4[] == "Perivascular macrophages" | lfile$col.attrs$TaxonomyRank4[] == "Astrocytes" | lfile$col.attrs$ClusterName[] == "MFOL1" | lfile$col.attrs$ClusterName[] == "MFOL2" | lfile$col.attrs$ClusterName[] == "MOL1" | lfile$col.attrs$ClusterName[] == "COP1", ]

# Rename rows and columns. Recall that genes are columns and cells are rows.
rownames(data.subset) <- CellID.names
colnames(data.subset) <- Gene.names

# Transpose matrix, such that rows become genes and columns are cells.
data.subset <- t(data.subset)

# Create list with metadat and matrix.
# To use with EWCE, must ensure matrix is named 'exp' and metadat is 'annot'.
Zeisel2018 <- list(data.subset, metadat_filtered_level1)
names(Zeisel2018) <- c("exp", "annot")

# Close loom object
lfile$close_all()

# Save file
save(Zeisel2018,file="./Zeisel2018_CNSClusters_DataForEWCE.Rda")

# 2. Generating ctd using EWCE ---###

load("./Zeisel2018_CNSClusters_DataForEWCE.Rda")

# Correct gene symbols
Zeisel2018$exp = fix_bad_mgi_symbols(Zeisel2018$exp)

save(Zeisel2018, file="./Zeisel2018_CNSClusters_DataForEWCE_MGIfixed.Rda")

load("./Zeisel2018_CNSClusters_DataForEWCE_MGIfixed.Rda")

# Calculate specificity matrices

## drops uninformative genes in order to reduce compute time and noise in subsequent steps. 
# It achieves this through several steps, each of which are optional: 1. Drop non-1:1 orthologs: Removes genes that don’t have 1:1 orthologs with the output_species 
# 2. Drop non-varying genes: Removes genes that don’t vary across cells based on variance deciles. 
#3. Drop non-differentially expressed genes (DEGs): Removes genes that are not significantly differentially expressed across cell-types. 

exp_DROPPED <- EWCE::drop_uninformative_genes(
  exp = Zeisel2018$exp,
  input_species = "mouse",
  output_species = "human",
  level2annot = Zeisel2018$annot$level2class,
  no_cores = 20) 

# Generate ctd
annotLevels <- list(level1class=Zeisel2018$annot$level1class,
                    level2class=Zeisel2018$annot$level2class)

fNames_CNS <- EWCE::generate_celltype_data(
  exp = exp_DROPPED,
  annotLevels = annotLevels,
  groupName = "CNS") 

ctd_CNS <- EWCE::load_rdata(fNames_CNS)

# Check working
EWCE::plot_ctd(ctd = ctd_CNS,
               level = 1,
               genes = c("Apoe","Gfap","Gapdh"),
               metric = "mean_exp") 


