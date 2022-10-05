#---Load Libraries and data--------------------------------------------------------------------------------------------------------------####
library(loomR)
library(tidyverse)
library(EWCE)
library(limma)
library(readxl)
library(tidyverse)
library(stringr)

setwd("~/projects/LDSC_Regression/Data_preLDSC/Zeisel2018/")

# Connect to the loom file in read/write mode
lfile <- connect(filename = "/home/rreynolds/projects/LDSC_Regression/Data_preLDSC/Zeisel2018/l5_all.loom", mode = "r+")

#---WARNING---###
# In order to maintain maximum efficiency and I/O speed, the HDF5 library transposes the access for the underlying data matrix.
# This means that genes are stored in columns, and cells are stored in rows. The boost in speed is well worth it, and we’ve
# taken precaution in developing loomR so that all methods take this into account, and properly transpose data for users.
# However, when interacting with the data matrices directly in loomR, please be sure to keep this in mind. Note also that,
# in order to maintain variable compatibility with loomR, the row.attrs slot still refers to gene metadata, and the
# col.attrs still refers to the cell metadata.

#---Tutorial-----------------------------------------------------------------------------------------------------------------------------####
# # For the sake of consistency within the single-cell community, we've
# # reversed the dimensions for the `shape` field.  As such, the number of
# # genes is stored in `lfile$shape[1]`; the number of cells is stored in the
# # second field
# lfile[["row_attrs/Gene"]]$dims == lfile$shape[1]
#
# # Pull gene expression data for all genes, for the first 5 cells Note that
# # we're using the row position for cells
# data.subset <- lfile[["matrix"]][1:5,1:5 ]
# dim(x = data.subset)
#
# # You can transpose this matrix if you wish to restore the standard
# # orientation
# data.subset <- t(x = data.subset)
# dim(x = data.subset)
#
# # Pull gene expression data. Note that we're using the column position for genes.
# data.gene <- lfile[["matrix"]][1:10, lfile$row.attrs$Gene[] == "Cbln2"| lfile$row.attrs$Gene[] == "Ptchd2"]
# datacell <- lfile[["matrix"]][lfile$col.attrs$CellID[] == "10X82_2_TCTCTCACCAGTTA-", lfile$row.attrs$Gene[] == "Cbln2"| lfile$row.attrs$Gene[] == "Ptchd2"]

#---Main---------------------------------------------------------------------------------------------------------------------------------####
# EWCE requires:
# 1. Matrix with data, with rows as genes and columns as cells.
# 2. Annotation dataframe with a minimum of cell__id, level1class (taxonomy rank 4) and level2class (ClusterName).

# Pull necessary metadata from the column attributes
cell_id <- lfile$col.attrs$CellID[]
Class <- lfile$col.attrs$Class[]
level1class <- lfile$col.attrs$TaxonomyRank4[]
level2class <- lfile$col.attrs$ClusterName[]
region <- lfile$col.attrs$Region[]


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
save(Zeisel2018,file="~/projects/LDSC_Regression/Data_preLDSC/Zeisel2018/Zeisel2018_CNSClusters_DataForEWCE.Rda")

# #---EWCE: Checking data------------------------------------------------------------------------------------------------------------------####
# load("~/projects/LDSC_Regression/Data_preLDSC/Zeisel2018/Zeisel2018_NeuronalClusters_DataForEWCE.Rda")
#
# print("Data now loaded.")

# Function to correct incorrect gene symbols.
Zeisel2018$exp = fix.bad.mgi.symbols(Zeisel2018$exp,
                                     mrk_file_path="/home/rreynolds/projects/LDSC_Regression/Data_preLDSC/Zeisel2015/MRK_List2.rpt")

print("Symbols now fixed.")

save(Zeisel2018, file="~/projects/LDSC_Regression/Data_preLDSC/Zeisel2018/Zeisel2018_NeuronalClusters_DataForEWCE_MGIfixed.Rda")

print("Data with symbols now saved.")

#---EWCE: Calculating specificity matrices-----------------------------------------------------------------------------------------------####
# Commented out, as this only requires performing once, after which generated files can simply be re-loaded.
# Drop genes which do not show significant evidence of varying between level 2 celltypes (based on ANOVA)
# load("~/projects/LDSC_Regression/Data_preLDSC/Zeisel2018/Zeisel2018_AllClusters_DataForEWCE_MGIfixed.Rda")
# print("Data loaded.")

exp_DROPPED = drop.uninformative.genes(exp=Zeisel2018$exp,level2annot = Zeisel2018$annot$level2class)
print("Uninformative genes now removed.")

annotLevels = list(level1class=Zeisel2018$annot$level1class,level2class=Zeisel2018$annot$level2class)
print("Annotation levels assigned.")

# Remove unneccessary files
rm(Zeisel2018)
print("Zeisel2018 file now removed.")

# Calculate cell type averages and specificity for each gene
fNames_Zeisel2018 = generate.celltype.data(exp=exp_DROPPED,annotLevels=annotLevels,groupName="kiZeisel2018_CNSClusters")
print("Cell type averages and specificity calculated.")

# Drop all genes which do not have 1:1 mouse:human orthologs. Only necessary if planning to compare to human data i.e. gene sets
# resulting from human genetics
fNames_Zeisel2018 = filter.genes.without.1to1.homolog(fNames_Zeisel2018)
print("All genes without 1:1 mouse:human orthologs now dropped.")
