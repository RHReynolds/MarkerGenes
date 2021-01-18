## code to prepare `mouse_to_human_orthologs` dataset

# Get all ensembl IDs, with corresponding HGNC_symbol and entrez ID. Using latest (GRCh38.p13).
library(biomaRt)

human <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
hgnc_symbols <- getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "entrezgene_id"), mart = human)

# 15/01/2021 - Downloaded and formatted list of homologues from ensembl BioMart, using query below.
# Query can be updated and saved at any time. Just remember to update download date.
mouse_to_human_orthologs <- getBM(attributes=c("ensembl_gene_id",
                                               "mmusculus_homolog_ensembl_gene",
                                               "mmusculus_homolog_associated_gene_name",
                                               "mmusculus_homolog_orthology_type",
                                               "mmusculus_homolog_orthology_confidence"),
                                  filters = c("with_mmusculus_homolog"),
                                  values = TRUE, mart = human) %>%
  dplyr::left_join(hgnc_symbols, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
  dplyr::select(ensembl_gene_id, hgnc_symbol, entrezgene = entrezgene_id, everything())

# Save data
usethis::use_data(mouse_to_human_orthologs, overwrite = TRUE)
