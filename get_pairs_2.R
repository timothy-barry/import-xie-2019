library(magrittr)
xie_offsite <- .get_config_path("LOCAL_XIE_2019_DATA_DIR")
intermediate_data_dir <- paste0(xie_offsite, "intermediate/")
raw_data_dir <- paste0(xie_offsite, "raw/")
aux_dir <- paste0(xie_offsite, "processed/aux/")

# load the cis pairs
cis_pairs <- readRDS(paste0(intermediate_data_dir, "select_gRNA_gene_pair.rds")) %>%
  dplyr::mutate(type = "cis")

# load the negative control pairs
neg_control_pairs <- readRDS(paste0(intermediate_data_dir, "neg_control_pair.rds")) %>%
  dplyr::mutate(type = "neg_control")
neg_control_pairs$gene.hgnc.id <- NULL

# load the list of protein-coding genes, and gene names
gene_df <- openxlsx::read.xlsx(xlsxFile = paste0(raw_data_dir, "/Genes.xlsx"), sheet = 1)
all_protein_coding_genes <- gene_df$Gene_Symbol

h5_loc <- paste0(raw_data_dir, "GSM3722729_K562-dCas9-KRAB_5K-sgRNAs_Batch-1_1_filtered_gene_bc_matrices_h5.h5")
gene_id <- rhdf5::h5read(file = h5_loc, name = "/refgenome_hg38_CROP-Guide-MS2-2.1.0/genes")
gene_name <- rhdf5::h5read(file = h5_loc, name = "/refgenome_hg38_CROP-Guide-MS2-2.1.0/gene_names")
protein_coding <- gene_name %in% all_protein_coding_genes

extra_gene_info <- data.frame(gene_id = gene_id, gene_name = gene_name, protein_coding = protein_coding)

# combine, convert to factor
pairs <- rbind(cis_pairs, neg_control_pairs) %>%
  dplyr::rename("gene_id" = "gene.id", "gRNA_id" = "gRNA.id")
pairs_plus <- dplyr::left_join(x = pairs, y = extra_gene_info, by = "gene_id") %>%
  dplyr::mutate_at(.tbl = ., .vars = c("gene_id", "gRNA_id", "type", "gene_name"), .funs = factor)
saveRDS(object = pairs_plus, file = paste0(aux_dir, "pairs_grouped.rds"))

# finally, add the spacer sequences to create the ungrouped pairs data frame
guide_seqs <- readRDS(paste0(aux_dir, "gRNA_sequence_dictionary.rds")) %>% 
  dplyr::rename("enh_region" = "hg38_enh_region") %>% dplyr::mutate_all(factor)
grouped_pairs <- pairs_plus %>% dplyr::rename("enh_region" = "gRNA_id")
ungrouped_pairs <- dplyr::inner_join(x = grouped_pairs, y = guide_seqs, by = "enh_region") %>%
  dplyr::rename("gRNA_id" = "spacer_sequence")
saveRDS(object = ungrouped_pairs, file = paste0(aux_dir, "pairs_ungrouped.rds"))
