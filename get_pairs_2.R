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

# combine, convert to factor
pairs <- rbind(cis_pairs, neg_control_pairs) %>%
  dplyr::mutate_all(factor) %>%
  dplyr::rename("gene_id" = "gene.id", "gRNA_id" = "gRNA.id")
saveRDS(object = pairs, file = paste0(aux_dir, "pairs_grouped.rds"))

# finally, add the spacer sequences to create the ungrouped pairs data frame
guide_seqs <- readRDS(paste0(aux_dir, "gRNA_sequence_dictionary.rds")) %>% 
  dplyr::rename("enh_region" = "hg38_enh_region") %>% dplyr::mutate_all(factor)
grouped_pairs <- pairs %>% dplyr::rename("enh_region" = "gRNA_id")
ungrouped_pairs <- dplyr::inner_join(x = grouped_pairs, y = guide_seqs, by = "enh_region")
saveRDS(object = ungrouped_pairs, file = paste0(aux_dir, "pairs_ungrouped.rds"))
