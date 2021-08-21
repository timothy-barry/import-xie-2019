################################################################
# 0. set file paths, load packages and gene IDs / gRNA positions
################################################################
library(magrittr)
xie_offsite <-.get_config_path("LOCAL_XIE_2019_DATA_DIR")
raw_data_dir <- paste0(xie_offsite, "raw/")
aux_data_dir <- paste0(xie_offsite, "processed/aux/")
gene_df <- openxlsx::read.xlsx(xlsxFile = paste0(raw_data_dir, "/Genes.xlsx"), sheet = 1)
all_protein_coding_genes <- gene_df$Gene_Symbol


#########################################
# 1. Load gRNA identification information
#########################################
enh_targets_df <- openxlsx::read.xlsx(xlsxFile = paste0(raw_data_dir, "/enh_targets.xlsx"), sheet = 1)
bulk_region_names <- dplyr::filter(enh_targets_df, gene_names %in% c("ARL15", "MYB")) %>%
  dplyr::select(region, region_name = Denoted.Region.Name.used.in.the.paper, targeted_gene = gene_names) %>%
  dplyr::filter(region_name %in% c("ARL15-enh", "MYB-enh-3"))
saveRDS(object = bulk_region_names, paste0(aux_data_dir, "/bulk_region_names.rds"))

##############
# Bulk RNA-seq
##############
bulk_info <- openxlsx::read.xlsx(xlsxFile = paste0(raw_data_dir, "/bulk_rna_info.xlsx"), sheet = 3) %>% dplyr::select(library_name = Library.Name, gRNA = sgRNA, region = Region, biological_duplicate = Biological.Duplicate)
bulk_info_arl15_enh <- dplyr::slice(bulk_info, 1:25)
bulk_info_myb_enh3 <- dplyr::slice(bulk_info, 26:49)

bulk_df <- suppressWarnings(readr::read_tsv(file = paste0(raw_data_dir, "/GSE129825_Libraries.FeatureCounts.ARL15_enhancer.txt"),
                                            col_types = "ccccciiiiiiiiiiiiiiiiiiiiiiiiii")) %>% dplyr::rename("PZ788" = "PZ778...10", "PZ778" = "PZ778...13")
bulk_df_arl15_enh <- dplyr::filter(bulk_df, Geneid %in% all_protein_coding_genes)
bulk_df <- suppressWarnings(readr::read_tsv(file = paste0(raw_data_dir, "/GSE129825_Libraries.FeatureCounts.MYB_enhancer.txt"),
                                            col_types = "ccccciiiiiiiiiiiiiiiiiiiiiiiiii"))
bulk_df_myb_enh3 <- dplyr::filter(bulk_df, Geneid %in% all_protein_coding_genes)

bulk_rnaseq <- list(data = list(arl15_enh = bulk_df_arl15_enh, myb_enh3 = bulk_df_myb_enh3), info = list(arl15_enh = bulk_info_arl15_enh, myb_enh3 = bulk_info_myb_enh3))
saveRDS(object = bulk_rnaseq, file = paste0(aux_data_dir, "bulk_RNAseq.rds"))
