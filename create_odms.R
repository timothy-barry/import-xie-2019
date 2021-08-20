#########################################
# 0. Load data, packages, set directories
#########################################
library(magrittr)
library(ondisc)
xie_offsite <-.get_config_path("LOCAL_XIE_2019_DATA_DIR")
raw_data_dir <- paste0(xie_offsite, "raw/")
intermediate_file_dir <- paste0(xie_offsite, "intermediate/")
processed_data_dir <- paste0(xie_offsite, "processed/")
gRNA_processed_data_dir <- paste0(processed_data_dir, "gRNA/")
gene_processed_data_dir <- paste0(processed_data_dir, "gene/")
aux_data_dir <- paste0(processed_data_dir, "aux/")
if (!dir.exists(gRNA_processed_data_dir)) dir.create(gRNA_processed_data_dir, recursive = TRUE)
if (!dir.exists(gene_processed_data_dir)) dir.create(gene_processed_data_dir, recursive = TRUE)
cells_in_use <- readRDS(paste0(intermediate_file_dir, "cells_in_use.rds"))

############################
# Create gene expression odm
############################
exp_mat <- readRDS(paste0(intermediate_file_dir, "gene_exp_mat.rds"))[,cells_in_use]
exp_mat <- as(exp_mat, "dgTMatrix")
features_df <- data.frame(gene_id = rownames(exp_mat))
odm_fp <- paste0(gene_processed_data_dir, "expression_matrix.odm")
metadata_fp <- paste0(gene_processed_data_dir, "metadata.rds")

odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = exp_mat,
                                                  barcodes = cells_in_use,
                                                  features_df = features_df,
                                                  odm_fp = odm_fp,
                                                  metadata_fp = metadata_fp)

# odm <- read_odm(odm_fp = odm_fp, metadata_fp = metadata_fp)
batch <- readRDS(file = paste0(intermediate_file_dir, "batch_vector.rds"))[cells_in_use]
batch_relevel <- factor(x = batch, levels = c("1_1", "1_2", "2_1", "2_2",
                                              "3_1", "3_2", "4_1", "4_2",
                                              "5_1", "5_2"),
                        labels = c("batch_1", "batch_1", "batch_2", "batch_2",
                                   "batch_3", "batch_3", "batch_4", "batch_4",
                                   "batch_5", "batch_5"))
odm_with_batch <- odm %>% mutate_cell_covariates(batch = batch_relevel)
save_odm(odm = odm_with_batch, metadata_fp = paste0(gene_processed_data_dir, "metadata.rds"))
rm(exp_mat, features_df, odm_fp, metadata_fp, odm, batch, batch_relevel, odm_with_batch)

#####################
# 2. create gRNA odms
#####################
# load gRNA information from raw directory
guide_seqs <- openxlsx::read.xlsx(xlsxFile = paste0(raw_data_dir, "/all_oligos.xlsx"),
                                  sheet = 1) %>% dplyr::select(spacer_sequence = "spacer.sequence",
                                                               hg38_enh_region = "region.pos.(hg38)")
# remove "chr11:5280670-5280820" and "chr11:61834748-61834898," which are weird.
guide_seqs <- dplyr::filter(guide_seqs, !(hg38_enh_region %in% c("chr11:5280670-5280820", "chr11:61834748-61834898")))

# 1. ungrouped-raw
ungrouped_gRNA_mat <- readRDS(paste0(intermediate_file_dir, "ungrouped_matrix.rds"))
ungrouped_gRNA_mat <- ungrouped_gRNA_mat[guide_seqs$spacer_sequence, cells_in_use]

barcodes <- colnames(ungrouped_gRNA_mat)
odm_fp <- paste0(gRNA_processed_data_dir, "raw_ungrouped.odm")
metadata_fp <- paste0(gRNA_processed_data_dir, "raw_ungrouped_metadata")
odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = ungrouped_gRNA_mat,
                                                  barcodes = barcodes,
                                                  features_df = guide_seqs,
                                                  odm_fp = odm_fp,
                                                  metadata_fp = metadata_fp)
rm(ungrouped_gRNA_mat, barcodes, spacer_seqs, odm_fp, metadata_fp, odm)

# 2. ungrouped-binary
ungrouped_gRNA_mat_bin <- readRDS(paste0(intermediate_file_dir, "binary_matrix_ungrouped.rds"))
ungrouped_gRNA_mat_bin <- ungrouped_gRNA_mat_bin[guide_seqs$spacer_sequence, cells_in_use]
barcodes <- colnames(ungrouped_gRNA_mat_bin)
odm_fp <- paste0(gRNA_processed_data_dir, "binary_ungrouped.odm")
metadata_fp <- paste0(gRNA_processed_data_dir, "binary_ungrouped_metadata.rds")
odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = ungrouped_gRNA_mat_bin,
                                                  barcodes = barcodes,
                                                  features_df = guide_seqs,
                                                  odm_fp = odm_fp,
                                                  metadata_fp = metadata_fp)
rm(ungrouped_gRNA_mat_bin, barcodes, spacer_seqs, odm_fp, metadata_fp, odm)

# 3. grouped-raw
grouped_raw <- readRDS(paste0(intermediate_file_dir, "summed_matrix.rds"))
grouped_raw <- grouped_raw[unique(guide_seqs$hg38_enh_region), cells_in_use]
barcodes <- colnames(grouped_raw)
chrom_ids <- rownames(grouped_raw)
odm_fp <- paste0(gRNA_processed_data_dir, "raw_grouped.odm")
metadata_fp <- paste0(gRNA_processed_data_dir, "raw_grouped_metadata.rds")
odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = grouped_raw,
                                                  barcodes = barcodes,
                                                  features_df = data.frame(chrom_ids),
                                                  odm_fp = odm_fp,
                                                  metadata_fp = metadata_fp)
rm(grouped_raw, barcodes, chrom_ids, odm_fp, metadata_fp, odm)

# 4. grouped-binary
grouped_binary <- readRDS(paste0(intermediate_file_dir, "binary_matrix_grouped.rds"))
grouped_binary <- grouped_binary[unique(guide_seqs$hg38_enh_region), cells_in_use]
barcodes <- colnames(grouped_binary)
chrom_ids <- rownames(grouped_binary)
odm_fp <- paste0(gRNA_processed_data_dir, "binary_grouped.odm")
metadata_fp <- paste0(gRNA_processed_data_dir, "binary_grouped_metadata.rds")
odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = grouped_binary,
                                                  barcodes = barcodes,
                                                  features_df = data.frame(chrom_ids),
                                                  odm_fp = odm_fp,
                                                  metadata_fp = metadata_fp)

###############################
# save gRNA sequence dictionary
###############################
saveRDS(object = guide_seqs, file = paste0(aux_data_dir, "gRNA_sequence_dictionary.rds"))
