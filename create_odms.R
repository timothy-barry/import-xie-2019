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

############################
# Create gene expression odm
############################
h5_loc <- paste0(raw_data_dir, "GSM3722729_K562-dCas9-KRAB_5K-sgRNAs_Batch-1_1_filtered_gene_bc_matrices_h5.h5")
gene_id <- rhdf5::h5read(file = h5_loc, name = "/refgenome_hg38_CROP-Guide-MS2-2.1.0/genes")
gene_name <- rhdf5::h5read(file = h5_loc, name = "/refgenome_hg38_CROP-Guide-MS2-2.1.0/gene_names")
cells_in_use <- readRDS(paste0(intermediate_file_dir, "cells_in_use.rds"))

exp_mat <- readRDS(paste0(intermediate_file_dir, "gene_exp_mat.rds"))[,cells_in_use]
exp_mat <- as(exp_mat, "dgTMatrix")
# check that gene IDs are correctly ordered
all(rownames(exp_mat) == gene_id)
features_df <- data.frame(gene_id = gene_id, gene_name = gene_name)
odm_fp <- paste0(gene_processed_data_dir, "expression_matrix.odm")
metadata_fp <- paste0(gene_processed_data_dir, "metadata.rds")

odm <- create_ondisc_matrix_from_R_matrix(r_matrix = exp_mat,
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
guide_seqs <- readRDS(paste0(intermediate_file_dir, "guide_seqs.rds")) %>%
  dplyr::select(spacer_seq, hg38_enh_region)

# 1. ungrouped-raw
ungrouped_gRNA_mat <- readRDS(paste0(intermediate_file_dir, "gRNA_matrix.rds"))[,cells_in_use]
ungrouped_gRNA_mat <- as(ungrouped_gRNA_mat, "dgTMatrix")
barcodes <- colnames(ungrouped_gRNA_mat)
odm_fp <- paste0(gRNA_processed_data_dir, "raw_ungrouped.odm")
metadata_fp <- paste0(gRNA_processed_data_dir, "raw_ungrouped_metadata")
odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = ungrouped_gRNA_mat,
                                                  barcodes = barcodes,
                                                  features_df = guide_seqs,
                                                  odm_fp = odm_fp,
                                                  metadata_fp = metadata_fp)
rm(ungrouped_gRNA_mat, barcodes, odm_fp, metadata_fp, odm)

# 2. ungrouped-binary
ungrouped_gRNA_mat_bin <- readRDS(paste0(intermediate_file_dir, "thresholded_ungrouped_mat.rds"))[,cells_in_use]
ungrouped_gRNA_mat_bin <- as.matrix(ungrouped_gRNA_mat_bin) == 1
barcodes <- colnames(ungrouped_gRNA_mat_bin)
odm_fp <- paste0(gRNA_processed_data_dir, "binary_ungrouped.odm")
metadata_fp <- paste0(gRNA_processed_data_dir, "binary_ungrouped_metadata.rds")
odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = ungrouped_gRNA_mat_bin,
                                                  barcodes = barcodes,
                                                  features_df = guide_seqs,
                                                  odm_fp = odm_fp,
                                                  metadata_fp = metadata_fp)
rm(ungrouped_gRNA_mat_bin, barcodes, odm_fp, metadata_fp, odm)

# 3. grouped-raw
grouped_raw <- readRDS(paste0(intermediate_file_dir, "unthresholded_grouped_gRNA.rds"))[,cells_in_use]
grouped_raw <- as(grouped_raw, "dgTMatrix")
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
grouped_binary <- readRDS(paste0(intermediate_file_dir, "thresholded_grouped_mat_gRNA.rds"))[,cells_in_use]
grouped_binary <- as.matrix(grouped_binary) == 1
barcodes <- colnames(grouped_binary)
chrom_ids <- rownames(grouped_binary)
odm_fp <- paste0(gRNA_processed_data_dir, "binary_grouped.odm")
metadata_fp <- paste0(gRNA_processed_data_dir, "binary_grouped_metadata.rds")
odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = grouped_binary,
                                                  barcodes = barcodes,
                                                  features_df = data.frame(chrom_ids),
                                                  odm_fp = odm_fp,
                                                  metadata_fp = metadata_fp)
