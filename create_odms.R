#########################################
# 0. Load data, packages, set directories
#########################################
library(magrittr)
xie_offsite <-.get_config_path("LOCAL_XIE_2019_DATA_DIR")
raw_data_dir <- paste0(xie_offsite, "raw/")
intermediate_file_dir <- paste0(xie_offsite, "intermediate/")
processed_data_dir <- paste0(xie_offsite, "processed/")
gRNA_processed_data_dir <- paste0(processed_data_dir, "gRNA/")
if (!dir.exists(gRNA_processed_data_dir)) dir.create(gRNA_processed_data_dir, recursive = TRUE)

#####################
# 1. create gRNA odms
#####################
# load gRNA information from raw directory
guide_seqs <- openxlsx::read.xlsx(xlsxFile = paste0(raw_data_dir, "/all_oligos.xlsx"),
                                  sheet = 1) %>% dplyr::select(spacer_sequence = "spacer.sequence",
                                                               hg38_enh_region = "region.pos.(hg38)")
# 1. ungrouped-raw
ungrouped_gRNA_mat <- readRDS(paste0(intermediate_file_dir, "ungrouped_matrix.rds"))
barcodes <- colnames(ungrouped_gRNA_mat)
spacer_seqs <- row.names(ungrouped_gRNA_mat)
# confirm spacer seq order in matrix matches that of data frame
all(spacer_seqs == guide_seqs$spacer_sequence)
odm_fp <- paste0(gRNA_processed_data_dir, "raw_ungrouped.odm")
metadata_fp <- paste0(gRNA_processed_data_dir, "raw_ungrouped_metadata")
odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = ungrouped_gRNA_mat,
                                                  barcodes = barcodes,
                                                  features_df = guide_seqs,
                                                  odm_fp = odm_fp,
                                                  metadata_fp = metadata_fp)
rm(ungrouped_gRNA_mat, barcodes, spacer_seqs, odm_fp, metadata_fp, odm)

# 2. ungrouped-binary
ungrouped_gRNA_mat <- readRDS(paste0(intermediate_file_dir, "binary_matrix_ungrouped.rds"))
barcodes <- colnames(ungrouped_gRNA_mat)
spacer_seqs <- rownames(ungrouped_gRNA_mat)
all(spacer_seqs == guide_seqs$spacer_sequence)
odm_fp <- paste0(gRNA_processed_data_dir, "binary_ungrouped.odm")
metadata_fp <- paste0(gRNA_processed_data_dir, "binary_ungrouped_metadata.rds")
odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = ungrouped_gRNA_mat,
                                                  barcodes = barcodes,
                                                  features_df = guide_seqs,
                                                  odm_fp = odm_fp,
                                                  metadata_fp = metadata_fp)
rm(ungrouped_gRNA_mat, barcodes, spacer_seqs, odm_fp, metadata_fp, odm)

# 3. grouped-raw
grouped_raw <- readRDS(paste0(intermediate_file_dir, "summed_matrix.rds"))
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
barcodes <- colnames(grouped_binary)
chrom_ids <- rownames(grouped_binary)
odm_fp <- paste0(gRNA_processed_data_dir, "binary_grouped.odm")
metadata_fp <- paste0(gRNA_processed_data_dir, "binary_grouped_metadata.rds")
odm <- ondisc::create_ondisc_matrix_from_R_matrix(r_matrix = grouped_binary,
                                                  barcodes = barcodes,
                                                  features_df = data.frame(chrom_ids),
                                                  odm_fp = odm_fp,
                                                  metadata_fp = metadata_fp)
