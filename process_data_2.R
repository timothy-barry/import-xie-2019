###################################
# 0. Load packages, set directories
###################################
xie_offsite <-.get_config_path("LOCAL_XIE_2019_DATA_DIR")
raw_data_dir <- paste0(xie_offsite, "raw/")
intermediate_file_dir <- paste0(xie_offsite, "intermediate/")
if (!dir.exists(intermediate_file_dir)) dir.create(intermediate_file_dir, recursive = TRUE)
library(magrittr)

gRNA_matrix_raw <- readRDS(paste0(intermediate_file_dir, "gRNA_matrix.rds"))
guide_seqs <- readRDS(paste0(intermediate_file_dir, "guide_seqs.rds"))
unique_regions <- unique(guide_seqs$hg38_enh_region)

########################
# ungrouped, thresholded
########################
thresholded_ungrouped_mat <- apply(X = gRNA_matrix_raw, MARGIN = 1, FUN = function(r) {
  r_nonzero <- r[r >= 1]
  as.integer(r > mean(r_nonzero))
}) %>% Matrix::t()
colnames(thresholded_ungrouped_mat) <- colnames(gRNA_matrix_raw)

#################################
# grouped, unthresholded (summed)
#################################
unthresholded_grouped_mat <- sapply(X = unique_regions, FUN = function(curr_region) {
  print(curr_region)
  curr_spacer_seqs <- dplyr::filter(guide_seqs, hg38_enh_region == curr_region) %>%
    dplyr::pull(spacer_seq)
  s <- gRNA_matrix_raw[curr_spacer_seqs,] %>% Matrix::colSums()
}) %>% t()

######################
# grouped, thresholded
######################
thresholded_grouped_mat <- sapply(X = unique_regions, FUN = function(curr_region) {
  print(curr_region)
  curr_spacer_seqs <- dplyr::filter(guide_seqs, hg38_enh_region == curr_region) %>%
    dplyr::pull(spacer_seq)
  curr_m <- thresholded_ungrouped_mat[curr_spacer_seqs,]
  apply(X = curr_m, MARGIN = 2, FUN = function(col) max(col))
}) %>% t()


######
# save
######
saveRDS(thresholded_ungrouped_mat, paste0(intermediate_file_dir, "thresholded_ungrouped_mat.rds"))
saveRDS(unthresholded_grouped_mat, paste0(intermediate_file_dir, "unthresholded_grouped_gRNA.rds"))
saveRDS(thresholded_grouped_mat, paste0(intermediate_file_dir, "thresholded_grouped_mat_gRNA.rds"))
