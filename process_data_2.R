#########################################
# 0. Load data, packages, set directories
#########################################
library(magrittr)
xie_offsite <-.get_config_path("LOCAL_XIE_2019_DATA_DIR")
raw_data_dir <- paste0(xie_offsite, "raw/")
intermediate_file_dir <- paste0(xie_offsite, "intermediate/")
gRNA_matrix_list <- readRDS(paste0(intermediate_file_dir, "gRNA_matrix_list.rds"))

#####################################################
# 1. create the grouped and ungrouped count matrices;
# grouped has types "summed" and "binary"
#####################################################
# 1. ungrouped matrix
ungrouped_mat_list <- lapply(seq(1L, length(gRNA_matrix_list)), function(i) {
  x <- do.call(what = "cbind", args = gRNA_matrix_list[[i]]$count_matrix_list) %>% Matrix::t()
  colnames(x) <- gRNA_matrix_list[[i]]$cell_barcodes
  return(x)
})
ungrouped_mat <- do.call(what = "cbind", args = ungrouped_mat_list)
ungrouped_mat <- as(ungrouped_mat, "dgTMatrix")
saveRDS(ungrouped_mat, paste0(intermediate_file_dir, "ungrouped_matrix.rds"))

# 2. summed matrix
summed_mat_list <- lapply(seq(1L, length(gRNA_matrix_list)), function(i) {
  x <- do.call(what = "rbind", args = lapply(gRNA_matrix_list[[i]]$count_matrix_list, Matrix::rowSums))
  colnames(x) <- gRNA_matrix_list[[i]]$cell_barcodes
  return(x)
})
summed_mat <- do.call(what = "cbind", args = summed_mat_list)
summed_mat <-  as(summed_mat, "dgTMatrix")
saveRDS(summed_mat, paste0(intermediate_file_dir, "summed_matrix.rds"))

# 3. binary matrix
target_regions <- gRNA_matrix_list[[1]]$count_matrix_list %>% names()
gRNA_count_matrix_list <- lapply(target_regions, function(region) {
  lapply(gRNA_matrix_list, function(x) x$count_matrix_list[[region]]) %>%
    do.call(what = "rbind", args = .)
})
cell_barcodes <- lapply(gRNA_matrix_list, function(x) x$cell_barcodes) %>%
  do.call(what = "c", args = .)
# We reduce each matrix in the gRNA_count_matrix_list to a single logical vector
combine_gRNAs_in_group <- function(gRNA_count_matrix) {
  apply(X = gRNA_count_matrix, MARGIN = 2, FUN = function(column) {
    v <- sum(column >= 1)
    U <- sum(column)
    column/U > 1/v
  }) %>% apply(MARGIN = 1, FUN = function(r) any(r))
}
gRNA_indic_matrix_list <- lapply(gRNA_count_matrix_list, combine_gRNAs_in_group)
gRNA_indic_matrix <- do.call(what = "rbind", args = gRNA_indic_matrix_list)
colnames(gRNA_indic_matrix) <- cell_barcodes
rownames(gRNA_indic_matrix) <- target_regions
  
gRNA_indic_matrix_sparse <- as(gRNA_indic_matrix, "dgTMatrix")
saveRDS(gRNA_indic_matrix_sparse, paste0(intermediate_file_dir, "binary_matrix.rds"))
