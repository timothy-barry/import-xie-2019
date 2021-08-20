# This script converts the gene expression data into a sparse matrix.
# It save the sparse matrix into the intermediate file directory.

###################################
# 0. Load packages, set directories
###################################
xie_offsite <-.get_config_path("LOCAL_XIE_2019_DATA_DIR")
raw_data_dir <- paste0(xie_offsite, "raw/")
intermediate_file_dir <- paste0(xie_offsite, "intermediate/")
library(magrittr)

#######################
# get fps to hdf5 files
#######################
fs <- list.files(raw_data_dir)
h5_fs <- grep(pattern = "*.h5$", x = fs, value = TRUE)
batch <- stringr::str_extract(string = h5_fs, pattern = "Batch-[0-9]_[0-9]") %>%
  gsub(pattern = "Batch-", replacement = "")
base_name <- "/refgenome_hg38_CROP-Guide-MS2-2.1.0/"

# loop through the files, converting the h5 files to sparse matrix format and obtaining the number of cells
matrix_list <- lapply(seq(1L, length(h5_fs)), function(i) {
  print(paste0("Processing ", i))
  h5_f <- paste0(raw_data_dir, h5_fs[i])
  # get the components of the sparse matrix
  rhdf5::h5ls(h5_f)
  barcodes <- rhdf5::h5read(file = h5_f, name = paste0(base_name, "barcodes"))
  # ensure barcodes are unique
  if (any(duplicated(barcodes))) stop("Barcodes duplicated within batch/run.")
  # process barcodes
  barcodes_annotated <- paste0(gsub(pattern = "1", replacement = "", x = barcodes),
                               batch[i])
  gene_ids <- rhdf5::h5read(file = h5_f, name = paste0(base_name, "genes"))
  x <- rhdf5::h5read(file = h5_f, name = paste0(base_name, "data"))
  row_idxs <- rhdf5::h5read(file = h5_f, name = paste0(base_name, "indices"))
  col_ptr <- rhdf5::h5read(file = h5_f, name = paste0(base_name, "indptr"))
  dim <- rhdf5::h5read(file = h5_f, name = paste0(base_name, "shape"))
  # construct the sparse matrix
  m <- Matrix::sparseMatrix(i = row_idxs, p = col_ptr, x = x, dims = dim)
  row.names(m) <- gene_ids
  colnames(m) <- barcodes_annotated
  return(m)
})

# create the gene expression matrix
gene_exp_m <- do.call(what = "cbind", args = matrix_list)
any(duplicated(colnames(gene_exp_m)))

# compute the batch factor vector
n_cells_per_file <- sapply(matrix_list, function(m) ncol(m))
batch <- rep(x = batch, times = n_cells_per_file)
names(batch) <- colnames(gene_exp_m)

saveRDS(object = batch, file = paste0(intermediate_file_dir, "batch_vector.rds"))
saveRDS(object = gene_exp_m, paste0(intermediate_file_dir, "gene_exp_mat.rds"))

# Finally, obtain the cells to use (taken as the intersection of the gene cells and the gRNA cells)
# gene_exp_m <- readRDS(paste0(intermediate_file_dir, "gene_exp_mat.rds"))
gRNA_matrix <- readRDS(paste0(intermediate_file_dir, "summed_matrix.rds"))
cells_to_use <- intersect(colnames(gene_exp_m), colnames(gRNA_matrix))
saveRDS(object = cells_to_use, paste0(intermediate_file_dir, "cells_in_use.rds"))
