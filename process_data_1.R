# Process data 1 & 2 convert the gRNA data into sparse matrices, stored in intermediate dir.

###################################
# 0. Load packages, set directories
###################################
xie_offsite <-.get_config_path("LOCAL_XIE_2019_DATA_DIR")
raw_data_dir <- paste0(xie_offsite, "raw/")
intermediate_file_dir <- paste0(xie_offsite, "intermediate/")
if (!dir.exists(intermediate_file_dir)) dir.create(intermediate_file_dir, recursive = TRUE)
library(magrittr)

####################
# 1. gRNA UMI counts
####################
# remove the weird spacers
guide_seqs <- openxlsx::read.xlsx(xlsxFile = paste0(raw_data_dir, "/all_oligos.xlsx"),
                                sheet = 1) %>% dplyr::select(hg38_enh_region = "region.pos.(hg38)",
                                                             spacer_seq = "spacer.sequence")
weird_spacers <- names(which(table(guide_seqs$spacer_seq) >= 2))
guide_seqs <- dplyr::filter(guide_seqs, !(spacer_seq %in% weird_spacers))
saveRDS(object = guide_seqs, file = paste0(intermediate_file_dir, "guide_seqs.rds"))
all(table(guide_seqs$spacer_seq) == 1); all(table(guide_seqs$hg38_enh_region) == 10)

# get the raw file names
raw_fs <- list.files(raw_data_dir)
gRNA_files <- paste0(raw_data_dir, "/", grep(pattern = "sgRNA-enrichment_5K-sgRNAs_Batch", x = raw_fs, value = TRUE))

# get the batch ID of each file
batch <- stringr::str_extract(string = gRNA_files, pattern = "Batch_[0-9]_[0-9]") %>%
  gsub(pattern = "Batch_", replacement = "")

# construct a sparse matrix of gRNA counts for each file
res <- lapply(X = seq(1L, length(gRNA_files)), FUN = function(i) {
  curr_file <- gRNA_files[i]
  curr_batch <- batch[i]
  print(paste("Working on file", curr_file))
  curr_gRNA_count_matrix <- readr::read_tsv(file = curr_file,
                                            col_names = c("cell_barcode", "total_read_count", "total_umi_count", "gRNA_spacer_seqs", "read_counts", "umi_counts"),
                                            col_types = c("cccccc")) %>% dplyr::select(cell_barcode, gRNA_spacer_seqs, umi_counts, total_umi_count)
  cell_barcodes <- curr_gRNA_count_matrix$cell_barcode
  # convert to sparse triplet matrix
  umi_count_list <- lapply(curr_gRNA_count_matrix$umi_counts, function(v) stringr::str_split(v, pattern = ";") %>% unlist() %>% as.integer())
  gRNA_spacer_seq_list <- lapply(curr_gRNA_count_matrix$gRNA_spacer_seqs, function(v) stringr::str_split(v, pattern = ";") %>% unlist())
  n_entries_per_cell <- sapply(X = umi_count_list, function(v) length(v))
  barcode_vector <- rep(x = cell_barcodes, times = n_entries_per_cell)
  gRNA_spacer_seq_vector <- unlist(gRNA_spacer_seq_list)
  umi_count_vector <- unlist(umi_count_list)
  ch_df <- data.frame(barcode = barcode_vector,
                      spacer_seq = gRNA_spacer_seq_vector,
                      umi_count = umi_count_vector)
  # remove all entries with spacer seqs not in the spacer seq list
  ch_df <- ch_df %>% dplyr::filter(gRNA_spacer_seq_vector %in% guide_seqs$spacer_seq)
  ordered_barcode_vector <- sort(unique(ch_df$barcode))
  
  spacer_seq_idx <- match(x = ch_df$spacer_seq, guide_seqs$spacer_seq)
  barcode_idx <- match(x = ch_df$barcode, table = ordered_barcode_vector)
  
  m <- Matrix::sparseMatrix(i = spacer_seq_idx,
                            j = barcode_idx,
                            x = ch_df$umi_count,
                            dims = c(length(guide_seqs$spacer_seq),
                                     length(ordered_barcode_vector)))
  rownames(m) <- guide_seqs$spacer_seq
  colnames(m) <- paste0(ordered_barcode_vector, "_", curr_batch)
  return(m)
})

combined_m <- do.call(what = cbind, args = res)
all(Matrix::colSums(combined_m) >= 1); all(Matrix::rowSums(combined_m) >= 1)

# save in intermediate file directory
saveRDS(object = combined_m, file = paste0(intermediate_file_dir, "gRNA_matrix.rds"))
