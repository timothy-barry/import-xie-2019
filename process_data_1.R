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
guide_seqs <- openxlsx::read.xlsx(xlsxFile = paste0(raw_data_dir, "/all_oligos.xlsx"),
                                sheet = 1) %>% dplyr::rename(hg38_enh_region = "region.pos.(hg38)")
target_regions <- guide_seqs %>% dplyr::pull(hg38_enh_region) %>% unique() # 518 target regions; all but two targeted by 10 gRNAs.
names(target_regions) <- target_regions

raw_fs <- list.files(raw_data_dir)
gRNA_files <- paste0(raw_data_dir, "/", grep(pattern = "sgRNA-enrichment_5K-sgRNAs_Batch", x = raw_fs, value = TRUE))

# get sparse matrix of ungrouped counts
future::plan(future::multisession())
res <- furrr::future_map(.x = gRNA_files, .f = function(curr_file) {
  print(paste("Working on file", curr_file))
  curr_gRNA_count_matrix <- readr::read_tsv(file = curr_file, 
                                            col_names = c("cell_barcode", "total_read_count", "total_umi_count", "gRNA_spacer_seqs", "read_counts", "umi_counts"),
                                            col_types = c("cccccc")) %>% dplyr::select(cell_barcode, gRNA_spacer_seqs, umi_counts, total_umi_count)
  cell_barcodes <- dplyr::pull(curr_gRNA_count_matrix, cell_barcode)
  
  count_matrix_list <- purrr:::map(target_regions, function(region) {
    cat(paste0("Working on region ", region, ".\n"))
    region_spacer_seqs <- dplyr::filter(guide_seqs, hg38_enh_region == region) %>% dplyr::pull(spacer.sequence)
    curr_batch_gRNA_umi_counts <- sapply(X = 1:nrow(curr_gRNA_count_matrix), FUN = function(row_id) {
      r <- curr_gRNA_count_matrix[row_id,]
      spacers <- stringr::str_split(r$gRNA_spacer_seqs, pattern = ";") %>% unlist()
      umi_counts <- stringr::str_split(r$umi_counts, pattern = ";") %>% unlist() %>% as.integer()
      umi_locs <- match(x = region_spacer_seqs, table = spacers)
      curr_counts <- sapply(umi_locs, function(i) if (is.na(i)) 0 else umi_counts[i])
      names(curr_counts) <- region_spacer_seqs
      curr_counts
    }) %>% t() %>% Matrix::Matrix(sparse = TRUE)
  })
  list(cell_barcodes = cell_barcodes, count_matrix_list = count_matrix_list, total_umis = as.integer(curr_gRNA_count_matrix$total_umi_count))
})

# save in intermediate file directory
saveRDS(object = res, file = paste0(intermediate_file_dir, "gRNA_matrix_list.rds"))

# gRNA_count_matrix_list <- purrr::map(target_regions, function(region) {
#  map(res, function(x) x$count_matrix_list[[region]]) %>% purrr::reduce(.f = rbind)
# })
# cell_barcodes <- purrr::map(.x = res, .f = function(x) x$cell_barcodes) %>% purrr::reduce(.f = c)
# cell_gRNA_umi_counts <- purrr::map(.x = res, function(x) x$total_umis) %>% purrr::reduce(.f = c)
# We reduce each matrix in the gRNA_count_matrix_list to a single logical vector
# combine_gRNAs_in_group <- function(gRNA_count_matrix) {
#  apply(X = gRNA_count_matrix, MARGIN = 2, FUN = function(column) {
#    v <- sum(column >= 1)
#    U <- sum(column)
#    column/U > 1/v
#  }) %>% apply(MARGIN = 1, FUN = function(r) any(r))
# }
# gRNA_indic_matrix <- map_dfr(gRNA_count_matrix_list, combine_gRNAs_in_group)
# Finally, confirm that the cell barcode order for the gRNA indicator matrix matches that of the cell-by-gene expression matrix and cell-specific covariate matrix. Also, append the gRNA UMI count to the cell covariate matrix.
# cell_covariate_matrix <- read.fst(paste0(processed_dir, "/cell_covariate_matrix.fst"))
# cell_barcodes_to_check <- pull(cell_covariate_matrix, ordered_cell_barcodes) %>% gsub(pattern = "-1", replacement = "")
# m <- match(x = cell_barcodes_to_check, table = cell_barcodes) # There will be some na's.
# gRNA_indic_matrix_ordered <- gRNA_indic_matrix[m,]
# cell_gRNA_umi_counts <- cell_gRNA_umi_counts[m]