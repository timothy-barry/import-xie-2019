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
guide_seqs <- openxlsx::read.xlsx(xlsxFile = paste0(raw_data_dir, "/all_oligos.xlsx"),
                                sheet = 1) %>% dplyr::rename(hg38_enh_region = "region.pos.(hg38)")
target_regions <- guide_seqs %>% dplyr::pull(hg38_enh_region) %>% unique() # 518 target regions; all but two targeted by 10 gRNAs.
names(target_regions) <- target_regions
raw_fs <- list.files(raw_data_dir)
gRNA_files <- paste0(raw_data_dir, "/", grep(pattern = "sgRNA-enrichment_5K-sgRNAs_Batch", x = raw_fs, value = TRUE))

# get the batch ID of each file
batch <- stringr::str_extract(string = gRNA_files, pattern = "Batch_[0-9]_[0-9]") %>%
  gsub(pattern = "Batch_", replacement = "")

# get sparse matrix of ungrouped counts
future::plan(future::multisession())
res <- furrr::future_map(.x = seq(1L, length(gRNA_files)), .f = function(i) {
  curr_file <- gRNA_files[i]
  curr_batch <- batch[i]
  print(paste("Working on file", curr_file))
  curr_gRNA_count_matrix <- readr::read_tsv(file = curr_file,
                                            col_names = c("cell_barcode", "total_read_count", "total_umi_count", "gRNA_spacer_seqs", "read_counts", "umi_counts"),
                                            col_types = c("cccccc")) %>% dplyr::select(cell_barcode, gRNA_spacer_seqs, umi_counts, total_umi_count)
  cell_barcodes <- dplyr::pull(curr_gRNA_count_matrix, cell_barcode)
  # check for duplicates
  if (any(duplicated(cell_barcodes))) stop("Duplicated barcodes")
  barcodes_annotated <- paste0(cell_barcodes, "-", curr_batch)
  # create the count matrix list
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
  list(cell_barcodes = barcodes_annotated, count_matrix_list = count_matrix_list, total_umis = as.integer(curr_gRNA_count_matrix$total_umi_count))
})

# save in intermediate file directory
saveRDS(object = res, file = paste0(intermediate_file_dir, "gRNA_matrix_list.rds"))
