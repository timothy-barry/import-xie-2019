#########################################
# 0. Load data, packages, set directories
#########################################
library(magrittr)
xie_offsite <-.get_config_path("LOCAL_XIE_2019_DATA_DIR")
raw_data_dir <- paste0(xie_offsite, "raw/")
intermediate_file_dir <- paste0(xie_offsite, "intermediate/")
processed_data_dir <- paste0(xie_offsite, "processed/")
if (!dir.exists(processed_data_dir)) dir.create(processed_data_dir)
