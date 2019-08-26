#fixing metadata EXAMPLE

library(tidyverse)

root_dir = "/Volumes/microbiome/analysis/danielsg/anguera_mont"

zwei_mp_fp <- file.path(root_dir, "anguera_run2_20190507", "anguera_run2_mapping_file.tsv")

metadata_causing_problems <- read_tsv(zwei_mp_fp)

new_metadata <- metadata_causing_problems %>%
  mutate(new_sample_id = paste0(`#SampleID`,"_run2"))

write_tsv(x = new_metadata, path = file.path(root_dir, "anguera_run2_20190507", "anguera_run2_mapping_file_pre.tsv"))

#then, in bash

###bash###
# conda activate qiime2-2019.4
# 
# qiime feature-table group --i-table denoising_results/table.qza --m-metadata-file anguera_run2_mapping_file_for_renaming_trick.tsv --m-metadata-column "new_sample_id" --p-mode "sum" --o-grouped-table denoising_results/new_table.qza --p-axis "sample"
# 
# cp anguera_run2_mapping_file_pre.tsv anguera_run2_mapping_file_post.tsv
# 
# cp denoising_results/table.qza denoising_results/old_table.qza
# 
# mv denoising_results/new_table.qza denoising_results/table.qza

new_new_metadata <- new_metadata %>%
  select(-`#SampleID`) %>%
  select(`#SampleID` = new_sample_id, everything())

write_tsv(x = new_new_metadata, path = file.path(root_dir, "anguera_run2_20190507", "anguera_run2_mapping_file_post.tsv"))

#and then back to qiime2_merge_runs.sh
