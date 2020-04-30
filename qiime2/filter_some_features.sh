#!/usr/bin/env bash

#$ -cwd
#$ -r n
#$ -j y
#$ -m ea
#$ -M danielsg@email.chop.edu
#$ -l h_vmem=50G
#$ -l m_mem_free=20G

set -x
set -e
#set -u

source ~/.bashrc
conda activate qiime2-2019.4

if [ $# -ne 1 ]; then
	echo "Usage: $0 MAPPING_FP"
    echo "MAPPING_FP is the mapping file for Qiime2"
    echo "Data files should be in a directory named "data_files""
	exit 1
fi

MAPPING_FP=$1
WORK_DIR="/mnt/isilon/microbiome/analysis/danielsg/sehgal_run_6"

# See https://www.ostricher.com/2014/10/the-right-way-to-get-the-directory-of-a-bash-script/
SOURCE_REL="${BASH_SOURCE[0]}"
SOURCE_ABS=$(readlink -f "${SOURCE_REL}")
SOURCE_DIR=$(dirname "${SOURCE_ABS}")

### PATH TO Ceylan's CODE TO COMBINE I1 and I2
INDEX1_INDEX2_COMBINE_SCRIPT="${SOURCE_DIR}/combine_barcodes.py"

## Taxonomy classifier setup. Two classifiers are currently available:
## classifiers trained on full length and on 515F/806R region of Greengenes 13_8 99% OTUs
## These can be downloaded from https://data.qiime2.org/2017.9/common/gg-13-8-99-nb-classifier.qza (full length)
## or https://data.qiime2.org/2017.9/common/gg-13-8-99-515-806-nb-classifier.qza (515F/806R region)

#CLASSIFIER_FP="${HOME}/gg-13-8-99-nb-classifier.qza"
#CLASSIFIER_FP="gg-13-8-99-515-806-nb-classifier.qza" ## used for V4 region
CLASSIFIER_FP="${WORK_DIR}/gg-13-8-99-27-338-nb-classifier.qza" ## trained for V1V2 region truncated at 350 bp

DATA_DIR="${WORK_DIR}/Data"

#EMP_PAIRED_END_SEQUENCES_DIR="${DATA_DIR}/emp_paired_end_sequences"
#DATA_DIR="${DATA_DIR}/data_files"
DEMUX_DIR="${DATA_DIR}/demux_results"
DENOISE_DIR="${DATA_DIR}/denoising_results"
METRIC_DIR="${DATA_DIR}/core_metrics_results"

qiime taxa filter-table \
  --i-table "${DENOISE_DIR}/table.qza" \
  --i-taxonomy "${DENOISE_DIR}/taxonomy.qza" \
  --p-exclude "Corynebacterium,Propionibacterium,Streptococcus,Acidovorax,Bradyrhizobium,Comamonas,Pelomonas,Pseudomonas,Ralstonia" \
  --o-filtered-table "${DENOISE_DIR}/table-no-contaminants.qza"
#
#qiime feature-table summarize \
#      --i-table "${DENOISE_DIR}/table-no-contaminants.qza" \
#      --o-visualization "${DENOISE_DIR}/table-no-contaminants.qzv" \
#      --m-sample-metadata-file subset_metadata.tsv
#
qiime tools export \
      --input-path "${DENOISE_DIR}/table-no-contaminants.qza" \
      --output-path "${DENOISE_DIR}/table-no-contaminants"

qiime diversity alpha-phylogenetic \
      --i-phylogeny "${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${DENOISE_DIR}/table-no-contaminants.qza" \
      --p-metric faith_pd \
      --o-alpha-diversity "${METRIC_DIR}/faith_pd_vector_no_contam.qza"

qiime tools export \
      --input-path "${METRIC_DIR}/faith_pd_vector_no_contam.qza" \
      --output-path "${METRIC_DIR}/faith_no_contam"

qiime diversity beta-phylogenetic \
      --i-phylogeny "${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${DENOISE_DIR}/table-no-contaminants.qza" \
      --p-metric weighted_unifrac \
      --o-distance-matrix "${METRIC_DIR}/weighted_unifrac_distance_matrix_no_contam.qza"

qiime tools export \
  --input-path "${METRIC_DIR}/weighted_unifrac_distance_matrix_no_contam.qza" \
  --output-path "${METRIC_DIR}/wu_no_contam"

qiime diversity beta-phylogenetic \
      --i-phylogeny "${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${DENOISE_DIR}/table-no-contaminants.qza" \
      --p-metric unweighted_unifrac \
      --o-distance-matrix "${METRIC_DIR}/unweighted_unifrac_distance_matrix_no_contam.qza"

qiime tools export \
  --input-path "${METRIC_DIR}/unweighted_unifrac_distance_matrix_no_contam.qza" \
  --output-path "${METRIC_DIR}/uu_no_contam"

###=====================
###  BIOM CONVERT
###=====================
if [ ! -e "${DENOISE_DIR}/table-no-contaminants/feature-table.tsv" ]; then
    biom convert \
      -i "${DENOISE_DIR}/table-no-contaminants/feature-table.biom" \
      -o "${DENOISE_DIR}/table-no-contaminants/feature-table.tsv" \
      --to-tsv
fi

cd "${DENOISE_DIR}"

cp taxonomy/taxonomy.tsv biom-taxonomy.tsv

# had to vim biom-taxonomy.tsv and add 
# "#OTUID	taxonomy	confidence" (without quotes)
# line to make it recognize the header and correctly add the taxonomy
# see here: https://forum.qiime2.org/t/exporting-and-modifying-biom-tables-e-g-adding-taxonomy-annotations/3630

biom add-metadata \
    -i table-no-contaminants/feature-table.biom \
    -o table-with-taxonomy.biom \
    --observation-metadata-fp biom-taxonomy.tsv \
    --sc-separated taxonomy

biom convert \
    -i table-with-taxonomy.biom \
    -o table-with-taxa.tsv \
    --to-tsv --header-key taxonomy

echo "****"
echo "Done!"
echo "If nothing happened,"
echo "You may already have result files, check your directories!"
