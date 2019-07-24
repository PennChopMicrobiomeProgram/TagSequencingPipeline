#!/bin/bash
set -x
set -e
#set -u

#source activate qiime2-2018.11

if [ $# -ne 2 ]; then
    echo "Usage: $0 MAPPING_FP SAMPLING_DEPTH"
    echo "MAPPING_FP is the mapping file for Qimme2"
    echo "SAMPLING_DEPTH is the number of features to rarefy to"
    echo "Required the original qiime2_pipeline.bash to be run first"
    exit 1
fi

MAPPING_FP=$1
SAMPLING_DEPTH=$2
WORK_DIR="$(dirname ${MAPPING_FP})"

DENOISE_DIR="${WORK_DIR}/denoising_results"
RARE_DIR="${WORK_DIR}/rarefied_seqs"

###=====================
### RAREFY TABLE
###=====================

if [ ! -e "${RARE_DIR}/rarefied_table.qza" ]; then
    qiime feature-table rarefy \
        --i-table "${DENOISE_DIR}/table.qza" \
        --p-sampling-depth $SAMPLING_DEPTH \
        --output-dir ${RARE_DIR}
fi

if [ -e "${RARE_DIR}/rarefied_table.qza" ]; then
    qiime feature-table summarize \
      --i-table "${RARE_DIR}/rarefied_table.qza" \
      --o-visualization "${RARE_DIR}/table.qzv" \
      --m-sample-metadata-file ${MAPPING_FP}
fi

if [ -e "${RARE_DIR}/rarefied_table.qza" ]; then
    qiime tools export \
      --input-path "${RARE_DIR}/table.qzv" \
      --output-path "${RARE_DIR}/table"
fi

if [ -e "${RARE_DIR}/rarefied_table.qza" ]; then
    qiime tools export \
      --input-path "${RARE_DIR}/rarefied_table.qza" \
      --output-path "${RARE_DIR}/table"
fi


###=====================
###  ALPHA AND BETA DIVERSITY
###=====================

if [ ! -e "${RARE_DIR}/faith_pd_vector.qza" ]; then
    qiime diversity alpha-phylogenetic \
      --i-phylogeny "${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${RARE_DIR}/rarefied_table.qza" \
      --p-metric faith_pd \
      --o-alpha-diversity "${RARE_DIR}/faith_pd_vector.qza"
fi

qiime tools export \
  --input-path "${RARE_DIR}/faith_pd_vector.qza" \
  --output-path "${RARE_DIR}/faith"

if [ ! -e "${RARE_DIR}/weighted_unifrac_distance_matrix.qza" ]; then
    qiime diversity beta-phylogenetic \
      --i-phylogeny "${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${RARE_DIR}/rarefied_table.qza" \
      --p-metric weighted_unifrac \
      --o-distance-matrix "${RARE_DIR}/weighted_unifrac_distance_matrix.qza"
fi

qiime tools export \
  --input-path "${RARE_DIR}/weighted_unifrac_distance_matrix.qza" \
  --output-path "${RARE_DIR}/wu"

if [ ! -e "${RARE_DIR}/unweighted_unifrac_distance_matrix.qza" ]; then
    qiime diversity beta-phylogenetic \
      --i-phylogeny "${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${RARE_DIR}/rarefied_table.qza" \
      --p-metric unweighted_unifrac \
      --o-distance-matrix "${RARE_DIR}/unweighted_unifrac_distance_matrix.qza"
fi

qiime tools export \
  --input-path "${RARE_DIR}/unweighted_unifrac_distance_matrix.qza" \
  --output-path "${RARE_DIR}/uu"

###=====================
###  BIOM CONVERT
###=====================

if [ ! -e "${RARE_DIR}/table/feature-table.tsv" ]; then
    biom convert \
      -i "${RARE_DIR}/table/feature-table.biom" \
      -o "${RARE_DIR}/table/feature-table.tsv" \
      --to-tsv
fi
