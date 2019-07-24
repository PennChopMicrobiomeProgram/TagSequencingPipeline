#!/bin/bash
set -x
set -e
#set -u

if [ $# -ne 3 ]; then
    echo "Usage: $0 DIRECTORY1 DIRECTORY2 DIRECTORY_OUT"
    exit 1
fi

DIR1=$1
DIR2=$2
DIR_OUT=$3

DENOISE_DIR="denoising_results"
METRIC_DIR="${DIR_OUT}/core_metrics_results"

mkdir -p "${DIR_OUT}"
mkdir -p "${DIR_OUT}/${DENOISE_DIR}"


## Taxonomy classifier setup. Two classifiers are currently available:
## classifiers trained on full length and on 515F/806R region of Greengenes 13_8 99% OTUs
## These can be downloaded from https://data.qiime2.org/2017.9/common/gg-13-8-99-nb-classifier.qza (full length)
## or https://data.qiime2.org/2017.9/common/gg-13-8-99-515-806-nb-classifier.qza (515F/806R region)

CLASSIFIER_FP="gg-13-8-99-nb-classifier.qza"
#CLASSIFIER_FP="gg-13-8-99-515-806-nb-classifier.qza" ## used for V4 region
#CLASSIFIER_FP="gg-13-8-99-27-338-nb-classifier.qza" ## trained for V1V2 region truncated at 350 bp



###===============
###  MERGE FILES 
###===============

#NOtE: if you have duplicate sample-ids you can use a trick here to rename them without running an entire pipeline again:
#https://forum.qiime2.org/t/change-sample-ids-after-running-dada2/3918/3

qiime feature-table merge \
      --i-tables "${DIR1}/${DENOISE_DIR}/table.qza" "${DIR2}/${DENOISE_DIR}/table.qza" \
      --o-merged-table "${DIR_OUT}/${DENOISE_DIR}/table.qza"


qiime tools export \
      --input-path "${DIR_OUT}/${DENOISE_DIR}/table.qza" \
      --output-path "${DIR_OUT}/${DENOISE_DIR}/table"


qiime feature-table merge-seqs \
      --i-data "${DIR1}/${DENOISE_DIR}/rep-seqs.qza" "${DIR2}/${DENOISE_DIR}/rep-seqs.qza" \
      --o-merged-data "${DIR_OUT}/${DENOISE_DIR}/rep-seqs.qza"


###=====================
###  TAXONOMIC ANALYSIS
###=====================

qiime feature-classifier classify-sklearn \
      --i-classifier "${CLASSIFIER_FP}" \
      --i-reads "${DIR_OUT}/${DENOISE_DIR}/rep-seqs.qza" \
      --o-classification "${DIR_OUT}/${DENOISE_DIR}/taxonomy.qza"

qiime metadata tabulate \
      --m-input-file "${DIR_OUT}/${DENOISE_DIR}/taxonomy.qza" \
      --o-visualization "${DIR_OUT}/${DENOISE_DIR}/taxonomy.qzv"

qiime tools export \
      --input-path "${DIR_OUT}/${DENOISE_DIR}/taxonomy.qza" \
      --output-path "${DIR_OUT}/${DENOISE_DIR}/taxonomy"



###=====================
###  GENERATE TREES
###=====================

qiime alignment mafft \
      --i-sequences "${DIR_OUT}/${DENOISE_DIR}/rep-seqs.qza" \
      --o-alignment "${DIR_OUT}/${DENOISE_DIR}/aligned-rep-seqs.qza"

qiime alignment mask \
      --i-alignment "${DIR_OUT}/${DENOISE_DIR}/aligned-rep-seqs.qza" \
      --o-masked-alignment "${DIR_OUT}/${DENOISE_DIR}/masked-aligned-rep-seqs.qza"

qiime phylogeny fasttree \
      --i-alignment "${DIR_OUT}/${DENOISE_DIR}/masked-aligned-rep-seqs.qza" \
      --o-tree "${DIR_OUT}/${DENOISE_DIR}/unrooted-tree.qza"

qiime phylogeny midpoint-root \
      --i-tree "${DIR_OUT}/${DENOISE_DIR}/unrooted-tree.qza" \
        --o-rooted-tree "${DIR_OUT}/${DENOISE_DIR}/rooted-tree.qza"



###=====================
###  ALPHA AND BETA DIVERSITY
###=====================

if [ ! -d ${METRIC_DIR} ]; then
    mkdir -p ${METRIC_DIR}
fi

qiime diversity alpha-phylogenetic \
      --i-phylogeny "${DIR_OUT}/${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${DIR_OUT}/${DENOISE_DIR}/table.qza" \
      --p-metric faith_pd \
      --o-alpha-diversity "${METRIC_DIR}/faith_pd_vector.qza"

qiime tools export \
      --input-path "${METRIC_DIR}/faith_pd_vector.qza" \
      --output-path "${METRIC_DIR}/faith"

qiime diversity beta-phylogenetic \
      --i-phylogeny "${DIR_OUT}/${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${DIR_OUT}/${DENOISE_DIR}/table.qza" \
      --p-metric weighted_unifrac \
      --o-distance-matrix "${METRIC_DIR}/weighted_unifrac_distance_matrix.qza"

qiime tools export \
      --input-path "${METRIC_DIR}/weighted_unifrac_distance_matrix.qza" \
      --output-path "${METRIC_DIR}/wu"

qiime diversity beta-phylogenetic \
      --i-phylogeny "${DIR_OUT}/${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${DIR_OUT}/${DENOISE_DIR}/table.qza" \
      --p-metric unweighted_unifrac \
      --o-distance-matrix "${METRIC_DIR}/unweighted_unifrac_distance_matrix.qza"

qiime tools export \
      --input-path "${METRIC_DIR}/unweighted_unifrac_distance_matrix.qza" \
      --output-path "${METRIC_DIR}/uu"




###=====================
###  BIOM CONVERT
###=====================

biom convert \
     -i "${DIR_OUT}/${DENOISE_DIR}/table/feature-table.biom" \
     -o "${DIR_OUT}/${DENOISE_DIR}/table/feature-table.tsv" \
     --to-tsv
