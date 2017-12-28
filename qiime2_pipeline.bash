#!/bin/bash

if [ $# -ne 3 ]; then
	echo "Usage: $0 WORK_DIR MAPPING_FP PROJECT_NAME"
	exit 1
fi

WORK_DIR=$1
MAPPING_FP=$2
PROJECT_NAME=$3

## Taxonomy classifier setup. Two classifiers are currently available:
## classifiers trained on full length and on 515F/806R region of Greengenes 13_8 99% OTUs
## These can be downloaded from https://data.qiime2.org/2017.9/common/gg-13-8-99-nb-classifier.qza (full length)
## or https://data.qiime2.org/2017.9/common/gg-13-8-99-515-806-nb-classifier.qza (515F/806R region)

CLASSIFIER_FP="${HOME}/gg-13-8-99-nb-classifier.qza"
#CLASSIFIER_FP="${HOME}/gg-13-8-99-515-806-nb-classifier.qza" ## used for V4 region

### PATH TO PYTHON2
PHYTON2="/home/leej39/miniconda3/bin/python"

### PATH TO Ceylan's CODE TO COMBINE I1 and I2
INDEX1_INDEX2_COMBINE_SCRIPT="${HOME}/combine_barcodes.py"

DATA_DIR="${WORK_DIR}/data_files"

PROJECT_DIR="${WORK_DIR}/${PROJECT_NAME}"
if [ ! -d ${PROJECT_DIR} ]; then
	mkdir ${PROJECT_DIR}
fi

###=====================
### gunzip INDEX1 AND INDEX2, IF NECESSARY
###=====================

if [ -e "${DATA_DIR}/Undetermined_S0_L001_I1_001.fastq.gz" ]; then
	gunzip "${DATA_DIR}/Undetermined_S0_L001_I1_001.fastq.gz"
fi

if [ -e "${DATA_DIR}/Undetermined_S0_L001_I2_001.fastq.gz" ]; then
        gunzip "${DATA_DIR}/Undetermined_S0_L001_I2_001.fastq.gz"
fi

###=====================
### gzip R1 AND R2, IF NECESSARY
###=====================

if [ -e "${DATA_DIR}/Undetermined_S0_L001_R1_001.fastq" ]; then
        gzip "${DATA_DIR}/Undetermined_S0_L001_R1_001.fastq"
fi

if [ -e "${DATA_DIR}/Undetermined_S0_L001_R2_001.fastq" ]; then
        gzip "${DATA_DIR}/Undetermined_S0_L001_R2_001.fastq"
fi

###=====================
### COMBINE INDEX1 AND INDEX2 AND gzip
###=====================

${PHYTON2} ${INDEX1_INDEX2_COMBINE_SCRIPT} ${DATA_DIR}
gzip "${DATA_DIR}/Undetermined_S0_L001_I12_001.fastq"

FWD="${DATA_DIR}/Undetermined_S0_L001_R1_001.fastq.gz"
REV="${DATA_DIR}/Undetermined_S0_L001_R2_001.fastq.gz"
IDX="${DATA_DIR}/Undetermined_S0_L001_I12_001.fastq.gz"

###=====================
### DATA IMPORT
###=====================

EMP_PAIRED_END_SEQUENCES_DIR="${PROJECT_DIR}/emp-paired-end-sequences"
if [ ! -d ${EMP_PAIRED_END_SEQUENCES_DIR} ]; then
        mkdir ${EMP_PAIRED_END_SEQUENCES_DIR}
fi

mv ${FWD} "${EMP_PAIRED_END_SEQUENCES_DIR}/forward.fastq.gz"
mv ${REV} "${EMP_PAIRED_END_SEQUENCES_DIR}/reverse.fastq.gz"
mv ${IDX} "${EMP_PAIRED_END_SEQUENCES_DIR}/barcodes.fastq.gz"

qiime tools import \
  --type EMPPairedEndSequences \
  --input-path ${EMP_PAIRED_END_SEQUENCES_DIR} \
  --output-path "${PROJECT_DIR}/emp-paired-end-sequences.qza"

###=====================
### DEMULTIPLEXING SEQUENCE
###=====================

qiime demux emp-paired \
  --m-barcodes-file ${MAPPING_FP} \
  --m-barcodes-category BarcodeSequence \
  --i-seqs "${PROJECT_DIR}/emp-paired-end-sequences.qza" \
  ## --p-rev-comp-mapping-barcodes \ ## use if mapping file BarcodeSequence nees to be reverse complemented
  --o-per-sample-sequences "${PROJECT_DIR}/demux.qza"

qiime demux summarize \
  --i-data "${PROJECT_DIR}/demux.qza" \
  --o-visualization "${PROJECT_DIR}/demux.qzv"

qiime tools export \
  "${PROJECT_DIR}/demux.qzv" \
  --output-dir "${PROJECT_DIR}/demux"

###=====================
###  SEQUENCE QC AND FEATURE TABLE
###=====================

## discussion needed for denosing parameters below

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${PROJECT_DIR}/demux.qza" \
  --p-trim-left-f 0 \
  --p-trunc-len-f 240 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 240 \
  --p-n-threads 8 \
  --o-representative-sequences "${PROJECT_DIR}/rep-seqs.qza" \
  --o-table "${PROJECT_DIR}/table.qza"

qiime feature-table summarize \
  --i-table "${PROJECT_DIR}/table.qza" \
  --o-visualization "${PROJECT_DIR}/table.qzv" \
  --m-sample-metadata-file ${MAPPING_FP}

qiime feature-table tabulate-seqs \
  --i-data "${PROJECT_DIR}/rep-seqs.qza" \
  --o-visualization "${PROJECT_DIR}/rep-seqs.qzv"

qiime tools export \
  "${PROJECT_DIR}/table.qza" \
  --output-dir "${PROJECT_DIR}/table"

###=====================
###  TAXONOMIC ANALYSIS
###=====================

qiime feature-classifier classify-sklearn \
  --i-classifier ${CLASSIFIER_FP} \
  --i-reads "${PROJECT_DIR}/rep-seqs.qza" \
  --o-classification "${PROJECT_DIR}/taxonomy.qza"

qiime metadata tabulate \
  --m-input-file "${PROJECT_DIR}/taxonomy.qza" \
  --o-visualization "${PROJECT_DIR}/taxonomy.qzv"

qiime tools export \
  "${PROJECT_DIR}/taxonomy.qza" \
  --output-dir "${PROJECT_DIR}/taxonomy"

###=====================
###  GENERATE TREES
###=====================

qiime alignment mafft \
  --i-sequences "${PROJECT_DIR}/rep-seqs.qza" \
  --o-alignment "${PROJECT_DIR}/aligned-rep-seqs.qza"

qiime alignment mask \
  --i-alignment "${PROJECT_DIR}/aligned-rep-seqs.qza" \
  --o-masked-alignment "${PROJECT_DIR}/masked-aligned-rep-seqs.qza"

qiime phylogeny fasttree \
  --i-alignment "${PROJECT_DIR}/masked-aligned-rep-seqs.qza" \
  --o-tree "${PROJECT_DIR}/unrooted-tree.qza"

qiime phylogeny midpoint-root \
  --i-tree "${PROJECT_DIR}/unrooted-tree.qza" \
  --o-rooted-tree "${PROJECT_DIR}/rooted-tree.qza"

###=====================
###  ALPHA AND BETA DIVERSITY
###=====================

### method "core-metrics-phylogenetic" was added in qiime2-2017.10

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "${PROJECT_DIR}/rooted-tree.qza" \
  --i-table "${PROJECT_DIR}/table.qza" \
  --p-sampling-depth 10000 \
  --m-metadata-file ${MAPPING_FP} \
  --output-dir "${PROJECT_DIR}/core-metrics-results"

qiime tools export \
  "${PROJECT_DIR}/core-metrics-results/faith_pd_vector.qza" \
  --output-dir "${PROJECT_DIR}/core-metrics-results/faith"

qiime tools export \
  "${PROJECT_DIR}/core-metrics-results/unweighted_unifrac_distance_matrix.qza" \
  --output-dir "${PROJECT_DIR}/core-metrics-results/uu"

qiime tools export \
  "${PROJECT_DIR}/core-metrics-results/weighted_unifrac_distance_matrix.qza" \
  --output-dir "${PROJECT_DIR}/core-metrics-results/wu"



###=====================
###  BIOM CONVERT
###=====================

biom convert \
  -i "${PROJECT_DIR}/table/feature-table.biom" \
  -o "${PROJECT_DIR}/table/feature-table.tsv" \
  --to-tsv
