#!/bin/bash
set -x
set -e
set -u

if [ $# -ne 2 ]; then
	echo "Usage: $0 WORK_DIR MAPPING_FP"
	exit 1
fi

WORK_DIR=$1
MAPPING_FP=$2

## Taxonomy classifier setup. Two classifiers are currently available:
## classifiers trained on full length and on 515F/806R region of Greengenes 13_8 99% OTUs
## These can be downloaded from https://data.qiime2.org/2017.9/common/gg-13-8-99-nb-classifier.qza (full length)
## or https://data.qiime2.org/2017.9/common/gg-13-8-99-515-806-nb-classifier.qza (515F/806R region)

#CLASSIFIER_FP="${HOME}/gg-13-8-99-nb-classifier.qza"
#CLASSIFIER_FP="${HOME}/gg-13-8-99-515-806-nb-classifier.qza" ## used for V4 region
CLASSIFIER_FP="gg-13-8-99-27-338-nb-classifier.qza" ## trained for V1V2 region truncated at 350 bp


### PATH TO Ceylan's CODE TO COMBINE I1 and I2

INDEX1_INDEX2_COMBINE_SCRIPT="${HOME}/combine_barcodes.py"

EMP_PAIRED_END_SEQUENCES_DIR="${WORK_DIR}/emp-paired-end-sequences"
DATA_DIR="${WORK_DIR}/data_files"
DEMUX_DIR="${WORK_DIR}/demux-results"
DENOISE_DIR="${WORK_DIR}/denoising-results"
METRIC_DIR="${WORK_DIR}/core-metric-results"

if [ ! -d ${EMP_PAIRED_END_SEQUENCES_DIR} ]; then
        mkdir ${EMP_PAIRED_END_SEQUENCES_DIR}
fi

if [ ! -d ${DEMUX_DIR} ]; then
        mkdir ${DEMUX_DIR}
fi

if [ ! -d ${DENOISE_DIR} ]; then
        mkdir ${DENOISE_DIR}
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

python ${INDEX1_INDEX2_COMBINE_SCRIPT} ${DATA_DIR}
gzip "${DATA_DIR}/Undetermined_S0_L001_I12_001.fastq"

FWD="${DATA_DIR}/Undetermined_S0_L001_R1_001.fastq.gz"
REV="${DATA_DIR}/Undetermined_S0_L001_R2_001.fastq.gz"
IDX="${DATA_DIR}/Undetermined_S0_L001_I12_001.fastq.gz"

###=====================
### DATA IMPORT
###=====================


if [ ! -e "${EMP_PAIRED_END_SEQUENCES_DIR}/forward.fastq.gz" ]; then
    mv ${FWD} "${EMP_PAIRED_END_SEQUENCES_DIR}/forward.fastq.gz"
    mv ${REV} "${EMP_PAIRED_END_SEQUENCES_DIR}/reverse.fastq.gz"
    mv ${IDX} "${EMP_PAIRED_END_SEQUENCES_DIR}/barcodes.fastq.gz"
fi

qiime tools import \
  --type EMPPairedEndSequences \
  --input-path ${EMP_PAIRED_END_SEQUENCES_DIR} \
  --output-path "${WORK_DIR}/emp-paired-end-sequences.qza"

###=====================
### DEMULTIPLEXING SEQUENCE
###=====================

qiime demux emp-paired \
  --m-barcodes-file ${MAPPING_FP} \
  --m-barcodes-category BarcodeSequence \
  --i-seqs "${WORK_DIR}/emp-paired-end-sequences.qza" \
  --p-rev-comp-mapping-barcodes \
  --o-per-sample-sequences "${DEMUX_DIR}/demux.qza"

qiime demux summarize \
  --i-data "${DEMUX_DIR}/demux.qza" \
  --o-visualization "${DEMUX_DIR}/demux.qzv"

qiime tools export \
  "${DEMUX_DIR}/demux.qzv" \
  --output-dir "${DEMUX_DIR}/demux-qzv"

###=====================
###  SEQUENCE QC AND FEATURE TABLE
###=====================

## discussion needed for denosing parameters below

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${DEMUX_DIR}/demux.qza" \
  --p-trim-left-f 0 \
  --p-trunc-len-f 230 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 230 \
  --p-n-threads 8 \
  --o-representative-sequences "${DENOISE_DIR}/rep-seqs.qza" \
  --o-table "${DENOISE_DIR}/table.qza"

qiime feature-table summarize \
  --i-table "${DENOISE_DIR}/table.qza" \
  --o-visualization "${DENOISE_DIR}/table.qzv" \
  --m-sample-metadata-file ${MAPPING_FP}

qiime feature-table tabulate-seqs \
  --i-data "${DENOISE_DIR}/rep-seqs.qza" \
  --o-visualization "${DENOISE_DIR}/rep-seqs.qzv"

qiime tools export \
  "${DENOISE_DIR}/table.qzv" \
  --output-dir "${DENOISE_DIR}/table-qzv"

qiime tools export \
  "${DENOISE_DIR}/table.qza" \
  --output-dir "${DENOISE_DIR}/table-qza"

###=====================
###  TAXONOMIC ANALYSIS
###=====================

qiime feature-classifier classify-sklearn \
  --i-classifier ${CLASSIFIER_FP} \
  --i-reads "${DENOISE_DIR}/rep-seqs.qza" \
  --o-classification "${DENOISE_DIR}/taxonomy.qza"

qiime metadata tabulate \
  --m-input-file "${DENOISE_DIR}/taxonomy.qza" \
  --o-visualization "${DENOISE_DIR}/taxonomy.qzv"

qiime tools export \
  "${DENOISE_DIR}/taxonomy.qza" \
  --output-dir "${DENOISE_DIR}/taxonomy-qza"

###=====================
###  GENERATE TREES
###=====================

qiime alignment mafft \
  --i-sequences "${DENOISE_DIR}/rep-seqs.qza" \
  --o-alignment "${DENOISE_DIR}/aligned-rep-seqs.qza"

qiime alignment mask \
  --i-alignment "${DENOISE_DIR}/aligned-rep-seqs.qza" \
  --o-masked-alignment "${DENOISE_DIR}/masked-aligned-rep-seqs.qza"

qiime phylogeny fasttree \
  --i-alignment "${DENOISE_DIR}/masked-aligned-rep-seqs.qza" \
  --o-tree "${DENOISE_DIR}/unrooted-tree.qza"

qiime phylogeny midpoint-root \
  --i-tree "${DENOISE_DIR}/unrooted-tree.qza" \
  --o-rooted-tree "${DENOISE_DIR}/rooted-tree.qza"

###=====================
###  ALPHA AND BETA DIVERSITY
###=====================

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "${DENOISE_DIR}/rooted-tree.qza" \
  --i-table "${DENOISE_DIR}/table.qza" \
  --p-sampling-depth 10 \
  --m-metadata-file ${MAPPING_FP} \
  --output-dir "${METRIC_DIR}"

qiime tools export \
  "${METRIC_DIR}/faith_pd_vector.qza" \
  --output-dir "${METRIC_DIR}/faith"

qiime tools export \
  "${METRIC_DIR}/unweighted_unifrac_distance_matrix.qza" \
  --output-dir "${METRIC_DIR}/uu"

qiime tools export \
  "${METRIC_DIR}/weighted_unifrac_distance_matrix.qza" \
  --output-dir "${METRIC_DIR}/wu"

###=====================
###  BIOM CONVERT
###=====================

biom convert \
  -i "${DENOISE_DIR}/table/feature-table.biom" \
  -o "${DENOISE_DIR}/table/feature-table.tsv" \
  --to-tsv

