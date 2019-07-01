#!/usr/bin/env bash

set -x
set -e
#set -u

#source ~/.bashrc
#conda activate qiime2-2018.11

if [ $# -ne 1 ]; then
	echo "Usage: $0 MAPPING_FP"
    echo "MAPPING_FP is the mapping file for Qiime2"
    echo "Data files should be in a directory named "data_files""
	exit 1
fi

MAPPING_FP=$1
WORK_DIR="$(dirname ${MAPPING_FP})"

# See https://www.ostricher.com/2014/10/the-right-way-to-get-the-directory-of-a-bash-script/
SOURCE_REL="${BASH_SOURCE[0]}"
SOURCE_ABS=$(readlink -f "${SOURCE_REL}")
SOURCE_DIR=$(dirname "${SOURCE_ABS}")

### PATH TO Ceylan's CODE TO COMBINE I1 and I2
INDEX1_INDEX2_COMBINE_SCRIPT="$(dirname ${SOURCE_DIR})/combine_barcodes.py"

## Taxonomy classifier setup. Two classifiers are currently available:
## classifiers trained on full length and on 515F/806R region of Greengenes 13_8 99% OTUs
## These can be downloaded from https://data.qiime2.org/2017.9/common/gg-13-8-99-nb-classifier.qza (full length)
## or https://data.qiime2.org/2017.9/common/gg-13-8-99-515-806-nb-classifier.qza (515F/806R region)

CLASSIFIER_FP="gg-13-8-99-nb-classifier.qza"
#CLASSIFIER_FP="gg-13-8-99-515-806-nb-classifier.qza" ## used for V4 region
#CLASSIFIER_FP="gg-13-8-99-27-338-nb-classifier.qza" ## trained for V1V2 region truncated at 350 bp

EMP_PAIRED_END_SEQUENCES_DIR="${WORK_DIR}/emp_paired_end_sequences"
DATA_DIR="${WORK_DIR}/data_files"
DEMUX_DIR="${WORK_DIR}/demux_results"
DENOISE_DIR="${WORK_DIR}/denoising_results"
METRIC_DIR="${WORK_DIR}/core_metrics_results"
EXTRAS_DIR="${WORK_DIR}/extras"

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

if [ ! -e "${DATA_DIR}/barcodes.fastq.gz" ]; then
    python ${INDEX1_INDEX2_COMBINE_SCRIPT} --work-dir "${DATA_DIR}"
    gzip "${DATA_DIR}/barcodes.fastq"
fi

FWD="${DATA_DIR}/Undetermined_S0_L001_R1_001.fastq.gz"
REV="${DATA_DIR}/Undetermined_S0_L001_R2_001.fastq.gz"
IDX="${DATA_DIR}/barcodes.fastq.gz"


###=====================
### DATA IMPORT
###=====================

if [ ! -d ${EMP_PAIRED_END_SEQUENCES_DIR} ]; then
        mkdir ${EMP_PAIRED_END_SEQUENCES_DIR}
fi

if [ ! -e "${EMP_PAIRED_END_SEQUENCES_DIR}/forward.fastq.gz" ]; then
    mv ${FWD} "${EMP_PAIRED_END_SEQUENCES_DIR}/forward.fastq.gz"
    mv ${REV} "${EMP_PAIRED_END_SEQUENCES_DIR}/reverse.fastq.gz"
    mv ${IDX} "${EMP_PAIRED_END_SEQUENCES_DIR}"
fi

if [ ! -e "${WORK_DIR}/emp-paired-end-sequences.qza" ]; then
    qiime tools import \
      --type EMPPairedEndSequences \
      --input-path ${EMP_PAIRED_END_SEQUENCES_DIR} \
      --output-path "${WORK_DIR}/emp-paired-end-sequences.qza"
fi

###=====================
### DEMULTIPLEXING SEQUENCE
###=====================

if [ ! -e "${DEMUX_DIR}/per_sample_sequences.qza" ]; then
    qiime demux emp-paired \
      --m-barcodes-file ${MAPPING_FP} \
      --m-barcodes-column BarcodeSequence \
      --i-seqs "${WORK_DIR}/emp-paired-end-sequences.qza" \
      --p-rev-comp-mapping-barcodes \
      --output-dir "${DEMUX_DIR}"
fi

if [ ! -e "${DEMUX_DIR}/per_sample_sequences.qzv" ]; then
    qiime demux summarize \
      --i-data "${DEMUX_DIR}/per_sample_sequences.qza" \
      --o-visualization "${DEMUX_DIR}/per_sample_sequences.qzv"
fi

if [[ ! -e "${DEMUX_DIR}/demux" && -e "${DEMUX_DIR}/per_sample_sequences.qzv" ]]; then
    qiime tools export \
      --input-path "${DEMUX_DIR}/per_sample_sequences.qzv" \
      --output-path "${DEMUX_DIR}/demux"
fi

###=====================
###  SEQUENCE QC AND FEATURE TABLE
###=====================

if [ ! -d ${DENOISE_DIR} ]; then
        mkdir ${DENOISE_DIR}
fi

## discussion needed for denosing parameters below

if [ ! -e "${DENOISE_DIR}/table_main.qza" ]; then
    qiime dada2 denoise-paired \
      --i-demultiplexed-seqs "${DEMUX_DIR}/per_sample_sequences.qza" \
      --p-trim-left-f 0 \
      --p-trunc-len-f 230 \
      --p-trim-left-r 0 \
      --p-trunc-len-r 230 \
      --p-n-threads 8 \
      --o-representative-sequences "${DENOISE_DIR}/rep-seqs.qza" \
      --o-denoising-stats "${DENOISE_DIR}/denoising-stats.qza" \
      --o-table "${DENOISE_DIR}/table_main.qza"
fi

# filter out the samples with no featured from the table
if [ ! -e "${DENOISE_DIR}/table.qza" ]; then
    qiime feature-table filter-samples \
      --i-table "${DENOISE_DIR}/table_main.qza" \
      --p-min-frequency 1 \
      --o-filtered-table "${DENOISE_DIR}/table.qza"
fi

if [[ ! -e  "${DENOISE_DIR}/table.qzv" && -e "${DENOISE_DIR}/table.qza" ]]; then
    qiime feature-table summarize \
      --i-table "${DENOISE_DIR}/table.qza" \
      --o-visualization "${DENOISE_DIR}/table.qzv" \
      --m-sample-metadata-file ${MAPPING_FP}
fi

if [[ ! -e "${DENOISE_DIR}/rep-seqs.qzv" && -e "${DENOISE_DIR}/rep-seqs.qza" ]]; then
    qiime feature-table tabulate-seqs \
      --i-data "${DENOISE_DIR}/rep-seqs.qza" \
      --o-visualization "${DENOISE_DIR}/rep-seqs.qzv"
fi

if [[ ! -e "${DENOISE_DIR}/table" && -e "${DENOISE_DIR}/table.qza" ]]; then
    qiime tools export \
      --input-path "${DENOISE_DIR}/table.qza" \
      --output-path "${DENOISE_DIR}/table"
fi

if [[ ! -e "${DENOISE_DIR}/denoising-stats" && -e "${DENOISE_DIR}/denoising-stats.qza" ]]; then
    qiime tools export \
	  --input-path "${DENOISE_DIR}/denoising-stats.qza" \
	  --output-path "${DENOISE_DIR}/denoising-stats"
fi


###=====================
###  TAXONOMIC ANALYSIS
###=====================

if [ ! -e "${DENOISE_DIR}/taxonomy.qza" ]; then
    qiime feature-classifier classify-sklearn \
      --i-classifier ${CLASSIFIER_FP} \
      --i-reads "${DENOISE_DIR}/rep-seqs.qza" \
      --o-classification "${DENOISE_DIR}/taxonomy.qza"
fi

if [[ ! -e "${DENOISE_DIR}/taxonomy.qzv" && -e "${DENOISE_DIR}/taxonomy.qza" ]]; then
    qiime metadata tabulate \
      --m-input-file "${DENOISE_DIR}/taxonomy.qza" \
      --o-visualization "${DENOISE_DIR}/taxonomy.qzv"
fi

if [[ ! -e  "${DENOISE_DIR}/taxonomy" && -e "${DENOISE_DIR}/taxonomy.qza" ]]; then
    qiime tools export \
      --input-path "${DENOISE_DIR}/taxonomy.qza" \
      --output-path "${DENOISE_DIR}/taxonomy"
fi

###=====================
###  GENERATE TREES
###=====================

if [ ! -e "${DENOISE_DIR}/aligned-rep-seqs.qza" ]; then
    qiime alignment mafft \
      --i-sequences "${DENOISE_DIR}/rep-seqs.qza" \
      --o-alignment "${DENOISE_DIR}/aligned-rep-seqs.qza"
fi

if [ ! -e "${DENOISE_DIR}/masked-aligned-rep-seqs.qza" ]; then
    qiime alignment mask \
      --i-alignment "${DENOISE_DIR}/aligned-rep-seqs.qza" \
      --o-masked-alignment "${DENOISE_DIR}/masked-aligned-rep-seqs.qza"
fi

if [ ! -e "${DENOISE_DIR}/unrooted-tree.qza" ]; then
    qiime phylogeny fasttree \
      --i-alignment "${DENOISE_DIR}/masked-aligned-rep-seqs.qza" \
      --o-tree "${DENOISE_DIR}/unrooted-tree.qza"
fi

if [ ! -e "${DENOISE_DIR}/rooted-tree.qza" ]; then
    qiime phylogeny midpoint-root \
      --i-tree "${DENOISE_DIR}/unrooted-tree.qza" \
      --o-rooted-tree "${DENOISE_DIR}/rooted-tree.qza"
fi

###=====================
###  ALPHA AND BETA DIVERSITY
###=====================

if [ ! -d ${METRIC_DIR} ]; then
        mkdir ${METRIC_DIR}
fi

if [ ! -e "${METRIC_DIR}/faith_pd_vector.qza" ]; then
    qiime diversity alpha-phylogenetic \
      --i-phylogeny "${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${DENOISE_DIR}/table.qza" \
      --p-metric faith_pd \
      --o-alpha-diversity "${METRIC_DIR}/faith_pd_vector.qza"
fi

if [[ ! -e "${METRIC_DIR}/faith" && -e "${METRIC_DIR}/faith_pd_vector.qza" ]]; then
    qiime tools export \
      --input-path "${METRIC_DIR}/faith_pd_vector.qza" \
      --output-path "${METRIC_DIR}/faith"
fi

if [ ! -e "${METRIC_DIR}/weighted_unifrac_distance_matrix.qza" ]; then
    qiime diversity beta-phylogenetic \
      --i-phylogeny "${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${DENOISE_DIR}/table.qza" \
      --p-metric weighted_unifrac \
      --o-distance-matrix "${METRIC_DIR}/weighted_unifrac_distance_matrix.qza"
fi

qiime tools export \
  --input-path "${METRIC_DIR}/weighted_unifrac_distance_matrix.qza" \
  --output-path "${METRIC_DIR}/wu"

if [ ! -e "${METRIC_DIR}/unweighted_unifrac_distance_matrix.qza" ]; then
    qiime diversity beta-phylogenetic \
      --i-phylogeny "${DENOISE_DIR}/rooted-tree.qza" \
      --i-table "${DENOISE_DIR}/table.qza" \
      --p-metric unweighted_unifrac \
      --o-distance-matrix "${METRIC_DIR}/unweighted_unifrac_distance_matrix.qza"
fi

qiime tools export \
  --input-path "${METRIC_DIR}/unweighted_unifrac_distance_matrix.qza" \
  --output-path "${METRIC_DIR}/uu"

###=====================
###  BIOM CONVERT
###=====================

if [ ! -e "${DENOISE_DIR}/table/feature-table.tsv" ]; then
    biom convert \
      -i "${DENOISE_DIR}/table/feature-table.biom" \
      -o "${DENOISE_DIR}/table/feature-table.tsv" \
      --to-tsv
fi

###=====================
###  ADDITIONAL FIGURES
###=====================


#if [ ! -d "${EXTRAS_DIR}/taxa_barplot" ]; then
#    qiime taxa barplot \
#      --i-table "${DENOISE_DIR}/table.qza" \
#      --i-taxonomy "${DENOISE_DIR}/taxonomy.qza" \
#      --m-metadata-file ${MAPPING_FP} \
#      --output-dir "${EXTRAS_DIR}/taxa_barplot"
#fi


#if [ ! -d "${EXTRAS_DIR}/core-metrics-phylogenetic" ]; then
#qiime diversity core-metrics-phylogenetic \
#      --p-sampling-depth 1000 \
#      --i-phylogeny "${DENOISE_DIR}/rooted-tree.qza" \
#      --i-table "${DENOISE_DIR}/table.qza" \
#      --m-metadata-file "${MAPPING_FP}" \
#      --p-n-jobs 8 \
#      --o-observed-otus-vector "${EXTRAS_DIR}/obs_otus.qza" \
#      --output-dir "${EXTRAS_DIR}/core-metrics-phylogenetic"
#fi

#if [ ! -d "${EXTRAS_DIR}/alpha-group-signf-otus" ]; then
#qiime diversity alpha-group-significance \
#      --i-alpha-diversity "${EXTRAS_DIR}/obs_otus.qza" \
#      --m-metadata-file "${MAPPING_FP}" \
#      --output-dir "${EXTRAS_DIR}/alpha-group-signf-otus"
#fi


echo "****"
echo "Done!"
echo "If nothing happened,"
echo "You may already have result files, check your directories!"



