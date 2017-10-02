#!/bin/bash

if [ $# -ne 4 ]; then
	echo "Usage: $0 WORK_DIR MAPPING_FP LANE_NUM PROJECT_NAME"
	exit 1 
fi

WORK_DIR=$1
MAPPING_FP=$2
LANE_NUM=$3
PROJECT_NAME=$4

## GreenGene Training Set artiffact FP
## This can be downloaded from https://data.qiime2.org/2017.8/common/gg-13-8-99-nb-classifier.qza
gg_FP=$HOME"/gg-13-8-99-nb-classifier.qza"

DATA_DIR=$WORK_DIR"/data_files"

PROJECT_DIR=$WORK_DIR"/"$PROJECT_NAME
if [ ! -d $PROJECT_DIR ]; then
	mkdir $PROJECT_DIR
fi

FWD=$DATA_DIR"/Undetermined_S0_L$(printf "%03d" $LANE_NUM)_R1_001.fastq.gz"
REV=$DATA_DIR"/Undetermined_S0_L$(printf "%03d" $LANE_NUM)_R2_001.fastq.gz"
IDX=$DATA_DIR"/Undetermined_S0_L$(printf "%03d" $LANE_NUM)_I1_001.fastq.gz"

###=====================
### Maybe need a step to check mapping file using Keemei?
###=====================

###=====================
### DATA IMPORT 
###=====================

EMP_PAIRED_END_SEQUENCES_DIR=$PROJECT_DIR"/emp-paired-end-sequences"
if [ ! -d $EMP_PAIRED_END_SEQUENCES_DIR ]; then
        mkdir $EMP_PAIRED_END_SEQUENCES_DIR
fi

cp $FWD $EMP_PAIRED_END_SEQUENCES_DIR/"forward.fastq.gz"
cp $REV $EMP_PAIRED_END_SEQUENCES_DIR/"reverse.fastq.gz"
cp $IDX $EMP_PAIRED_END_SEQUENCES_DIR/"barcodes.fastq.gz"

qiime tools import \
  --type EMPPairedEndSequences \
  --input-path $EMP_PAIRED_END_SEQUENCES_DIR \
  --output-path $PROJECT_DIR"/emp-paired-end-sequences.qza"

###=====================
### DEMULTIPLEXING SEQUENCE
###=====================

qiime demux emp-paired \
  --m-barcodes-file $MAPPING_FP \
  --m-barcodes-category BarcodeSequence \
  --i-seqs $PROJECT_DIR"/emp-paired-end-sequences.qza" \
  --o-per-sample-sequences $PROJECT_DIR"/demux" \
  --p-rev-comp-mapping-barcodes 

qiime demux summarize \
  --i-data $PROJECT_DIR"/demux.qza" \
  --o-visualization $PROJECT_DIR"/demux.qzv" 

qiime tools export \
  $PROJECT_DIR"/demux.qzv" \
  --output-dir $PROJECT_DIR"/demux"  

###=====================
###  SEQUENCE QC AND FEATURE TABLE
###=====================

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $PROJECT_DIR"/demux.qza" \
  --p-trim-left-f 0 \
  --p-trunc-len-f 230 \
  --p-trim-left-r 0 \
  --p-trunc-len-r 230 \
  --p-n-threads 8 \
  --o-representative-sequences $PROJECT_DIR"/rep-seqs.qza" \
  --o-table $PROJECT_DIR"/table.qza"

qiime feature-table summarize \
  --i-table $PROJECT_DIR"/table.qza" \
  --o-visualization $PROJECT_DIR"/table.qzv" \
  --m-sample-metadata-file $MAPPING_FP 

qiime feature-table tabulate-seqs \
  --i-data $PROJECT_DIR"/rep-seqs.qza" \
  --o-visualization $PROJECT_DIR"/rep-seqs.qzv"

qiime tools export \
  $PROJECT_DIR"/table.qzv" \
  --output-dir $PROJECT_DIR"/table"

qiime tools export \
  $PROJECT_DIR"/table.qza" \
  --output-dir $PROJECT_DIR"/table"

###=====================
###  TAXONOMIC ANALYSIS
###=====================

qiime feature-classifier classify-sklearn \
  --i-classifier $gg_FP \
  --i-reads $PROJECT_DIR"/rep-seqs.qza" \
  --o-classification $PROJECT_DIR"/taxonomy.qza"

qiime metadata tabulate \
  --m-input-file $PROJECT_DIR"/taxonomy.qza" \
  --o-visualization $PROJECT_DIR"/taxonomy.qzv"

qiime tools export \
  $PROJECT_DIR"/taxonomy.qza" \
  --output-dir $PROJECT_DIR"/taxonomy"

###=====================
###  GENERATE TREES
###=====================

qiime alignment mafft \
  --i-sequences $PROJECT_DIR"/rep-seqs.qza" \
  --o-alignment $PROJECT_DIR"/aligned-rep-seqs.qza"

qiime alignment mask \
  --i-alignment $PROJECT_DIR"/aligned-rep-seqs.qza" \
  --o-masked-alignment $PROJECT_DIR"/masked-aligned-rep-seqs.qza"

qiime phylogeny fasttree \
  --i-alignment $PROJECT_DIR"/masked-aligned-rep-seqs.qza" \
  --o-tree $PROJECT_DIR"/unrooted-tree.qza"

qiime phylogeny midpoint-root \
  --i-tree $PROJECT_DIR"/unrooted-tree.qza" \
  --o-rooted-tree $PROJECT_DIR"/rooted-tree.qza"

###=====================
###  ALPHA AND BETA DIVERSITY
###=====================

qiime diversity core-metrics \
  --i-phylogeny $PROJECT_DIR"/rooted-tree.qza" \
  --i-table $PROJECT_DIR"/table.qza" \
  --p-sampling-depth 1000 \
  --output-dir $PROJECT_DIR"/core-metrics-results"

qiime tools export \
  $PROJECT_DIR"/core-metrics-results/observed_otus_vector.qza" \
  --output-dir $PROJECT_DIR"/core-metrics-results/observed_otu"

qiime tools export \
  $PROJECT_DIR"/core-metrics-results/shannon_vector.qza" \
  --output-dir $PROJECT_DIR"/core-metrics-results/shannon"

qiime tools export \
  $PROJECT_DIR"/core-metrics-results/unweighted_unifrac_distance_matrix.qza" \
  --output-dir $PROJECT_DIR"/core-metrics-results/uu"

qiime tools export \
  $PROJECT_DIR"/core-metrics-results/weighted_unifrac_distance_matrix.qza" \
  --output-dir $PROJECT_DIR"/core-metrics-results/wu"


###=====================
###  BIOM CONVERT
###=====================

biom convert \
  -i $PROJECT_DIR"/table/feature-table.biom" \
  -o $PROJECT_DIR"/table/feature-table.tsv" \
  --to-tsv \
  -m $MAPPING_FP \
  --header-key=taxonomy \
  --process-obs-metadata=taxonomy






