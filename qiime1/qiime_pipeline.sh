#!/bin/bash
set -x
set -e
set -u

if [ $# -ne 4 ]; then
    echo "Usage: $0 WORK_DIR MAPPING_FP LANE_NUM PROJECT_NAME"
    exit 1
fi

WORK_DIR="$1"
MAPPING_FP="$2"
LANE_NUM="$3"
PROJECT_NAME="$4"

DATA_DIR="${WORK_DIR}/data_files"
JOINED_DIR="${WORK_DIR}/qiime_joined"

PROJECT_DIR="${WORK_DIR}/${PROJECT_NAME}"
if [ ! -d $PROJECT_DIR ]; then
    mkdir $PROJECT_DIR
fi
cp $MAPPING_FP $PROJECT_DIR

LIBRARY_DIR="${PROJECT_DIR}/library"
OTU_DIR="${PROJECT_DIR}/otu"
BETA_DIR="${PROJECT_DIR}/beta_diversity"

QIIME_PARAMS_FP="${HOME}/qiime_parameters.txt"

FWD="${DATA_DIR}/Undetermined_S0_L$(printf "%03d" $LANE_NUM)_R1_001.fastq"
REV="${DATA_DIR}/Undetermined_S0_L$(printf "%03d" $LANE_NUM)_R2_001.fastq"
IDX="${DATA_DIR}/Undetermined_S0_L$(printf "%03d" $LANE_NUM)_I1_001.fastq"

# check if joined_dir exists, if it doesn't, execute
if [ ! -d $JOINED_DIR ]; then
    join_paired_ends.py \
	-f $FWD \
	-r $REV \
	-b $IDX \
	--min_overlap 35 \
	--perc_max_diff 15 \
	-o $JOINED_DIR
fi

split_libraries_fastq.py \
    -i "${JOINED_DIR}/fastqjoin.join.fastq" \
    -b "${JOINED_DIR}/fastqjoin.join_barcodes.fastq" \
    -m $MAPPING_FP \
    -q19 \
    -o $LIBRARY_DIR \
    --rev_comp_mapping_barcodes \
    --barcode_type golay_12 \
    --max_barcode_errors=0

pick_de_novo_otus.py \
    --input_fp "${LIBRARY_DIR}/seqs.fna" \
    --output_dir $OTU_DIR \
    -p $QIIME_PARAMS_FP

biom convert \
    -i "${OTU_DIR}/otu_table.biom"  \
    -o "${OTU_DIR}/otu_table.txt" \
    --to-tsv \
    -m $MAPPING_FP \
    --header-key=taxonomy \
    --output-metadata-id="Consensus Lineage" \
    --process-obs-metadata=taxonomy

beta_diversity_through_plots.py \
    --otu_table_fp "${OTU_DIR}/otu_table.biom" \
    --mapping_fp $MAPPING_FP \
    --parameter_fp $QIIME_PARAMS_FP \
    --output $BETA_DIR \
    --tree_fp "${OTU_DIR}/rep_set.tre"
