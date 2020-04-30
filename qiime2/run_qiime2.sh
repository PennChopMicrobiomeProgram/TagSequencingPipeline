#!/usr/bin/env bash

#$ -cwd
#$ -r n
#$ -j y
#$ -m ea
#$ -M danielsg@email.chop.edu
#$ -l h_vmem=50G
#$ -l m_mem_free=20G

source ~/.bashrc
conda activate qiime2-2019.4

export PRJ_DIR="/mnt/isilon/microbiome/analysis/danielsg/sehgal_run_6"

export METADATA="${PRJ_DIR}/Data/sehgal_run_6_mapping_file.tsv"

./qiime2_pipeline.bash $METADATA
