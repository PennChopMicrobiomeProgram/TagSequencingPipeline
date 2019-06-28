#!/usr/bin/env bash

#$ -cwd
#$ -r n
#$ -j y
#$ -m ea
#$ -M danielsg@email.chop.edu
#$ -l h_vmem=50G
#$ -l m_mem_free=20G

source ~/.bashrc
conda activate qiime2-2018.11

./qiime2_pipeline.bash Robinson_Run_1_Metadata.tsv
