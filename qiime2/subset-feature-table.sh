#!/usr/bin/env bash

set -x

qiime feature-table filter-samples \
    --i-table ../../denoising-results/table.qza \
    --m-metadata-file ./for_forum/subset_metadata2.txt \
    --o-filtered-table ./for_forum/subset_table2.qza
    
# the first subset:
#    --p-where '"#SampleID" IN ("geneblock5","Feed","Extractemptywell6","DNAfreewater1","Day36.TC2.L","Day36.TC1.L","Day36.DC2.L","Day36.DC1.L","Day36.AC2.L")' \
# which actually worked!
