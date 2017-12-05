#! /bin/bash
set -x 
set -e 
set -u 

if [ $# -ne 4 ]; then
    echo "Usage: $0 PROJECT_DIR OTU_ID TAXON_NAME OLIGO_NAME"
    exit 1
fi

PROJECT_DIR="$1"
OTU_ID="$2"
TAXON_NAME="$3"
OLIGO_NAME="$4"

# necessary QIIME directory 
LIBRARY_DIR="${PROJECT_DIR}/library"
OTU_DIR="${PROJECT_DIR}/otu"

# oligotyping result directory
OLIGO_DIR="${PROJECT_DIR}/${OLIGO_NAME}_${TAXON_NAME}"

# helper package
QIIME_PARAMS_FP="${HOME}/qiime_parameters.txt"
#Q2Oligo_FP="${HOME}/oligotyping/q2oligo"

# make the oligotyping result directory
mkdir -p ${OLIGO_DIR}

# get otu map
SEQS_MAP="${OTU_DIR}/uclust_picked_otus/seqs_otus.txt"
OTU_MAP="${OLIGO_DIR}/${OTU_ID}_otu_map.txt"
grep -w ${OTU_ID} ${SEQS_MAP} > ${OTU_MAP}

# get sequences
OTU_RAW_SEQS="${OLIGO_DIR}/${OTU_ID}_raw.fasta"
filter_fasta.py -f "${LIBRARY_DIR}/seqs.fna" -m "${OTU_MAP}" -o "${OTU_RAW_SEQS}"

# strip the header
OTU_STRIP_SEQS="${OLIGO_DIR}/${OTU_ID}_stripped.fasta"
sed 's/\s.*$//' "${OTU_RAW_SEQS}" "${OTU_STRIP_SEQS}"
#python "${Q2Oligo_FP}/stripMeta.py" "${OTU_RAW_SEQS}" "${OTU_STRIP_SEQS}"

## filter samples to keep
OTU_SEQS="${OLIGO_DIR}/${OTU_ID}.fasta"
SAMPLES_FP="${PROJECT_DIR}/samples_tokeep.txt"
filter_fasta.py -f "${OTU_STRIP_SEQS}" --sample_id_fp "${SAMPLES_FP}" -o "${OTU_SEQS}"

## all living tree project
LTP_FP="${PROJECT_DIR}/LTPs128/strain_aligned.fasta-TRIMMED"
TEMPLATE_FP="${OLIGO_DIR}/strain_aligned.fasta-TRIMMED"
cp "${LTP_FP}" "${OLIGO_DIR}"

## pynast alignment
OTU_PYNAST_SEQS="${OLIGO_DIR}/${OTU_ID}_pynast_aligned.fasta"
pynast -i "${OTU_SEQS}" -t "${TEMPLATE_FP}" -l 250

## trim uninformatics column
OTU_PYNAST_TRIM_SEQS="${OTU_PYNAST_SEQS}-TRIMMED"
o-trim-uninformative-columns-from-alignment "${OTU_PYNAST_SEQS}"

## position selection 
entropy-analysis "${OTU_PYNAST_TRIM_SEQS}"

## finally 
oligotype "${OTU_PYNAST_TRIM_SEQS}" "${OTU_PYNAST_TRIM_SEQS}-ENTROPY" -c 2 -M 1000
