#!/bin/bash
# ==============================================================================
# Software:
    # pfam_scan (v1.6)
    # bioawk (v1.0)
# ==============================================================================

# Define function
function do_pfamscan {

    # input
    fas=$1
    pfam_db=$2
    n_cpu=$3

    # pfamscan
    pfam_scan.pl  -fasta ${fas}  -dir ${pfam_db} -cpu ${n_cpu} > ${fas%%.fasta}.pfamscan.csv

    # list of proteins
    bioawk -c fastx '{ print $1 }' ${fas} > ${fas}_TMP_list

    # table with pfam architectures
    while read seqid ; do
    grep -w "$seqid" ${fas%%.fasta}.pfamscan.csv \
    | tr -s  ' ' '\t' | cut -f 7 | xargs
    done < ${fas}_TMP_list > ${fas}_TMP_architectures
    paste ${fas}_TMP_list ${fas}_TMP_architectures > ${fas%%.fasta}.pfamscan_archs.csv

    # clean
    rm ${fas}_TMP_architectures ${fas}_TMP_list

}

# Run analyses
do_pfamscan $1 $2 $3
