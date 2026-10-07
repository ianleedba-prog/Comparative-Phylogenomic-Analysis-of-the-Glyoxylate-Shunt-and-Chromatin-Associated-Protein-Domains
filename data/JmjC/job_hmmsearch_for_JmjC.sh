#!/bin/bash
# ==============================================================================
# Software:
    # hmmer (v3.4)
    # diamond (v2.1.10)
    # mcl (v22.282)
    # bedtools (v2.31.1)
    # samtools (v1.21)
# ==============================================================================

# Directories
dir="/Directory/to/data/JmjC"
pfam_db="$dir/Pfam"

# Additional scripts
job_pfamscan="$dir/job_pfamscan.sh"
job_alignment="$dir/job_alignment.sh"

# Number of cpu
n_cpu="8"

# Check input status
if [ -z ${10} ] ; then
    echo "missing arguments:"
    echo "1. family id. e.g. homeobox"
    echo "2. comma-separated list of HMM files, eg. Homeodomain,Heomeobox_KN"
    echo "3. inflation for MCL, e.g. 1.1"
    echo "4. min alignment size for phhylogenetic analysis, e.g. 10"
    echo "5. Homology Cluster? na/HC"
    echo "6. database to search e.g. Compiled.All.fasta"
    echo "7. Do MCLclust+MSA+phylogeneis? Y/N"
    echo "8. All species (All) or Representative species (Rep)"
    echo "9. Out directory"
    echo "10. Log file"
    exit
fi

# Variables
fam_id=$1
hmm_syn=$(echo $2 | tr ',' ' ')
inflation=$3
min_alignment=$4
homology_cluster=$5
dbfile=$6
dophy=$7
prefix=$8
outdir=$9
log_file=$10

alignments="${outdir}/alignments/"; mkdir -p ${alignments}
searches="${outdir}/searches/"; mkdir -p ${searches}

# Define function
function do_hmmsearch {

    local hmm=$1
    local out=$2
    local fas=$3
    local cpu=$4
    local thr=$5

    if [ $thr == "GA" ] ; then

    hmmsearch \
        --domtblout ${out}.domtable \
        --cut_ga \
        --cpu ${cpu} \
        ${hmm} \
        ${fas} 1> /dev/null
    else

     hmmsearch \
        --domtblout ${out}.domtable \
        --domE ${thr} \
        --cpu ${cpu} \
        ${hmm} \
        ${fas} 1> /dev/null
    fi

    grep -v "^#" ${out}.domtable \
    | awk 'BEGIN { OFS="\t" } { print $1,$18,$19,$4,$5,$12 }' \
    > ${out}.domtable.csv
    
}

function do_homology_clusters {

    local fas=$1
    local out=$2
    local cpu=$3
    local inf=$4
    local min=$5

    #pairwise alignments
    diamond makedb --in ${fas} -d ${fas} --quiet
    diamond blastp \
        --more-sensitive \
        --max-target-seqs 100 \
        -d ${fas} \
        -q ${fas} \
        -o ${out}.diamond.csv \
        --quiet \
        --threads ${cpu}

    #partition sequences with MCL(Markov Cluster Algorithm) (edge weights are alignment bitscores)
    #low inflation value gives large, non-granular clusters
    awk '{ print $1,$2,$12 }' ${out}.diamond.csv > ${out}.diamond.abc
    mcl ${out}.diamond.abc --abc -I ${inf} -o ${out}.mcl.csv 2> /dev/null

    #compress stuff
    gzip -f ${out}.diamond.abc
    gzip -f ${out}.diamond.csv

    #ignore homology clusters with <X seqs
    awk 'NF > 2' ${out}.mcl.csv > ${out}.mcl_filtered.csv

    #remove clusters with <X seqs
    awk '{ for(i =1; i <= NF; i++) { print "HG"NR,$i } }' ${out}.mcl_filtered.csv > ${out}.mcl_filtered.txt

}

# Run analyses
echo "## hmmsearch & alignment | $(basename ${dbfile}): ${fam_id} | HMM search" | tee -a $log_file
for hmm in ${hmm_syn} ; do
    do_hmmsearch \
        ${pfam_db}/hmm/${hmm}.hmm \
        ${searches}/${prefix}.hmmsearch.domain_${hmm} \
        ${dbfile} \
        ${n_cpu} \
        GA
    cat ${searches}/${prefix}.hmmsearch.domain_${hmm}.domtable.csv \
    >> ${searches}/${prefix}.${fam_id}.domains.csv.tmp
done

#any hits?
cut -f 1 ${searches}/${prefix}.${fam_id}.domains.csv.tmp | sort -u > ${searches}/${prefix}.${fam_id}.genes.list
num_hit_genes=$(cat ${searches}/${prefix}.${fam_id}.genes.list | wc -l)
echo "# $(basename ${dbfile}): ${fam_id} | #genes found = ${num_hit_genes}" | tee -a $log_file

if [ ${num_hit_genes} -eq 0 ] ; then

    echo "## hmmsearch & alignment | $(basename ${dbfile}): ${fam_id} | Omit downstream analyses" | tee -a $log_file

else

    # merge overlapping domains
    bedtools merge \
    -i <(sort -k1,1 -k 2,2n ${searches}/${prefix}.${fam_id}.domains.csv.tmp) \
    -c 4 -o collapse \
    -d 0 \
    > ${searches}/${prefix}.${fam_id}.domains.csv.tmp2 \
    && rm -f ${searches}/${prefix}.${fam_id}.domains.csv.tmp

    # report
    num_unique_domns=$(cut -f1 ${searches}/${prefix}.${fam_id}.domains.csv.tmp2 | wc -l)
    echo "## hmmsearch & alignment | $(basename ${dbfile}): ${fam_id} | # unique domains = $num_unique_domns" | tee -a $log_file

    # extract complete sequences
    esl-sfetch -f ${dbfile} ${searches}/${prefix}.${fam_id}.genes.list \
    > ${searches}/${prefix}.${fam_id}.seqs.fasta 2>> $log_file

    #dict sequence lengths
    samtools faidx ${searches}/${prefix}.${fam_id}.seqs.fasta

    #expand domain region by a fixed amount of aa
    bedtools slop \
        -i ${searches}/${prefix}.${fam_id}.domains.csv.tmp2 \
        -g ${searches}/${prefix}.${fam_id}.seqs.fasta.fai \
        -b 10 \
        | awk '{ print $1, $2+1, $3, $4 }' \
        > ${searches}/${prefix}.${fam_id}.domains.csv \
        && rm ${searches}/${prefix}.${fam_id}.domains.csv.tmp2

    #extract domain region
    esl-sfetch --index ${searches}/${prefix}.${fam_id}.seqs.fasta 1>> $log_file 2>> $log_file
    esl-sfetch -C -f ${searches}/${prefix}.${fam_id}.seqs.fasta \
        <(awk '{ print $1"_"$2"-"$3, $2, $3, $1 }' ${searches}/${prefix}.${fam_id}.domains.csv) \
        > ${searches}/${prefix}.${fam_id}.domains.fasta 2>> $log_file
   
    if [ "${dophy}" == "Y" ] && [ "${homology_cluster}" == "HC" ]; then

    #find clusters of homology
    echo "## hmmsearch & alignment | $(basename ${dbfile}): ${fam_id} | pairwise alignements + MCL..." | tee -a $log_file
    do_homology_clusters \
        ${searches}/${prefix}.${fam_id}.seqs.fasta \
        ${alignments}/${prefix}.${fam_id} \
        ${n_cpu} \
        ${inflation} \
        ${min_alignment_size} 2>> $log_file

    #report
    num_homgroups=$( awk 'END { print NR }' ${alignments}/${prefix}.${fam_id}.mcl_filtered.csv)
    num_seqs_in_hg=$( awk 'END { print NR }' ${alignments}/${prefix}.${fam_id}.mcl_filtered.txt)
    echo "## hmmsearch & alignment | $(basename ${dbfile}): ${fam_id} | homology groups = ${num_homgroups} (${num_seqs_in_hg} / ${num_hit_genes})" | tee -a $log_file

    #pfamscan
    echo "## hmmsearch & alignment | $(basename ${dbfile}): ${fam_id} | Submit pfamscan job | pf.${fam_id}" | tee -a $log_file
    bash $job_pfamscan ${searches}/${prefix}.${fam_id}.seqs.fasta ${pfam_db} ${n_cpu} 1>> $log_file 2>> $log_file

    # obtain lists of clusters
    clu_list=$(awk '{ print $1 }' ${alignments}/${prefix}.${fam_id}.mcl_filtered.txt | sort -u -V)

        for clu in $clu_list; do

            # retrieve sequence list
            awk -v clu="$clu" '$1 == clu { print $2 }' ${alignments}/${prefix}.${fam_id}.mcl_filtered.txt \
            > ${alignments}/${prefix}.${fam_id}.${clu}.seqs.list

            # retrieve sequences from original fasta
            xargs -n 1 samtools faidx ${searches}/${prefix}.${fam_id}.seqs.fasta \
            < ${alignments}/${prefix}.${fam_id}.${clu}.seqs.list \
            > ${alignments}/${prefix}.${fam_id}.${clu}.seqs.fasta \
            2>> $log_file

            # report
            num_seqs_in_ali=$(grep -c ">" ${alignments}/${prefix}.${fam_id}.${clu}.seqs.fasta)
            echo "## hmmsearch & alignment | $(basename ${dbfile}): ${fam_id} | $clu | create fasta n = $num_seqs_in_ali" | tee -a $log_file

            # submit alignment job
            echo "## hmmsearch & alignment | $(basename ${dbfile}): ${fam_id} | $clu | submit alignment job ${fam_id}.${clu}" | tee -a $log_file
            bash $job_alignment ${alignments}/${prefix}.${fam_id}.${clu}.seqs.fasta ${n_cpu} 1>> $log_file 2>> $log_file
            
        done # end FOR loop that iterates over homology groups


    fi # end IF statement that allows you to skip phylogenies
    fi # end IF statement that excludes empty searches

exit 0
