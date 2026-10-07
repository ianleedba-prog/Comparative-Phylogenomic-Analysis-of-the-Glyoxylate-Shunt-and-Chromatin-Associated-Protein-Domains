#!/bin/bash
# ==============================================================================
# Software:
    # ete3 (ETE ToolKit) (v3.1.3)
    # possvm
# ==============================================================================

# Directories
DIR="/Directory/to/data/JmjC"
ANALYSIS="${DIR}/analysis"
POSSVM="${DIR}/possvm-orthology-master/possvm.py" # possvm script downloaded from https://github.com/xgrau/possvm-orthology

# Define functions
possvm() {
    local REF=$1
    local DOM=$2
    local GEN=$3
    
    if [ "$REF" == "HsapAtha" ]; then
        REFSPS="Hsap,Atha"
    else
        REFSPS="$REF"
    fi
    
    # 1. possvm
    for i in ${ANALYSIS}/gene_trees/*${DOM}*treefile ; do
        python $POSSVM \
        -i $i \
        -o ${ANALYSIS}/gene_trees/ \
        -r ${GEN}  \
        -refsps ${REFSPS} \
        --skipprint \
        -itermidroot 10 \
        -min_support_transfer 50 \
        -cut_gene_names 100 \
        -ogprefix $(basename $i | cut -f2,3 -d '.'). \
        -p $(basename $i | sed "s/.seqs.iqtree.treefile/.possvm.${REF}/")
    done
    
    # 2. process result
    process_orthologs $REF $DOM
}

process_orthologs() {
    local REF=$1
    local DOM=$2

    DOM_DIR="${ANALYSIS}/orthologs/${REF}/${DOM}"
    
    mkdir -p ${DOM_DIR}
    COMPILED_ORT="${DOM_DIR}/${DOM}.${REF}.orthogroups.csv"
    
    echo -e "gene\torthogroup\torthologous_to" > $COMPILED_ORT

    # compile data
    for i in ${ANALYSIS}/gene_trees/Rep.${DOM}*${REF}.ortholog_groups.csv; do
        awk 'NR > 1' "$i"
    done >> $COMPILED_ORT

    # unique names of orthogroups
    awk -F'\t' 'NR > 1 {print $2}' "$COMPILED_ORT" | sort -u > ${DOM_DIR}/${DOM}.${REF}.orthogroups.unique.txt
    
    # store only significant orthogroups
    FILTERED_TYPES="${DOM_DIR}/${DOM}.${REF}.orthogroups.unique.filt.txt"
    
    if [ "$DOM" == "JmjC" ]; then
        # Remove that are not histone demethylase
        grep -vE "HSPBAP1|JARID2|HIF1A|TYW5|JMJD4|JMJD8|RIOX1|RIOX2|NA" \
        ${DOM_DIR}/${DOM}.${REF}.orthogroups.unique.txt \
        > ${FILTERED_TYPES}
    elif [ "$DOM" == "AOD" ]; then
        # Store only histone demethylases
        grep -E "KDM1A|KDM1B|FLD|LDL1|LDL2|LDL3" \
        ${DOM_DIR}/${DOM}.${REF}.orthogroups.unique.txt \
        > ${FILTERED_TYPES}
    elif [ "$DOM" == "Tet_JBP" ]; then
        cp ${DOM_DIR}/${DOM}.${REF}.orthogroups.unique.txt \
        ${FILTERED_TYPES}
    fi
    
    FILTERED_COMPILED_ORT="${DOM_DIR}/${DOM}.${REF}.orthogroups.filt.csv"
    while read clu; do awk -v clu="$clu" -F "\t" '$2==clu' "$COMPILED_ORT" ; done < "$FILTERED_TYPES" > "$FILTERED_COMPILED_ORT"
}

# Run analyses
possvm $1 $2 $3
