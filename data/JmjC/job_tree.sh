#!/bin/bash
# ==============================================================================
# Software:
    # diamond (v2.1.10)
    # mcl (v22.282)
    # mafft (v7.525)
    # clipkit (v1.4.1)
    # iqtree (v2.4.0)
# ==============================================================================

# Directories
DIR="/Directory/to/data/JmjC"
ANALYSIS="${DIR}/analysis"

# Define function
tree() {

    # local variables
    local GRP=$1; local DOM=$2;
    local REF=$([ "$GRP" == "PLT" ] && echo "HsapAtha" || echo "Hsap")
    local DIR_GRP="$ANALYSIS/orthologs/$REF/$DOM/$GRP"
    local FASTA_UNIQ="$DIR_GRP/${GRP}.${DOM}.seqs.fasta"
    
    # sequence clustering
    echo "## Tree construction | $GRP | $DOM | sequence clustering (diamond, mcl)"
    diamond makedb --in $FASTA_UNIQ -d $FASTA_UNIQ --quiet
    diamond blastp --more-sensitive --max-target-seqs 100 -d $FASTA_UNIQ -q $FASTA_UNIQ -o $DIR_GRP/${GRP}.${DOM}.diamond.csv --quiet --threads 8
    
    # generate abc format using bitscore ($12)
    awk '{print $1, $2, $12}' $DIR_GRP/${GRP}.${DOM}.diamond.csv > $DIR_GRP/${GRP}.${DOM}.diamond.abc
    mcl $DIR_GRP/${GRP}.${DOM}.diamond.abc --abc -I 1.5 -o $DIR_GRP/${GRP}.${DOM}.mcl.csv
    
    # multiple sequence alignment
    echo "## Tree construction | $GRP | $DOM | multiple sequence alignment (mafft)"
    mafft --genafpair --thread 8 --reorder --maxiterate 10000 $FASTA_UNIQ > $DIR_GRP/${GRP}.${DOM}.seqs.l.fasta
    
    # alignment trimming
    echo "## Tree construction | $GRP | $DOM | alignment trimming (clipkit)"
    clipkit $DIR_GRP/${GRP}.${DOM}.seqs.l.fasta -m kpic-gappy -o $DIR_GRP/${GRP}.${DOM}.seqs.lt.fasta -g 0.7

    # tree construction
    echo "## Tree construction | $GRP | $DOM | tree construction (iqtree2)"
    
    rm -f "$DIR_GRP/${GRP}.${DOM}.iqtree".* # old files will halt following processes. delete them if they are present.
    
    iqtree2 \
        -s $DIR_GRP/${GRP}.${DOM}.seqs.lt.fasta \
        -m TEST \
        -mset LG,WAG,JTT \
        -bb 1000 \
        -nt AUTO \
        -ntmax 8 \
        -pre $DIR_GRP/${GRP}.${DOM}.iqtree \
        -nstop 200 \
        -cptime 1800 \
        --quiet
    
    # tree shrink
    echo "## Tree construction | $GRP | $DOM | treeshrink"
    run_treeshrink.py -c -t $DIR_GRP/${GRP}.${DOM}.iqtree.treefile -m per-gene -q 0.05 -s 10,1 -f
    
}

# Run analyses
tree $1 $2

