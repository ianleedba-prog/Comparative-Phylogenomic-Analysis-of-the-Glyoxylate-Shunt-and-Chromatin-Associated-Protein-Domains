#!/bin/bash
# ==============================================================================
# Software:
    # hmmer (v3.4)
    # bioawk
# ==============================================================================

# Directories
DIR="/Directory/to/data/"
JmjC="${DIR}/JmjC"
FASTA="${DIR}/protein_fa"
COMPILED="${JmjC}/compiled_fa"
ANN="${DIR}/annotation"
ANALYSIS="${JmjC}/analysis"

# Metadata for representative species
METADATA="${JmjC}/Metadata_rep.species_info.txt"

# ==============================================================================
# 1. FASTA compilation (proceed when all files are prepared. (NCBI/WormBase/Ensembl/+other.publications)
# ==============================================================================
mkdir -p ${COMPILED}

# log file
COMPILE_LOG="$JmjC/s01_compile_rep.log"
> $COMPILE_LOG

# failed list
COMPILE_MISSING_LOG="$JmjC/s01_compile_rep.failed.list.txt"
echo -e "Taxonomy_name\tAbbreviation\tNCBI_Taxonomy_ID\tSource\tFile\tAnnotation" > $COMPILE_MISSING_LOG

# checking list
echo "## FASTA compilation | Rep. species (n = 118) | Checking all FASTA files within directory." | tee -a $COMPILE_LOG
tail -n +2 $METADATA | while IFS=$'\t' read -r tax_name abbr ncbi_id file_source file_name annotation; do
    faa=$(find $FASTA -maxdepth 1 -name "*${file_name}*.gz" | head -n 1)
    if [[ ! -f $faa ]]; then
        echo -e "$tax_name\t$abbr\t$ncbi_id\t$file_source\t$file_name\t$annotation" >> $COMPILE_MISSING_LOG
    fi
done

# halt process if missing files are detected
if [ $(wc -l < $COMPILE_MISSING_LOG) -gt 1 ]; then
    echo "## FASTA compilation | Rep. species (n = 118) | Missing files detected ❌" | tee -a $COMPILE_LOG
    echo "## FASTA compilation | Rep. species (n = 118) | Please check missing list & refine the data. | Exiting Script." | tee -a $COMPILE_LOG
    exit 1
fi

# compilation begin
echo "## FASTA compilation | Rep. species (n = 118) | All files verified." | tee -a $COMPILE_LOG
echo "## FASTA compilation | Rep. species (n = 118) | Starting compilation." | tee -a $COMPILE_LOG

> ${COMPILED}/Compiled.Rep.fasta
tail -n +2 $METADATA | while IFS=$'\t' read -r tax_name abbr ncbi_id file_source file_name annotation; do
    faa=$(find $FASTA -maxdepth 1 -name "*${file_name}*.gz" | head -n 1)
    gzcat $faa >> ${COMPILED}/Compiled.Rep.fasta
    echo "## FASTA compilation | Rep. species (n = 118) | $tax_name | DONE ✅" >> $COMPILE_LOG
done

echo "## FASTA compilation | Rep. species (n = 118) | Compilation done." | tee -a $COMPILE_LOG

# index newly compiled fasta
echo "## FASTA compilation | Rep. species (n = 118) | Indexing compiled FASTA." | tee -a $COMPILE_LOG
esl-sfetch --index ${COMPILED}/Compiled.Rep.fasta

# remove failed list if no record
if [ -s $COMPILE_MISSING_LOG ] && [ $(wc -l < $COMPILE_MISSING_LOG) -eq 1 ]; then
    rm $COMPILE_MISSING_LOG
fi

echo "## FASTA compilation | Rep. species (n = 118) | Completed ✅" | tee -a $COMPILE_LOG

# ==============================================================================
# 2. Rename annotation (proceed when all files are prepared. (NCBI/WormBase/Ensembl/+other.publications)
# ==============================================================================
# log file
RENAME_LOG="$JmjC/s01_rename_annotation.log"
> $RENAME_LOG

# failed list
RENAME_MISSING_LOG="$JmjC/s01_rename_annotation.missing.list.txt"
echo -e "Taxonomy_name\tAbbreviation\tNCBI_Taxonomy_ID\tSource\tFile\tAnnotation" > $RENAME_MISSING_LOG

# checking list
echo "## Rename Annotation File | Rep. species (n = 118) | Checking all annotation files within directory." | tee -a $RENAME_LOG
tail -n +2 $METADATA | while IFS=$'\t' read -r tax_name abbr ncbi_id file_source file_name annotation; do
    
    original_file=$ANN/$annotation
    renamed_file=$ANN/${abbr}_${annotation}

    # skip na (Tax_id: 749232, 433461, 749231)
    if [[ $annotation == "na" ]]; then
        continue
    fi
    
    # record missing files
    if [[ ! -f $original_file && ! -f $renamed_file ]]; then
        echo -e "$tax_name\t$abbr\t$ncbi_id\t$file_source\t$file_name\t$annotation" >> $RENAME_MISSING_LOG
    fi
done

# halt process if missing files are detected
if [ $(wc -l < $RENAME_MISSING_LOG) -gt 1 ]; then
    echo "## Rename Annotation File | Rep. species (n = 118) | Missing files detected ❌" | tee -a $RENAME_LOG
    echo "## Rename Annotation File | Rep. species (n = 118) | Please check missing list & refine the data. | Exiting Script." | tee -a $RENAME_LOG
    exit 1
fi

# rename begin
echo "## Rename Annotation File | Rep. species (n = 118) | All files verified." | tee -a $RENAME_LOG
echo "## Rename Annotation File | Rep. species (n = 118) | Starting to rename annotations." | tee -a $RENAME_LOG

tail -n +2 $METADATA | while IFS=$'\t' read -r tax_name abbr ncbi_id file_source file_name annotation; do

    original_file=$ANN/$annotation
    renamed_file=$ANN/${abbr}_${annotation}

    # skip na (Tax_id: 749232, 433461, 749231)
    if [[ $annotation == "na" ]]; then
        continue
    fi

    # rename
    if [[ -f $original_file ]]; then
        mv $original_file $renamed_file
        echo "## Rename Annotation File | Rep. species (n = 118) | $tax_name | $annotation -> ${abbr}_${annotation} | DONE ✅" >> $RENAME_LOG
    elif [[ -f $renamed_file ]]; then
        echo "## Rename Annotation File | Rep. species (n = 118) | $tax_name | ${abbr}_${annotation}: File exists | SKIP ⏭️" >> $RENAME_LOG
    fi
done

echo "## Rename Annotation File | Rep. species (n = 118) | Completed ✅" | tee -a $RENAME_LOG

# remove missing list with no record
if [ -s $RENAME_MISSING_LOG ] && [ $(wc -l < $RENAME_MISSING_LOG) -eq 1 ]; then
    rm $RENAME_MISSING_LOG
fi

# ==============================================================================
# 3. hmmsearch & alignment
# ==============================================================================
mkdir -p ${ANALYSIS}

# protein domain list
DOM_LIST="${DIR}/protein_domains_hmm.csv"

# log file
HMM_LOG="$JmjC/s01_hmmsearch_alignment.Rep.log"
> $HMM_LOG

# hmmsearch
echo "## hmmsearch & alignment | Rep. species (n = 118) | Begin" | tee -a $HMM_LOG
while read -a i ; do
    bash ${JmjC}/job_hmmsearch_for_JmjC.sh ${i[2]} ${i[3]} ${i[5]} ${i[6]} ${i[4]} ${COMPILED}/Compiled.Rep.fasta Y Rep ${ANALYSIS} $HMM_LOG
done < ${DOM_LIST}
echo "## hmmsearch & alignment | Rep. species (n = 118) | Completed ✅" | tee -a $HMM_LOG

# ==============================================================================
# 4. Gene counts
# ==============================================================================
mkdir -p ${ANALYSIS}/gene_counts

# summarize gene counts into single file
while read -a i ; do
    cat ${ANALYSIS}/searches/Rep.${i[2]}.genes.list | sed "s/$/\t${i[2]}/"
done < ${DOM_LIST} > ${ANALYSIS}/gene_counts/Rep.genecounts.csv
echo "## Gene counts | Rep. species (n = 118) | Completed ✅"

# ==============================================================================
# 5. POSSVM (identify pairs and clusters of orthologous genes)
# ==============================================================================
mkdir -p ${ANALYSIS}/gene_trees
mkdir -p ${ANALYSIS}/gene_sequences

cp ${ANALYSIS}/alignments/*treefile ${ANALYSIS}/gene_trees/
cp ${ANALYSIS}/alignments/*.seqs.fasta ${ANALYSIS}/gene_sequences/

# ortholog idendification using POSSVM
REFS=("Hsap" "HsapAtha"); DOMS=("Tet_JBP" "JmjC" "AOD")
for ref in "${REFS[@]}"; do
    for dom in "${DOMS[@]}"; do
        echo "## Ortholog identification | ref: $ref | $dom"
        bash ${JmjC}/job_possvm.sh $ref $dom $JmjC/${ref}_gene_names.txt
    done
done

# ==============================================================================
# 6. Shrink isoforms
# ==============================================================================

# PEM = Premetazoa vs Early metazoa
# MET = Metazoans with vs without glyoxylate cycle
# PLT = Metazoa vs Embryophyta

# shrink isoforms into a single gene
GRPS=("PEM" "MET" "PLT"); DOMS=("Tet_JBP" "JmjC" "AOD")
for grp in "${GRPS[@]}"; do
    for dom in "${DOMS[@]}"; do
        echo "## Processing isoforms | $grp | $dom |"
        bash ${JmjC}/job_isoformshrink.sh $grp $dom $JmjC/Sps.prefix_${grp}.txt
    done
done

# ==============================================================================
# 7. Align & Tree construction
# ==============================================================================

# align orthologs and construct tree files
GRPS=("PEM" "MET" "PLT"); DOMS=("Tet_JBP" "JmjC" "AOD")
for grp in "${GRPS[@]}"; do
    for dom in "${DOMS[@]}"; do
        echo "## Tree construction | $grp | $dom |"
        bash ${JmjC}/job_tree.sh $grp $dom
    done
done

# ==============================================================================
# 8. Find conserved sites for JmjC-Histone demethylases
# ==============================================================================

# profiling conserved catalytic sites within KDM (JmjC-containing) orthologs
GRPS=("PEM" "MET"); DOMS=("JmjC")
for grp in "${GRPS[@]}"; do
    for dom in "${DOMS[@]}"; do
        echo "## Conserved sites identification (ref: Human) | $grp | $dom |"
        bash ${JmjC}/job_tree.sh $grp $dom
    done
done
