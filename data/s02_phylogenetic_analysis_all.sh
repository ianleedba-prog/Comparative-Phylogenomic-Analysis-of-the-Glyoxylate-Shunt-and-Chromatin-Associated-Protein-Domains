#!/bin/bash
# ==============================================================================
# Software:
    # busco (v6.0.0)
    # hmmer (v3.4)
# ==============================================================================

# Directories
DIR="/Directory/to/data"
FASTA="${DIR}/protein_fa"
COMPILED="${DIR}/compiled_fa"
ANALYSIS="${DIR}/analysis"

# Species metadata
METADATA="${DIR}/Metadata_species_info.txt"

# ==============================================================================
# 1. Validate assemblies (BUSCO) & Store only qualified or significant data
# ==============================================================================
echo "## BUSCO | Downloading eukaryota_odb10"
busco --download eukaryota_odb10 --quiet --download_path $DIR > /dev/null 2>&1

echo "## BUSCO | All species (n = 1609) | Calculating completeness scores of each assembly"
bash ${DIR}/job_busco.sh
echo "## BUSCO | All species (n = 1609) | Completed ✅"

# ==============================================================================
# 2. Prepare taxonomy database ** The result file (Taxonomy_info.txt) is provided by author.
# ==============================================================================
#echo "## Prepare taxonomy database |"
#TAX_RANK="${DIR}/taxonomy/taxonomy_ranked.tsv"
#sed "s/\t|//g" ${DIR}/taxonomy/rankedlineage.dmp > ${TAX_RANK}

# Match taxonomy rank to metadata species
#METADATA_TAX="${DIR}/Taxonomy_info.txt"
#UNMATCHED="${DIR}/Unmatched.txt"

#> "$METADATA_TAX"
#> "${DIR}/Unmatched.txt"

#echo -e "Taxonomy_name\tNCBI_Taxonomy_ID\tSpecies\tGenus\tFamily\tOrder\tClass\tPhylum\tKingdom\tSuperkingdom" > $METADATA_TAX
#awk -F'\t' -v matched="$METADATA_TAX" -v unmatched="$UNMATCHED" '
#    # 1. taxonomy_ranked.tsv
#    NR==FNR {
#        # tax information ($3~$10)
#        val = $3
#        for(i=4; i<=10; i++) {
#            val = val "\t" $i
#        }
#        target[$1] = val
#        next
#    }
#    # 2. Metadata
#    FNR == 1 { next }
#    {
#        if ($3 in target) {
#            print $1 "\t" $3 "\t" target[$3] >> matched
#        } else {
#            print $0 >> unmatched
#        }
#    }
#' "$TAX_RANK" "$METADATA"

# Check unmatched species. Remove if all species are matched.
#[ -f "$UNMATCHED" ] && [ ! -s "$UNMATCHED" ] && rm "$UNMATCHED"
#echo "## Prepare taxonomy database | All species matched | ✅"

# ==============================================================================
# 3. FASTA compilation
# ==============================================================================
mkdir -p ${COMPILED}

# log file
COMPILE_LOG="$DIR/s02a_compile_all.log"
> $COMPILE_LOG

# failed list
COMPILE_MISSING_LOG="$DIR/s02a_compile_all.failed.list.txt"
echo -e "Taxonomy_name\tAbbreviation\tNCBI_Taxonomy_ID\tSource\tFile" > $COMPILE_MISSING_LOG

# checking list
echo "## FASTA compilation | All species (n = 1609) | Checking all FASTA files within directory." | tee -a $COMPILE_LOG
tail -n +2 $METADATA | while IFS=$'\t' read -r tax_name abbr ncbi_id file_source file_name; do
    faa=$(find $FASTA -maxdepth 1 -name "*${file_name}*.gz" | head -n 1)
    if [[ ! -f $faa ]]; then
        echo -e "$tax_name\t$abbr\t$ncbi_id\t$file_source\t$file_name" >> $COMPILE_MISSING_LOG
    fi
done

# halt process if missing files are detected
if [ $(wc -l < $COMPILE_MISSING_LOG) -gt 1 ]; then
    echo "## FASTA compilation | All species (n = 1609) | Missing files detected ❌" | tee -a $COMPILE_LOG
    echo "## FASTA compilation | All species (n = 1609) | Please check missing list & refine the data. | Exiting Script." | tee -a $COMPILE_LOG
    exit 1
fi

# compilation begin
echo "## FASTA compilation | All species (n = 1609) | All files verified." | tee -a $COMPILE_LOG
echo "## FASTA compilation | All species (n = 1609) | Starting compilation." | tee -a $COMPILE_LOG

> ${COMPILED}/Compiled.All.fasta
tail -n +2 $METADATA | while IFS=$'\t' read -r tax_name abbr ncbi_id file_source file_name; do
    faa=$(find $FASTA -maxdepth 1 -name "*${file_name}*.gz" | head -n 1)
    gzcat $faa >> ${COMPILED}/Compiled.All.fasta
    echo "## FASTA compilation | All species (n = 1609) | $tax_name | DONE ✅" >> $COMPILE_LOG
done

echo "## FASTA compilation | All species (n = 1609) | Compilation done." | tee -a $COMPILE_LOG

# index newly compiled fasta
echo "## FASTA compilation | All species (n = 1609) | Indexing compiled FASTA" | tee -a $COMPILE_LOG
esl-sfetch --index ${COMPILED}/Compiled.All.fasta

# remove missing list with no record
if [ -s $COMPILE_MISSING_LOG ] && [ $(wc -l < $COMPILE_MISSING_LOG) -eq 1 ]; then
    rm $COMPILE_MISSING_LOG
fi

echo "## FASTA compilation | All species (n = 1609) | Completed ✅" | tee -a $COMPILE_LOG

# ==============================================================================
# 4. hmmsearch & alignment
# ==============================================================================
mkdir -p ${ANALYSIS}

# protein domain list
DOM_LIST="${DIR}/protein_domains_hmm.csv"

# log file
HMM_LOG="$DIR/s02a_hmmsearch_alignment.All.log"
> $HMM_LOG

# hmmsearch
echo "## hmmsearch & alignment | All Species (n = 1609) | Begin" | tee -a $HMM_LOG
while read -a i ; do
    bash ${DIR}/job_hmmsearch.sh ${i[2]} ${i[3]} ${i[5]} ${i[6]} ${i[4]} ${COMPILED}/Compiled.All.fasta N All ${ANALYSIS} $HMM_LOG
done < ${DOM_LIST}
echo "## hmmsearch & alignment | All Species (n = 1609) | Completed ✅" | tee -a $HMM_LOG

# ==============================================================================
# 5. Gene counts
# ==============================================================================
mkdir -p ${ANALYSIS}/gene_counts

# summarize gene counts into single file
while read -a i ; do
    cat ${ANALYSIS}/searches/All.${i[2]}.genes.list | sed "s/$/\t${i[2]}/"
done < ${DOM_LIST} > ${ANALYSIS}/gene_counts/1609.genecounts.csv
echo "## Gene counts | All Species (n = 1609) | Completed ✅"
