#!/bin/bash
# ==============================================================================
# Software:
    # ncbi-datasets-cli (v16.40.1)
# ==============================================================================

# Directories
DIR="/Directory/to/data"
FASTA="$DIR/protein_fa"; mkdir -p $FASTA
ENS_SOURCE="$DIR/ensembl_source" # folder where species_Ensembl(Group).txt (downloaded from Ensembl FTP sites) is located.

# Species metadata
METADATA="${DIR}/Metadata_species_info.txt"

# Log file
LOG_FILE="$DIR/s01a_download_fasta.sh.log"
> $LOG_FILE

# Failed Log
FAILED_LOG="$DIR/s01a.failed.list.txt"
echo -e "Tax_name\tSource\tFile_name\tURL" > $FAILED_LOG

# Multi jobs
MAX_JOBS=6

# Define function
download_and_process() {
    local tax_name=$1
    local abbr=$2
    local ncbi_id=$3
    local file_source=$4
    local file_name=$5

    # ==========================================
    # 1. Source = NCBI
    # ==========================================
    if [[ $file_source == NCBI* ]]; then
        zip_file=$FASTA/${ncbi_id}.zip
        tmp_dir=$FASTA/tmp_${ncbi_id}
        
        if datasets download genome accession "$file_name" --include protein --filename "$zip_file" >/dev/null 2>&1; then
            unzip -q $zip_file -d $tmp_dir
            target_path=$(grep ${file_name}/protein.faa $tmp_dir/md5sum.txt | awk '{print $2}')

            if [[ -n $target_path && -f $tmp_dir/$target_path ]]; then
                mv $tmp_dir/$target_path $FASTA/${file_name}.faa
                echo "## wget: FASTA | $file_source | $tax_name | DONE ✅" | tee -a $LOG_FILE
                sed "s/^>/>${abbr}_/" $FASTA/${file_name}.faa | gzip > $FASTA/${file_name}.faa.gz && rm -f $FASTA/${file_name}.faa
                echo "## wget: FASTA | $file_source | $tax_name | Added species prefix & Compressed to gz | DONE ✅" | tee -a $LOG_FILE
            else
                echo "## wget: FASTA | $file_source | $tax_name | Downloaded success, but no file/directory inside-FAILED ❌" | tee -a $LOG_FILE
                echo -e "${tax_name}\t${file_source}\t${file_name}\tdatasets_download_success_but_no_file_inside" >> $FAILED_LOG
            fi
            # clean
            rm -rf $tmp_dir $zip_file
        else
            echo "## wget: FASTA | $file_source | $tax_name | FAILED ❌" | tee -a $LOG_FILE
            echo -e "${tax_name}\t${file_source}\t${file_name}\tdatasets_download_fail" >> $FAILED_LOG
            rm -rf $tmp_dir $zip_file
        fi

    # ==========================================
    # 2. Source = WormBase
    # ==========================================
    elif [[ $file_source == WormBase* ]]; then
        
        release_version=$(echo $file_source | cut -d'_' -f2)
        species_name=$(echo $file_name | cut -d'.' -f1)
        assembly_id=$(echo $file_name | cut -d'.' -f2)
        
        # define url
        url="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS${release_version}/species/${species_name}/${assembly_id}/${file_name}.gz"
        
        # download data
        if wget -qO- $url | zcat 2>/dev/null > $FASTA/${file_name}; then
            echo "## wget: FASTA | $file_source | $tax_name | DONE ✅" | tee -a $LOG_FILE
            sed "s/^>/>${abbr}_/" $FASTA/$file_name | gzip > $FASTA/${file_name}.gz && rm -f $FASTA/$file_name
            echo "## wget: FASTA | $file_source | $tax_name | Added species prefix & Compressed to gz | DONE ✅" | tee -a $LOG_FILE
        else
            echo "## wget: FASTA | $file_source | $tax_name | FAILED ❌" | tee -a $LOG_FILE
            echo -e "${tax_name}\t${file_source}\t${file_name}\t${url}" >> $FAILED_LOG
            rm -f $FASTA/${file_name}
        fi

    # ==========================================
    # 3. Source = Ensembl
    # ==========================================
    elif [[ $file_source == Ensembl* ]]; then

        col2=$(echo $file_source | cut -d'_' -f2)
        col3=$(echo $file_source | cut -d'_' -f3)
        
        if [[ -z $col3 ]]; then
            tax_group="vertebrates" # case for Ensembl_115
            release_version=$col2
        else
            tax_group=$col2         # case for Ensembl_(metazoa/fungi/protists/plants)_61
            release_version=$col3
        fi

        ens_source_file=${ENS_SOURCE}/species_Ensembl${tax_group}.txt # downloaded from Ensembl FTP sites

        # get species name from file name, then convert all letters to lower case
        species_name=$(echo $file_name | cut -d'.' -f1)
        species_name=$(echo $species_name | tr '[:upper:]' '[:lower:]')
        
        # match species name within ens_source_file, then extract collection folder name
        core_db=$(awk -F'\t' -v sp="$species_name" '$2 == sp {print $14; exit}' $ens_source_file)
        collection_name=$(echo $core_db | sed 's/_core.*//')

        # define url
        if [[ ${tax_group} != "vertebrates" ]]; then    # case for Ensembl_(metazoa/fungi/protists/plants)_61
            url1="http://ftp.ensemblgenomes.org/pub/${tax_group}/release-${release_version}/fasta/${collection_name}/${species_name}/pep/${file_name}.gz" # Plan A
            url2="http://ftp.ensemblgenomes.org/pub/${tax_group}/release-${release_version}/fasta/${collection_name}/pep/${file_name}.gz" # Plan B
        else                                            # case for Ensembl_115
            url1="http://ftp.ensembl.org/pub/release-${release_version}/fasta/${species_name}/pep/${file_name}.gz" # Plan A
            url2="http://ftp.ensembl.org/pub/release-${release_version}/fasta/pep/${file_name}.gz" # Plan B
        fi

        # download data ( Plan A -> if url failed -> Plan B )
        if wget -qO- $url1 | zcat 2>/dev/null > "$FASTA/${file_name}"; then # Plan A
            echo "## wget: FASTA | $file_source | $tax_name | DONE ✅" | tee -a $LOG_FILE
            sed "s/^>/>${abbr}_/" $FASTA/$file_name | gzip > $FASTA/${file_name}.gz && rm -f $FASTA/$file_name
            echo "## wget: FASTA | $file_source | $tax_name | Added species prefix & Compressed to gz | DONE ✅" | tee -a $LOG_FILE
            
        else
            rm -f $FASTA/${file_name} # clean file
            if wget -qO- $url2 | zcat 2>/dev/null > $FASTA/${file_name}; then # Plan B
                echo "## wget: FASTA | $file_source | $tax_name | DONE ✅" | tee -a $LOG_FILE
                sed "s/^>/>${abbr}_/" $FASTA/$file_name | gzip > $FASTA/${file_name}.gz && rm -f $FASTA/$file_name
                echo "## wget: FASTA | $file_source | $tax_name | Added species prefix & Compressed to gz | DONE ✅" | tee -a $LOG_FILE
                
            else
                # if both urls failed, record into failed log file
                echo "## wget: FASTA | $file_source | $tax_name | FAILED ❌" | tee -a $LOG_FILE
                echo -e "${tax_name}\t${file_source}\t${file_name}\t${url1} (or ${url2})" >> $FAILED_LOG
                rm -f $FASTA/${file_name}
            fi
        fi
    fi
}

# Download data (First trial)
echo "## wget data | Starting parallel downloads (Max concurrent jobs: $MAX_JOBS)" | tee -a $LOG_FILE
while IFS=$'\t' read -r tax_name abbr ncbi_id file_source file_name; do
    
    ( download_and_process "$tax_name" "$abbr" "$ncbi_id" "$file_source" "$file_name" ) &
    
    while [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; do
        sleep 1
    done
done < <(tail -n +2 $METADATA)

wait

# Download data (Second trial if there are failed cases due to traffic issues)
if [ -s $FAILED_LOG ] && [ $(wc -l < $FAILED_LOG) -gt 1 ]; then
echo "## wget: FASTA | Download failure detected: ** Traffic issues **. Retrial one-by-one for those failure." | tee -a $LOG_FILE
    
    mv $FAILED_LOG ${FAILED_LOG}.bak # rename failed log as backup
    echo -e "Tax_name\tSource\tFile_name\tURL" > $FAILED_LOG # recreate failed log file
    
    # create tmp_metadata for failed list
    TMP_METADATA="$DIR/Tmp_failed_metadata.txt"
    awk -F'\t' 'NR==FNR {if(NR>1) failed[$1]; next} $1 in failed' ${FAILED_LOG}.bak $METADATA > $TMP_METADATA
    
    while IFS=$'\t' read -r tax_name abbr ncbi_id file_source file_name; do
        echo "## wget: FASTA | $file_source | $tax_name | RETRIAL 🔁" | tee -a $LOG_FILE
        download_and_process "$tax_name" "$abbr" "$ncbi_id" "$file_source" "$file_name"
    done < $TMP_METADATA
fi

# WARNING if there still is failure.
if [ -s $FAILED_LOG ] && [ $(wc -l < $FAILED_LOG) -gt 1 ]; then
    echo "## wget: FASTA | Download failure unresolved. Please check s01a.failed.list.txt." | tee -a $LOG_FILE
else
    echo "## wget: FASTA | Download failure all resolved." | tee -a $LOG_FILE
    echo "## wget: FASTA | All downloads successfully completed (NCBI, WormBase, Ensembl) ✅" | tee -a $LOG_FILE
    echo "## wget: FASTA | Proceed to those derived from ** other publications **. Please download manually according to the Metadata, then save as .gz." | tee -a $LOG_FILE
fi

# remove backup, tmp files
rm -f ${FAILED_LOG}.bak $TMP_METADATA

# remove failed list if no record
if [ -s $FAILED_LOG ] && [ $(wc -l < $FAILED_LOG) -eq 1 ]; then
    rm $FAILED_LOG
fi

## Please download data from ** other publication **. Store in protein_fa/, prepend species abbreviation to all sequence headers ##
# PUB_METADATA="$DIR/publication_metadata.txt"
# awk 'NR==1 || $2=="apple"' $METADATA > $PUB_METADATA
# while IFS=$'\t' read -r tax_name abbr ncbi_id file_source file_name; do
    # sed "s/^>/>${abbr}_/" $FASTA/$file_name | gzip > $FASTA/${file_name}.gz && rm -f $FASTA/$file_name
# done < $TMP_METADATA
