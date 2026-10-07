#!/bin/bash
# ==============================================================================
# Software:
    # ncbi-datasets-cli (v16.40.1)
# ==============================================================================

# Directories
DIR="/Directory/to/data"
ANN="$DIR/annotation"; mkdir -p $ANN
ENS_SOURCE="$DIR/ensembl_source" # folder where species_Ensembl(Group).txt (downloaded from Ensembl FTP sites) is located.

# Representative species metadata
METADATA="${DIR}/Metadata_rep.species_info.txt"

# Log file
LOG_FILE="$DIR/s01b_download_annotation.sh.log"
> $LOG_FILE

# Failed Log
FAILED_LOG="$DIR/s01b.failed.list.txt"
echo -e "Tax_name\tSource\tAnnotation\tURL" > $FAILED_LOG

# Multi jobs
MAX_JOBS=6

# Define function
download_and_process() {
    local tax_name=$1
    local abbr=$2
    local ncbi_id=$3
    local file_source=$4
    local file_name=$5
    local file_annotation=$6

    # ==========================================
    # 1. Source = NCBI
    # ==========================================
    if [[ $file_source == NCBI* ]]; then
        zip_file=$ANN/${ncbi_id}.zip
        tmp_dir=$ANN/tmp_${ncbi_id}
        
        if datasets download genome accession "$file_name" --include gtf --filename "$zip_file" >/dev/null 2>&1; then
            unzip -q $zip_file -d $tmp_dir
            target_path=$(grep ${file_name}/genomic.gtf $tmp_dir/md5sum.txt | awk '{print $2}')

            if [[ -n $target_path && -f $tmp_dir/$target_path ]]; then
                mv $tmp_dir/$target_path $ANN/${file_name}.gtf
                echo "## wget: Annotation File | $file_source | $tax_name | DONE ✅" | tee -a $LOG_FILE
            else
                echo "## wget: Annotation File | $file_source | $tax_name | Downloaded success, but no file/directory inside-FAILED ❌" | tee -a $LOG_FILE
                echo -e "${tax_name}\t${file_source}\t${file_annotation}\tdatasets_download_success_but_no_file_inside" >> $FAILED_LOG
            fi
            # clean
            rm -rf $tmp_dir $zip_file
        else
            echo "## wget: Annotation File | $file_source | $tax_name | FAILED ❌" | tee -a $LOG_FILE
            echo -e "${tax_name}\t${file_source}\t${file_annotation}\tdatasets_download_fail" >> $FAILED_LOG
            rm -rf $tmp_dir $zip_file
        fi

    # ==========================================
    # 2. Source = WormBase
    # ==========================================
    elif [[ $file_source == WormBase* ]]; then
        
        release_version=$(echo $file_source | cut -d'_' -f2)
        species_name=$(echo $file_annotation | cut -d'.' -f1)
        assembly_id=$(echo $file_annotation | cut -d'.' -f2)
        
        # define url
        url="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS${release_version}/species/${species_name}/${assembly_id}/${file_annotation}.gz"
        
        # download data
        if wget -qO- $url | zcat 2>/dev/null > $ANN/${file_annotation}; then
            echo "## wget: Annotation File | $file_source | $tax_name | DONE ✅" | tee -a $LOG_FILE
        else
            echo "## wget: Annotation File | $file_source | $tax_name | FAILED ❌" | tee -a $LOG_FILE
            echo -e "${tax_name}\t${file_source}\t${file_annotation}\t${url}" >> $FAILED_LOG
            rm -f $ANN/${file_annotation}
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
        species_name=$(echo $file_annotation | cut -d'.' -f1)
        species_name=$(echo $species_name | tr '[:upper:]' '[:lower:]')
        
        # match species name within ens_source_file, then extract collection folder name
        core_db=$(awk -F'\t' -v sp="$species_name" '$2 == sp {print $14; exit}' $ens_source_file)
        collection_name=$(echo $core_db | sed 's/_core.*//')

        # define url
        if [[ ${tax_group} != "vertebrates" ]]; then    # case for Ensembl_(metazoa/fungi/protists/plants)_61
            url1="http://ftp.ensemblgenomes.org/pub/${tax_group}/release-${release_version}/gtf/${collection_name}/${species_name}/${file_annotation}.gz" # Plan A
            url2="http://ftp.ensemblgenomes.org/pub/${tax_group}/release-${release_version}/gtf/${collection_name}/${file_annotation}.gz" # Plan B
        else                                            # case for Ensembl_115
            url1="http://ftp.ensembl.org/pub/release-${release_version}/gtf/${species_name}/${file_annotation}.gz" # Plan A
            url2="http://ftp.ensembl.org/pub/release-${release_version}/gtf/${file_annotation}.gz" # Plan B
        fi

        # download data ( Plan A -> if url failed -> Plan B )
        if wget -qO- $url1 | zcat 2>/dev/null > "$ANN/${file_annotation}"; then # Plan A
            echo "## wget: Annotation File | $file_source | $tax_name | DONE ✅" | tee -a $LOG_FILE
        else
            rm -f $ANN/${file_annotation} # clean file
            if wget -qO- $url2 | zcat 2>/dev/null > $ANN/${file_annotation}; then # Plan B
                echo "## wget: Annotation File | $file_source | $tax_name | DONE ✅" | tee -a $LOG_FILE
            else
                # if both urls failed, record into failed log file
                echo "## wget: Annotation File | $file_source | $tax_name | FAILED ❌" | tee -a $LOG_FILE
                echo -e "${tax_name}\t${file_source}\t${file_annotation}\t${url1} (or ${url2})" >> $FAILED_LOG
                rm -f $ANN/${file_annotation}
            fi
        fi
    fi
}

# Download data (First trial)
echo "## wget data | Starting parallel downloads (Max concurrent jobs: $MAX_JOBS)" | tee -a $LOG_FILE
while IFS=$'\t' read -r tax_name abbr ncbi_id file_source file_name file_annotation; do
    
    ( download_and_process "$tax_name" "$abbr" "$ncbi_id" "$file_source" "$file_name" "$file_annotation" ) &
    
    while [[ $(jobs -r -p | wc -l) -ge $MAX_JOBS ]]; do
        sleep 1
    done
done < <(tail -n +2 $METADATA)

wait

# Download data (Second trial if there are failed cases due to traffic issues)
if [ -s $FAILED_LOG ] && [ $(wc -l < $FAILED_LOG) -gt 1 ]; then
echo "## wget: Annotation File | Download failure detected: ** Traffic issues **. Retrial one-by-one for those failure." | tee -a $LOG_FILE
    
    mv $FAILED_LOG ${FAILED_LOG}.bak # rename failed log as backup
    echo -e "Tax_name\tSource\tAnnotation\tURL" > $FAILED_LOG # recreate failed log file
    
    # create tmp_metadata for failed list
    TMP_METADATA="$DIR/Tmp_failed_metadata.txt"
    awk -F'\t' 'NR==FNR {if(NR>1) failed[$1]; next} $1 in failed' ${FAILED_LOG}.bak $METADATA > $TMP_METADATA
    
    while IFS=$'\t' read -r tax_name abbr ncbi_id file_source file_name file_annotation; do
        echo "## wget: Annotation File | $file_source | $tax_name | RETRIAL 🔁" | tee -a $LOG_FILE
        download_and_process "$tax_name" "$abbr" "$ncbi_id" "$file_source" "$file_name" "$file_annotation"
    done < $TMP_METADATA
fi

# WARNING if there still is failure.
if [ -s $FAILED_LOG ] && [ $(wc -l < $FAILED_LOG) -gt 1 ]; then
    echo "## wget: Annotation File | Download failure unresolved. Please check s01b.failed.list.txt." | tee -a $LOG_FILE
else
    echo "## wget: Annotation File | Download failure all resolved." | tee -a $LOG_FILE
    echo "## wget: Annotation File | All downloads complete (NCBI, WormBase, Ensembl) ✅" | tee -a $LOG_FILE
    echo "## wget: Annotation File | Proceed to those derived from ** other publications **. Please download manually according to the Metadata, then save UNZIPPED." | tee -a $LOG_FILE
fi

# remove backup, tmp files
rm -f ${FAILED_LOG}.bak $TMP_METADATA

# remove failed list if no record
if [ -s $FAILED_LOG ] && [ $(wc -l < $FAILED_LOG) -eq 1 ]; then
    rm $FAILED_LOG
fi
