#!/bin/bash
# ==============================================================================
# Software:
    # busco (v6.0.0)
    # GNU parallel
# ==============================================================================

# Directories
DIR="/Directory/to/data"
FASTA="${DIR}/protein_fa"
BUSCO="${DIR}/busco"

# remove existing files if restart is necessary
rm -rf $BUSCO; mkdir -p $BUSCO

# Species metadata
METADATA="${DIR}/Metadata_species_info.txt"

# Run analyses
# 1. BUSCO analysis
parallel --colsep '\t' --header : -j 8 \
    "cd ${BUSCO}; \
     matched_file=\$(find ${FASTA} -maxdepth 1 -name \"*{File}*.gz\" | head -n 1); \
     busco -i \"\$matched_file\" \
           -o {Abbreviation}_busco \
           -m prot \
           -l ${DIR}/lineages/eukaryota_odb10 \
           --cpu 1 \
           --out_path ${BUSCO} \
           --download_path ${DIR} \
           --offline \
           --quiet > /dev/null 2>&1;" \
    < $METADATA
    
rm -rf ${BUSCO}/busco_downloads ${BUSCO}/busco_*.log ${DIR}/file_versions.tsv # clean

# 3. Summary
echo -e "Sps_abb\tComplete\tFragmented\tMissing\tComplete_percentage" > ${BUSCO}/summary_busco.txt

find ${BUSCO} -name "short_summary.*.json" | sort | while read json_file; do
    # Extract species abbreviation
    parent_dir=$(dirname "$json_file")
    folder_name=$(basename "$parent_dir")
    abbr="${folder_name%_busco}"
    
    # Complete, Fragmented, Missing BUSCOs, Complete percentage
    C=$(grep '"Complete BUSCOs":' "$json_file" | tr -dc '0-9')
    F=$(grep '"Fragmented BUSCOs":' "$json_file" | tr -dc '0-9')
    M=$(grep '"Missing BUSCOs":' "$json_file" | tr -dc '0-9')
    Comp_Pct=$(grep '"Complete percentage":' "$json_file" | tr -dc '0-9.')

    # Record in summary file
    echo -e "${abbr}\t${C}\t${F}\t${M}\t${Comp_Pct}" >> ${BUSCO}/summary_busco.txt
done

# 4. Add busco results in Metadata
awk -F'\t' 'BEGIN {OFS="\t"}
    # 1. summary_busco.txt
    NR == FNR {
        if (FNR == 1) {
            header_name = $5
        } else {
            # species abbreviation, Complete percentage
            val_map[$1] = $5
        }
        next
    }
    # 2. Metadata
    FNR == 1 {
        # Add "Complete_percentage" to header of Metadata
        $(NF+1) = header_name
        print
        next
    }
    {
        if ($2 in val_map) {
            $(NF+1) = val_map[$2]
        } else {
            $(NF+1) = "NA"
        }
        print
    }
' ${BUSCO}/summary_busco.txt $METADATA | cut -f1,3,6 > ${DIR}/busco_scores.txt # Recorded in Source_data Supple_Fig.2A
