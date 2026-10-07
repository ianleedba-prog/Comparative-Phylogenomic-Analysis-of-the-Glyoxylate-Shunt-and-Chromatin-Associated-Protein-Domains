#!/bin/bash
# ==============================================================================
# Software:
    # bioawk (v1.0)
    # seqkit (v2.10.0)
    # mafft (v7.525)
    # clipkit (v1.4.1)
# ==============================================================================

# Directories
DIR="/Directory/to/data/JmjC"
ANALYSIS="${DIR}/analysis"
SEARCHES="${ANALYSIS}/searches"

# Define function
find_conserved_sites() {
    local GRP=$1
    local DOM=$2
    local REF=$([ "$GRP" == "PLT" ] && echo "HsapAtha" || echo "Hsap")
    
    local DIR_GRP="$ANALYSIS/orthologs/$REF/$DOM/$GRP"
    local DIR_CON="$DIR_GRP/conserved_sites"; mkdir -p $DIR_CON
    
    local ORT_UNIQ="$DIR_GRP/${GRP}.${DOM}.orthogroups.uniq.csv"
    local FASTA="$SEARCHES/Rep.${DOM}.seqs.fasta"
    local LOG_FILE="$DIR_CON/${GRP}.${DOM}.Conserved.sites.log"
    
    # Reset log file
    > $LOG_FILE
    
    # --------------------------------------------------------------------------
    # 1. Retrieve Domain Information & Expand Regions
    # --------------------------------------------------------------------------
    echo "## Conserved sites identification (ref: Human) | $GRP | $DOM | Retrieving domain info & Expand regions (-100 +50)" | tee -a $LOG_FILE
    # get domain ids
    cut -d$'\t' -f1 $ORT_UNIQ > $DIR_CON/id.tmp
    grep -F -f $DIR_CON/id.tmp $SEARCHES/Rep.${DOM}.domains.csv > $DIR_CON/${GRP}.${DOM}.domains.tmp
    
    # get length info
    > $DIR_CON/length.tmp
    while read -r id; do
        bioawk -c fastx -v id="$id" '$name == id {print $name, length($seq)}' $FASTA >> $DIR_CON/length.tmp
    done < <(cut -f1 -d ' ' $DIR_CON/${GRP}.${DOM}.domains.tmp)
    
    paste -d " " $DIR_CON/${GRP}.${DOM}.domains.tmp <(cut -f2 -d$'\t' $DIR_CON/length.tmp) > $DIR_CON/${GRP}.${DOM}.domains.tmp2
    
    # expand domain regions (-100, +50)
    awk -F' ' '{
        adjusted_start = ($2 - 100 > 0) ? $2 - 100 : 1
        adjusted_end = ($3 + 50 > $5) ? $5 : $3 + 50
        print $1, adjusted_start, adjusted_end
        }
    ' $DIR_CON/${GRP}.${DOM}.domains.tmp2 > $DIR_CON/${GRP}.${DOM}.domains.tmp3
    
    # --------------------------------------------------------------------------
    # 2. Organize domain info (+Homology Groups) & get fasta sequences
    # --------------------------------------------------------------------------
    echo "## Conserved sites identification (ref: Human) | $GRP | $DOM | Organize domain info & get fasta sequences" | tee -a $LOG_FILE
    local DOM_INFO=$DIR_CON/${GRP}.${DOM}.domains.csv # summary of domains (id, start (-100), end (+50), homology group)
    local DOM_FASTA=$DIR_CON/${GRP}.${DOM}.domains.fasta # fasta sequence
    
    # summary of domains
    awk '
        BEGIN {
            OFS = " ";
        }
        FNR == NR {
            n = split($2, arr, ":");
            result = "";
            for (i = 2; i <= n; i++) {
                if (i > 2) result = result ".";
                gsub("/", ".", arr[i]);
                result = result arr[i];
            }
            a[$1] = result;
            next;
        }
        {
            if ($1 in a) {
                print $1, $2, $3, a[$1];
            }
        }
    ' FS="\t" $ORT_UNIQ FS=" " $DIR_CON/${GRP}.${DOM}.domains.tmp3 > $DOM_INFO
    
    # get fasta
    esl-sfetch -C -f $FASTA <(awk '{ print $1"_"$2"-"$3, $2, $3, $1 }' $DOM_INFO) > $DOM_FASTA
    
    # --------------------------------------------------------------------------
    # 3. Find conserved sites using Human reference sequences
    # --------------------------------------------------------------------------
    
    # final output file
    local FIN_OUT=$DIR_CON/${GRP}.${DOM}.conserved.sites.csv
    > $FIN_OUT
    
    for HG in $(cut -d " " -f4 $DOM_INFO | sort | uniq -d); do
        echo "## Conserved sites identification (ref: Human) | $GRP | $DOM | Processing each HG | $HG | Preparing data" | tee -a $LOG_FILE
        
        local DOM_hFASTA=$DIR_CON/${GRP}.${DOM}.domains_${HG}.fasta
        local DOM_hFASTA_dedup=$DIR_CON/${GRP}.${DOM}.domains_${HG}.dedup.fasta
        local ALIGNED_DOM_FASTA=$DIR_CON/${GRP}.${DOM}.domains_${HG}.dedup.l.fasta
        local TRIMMED_DOM_FASTA=$DIR_CON/${GRP}.${DOM}.domains_${HG}.dedup.lt.fasta
        
        # a. Fetch domain sequences for current HG
        awk -v hg="$HG" '$4 == hg' "$DOM_INFO" > $DIR_CON/${GRP}.${DOM}.domains_${HG}.tmp
        esl-sfetch -C -f $FASTA <(awk '{ print $1"_"$2"-"$3, $2, $3, $1 }' $DIR_CON/${GRP}.${DOM}.domains_${HG}.tmp) > $DOM_hFASTA
        esl-sfetch --index $DOM_hFASTA >> $LOG_FILE 2>&1
        
        # b. Define REF_ID and binding sites
        if [[ "$HG" == *KDM2* || "$HG" == *KDM7* ]]; then
            REF_DOM="KDM2B"
            REF_ID=$(grep -F Hsap_NP_115979.3 "$DOM_INFO" | awk '{ print $1"_"$2"-"$3 }')
            binding_sites=(208 211 213 228 283)
        elif [[ "$HG" == *KDM3* ]]; then
            REF_DOM="KDM3B"
            REF_ID=$(grep -F Hsap_NP_057688.3 "$DOM_INFO" | awk '{ print $1"_"$2"-"$3 }')
            binding_sites=(52 68 71 73 79 200)
        elif [[ "$HG" == *KDM4* || "$HG" == *KDM5* ]]; then
            REF_DOM="KDM4A"
            REF_ID=$(grep -F Hsap_NP_055478.2 "$DOM_INFO" | awk '{ print $1"_"$2"-"$3 }')
            binding_sites=(67 123 125 133 141 176 211)
        elif [[ "$HG" == *KDM6* ]]; then
            REF_DOM="KDM6A"
            REF_ID=$(grep -F Hsap_NP_001406738.1 "$DOM_INFO" | awk '{ print $1"_"$2"-"$3 }')
            binding_sites=(114 120 123 125 131 133 203)
        elif [[ "$HG" == *KDM8* || "$HG" == *JMJD6* || "$HG" == *JMJD7* ]]; then
            REF_DOM="KDM8"
            REF_ID=$(grep -F Hsap_NP_001138820.1 "$DOM_INFO" | awk '{ print $1"_"$2"-"$3 }')
            binding_sites=(187 233 236 238 242 251 315 329)
        fi
        
        # c. Add Human reference sequence to each fasta & Deduplicate (seqkit)
        seqkit grep -p $REF_ID $DOM_FASTA > $DIR_CON/${REF_DOM}.tmp.fasta
        cat $DIR_CON/${REF_DOM}.tmp.fasta $DOM_hFASTA > ${DOM_hFASTA/.fasta/.tmp.fasta}
        seqkit rmdup -n ${DOM_hFASTA/.fasta/.tmp.fasta} -o $DOM_hFASTA_dedup >> $LOG_FILE 2>&1
        
        # d. MSA + Trim
        echo "## Conserved sites identification (ref: Human) | $GRP | $DOM | Processing each HG | $HG | Multiple sequence alignment" | tee -a $LOG_FILE
        mafft --quiet --genafpair --thread 8 --reorder --maxiterate 1000 $DOM_hFASTA_dedup > $ALIGNED_DOM_FASTA 2>> $LOG_FILE
        clipkit $ALIGNED_DOM_FASTA -m kpic-gappy -o $TRIMMED_DOM_FASTA -g 0.7 >> $LOG_FILE 2>&1
        
        # e. Conserved sites position in aligned reference sequence
        echo "## Conserved sites identification (ref: Human) | $GRP | $DOM | Processing each HG | $HG | Sites Identification" | tee -a $LOG_FILE
        local REF_SEQ=$(awk -v id="$REF_ID" '
            $0 ~ id {flag=1; next}
            /^>/ && flag {exit}
            flag {seq = seq $0}
            END {gsub(/[ \t]/, "", seq); print seq}' $ALIGNED_DOM_FASTA)
            
        local position_count=0
        local msa_positions=()

        for ((i = 0; i < ${#REF_SEQ}; i++)); do
            local char=${REF_SEQ:i:1}
            if [[ "$char" != "-" ]]; then
                position_count=$((position_count + 1))
            fi
            
            for site in "${binding_sites[@]}"; do
                if [[ "$position_count" -eq "$site" && "$char" != "-" ]]; then
                    if [[ ! " ${msa_positions[@]} " =~ " $((i + 1)) " ]]; then
                        msa_positions+=($((i + 1)))
                    fi
                fi
            done
        done
        
        local positions_str=$(IFS=" "; echo "${msa_positions[*]}") # conserved sites profile
        
        # f. Extract conserved sites for query domains
        local HG_OUT=$DIR_CON/${GRP}.${DOM}.conserved.sites_${HG}.csv # output per homology group
        > $HG_OUT
        
        awk -v id="$REF_ID" -v positions="$positions_str" -v msa_positions="${msa_positions[*]}" -v output="$HG_OUT" '
        BEGIN {
            split(positions, pos_arr, " ")
            split(msa_positions, msa_arr, " ")
            seq_id = ""; seq = ""; OFS = "\t"
        }
        /^>/ {
            if (seq_id != "") {
                binding_residues = ""; msa_location = ""
                for (p in pos_arr) {
                    binding_residues = binding_residues substr(seq, pos_arr[p], 1) " "
                    msa_location = msa_location msa_arr[p] " "
                }
                print seq_id, binding_residues, msa_location >> output
            }
            seq_id = $1; gsub(/^>/, "", seq_id); seq = ""; next
        }
        { gsub(/[ \t]/, "", $0); seq = seq $0 }
        END {
            binding_residues = ""; msa_location = ""
            for (p in pos_arr) {
                binding_residues = binding_residues substr(seq, pos_arr[p], 1) " "
                msa_location = msa_location msa_arr[p] " "
            }
            if (seq_id != "") print seq_id, binding_residues, msa_location >> output
        }' $ALIGNED_DOM_FASTA
        
        # g. Record into final output file
        if [[ "$HG" == *${REF_DOM}* ]]; then
            cut -d$'\t' -f1,2 $HG_OUT >> $FIN_OUT
        else
            grep -v "^${REF_ID}" $HG_OUT | cut -d$'\t' -f1,2 >> $FIN_OUT
        fi
        
        rm $HG_OUT $DOM_hFASTA ${DOM_hFASTA}.ssi
        
    done
    
    # --------------------------------------------------------------------------
    # 4. Sort final output & remove tmp files
    # --------------------------------------------------------------------------
    sort $FIN_OUT > ${FIN_OUT/.csv/.sorted.csv}
    
    rm $DIR_CON/*.tmp* $FIN_OUT
    
    echo "## Conserved sites identification (ref: Human) | $GRP | $DOM | Processing each HG | $HG | Sites Identification" | tee -a $LOG_FILE

}

find_conserved_sites $1 $2
