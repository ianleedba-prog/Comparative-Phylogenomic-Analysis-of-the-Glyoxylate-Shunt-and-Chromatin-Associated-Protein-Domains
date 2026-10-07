#!/bin/bash
# ==============================================================================
# Software:
    # bioawk (v1.0)
# ==============================================================================

# Directories
DIR="/Directory/to/data"
ANALYSIS="${DIR}/analysis"
DIR_ANNO="${DIR}/annotation"

# Define functions
process_isoforms() {
    # local variables
    local GRP=$1; local DOM=$2; local SPS=$3; local MODE=$4
    local REF=$([ "$GRP" == "PLT" ] && echo "HsapAtha" || echo "Hsap")
    
    local DIR_ORT="$ANALYSIS/orthologs/$REF/$DOM"
    local DIR_GRP="$DIR_ORT/$GRP"; mkdir -p $DIR_GRP
    
    local ORT_LIST="$DIR_ORT/$DOM.$REF.orthogroups.filt.csv"
    local ISO_LIST="$DIR_GRP/${GRP}.${DOM}.isoforms.${MODE}.csv"
    local EXT=$([ "$MODE" == "gtf" ] && echo "gtf" || echo "gff*")
    local PATTERN1=$([ "$MODE" == "gtf" ] && echo '(protein_id) "([^"]+)"' || echo '(Parent|Name)=([^;]+)')
    local PATTERN2=$([ "$MODE" == "gtf" ] && echo '(gene_id) "([^"]+)"' || echo '(ID|gene|Name)=([^;]+)')

    # species list
    local SPS_LIST=$(ls $DIR_ANNO/*.$EXT | xargs -n 1 basename | cut -d '_' -f1 | grep -Fxf "$SPS")
    
    # skip function when gff mode is not needed
    if [[ "$MODE" == "gff" && -z "$SPS_LIST" ]]; then
        return 0
    fi

    # find gene id for each ortholog from annotation file
    > $ISO_LIST
    while read -r species; do
        ANN=$(ls $DIR_ANNO/${species}_*.$EXT | head -n 1)
        [ -z "$ANN" ] && continue

        gawk -v mode="$MODE" -v pat1="$PATTERN1" -v pat2="$PATTERN2" '
        BEGIN { FS=OFS="\t" }
        
        # 1. load ortholog ids
        NR == FNR { 
            full_id = $1
            clean_id = full_id; sub(/\|.*$/, "", clean_id)

            # Case for Twil gff
            if (mode == "gff" && match(clean_id, /c[0-9]+_g[0-9]+_i[0-9]+/)) {
                core_key = substr(clean_id, RSTART, RLENGTH)
            } 
            # Case for Pigvie, Pigchi gff
            else if (mode == "gff" && match(clean_id, /_g[0-9]+/)) {
                core_key = substr(clean_id, RSTART + 1)
            } 
            # Case for the rest of annotation files
            else {
                core_key = substr(clean_id, index(clean_id, "_") + 1)
            }
            
            query[core_key] = full_id
            next 
        }
        
        # save pos for each gene listed in annotation file
        ARGIND == 2 && $3 == "gene" {
            n = split($9, words, /[;=" ]+/)
            for (i=1; i<=n; i++) {
                if (words[i] != "") {
                    g_key = words[i]
                    gsub(/[.](path|mrna|t|cds|exon)[0-9]*$/, "", g_key)
                    gene_coords[g_key] = $4 "_" $5
                }
            }
            next
        }
        
        # 3. match ortholog ids with protein id listed in annotation file
        ARGIND == 3 {
            is_target = 0
            if (has_cds == 1) {
                if ($3 == "CDS") is_target = 1
            } 
            else {
                if ($3 !~ /^(gene|exon|transcript|start_codon|stop_codon)$/) is_target = 1
            }

            if (is_target == 1) {
                if (match($9, pat1, arr1)) {
                    p_id = arr1[2]; if (p_id == "") next

                    if (mode == "gff") {
                        gsub(/^cds[.]|[.](path|mrna|t|cds|exon)[0-9]*$/, "", p_id)
                    }
                    
                    if (p_id in query) {
                        if (match($9, pat2, arr2)) {
                            g_id = arr2[2]
                            
                            if (mode == "gff") {
                                gsub(/^cds[.]|[.](path|mrna|t|cds|exon)[0-9]*$/, "", g_id)
                            }
                            
                            if (!(g_id in gene_coords)) { sub(/[.][^.]+$/, "", g_id) }

                            if (g_id in gene_coords) {
                                print query[p_id], g_id, gene_coords[g_id]
                                delete query[p_id]
                            }                    
                        }
                    }
                }
            }
        }
        ' <(gawk -v sp="${species}_" '$1 ~ sp {print $1}' $ORT_LIST) "$ANN" "$ANN" >> $ISO_LIST
    done <<< "$SPS_LIST"

    # add protein length
    local FASTA="$ANALYSIS/searches/Rep.${DOM}.seqs.fasta"
    local ISO_LENGTH="$DIR_GRP/${GRP}.${DOM}.isoforms.unshrinked.${MODE}.csv"

    bioawk -c fastx '{print $name, length($seq)}' $FASTA | \
    awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR {len[$1]=$2; next} ($1 in len) {print $0, len[$1]}' - $ISO_LIST > $ISO_LENGTH

    # shrink isoforms (maintain the longest isoform; best priority given to the id that starts with "NP_")
    local ISO_SHRINK="$DIR_GRP/${GRP}.${DOM}.isoforms.shrinked.${MODE}.csv"
    gawk -F'\t' '
        {
            # variants
            oid=$1; gid=$2; pos=$3; slen=$4;
            split(oid, a, "_"); sp=a[1]; key=sp"_"gid;
            
            # priority given to "NP_"
            is_np=(oid ~ /_NP_/) ? 1 : 0;
            score=(is_np * 1000000) + slen;
            
            # record the prioritized isoform
            if (!(key in b_score) || score > b_score[key]) { b_score[key]=score; b_rec[key]=$0 }
        }
        END {
            for (k in b_rec) {
            print b_rec[k];
            }
        }' $ISO_LENGTH | sort -k1,1 > $ISO_SHRINK
    
    rm $ISO_LIST
}

# Define function #2
reformat_data() {
    # local variables
    local GRP=$1; local DOM=$2; local SPS=$3
    local REF=$([ "$GRP" == "PLT" ] && echo "HsapAtha" || echo "Hsap")
    local DIR_GRP="${DIR}/analysis/orthologs/$REF/$DOM/$GRP"
    
    local ISO_GTF="$DIR_GRP/${GRP}.${DOM}.isoforms.shrinked.gtf.csv"
    local ISO_GFF="$DIR_GRP/${GRP}.${DOM}.isoforms.shrinked.gff.csv"
    local ISO_FIN="$DIR_GRP/${GRP}.${DOM}.isoforms.shrinked.total.csv"

    # concatenate isoform lists (from gtf, gff annotations)
    cat $ISO_GTF $ISO_GFF 2>/dev/null > $ISO_FIN
    [ ! -s $ISO_FIN ] && cp $ISO_GTF $ISO_FIN # if there is no file generated from gff annotation, use only that from gtf annotation
    
    # fill gaps for species lacking annotation files
    local PROCESSED_SPS=$(cut -f1 "$ISO_FIN" 2>/dev/null | cut -d'_' -f1 | sort -u)
    local MISSING_SPS=$(grep -vxf <(echo "$PROCESSED_SPS") "$SPS")

    if [ ! -z "$MISSING_SPS" ]; then
        while read -r ms; do
            gawk -v sp="${ms}_" -F'\t' 'BEGIN{OFS="\t"} $1 ~ sp {print $1, "NA", "NA", "0"}' \
                "$ANALYSIS/orthologs/$REF/$DOM/$DOM.$REF.orthogroups.filt.csv" >> "$ISO_FIN"
        done <<< "$MISSING_SPS"
    fi

    sort -u -k1,1 "$ISO_FIN" -o "$ISO_FIN"

    # reformat orthogroup results
    local ORT_LIST="$ANALYSIS/orthologs/$REF/$DOM/$DOM.$REF.orthogroups.filt.csv"
    local ORT_UNIQ="$DIR_GRP/${GRP}.${DOM}.orthogroups.uniq.csv"
    awk -F'\t' 'NR==FNR {b[$1]; next} ($1 in b)' $ISO_FIN $ORT_LIST > $ORT_UNIQ

    # reformat fasta
    local FASTA="$ANALYSIS/searches/Rep.${DOM}.seqs.fasta"
    local FASTA_UNIQ="$DIR_GRP/${GRP}.${DOM}.seqs.fasta"
    
    > $FASTA_UNIQ
    while read id; do bioawk -c fastx -v id="$id" '$name == id {print ">"$name"\n"$seq}' $FASTA >> $FASTA_UNIQ ; done < <(cut -f1 $ORT_UNIQ)
}

# Run analyses
process_isoforms $1 $2 $3 "gtf"
process_isoforms $1 $2 $3 "gff"
reformat_data $1 $2 $3
