#!/bin/bash
# ==============================================================================
# Software:
    # BayesTraitsV4.1.3
# ==============================================================================

# Directories
DIR="/directory/to/data"
BAY="$DIR/Bayesian_traits"
PRM="$BAY/parameters"
RES="$BAY/results"

PRG="/directory/to/BayesTraitsV4.1.3"

# Tables written by the R scripts; s06 has to run after both of them
## RTAB="$DIR/analysis/results/Binary_table_eight_characters_259.csv"   # s04_Dstatistic.R
## RCNT="$DIR/analysis/results/Count_matrix_by_species_259.csv"         # s03_gene_counts.R

# Prepare input data
## prepare_table() {
##
##     local tax="$1"
##     local out="$BAY/binary.table.GS_TET.species_${tax}.txt"
##
##         awk -F, -v tax="$tax" -v OFS='\t' '
##         function trim(s) { sub(/\r$/, "", s); return s }
##
##         function key(s) { sub(/\r$/, "", s); gsub(/[\047"]/, "", s)
##                           gsub(/_/, " ", s); return s }
##
##         ## the count matrix, read first: species -> kingdom
##         NR == FNR { if (FNR > 1) kingdom[key($2)] = trim($3); next }
##
##         FNR == 1 {
##             off = ($1 == "" ? 0 : 1)
##             for (i = 1; i <= NF; i++) col[trim($i)] = i + off
##             if (!("GS" in col) || !("TET" in col)) {
##                 print "GS or TET is missing from " FILENAME > "/dev/stderr"
##                 exit 1
##             }
##             next
##         }
##
##         {
##             sp = key($1)
##             k  = kingdom[sp]
##             if (k == "") { print "no kingdom for " sp > "/dev/stderr"; exit 1 }
##
##             if (tax == "Metazoa") { if (k != "Metazoa") next }
##             else if (k != "Metazoa" && k != "Premetazoa") next
##
##             gsub(/ /, "_", sp)
##             print sp, trim($(col["GS"])), trim($(col["TET"]))
##         }' "$RCNT" "$RTAB" > "$out"
##
##     printf '%-8s %3d species | GS+ %3d | TET+ %3d\n' "$tax" \
##         "$(wc -l < "$out")" "$(awk -F'\t' '$2 == 1' "$out" | wc -l)" \
##         "$(awk -F'\t' '$3 == 1' "$out" | wc -l)"
##
## }

# Taxon set
TAX=("Metazoa" "Holozoa")

# Parameter
Parameter=()
for file in $PRM/parameter_*.txt; do
    base_name=$(basename "$file" .txt); prm_name="${base_name#parameter_}"
    Parameter+=("$prm_name")
done

# Run analysis
for tax in "${TAX[@]}"; do

##  prepare_table "$tax"

    TABLE="${BAY}/binary.table.GS_TET.species_${tax}.txt"
    NEX="${BAY}/${tax}.nex"

    # output basename
    basename=$(basename "$TABLE" .txt)

    # RUN
    for prm in "${Parameter[@]}"; do

        echo "BayesTraits | Tax_range: $tax | Mode: $prm | Starting to analyze ..."
        # run program
        $PRG/BayesTraitsV4 $NEX $TABLE < $PRM/parameter_${prm}.txt > /dev/null

        # move output files to result folder
        mv -f $BAY/*.Log.txt $RES/${prm}.${basename}.Log.txt
        mv -f $BAY/*.Schedule.txt $RES/${prm}.${basename}.Schedule.txt
        mv -f $BAY/*.Stones.txt $RES/${prm}.${basename}.Stones.txt 2>/dev/null

        echo "BayesTraits | Tax_range: $tax | Mode: $prm | End"
    
    done

done
