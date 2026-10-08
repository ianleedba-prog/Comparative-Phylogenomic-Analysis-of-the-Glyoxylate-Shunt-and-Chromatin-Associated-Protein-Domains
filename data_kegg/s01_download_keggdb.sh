#!/bin/bash
# ==============================================================================
# Software:
# ==============================================================================

# Directories
DIR_kegg = "/Directory/to/data_kegg"

# KEGG-Module list
mod_list = "$DIR_kegg/rest.kegg.jp_list_module.txt"

# Retrieve data from KEGG taxonomy db (module presence per species): Provided in submission
mkdir -p $DIR_kegg/species_abb
while read module; do
    curl -s "https://www.kegg.jp/kegg-bin/show_brite?htext=br08611&pruning=taxmap&brel=off&hier=20&highlight=join_brite&mapper=taxmap%2ds%20${module}" -o page.html
    grep "10003" page.html | cut -f1 > $DIR_kegg/species_abb/${module}_species_abb.txt
    rm page.html
done < <(cut -f1 $mod_list)
