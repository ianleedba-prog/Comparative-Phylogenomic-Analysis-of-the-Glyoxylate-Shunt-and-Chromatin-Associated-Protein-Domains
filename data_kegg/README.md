# KEGG Metabolic Module Presence Analysis

**Associated manuscript: NCOMMS-26-022536 (revised)**

*Resurrecting the glyoxylate cycle constrains TET epigenetics via governing 2-oxoglutarate and its competitive inhibitors*  
Il-Hwan Lee, Jinmi Choi, So-Yeon Kim, Young Ah Kim, Yufei Li, Joo-Youn Cho, Eun-Jung Cho, Hong-Duk Youn

**Corresponding author:** Hong-Duk Youn (hdyoun@snu.ac.kr)  

**Figure produced by this workflow:** Supplementary Figure 1A

---

## Overview

This workflow surveys the presence of all 477 KEGG metabolic modules across the full set of species registered in the KEGG taxonomy database. Module presence data are retrieved from the KEGG web service, aggregated by taxonomic group, and visualized as a heatmap with hierarchical column clustering. The resulting figure (Supplementary Figure 1A) identifies three major clusters of metabolic modules distinguished by their phylogenetic distribution patterns across eukaryotic and prokaryotic groups — with the glyoxylate cycle (M00012) highlighted as a module in the cluster that shows pronounced loss at the metazoan boundary.

The workflow consists of two sequential steps:

**Step 1 — Data retrieval (`s01_download_keggdb.sh`):**  
Queries the KEGG taxonomy browser for each of the 477 modules to retrieve the list of KEGG species abbreviations in which the module is annotated as present. Output is one text file per module stored in `$DIR_kegg/species_abb/`.

**Step 2 — Data processing and visualization (`s02_module_presence_sh.R`):**  
Integrates all per-module species lists with the KEGG taxonomy table, computes module presence fractions at multiple taxonomic levels, filters out low-coverage modules and groups, performs hierarchical clustering, and generates the heatmap and cluster assignment tables.

---

## Directory Structure

All files must be placed under a single root directory referred to as `$DIR_kegg`. Update the hardcoded path (`/directory/to/data_kegg`) in each script to match your local environment before running.

```
 $DIR_kegg/
 │
 ├── README.md
 │
 ├── ── Main pipeline scripts ────────────────────────────────────────────────────
 ├── s01_download_keggdb.sh              ← Step 1: Download per-module species lists from KEGG
 ├── s02_module_presence_sh.R            ← Step 2: Aggregate, compute fractions, heatmap
 │
 ├── ── Provided input data files ────────────────────────────────────────────────
 ├── KEGG_Taxonomy_list.csv              ← KEGG species taxonomy table (~10,554 entries)
 ├── rest.kegg.jp_list_module.txt        ← Complete list of 477 KEGG metabolic modules
 │
 ├── ── Directory created during Step 1 ─────────────────────────────────────────
 └── species_abb/                        ← Per-module species abbreviation files
     ├── M00001_species_abb.txt
     ├── M00002_species_abb.txt
     └── ... (477 files total)
 │
 └── ── Results directory (created by Step 2) ────────────────────────────────────
     └── results/
         ├── module.presence.species.csv     ← Module presence per species (binary)
         ├── module.presence.group1.csv      ← Module presence fraction per Group1 taxon
         ├── module.presence.group2.csv      ← Module presence fraction per Group2 taxon
         ├── module_presence_final.csv       ← Filtered final matrix (input to heatmap)
         ├── module.presence.final.pdf       ← Heatmap figure (Supplementary Figure 1A)
         ├── cluster1_modules.csv            ← Modules in Cluster 1 (mixed distribution)
         ├── cluster2_modules.csv            ← Modules in Cluster 2 (gained at metazoan boundary)
         ├── cluster3_modules.csv            ← Modules in Cluster 3 (lost at metazoan boundary)
         └── clustered_modules.xlsx          ← All three cluster CSVs merged into one workbook
```

---

## Software Requirements

### R packages

| Package | Version | Purpose |
|---|---|---|
| R base | v4.4.3 (arm64) | Runtime environment |
| data.table | v1.18.2.1 | Fast tabular data loading |
| pheatmap | v1.0.13 | Heatmap visualization with hierarchical clustering |
| seriation | v1.5.8 | OLO (Optimal Leaf Ordering) reordering of dendrogram |

---

## Provided Data Files

### KEGG_Taxonomy_list.csv
Taxonomy classification table for all species registered in the KEGG database (~10,554 entries). Obtained by downloading and reformatting the KEGG organism list from https://www.kegg.jp/kegg-bin/show_brite?htext=br08601. Each row represents one KEGG species entry. Columns:
- `Supergroup`: Broad classification; `Eukaryotes` or `Prokaryotes`
- `Group1`: Major taxonomic group (e.g., `Animals`, `Plants`, `Fungi`, `Bacteria`, `Archaea`)
- `Group2`: Finer taxonomic group used as the primary display unit in the heatmap (e.g., `Mammals`, `Insects`, `Ascomycetes`, `Eudicots`)
- `Group3`: Most granular classification available (e.g., `Primates`, `Rodents`)
- `Abbreviation`: KEGG 3–5 character species abbreviation (e.g., `hsa` for *Homo sapiens*)
- `Species`: Full species name

Rows with missing `Abbreviation` values (non-species entries) are excluded during data processing.

### rest.kegg.jp_list_module.txt
Complete list of 477 KEGG metabolic modules retrieved from the KEGG REST API (`https://rest.kegg.jp/list/module`). Tab-delimited, two columns, no header:
- Column 1: KEGG module identifier (e.g., `M00001`)
- Column 2: Module name and description (e.g., `Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate`)

This file defines the complete set of modules queried in Step 1 and serves as the reference for module name annotation in the cluster output tables.

---

## Provided Result Files

The following output files from Step 2 are provided in this submission for reference and reproducibility verification. They are located in `$DIR_kegg/results/`.

| File | Description |
|---|---|
| `module.presence.species.csv` | Binary module presence per KEGG species; rows = species, columns = modules |
| `module.presence.group1.csv` | Module presence fraction (0–100%) aggregated at Group1 level |
| `module.presence.group2.csv` | Module presence fraction (0–100%) aggregated at Group2 level |
| `module_presence_final.csv` | Filtered final matrix used as direct input to the heatmap; rows = taxonomic display groups (Group2 for eukaryotes; Group1 for prokaryotes), columns = modules passing the 10% presence threshold in at least one group |
| `module.presence.final.pdf` | Heatmap output (Supplementary Figure 1A) |
| `clustered_modules.xlsx` | Excel workbook containing the cluster assignment table for all modules in the heatmap, with each sheet corresponding to one of the three column clusters |

---

## Pipeline Description

### Step 1: Per-module species list retrieval — `s01_download_keggdb.sh`

**Output:** `$DIR_kegg/species_abb/{module_id}_species_abb.txt` — one file per module (477 total)

---

### Step 2: Module presence aggregation and heatmap — `s02_module_presence_sh.R`

**Output:** Presence matrices (CSV), heatmap (PDF), cluster assignment tables (CSV, XLSX)  
**Figure:** Supplementary Figure 1A

**1. Per-species presence matrix construction**  
For each module, reads the corresponding `{module_id}_species_abb.txt` file and marks presence (`1`) for each KEGG species abbreviation found. Species not in a module's list are assigned absence (`0`). The result is a species × module binary matrix merged with the taxonomy table.  
**2. Fraction calculation by taxonomic group + Final matrix assembly**  
Module presence fractions (0–100%) are computed at two levels by dividing the count of species with a module by the total number of species in the group:  
- **Group1** (e.g., Animals, Plants, Fungi, Bacteria, Archaea)
- **Group2** (e.g., Mammals, Insects, Ascomycetes, Eudicots)
**3. Clusters**  
- **Cluster 1:** Modules with mixed taxonomic distribution (present broadly across both eukaryotes and prokaryotes)
- **Cluster 2:** Modules gained or expanded around the metazoan boundary (enriched in animals relative to other groups)
- **Cluster 3:** Modules showing pronounced loss at the metazoan boundary (present in most eukaryotes and prokaryotes but absent in most Metazoa); this cluster contains the glyoxylate cycle (M00012)

Per-cluster module lists with KEGG IDs and descriptions are saved as individual CSVs and merged into `clustered_modules.xlsx`.

---

## Contact

**Il-Hwan Lee** (first author, code maintainer) - ianlee.dba@snu.ac.kr  
Department of Biomedical Sciences, Seoul National University College of Medicine, Seoul 03080, Republic of Korea

**Hong-Duk Youn** (corresponding author) - hdyoun@snu.ac.kr  
Department of Biomedical Sciences, Seoul National University College of Medicine, Seoul 03080, Republic of Korea
