# README — Comparative Phylogenomic Analysis of the Glyoxylate Shunt and Chromatin-Associated Protein Domains

**Associated Manuscript: NCOMMS-26-022536 / Submission Date: 17th March 26 **
*Resurrecting the glyoxylate cycle constrains TET epigenetics through 2-oxoglutarate limitation*
Il-Hwan Lee, Jinmi Choi, So-Yeon Kim, Young Ah Kim, Yufei Li, Joo-Youn Cho, Eun-Jung Cho, Hong-Duk Youn

**Corresponding author:** Hong-Duk Youn (hdyoun@snu.ac.kr)
**Figures produced by this workflow:** Figure 1, Figure 2, Supplementary Figure 1A, Supplementary Figure 2A, Supplementary Figure 3A–C

---

## Table of Contents

1. [Overview](#overview)
2. [Directory Structure](#directory-structure)
3. [Software Requirements](#software-requirements)
4. [Provided Data Files](#provided-data-files)
5. [Evolutionary Boundary Groups](#evolutionary-boundary-groups)
6. [Pipeline Description](#pipeline-description)
7. [How to Run](#how-to-run)

---

## Overview

This workflow performs a comparative phylogenomic survey of glyoxylate shunt enzymes and chromatin-associated protein domains across eukaryotes, and identifies orthologous groups of histone lysine demethylases (KDMs) across three evolutionary boundaries. The pipeline proceeds through two parallel tracks:

**Track A — Broad phylogenomic survey (1,610 species) → Figure 1**
Determines presence/absence of 42 protein domain families across the full eukaryotic dataset. Produces a phylum-level fraction matrix visualized as a heatmap, and computes Jaccard overlap scores quantifying co-occurrence between isocitrate lyase (ICL) and each chromatin-associated domain across taxa.

**Track B — Representative species analysis (112 species) → Figure 2, Supplementary Figure 3**
Performs ortholog identification and subfamily classification of JmjC-containing KDMs, AOD-containing KDMs (KDM1 family), and Tet_JBP dioxygenases across three evolutionary boundaries: (i) premetazoans vs. early metazoans (PEM), (ii) metazoans with vs. without the glyoxylate shunt (MET), and (iii) metazoans vs. embryophytes (PLT). Additionally identifies conserved catalytic and cofactor-binding residues within JmjC-containing KDM orthologs.

---

## Directory Structure

All scripts and data files must be placed under a single root directory referred to as `$DIR` throughout this document. Update all hardcoded directory paths (`/Directory/to/data`) in each script to match your local environment before running.

```
$DIR/
│
├── README.md
│
├── ── Main pipeline scripts ──────────────────────────────────────────────────
├── s01a_download_fasta.sh                 ← Step 1a: Download protein FASTAs
├── s01b_download_annotation.sh            ← Step 1b: Download annotation files (GTF/GFF)
├── s02a_phylogenetic_analysis_all.sh      ← Step 2a: Full-scale analysis (1,610 species)
├── s02b_phylogenetic_analysis_rep.sh      ← Step 2b: Representative species analysis (112 species)
├── s03_gene_counts.R                      ← Step 3: Heatmap + Jaccard index (Figure 1)
├── s04_ortholog_counts.R                  ← Step 4: KDM ortholog counts (Figure 2)
├── s05_conserved_sites.R                  ← Step 5: Catalytic site conservation (Suppl. Fig. 3B, C)
│
├── ── Job module scripts (called by Step 2) ───────────────────────────────────
├── job_busco.sh                           ← Assembly completeness validation (BUSCO)
├── job_hmmsearch.sh                       ← HMM search + optional MCL clustering + alignment
├── job_alignment.sh                       ← MSA + trimming + IQ-TREE + TreeShrink
├── job_pfamscan.sh                        ← Pfam domain architecture annotation
├── job_possvm.sh                          ← Ortholog group identification (POSSVM)
├── job_isoformshrink.sh                   ← Isoform deduplication per genomic locus
├── job_tree.sh                            ← Final tree construction per group/domain
├── job_conserved_sites.sh                 ← Catalytic/cofactor binding site identification
│
├── ── Provided data files ─────────────────────────────────────────────────────
├── Metadata_species_info.txt              ← Full species metadata (1,610 species)
├── Metadata_rep.species_info.txt          ← Representative species metadata (112 species)
├── protein_domains_hmm.csv                ← Domain class definitions + HMM search parameters
├── Taxonomy_info.txt                      ← Pre-computed NCBI taxonomy per species (provided)
├── Hsap_gene_names.txt                    ← Human reference gene names for POSSVM
├── HsapAtha_gene_names.txt                ← Human + Arabidopsis reference gene names (PLT)
├── Sps.prefix_PEM.txt                     ← Species list: Premetazoa vs. Early metazoa (20 spp.)
├── Sps.prefix_MET.txt                     ← Species list: GS+ vs. GS- Metazoa (64 spp.)
├── Sps.prefix_PLT.txt                     ← Species list: Metazoa vs. Embryophyta (47 spp.)
├── Sps.tree_PEM.newick.txt                ← Species tree for PEM group (Newick)
├── Sps.tree_MET.newick.txt                ← Species tree for MET group (Newick)
├── Sps.tree_PLT.newick.txt                ← Species tree for PLT group (Newick)
├── Sps.tree.newick.txt                    ← Full representative species tree (Newick)
│
├── ── Directories created during runtime ──────────────────────────────────────
├── Pfam/                                  ← Pfam HMM database
│   ├── Pfam-A.hmm
│   ├── Pfam-A.hmm.dat
│   └── hmm/                               ← Per-domain .hmm files (provided in submission)
├── taxonomy/                              ← NCBI taxonomy dump
├── ensembl_source/                        ← Ensembl species registry files (provided in submission)
├── lineages/                              ← BUSCO lineage database (eukaryota_odb10) (will be downloaded automatically when job_busco.sh is run)
├── protein_fa/                            ← Per-species protein FASTA files (directory where protein fasta should be stored)
├── compiled_fa/                           ← Compiled multi-species FASTA
├── annotation/                            ← Per-species GTF/GFF annotation files (directory where annotation files should be stored)
├── busco/                                 ← BUSCO output per species
└── analysis/                              ← All analysis outputs
    ├── searches/                          ← HMM search results per domain
    ├── alignments/                        ← Per-cluster MSA and tree files
    ├── gene_counts/                       ← Gene count tables (All.genecounts.csv, Rep.genecounts.csv)
    ├── gene_trees/                        ← Compiled gene tree files (.treefile) for POSSVM
    ├── gene_sequences/                    ← Compiled gene sequence FASTA files
    ├── orthologs/                         ← POSSVM ortholog results
    │   ├── Hsap/                          ← PEM and MET groups (Human reference)
    │   │   ├── Tet_JBP/
    │   │   │   ├── PEM/                   ← Isoform-deduplicated data and trees
    │   │   │   └── MET/
    │   │   ├── JmjC/
    │   │   └── AOD/
    │   └── HsapAtha/                      ← PLT group (Human + Arabidopsis reference)
    │       ├── Tet_JBP/
    │       ├── JmjC/
    │       └── AOD/
    └── results/                           ← *** Final output tables and figures ***
```

---

## Software Requirements

All versions listed below match those used in the paper.

### Shell-based tools

| Software | Version | Purpose |
|---|---|---|
| ncbi-datasets-cli | v16.40.1 | Proteome and annotation download from NCBI |
| BUSCO | v6.0.0 | Assembly completeness validation |
| HMMER (hmmsearch, hmmpress, esl-sfetch) | v3.4 | Protein domain HMM search and sequence extraction |
| DIAMOND | v2.1.10 | All-vs-all pairwise alignment for MCL clustering |
| MCL | v22.282 | Markov Cluster Algorithm for homology grouping |
| MAFFT | v7.525 | Multiple sequence alignment |
| ClipKIT | v1.4.1 | Alignment trimming |
| IQ-TREE2 | v2.4.0 | Maximum-likelihood phylogenetic tree construction |
| TreeShrink | v1.3.9 | Long-branch outlier sequence detection and removal |
| POSSVM | — | Ortholog group identification from gene trees |
| pfam_scan (pfam_scan.pl) | v1.6 | Pfam domain architecture annotation |
| seqkit | v2.10.0 | FASTA deduplication and sequence grep |
| Samtools | v1.21 | FASTA indexing |
| bedtools | v2.31.1 | Domain coordinate merging and region expansion |
| bioawk | v1.0 | FASTA parsing and sequence length operations |
| GNU parallel | — | Parallel BUSCO execution |
| ETE Toolkit (ete3) | v3.1.3 | Tree manipulation (dependency of POSSVM) |

POSSVM source code is available at: https://github.com/xgrau/possvm-orthology; Clone the repository and place it at `$DIR/possvm-orthology-master/`
cite:   (1) Huerta-Cepas, J., Dopazo, H., Dopazo, J. et al. The human phylome. Genome Biol 8, R109 (2007). https://doi.org/10.1186/gb-2007-8-6-r109
        (2) Jaime Huerta-Cepas, François Serra, Peer Bork, ETE 3: Reconstruction, Analysis, and Visualization of Phylogenomic Data, Molecular Biology and Evolution, Volume 33, Issue 6, June 2016, Pages 1635–1638, https://doi.org/10.1093/molbev/msw046
        (3) Xavier Grau-Bové, Arnau Sebé-Pedrós, Orthology Clusters from Gene Trees with Possvm, Molecular Biology and Evolution, Volume 38, Issue 11, November 2021, Pages 5204–5208, https://doi.org/10.1093/molbev/msab234
        (4) A. J. Enright, S. Van Dongen, C. A. Ouzounis, An efficient algorithm for large-scale detection of protein families, Nucleic Acids Research, Volume 30, Issue 7, 1 April 2002, Pages 1575–1584, https://doi.org/10.1093/nar/30.7.1575

### R packages

| Package | Version | Purpose |
|---|---|---|
| R base | v4.4.3 (arm64) | Runtime environment |
| dplyr | v1.2.0 | Data manipulation |
| tidyr | v1.3.2 | Data reshaping |
| data.table | v1.18.2.1 | Fast tabular data loading |
| pheatmap | v1.0.13 | Heatmap visualization |
| ggplot2 | v4.0.2 | General plotting |
| ggtree | v3.14.0 | Phylogenetic tree visualization |
| ape | 5.8-1 | Newick tree reading |
| phytools | 2.5-2 | Phylogenetic utility functions |

---

## Provided Data Files

The following files are included in this code submission and must be present in `$DIR` before running.

### Metadata_species_info.txt
Full species metadata for 1,610 eukaryotic species (Track A). Tab-delimited, with a header row. Columns:
- `Taxonomy_name`: Full species name
- `Abbreviation`: 4–6 character prefix used in all FASTA headers and output file names
- `NCBI_Taxonomy_ID`: NCBI numerical taxonomy ID
- `Source`: Download source; one of `NCBI`, `WormBase_18`, `Ensembl_115`, `Ensembl_metazoa_61`, `Ensembl_fungi_61`, `Ensembl_plants_61`, `Ensembl_protists_61`, or `Other` (requires manual download)
- `File`: File identifier for download; NCBI genome accession (e.g., `GCF_000001405.40`) for NCBI entries, or filename stem for WormBase/Ensembl entries

### Metadata_rep.species_info.txt
Metadata for 112 representative species (Track B). Identical column structure to `Metadata_species_info.txt` with one additional column:
- `Annotation`: GTF or GFF annotation filename; set to `na` for three species (NCBI Taxonomy IDs: 749232, 433461, 749231) lacking available annotation files

### protein_domains_hmm.csv
Defines the 42 protein domain classes surveyed in this study. No header. Tab-delimited columns:
- Col 1 (`Class`): Broad functional class (e.g., `Metabolism`, `DNA_modification`, `Histone_modification`, `Readers`)
- Col 2 (`Type`): Sub-type within class (e.g., `TCAcycle`, `GlyoxylateShunt`, `Demethylation`)
- Col 3 (`Family`): Domain family name used as the column identifier in all output matrices
- Col 4 (`Domains`): Comma-separated list of Pfam HMM profile names (e.g., `Tet_JBP,2OG-FeII_Oxy_6`)
- Col 5 (`homology_cluster`): `HC` = apply MCL-based homology clustering; `na` = skip. Only `Tet_JBP`, `JmjC`, and `AOD` use `HC`
- Col 6 (`inflation`): MCL inflation value; `1.5` for Tet_JBP, `1.3` for JmjC and AOD, `1.1` for all others
- Col 7 (`min_phylo_size`): Minimum sequences per cluster to proceed to phylogenetic analysis; `2` for all entries

### hmm files (from Pfam release 37.2 / InterPro) (provided in submission)
Source data was downloaded as following steps:

    wget -qO- https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.2/Pfam-A.hmm.gz > $PFAM/Pfam-A.hmm.gz
    wget -qO- https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.2/Pfam-A.hmm.dat.gz > $PFAM/Pfam-A.hmm.dat.gz
    gunzip $PFAM/Pfam-A.hmm.gz
    gunzip $PFAM/Pfam-A.hmm.dat.gz
    hmmpress $PFAM/Pfam-A.hmm

    # Note: Individual per-domain .hmm files used by job_hmmsearch.sh are provided in $PFAM/hmm/ as part of this code submission. 
    # These were extracted from Pfam-A.hmm using hmmfetch and correspond to the Pfam accessions listed in column 4 of protein_domains_hmm.csv.

### Taxonomy_info.txt
Pre-computed NCBI taxonomy table mapping each of the 1,610 species to its full taxonomic lineage (Species → Genus → Family → Order → Class → Phylum → Kingdom → Superkingdom). 
Provided directly by the author, which was generated by using new_taxdump.tar.gz downloaded from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/.

### Species prefix lists
Plain-text files defining the species included in each evolutionary boundary group. The first entry is always `Hsap` (human), which serves as the reference anchor in downstream analyses.

### Reference gene name files
`Hsap_gene_names.txt`: Human gene symbols used by POSSVM to name ortholog groups.
`HsapAtha_gene_names.txt`: Combined human and *Arabidopsis thaliana* gene symbols for the PLT group, enabling cross-kingdom ortholog annotation.

### Species tree files
Pre-constructed reference species trees in Newick format, used for visualization in `s04_ortholog_counts.R` and `s05_conserved_sites.R`.

---

## Evolutionary Boundary Groups

| Code | Comparison | Species file | n species | Tree file | POSSVM ref |
|---|---|---|---|---|---|---|
| **PEM** | Premetazoans vs. Early Metazoans | `Sps.prefix_PEM.txt` | 20 | `Sps.tree_PEM.newick.txt` | Hsap |
| **MET** | Metazoans with GS vs. without GS | `Sps.prefix_MET.txt` | 64 | `Sps.tree_MET.newick.txt` | Hsap |
| **PLT** | Metazoans vs. Embryophytes | `Sps.prefix_PLT.txt` | 47 | `Sps.tree_PLT.newick.txt` | HsapAtha |

---

## Pipeline Description

Scripts must be executed in order. Steps 1–2 are shell scripts; Steps 3–5 are R scripts.

---

### Step 1a-b: Download protein FASTAs, annotation files

s01a_downlaod_fasta.sh: Downloads predicted proteome FASTA files for all species (1,610) from major three sources → species abbreviation is prepended to all sequence headers:
**NCBI** — uses `ncbi-datasets-cli` (`datasets download genome --include protein`). Accession numbers are taken from the `File` column.
**WormBase** — WBPS18
**Ensembl** — vertebrates release 115; non-vertebrates release 6
**Other publication** — individual genome publications. → → → → Please download and process manually as described in s01a_downlaod_fasta.sh (end of the script)


s01b_download_annotation.sh: Downloads annotation files (GTF/GFF) for representative species (112).

---

### Step 2a: Full-scale phylogenomic survey — `s02a_phylogenetic_analysis_all.sh`

**Output:** `$DIR/busco_scores.txt`, `$DIR/analysis/searches/*`, `$DIR/analysis/gene_counts/All.genecounts.csv`
**Figures:** Supplementary Figure 2A-B, Figure 1A–B

**Substep 1 — BUSCO validation (`job_busco.sh`)**
Runs BUSCO v6.0.0 on all 1,610 assemblies using the `eukaryota_odb10` lineage dataset in protein mode. Per-species completeness scores are summarized in `$DIR/busco_scores.txt`. 
The median completeness across all assemblies was 96.1%.
**Substep 2 — FASTA compilation**
Concatenates all species FASTAs into `Compiled.All.fasta`.
**Substep 3 — HMM search (all 1,610 species)**
For each of the 42 domain classes in `protein_domains_hmm.csv`, calls `job_hmmsearch.sh`.
Results are written to `$DIR/analysis/searches/All.*.genes.list`.
**Substep 4 — Gene count compilation**
Concatenates all per-domain gene lists into `All.genecounts.csv`, a two-column file (protein ID, domain family name).
This is the primary input for `s03_gene_counts.R`.

---

### Step 2b: Representative species analysis — `s02b_phylogenetic_analysis_rep.sh`

**Output:** `$DIR/analysis/alignment/*`, `$DIR/analysis/gene_sequences/*`, `$DIR/analysis/gene_trees/*`, `$DIR/analysis/orthologs/*`.
**Figures:** Figure 2A–F, Supplementary Figure 3A–D

**Substep 1 — FASTA compilation (112 species)**
**Substep 2 — HMM search (112 representative species)**
For the three domains (Tet_JBP, JmjC, AOD), this additionally performs: all-vs-all `DIAMOND blastp`, `MCL` clustering at the domain-specific inflation value, `pfam_scan.pl` domain architecture annotation, and per-cluster MSA + tree construction via `job_alignment.sh`.
**Substep 3 — Gene count compilation (112 species)**
`Rep.genecounts.csv`.
**Substep 4 — Ortholog group identification (`job_possvm.sh`)**
For each domain (Tet_JBP, JmjC, AOD) and reference set (Hsap; HsapAtha for PLT group), runs POSSVM on all gene trees to assign ortholog group names anchored to human (and Arabidopsis) gene symbols. 
**Substep 5 — Isoform deduplication (`job_isoformshrink.sh`)**
Maps each ortholog protein ID to its genomic locus using annotation files. Among all transcript isoforms at the same gneomic locus, retains only the single longest protein. 
**Substep 6 — Final tree construction (`job_tree.sh`)**
Constructs maximum-likelihood trees from the isoform-deduplicated sequence sets.
**Substep 7 — Conserved site identification (`job_conserved_sites.sh`)**
Identifies amino acid identities at catalytic and cofactor-binding positions for JmjC-containing KDMs in PEM and MET groups.

---

### Step 3: Gene counts and Jaccard index — `s03_gene_counts.R`

**Output:** Gene presence/absence matrix-related Tables & Heatmap, Jaccard index-related Tables & Bargraphs
**Figures:** Figure 1A-H, Supplementary Figure 2A

---

### Step 4: KDM ortholog counts — `s04_ortholog_counts.R`

**Output:** Count matrices, fraction matrices, heatmaps, bar graphs, species tree PDFs
**Figures:** Figure 2A–F, Supplementary Figure 3A-D

Iterates over all combinations of evolutionary boundary group (PEM, MET, PLT) × domain (Tet_JBP, JmjC, AOD). For each:

---

### Step 5: Conserved catalytic site visualization — `s05_conserved_sites.R`

**Output:** Annotated tree PDFs, tip data CSVs
**Figures:** Supplementary Figure 3B-C


---

## How to Run

```bash
# 1. Data acquisition
# --- KEGG db, Protein fasta, annotation files ---

# 2. Phylogenomic analysis
bash s02a_phylogenetic_analysis_all.sh   # 1~2 hours; BUSCO + HMM survey
bash s02b_phylogenetic_analysis_rep.sh   # 2~3 days; ortholog analysis + trees

# 3–5. R-based visualization
Rscript s03_gene_counts.R
Rscript s04_ortholog_counts.R
Rscript s05_conserved_sites.R
```

---

## Contact

**Il-Hwan Lee** (first author, primary code maintainer)
Department of Biomedical Sciences, Seoul National University College of Medicine, Seoul 03080, Republic of Korea

**Hong-Duk Youn** (corresponding author)
hdyoun@snu.ac.kr
Department of Biomedical Sciences, Seoul National University College of Medicine, Seoul 03080, Republic of Korea
