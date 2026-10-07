[README.md](https://github.com/user-attachments/files/33150156/README.md)
# Comparative Phylogenomic Analysis of the Glyoxylate Shunt and Tet_JBP

**Associated manuscript: NCOMMS-26-022536 (revised)**

*Resurrecting the glyoxylate cycle constrains TET epigenetics via governing 2-oxoglutarate and its competitive inhibitors*
Il-Hwan Lee, Jinmi Choi, So-Yeon Kim, Young Ah Kim, Yufei Li, Joo-Youn Cho, Eun-Jung Cho, Hong-Duk Youn

**Corresponding author:** Hong-Duk Youn (hdyoun@snu.ac.kr)

**Figures produced by this workflow:** Figure 1, Figure 2, Supplementary Figures 1–3, Supplementary Figure 9

---

## Overview

**1,609 eukaryotes.** Presence and absence of 42 protein domain families, compiled into phylum- and species-level matrices, and the co-occurrence between the shunt and each chromatin-associated domain. → **Figure 1**, Supplementary Figures 1–2

**259 species**, a subset of the 1,609 carrying a species tree. Within it, the 125 holozoans and the 112 metazoans are the trees on which phylogenetic signal, ancestral states and Bayesian model comparison are estimated. → **Figure 2**, Supplementary Figures 2–3

A separate sub-pipeline under `JmjC/` covers the histone demethylase survey on 112 species. It is self-contained, does not feed the main pipeline, and supports **Supplementary Figure 9** only.

---

## Directory structure

All scripts take a single root directory `$DIR`. Update the hardcoded path (`/Directory/to/data`) at the top of each script before running.

```
$DIR/
│
├── ── Main pipeline ───────────────────────────────────────────────────────────
├── s01a_download_fasta.sh             ← protein FASTAs (1,609 species)
├── s01b_download_annotation.sh        ← GTF/GFF annotations
├── s02_phylogenetic_analysis_all.sh   ← BUSCO validation + HMM domain search
├── s03_gene_counts.R                  ← domain matrices, co-occurrence, 259-species subset
├── s04_Dstatistic.R                   ← phylogenetic signal of each character
├── s05_ancestral_state_reconstruction.R ← four-state model, ancestral states, branch events
├── s06_Bayesian.sh                    ← BayesTraits runs
├── s07_bayesian_traits.R              ← posteriors, Bayes factors, clade decomposition
│
├── job_busco.sh                       ← called by s02
├── job_hmmsearch.sh                   ← called by s02
│
├── ── Provided data ───────────────────────────────────────────────────────────
├── Metadata_species_info.txt          ← 1,609 species: source, accession, abbreviation
├── Taxonomy_info.txt                  ← NCBI lineage, 1,609 species
├── Taxonomy_info_259.txt              ← NCBI lineage, 259-species subset
├── protein_domains_hmm.csv            ← 42 domain classes and their HMM search parameters
├── Eukaryota_259.nwk                  ← species tree of the 259-species subset
├── Eukaryota_259_dichotomous.nwk      ← the same tree with polytomies resolved (see s04)
├── ensembl_source/                    ← Ensembl species registries
│
├── analysis/                          ← outputs of the main pipeline (provided)
│   ├── searches/                      ← per-domain HMM hits (gene lists, domain tables)
│   ├── gene_counts/                   ← 1609.genecounts.csv
│   └── results/                       ← final tables and figures
│
├── Bayesian_traits/                   ← BayesTraits inputs and outputs (provided)
│   ├── Holozoa.nex, Metazoa.nex       ← trees handed to BayesTraits
│   ├── binary.table.GS_TET.species_*.txt   ← GS and TET per species
│   ├── parameters/                    ← one file per model (DP, DP with q13 = q24, IDP)
│   └── results/                       ← .Log.txt, .Schedule.txt, .Stones.txt per run
│
└── JmjC/                              ← KDM sub-pipeline, self-contained (see below)
    ├── s01_phylogenetic_analysis_rep.sh, s02_ortholog_counts.R, s03_conserved_sites.R
    ├── job_*.sh                       ← alignment, Pfam scan, POSSVM, isoform collapse, trees
    ├── Metadata_rep.species_info.txt, Sps.prefix_*.txt, Sps.tree*.newick.txt
    ├── Hsap_gene_names.txt, HsapAtha_gene_names.txt
    ├── Pfam/
    └── alignments/, gene_trees/, gene_sequences/, orthologs/, searches/, results/
```
**Not included.** 
The raw hmmsearch output under `analysis/searches/` (`*.seqs.fasta`, `*.domains.fasta`, `*.domtable`) is left out to keep the repository within GitHub's size limits;
`s02` regenerates it, and the gene lists and domain tables the R scripts read are included.

---

## Pipeline

### Step 1 — Data acquisition (`s01a`, `s01b`)
Downloads proteomes for the 1,609 species from NCBI, WormBase (WBPS18) and Ensembl (vertebrates 115, non-vertebrates 61), prepending the species abbreviation to every FASTA header, and the annotation files used by the `JmjC/` sub-pipeline. A handful of species from individual genome publications are downloaded by hand, as described at the end of `s01a`.

### Step 2 — Domain survey (`s02_phylogenetic_analysis_all.sh`)
BUSCO v6.0.0 (`eukaryota_odb10`, protein mode) on every assembly, then an `hmmsearch` for each of the 42 domain classes in `protein_domains_hmm.csv` across the compiled proteome. Per-domain hits land in `analysis/searches/`, and are compiled into `analysis/gene_counts/1609.genecounts.csv`, the input to Step 3.

### Step 3 — Domain matrices and co-occurrence (`s03_gene_counts.R`)
Builds the species- and phylum-level count and fraction matrices, the clustered heatmap, and Jaccard overlap between the shunt and each chromatin-associated domain at four taxonomic ranks. Defines the two binary characters and the four states, draws the 259-species subset, and computes the odds ratio of carrying Tet_JBP given the shunt (with the h-module as a second trait) for Eukaryota, Opisthokonta, Holomycota, Holozoa and Metazoa.

### Step 4 — Phylogenetic signal (`s04_Dstatistic.R`)
Fritz and Purvis' *D* for each of eight binary characters in each of five nested taxon sets, with a leave-one-out jackknife interval and the two null distributions rescaled onto the *D* axis. `caper` needs a dichotomous tree, so the polytomies of the 259-species tree are resolved once and stored as `Eukaryota_259_dichotomous.nwk`.

### Step 5 — Ancestral state reconstruction (`s05_ancestral_state_reconstruction.R`)
Fits a dependent four-state model with all eight single-character transitions to the 125 holozoans: marginal states at every node, the crown and parent state of each named clade, and the expected number of each transition on every branch.

### Step 6 — BayesTraits runs (`s06_Bayesian.sh`)
Runs BayesTraits v4.1.3 on the holozoan (125) and metazoan (112) trees under three models — the full dependent model, the dependent model with the two rates of shunt gain held equal (`q13 = q24`), and the independent model.

### Step 7 — Posteriors, Bayes factors, clade decomposition (`s07_bayesian_traits.R`)
Reads the BayesTraits logs for the posterior of the two rates of shunt gain, the stone files for the marginal likelihoods and the log Bayes factors, and refits the maximum-likelihood model with each metazoan clade removed in turn to see which clades carry the difference between the two rates.

### `JmjC/` — Histone demethylase survey (supporting)
Independent of Steps 3–7 and run separately. On 112 species it identifies ortholog groups of JmjC- and AOD-containing KDMs and of Tet_JBP with POSSVM, across three evolutionary boundaries — premetazoans against early metazoans (PEM, 20 spp.), metazoans with against without the shunt (MET, 64 spp.), and metazoans against embryophytes (PLT, 47 spp.) — and scores the amino acid at the catalytic and cofactor-binding positions of JmjC KDMs. `s01` builds the trees, `s02` counts orthologs, `s03` reads off the conserved sites.

---

## Software

| Software | Version | Used in |
|---|---|---|
| ncbi-datasets-cli | v16.40.1 | s01 |
| BUSCO | v6.0.0 | s02 |
| HMMER | v3.4 | s02, JmjC |
| BayesTraits | v4.1.3 | s06 |
| DIAMOND / MCL | v2.1.10 / v22.282 | JmjC |
| MAFFT / ClipKIT | v7.525 / v1.4.1 | JmjC |
| IQ-TREE2 / TreeShrink | v2.4.0 / v1.3.9 | JmjC |
| POSSVM / ETE3 | — / v3.1.3 | JmjC |
| pfam_scan.pl | v1.6 | JmjC |
| seqkit / samtools / bedtools / bioawk | v2.10.0 / v1.21 / v2.31.1 / v1.0 | JmjC |

| R package | Version | Used in |
|---|---|---|
| R base | v4.4.3 (arm64) | all |
| ape | v5.8-1 | s03, s04, s05, s07, JmjC |
| phytools | v2.4-4 (v2.5-2 in JmjC) | s05, s07, JmjC |
| caper | v1.0.4 | s04 |
| expm | v1.0-0 | s05 |
| dplyr / data.table | v1.2.0 / v1.18.2.1 | s03, JmjC |
| pheatmap / dendextend | v1.0.13 / v1.19.1 | s03, JmjC |
| tidyr / ggplot2 / ggtree | v1.3.2 / v4.0.2 / v3.14.0 | JmjC |

POSSVM is available at https://github.com/xgrau/possvm-orthology; clone it to `$DIR/JmjC/possvm-orthology-master/`.

Pfam HMM profiles are from Pfam release 37.2. Per-domain `.hmm` files are provided in `JmjC/Pfam/hmm/`; the full database is rebuilt with:

```bash
wget -qO- https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.2/Pfam-A.hmm.gz | gunzip > Pfam-A.hmm
wget -qO- https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.2/Pfam-A.hmm.dat.gz | gunzip > Pfam-A.hmm.dat
hmmpress Pfam-A.hmm
```

`Taxonomy_info.txt` was generated from the NCBI `new_taxdump` (https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/) and is provided.

---

## How to run

```bash
# 1. Data acquisition
bash s01a_download_fasta.sh
bash s01b_download_annotation.sh

# 2. Domain survey across 1,609 species          (~1-2 h)
bash s02_phylogenetic_analysis_all.sh

# 3-5. Matrices, phylogenetic signal, ancestral states
Rscript s03_gene_counts.R
Rscript s04_Dstatistic.R
Rscript s05_ancestral_state_reconstruction.R

# 6-7. Bayesian analysis
bash s06_Bayesian.sh
Rscript s07_bayesian_traits.R

# Supporting: histone demethylase survey          (~2-3 days)
cd JmjC
bash s01_phylogenetic_analysis_rep.sh
Rscript s02_ortholog_counts.R
Rscript s03_conserved_sites.R
```

---

## Contact

**Il-Hwan Lee** (first author, code maintainer) — ianlee.dba@snu.ac.kr
Department of Biomedical Sciences, Seoul National University College of Medicine, Seoul 03080, Republic of Korea

**Hong-Duk Youn** (corresponding author) — hdyoun@snu.ac.kr
Department of Biomedical Sciences, Seoul National University College of Medicine, Seoul 03080, Republic of Korea
