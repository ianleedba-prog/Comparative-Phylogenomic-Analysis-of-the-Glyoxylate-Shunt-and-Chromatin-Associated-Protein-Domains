# ==============================================================================
# R 4.4.3 (arm64)
# ==============================================================================
# Load Libraries
library(dplyr) # v1.2.0
library(data.table) # v1.18.2.1
library(pheatmap) # v1.0.13
library(dendextend) #v1.19.1
library(ape) #v5.8-1

# ==============================================================================
# Define functions
# ==============================================================================

# 1. Function to modify uncertain tax levels-written as incertae sedis
fix_incertae_by_phylum <- function(df) {
  tax_levels <- c("Species", "Genus", "Family", "Order", "Class")
  
  df_final <- df %>%
    mutate(across(all_of(tax_levels), 
                  ~ ifelse(tolower(trimws(as.character(.))) == "incertae sedis",
                           paste0("Uncertain_", cur_column(), "_", as.character(Phylum)), 
                           as.character(.))))
  
  return(df_final)
}

# 2. Function to calculate jaccard index
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

# 3. Function to build a species by protein family presence matrix
presence_matrix <- function(gen_tab) {
  mat <- xtabs(formula = ~ Abbreviation + Protein_family, data = gen_tab, drop.unused.levels = F)
  as.data.frame.matrix(mat) > 0
}

# 4. Function to score every species for the glyoxylate shunt and the Tet_JBP characters
#    (characters are returned as 1 for present and 0 for absent)
score_traits <- function(gen_tab, tax_tab, gs_doms, tet_dom, acc_doms) {
  pres <- presence_matrix(gen_tab)
  ix   <- match(rownames(pres), tax_tab$Abbreviation)
  
  data.frame(
    Abbreviation = rownames(pres),
    Kingdom      = tax_tab$Kingdom[ix],
    Clade        = tax_tab$Clade[ix],
    Phylum       = tax_tab$Phylum[ix],
    GS           = as.integer(rowSums(pres[, gs_doms, drop = F]) == length(gs_doms)),
    TET          = as.integer(pres[, tet_dom]),
    h_module     = as.integer(pres[, tet_dom] & rowSums(pres[, acc_doms, drop = F]) > 0),
    row.names    = NULL, stringsAsFactors = F
  )
}

# 5. Function to estimate an odds ratio and its confidence interval
#    (Haldane-Anscombe corrected sample estimate with a Woolf log interval,
#     which stays finite when a cell of the table is empty)
odds_ratio_ci <- function(tab, conf = 0.95) {
  v  <- as.numeric(t(tab)) + 0.5 # a, b, c, d with the correction
  or <- (v[1] * v[4]) / (v[2] * v[3])
  se <- sqrt(sum(1 / v))
  c(OR = or, or * exp(c(CI_low = -1, CI_high = 1) * qnorm((1 + conf) / 2) * se))
}

# 6. Function to test GS presence against each predictor within a set of taxa
#    (P values from the two-tailed Fisher's exact test on the uncorrected table)
odds_ratio_table <- function(traits, taxon_sets) {
  
  res <- data.frame()
  
  for (trait in c("h_module", "TET")) {
    for (set_name in names(taxon_sets)) {
      
      sub <- traits[taxon_sets[[set_name]](traits), ]
      trait_label <- ifelse(trait == "h_module", "h-module", "TET")
      
      # no fungus in the collection carries an accessory domain,
      # so the h-module test is not applied to Holomycota
      if (trait == "h_module" && set_name == "Holomycota") {
        res <- rbind(res, data.frame(Trait = trait_label, Clade = set_name, n = nrow(sub),
                                     OR = NA, CI_low = NA, CI_high = NA, P = NA))
        next
      }
      
      # 2 x 2 table: rows GS present/absent, columns predictor present/absent
      tab <- table(factor(sub$GS,       levels = c(1, 0)),
                   factor(sub[[trait]], levels = c(1, 0)))
      
      est <- odds_ratio_ci(tab)
      
      res <- rbind(res, data.frame(Trait = trait_label, Clade = set_name, n = nrow(sub),
                                   OR = est[["OR"]], CI_low = est[["CI_low"]], CI_high = est[["CI_high"]],
                                   P = fisher.test(tab)$p.value))
    }
  }
  res
}

# 7. Function to combine the results of the two collections into the reported table
format_odds_ratio_table <- function(or_a, or_b, label_a = "1609", label_b = "259") {
  
  fmt_or <- function(or, lo, hi) ifelse(is.na(or), "", sprintf("%.4f [%.4f, %.4f]", or, lo, hi))
  fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-4, formatC(p, format = "e", digits = 2),
                                                   trimws(formatC(p, format = "g", digits = 3))))
  
  out <- data.frame(Trait = or_a$Trait, Clade = or_a$Clade,
                    or_a$n, fmt_or(or_a$OR, or_a$CI_low, or_a$CI_high), fmt_p(or_a$P),
                    or_b$n, fmt_or(or_b$OR, or_b$CI_low, or_b$CI_high), fmt_p(or_b$P),
                    check.names = F, stringsAsFactors = F)
  
  colnames(out) <- c("Trait", "Clade",
                     sprintf("n_%s", label_a), sprintf("OR_%s [95%% CI]", label_a), sprintf("P_%s", label_a),
                     sprintf("n_%s", label_b), sprintf("OR_%s [95%% CI]", label_b), sprintf("P_%s", label_b))
  out
}

# 8. Function to tabulate metazoan species by phylum and combined character state
state_crosstab <- function(traits, phylum_order, state_levels) {
  mz <- traits[traits$Kingdom == "Metazoa", ]
  
  gs  <- ifelse(mz$GS == 1, "GS+", "GS-")
  tet <- ifelse(mz$h_module == 1, "h-module+", ifelse(mz$TET == 1, "TET-only", "TET-"))
  
  mz$State  <- factor(paste(gs, tet, sep = "/"), levels = state_levels)
  mz$Phylum <- factor(mz$Phylum, levels = phylum_order)
  
  as.data.frame.matrix(xtabs(formula = ~ Phylum + State, data = mz, drop.unused.levels = F))
}

# 9. Function to compress a count axis beyond a break point
#    (counts up to brk are linear; beyond it they are compressed by cmp)
xmap <- function(v, brk, cmp) ifelse(abs(v) <= brk, v, sign(v) * (brk + (abs(v) - brk) / cmp))

# 10. Function to draw the odds ratio forest plot
plot_odds_ratio <- function(or_res, fp, main_lab, set_order, trait_pch, trait_col) {
  
  or_res <- or_res[!is.na(or_res$OR), ]
  or_res$y <- match(or_res$Clade, set_order) + ifelse(or_res$Trait == "TET", -0.14, 0.14) # Tet_JBP above the h-module
  
  xlim_log <- c(-3.3, 1.3) # 10^-3 to 10^1, as in the published axis
  
  pdf(file = fp, height = 3.6, width = 3.9, pointsize = 7)
  par(mar = c(5.4, 6.8, 1.2, 0.8), mgp = c(1.9, 0.5, 0), tcl = -0.25, las = 1, xpd = NA)
  plot(NA, xlim = xlim_log, ylim = c(length(set_order) + 0.6, 0.4),
       axes = F, xlab = "Odds ratio", ylab = "", main = main_lab, cex.main = 1)
  
    ## axes
    at_log <- -3:1
    axis(1, at = at_log, labels = parse(text = sprintf("10^%d", at_log)), lwd = 0, lwd.ticks = 0.6)
    axis(1, at = xlim_log, labels = F, lwd.ticks = 0, lwd = 0.6)
    axis(2, at = seq_along(set_order), labels = set_order, lwd = 0, lwd.ticks = 0, cex.axis = 1)
  
    ## reference line at an odds ratio of 1
    segments(0, 0.4, 0, length(set_order) + 0.6, lty = 2, lwd = 0.6, col = "grey50")
  
    ## confidence intervals and point estimates
    segments(log10(or_res$CI_low), or_res$y, log10(or_res$CI_high), or_res$y,
             lwd = 0.8, col = trait_col[or_res$Trait])
    points(log10(or_res$OR), or_res$y, pch = trait_pch[or_res$Trait], cex = 0.8,
           bg = trait_col[or_res$Trait], col = trait_col[or_res$Trait])
  
    ## separator: every comparison above it is significant, the set below it is not
    segments(xlim_log[1], length(set_order) - 0.5, xlim_log[2], length(set_order) - 0.5,
             lty = 3, lwd = 0.6, col = "grey50")
  
    ## legend, placed under the axis in user coordinates (y increases downwards here)
    usr <- par("usr")
    legend(usr[1] - 0.22 * diff(usr[1:2]), usr[3] + 0.55,
           bty = "n", cex = 0.93, pt.cex = 0.8, y.intersp = 1.1,
           legend = c("Tet_JBP", "h-module+ (Tet_JBP+) & (other domain >= 1)"),
           pch = trait_pch[c("TET", "h-module")],
           pt.bg = trait_col[c("TET", "h-module")], col = trait_col[c("TET", "h-module")])
  dev.off()
}

# 11. Function to draw the metazoan state composition chart
#     (GS-positive states extend to the left of zero, GS-negative states to the right)
plot_state_composition <- function(state_tab, fp, x_ticks, main_lab, state_levels, state_col, brk, cmp) {
  
  state_side <- setNames(ifelse(grepl("^GS\\+", state_levels), -1, 1), state_levels)
  left_states  <- state_levels[state_side == -1]
  right_states <- state_levels[state_side ==  1]
  
  n_phy <- nrow(state_tab); bar_h <- 0.7
  y     <- seq_len(n_phy)
  
  xlim <- range(xmap(c(-max(rowSums(state_tab[, left_states,  drop = F])),
                        max(rowSums(state_tab[, right_states, drop = F])), 0), brk, cmp)) * 1.04
  
  pdf(file = fp, height = 5.1, width = 5.6, pointsize = 7)
  par(mar = c(6.9, 7.4, 1.2, 1.4), mgp = c(2.4, 0.5, 0), tcl = -0.25, las = 1, xpd = NA)
  plot(NA, xlim = xlim, ylim = c(n_phy + 0.6, 0.4), axes = F,
       xlab = "Number of species", ylab = "", main = main_lab, cex.main = 1)
  
    ## stacked bars, drawn outward from zero on each side
    for (side_states in list(left_states, right_states)) {
      for (i in y) {
        cum <- 0
        for (st in side_states) {
          n <- state_tab[i, st]
          if (n == 0) next
          rect(xmap(state_side[st] * cum, brk, cmp), y[i] - bar_h/2,
               xmap(state_side[st] * (cum + n), brk, cmp), y[i] + bar_h/2,
               col = state_col[st], border = NA)
          cum <- cum + n
        }
        ## total on each side, printed beyond the end of the stack
        if (cum > 0) {
          side <- state_side[side_states[1]]
          text(xmap(side * cum, brk, cmp) + side * 0.5, y[i], cum,
               adj = c(ifelse(side < 0, 1, 0), 0.5), cex = 0.9)
        }
      }
    }
  
    ## axes
    axis(1, at = xmap(x_ticks, brk, cmp), labels = abs(x_ticks), lwd = 0, lwd.ticks = 0.6)
    axis(1, at = xlim, labels = F, lwd.ticks = 0, lwd = 0.6)
    axis(2, at = y, labels = rownames(state_tab), lwd = 0, lwd.ticks = 0, cex.axis = 1)
    segments(0, 0.4, 0, n_phy + 0.6, lwd = 0.6, col = "grey30")
  
    ## legend
    usr <- par("usr")
    legend(usr[1] - 0.13 * diff(usr[1:2]), usr[3] + 1.7,
           bty = "n", ncol = 2, cex = 0.93, y.intersp = 1.1,
           legend = state_levels, fill = state_col[state_levels], border = NA)
  dev.off()
}

# ==============================================================================
# Environment & Data Setup
# ==============================================================================

# Setup directories
dir <- "/directory/to/data"
dir_analysis <- file.path(dir, "analysis")
dir_res <- file.path(dir_analysis, "results")
if(!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE, showWarnings = FALSE)

# Domain list
dom_list_fp <- file.path(dir, "protein_domains_hmm.csv")
dom_list <- read.table(dom_list_fp, header = F, stringsAsFactors = T)
colnames(dom_list) <- c("Class","Type","Family","Domains","homology_cluster","inflation","min_phylo_size")
dom_list_gap_ixs <- c(1,1+which(diff(as.numeric(dom_list$Type))!=0)) - 1

# Taxonomy information (taxonomy matched from new_taxdump.tar.gz)
sps_tax_fp <- file.path(dir, "Taxonomy_info.txt")
sps_tax <- read.table(sps_tax_fp, sep="\t", quote="\"", fill = T, header = T, comment.char = "", stringsAsFactors = F)
sps_tax <- fix_incertae_by_phylum(sps_tax) # fix uncertain information (Incertae sedis)

sps_tax <- sps_tax %>% # add Clade column
  mutate(Clade = case_when(
    Kingdom %in% c("Metazoa", "Premetazoa") ~ "Holozoa",
    Kingdom %in% c("Fungi", "Rotosphaerida") ~ "Holomycota",
    Kingdom %in% c("Viridiplantae", "Rhodophyta") ~ "Plantae",
    Kingdom %in% c("Malawimonada", "Ancyromonadida") ~ "Mal_Anc",
    Kingdom %in% c("Chromista") ~ "Chromista",
    Kingdom %in% c("Jakobea", "Heterolobosea", "Euglenozoa", "Parabasalia") ~ "Excavata",
    TRUE ~ Kingdom
  )) %>%
  relocate(Clade, .after = Kingdom)

# Gene counts
gen <- fread(input = sprintf("%s/gene_counts/1609.genecounts.csv", dir_analysis), data.table = F, header = F, col.names = c("Protein_id", "Protein_family"))
gen$Abbreviation <- sub("_.*", "", gen[, 1])
gen$Protein_family <- factor(gen$Protein_family, levels = unique(dom_list$Family))

# ==============================================================================
# Main - 1. Gene counts
# ==============================================================================

# Match gene counts with species taxonomy
gen_m <- merge(gen, sps_tax, by.x = "Abbreviation", by.y = "Abbreviation", all.x = T)

# Convert to gene presence table (one entry per species)
gen_u <- within(gen_m, rm("Protein_id"))
gen_u <- unique(gen_u)

# 1-1. Create "counts" matrix (Phylum by Protein family)
phylum_list <- unique(data.frame(sps_tax$Phylum, sps_tax$Kingdom))
phylum_order <- phylum_list[,1]
phylum_gap_ixs <- c(1,1+which(diff(as.numeric(factor(phylum_list$sps_tax.Kingdom)))!=0)) - 1

gen_crosstab <- xtabs(formula = ~ Phylum + Protein_family, data = gen_u, drop.unused.levels	= F)
gen_crosstab <- as.data.frame.matrix(gen_crosstab)
gen_crosstab <- gen_crosstab[phylum_order, ] # row sort by phylum order

# 1-2. Create "fraction" matrix
sps_counts_by_phylum <- table(factor(sps_tax$Phylum, levels=unique(sps_tax$Phylum))) # number of species per phylum
gen_crosstab_frac <- sweep(gen_crosstab, MARGIN = 1, sps_counts_by_phylum, "/")

# 1-3. Column (protein family) sort by clustered pattern

TCA_names <- c("Citrate_synt", "Aconitase", "Iso_dh", "2-oxogl_dehyd", "FumaraseC_C")
GS_names <- c("Malate_synthase", "ICL")

TCA_pos <- match(TCA_names, names(gen_crosstab_frac)) # 1. fixed positions for TCA enzymes
GS_pos <- match(GS_names, names(gen_crosstab_frac)) # 2. fixed positions for GS enzymes

  ## --- a-1. cluster protein families ---
  doms_to_cluster <- gen_crosstab_frac[, setdiff(names(gen_crosstab_frac), c(TCA_names, GS_names))]
  
  dist_cols <- dist(t(doms_to_cluster))
  hc_cols <- hclust(dist_cols)
  
  ## --- a-2. generate dendogram to determine clusters ---
  dend <- as.dendrogram(hc_cols) # dendogram
  node_leaves <- partition_leaves(dend)
  node_pos <- get_nodes_xy(dend)
  nn <- nnodes(dend)
  is_leaf <- is.leaf(dend)
  
  plot(dend)
  text(node_pos, 
       labels = 1:nn, 
       adj = c(1.2, 1.2), 
       col = "blue", cex = 0.7)
  
  cluster_1 <- c(11,13,15,16,9,6,7,3)
  cluster_2 <- c(64,66,68,69)
  cluster_3 <- c(61,62,59,58,55)
  cluster_4 <- c(21,23,25,26)
  cluster_5 <- c(28,31,32,34,36,38,40,42,44,46,48,50,52,53)
  
  target_nodes <- c(cluster_1, cluster_2, cluster_3, cluster_4, cluster_5)
  target_names <- unlist(node_leaves[target_nodes])

reordered_pos <- match(target_names, names(gen_crosstab_frac)) # 3. reordered positions for other protein families

gen_crosstab <- gen_crosstab[,c(TCA_pos, GS_pos, reordered_pos)] # reorder matrix
gen_crosstab_frac <- gen_crosstab_frac[,c(TCA_pos, GS_pos, reordered_pos)] # reorder matrix

    # Save matrix (CSVs)
    write.csv(gen_crosstab, file=sprintf("%s/Count_matrix_by_phylum_1609.csv", dir_res), row.names = T, quote = F)
    write.csv(gen_crosstab_frac, file=sprintf("%s/Fraction_matrix_by_phylum_1609.csv", dir_res), row.names = T, quote = F)

# 1-4. Create per-species "counts" matrix (Species by Protein family, same column order as above)
gen_m$Abbreviation <- factor(gen_m$Abbreviation, levels = sps_tax$Abbreviation) # row order follows the taxonomy table

sps_crosstab <- xtabs(formula = ~ Abbreviation + Protein_family, data = gen_m, drop.unused.levels = F)
sps_crosstab <- as.data.frame.matrix(sps_crosstab)
sps_crosstab <- sps_crosstab[, colnames(gen_crosstab)] # column sort as in the phylum-level matrix

sps_count_tab <- cbind(sps_tax[match(rownames(sps_crosstab), sps_tax$Abbreviation),
                               c("Abbreviation", "Species", "Kingdom", "Phylum", "Class")],
                       sps_crosstab)

    # Save matrix (CSV)
    write.csv(sps_count_tab, file=sprintf("%s/Count_matrix_by_species_1609.csv", dir_res), row.names = F, quote = F)

# 1-5. Heatmap
col_heatmap <- colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue2","dodgerblue4"))

pdf(file=sprintf("%s/Fraction_matrix_by_phylum_heatmap.pdf", dir_res), height=10,width=10)
pheatmap(gen_crosstab_frac, color = col_heatmap(20), breaks = c(0,0.0001,seq(0.05,1,length.out = 19)),
         gaps_col = c(length(TCA_pos), length(c(TCA_pos, GS_pos))),
         gaps_row = phylum_gap_ixs,
         cellwidth = 5, cellheight = 5, fontsize = 5,
         border_color = "white", cluster_cols=F, cluster_rows=F, display_numbers = F)
dev.off()

# Save number of species per phylum (CSVs, Bar Graph)
write.csv(sps_counts_by_phylum, file=sprintf("%s/Species_count_per_phylum.csv", dir_res), row.names = F, quote = F)

pdf(file=sprintf("%s/Species_count_per_group_1609.pdf", dir_res), height=5,width=10)
barplot(sps_counts_by_phylum, 
        las=2, border = NA, cex.names = 0.5,
        main="# of species per phylum")
dev.off()

# ==============================================================================
# Main - 2. Jaccard index
# ==============================================================================

tax_rank <- c("Phylum", "Kingdom", "Clade", "Superkingdom")
base_domains <- c("ICL", "Tet_JBP")
target_list <- as.character(dom_list$Family)

for (rank in tax_rank) {
  name_list = as.character(unique(sps_tax[,rank]))
  
  # matrix to store result
  target_mat_list <- list(
    ICL = matrix(ncol = length(target_list) - 1, nrow = 0), # exclude itself -1
    Tet_JBP = matrix(ncol = length(target_list) - 1, nrow = 0)
  )
  rowname_list <- list(ICL = c(), Tet_JBP = c())
  
  for (name in name_list) {
    sub_sps_tax <- sps_tax[sps_tax[,rank] == name, ]
    sub_gen_m <- gen_m[gen_m[,rank] == name, ]
    
    if (nrow(sub_gen_m) == 0) next
    
    sub_gen_m$Taxonomy_name_factor <- factor(sub_gen_m$Taxonomy_name, levels = unique(sub_sps_tax$Taxonomy_name))
    sub_gen_u <- unique(within(sub_gen_m, rm("Protein_id")))
    
    gen_crosstab_jacc <- xtabs(formula = ~ Protein_family + Taxonomy_name_factor, data = sub_gen_u, drop.unused.levels = F)
    gen_crosstab_jacc[gen_crosstab_jacc > 1] <- 1
    
    # calculate jaccard overlap score
    for (base_dom in base_domains) {
      
      # remove base_domain from target list
      curr_targets <- setdiff(target_list, base_dom) 
      
      if (base_dom %in% rownames(gen_crosstab_jacc)) {
        tax_with_base <- colnames(gen_crosstab_jacc)[ gen_crosstab_jacc[base_dom, ] > 0 ]
        
        jac_scores <- unlist(lapply(curr_targets, function(target) { 
          if (target %in% rownames(gen_crosstab_jacc)) {
            tax_with_target <- colnames(gen_crosstab_jacc)[ gen_crosstab_jacc[target, ] > 0 ]
            return(jaccard(tax_with_base, tax_with_target))
          } else {
            return(0)
          }
        }))
        names(jac_scores) <- curr_targets
        
        # store results
        rowname_list[[base_dom]] <- c(rowname_list[[base_dom]], name)
        target_mat_list[[base_dom]] <- rbind(target_mat_list[[base_dom]], jac_scores)
      }
    }
  }
  
  # Save results (CSVs, Bar Graphs)
  for (base_dom in base_domains) {
    
    mat <- target_mat_list[[base_dom]]
    
    if(nrow(mat) == 0) next
    
    rownames(mat) <- rowname_list[[base_dom]]
    colnames(mat) <- setdiff(target_list, base_dom)
    mat <- na.omit(mat)
    
    write.table(mat, file=sprintf("%s/Jaccard_index_%s_%s_matrix.csv", dir_res, rank, base_dom), row.names = T, col.names = NA, quote = F, sep = "\t")
    
    # Barplots
    pdf(file=sprintf("%s/Jaccard_index_%s_%s_barplot.pdf", dir_res, rank, base_dom), width = 4, height = 6)
    for (tax_name in row.names(mat)) {
      
      sub_vec <- mat[tax_name, ]
      sub_vec <- sort(sub_vec, decreasing = FALSE) 
      
      par(mgp=c(1,0.5,0), mar=c(5,4,4,2)+0.1)
      barplot(sub_vec, horiz = T, las=2, col = rainbow(n = 1, start = 0.1, end = 0.85, v = 0.9), cex.names=0.3, cex.axis = 0.6, xlab = "Jaccard index", cex.lab=0.6)
      title(main=sprintf("Jaccard overlap scores (%s)\nSimilarity between the %s-level distribution of %s with other domains", tax_name, rank, base_dom), cex.main=0.6)
    }
    dev.off()
  }
}

# ==============================================================================
# Main - 3. Subset species (n = 259)
# ==============================================================================

# Taxonomy information of "representative 259 species" across eukaryotes
rep_tree <- read.tree(file.path(dir,"Eukaryota_259.nwk"))
rep_sps <- gsub("lucimarinus", "'lucimarinus'", gsub("_", " ", rep_tree$tip.label))

# Subset species (n = 259) from sps_tax, gen_m
rep_tax <- sps_tax[match(rep_sps, sps_tax$Species), ]
rep_gen_m <- gen_m[gen_m$Species %in% rep_sps, ]

# 3-1. Plot tree (rectangular phylogram with branch lengths, abbreviated tip labels, scale bar)
rep_tree_plot <- rotateConstr(rep_tree, rev(rep_tree$tip.label)) # newick order runs from top to bottom
rep_sps_abb <- setNames(sps_tax$Abbreviation[match(rep_sps, sps_tax$Species)], rep_tree$tip.label)
rep_tree_plot$tip.label <- unname(rep_sps_abb[rep_tree_plot$tip.label]) # tip labels as species abbreviations

pdf(file=sprintf("%s/Eukaryota_259.newick.tree.pdf", dir_res), height=20, width=8)
par(xpd = NA)
plot.phylo(rep_tree_plot, type = "phylogram", direction = "rightwards",
           show.tip.label = TRUE, label.offset = 0.02, cex = 0.4, font = 3,
           edge.width = 0.4, no.margin = TRUE)

  ## scale bar (expected domain changes per branch)
  scale_bar <- 0.5
  usr <- par("usr")
  x0 <- usr[1] + 0.02 * diff(usr[1:2]); y0 <- usr[3] + 0.01 * diff(usr[3:4])
  segments(x0, y0, x0 + scale_bar, y0, lwd = 1)
  text(x0, y0 - 0.008 * diff(usr[3:4]), sprintf("%.1f expected domain changes", scale_bar), adj = c(0, 1), cex = 0.6)
dev.off()

# 3-2. Create per-species "counts" matrix (same format as the 1609-species matrix)
rep_gen_m$Abbreviation <- factor(as.character(rep_gen_m$Abbreviation), levels = rep_tax$Abbreviation)

rep_sps_crosstab <- xtabs(formula = ~ Abbreviation + Protein_family, data = rep_gen_m, drop.unused.levels = F)
rep_sps_crosstab <- as.data.frame.matrix(rep_sps_crosstab)
rep_sps_crosstab <- rep_sps_crosstab[, colnames(gen_crosstab)]

rep_sps_count_tab <- cbind(rep_tax[match(rownames(rep_sps_crosstab), rep_tax$Abbreviation),
                                   c("Abbreviation", "Species", "Kingdom", "Phylum", "Class")],
                           rep_sps_crosstab)
    # Save matrix (CSV)
    write.csv(rep_sps_count_tab, file=sprintf("%s/Count_matrix_by_species_259.csv", dir_res), row.names = F, quote = F)

# 3-3. Number of species per group (259)
  ## plotting order
  pool <- c(Dinoflagellata="Other Alveolata", Chromerida="Other Alveolata",
          Bacillariophyta="Other Stramenopiles", Ochrophyta="Other Stramenopiles", Bigyra="Other Stramenopiles",
          Cercozoa="Rhizaria", Endomyxa="Rhizaria",
          Haptophyta="Hacrobia", Cryptophyceae="Hacrobia")
  rep_tax$Group_tmp <- ifelse(rep_tax$Phylum %in% names(pool), pool[rep_tax$Phylum], rep_tax$Phylum)

  layout_grp <- list(
    "Metazoa"        = c("Chordata","Hemichordata","Echinodermata","Arthropoda","Tardigrada","Nematoda","Priapulida",
                         "Mollusca","Annelida","Brachiopoda","Platyhelminthes","Rotifera","Acanthocephala",
                         "Xenacoelomorpha","Cnidaria","Placozoa","Ctenophora","Porifera"),
    "Premetazoa"     = c("Choanoflagellata","Filasterea","Corallochytrea","Ichthyosporea"),
    "Holomycota"     = c("Basidiomycota","Ascomycota","Mucoromycota","Zoopagomycota","Blastocladiomycota",
                         "Chytridiomycota","Rotosphaerida"),
    " "              = c("Apusomonadida"),
    "Amoebozoa"      = c("Evosea","Discosea"),
    "  "             = c("Malawimonada","Ancyromonadida"),
    "Plantae"        = c("Streptophyta","Chlorophyta","Rhodophyta"),
    "SAR"            = c("Ciliophora","Apicomplexa","Other Alveolata","Oomycota","Other Stramenopiles","Rhizaria"),
    "   "            = c("Hacrobia"),
    "Discoba"        = c("Jakobea","Heterolobosea","Euglenozoa"),
    "    "           = c("Parabasalia"))
  
  rep_grp_order <- unlist(layout_grp, use.names = F)
  rep_sps_counts <- table(factor(rep_tax$Group_tmp, levels = rep_grp_order))

  ## Bar plot
  n_grp <- length(rep_grp_order); bar_w <- 0.8
  cnt_v <- as.numeric(rep_sps_counts)
  med <- median(cnt_v)
  
  pdf(file=sprintf("%s/Species_count_per_group_259.pdf", dir_res), height=3.6, width=7.2, pointsize=7)
  par(mar = c(9.5, 3, 0.6, 9.5), mgp = c(1.8, 0.5, 0), tcl = -0.25, las = 1, xpd = NA)
  plot(NA, xlim = c(-0.9, n_grp - 0.4), ylim = c(0, 31.5), xaxs = "i", yaxs = "i",
       axes = F, xlab = "", ylab = "Number of species")
  axis(2, at = seq(0, 30, 5), lwd = 0, lwd.ticks = 0.6)
  axis(2, at = c(0, 31.5), labels = F, lwd.ticks = 0, lwd = 0.6)
  segments(-0.9, 0, n_grp - 0.4, 0, lwd = 0.6) # baseline
  
    ## bars and the count above each of them
    rect(0:(n_grp-1) - bar_w/2, 0, 0:(n_grp-1) + bar_w/2, cnt_v, col = "grey83", border = NA)
    text(0:(n_grp-1), cnt_v + 0.3, cnt_v, adj = c(0.5, 0), cex = 0.93)
  
    ## median
    abline(h = med, lty = 2, lwd = 0.6, col = "blue")
    text(n_grp - 0.3, med, sprintf("Median number of\nspecies per group (%g)", med),
         adj = c(0, 0.5), col = "blue", cex = 0.93)
  
    ## group names
    text(0:(n_grp-1), -0.3, rep_grp_order, srt = 90, adj = c(1, 0.5), cex = 1)
  dev.off()
  
# ==============================================================================
# Main - 4. Odds ratio (1609 & 259) - Opisthokont, Holozoa, Holomycota, Metazoa
# ==============================================================================

# Character definitions
#   GS       : the shunt, coexistence of ICL and Malate_synthase
#   TET      : Tet_JBP alone
#   h-module : Tet_JBP together with at least one of ADD_ATRX, ADD_DNMT3 and TTD
GS_doms  <- c("ICL", "Malate_synthase")
TET_dom  <- "Tet_JBP"
ACC_doms <- c("ADD_ATRX", "ADD_DNMT3", "TTD")

# Taxon sets
taxon_sets <- list(
  Opisthokonta = function(d) d$Clade %in% c("Holozoa", "Holomycota"),
  Holomycota   = function(d) d$Clade == "Holomycota",
  Holozoa      = function(d) d$Clade == "Holozoa",
  Metazoa      = function(d) d$Kingdom == "Metazoa"
)

# 4-1. Score every species for the two characters
sps_traits <- score_traits(gen_m, sps_tax, GS_doms, TET_dom, ACC_doms) # 1,609 species
rep_traits <- sps_traits[match(rep_tax$Abbreviation, sps_traits$Abbreviation), ] # 259 species
rownames(rep_traits) <- NULL

    # Save matrix (CSV)
    write.csv(sps_traits, file=sprintf("%s/Character_states_by_species_1609.csv", dir_res), row.names = F, quote = F)
    write.csv(rep_traits, file=sprintf("%s/Character_states_by_species_259.csv", dir_res), row.names = F, quote = F)

# 4-2. Association between GS presence and each predictor
or_1609 <- odds_ratio_table(sps_traits, taxon_sets)
or_259  <- odds_ratio_table(rep_traits, taxon_sets)

or_tab <- format_odds_ratio_table(or_1609, or_259, label_a = "1609", label_b = "259")

    # Save matrix (CSV)
    write.csv(or_tab, file=sprintf("%s/Odds_ratio_GS_vs_TET_hmodule.csv", dir_res), row.names = F, quote = F)

# 4-3. Forest plot
set_order <- c("Opisthokonta", "Holozoa", "Metazoa", "Holomycota") # top to bottom; Holomycota below the separator
trait_pch <- c("TET" = 21, "h-module" = 22) # circles for Tet_JBP alone, squares for the h-module
trait_col <- c("TET" = "grey25", "h-module" = "dodgerblue4")

plot_odds_ratio(or_1609, sprintf("%s/Odds_ratio_forest_1609.pdf", dir_res), "1,609-species collection",
                set_order, trait_pch, trait_col)
plot_odds_ratio(or_259,  sprintf("%s/Odds_ratio_forest_259.pdf", dir_res), "259-species collection",
                set_order, trait_pch, trait_col)

# ==============================================================================
# Main - 5. Metazoan species state composition (1609 & 259) - GS, h-module
# ==============================================================================

# Six combined states
state_levels <- c("GS-/TET-", "GS-/TET-only", "GS-/h-module+",
                  "GS+/TET-", "GS+/TET-only", "GS+/h-module+")
state_col <- c("GS-/TET-" = "#6E6E6E", "GS-/TET-only" = "#A6A6A6", "GS-/h-module+" = "#D3D3D3",
               "GS+/TET-" = "#0A22A0", "GS+/TET-only" = "#2F63D6", "GS+/h-module+" = "#5DA4FF")

x_break <- 15
x_cmp   <- 8

# 5-1. Species composition per metazoan phylum
mz_phylum_order <- layout_grp[["Metazoa"]] # 18 metazoan phyla, in the order used throughout

state_1609 <- state_crosstab(sps_traits, mz_phylum_order, state_levels)
state_259  <- state_crosstab(rep_traits, mz_phylum_order, state_levels)

    # Save matrices (CSVs)
    write.csv(state_1609, file=sprintf("%s/Metazoan_state_composition_1609.csv", dir_res), row.names = T, quote = F)
    write.csv(state_259,  file=sprintf("%s/Metazoan_state_composition_259.csv", dir_res), row.names = T, quote = F)

# 5-2. Bar plot
plot_state_composition(state_1609, sprintf("%s/Metazoan_state_composition_1609.pdf", dir_res),
                       x_ticks = c(-100, -15, 0, 15, 100, 150, 200, 250),
                       main_lab = "Metazoa (n = 570 species, 18 phyla)",
                       state_levels, state_col, x_break, x_cmp)

plot_state_composition(state_259, sprintf("%s/Metazoan_state_composition_259.pdf", dir_res),
                       x_ticks = c(-15, -10, -5, 0, 5, 10, 15),
                       main_lab = "Metazoa (n = 112 species, 18 phyla)",
                       state_levels, state_col, x_break, x_cmp)
