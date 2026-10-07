# ==============================================================================
# Title: Ortholog Count/Fraction Analysis & Tree visualization
# Author: Il-Hwan LEE
# ==============================================================================
# R 4.4.3 (arm64)
# ==============================================================================
# Load Libraries
library(phytools) # v2.5-2
library(ape) # v5.8-1
library(ggtree) # v3.14.0
library(dplyr) # v1.2.0
library(tidyr) # v1.3.2
library(data.table) # v1.18.2.1
library(ggplot2) # v4.0.2
library(pheatmap) # v1.0.13

# ==============================================================================
# Define Functions
# ==============================================================================

plot_taxonomy_tree <- function(tree_file, xlim_max = 15, label_size = 3, node_label_size = 2) {
  tree <- read.tree(tree_file)
  tip_data <- data.frame(label = tree$tip.label)
  
  tree_plot <- suppressWarnings({
    ggtree(tree, layout = "rectangular", ladderize = FALSE, branch.length = "none", linewidth = 0.5) %<+% tip_data +
      geom_tiplab(aes(label = label), size = label_size) +
      xlim(0, xlim_max)
  })
  
  highest_node <- tree$edge[1, 1]
  tree_plot <- ggtree::rotate(tree_plot, highest_node)
  
  tree_plot <- tree_plot + 
    geom_text2(aes(label = label), size = node_label_size, hjust = 1.5, data = ~subset(.x, !isTip))
  
  return(tree_plot)
}

get_orthogroup_file <- function(base_dir, ref, domain, group) {
  file.path(base_dir, "orthologs", ref, domain, group, 
            sprintf("%s.%s.orthogroups.uniq.csv", group, domain))
}

count_genes_by_subtypes <- function(orthogroups, subtypes, sps_list) {
  # Create a complete grid to ensure NO species or subtype is left out (filled with 0)
  expand_grid(abbreviation = sps_list, ort_name_broad = subtypes) %>%
    left_join(orthogroups, by = c("abbreviation", "ort_name_broad")) %>%
    group_by(abbreviation, ort_name_broad) %>%
    summarise(n = sum(!is.na(sequence_id)), .groups = "drop") %>%
    pivot_wider(names_from = ort_name_broad, values_from = n, values_fill = list(n = 0)) %>%
    arrange(match(abbreviation, sps_list)) %>%
    tibble::column_to_rownames("abbreviation") %>%
    as.matrix()
}

calculate_subtype_fraction_by_taxonomy <- function(subtype_mat, subgroup_list) {
  frac_mat <- matrix(nrow = length(subgroup_list), ncol = ncol(subtype_mat))
  colnames(frac_mat) <- colnames(subtype_mat)
  rownames(frac_mat) <- names(subgroup_list)
  
  for (i in seq_along(subgroup_list)) {
    group_sps <- subgroup_list[[i]]
    # Intersect to avoid subscript out of bounds errors
    valid_sps <- intersect(group_sps, rownames(subtype_mat)) 
    
    if (length(valid_sps) > 0) {
      frac_mat_tmp <- subtype_mat[valid_sps, , drop = FALSE]
      frac_mat_tmp[frac_mat_tmp > 0] <- 1
      frac_mat[i, ] <- colSums(frac_mat_tmp) / length(group_sps)
    } else {
      frac_mat[i, ] <- 0
    }
  }
  return(frac_mat)
}

plot_count <- function(data, color_vec, palette, main_title, xlab_text) {
  for (i in 1:ncol(data)) {
    
    y_max <- max(data[, i], na.rm=TRUE)
    y_limit <- if(y_max == 0) 1 else y_max * 1.25
    
    bp <- barplot(data[, i], col = "lightgrey", border = color_vec, las = 1, 
                  ylim = c(0, y_limit),
                  main = sprintf("%s | %s: gene counts", group, colnames(data)[i]),
                  xlab = xlab_text, ylab = "Number of genes", names.arg = rep("", nrow(data)))
    mtext(sps_df$abbreviation, side = 1, at = bp, col = color_vec, line = 1, cex = 0.7, las = 2)
    
    legend("topright", 
           legend = names(palette),      
           fill = "lightgrey",           
           border = unname(palette),     
           bty = "n",                    
           cex = 0.8)
  }
}

plot_frac <- function(mat, colors) {
  for (i in 1:nrow(mat)) {
    bp <- barplot(mat[i, ], col = "lightgrey", border = colors[i+1], las = 2, ylim = c(0, 1.2), 
                  main = sprintf("%s | %s: fraction of taxa with domains", group, rownames(mat)[i]),
                  ylab = "# of Taxa / Total # of taxa", names.arg = rep("", ncol(mat)))
    text(x = bp, y = par("usr")[3] - 0.02, labels = colnames(mat), srt = 45, adj = 1, xpd = TRUE, cex = 0.6)
    
    legend("topright", 
           legend = rownames(mat)[i], 
           fill = "lightgrey", 
           border = colors[i+1], 
           bty = "n", 
           cex = 0.8)
  }
}

get_tree_file <- function(base_dir, ref, domain, group) {
  file.path(base_dir, "orthologs", ref, domain, group,
            sprintf("%s.%s.iqtree_treeshrink", group, domain),
            "output.treefile")
}

# ==============================================================================
# Environment & Data Setup
# ==============================================================================

# Directories
dir <- "/directory/to/data/JmjC"
dir_analysis <- file.path(dir, "analysis")
dir_res <- file.path(dir_analysis, "results")
if(!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE, showWarnings = FALSE)

# Load Gene counts
gen_path <- file.path(dir_analysis, "gene_counts", "Rep.genecounts.csv")
gen <- fread(input = gen_path, data.table = FALSE, header = FALSE, col.names = c("Protein_id", "Protein_family"))

# Filter Glyoxylate Shunt (GS) related genes
gen_gs <- gen %>% 
  filter(Protein_family %in% c("ICL", "Malate_synthase")) %>%
  mutate(Abbreviation = sapply(strsplit(Protein_id, "_"), `[`, 1))

# Identify species with complete GS
species_gs <- gen_gs %>%
  group_by(Abbreviation) %>%
  filter(all(c("ICL", "Malate_synthase") %in% Protein_family)) %>%
  pull(Abbreviation) %>%
  unique()

# ==============================================================================
# Main
# ==============================================================================

Group_list <- c("PEM", "MET", "PLT")
Domain_list <- c("Tet_JBP", "JmjC", "AOD")

for (group in Group_list) {
  
  REF <- ifelse(group == "PLT", "HsapAtha", "Hsap")
  
  # Load Species List
  sps_file <- file.path(dir, sprintf("Sps.prefix_%s.txt", group))
  sps_list <- read.table(sps_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
  
  # --- a. Plot tree ---
  tree_file <- file.path(dir, sprintf("Sps.tree_%s.newick.txt", group))
  tree_plot <- plot_taxonomy_tree(tree_file)
  
  ggsave(filename = file.path(dir_res, sprintf("Sps.tree_%s.newick.pdf", group)), 
         plot = tree_plot, height = 7, width = 5)
  
  for (domain in Domain_list) {
    
    # --- b. Species division ---
    if (group == "PEM") {
      sps_grp <- list(
        Earlymetazoa = sps_list[which(sps_list == "Tadh"):which(sps_list == "Lcom")],
        Premetazoa   = sps_list[which(sps_list == "Mbre"):which(sps_list == "Clim")]
      )
    } else if (group == "MET") {
      sps_grp <- list(
        Panarthropoda = sps_list[which(sps_list == "Bmor"):which(sps_list == "Rvar")],
        Nematoda      = sps_list[which(sps_list == "Cele"):which(sps_list == "Tbri")],
        Spiralia      = sps_list[which(sps_list == "Aplcal"):which(sps_list == "Avag")]
      )
    } else {
      sps_grp <- list(
        Metazoa = sps_list[which(sps_list == "Hsap"):which(sps_list == "Scil")],
        Plantae = sps_list[which(sps_list == "Psat"):which(sps_list == "Ccri")]
      )
    }
    
    sps_gs <- list(
      gs_neg = sps_list[!sps_list %in% species_gs & sps_list != "Hsap"],
      gs_pos = sps_list[sps_list %in% species_gs]
    )
    
    # --- c. Ortholog result dataframe ---
    ort_fp <- get_orthogroup_file(dir_analysis, REF, domain, group)
    ort <- fread(ort_fp, header = FALSE, select = 1:2, col.names = c("sequence_id", "ort_fullname"), data.table = FALSE)
    
    # Clean ortholog names
    ort <- ort %>%
      mutate(
        abbreviation = sapply(strsplit(sequence_id, "_"), `[`, 1),
        ort_name = sub("^[^:]*:(.*)", "\\1", ort_fullname),
        ort_name_broad = gsub("^like:", "", ort_name)
      )
    
    subtypes <- unique(ort$ort_name_broad)
    
    # --- d. Species group profile creation ---
    grp_mapping <- stack(sps_grp) %>% setNames(c("abbreviation", "tax_group")) %>% mutate(abbreviation = as.character(abbreviation))
    gs_mapping <- stack(sps_gs) %>% setNames(c("abbreviation", "gs_group")) %>% mutate(abbreviation = as.character(abbreviation))
    
    sps_df <- data.frame(abbreviation = sps_list, stringsAsFactors = FALSE) %>%
      left_join(grp_mapping, by = "abbreviation") %>%
      left_join(gs_mapping, by = "abbreviation") %>%
      mutate(
        tax_group = replace_na(as.character(tax_group), "Human"),
        gs_group = replace_na(as.character(gs_group), "Human")
      )
    
      ## colors
      tax_colors <- c("deepskyblue1", "goldenrod1", "firebrick2", "darkred")
      gs_colors <- c("deepskyblue1", "magenta1", "chartreuse1")
      
        ## color defined by sps_grp
        if (group == "PLT") {
          tax_color_palette <- setNames(tax_colors[2:(length(sps_grp) + 1)], names(sps_grp))
        } else {
          tax_color_palette <- setNames(tax_colors[1:(length(sps_grp) + 1)], c("Human", names(sps_grp)))
        }
        
        ## color defined by sps_gs
        gs_color_palette  <- setNames(gs_colors[1:(length(sps_gs) + 1)], c("Human", names(sps_gs)))
      
      sps_df <- sps_df %>%
        mutate(
          tax_color = tax_color_palette[tax_group],
          gs_color = gs_color_palette[gs_group]
        )

    # --- e. Counts & Fraction Calculation ---
    subtype_mat <- count_genes_by_subtypes(ort, subtypes, sps_list)
    frac_mat_grp <- calculate_subtype_fraction_by_taxonomy(subtype_mat, sps_grp)
    frac_mat_gs <- calculate_subtype_fraction_by_taxonomy(subtype_mat, sps_gs)
    
    # Remove empty columns
    valid_cols <- colSums(subtype_mat) != 0
    subtype_mat <- subtype_mat[, valid_cols, drop = FALSE]
    
    # Reorder columns
    gap_col_ixs <- NULL
    if (group == "PLT" && domain == "JmjC") {
      cluster_1 <- c("KDM2A/KDM2B", "KDM6A/KDM6B/UTY", "KDM7A/PHF2/PHF8")
      cluster_2 <- c("HR/JMJD1C/KDM3A/KDM3B", "B160/IBM1/JMJ24/JMJ26/JMJ27/JMJ29")
      cluster_3 <- c("KDM4A/KDM4B/KDM4C/KDM4D/KDM4E/KDM4F", "ELF6/JMJ13/REF6")
      cluster_4 <- c("JMJ17/KDM5A/KDM5B/KDM5C/KDM5D", "JMJ14/JMJ18/MEE27", "JMJ19", "PKDM7D", "JMJ14/JMJ18/MEE27/PKDM7D")
      cluster_5 <- c("JMJD5/KDM8", "JMJ21", "JMJ22", "JMJ21/JMJ22/JMJD6", "JMJ31", "JMJD6", "JMJ32/JMJD7/JMJD7-PLA2G4B")
      
      col_order <- c(cluster_1, cluster_2, cluster_3, cluster_4, cluster_5)
      gap_col_ixs <- cumsum(c(length(cluster_1), length(cluster_2), length(cluster_3), length(cluster_4)))
    } else if (group != "PLT" && domain == "JmjC") {
      col_order <- c("KDM2A/KDM2B", "HR/JMJD1C/KDM3A/KDM3B", "KDM4A/KDM4B/KDM4C/KDM4D/KDM4E/KDM4F", 
                     "KDM5A/KDM5B/KDM5C/KDM5D", "KDM6A/KDM6B/UTY", "KDM7A/PHF2/PHF8", "KDM8", 
                     "JMJD6", "JMJD7/JMJD7-PLA2G4B")
    } else {
      col_order <- colnames(subtype_mat)
    }
    
    subtype_mat <- subtype_mat[, col_order, drop = FALSE]
    frac_mat_grp <- frac_mat_grp[, col_order, drop = FALSE]
    frac_mat_gs <- frac_mat_gs[, col_order, drop = FALSE]
    
    # Save CSVs
    write.csv(sps_df, file.path(dir_res, sprintf("%s.species_group_profile.csv", group)), quote = F, row.names = F)
    write.csv(subtype_mat, file.path(dir_res, sprintf("%s.%s.counts.csv", group, domain)), quote = F)
    write.csv(frac_mat_grp, file.path(dir_res, sprintf("%s.%s.fraction_grouped_by_taxinfo.csv", group, domain)), quote = F)
    write.csv(frac_mat_gs, file.path(dir_res, sprintf("%s.%s.fraction.grouped_by_gspresence.csv", group, domain)), quote = F)
    
    # --- f. Plots ---
    
    # Heatmap Gaps
    gap_row_ixs <- if (group == "PLT") {
      length(sps_grp[[1]])
    } else {
      head(cumsum(c(1, sapply(sps_grp, length))), -1)
    }
    
    # --- f-1. Heatmap ---
    pdf(file.path(dir_res, sprintf("%s.%s.counts.heatmap.pdf", group, domain)), height=5, width=5)
    pheatmap(subtype_mat, color = c("gray90", "grey13"), breaks = c(-0.001, 0.001, max(subtype_mat, na.rm=T)),
             gaps_row = gap_row_ixs, gaps_col = gap_col_ixs, cellwidth = 4, cellheight = 4, 
             number_color = "aliceblue", fontsize = 4, border_color = "white", 
             cluster_cols=FALSE, cluster_rows=FALSE, display_numbers = TRUE, legend = FALSE, number_format = "%.0f")
    dev.off()
    
    # --- f-2. Bar Graphs ---
    pdf(file.path(dir_res, sprintf("%s.%s.fraction.bargraph.pdf", group, domain)), height=5, width=8)
    plot_count(subtype_mat, sps_df$tax_color, tax_color_palette, "taxonomy groups", "Colored by taxonomy groups")
    plot_count(subtype_mat, sps_df$gs_color, gs_color_palette, "glyoxylate shunt existence", "Colored by glyoxylate shunt existence")
    plot_frac(frac_mat_grp, tax_colors)
    plot_frac(frac_mat_gs, gs_colors)
    dev.off()
    
    # --- f-3. Tree Plots for KDM1A/B, Tet_JBP --- ! JmjC-containing KDMs will be analyzed in s05_conserved_sites.R
    if (domain == "AOD" | domain == "Tet_JBP") {
      tree_file <- get_tree_file(dir_analysis, REF, domain, group)
      tree <- read.tree(tree_file)
      
      sps_abb <- sapply(strsplit(tree$tip.label, "_"), `[`, 1) # extract species prefix per tip
      
      tip_data <- data.frame(sequence_id = tree$tip.label, abbreviation = sps_abb, stringsAsFactors = FALSE) # dataframe with tip info
      tip_data <- merge(tip_data, sps_df, by.x = "abbreviation", by.y = "abbreviation", all.x = T) # merge with sps_df
      tip_data <- merge(tip_data, ort[,-c(2,3)], by.x = "sequence_id", by.y = "sequence_id", all.x = T) # merge with homology group info
      tip_data <- tip_data %>% # add column that states "sequence_id | homology group"
        mutate(seqid_hg = paste(sequence_id, " | ", ort_name_broad)) 

        ## rectangular tree
        tree_rec <- ggtree(tree, layout = "rectangular") %<+% tip_data +
          geom_tiplab(aes(label = seqid_hg, color = gs_group), size = 1.0, hjust = -0.2, key_glyph = "point") +
          geom_point2(aes(subset = isTip, fill = tax_group), shape = 21, size = 2, stroke = 0.1, alpha = 0.5) +
          scale_color_manual(values = gs_color_palette, breaks = names(gs_color_palette), na.translate = FALSE, na.value = "transparent") +
          scale_fill_manual(values = tax_color_palette, breaks = names(tax_color_palette), na.translate = FALSE, na.value = "transparent") +
          theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top")) +
          guides(color = guide_legend(override.aes = list(size = 3)), fill = guide_legend(override.aes = list(size = 3))) +
          xlim(0, 10)
        
        ## circular tree
        tree_cir <- ggtree(tree, layout = "circular") %<+% tip_data +
          geom_tiplab(aes(label = seqid_hg, color = gs_group), size = 2, hjust = -0.5, key_glyph = "point") +
          geom_point2(aes(subset = isTip, fill = tax_group), shape = 21, size = 2, stroke = 0.1, alpha = 0.5) +
          scale_color_manual(values = gs_color_palette, breaks = names(gs_color_palette), na.translate = FALSE, na.value = "transparent") +
          scale_fill_manual(values = tax_color_palette, breaks = names(tax_color_palette), na.translate = FALSE, na.value = "transparent") +
          theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top")) +
          guides(color = guide_legend(override.aes = list(size = 3)), fill = guide_legend(override.aes = list(size = 3)))
        
      pdf(file.path(dir_res, sprintf("%s.%s.tree.pdf", group, domain)), width = 12, height = 12)
      print(tree_rec) # ignore "In fortify" WARNINGS
      print(tree_cir) # ignore "In fortify" WARNINGS
      dev.off()
    }
  }
}
