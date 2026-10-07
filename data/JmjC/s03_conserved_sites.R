# ==============================================================================
# Title: Analysis of Conserved Sites for Catalytic Activity in KDM Orthologs
# Author: Il-Hwan LEE
# ==============================================================================
# R 4.4.3 (arm64)
# ==============================================================================
# Load Libraries
library(ggtree) # v3.14.0
library(dplyr) # v1.2.0
library(tidyr) # v1.3.2
library(data.table) # v1.18.2.1
library(ggplot2) # v4.0.2

# ==============================================================================
# Define Functions
# ==============================================================================

get_orthogroup_file <- function(base_dir, ref, domain, group) {
  file.path(base_dir, "orthologs", ref, domain, group, 
            sprintf("%s.%s.orthogroups.uniq.csv", group, domain))
}

get_cof_file <- function(base_dir, ref, domain, group) {
  file.path(base_dir, "orthologs", ref, domain, group, "conserved_sites",
            sprintf("%s.%s.conserved.sites.sorted.csv", group, domain))
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

# ==============================================================================
# Main
# ==============================================================================

Group_list <- c("PEM", "MET")
Domain_list <- c("JmjC")

for (group in Group_list) {
  
  REF <- ifelse(group == "PLT", "HsapAtha", "Hsap")
  
  # Load species_group_profile
  sps_df_path <- file.path(dir_res, sprintf("%s.species_group_profile.csv", group))
  sps_df <- fread(input = sps_df_path, data.table = FALSE, header = T)

  for (domain in Domain_list) {
    
    # --- a. get domain tree & merge with species info ---
    tree_file <- get_tree_file(dir_analysis, REF, domain, group)
    tree <- read.tree(tree_file)
    
    sps_abb <- sapply(strsplit(tree$tip.label, "_"), `[`, 1) # extract species prefix per tip
    
    tip_data <- data.frame(gen_id = tree$tip.label, abbreviation = sps_abb, stringsAsFactors = FALSE)
    tip_data <- merge(tip_data, sps_df, by.x = "abbreviation", by.y = "abbreviation", all.x = T)
    
    # --- b. get conserved sites information & merge with tip_data ---
    cof_file <- get_cof_file(dir_analysis, REF, domain, group)
    cof_df <- fread(input=cof_file, data.table = F, header = F, col.names = c("id", "sites"))
    cof_df$id <- sub("_(?!.*_).*$", "", cof_df$id, perl = TRUE) # remove domain_range info
    
    tip_data$cof_sites <- cof_df$sites[match(tip_data$gen_id, cof_df$id)]
    
    # --- c. match with homology group & sort ---
    ort_fp <- get_orthogroup_file(dir_analysis, REF, domain, group)
    ort <- fread(ort_fp, header = FALSE, select = 1:2, col.names = c("gen_id", "ort_fullname"), data.table = FALSE)
    
    ort <- ort %>% # clean names
      mutate(
        abbreviation = sapply(strsplit(gen_id, "_"), `[`, 1),
        ort_name = sub("^[^:]*:(.*)", "\\1", ort_fullname),
        ort_name_broad = gsub("^like:", "", ort_name)
      )
    
    h_order <- c("KDM2A/KDM2B", "HR/JMJD1C/KDM3A/KDM3B", "KDM4A/KDM4B/KDM4C/KDM4D/KDM4E/KDM4F", 
                   "KDM5A/KDM5B/KDM5C/KDM5D", "KDM6A/KDM6B/UTY", "KDM7A/PHF2/PHF8", "KDM8", 
                   "JMJD6", "JMJD7/JMJD7-PLA2G4B")
    
    tip_data <- merge(tip_data, ort[,-c(2,3)], by.x = "gen_id", by.y = "gen_id", all.x = T) %>% 
      mutate(ort_name_broad = factor(ort_name_broad, levels = h_order)) %>%
      arrange(ort_name_broad)
    
    
    tip_data <- tip_data %>% # add column that states "ort_name_broad | id | sites"
      mutate(name_id_sites = paste(ort_name_broad, " | ", gen_id, " | ", cof_sites))
    
    # --- d. save CSVs, Tree Plots ---
    
    # --- d-0. Colors ---
    tax_colors <- c("deepskyblue1", "goldenrod1", "firebrick2", "darkred")
    gs_colors <- c("deepskyblue1", "magenta1", "chartreuse1")

    tax_color_palette <- setNames(tax_colors[1:(length(unique(sps_df$tax_group)))], c(unique(sps_df$tax_group)))
    gs_color_palette  <- setNames(gs_colors[1:(length(unique(sps_df$gs_group)))], c(unique(sps_df$gs_group)))
    
    # --- d-1. CSVs ---
    write.csv(tip_data, file.path(dir_res, sprintf("%s.%s.tree_conservedsites_tipdata.csv", group, domain)), quote = F, row.names = F)
    
    # --- d-2. Tree Plots ---
    p_tree <- ggtree(tree, layout = "rectangular") %<+% tip_data + 
      geom_tiplab(aes(label = name_id_sites, color = gs_group), size = 1.5, hjust = -0.2, key_glyph = "point") +
      geom_point2(aes(subset = isTip, fill = tax_group), shape = 21, size = 2, stroke = 0.1, alpha = 0.5) +
      scale_color_manual(values = gs_color_palette, breaks = names(gs_color_palette), na.translate = FALSE, na.value = "transparent") +
      scale_fill_manual(values = tax_color_palette, breaks = names(tax_color_palette), na.translate = FALSE, na.value = "transparent") +
      theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top")) +
      guides(color = guide_legend(override.aes = list(size = 3)), fill = guide_legend(override.aes = list(size = 3))) +
      xlim(0, 10)
    
    pdf(file.path(dir_res, sprintf("%s.%s.tree_conservedsites.pdf", group, domain)), width = 10, height = 40)
    print(p_tree) # ignore "In fortify" WARNINGS
    dev.off()
  }
}
  

