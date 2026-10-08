# ==============================================================================
# Title: KEGG Metabolic Modules Investigation with KEGG Taxonomy DB
# Author: Il-Hwan LEE
# ==============================================================================
# R 4.4.3 (arm64)
# ==============================================================================
# Load Libraries
library(data.table) # v1.18.2.1
library(pheatmap) # v1.0.13
library(seriation) # v1.5.8

# ==============================================================================
# Environment & Data Setup
# ==============================================================================

# Setup directories
dir = "/directory/to/data_kegg"
dir_res <- file.path(dir, "results")
if(!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE, showWarnings = FALSE)

# KEGG Taxonomy database
Tax_df = file.path(dir, "KEGG_Taxonomy_list.csv")
Tax_df = read.csv(Tax_df, header = T, stringsAsFactors = T)
Tax_df = Tax_df[!is.na(Tax_df$Abbreviation), ]  #Remove rows without species information

# KEGG metabolic module list
Module_list = file.path(dir, "rest.kegg.jp_list_module.txt")
Module_list = read.table(Module_list, sep = "\t", header = F, stringsAsFactors = F)
Module_id_list = Module_list[,1]

# ==============================================================================
# Main - 1. KEGG module presence across species grouped by Group2 tax list
# ==============================================================================

# Merge all results by species
res_df = Tax_df

for ( id in Module_id_list ) {
  
  res = read.table(sprintf("%s/species_abb/%s_species_abb.txt", dir, id), sep=" ", header = F, stringsAsFactors = F)
  res[[id]] = 1  #presence = 1
  colnames(res) = c("Abbreviation", id)
  
  res_df[[id]] <- res[ match(res_df$Abbreviation, res$Abbreviation), id ]
  
}

res_df[is.na(res_df)] = 0 #NA->0

# Calculate module presence fraction by Group 1,2
tax_count_group1 = table(res_df$Group1)
tax_count_group2 = table(res_df$Group2)

  ## --- a-1. Group1 ---
  res_df_group1 <- rowsum(res_df[ , -c(1:6)], group = res_df$Group1)
  res_df_group1 <- sweep(res_df_group1, 1, tax_count_group1, FUN = "/") #calculate fraction
  res_df_group1 = 100*res_df_group1 #convert to percentage
  
  res_df_group1 = res_df_group1[unique(Tax_df$Group1),] #re-order groups
  
  ## --- a-2. Group2 ---
  res_df_group2 <- rowsum(res_df[ , -c(1:6)], group = res_df$Group2)
  res_df_group2 <- sweep(res_df_group2, 1, tax_count_group2, FUN = "/") #calculate fraction
  res_df_group2 = 100*res_df_group2 #convert to percentage

  res_df_group2 = res_df_group2[unique(Tax_df$Group2),] #re-order groups

  ## --- a-3. Prokaryotes + Eukaryotes (Group 2) ---
  Tax_df_group1 = unique(Tax_df[,c("Supergroup","Group1")])
  Tax_df_group2 = unique(Tax_df[,c("Supergroup","Group1","Group2")])
  Pro_idx_group1 = which(Tax_df_group1$Supergroup == "Prokaryotes")
  Pro_idx_group2 = which(Tax_df_group2$Supergroup == "Prokaryotes")
  
  res_df_final = rbind(res_df_group2[-Pro_idx_group2,], res_df_group1[Pro_idx_group1,])
  res_df_final = res_df_final[,!colSums(res_df_final)<10*nrow(res_df_final)]  #filter-out unsignificant modules
  res_df_final = res_df_final[!rowSums(res_df_final)<10*ncol(res_df_final),]  #filter-out unsignificant tax groups

# Save CSVs
write.csv(res_df, sprintf("%s/module.presence.species.csv", dir_res))
write.csv(res_df_group1, sprintf("%s/module.presence.group1.csv", dir_res))
write.csv(res_df_group2, sprintf("%s/module.presence.group2.csv", dir_res))
write.csv(res_df_final, sprintf("%s/module_presence_final.csv", dir_res))


# Heatmap
group_list_gap = Tax_df_group2[-Pro_idx_group2,]
group_list_gap_idx = c(1,1+which(diff(as.numeric(group_list_gap$Group1))!=0)) - 1
group_list_gap_idx = c(group_list_gap_idx, nrow(res_df_final)-2)

col_red = colorRampPalette(interpolate="l",c("gray90", "rosybrown1","maroon3","deeppink4"))

  ## --- b. Define clustering function ---
  cl_cb <- function(hcl, mat){
    # Recalculate manhattan distances for reorder method
    dists <- dist(mat, method = "manhattan")
    
    # Perform reordering according to OLO method
    hclust_olo <- reorder(hcl, dists)
    return(hclust_olo)
  }

heatmap_p = pheatmap(res_df_final, color = col_red(10), 
           gaps_row = group_list_gap_idx,
           cellwidth = 5, cellheight = 5, fontsize = 5,
           border_color = "white", cluster_cols=T, cluster_rows=F, cutree_cols = 3, display_numbers = F, 
           clustering_distance_cols = "manhattan", clustering_callback = cl_cb,
           main="KEGG metablic module presence (fraction)")

pdf(file=sprintf("%s/module.presence.final.pdf", dir_res),height=5,width=12)
heatmap_p
dev.off()

# Save Tables for cluster information (each three CSVs were later merged into a single Xlsx file: results/clustered_modules.xlsx)
hc_col <- heatmap_p$tree_col #Clustering result
clusters <- cutree(hc_col, k = 3) #Clustering into three groups as done with heatmap

cluster_list = c(1, 2, 3)
for ( clu in cluster_list ) {
  cluster_module_names <- names(clusters)[clusters == clu]
  cluster_module_df = Module_list[Module_list$V1 %in% cluster_module_names, ]
  write.csv(cluster_module_df, sprintf("%s/cluster%s_modules.csv", dir_res, clu), row.names = F, col.names = F)
}


