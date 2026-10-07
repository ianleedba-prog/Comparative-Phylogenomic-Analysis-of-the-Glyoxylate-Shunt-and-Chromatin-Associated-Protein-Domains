# ==============================================================================
# R 4.4.3 (arm64)
# ==============================================================================
# Load Libraries
library(ape) #v5.8-1
library(caper) # v1.0.4

# ==============================================================================
# Define functions
# ==============================================================================

# 1. Function to bring tip labels and species names to one form (no quotes, spaces)
norm_name <- function(v) gsub("_", " ", gsub("['\"]", "", v))

# 2. Function to read a newick tree and bring its labels to the common form
read_collection_tree <- function(fp) {
  tr <- read.tree(fp)
  tr$node.label <- NULL
  tr$tip.label  <- norm_name(tr$tip.label)
  tr
}

# 3. Function to make a tree fully dichotomous
resolve_polytomies <- function(tr, min_len = 1e-6) {
  tr <- multi2di(tr)
  tr$edge.length[tr$edge.length == 0] <- min_len
  tr
}

# 4. Function to cut one taxon set out of the collection tree
taxon_tree <- function(species, tree) keep.tip(tree, intersect(tree$tip.label, species))

# 5. Function to place a leave-one-out jackknife confidence interval on D
#    (refit on every n-1 subset, then mean +/- 1.96 * jackknife standard error)
jackknife_D <- function(tax_bin_dat, tax_tree, varname, permut = 100) {

  species <- tax_bin_dat$Species
  n <- length(species)
  D_jack <- numeric(n)

  for (i in seq_len(n)) {
    dat_sub  <- tax_bin_dat[tax_bin_dat$Species %in% species[-i], ]
    tree_sub <- drop.tip(tax_tree, species[i])

    tryCatch({
      cd_sub <- comparative.data(tree_sub, dat_sub, names.col = "Species")
      res <- eval(bquote(phylo.d(cd_sub, binvar = .(as.name(varname)), permut = permut)))
      D_jack[i] <- res$DEstimate
    }, error = function(e) {
      D_jack[i] <<- NA
    })
  }

  D_jack  <- D_jack[!is.na(D_jack)]
  valid_n <- length(D_jack)

  if (valid_n < 2) return(list(D_CI_lo = NA, D_CI_hi = NA, n_jack = valid_n))

  D_mean <- mean(D_jack)
  D_SE   <- sqrt(((valid_n - 1) / valid_n) * sum((D_jack - D_mean)^2))

  list(D_CI_lo = D_mean - 1.96 * D_SE,
       D_CI_hi = D_mean + 1.96 * D_SE,
       n_jack  = valid_n)
}

# 6. Function to run phylo.d for every character in every taxon set
dstat_table <- function(bin_dat, trees, vars, permut = 10000, permut_jack = 100) {

  res_all <- data.frame()

  for (taxon in names(trees)) {

    tax_tree <- trees[[taxon]]

    tax_bin_dat <- bin_dat[rownames(bin_dat) %in% tax_tree$tip.label, , drop = F]
    tax_bin_dat$Species <- rownames(tax_bin_dat)
    tax_bin_dat <- tax_bin_dat[, c("Species", vars)]

    ## every tip must carry a character value, or comparative.data drops it silently
    stopifnot(setequal(tax_bin_dat$Species, tax_tree$tip.label))

    comp_data <- comparative.data(tax_tree, tax_bin_dat, names.col = "Species")
    cat(sprintf("\n[%s] n = %d tips\n", taxon, Ntip(tax_tree)))

    for (varname in vars) {

      v  <- tax_bin_dat[[varname]]
      n0 <- sum(v == 0); n1 <- sum(v == 1)

      if (n0 == 0 || n1 == 0) {
        res_all <- rbind(res_all, data.frame(
          Taxon_Set = taxon, Variable = varname, State_0 = n0, State_1 = n1,
          Estimated_D = NA, D_CI_lo = NA, D_CI_hi = NA,
          D_random_lo = NA, D_random_hi = NA, D_brownian_lo = NA, D_brownian_hi = NA,
          Pval_Random = NA, Pval_Brownian = NA, N_species = Ntip(tax_tree),
          Note = sprintf("invariant (all %d): D not defined", ifelse(n0 == 0, 1, 0))))
        cat(sprintf("  %-10s invariant (%d / %d) - skipped\n", varname, n0, n1))
        next
      }

      res <- eval(bquote(phylo.d(comp_data, binvar = .(as.name(varname)), permut = permut)))

      ## the two null distributions rescaled onto the D axis, on which
      ## the Brownian mean is 0 and the random mean is 1
      mean_rand  <- mean(res$Permutations$random,   na.rm = T)
      mean_brown <- mean(res$Permutations$brownian, na.rm = T)
      denom      <- mean_rand - mean_brown

      D_rand_ci  <- unname(quantile((res$Permutations$random   - mean_brown) / denom,
                                    c(0.025, 0.975), na.rm = T))
      D_brown_ci <- unname(quantile((res$Permutations$brownian - mean_brown) / denom,
                                    c(0.025, 0.975), na.rm = T))

      jack <- jackknife_D(tax_bin_dat, tax_tree, varname, permut = permut_jack)

      res_all <- rbind(res_all, data.frame(
        Taxon_Set = taxon, Variable = varname,
        State_0 = unname(res$States["0"]), State_1 = unname(res$States["1"]),
        Estimated_D   = res$DEstimate,
        D_CI_lo       = jack$D_CI_lo,       D_CI_hi       = jack$D_CI_hi,
        D_random_lo   = D_rand_ci[1],       D_random_hi   = D_rand_ci[2],
        D_brownian_lo = D_brown_ci[1],      D_brownian_hi = D_brown_ci[2],
        Pval_Random   = res$Pval1,          # against random association (D = 1)
        Pval_Brownian = 1 - res$Pval0,      # against the Brownian expectation (D = 0)
        N_species     = Ntip(tax_tree),
        Note = if (jack$n_jack < Ntip(tax_tree))
                 sprintf("jackknife on %d of %d refits", jack$n_jack, Ntip(tax_tree)) else ""))

      cat(sprintf("  %-10s D = %7.4f  [%7.4f, %7.4f]  P_rand = %.4f  P_brown = %.4f\n",
                  varname, res$DEstimate, jack$D_CI_lo, jack$D_CI_hi,
                  res$Pval1, 1 - res$Pval0))
    }
  }
  res_all
}

# ==============================================================================
# Environment & Data Setup
# ==============================================================================

# Setup directories
dir <- "/directory/to/data"
dir_analysis <- file.path(dir, "analysis")
dir_res <- file.path(dir_analysis, "results")
if(!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE, showWarnings = FALSE)

# Per-species domain counts of the 259-species collection (written by s03_gene_counts.R)
rep_sps_count_tab <- read.csv(sprintf("%s/Count_matrix_by_species_259.csv", dir_res),
                              header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# ==============================================================================
# Main - 1. Fritz and Purvis' D statistic
# ==============================================================================

# Eight binary characters scored on the 259-species collection
#   TCA      : all five components of the cycle present
#   ICL, MLS : the two shunt enzymes on their own
#   GS       : the shunt, ICL together with Malate_synthase
#   TET      : Tet_JBP
#   the three reading domains that co-distribute with Tet_JBP
cm <- rep_sps_count_tab
bin_dat <- data.frame(
  TCA       = as.integer(cm$Citrate_synt >= 1 & cm$Aconitase >= 1 & cm$Iso_dh >= 1 &
                         cm$`2-oxogl_dehyd` >= 1 & cm$FumaraseC_C >= 1),
  ICL       = as.integer(cm$ICL >= 1),
  MLS       = as.integer(cm$Malate_synthase >= 1),
  GS        = as.integer(cm$ICL >= 1 & cm$Malate_synthase >= 1),
  TET       = as.integer(cm$Tet_JBP >= 1),
  ADD_ATRX  = as.integer(cm$ADD_ATRX >= 1),
  ADD_DNMT3 = as.integer(cm$ADD_DNMT3 >= 1),
  TTD       = as.integer(cm$TTD >= 1),
  row.names = norm_name(cm$Species)
)

VARS        <- c("TCA", "ICL", "MLS", "GS", "TET", "ADD_ATRX", "ADD_DNMT3", "TTD")
N_PERMUT    <- 10000 # permutations of the two null distributions
N_PERMUT_JK <- 100   # permutations inside each jackknife refit

    # Save binary table (CSV)
    write.csv(bin_dat, file=sprintf("%s/Binary_table_eight_characters_259.csv", dir_res), row.names = T, quote = F)

# 1-1. The dichotomous collection tree
fp_tree <- file.path(dir, "Eukaryota_259.nwk")
fp_di   <- file.path(dir, "Eukaryota_259_dichotomous.nwk")

if (!file.exists(fp_di)) {
  set.seed(20260914)
  write.tree(resolve_polytomies(read_collection_tree(fp_tree)), file = fp_di)
}

dic_tree <- read_collection_tree(fp_di)

# 1-2. Taxon sets
set_order <- c("Eukaryota", "Opisthokonta", "Holomycota", "Holozoa", "Metazoa")

clade <- ifelse(cm$Kingdom %in% c("Metazoa", "Premetazoa"), "Holozoa",
         ifelse(cm$Kingdom %in% c("Fungi", "Rotosphaerida"), "Holomycota", NA))
sp    <- norm_name(cm$Species)

# return one NA per non-opisthokont species alongside the names wanted
taxon_species <- list(
  Eukaryota    = sp,
  Opisthokonta = sp[which(clade %in% c("Holozoa", "Holomycota"))],
  Holomycota   = sp[which(clade == "Holomycota")],
  Holozoa      = sp[which(clade == "Holozoa")],
  Metazoa      = sp[which(cm$Kingdom == "Metazoa")]
)[set_order]

trees <- lapply(taxon_species, taxon_tree, tree = dic_tree)
cat(sprintf("%-13s %3d species\n", names(trees), sapply(trees, Ntip)), sep = "")

# 1-3. D statistic for every character in every taxon set
set.seed(20260914)
dST_res <- dstat_table(bin_dat, trees, VARS, permut = N_PERMUT, permut_jack = N_PERMUT_JK)

    # Save results (CSV)
    write.csv(dST_res, file=sprintf("%s/Dstatistic_eight_characters_259.csv", dir_res), row.names = F, quote = F)
