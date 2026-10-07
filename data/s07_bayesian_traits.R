# ==============================================================================
# R 4.4.3 (arm64)
# ==============================================================================
# Load Libraries
library(ape)      # v5.8-1
library(phytools) # v2.4-4

# ==============================================================================
# Define functions
# ==============================================================================

# 1. Function to bring tip labels and species names to one form (no quotes, spaces)
norm_name <- function(v) gsub("_", " ", gsub("['\"]", "", v))

# 2. Function to read a newick tree
read_collection_tree <- function(fp, drop_node_labels = FALSE) {
  tr <- read.tree(fp)
  if (drop_node_labels) tr$node.label <- NULL
  tr$tip.label <- norm_name(tr$tip.label)
  tr
}

# 3. Function to cut one taxon set out of the collection tree
taxon_tree <- function(species, tree) keep.tip(tree, intersect(tree$tip.label, species))

# 4. Function to write the root prior that fixes the root to one state
root_pi <- function(state) replace(numeric(4), state, 1)

# 5. Function to read the posterior sample of a BayesTraits run
read_bayestraits_log <- function(fp) {
  txt <- readLines(fp)
  i   <- grep("^Iteration\t", txt)[1]
  if (is.na(i)) stop("no sample table in ", fp)
  tab <- read.delim(text = paste(txt[i:length(txt)], collapse = "\n"),
                    check.names = FALSE, stringsAsFactors = FALSE)
  tab[, trimws(names(tab)) != "", drop = FALSE]
}

# 6. Function to read the marginal likelihood a stepping-stone run ended on
read_stones <- function(fp) {
  line <- grep("^Log marginal likelihood", readLines(fp), value = TRUE)[1]
  if (is.na(line)) stop("no marginal likelihood in ", fp)
  as.numeric(sub(".*:[[:space:]]*", "", line))
}

# 7. Function to summarise a posterior sample by its median and central 95%
posterior_summary <- function(v) {
  q <- unname(quantile(v, c(0.5, 0.025, 0.975)))
  c(Median = q[1], CrI_lower = q[2], CrI_upper = q[3])
}

# 8. Function to take one rate out of a fit by the cell of the model it sits in
rate_of <- function(fit, model, from, to) unname(fit$rates[model[from, to]])

# 9. Function to fit the model with the root fixed to one state
fit_fixed_root <- function(tr, states, state, model) {
  fitMk(tr, states, model = model, pi = root_pi(state), quiet = TRUE)
}

# 10. Function to refit the model on the tree with one clade removed (clade decomposition)
leave_clade_out <- function(tr, states, state, drop, model, model_eq, rate_zero) {

  keep <- setdiff(tr$tip.label, drop)
  sub  <- if (length(drop)) keep.tip(tr, keep) else tr
  st   <- states[sub$tip.label]

  free <- fit_fixed_root(sub, st, state, model)
  eq   <- fit_fixed_root(sub, st, state, model_eq)

  q13 <- rate_of(free, model, 1, 3)
  q24 <- rate_of(free, model, 2, 4)
  lr  <- 2 * (as.numeric(logLik(free)) - as.numeric(logLik(eq)))

  ## a rate driven onto the boundary leaves the ratio undefined, and the figure
  ## reports it as such rather than as a very large number
  ratio <- if (q24 > rate_zero) q13 / q24 else NA_real_

  data.frame(
    Species_removed   = length(drop),
    Species_remaining = Ntip(sub),
    q13 = q13, q24 = q24,
    ratio = ratio, log2_ratio = log2(ratio),
    logLik_free_rates = as.numeric(logLik(free)),
    logLik_q13_eq_q24 = as.numeric(logLik(eq)),
    LR = lr, P = pchisq(lr, df = 1, lower.tail = FALSE),
    stringsAsFactors = FALSE)
}

# ==============================================================================
# Environment & Data Setup
# ==============================================================================

# Setup directories
dir <- "/directory/to/data"
dir_analysis <- file.path(dir, "analysis")
dir_res <- file.path(dir_analysis, "results")
if(!dir.exists(dir_res)) dir.create(dir_res, recursive = TRUE, showWarnings = FALSE)

# BayesTraits runs, laid out by s06_Bayesian.sh
dir_bayes     <- file.path(dir, "Bayesian_traits")
dir_bayes_res <- file.path(dir_bayes, "results")

# Per-species domain counts of the 259-species collection (written by s03_gene_counts.R)
rep_sps_count_tab <- read.csv(sprintf("%s/Count_matrix_by_species_259.csv", dir_res),
                              header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# The four states and the eight transitions between them
STATE_LAB <- c("GS-/TET-", "GS-/TET+", "GS+/TET-", "GS+/TET+")

# The index of each free rate in the model matrix passed to fitMk.
MODEL <- matrix(c(0, 1, 2, 0,
                  3, 0, 0, 4,
                  5, 0, 0, 6,
                  0, 7, 8, 0), 4, 4, byrow = TRUE, dimnames = list(1:4, 1:4))

MODEL_EQ <- matrix(c(0, 1, 2, 0,
                     3, 0, 0, 2,
                     4, 0, 0, 5,
                     0, 6, 7, 0), 4, 4, byrow = TRUE, dimnames = list(1:4, 1:4))

BAYES_MODEL <- c("Full dependent (DP)"      = "DP_full",
                 "Dependent, q13 = q24"     = "DP_res.q13.q24",
                 "Independent (IDP)"        = "IDP_full")
BAYES_TAXON <- c("Metazoa", "Holozoa")

bayes_file <- function(prm, tax, ext)
  sprintf("%s/%s.binary.table.GS_TET.species_%s.%s", dir_bayes_res, prm, tax, ext)

# ==============================================================================
# Main - 1. The four-state character and the holozoan tree
# ==============================================================================

# 1-1. State of every species
cm <- rep_sps_count_tab
holozoa <- cm$Kingdom %in% c("Metazoa", "Premetazoa")

sps_state <- data.frame(
  Species = norm_name(cm$Species[holozoa]),
  Kingdom = cm$Kingdom[holozoa],
  Phylum  = cm$Phylum[holozoa],
  GS      = as.integer(cm$ICL[holozoa] >= 1 & cm$Malate_synthase[holozoa] >= 1),
  TET     = as.integer(cm$Tet_JBP[holozoa] >= 1),
  stringsAsFactors = FALSE)

sps_state$State <- 1L + 2L * sps_state$GS + sps_state$TET
sps_state$Clade <- ifelse(sps_state$Kingdom == "Premetazoa", "Premetazoa", sps_state$Phylum)

states <- setNames(as.character(sps_state$State), sps_state$Species)

# 1-2. The holozoan tree
tree <- taxon_tree(sps_state$Species, read_collection_tree(file.path(dir, "Eukaryota_259.nwk")))

# ==============================================================================
# Main - 2. Posterior rates of gain of the shunt
# ==============================================================================

# 2-1. The posterior sample of each run
post <- do.call(rbind, lapply(BAYES_TAXON, function(tax) {
  s <- read_bayestraits_log(bayes_file("DP_full", tax, "Log.txt"))
  data.frame(Taxon = tax, Iteration = s$Iteration, lnL = s$Lh,
             q13 = s$q13, q24 = s$q24, ratio = s$q13 / s$q24,
             stringsAsFactors = FALSE)
}))

cat(sprintf("%-8s %5d posterior samples, iterations %d to %d\n",
            BAYES_TAXON,
            tapply(post$Iteration, post$Taxon, length)[BAYES_TAXON],
            tapply(post$Iteration, post$Taxon, min)[BAYES_TAXON],
            tapply(post$Iteration, post$Taxon, max)[BAYES_TAXON]), sep = "")

    ### Save the sample (CSV)
    write.csv(post, file = sprintf("%s/Posterior_rates_GS_gain.csv", dir_res),
              row.names = FALSE, quote = FALSE)

# 2-2. Median and central 95% of each quantity
post_summary <- do.call(rbind, lapply(BAYES_TAXON, function(tax) {
  p <- post[post$Taxon == tax, ]
  s <- rbind(posterior_summary(p$q13),
             posterior_summary(p$q24),
             posterior_summary(p$ratio))
  data.frame(Taxon = tax, N_samples = nrow(p),
             Quantity = c("q13 (GS gain, TET- background)",
                          "q24 (GS gain, TET+ background)",
                          "q13 / q24"),
             s, P_q13_gt_q24 = c(mean(p$q13 > p$q24), NA, NA),
             row.names = NULL, stringsAsFactors = FALSE)
}))

print(format(post_summary, digits = 4), row.names = FALSE)

    ### Save results (CSV)
    write.csv(post_summary, file = sprintf("%s/Posterior_summary_GS_gain.csv", dir_res),
              row.names = FALSE, quote = which(names(post_summary) == "Quantity"))

# ==============================================================================
# Main - 3. Marginal likelihoods and Bayes factors
# ==============================================================================

# 3-1. The marginal likelihood each stepping-stone run ended on
marginal <- do.call(rbind, lapply(BAYES_TAXON, function(tax)
  data.frame(Taxon_set = tax,
             Species   = if (tax == "Metazoa") sum(sps_state$Kingdom == "Metazoa") else Ntip(tree),
             Model     = names(BAYES_MODEL),
             log_marginal_likelihood =
               sapply(BAYES_MODEL, function(prm) read_stones(bayes_file(prm, tax, "Stones.txt"))),
             row.names = NULL, stringsAsFactors = FALSE)))

print(format(marginal, digits = 6), row.names = FALSE)

    ### Save results (CSV)
    write.csv(marginal, file = sprintf("%s/BayesTraits_marginal_likelihoods.csv", dir_res),
              row.names = FALSE, quote = which(names(marginal) == "Model"))

# 3-2. Each restriction against the full dependent model
bayes_factor <- do.call(rbind, lapply(BAYES_TAXON, function(tax) {
  m    <- marginal[marginal$Taxon_set == tax, ]
  full <- m$log_marginal_likelihood[m$Model == "Full dependent (DP)"]
  rest <- m[m$Model != "Full dependent (DP)", ]
  data.frame(Taxon_set = tax, Species = rest$Species,
             Comparison = paste(rest$Model, "vs full DP"),
             log_Bayes_factor = 2 * (full - rest$log_marginal_likelihood),
             row.names = NULL, stringsAsFactors = FALSE)
}))

print(format(bayes_factor, digits = 4), row.names = FALSE)

    ### Save results (CSV)
    write.csv(bayes_factor, file = sprintf("%s/BayesTraits_bayes_factors.csv", dir_res),
              row.names = FALSE, quote = which(names(bayes_factor) == "Comparison"))

# ==============================================================================
# Main - 4. Clade decomposition of the two rates of gain (q13, q24)
# ==============================================================================

# 4-1. The clades removed
CLADE_MIN <- 4
RATE_ZERO <- 1e-8

clade_n <- table(sps_state$Clade[sps_state$Kingdom == "Metazoa"])
CLADES  <- names(clade_n)[clade_n >= CLADE_MIN]

# 4-2. Refit under each root with each clade left out (This process takes time)
ROOTS <- c(S1 = 1L, S3 = 3L)

decomposition <- do.call(rbind, lapply(names(ROOTS), function(r) {

  rows <- do.call(rbind, lapply(c("(none: full tree)", CLADES), function(cl) {
    drop <- if (cl == "(none: full tree)") character(0) else
            sps_state$Species[sps_state$Clade == cl]
    cbind(Root = r, Clade_removed = cl,
          leave_clade_out(tree, states, ROOTS[[r]], drop, MODEL, MODEL_EQ, RATE_ZERO))
  }))

  ## every change is taken against the full-tree fit under the same root
  rows$change_in_log2_ratio <- rows$log2_ratio - rows$log2_ratio[rows$Clade_removed == "(none: full tree)"]

  full <- rows[rows$Clade_removed == "(none: full tree)", ]
  rest <- rows[rows$Clade_removed != "(none: full tree)", ]
  rbind(full, rest[order(!is.na(rest$change_in_log2_ratio),
                         -xtfrm(rest$change_in_log2_ratio)), ])
}))

decomposition <- decomposition[, c("Root", "Clade_removed", "Species_removed",
  "Species_remaining", "q13", "q24", "ratio", "log2_ratio", "change_in_log2_ratio",
  "logLik_free_rates", "logLik_q13_eq_q24", "LR", "P")]

print(format(decomposition, digits = 4), row.names = FALSE)

    ### Save results (CSV)
    write.csv(decomposition, file = sprintf("%s/Clade_decomposition.csv", dir_res),
              row.names = FALSE, quote = FALSE)

