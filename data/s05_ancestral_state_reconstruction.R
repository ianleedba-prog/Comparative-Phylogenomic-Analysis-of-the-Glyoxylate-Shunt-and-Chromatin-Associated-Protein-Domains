# ==============================================================================
# R 4.4.3 (arm64)
# ==============================================================================
# Load Libraries
library(ape)      # v5.8-1
library(phytools) # v2.4-4
library(expm)     # v1.0-0

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

# 4. Function to order the nodes of a tree
tree_order <- function(tr) {
  n_tip <- length(tr$tip.label)
  root  <- n_tip + 1L
  kids  <- split(tr$edge[, 2], tr$edge[, 1])

  preorder <- integer(0)
  stack    <- root
  while (length(stack)) {
    nd    <- stack[1L]
    stack <- stack[-1L]
    preorder <- c(preorder, nd)
    k <- kids[[as.character(nd)]]
    if (!is.null(k)) stack <- c(k, stack)
  }

  edge_of <- integer(n_tip + tr$Nnode)
  edge_of[tr$edge[, 2]] <- seq_len(nrow(tr$edge))

  list(children = kids, preorder = preorder, postorder = rev(preorder),
       edge_of = edge_of, root = root, n_tip = n_tip, n_node = n_tip + tr$Nnode)
}

# 5. Function to list the tips below every node
descendant_tips <- function(tr) {
  ord <- tree_order(tr)
  out <- vector("list", ord$n_node)
  for (nd in ord$postorder) {
    k <- ord$children[[as.character(nd)]]
    out[[nd]] <- if (is.null(k)) tr$tip.label[nd] else unlist(out[k], use.names = FALSE)
  }
  out
}

# 6. Function to name a branch by the node it subtends
branch_label <- function(tr, nd, n_desc) {
  if (nd <= length(tr$tip.label)) return(tr$tip.label[nd])
  lab <- if (is.null(tr$node.label)) "" else tr$node.label[nd - length(tr$tip.label)]
  if (is.na(lab) || lab == "") sprintf("internal node (%d spp.)", n_desc) else lab
}

# 7. Function to read off the crown ancestor of a clade and the node above it
clade_nodes <- function(tr, tips) {
  crown  <- unname(if (length(tips) == 1) match(tips, tr$tip.label) else getMRCA(tr, tips))
  parent <- unname(tr$edge[tr$edge[, 2] == crown, 1])
  c(crown = crown, parent = parent)
}

# 8. Function to write the root prior that fixes the root to one state
root_pi <- function(state) replace(numeric(4), state, 1)

# 9. Function to build the rate matrix of the dependent four-state model
#    1 = GS-/TET-, 2 = GS-/TET+, 3 = GS+/TET-, 4 = GS+/TET+. One character
#    changes at a time, so the four double changes stay at zero.
make_Q <- function(rates) {
  Q <- matrix(0, 4, 4, dimnames = list(1:4, 1:4))
  Q[1, 2] <- rates[["q12"]]; Q[1, 3] <- rates[["q13"]]
  Q[2, 1] <- rates[["q21"]]; Q[2, 4] <- rates[["q24"]]
  Q[3, 1] <- rates[["q31"]]; Q[3, 4] <- rates[["q34"]]
  Q[4, 2] <- rates[["q42"]]; Q[4, 3] <- rates[["q43"]]
  diag(Q) <- -rowSums(Q)
  Q
}

# 10. Function to run the downward and upward passes of the pruning algorithm
#     L  : likelihood of the subtree below a node, given its state
#     U  : message from the rest of the tree, taken at the lower end of the branch
#          above the node
#     Utop : the same message carried down that branch
tree_passes <- function(tr, states, Q, pi) {
  ord <- tree_order(tr)
  P   <- lapply(tr$edge.length, function(len) expm(Q * len))

  L         <- matrix(0, ord$n_node, 4)
  log_scale <- 0
  for (nd in ord$postorder) {
    k <- ord$children[[as.character(nd)]]
    ## states is the vector fitMk is given, whose entries are the state labels
    if (is.null(k)) { L[nd, as.integer(states[[tr$tip.label[nd]]])] <- 1; next }
    v <- rep(1, 4)
    for (ch in k) v <- v * as.vector(P[[ord$edge_of[ch]]] %*% L[ch, ])
    s_node    <- sum(v)
    log_scale <- log_scale + log(s_node)
    L[nd, ]   <- v / s_node
  }

  U    <- matrix(0, ord$n_node, 4)
  Utop <- matrix(0, ord$n_node, 4)
  Utop[ord$root, ] <- pi
  for (nd in ord$preorder) {
    k <- ord$children[[as.character(nd)]]
    if (is.null(k)) next
    for (ch in k) {
      v <- Utop[nd, ]
      for (sib in k[k != ch]) v <- v * as.vector(P[[ord$edge_of[sib]]] %*% L[sib, ])
      U[ch, ]    <- v
      Utop[ch, ] <- as.vector(v %*% P[[ord$edge_of[ch]]])
    }
  }

  list(L = L, U = U, Utop = Utop, P = P, order = ord,
       logLik = log(sum(pi * L[ord$root, ])) + log_scale)
}

# 11. Function to take the marginal posterior of every node from those passes
node_marginals <- function(pass) {
  M <- pass$Utop * pass$L
  M / rowSums(M)
}

# 12. Function to place the nodes and weights of an n-point Gauss-Legendre rule on [-1, 1], by eigen decomposition of the Jacobi matrix (Golub and Welsch)
gauss_legendre <- function(n) {
  k <- seq_len(n - 1)
  b <- k / sqrt(4 * k^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1)] <- b
  J[cbind(k + 1, k)] <- b
  ev <- eigen(J, symmetric = TRUE)
  o  <- order(ev$values)
  list(x = ev$values[o], w = 2 * ev$vectors[1, o]^2)
}

# 13. Function to count the expected number of each transition on every branch
#     On a branch of length t whose ends are in states a and b, the expected number of i -> j changes is q_ij * integral_0^t P_s[a, i] P_(t-s)[j, b] ds divided by P_t[a, b];
#     the integral is taken by Gauss-Legendre quadrature and the result averaged over the joint posterior of (a, b) on that branch.
#     The counts are therefore analytic rather than averages over stochastic maps, and the table is a deterministic function of the tree, the states and the rates.
expected_changes <- function(tr, Q, pass, pairs, n_quad = 30) {

  gl  <- gauss_legendre(n_quad)
  out <- matrix(0, nrow(tr$edge), length(pairs), dimnames = list(NULL, names(pairs)))

  for (e in seq_len(nrow(tr$edge))) {

    ch    <- tr$edge[e, 2]
    len   <- tr$edge.length[e]
    P_len <- pass$P[[e]]

    ## joint posterior of the states at the two ends of the branch
    J <- (pass$U[ch, ] %o% pass$L[ch, ]) * P_len
    J <- J / sum(J)

    if (len <= 0) next

    s  <- (gl$x + 1) / 2 * len
    w  <- gl$w * len / 2
    P_a <- lapply(s, function(si) expm(Q * si))
    P_b <- lapply(s, function(si) expm(Q * (len - si)))

    for (r in names(pairs)) {
      i <- pairs[[r]][1]; j <- pairs[[r]][2]
      Int <- matrix(0, 4, 4)
      for (k in seq_along(s)) Int <- Int + w[k] * (P_a[[k]][, i] %o% P_b[[k]][j, ])
      E <- Q[i, j] * Int / P_len
      E[!is.finite(E)] <- 0
      out[e, r] <- sum(J * E)
    }
  }
  out
}

# 14. Function to fit the model with the root fixed to one state
fit_fixed_root <- function(tr, states, state, model) {
  fitMk(tr, states, model = model, pi = root_pi(state))
}

# 15. Functions to draw circular tree plot (Fig. 2A)
  ## 15-1. Function to place the tree on a disc
  fan_layout <- function(tr, open_angle = 20, rotate = 100) {
    
    ord <- tree_order(tr)
    
    d <- numeric(ord$n_node)
    for (nd in ord$preorder) {
      k <- ord$children[[as.character(nd)]]
      if (!is.null(k)) d[k] <- d[nd] + 1
    }
    m <- numeric(ord$n_node)
    for (nd in ord$postorder) {
      k <- ord$children[[as.character(nd)]]
      if (!is.null(k)) m[nd] <- max(m[k] + 1)
    }
    r <- d / (d + m)
    stopifnot(all(r[tr$edge[, 2]] > r[tr$edge[, 1]]))
    
    tip_order <- ord$preorder[ord$preorder <= ord$n_tip]
    step      <- 2 * pi * (1 - open_angle / 360) / (ord$n_tip - 1)
    theta     <- numeric(ord$n_node)
    theta[tip_order] <- rotate * pi / 180 + (seq_along(tip_order) - 1) * step
    for (nd in ord$postorder) {
      k <- ord$children[[as.character(nd)]]
      if (!is.null(k)) theta[nd] <- mean(theta[k])
    }
    
    list(r = r, theta = theta, tip_order = tip_order, tip_gap = step,
         open_mid = (rotate - open_angle / 2) * pi / 180)
  }
  
  ## 15-2. Function to draw one branch of that disc
  draw_fan_edge <- function(lay, parent, child, col, lwd) {
    r <- lay$r; th <- lay$theta
    if (r[parent] > 1e-9 && abs(th[child] - th[parent]) > 1e-9) {
      a <- seq(th[parent], th[child], length.out = 30)
      lines(r[parent] * cos(a), r[parent] * sin(a), col = col, lwd = lwd, lend = 1)
    }
    lines(c(r[parent], r[child]) * cos(th[child]),
          c(r[parent], r[child]) * sin(th[child]), col = col, lwd = lwd, lend = 1)
  }
  
  ## 15-3. Function to draw the circular tree
  plot_state_fan <- function(tr, lay, node_state, tip_state, fp, groups, bracket,
                             root_state, state_col, highlight = integer(0)) {
    
    r <- lay$r; th <- lay$theta; gap <- lay$tip_gap
    R_RING <- 1.035; RING_T <- 0.058
    R_ARC  <- R_RING + RING_T + 0.028
    R_BRK  <- R_ARC + 0.052
    R_LAB  <- R_BRK + 0.030
    
    span <- function(tips) range(th[match(tips, tr$tip.label)]) + c(-1, 1) * gap * 0.32
    
    radial_label <- function(radius, angle, txt, ...) {
      flip <- cos(angle) < 0
      text(radius * cos(angle), radius * sin(angle), txt,
           srt = angle * 180 / pi + ifelse(flip, 180, 0),
           adj = c(ifelse(flip, 1, 0), 0.5), ...)
    }
    
    pdf(file = fp, height = 5.6, width = 5.6, pointsize = 7)
    par(mar = c(0.2, 0.2, 0.2, 0.2), xpd = NA)
    plot(NA, xlim = c(-1.55, 1.55), ylim = c(-1.55, 1.55), asp = 1,
         axes = F, xlab = "", ylab = "")
    
    ## branches, then the ones the text names drawn again over a white halo
    for (e in seq_len(nrow(tr$edge)))
      draw_fan_edge(lay, tr$edge[e, 1], tr$edge[e, 2],
                    state_col[node_state[tr$edge[e, 2]]], 0.95)
    for (nd in highlight) {
      pa <- tr$edge[tr$edge[, 2] == nd, 1]
      draw_fan_edge(lay, pa, nd, "white", 3.6)
      draw_fan_edge(lay, pa, nd, state_col[node_state[nd]], 2.2)
    }
    
    ## observed state of every species
    for (k in lay$tip_order) {
      aa <- seq(th[k] - gap * 0.42, th[k] + gap * 0.42, length.out = 6)
      polygon(c(R_RING * cos(aa), rev((R_RING + RING_T) * cos(aa))),
              c(R_RING * sin(aa), rev((R_RING + RING_T) * sin(aa))),
              col = state_col[tip_state[tr$tip.label[k]]], border = NA)
    }
    
    ## clade arcs and their labels
    for (g in names(groups)) {
      s  <- span(groups[[g]])
      aa <- seq(s[1], s[2], length.out = 60)
      lines(R_ARC * cos(aa), R_ARC * sin(aa), col = "black", lwd = 1.3, lend = 1)
      radial_label(R_LAB, mean(s), sprintf("%s  %d", g, length(groups[[g]])), cex = 0.78)
    }
    
    ## the bracket, one tier further out
    s  <- span(bracket$tips)
    aa <- seq(s[1], s[2], length.out = 400)
    lines(R_BRK * cos(aa), R_BRK * sin(aa), col = "#E08A1E", lwd = 1.6, lend = 1)
    for (a_end in s)
      lines(c(R_BRK, R_BRK + 0.03) * cos(a_end), c(R_BRK, R_BRK + 0.03) * sin(a_end),
            col = "#E08A1E", lwd = 1.6, lend = 1)
    ## the two labels sit on opposite sides of the open wedge, the bracket by the
    ## end of its own arc and the root by the other
    radial_label(R_BRK, s[1] - gap * 1.6,
                 sprintf("%s  %d", bracket$name, length(bracket$tips)),
                 cex = 0.95, font = 2)
    
    ## the root, with the state it was fixed to; its leader runs out through the
    ## open wedge, which is the one radius that crosses no branch
    a_root <- lay$open_mid
    lines(c(0.06, R_ARC) * cos(a_root), c(0.06, R_ARC) * sin(a_root),
          col = "grey55", lwd = 0.6, lty = 2)
    points(0, 0, pch = 21, cex = 1.5, lwd = 0.8, bg = state_col[root_state], col = "black")
    radial_label(R_LAB, a_root,
                 sprintf("Holozoa root (%s)", names(state_col)[root_state]), cex = 0.8)
    
    legend(-1.55, 1.55, bty = "n", cex = 0.85, pt.cex = 1.1, y.intersp = 1.15,
           title = "Branch, reconstructed state\nRing, observed state", title.adj = 0,
           legend = sprintf("S%d  %s", seq_along(state_col), names(state_col)),
           pch = 22, pt.bg = state_col, col = state_col)
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

# Per-species domain counts of the 259-species collection (written by s03_gene_counts.R)
rep_sps_count_tab <- read.csv(sprintf("%s/Count_matrix_by_species_259.csv", dir_res),
                              header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# The four states and the eight transitions between them
STATE_LAB <- c("GS-/TET-", "GS-/TET+", "GS+/TET-", "GS+/TET+")

RATE_PAIRS <- list(q12 = c(1, 2), q13 = c(1, 3), q21 = c(2, 1), q24 = c(2, 4),
                   q31 = c(3, 1), q34 = c(3, 4), q42 = c(4, 2), q43 = c(4, 3))

# The index of each free rate in the model matrix passed to fitMk
MODEL <- matrix(c(0, 1, 2, 0,
                  3, 0, 0, 4,
                  5, 0, 0, 6,
                  0, 7, 8, 0), 4, 4, byrow = TRUE, dimnames = list(1:4, 1:4))

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

sps_state$State <- 1L + 2L * sps_state$GS + sps_state$TET # 1 = GS-/TET-, 2 = GS-/TET+, 3 = GS+/TET-, 4 = GS+/TET+
sps_state$Clade <- ifelse(sps_state$Kingdom == "Premetazoa", "Premetazoa", sps_state$Phylum)

states <- setNames(as.character(sps_state$State), sps_state$Species)

  ## summary states
  cat(sprintf("%-10s %3d species\n", STATE_LAB, tabulate(sps_state$State, 4)), sep = "")

# 1-2. The holozoan tree
tree <- taxon_tree(sps_state$Species, read_collection_tree(file.path(dir, "Eukaryota_259.nwk")))

desc  <- descendant_tips(tree)
n_tip <- Ntip(tree)
n_node <- n_tip + tree$Nnode

# 1-3. The nodes the reconstruction is read at
phylum_order <- c("Chordata", "Mollusca", "Echinodermata", "Annelida", "Arthropoda",
                  "Porifera", "Nematoda", "Platyhelminthes", "Tardigrada", "Cnidaria")
phylum_n <- table(sps_state$Phylum[sps_state$Kingdom == "Metazoa"])

clade_tips <- c(list(Metazoa = sps_state$Species[sps_state$Kingdom == "Metazoa"]),
                split(sps_state$Species, sps_state$Phylum)[phylum_order])
clade_nd <- t(sapply(clade_tips, clade_nodes, tr = tree))
mrca_metazoa <- clade_nd["Metazoa", "crown"]

  ## summary metazoan crown node and descendents
  cat(sprintf("Metazoan crown ancestor: node %d, %d descendant species\n",
              mrca_metazoa, length(desc[[mrca_metazoa]])))

# ==============================================================================
# Main - 2. The model fitted under each root in turn
# ==============================================================================

# 2-1. Fit the model four times, with the root fixed to each state
fits <- lapply(1:4, function(k) fit_fixed_root(tree, states, k, MODEL))

rates_tab <- as.data.frame(t(sapply(fits, function(f) f$rates)))
names(rates_tab)    <- names(RATE_PAIRS)
rownames(rates_tab) <- paste0("S", 1:4)
fit_lnL <- sapply(fits, function(f) as.numeric(logLik(f)))

    ### Save the fits (CSV)
    write.csv(cbind(Root = paste0("S", 1:4), logLik = fit_lnL, rates_tab),
              file = sprintf("%s/Model_fits_four_roots_Holozoa.csv", dir_res),
              row.names = FALSE, quote = FALSE)

# 2-2. The likelihoods and the state of the metazoan crown ancestor under each root
passes <- lapply(1:4, function(k)
  tree_passes(tree, states, make_Q(rates_tab[k, ]), root_pi(k)))
margins <- lapply(passes, node_marginals)

lnL   <- sapply(passes, function(p) p$logLik)
d_lnL <- max(lnL) - lnL
crown <- t(sapply(margins, function(M) M[mrca_metazoa, ]))

root_state_tab <- data.frame(
  Root_fixed_to         = STATE_LAB,
  State                 = paste0("S", 1:4),
  logLik                = lnL,
  dlnL_from_best        = d_lnL,
  Relative_likelihood   = exp(-d_lnL),
  P_S1_crown            = crown[, 1],
  P_S2_crown            = crown[, 2],
  P_S3_crown            = crown[, 3],
  P_S4_crown            = crown[, 4],
  P_TETpos_crown        = crown[, 2] + crown[, 4],
  P_GSpos_crown         = crown[, 3] + crown[, 4],
  stringsAsFactors = FALSE)

    ### Save results (CSV)
    write.csv(root_state_tab, file = sprintf("%s/Root_state_fixed_roots_Holozoa.csv", dir_res),
              row.names = FALSE, quote = FALSE)

# ==============================================================================
# Main - 3. The reconstruction on every node
# ==============================================================================

# 3-1.
# Only the two candidate roots are carried forward: S1 (GS-/TET-) and S3 (GS+/TET-)
ROOTS <- c(S1 = 1L, S3 = 3L)

node_state_tab <- data.frame(
  Node     = seq_len(n_node),
  Label    = sapply(seq_len(n_node), function(nd) branch_label(tree, nd, length(desc[[nd]]))),
  Is_tip   = seq_len(n_node) <= n_tip,
  N_descendant_species = sapply(desc, length),
  Observed_state = c(STATE_LAB[sps_state$State[match(tree$tip.label, sps_state$Species)]],
                     rep(NA_character_, tree$Nnode)),
  stringsAsFactors = FALSE)

for (r in names(ROOTS)) {
  M <- margins[[ROOTS[[r]]]]
  node_state_tab[[paste0(r, " state")]] <- STATE_LAB[max.col(M, ties.method = "first")]
  node_state_tab[[paste0(r, " PP")]]    <- apply(M, 1, max)
  node_state_tab[[paste0(r, " P(GS+)")]]  <- M[, 3] + M[, 4]
  node_state_tab[[paste0(r, " P(TET+)")]] <- M[, 2] + M[, 4]
}

    ### Save results (CSV)
    write.csv(node_state_tab, file = sprintf("%s/Node_states_S1_S3_Holozoa.csv", dir_res),
              row.names = FALSE, quote = FALSE)
    
# 3-2. Plot tree
fan <- fan_layout(tree, open_angle = 20, rotate = 100)

fan_groups <- split(sps_state$Species, sps_state$Clade)
fan_groups <- fan_groups[sapply(fan_groups, length) >= 2]

fan_bracket <- list(name = "Metazoa",
                    tips = sps_state$Species[sps_state$Kingdom == "Metazoa"])

STATE_COL <- setNames(c("#BEBEBE", "#54CDF5", "#1771C9", "#093451"), STATE_LAB)

for (r in names(ROOTS))
  plot_state_fan(tree, fan,
                 node_state = max.col(margins[[ROOTS[[r]]]], ties.method = "first"),
                 tip_state  = setNames(sps_state$State, sps_state$Species),
                 fp         = sprintf("%s/Circular_tree_root_%s_Holozoa.pdf", dir_res, r),
                 groups     = fan_groups, bracket = fan_bracket,
                 root_state = ROOTS[[r]], state_col = STATE_COL)

# ==============================================================================
# Main - 4. The crown ancestor of each clade and the node above it
# ==============================================================================

clade_node_tab <- data.frame(Phylum = rownames(clade_nd), stringsAsFactors = FALSE)

for (r in names(ROOTS)) {
  M <- margins[[ROOTS[[r]]]]
  for (nd_type in c("crown", "parent")) {
    nd <- clade_nd[, nd_type]
    clade_node_tab[[sprintf("%s %s state", r, nd_type)]] <-
      STATE_LAB[max.col(M[nd, , drop = FALSE], ties.method = "first")]
    clade_node_tab[[sprintf("%s %s PP", r, nd_type)]] <-
      round(apply(M[nd, , drop = FALSE], 1, max), 4)
  }
}

## the published column order puts the parent before the crown within each root
clade_node_tab <- clade_node_tab[, c("Phylum",
  "S1 crown state", "S1 crown PP", "S1 parent state", "S1 parent PP",
  "S3 crown state", "S3 crown PP", "S3 parent state", "S3 parent PP")]

    ### Save results (CSV)
    write.csv(clade_node_tab, file = sprintf("%s/Clade_nodes_fixedroots_Holozoa.csv", dir_res),
              row.names = FALSE, quote = FALSE)

# ==============================================================================
# Main - 5. The expected changes on every branch
# ==============================================================================

# 5-1. Expected counts of the eight transitions on each of the branches
N_QUAD <- 30
    
ec <- lapply(ROOTS, function(k)
  expected_changes(tree, make_Q(rates_tab[k, ]), passes[[k]], RATE_PAIRS, n_quad = N_QUAD))

branch_nd <- tree$edge[, 2]

branch_tab <- data.frame(
  Branch  = sapply(branch_nd, function(nd) branch_label(tree, nd, length(desc[[nd]]))),
  Descendant_species = sapply(desc[branch_nd], length),
  Clades  = sapply(desc[branch_nd], function(tips)
              paste(sort(unique(sps_state$Clade[match(tips, sps_state$Species)])), collapse = ", ")),
  Branch_length = tree$edge.length,
  stringsAsFactors = FALSE)

  ## one column per event, named as the figure names it rather than by rate
  EVENT_LAB <- c(q13 = "GS gain (TET-)",  q24 = "GS gain (TET+)",
                 q31 = "GS loss (TET-)",  q42 = "GS loss (TET+)",
                 q12 = "TET gain (GS-)",  q34 = "TET gain (GS+)",
                 q21 = "TET loss (GS-)",  q43 = "TET loss (GS+)")
  
  for (r in names(ROOTS)) {
    for (q in names(EVENT_LAB))
      branch_tab[[sprintf("%s %s", r, EVENT_LAB[[q]])]] <- round(ec[[r]][, q], 3)
  }
  for (r in names(ROOTS)) {
    M <- margins[[ROOTS[[r]]]]
    branch_tab[[sprintf("%s P(GS+) at node",  r)]] <- round(M[branch_nd, 3] + M[branch_nd, 4], 4)
    branch_tab[[sprintf("%s P(TET+) at node", r)]] <- round(M[branch_nd, 2] + M[branch_nd, 4], 4)
  }
  
  ## cut-off value for expected count of each branch = EVENT_MIN
  EVENT_MIN <- 0.15 
  
  largest <- apply(cbind(ec$S1, ec$S3), 1, max)
  branch_tab <- branch_tab[largest >= EVENT_MIN, ]
  branch_tab <- branch_tab[order(-branch_tab$Descendant_species, branch_tab$Branch), ]
  
      ### Save results (CSV)
      write.csv(branch_tab, file = sprintf("%s/Branch_expected_changes_Holozoa.csv", dir_res),
                row.names = FALSE, quote = which(names(branch_tab) == "Clades"))

# 5-2. Gains of Tet_JBP across Metazoa
gain <- sapply(ec, function(m) m[, "q12"] + m[, "q34"])
rownames(gain) <- branch_nd

i_stem <- which(branch_nd == mrca_metazoa)
i_psam <- which(branch_nd == match("Plectus sambesii", tree$tip.label))
i_rest <- setdiff(seq_len(nrow(gain)), c(i_stem, i_psam))

gain_tab <- rbind(
  "Total expected gains"            = colSums(gain),
  "Metazoan stem"                   = gain[i_stem, ],
  "Plectus sambesii (terminal)"     = gain[i_psam, ],
  "All other branches"              = colSums(gain[i_rest, ]),
  "Share, metazoan stem"            = gain[i_stem, ] / colSums(gain),
  "Share, P. sambesii"              = gain[i_psam, ] / colSums(gain),
  "Share, two branches combined"    = (gain[i_stem, ] + gain[i_psam, ]) / colSums(gain),
  "Share, all other branches"       = colSums(gain[i_rest, ]) / colSums(gain),
  "Largest single other branch"     = apply(gain[i_rest, ], 2, max))

gain_tab <- data.frame(Quantity = rownames(gain_tab),
                       `Root S1 (GS-/TET-)` = gain_tab[, "S1"],
                       `Root S3 (GS+/TET-)` = gain_tab[, "S3"],
                       check.names = FALSE, stringsAsFactors = FALSE)

  ## summary of Tet_JBP gains
  print(format(gain_tab, digits = 4), row.names = FALSE)

    ### Save results (CSV)
    write.csv(gain_tab, file = sprintf("%s/Expected_gains_Tet_JBP_Metazoa.csv", dir_res),
              row.names = FALSE, quote = which(names(gain_tab) == "Quantity"))
