# Install needed packages
install.packages(pkgs = c("ape", "ade4", "distory", "gplots", "ggplot2", "phangorn", "phytools"), repos = "https://mirrors.nic.cz/R/", dependencies = "Imports")
# Install kdetrees package (removed from CRAN)
# Ensure package 'devtools' is installed
if (!"devtools" %in% installed.packages()) {
  install.packages("devtools")
}
# Install 'kdetrees' from https://github.com/V-Z/kdetrees Git repository
devtools::install_github("V-Z/kdetrees")

# Load libraries
library(ape)
library(ade4)
library(distory)
library(gplots)
library(ggplot2)
library(kdetrees)
library(phangorn)
library(phytools)

# Set working directory
setwd("~/dokumenty/vyuka/hybseq/")

# Load the list of trees
trees <- read.tree(file = "trees_ml_exons.nwk")
trees
print(trees, details = TRUE)

# Compute distance of topological similarities
trees.d <- dist.topo(x = trees, method = "score", mc.cores = 4) # Set number of cores according to your computer

# Plot the heatmap (package gplots)
png(filename = "trees_dist.png", width = 10000, height = 10000)
heatmap.2(x = as.matrix(trees.d), Rowv = FALSE, Colv = "Rowv", dendrogram = "none", symm = TRUE, scale = "none", na.rm = TRUE, revC = FALSE, col = rainbow(15), cellnote = as.matrix(trees.d), notecex = 1, notecol = "white", trace = "none", labRow = rownames(as.matrix(trees.d)), labCol = colnames(as.matrix(trees.d)), key = FALSE, main = "Correlation matrix of topographical distances")
dev.off() # Saves the image

# Test if the distance matrix is Euclidean
is.euclid(distmat = as.dist(trees.d), plot = TRUE, tol = 1e-05)

# PCoA
trees.pcoa <- dudi.pco(d = trees.d, scannf = FALSE, nf = 5)
trees.pcoa

# Plot PCoA
s.label(dfxy = trees.pcoa$li)
s.kde2d(dfxy = trees.pcoa$li, cpoint = 0, add.plot = TRUE)
add.scatter.eig(trees.pcoa[["eig"]], 3, 1, 2, posi = "topright")
title("PCoA of matrix of pairwise trees distances")

# Remove outlying trees
trees
trees[c("Assembly_1033", "Assembly_10103", "Assembly_10222")] <- NULL
trees

# Now you can repeat recalculation of distance matrix and PCoA and possibly remove more trees...

# # Possibly remove trees with too few tips
# print(trees, details = TRUE)
# trees[c(1, 2, 3, 4)] <- NULL
# trees

# # Possibly remove rare tips
# trees <- lapply(X = trees, FUN = drop.tip, tip = c("Amomum-sp7_S308_L001", "Amomum-trilobum_S12_L001"))
# class(trees) <- "multiPhylo" # Use after usage of lapply to multiPhylo

# Run kdetrees to detect outliers - play with k
?kdetrees # See options for kdetrees
trees.kde <- kdetrees(trees = trees, k = 0.9, distance = "dissimilarity", topo.only = FALSE, greedy = TRUE)
# See text results with list of outlying trees
trees.kde
# See graphical results
plot(x = trees.kde)
hist(x = trees.kde)
# See removed trees
plot.multiPhylo(trees.kde[["outliers"]])
# Save removed trees
write.tree(phy = trees.kde[["outliers"]], file = "trees_outliers.nwk")
# Save kdetrees report
write.table(x = as.data.frame(x = trees.kde), file = "trees_scores.tsv", quote = FALSE, sep = "\t")
# Extract passing trees
trees.good <- trees[names(trees) %in% names(trees.kde[["outliers"]]) == FALSE]
trees.good
# Save passing trees
write.tree(phy = trees.good, file = "trees_good.nwk")

# Compute parsimony super tree
?superTree # See help first...
tree.sp <- superTree(tree = trees.good, method = "NNI", rooted = TRUE, trace = 2, start = NULL, multicore = TRUE)
tree.sp # See details
# Root it
tree.sp <- root(phy = tree.sp, outgroup = c("Riedelia-arfakensis_S49_L001", "Zingiber-officinale_S242_L001"), resolve.root = TRUE)
# Save parsimony super tree
write.tree(phy = tree.sp, file = "parsimony_sp_tree.nwk")
# Plot parsimony super tree
plot.phylo(x = tree.sp, type = "phylogram", edge.width = 2, label.offset = 0.01, cex = 1.2)
add.scale.bar()
# Tune display of the tree...

# See help...
?ape::speciesTree
?phytools::mrp.supertree
?phangorn::coalSpeciesTree
# All trees must be ultrametric - chronos scale them
trees.ultra <- lapply(X = trees.good, FUN = chronos, model = "correlated")
class(trees.ultra) <- "multiPhylo"
# Calculate the species tree
# tree.sp.mean <- speciesTree(x=trees.ultra, FUN=mean)
tree.sp2 <- mrp.supertree(tree = trees.good, method = "optim.parsimony", rooted = TRUE)
tree.sp2 <- root(phy = tree.sp2, outgroup = c("Riedelia-arfakensis_S49_L001", "Zingiber-officinale_S242_L001"), resolve.root = TRUE)
plot.phylo(x = tree.sp2, type = "phylogram", edge.width = 2, label.offset = 0.01, cex = 1.2)

# # Consensus networks
# # Requires all trees to have same set of tips (no missing data)
# # ?consensusNet
# # Compute consensus network
# # tree.net <- consensusNet(obj = trees.good, prob = 0.25)
# # Plot 2D or 3D
# # plot(x = tree.net, planar = FALSE, type = "2D", use.edge.length = TRUE, show.tip.label = TRUE, show.edge.label = TRUE, show.node.label = TRUE, show.nodes = TRUE, edge.color = "black", tip.color = "blue") # 2D
# # plot(x = tree.net, planar = FALSE, type = "3D", use.edge.length = TRUE, show.tip.label = TRUE, show.edge.label = TRUE, show.node.label = TRUE, show.nodes = TRUE, edge.color = "black", tip.color = "blue") # 3D

# Save trees.good in NEXUS for PhyloNet
write.nexus(trees.good, file = "trees_good.nex", translate = FALSE)

# Cophyloplots - comparing 2 phylogenetic trees
# We need 2 column matrix with tip labels
tips.labels <- matrix(data = c(sort(tree.sp[["tip.label"]]), sort(tree.sp2[["tip.label"]])), nrow = length(tree.sp[["tip.label"]]), ncol = 2)
# Draw the tree, play with graphical parameters
# Click to nodes to rotate them to get better display
cophyloplot(x = tree.sp, y = tree.sp2, assoc = tips.labels, use.edge.length = FALSE, space = 60, length.line = 1, gap = 2, type = "phylogram", rotate = TRUE, col = "red", lwd = 1.5, lty = 2)
# Slihtly better display in phytools::cophylo
trees.cophylo <- cophylo(tr1 = tree.sp, tr2 = tree.sp2, assoc = tips.labels, rotate = TRUE)
plot.cophylo(x = trees.cophylo, lwd = 2, link.type = "curved")

# Density trees
is.rooted.multiPhylo(trees.ultra) # rooted
is.ultrametric.multiPhylo(trees.ultra) # ultrametric
is.binary.multiPhylo(trees.ultra) # binary bifurcating
# See help page
?phangorn::densiTree
# Plotting density trees
densiTree(x = trees.ultra, scaleX = TRUE, col = rainbow(6), width = 5, cex = 1.5)
densiTree(x = trees.ultra, direction = "upwards", scaleX = TRUE, width = 5)
densiTree(x = trees.ultra, scaleX = TRUE, width = 5, cex = 1.5)
densiTree(x = trees.ultra[1:10], scaleX = TRUE, width = 5, cex = 1.25)
# See help page
?phytools::densityTree
# Plotting density trees
densityTree(trees = c(tree.sp, tree.sp2), fix.depth = TRUE, lwd = 4)
densityTree(trees = trees.ultra, fix.depth = TRUE, use.gradient = TRUE, alpha = 0.5, lwd = 4)
# densityTree(trees = trees.ultra[1:3], fix.depth = TRUE, use.gradient = TRUE, alpha = 0.5, lwd = 4)
# densityTree(trees = trees.ultra[c(2, 4, 6)], fix.depth = TRUE, use.gradient = TRUE, alpha = 0.5, lwd = 4)
