if (FALSE) {
  # environment setup:
  devtools::install_github("https://github.com/Nanostring-Biostats/insitutype")
  install("InSituCor")
}
library(tictoc)
library(igraph)
library(pheatmap)
library(scales)
library(Matrix)
library(viridis)
library(InSituType)
library(SeuratObject)
library(SeuratDisk)
library(InSituCor)


#### load and parse data ---------------------------------

# first, download the data from figshare:
# https://figshare.com/articles/dataset/CosMx_6000-plex_colon_cancer_dataset/24164388


# load the data:
load("colon cancer dataset.RData")

# normalize:
norm = sweep(raw, 1, rowSums(raw), "/") * mean(rowSums(raw))

#### load cellKlatch results ----------------------------------

step1 = readRDS("processed_data/lrstep1.RDS")
modules = readRDS("processed_data/modules.RDS")
scores = readRDS("processed_data/scores.RDS")
attribution = readRDS("processed_data/attribution.RDS")
rawcor = readRDS("processed_data/rawcor.RDS")

#### explore ligands ------------------------------

load("data/CellChatDB.human.rda") # (from https://github.com/sqjin/CellChat/blob/master/data/CellChatDB.human.rda)
ligands = intersect(unique(CellChatDB.human$interaction$ligand), colnames(raw))
receptors = intersect(unique(CellChatDB.human$interaction$receptor), colnames(raw))
lrpairs = CellChatDB.human$interaction[, c("ligand", "receptor")]
rownames(lrpairs) = CellChatDB.human$interaction$interaction_name
lrpairs = lrpairs[is.element(lrpairs[,1], colnames(raw)) & is.element(lrpairs[,2], colnames(raw)), ]
lrpairs = lrpairs[lrpairs[,1] != lrpairs[,2], ]

## get ligand modules -------------------------

# get condcor
ligstep1 <- calcSpatialCor(counts = norm[, ligands], 
                           conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust),
                           neighbors = NULL,
                           xy = xy,
                           k = 50,
                           radius = NULL,
                           tissue = NULL, # example of using tissue in neighbor definition
                           roundcortozero = 1e-4,
                           max_cells = 10000,
                           verbose = FALSE)
#saveRDS(ligstep1, file = "processed_data/ligstep1.RDS")
#ligstep1 = readRDS(file = "processed_data/ligstep1.RDS")

# get modules:
lmods <- defineModules(condcor = ligstep1$condcor, 
                          env = ligstep1$env, 
                          min_module_size = 2,
                          max_module_size = 20,
                          resolution = 0.02,
                          corthresh = 0.1, 
                          min_module_cor = 0.1,
                          gene_weighting_rule = "inverse_sqrt")
lmods$modules

lscores =  scoreModules(counts = norm,
                        weights = lmods$weights,
                        neighbors = ligstep1$neighbors)

pheatmap(cor(lscores$scores_env))

# attribution:
ligattribution <- celltype_attribution(
  modulescores = lscores$scores_env,
  weights = lmods$weights,
  counts = norm,
  celltype = clust,
  neighbors = ligstep1$neighbors,
  nsub = 1000,
  verbose = TRUE)
#saveRDS(ligattribution, file = "processed_data/ligattribution.RDS")
#ligattribution = readRDS("processed_data/ligattribution.RDS")


p1 = pheatmap(ligattribution$involvescores,
         col = colorRampPalette(c("white", "darkblue"))(100),
         breaks = seq(0,0.8,length.out=101))
svg("results/fig 2b - attrib heatmap.svg", width = 8, height = 6)
pheatmap(ligattribution$involvescores[p1$tree_row$order, p1$tree_col$order],
         cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("white", "darkblue"))(100),
         breaks = seq(0,0.8,length.out=101), legend = FALSE)
dev.off()

svg("results/fig 2a - network.svg")
par(mar = c(0,0,0,0))
plotCorrelationNetwork(x = ligstep1$condcor,
                       modules = lmods$weightsdf, show_gene_names = T)
dev.off()

# spatial plots of ligand env scores:
for (name in colnames(lscores$scores_env)) {
  
  png(paste0("results/ligand modules env scores - ", name, ".png"), width = 3, height = 5.5, units= "in", res = 500)
  par(mar = c(0,0,2,0))
  plot(xy, pch =16, asp = 1, cex = 0.2, ylim = c(-0.51,4.61), 
       main = paste0(lmods$modules[[name]], collapse = ", "), cex.main = 0.65,
       col = viridis_pal(option = "B")(101)[
         1 + round(100 * pmin(lscores$scores_env[, name] / quantile(lscores$scores_env[, name], 0.95), 1))],
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  lines(c(2.1,3.1), c(-0.2,-0.2))
  text(2.6,-0.3,"1 mm", cex = 0.8)
  dev.off()
}

#### LR pair analysis -----------------------------------
# get condcor for just the lr genes:

# look at cor values for all LR pairs:
lrcors = lrrawcors = c()
for (name in rownames(lrpairs)) {
  lrcors[name] = step1$condcor[lrpairs[name, 1], lrpairs[name, 2]]
  lrrawcors[name] = rawcor[lrpairs[name, 1], lrpairs[name, 2]]
}

plot(lrrawcors, lrcors, asp = 1, col = 0)
abline(0,1)
abline(h = 0)
text(lrrawcors, jitter(lrcors, amount = 0.106), names(lrcors), cex = 0.5)

# histogram of cor values for all LR pairs, with lines at cor values for selected pairs
svg("results/fig2f - LR cor histogram.svg", height = 3, width = 4)
par(mar = c(4,4,0,1))
hist(lrcors, breaks = 50, main = "", xlab = "Conditional correlation",
     ylab = "Number of L-R pairs", col = alpha("dodgerblue4", 0.75))
text(0.13, 12, c("FCER2 - CR2"))
lines(rep(lrcors["FCER2A_CR2"], 2), c(0, 10))
dev.off()

lrcors[order(lrcors, decreasing = T)[1:10]]

# pick 1-2 pairs, and show transcript plots:

# make pseudomodules, then score:
lrmods = list()
for (i in 1:5) {
  name = names(lrcors)[order(lrcors, decreasing = T)[i]]
  lrmods[[name]] = c(lrpairs[name, 1], lrpairs[name, 2])
}

#### correlation network around a LR pair: 

usegenes = names(which(colSums(step1$condcor[c(gl, gr), ] > 0) > 0))
gr0 <- cellKlatch:::get_coregulation_network(cormat = step1$condcor[usegenes, usegenes], 
                                corthresh = 0.1)
svg("results/fig 2h - network around LR pair.svg", width = 3, height = 3)
par(mar = c(0,0,0,0))
igraph::plot.igraph(gr0, 
                    vertex.label = unlist(list(NA, names(igraph::V(gr0)))[1 + TRUE]), 
                    vertex.size = 1, 
                    vertex.color = "grey",
                    vertex.label.color = c("dodgerblue4", "red")[1 + is.element(usegenes, c(gl, gr))]) 
dev.off()


# make a module score:
gcscore = rowMeans(norm[, setdiff(usegenes, c("FCER2", "CR2"))])
png(paste0("results/ligand modules env scores - ", name, ".png"), width = 3, height = 5, units= "in", res = 300)
par(mar = c(0,0,0,0))
plot(xy, pch =16, asp = 1, cex = 0.2, ylim = c(2,6.5),
     col = viridis_pal(option = "B")(101)[
       1 + round(100 * pmin(gcscore / quantile(gcscore, 0.995), 1))],
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
dev.off()


#### find TLS's and see how B-cells in them are different:

name = "FCER2A_CR2"
gl = lrpairs[name, "ligand"]
gr = lrpairs[name, "receptor"]

# define TLS's by finding spatial clusters of B-cells:
db = dbscan::dbscan(xy[clust == "B-cell", ], eps = 0.02, minPts = 5)$cluster
table(table(db))

# look at B-cell spatial clusters:
plot(xy[use, ], pch = 16, cex = 0.2, 
     main = name, col = c("grey80", rep(brewer_pal(9, "Set3")(9), 3))[1+db])
for (id in unique(db)) {
  text(median(xy[use, 1][db == id]),
       median(xy[use, 2][db == id]), id)
}

# only keep the big b-cell clusters:
tlsclusts = setdiff(names(table(db))[table(db) > 100], 0)
in.tls = is.element(db, tlsclusts)

t.test(norm[use, gl] ~ in.tls)
0.566/0.22
t.test(norm[use, gr] ~ in.tls)
0.56/0.163
# delta method to get CI on the ratio:
ests = t.test(norm[use, gl] ~ in.tls)$est[c("mean in group TRUE", "mean in group FALSE")]

# function to get Conf interval on expression ratio:
getCI = function(est1, est2, se1, se2) {
  grad = matrix(c(1/est1, -est1/est2^2), 2)
  se = as.vector(sqrt(t(grad) %*% diag(c(se1, se2)^2) %*% grad))
  c(est1/est2, (est1 / est2) + c(-1.96, 1.96) * se)
}
# ligand result:
tt.l.yes = t.test(norm[use, gl][in.tls])
tt.l.no = t.test(norm[use, gl][!in.tls])
getCI(est1 = tt.l.yes$estimate, est2 = tt.l.no$estimate, 
      se1 = tt.l.yes$stderr, se2 = tt.l.no$stderr)
# receptor result:
tt.r.yes = t.test(norm[use, gr][in.tls])
tt.r.no = t.test(norm[use, gr][!in.tls])
getCI(est1 = tt.r.yes$estimate, est2 = tt.r.no$estimate, 
      se1 = tt.r.yes$stderr, se2 = tt.r.no$stderr)



# B-cell GC polygons:
gcpolys = list()
for (id in setdiff(names(which(table(db) > 100)), 0)) {
  gcpolys[[id]] = xy[clust == "B-cell", ][db == id, ][chull(xy[clust == "B-cell", ][db == id, ]), ]
}
name = "FCER2A_CR2"
gl = lrpairs[name, "ligand"]
gr = lrpairs[name, "receptor"]

png("results/fig2g - spatial plot of 2 genes.png", width = 3, height = 4, units= "in", res = 500)
par(mar = c(0,0,0,0))
par(bg = "black")
show = rowSums(raw[, c(gl, gr)]) > 0
o = order(rowSums(norm[, c(gl, gr)]))
plot(xy[o, ], 
     col = c("grey20", "yellow", "cyan", "red")[
       1 + (norm[, gl] > 2) + 2 * (norm[, gr] > 2)
     ][o],
     pch = 16, cex = 0.1 + 0.05 * (apply(norm[, c(gl, gr)], 1, max) > 2), asp = 1,
     ylim = c(3, 6), #ylim = c(-0.51,4.61),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for (i in 1:length(gcpolys)) {
  polygon(gcpolys[[i]], col = NULL, border = "white")
}
legend(1.725,6.3, pch = c(16,16,16,NA), lty = c(NA,NA,NA,1), cex = 0.7,  
       #col = c(rgb(1,0,0,1), rgb(0,0,1,1), rgb(1,0,1,1), "white"),
       col = c("yellow", "cyan", "red", "white"),
       legend = c(paste0("high ", gl), paste0("high ", gr), "high both", "Lymphoid structure"), text.col = "white")
lines(c(1.8, 2.8), rep(2.7,2), col = "white")
text(2.3,2.6,"1 mm", cex = 0.8, col = "white")
dev.off()

