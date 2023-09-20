if (FALSE) {
  # environment setup:
  install.packages("Rcpp")
  devtools::install_github("https://github.com/Nanostring-Biostats/insitutype")
  devtools::install("../InSituCor")
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

#### load data ---------------------------------

# first, download the data from figshare:
# https://figshare.com/articles/dataset/CosMx_6000-plex_colon_cancer_dataset/24164388


# load the data:
load("colon cancer dataset.RData")

# normalize:
norm = sweep(raw, 1, rowSums(raw), "/") * mean(rowSums(raw))

#### run whole SPARC pipeline at once to test runtime ---------------------------

tic()
res <- insitucor(counts = norm, 
                 conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust),
                 celltype = clust,
                 neighbors = NULL,
                 xy = xy,
                 k = 50)
toc()  #2.4945

#### run SPARC step by step ---------------------------------

step1 <- calcSpatialCor(counts = norm, 
                        conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust),
                        neighbors = NULL,
                        xy = xy,
                        k = 50,
                        radius = NULL,
                        tissue = NULL, # example of using tissue in neighbor definition
                        roundcortozero = 0.1,
                        max_cells = 10000,
                        verbose = FALSE)
#saveRDS(step1, file = "processed_data/step1.RDS")
#step1 = readRDS(file = "processed_data/step1.RDS")

# derive Modules
modules <- defineModules(condcor = step1$condcor,
                         env = step1$env,
                         min_module_size = 2,
                         max_module_size = 20,
                         resolution = 0.02,
                         corthresh = 0.1, 
                         min_module_cor = 0.1,
                         gene_weighting_rule = "inverse_sqrt")
#saveRDS(modules, file = "processed_data/modules.RDS")
#modules <- readRDS("processed_data/modules.RDS")

# get scores:
scores <- scoreModules(counts = norm,
                       weights = modules$weights,
                       neighbors = step1$neighbors)
#saveRDS(scores, file = "processed_data/scores.RDS")
#scores <- readRDS("processed_data/scores.RDS")

#### attribution analysis -----------------------------------------------

tic()
attribution <- celltype_attribution(
  modulescores = scores$scores_env,
  weights = modules$weights,
  counts = norm,
  celltype = clust,
  neighbors = step1$neighbors,
  nsub = 1000,
  verbose = TRUE)
toc()
#saveRDS(attribution, file = "processed_data/attribution.RDS")
#attribution = readRDS("processed_data/attribution.RDS")



#### plots for manuscript --------------------------------------------

# cell type map:
cols = InSituType::colorCellTypes(freqs = table(clust), palette = "brewers")

png("results/fig 1a - cell type map.png", width = 3, height = 5, units= "in", res = 300)
par(mar = c(0,0,0,0))
plot(xy, pch = 16, col = cols[clust], cex = 0.1, asp = 1, ylim = c(-0.51,4.61),
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(c(2.1,3.1), c(-0.2,-0.2))
text(2.6,-0.3,"1 mm", cex = 0.8)
dev.off()

png("results/fig 1a - tall cell type map.png", width = 3, height = 6.5, units= "in", res = 300)
par(mar = c(0,0,0,0))
plot(xy, pch = 16, col = cols[clust], cex = 0.1, asp = 1, 
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(c(2.1,3.1), c(-0.2,-0.2))
text(2.6,-0.3,"1 mm", cex = 0.8)
dev.off()

svg("results/fig 1a - cell type legend.svg", width = 3, height = 3)
par(mar = c(0,0,0,0))
frame()
o = order(names(cols), decreasing = F)
legend("center", pch = 16, col = cols[o], legend = gsub(" naive", "", gsub(" memory", "", names(cols)[o])), cex = 0.55, ncol = 2)
dev.off()

# cartoon of neighbors:
ind = 10000
use = ((xy[, 1] > xy[ind,1] - 0.1) & (xy[, 1] < xy[ind,1] + 0.1)) & 
  ((xy[, 2] > xy[ind,2] - 0.1) & (xy[, 2] < xy[ind,2] + 0.1))
ind.neighbors = names(which(step1$neighbors[ind, ]> 1e-6))

svg("results/fig 1b - nearest neighbors.svg", width = 2, height = 2)
par(mar = c(0,0,0,0))
plot(xy[use, ], pch = 16, col = cols[clust[use]], cex = 0.8, asp = 1,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
points(xy[ind, 1], xy[ind, 2], cex = 1, pch = 16)
points(xy[ind.neighbors, 1], xy[ind.neighbors, 2], cex = 1, pch = 1, col = "blue")
legend("topright", pch = c(16, 1), col = c("black", "blue"), 
       legend = c("given cell", "nearest neighbors"), cex = 0.5)
dev.off()

# heatmap of environment matrix:
set.seed(0)
subr = sample(1:nrow(step1$env), 500)
subc = order(colMeans(step1$env), decreasing = T)[1:100]

mat = step1$env[subr, subc]
mat = sweep(mat,2,apply(mat,2,max), "/")
p1 = pheatmap(mat)
dev.off()
png("results/fig 1c - env matrix.png", width = 3, height = 3, units ="in",res = 400)
pheatmap(mat[p1$tree_row$order, p1$tree_col$order], 
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F, legend = FALSE,
         col = viridis_pal(option = "B")(100))
dev.off()

# heatmap of cormat
rawcor = cor(step1$env)
# saveRDS(rawcor, file = "processed_data/rawcor.RDS")
#rawcor = readRDS("processed_data/rawcor.RDS")

set.seed(0)
inds = sample(1:nrow(rawcor), 300)
  
hc1 = hclust(dist(rawcor[inds, inds]))
png("results/fig 1e - raw cormat.png", width = 3, height = 3, units ="in",res = 400)
pheatmap(rawcor[inds, inds][hc1$order, hc1$order], cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("darkblue",'blue', "white","red","darkred"))(100),
         breaks = seq(-0.6,0.6,length.out = 101),
         #col = colorRampPalette(c('blue', "white","red"))(100),
        #breaks = seq(-0.4,0.4,length.out = 101),
        show_rownames = F, show_colnames = F, legend = FALSE)
dev.off()

# heatmap of cond cor:
hc2 = hclust(dist(step1$condcor[inds, inds]))
png("results/fig 1f - conditional cormat.png", width = 3, height = 3, units ="in",res = 400)
pheatmap(step1$condcor[inds, inds][hc2$order, hc2$order], cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("darkblue",'blue', "white","red","darkred"))(100),
         breaks = seq(-0.6,0.6,length.out = 101),
         show_rownames = F, show_colnames = F, legend = FALSE)
dev.off()

# heatmap of env confounders
condmat = cellKlatch:::build_conditional_matrix(cbind(annot[, c("nCount_RNA", "neg")], clust))
colnames(condmat)[1:2] = c("total counts", "negprobe counts")
condmat[, "total counts"] = condmat[, "total counts"] / mean(condmat[, "total counts"])

mat = condmat[1:500, order(colMeans(condmat[1:500, ]), decreasing = T)[1:20]]
p1 = pheatmap(mat)
dev.off()
png("results/fig 1d - conditioning matrix.png", width = 3, height = 3, units ="in",res = 400)
pheatmap(mat[p1$tree_row$order, p1$tree_col$order], cluster_rows = F, cluster_cols = F,
         #fontsize_col = 6,
         show_colnames = F,
         col = colorRampPalette(c("grey80","darkviolet"))(100),
         show_rownames = F, legend = FALSE)
dev.off()


## show examples of marker genes losing conditional cor:

# find marker genes:
profiles = InSituType::Estep(counts = raw, clust = clust, neg = annot$neg)
markers = list()
for (cell in unique(clust)) {
  markers[[cell]] = rownames(profiles)[
    order(profiles[, cell] / apply(profiles[, setdiff(colnames(profiles), cell)], 1, max) * sqrt(profiles[, cell]),
                          decreasing = T)[1:4]]
}


# get non-rounded conditional cor of only marker genes:
markerstep1 <- spatialcor(counts = norm[, unlist(markers)], 
                          conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust),
                          neighbors = NULL,
                          xy = xy,
                          k = 50,
                          radius = NULL,
                          tissue = NULL, # example of using tissue in neighbor definition
                          roundcortozero = 1e-4,
                          max_cells = 10000,
                          verbose = FALSE)

unroundedstep1 <- spatialcor(counts = norm, 
                                conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust),
                                neighbors = NULL,
                                xy = xy,
                                k = 50,
                                radius = NULL,
                                tissue = NULL, # example of using tissue in neighbor definition
                                roundcortozero = 1e-2,
                                max_cells = 10000,
                                verbose = FALSE)
saveRDS(unroundedstep1, file = "processed_data/unroundedstep1.RDS")
#unroundedstep1 = readRDS(file = "processed_data/unroundedstep1.RDS")


# look at cors:
gpairs = list(c("MS4A1", "CD19"),   #, "BLK", "BCL2"
              c("CD68", "CD163"),
              c("CD3E", "CD3D"), 
              c("COL5A1", "COL5A2"))

for (i in 1:length(gpairs)) {
  print(c(rawcor[gpairs[[i]][1], gpairs[[i]][2]],
          uncorrelatedstep1$condcor[gpairs[[i]][1], gpairs[[i]][2]]))
}
for (i in 1:length(gpairs)) {
  print(rawcor[gpairs[[i]], gpairs[[i]]])
}

for (i in 1:length(gpairs)) {
  print(step1$condcor[gpairs[[i]], gpairs[[i]]])
}

for (i in 1:length(gpairs)) {
  print(markerstep1$condcor[gpairs[[i]], gpairs[[i]]])
}


## summary of cors before and after:

# summarize how many cors are lost after conditioning:
top5kthresh = quantile(rawcor[upper.tri(rawcor)], 1 - 5000/sum(upper.tri(rawcor)))
top5k = (rawcor > top5kthresh) & upper.tri(rawcor)
sum(step1$condcor[top5k] > 0.2)
sum(step1$condcor[top5k] <= 0.2)
sum(top5k)


# plot cors before and after:
ynudge = c(0.035,-0.05,0.05,0.035)
png("results/fig 1g - scatterplot of cors before and after.png", width = 4, height = 5, units = "in", res = 300)
par(mar = c(5,4.5,0,0))
plot(rawcor[upper.tri(rawcor)], uncorrelatedstep1$condcor[upper.tri(step1$condcor)], 
     pch = 16, cex = 0.5,
     xlab = "Correlation of gene pairs\nin environment matrix", 
     ylab = "Correlation conditional on confounders",
     cex.lab = 1.1,
     xlim = c(0,1), ylim = c(-0.1,1),
     col = alpha("dodgerblue4", 0.25))
abline(0,1)
abline(h = 0, col = 1)
abline(v = 0, col = 1)
for (i in 1:length(gpairs)) {
  points(rawcor[gpairs[[i]][1], gpairs[[i]][2]],
         step1$condcor[gpairs[[i]][1], gpairs[[i]][2]],
         pch = 16, col = "darkred")
  text(rawcor[gpairs[[i]][1], gpairs[[i]][2]],
       step1$condcor[gpairs[[i]][1], gpairs[[i]][2]] + ynudge[i],
       labels = paste0(gpairs[[i]][1], "/", gpairs[[i]][2]),
       col = "darkred", cex = 0.75)
}
rect(top5kthresh, 0.205, 1, 1, border = "orange")
rect(top5kthresh, -0.1, 1, 0.195, border = "darkviolet")
text(0.85,-0.05,paste0(sum(step1$condcor[top5k] <= 0.2), " pairs"), col = "darkviolet")
text(0.85,0.85,paste0(sum(step1$condcor[top5k] > 0.2), " pairs"), col = "orange")
dev.off()





#### spatial plots -----------------------
name = "SFRP2_CCL18_17"
png(paste0("results/fig 1h - ", name, "env scores.png"), width = 3, height = 5, units= "in", res = 300)
par(mar = c(0,0,0,0))
plot(xy, pch =16, asp = 1, cex = 0.2, ylim = c(-0.51,4.61),
     col = viridis_pal(option = "B")(101)[
       1 + round(100 * pmin(scores$scores_env[, name] / quantile(scores$scores_env[, name], 0.95), 1))],
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(c(2.1,3.1), c(-0.2,-0.2))
text(2.6,-0.3,"1 mm", cex = 0.8)
#rect(0.2, 2.5, 1, 3.5, border = "cyan", lwd = 2)
rect(1.4, 2.5, 3, 4.5, border = "dodgerblue1", lwd = 2)
dev.off()

png(paste0("results/fig 1h short - ", name, "env scores.png"), width = 3, height = 4.25, units= "in", res = 300)
par(mar = c(0,0,0,0))
plot(xy, pch =16, asp = 1, cex = 0.2, ylim = c(0,4.61),
     col = viridis_pal(option = "B")(101)[
       1 + round(100 * pmin(scores$scores_env[, name] / quantile(scores$scores_env[, name], 0.95), 1))],
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(c(0.1,1.1), rep(4.5,2), col = "grey50")
text(0.6,4.45,"1 mm", cex = 0.8, col = "grey50")
rect(1.4, 2.5, 3, 4.5, border = "dodgerblue1", lwd = 2)
dev.off()


png(paste0("results/fig 1i - ", name, "env scores zoom.png"), width = 4, height = 5, units= "in", res = 300)
par(mar = c(0,0,0,0))
plot(xy, pch =16, asp = 1, 
     xlim = c(0.2,1), ylim = c(2.5,3.5),
     cex = .5,
     col = viridis_pal(option = "B")(101)[
       1 + round(100 * pmin(scores$scores_env[, name] / quantile(scores$scores_env[, name], 0.95), 1))],     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
dev.off()


png(paste0("results/fig 1j - ", name, "sc scores zoom OLD.png"), width = 4, height = 5, units= "in", res = 300)
par(mar = c(0,0,0,0))
plot(xy, pch =16, asp = 1, 
     xlim = c(0.2,1), ylim = c(2.5,3.5),
     cex = .2 + 1.5 * pmin(scores$scores_sc[, name] / quantile(scores$scores_sc[, name],0.99),1), 
     col = cols[clust],
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(c(0.8,1), rep(2.6,2))
text(0.9,2.57, "0.2 mm")
dev.off()

png(paste0("results/fig 1j - ", name, "sc scores zoom.png"), width = 4, height = 5, units= "in", res = 300)
par(mar = c(0,0,0,0))
plot(xy, pch =16, asp = 1, 
     xlim = c(1.4,3), ylim = c(2.5,4.5),
     cex = .2 + 1.0 * pmin(scores$scores_sc[, name] / quantile(scores$scores_sc[, name],0.99),1), 
     col = cols[clust],
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
lines(c(2.6,2.8), rep(2.6,2))
text(2.7,2.57, "0.2 mm")
showcells = c("CAF", "macrophage", "mast", "smooth.muscle.cell", "stromal.cell.type.4")
legend("bottomleft", pch = 16, col = cols[showcells], legend = showcells, cex = 0.9)
rect(1.55, 3.55, 1.95, 4.05, border = "dodgerblue1", lwd = 2)
dev.off()

### network diagram ------------------

svg("results/fig 1k - network.svg")
par(mar = c(0,0,0,0))
plotCorrelationNetwork(x = step1$condcor,
                       modules = modules$weightsdf, show_gene_names = F)
dev.off()

svg("results/fig 1k - network v2.svg", width = 8, height = 8)
par(mar = c(0,0,0,0))
plotCorrelationNetwork(x = step1$condcor,
                       modules = modules$weightsdf, show_gene_names = F, vertex_size = 3.5)
dev.off()

pdf("results/fig 1k - network w names.pdf")
par(mar = c(0,0,0,0))
plotCorrelationNetwork(x = step1$condcor,
                       modules = modules$weightsdf, show_gene_names = T)
dev.off()


### attribution plots --------------

svg("results/attribution matrix transposed.svg", height = 8, width = 14)
pheatmap(t(attribution$involvescores),
         border_color = NA, legend = F,
         col = colorRampPalette(c("white","darkblue"))(100),
         breaks = seq(0,0.9,length.out = 101))
         #fontsize_col = 6)
dev.off()

p1 = pheatmap(attribution$attributionmats[[name]], main = name,
              col = colorRampPalette(c("white","darkblue"))(100),
              breaks = seq(0,1,length.out=101))

svg(paste0("results/attribution for ", name, ".svg"), width = 6, height = 7)
pheatmap(attribution$attributionmats[[name]][p1$tree_row$order, p1$tree_col$order],
         border_color = NA, 
         cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("white","darkblue"))(100),
         breaks = seq(0,0.9,length.out = 101),
         fontsize_row = 10, fontsize_col = 8, legend = F)
dev.off()



pdf("results/attributionmats.pdf")
for (name in names(attribution$attributionmats)) {
  pheatmap(attribution$attributionmats[[name]],
           main = name,
           col = colorRampPalette(c("white","darkblue"))(100),
           breaks = seq(0,1,length.out=101))
}
dev.off()

