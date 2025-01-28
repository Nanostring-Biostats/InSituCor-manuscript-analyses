# test impact of tuning parameters: subset size, 

library(pheatmap)
library(scales)
library(Matrix)
library(viridis)
library(InSituType)
library(InSituCor)

#### load data ---------------------------------

# first, download the data from figshare:
# https://figshare.com/articles/dataset/CosMx_6000-plex_colon_cancer_dataset/24164388

# load the data:
load("colon cancer dataset.RData")

# normalize:
norm = sweep(raw, 1, rowSums(raw), "/") * mean(rowSums(raw))
# identify genes with strongest correlations:
temp <- InSituCor:::calcSpatialCor(
  counts = norm, 
  conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust),
  neighbors = NULL,
  xy = xy,
  k = 50,
  radius = NULL,
  tissue = NULL, # example of using tissue in neighbor definition
  roundcortozero = 0,
  max_cells = nsub,
  verbose = FALSE)$condcor
diag(temp) <- 0
genes <- order(apply(temp, 1, max), decreasing = T)[1:50]
#set.seed(0)
#genes <- sample(1:ncol(norm), 100)

#### explore impact of subset size: ----------------------------

resultslist <- list()
nsubs <- c(1000, 2500, 5000, 10000, 20000, 50000)
for (nsub in nsubs) {
  print(nsub)
  for (i in 1:10) {
    shuffle <- sample(1:nrow(norm), nrow(norm))
    resultslist[[paste0(nsub, "_", i)]] <- InSituCor:::calcSpatialCor(
      counts = norm[shuffle, genes], 
      conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust)[shuffle, ],
      neighbors = NULL,
      xy = xy[shuffle, ],
      k = 50,
      radius = NULL,
      tissue = NULL, # example of using tissue in neighbor definition
      roundcortozero = 0,
      max_cells = nsub,
      verbose = FALSE)$condcor
    diag(resultslist[[paste0(nsub, "_", i)]]) = 0
  }
}

# gold standard: all cells:
goldcor <- InSituCor:::calcSpatialCor(
  counts = norm[, genes], 
  conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust),
  neighbors = NULL,
  xy = xy,
  k = 50,
  radius = NULL,
  tissue = NULL, # example of using tissue in neighbor definition
  roundcortozero = 0,
  max_cells = nrow(norm),
  verbose = FALSE)$condcor
diag(goldcor) = 0

save(goldcor, resultslist, file = "results/guardbanding - condcors under varied max_cells.RData")

# get MSE of estimates:
rmses <- c()
for (nsub in nsubs) {
  inds <- which(grepl(paste0(nsub, "_"), names(resultslist)))
  temp <- sapply(inds, function(ind) {
    resultslist[[ind]][upper.tri(resultslist[[ind]])]
  })
  rmses <- cbind(rmses, sqrt(rowMeans(sweep(temp, 1, goldcor[upper.tri(goldcor)], "-")^2)))
}
colnames(rmses) <- as.character(nsubs)
colMeans(rmses)
mean(abs(goldcor[upper.tri(goldcor)] - 0.2) < 0.01375*2)
tab = 0
for (i in 1:10) {
  ind = paste0("5000_", i)
  tab = tab + table((goldcor[upper.tri(goldcor)] < 0.1), 
                    (resultslist[[ind]][upper.tri(resultslist[[ind]])] < 0.1))
}
(sum(tab)-sum(diag(tab))) / sum(tab)

svg("results/guardbanding max_cells.svg")
boxplot(rmses, ylab = "rMSE of conditional correlation estimates", xlab = "Subset size (max_cells argument)")
dev.off()

#### explore impact of neighborhood calc: ----------------------------



temp <- calcSpatialCor(counts = norm[, 1:10], 
                        conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust),
                        neighbors = NULL,
                        xy = xy,
                        k = NULL,
                        radius = 0.1,
                        tissue = NULL, # example of using tissue in neighbor definition
                        roundcortozero = 0.1,
                        max_cells = 1000,
                        verbose = FALSE)$neighbors

nresultslist <- list()
for (k in c(10,25, 50, 100, 200)) {
  nresultslist[[paste0("k", k)]] <- InSituCor:::calcSpatialCor(
    counts = norm[, genes], 
    conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust),
    neighbors = NULL,
    xy = xy,
    k = k,
    radius = NULL,
    roundcortozero = 0,
    max_cells = 10000,
    verbose = FALSE)$condcor
}
for (rad in c(0.02, 0.03, 0.05, 0.1, 0.2, 0.3)) {
  nresultslist[[paste0("Radius = ", rad)]] <- InSituCor:::calcSpatialCor(
    counts = norm[, genes], 
    conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust),
    neighbors = NULL,
    xy = xy,
    k = NULL,
    radius = rad,
    roundcortozero = 0,
    max_cells = 10000,
    verbose = FALSE)$condcor
}
for (i in 1:length(nresultslist)) {
  diag(nresultslist[[i]]) = 0
}
save(nresultslist, file = "results/guardbanding neighborhoods.RData")


pairwise <- sapply(nresultslist, function(mat){
  mat[upper.tri(mat)]
})

dev.off()
svg("results/neighborhood heatmap unthresholded.svg")
pheatmap(pairwise, cluster_cols = F, 
         col = colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(100),
         breaks = seq(-0.5,0.5,length.out=101))
dev.off()
#svg("results/neighborhood heatmap thresholded.svg")
png("results/neighborhood heatmap thresholded.png", units = "in", width = 7, height = 4, res = 400)
pheatmap(t(pairwise * (abs(pairwise) > 0.1)), cluster_rows = F, 
         col = colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(100),
         breaks = seq(-0.5,0.5,length.out=101))
dev.off()
#### 


#### explore impact of corthresh: --------------------------------

# get condcor:
step1 <- calcSpatialCor(counts = norm, 
                        conditionon = cbind(annot[, c("nCount_RNA", "neg")], clust),
                        xy = xy,
                        k = 50,
                        roundcortozero = 0.025,
                        max_cells = 10000,
                        verbose = FALSE)

modlist = list()
threshs <- c(0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25)
for (thresh in threshs) {
  modlist[[as.character(thresh)]] <- defineModules(
    condcor = step1$condcor,
    env = step1$env,
    min_module_size = 2,
    max_module_size = 20,
    resolution = 0.02,
    corthresh = thresh, 
    min_module_cor = 0.1,
    gene_weighting_rule = "inverse_sqrt")
}

# par(mfrow = c(3,3))
# for (name in names(modlist)) {
#   barplot(table(table(modlist[[name]]$weightsdf$module)), main = name)
# }

# par(mfrow = c(2,2))
# for (name in names(modlist)) {
#   set.seed(0)
#   obj <- plotCorrelationNetwork(x = step1$condcor,
#                          modules = modlist[[name]]$weightsdf, show_gene_names = F)
# }

for (i in 1:length(modlist)) {
  modlist[[i]]$membership <- modlist[[i]]$weightsdf$module[match(colnames(norm), modlist[[i]]$weightsdf$gene)]
}


for (i in 2:length(modlist)) {
  png(paste0("results/guardbanding corthresh - ", names(modlist)[i], ".png"), width = 6, height = 6, units = "in", res = 300)
  pheatmap(pmin(table(modlist[[1]]$membership, modlist[[i]]$membership), 10), 
           col = colorRampPalette(c("grey98", "darkorchid4", "black"))(13)[-(2:3)],
           breaks = 0:11, show_rownames = F, show_colnames = F, 
           main = paste0(
           "Number of genes shared between clusters from
   corthresh = 0.025 (rows) and ", names(modlist)[i], " (columns)"))
  dev.off()
}

#### OLD BELOW

# look at SD of estimates:
sds <- means <- c()
for (nsub in nsubs) {
  inds <- which(grepl(paste0(nsub, "_"), names(resultslist)))
  temp <- sapply(inds, function(ind) {
    resultslist[[ind]][upper.tri(resultslist[[ind]])]
  })
  sds <- cbind(sds, apply(temp, 1, sd))
  means <- cbind(means, apply(temp, 1, mean))
}
colnames(sds) <- as.character(nsubs)
pheatmap(sds, cluster_cols = F)

plot(means[,1], sds[,1])
plot(means[,6], sds[,6])

svg("results/guardbanding max_cells.svg")
boxplot(sds, ylab = "SD of conditional correlation estimates", xlab = "Subset size (max_cells argument)")
dev.off()



#### old below
# gold standard: mean of n = 50000
goldcor <- (resultslist$`50000_1` + resultslist$`50000_2` + resultslist$`50000_3`) / 3

mat <- goldcor[upper.tri(goldcor)]
for (i in 1:length(resultslist)) {
  mat <- cbind(mat, resultslist[[i]][upper.tri(goldcor)])
}
pheatmap(mat, cluster_cols = F)

