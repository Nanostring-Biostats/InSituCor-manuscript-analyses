library(viridis)
library(igraph)
library(uwot)
library(hexbin)
library(pheatmap)

#### load dataset -----------------------------

# (data available upon request)
load("~//NAS_data/lwu/testRun/SpatialTest/giotto_test/SMITAP_proj/MGH_WillHwang_PDAC_2023Jan/rnaRes/rna_postQC_fundamental_data.RData")
metadata <- read.csv("~//NAS_data/lwu/testRun/SpatialTest/giotto_test/SMITAP_proj/MGH_WillHwang_PDAC_2023Jan/rnaRes/cellType_summary/all_celltyping_final_clus.csv")

# how many cells per core?
table(metadata$TMA_coreID)
# cores to remove:
# 2B-2, 2B-3, 2B-8

# remove cores with low cell counts:
use <- !is.element(metadata$TMA_coreID, c("2B-2", "2B-3", "2B-8"))
spatCoord <- spatCoord[use, ]
metadata <- metadata[use, ]
rna_counts <- rna_counts[, use]

# arrange xy space:
source("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/refs/heads/Main/_code/condensing%20xy%20space/CondenseXY.R")
xy <- cbind(spatCoord$x_slide_mm, spatCoord$y_slide_mm)
set.seed(0)
xy <- condenseXY(xy, fov = paste0(metadata$fov, metadata$slide_ID_numeric), 
                 tissue = metadata$TMA_coreID, 
                  condensefovs = TRUE, condensetissues = TRUE, # which operations to perform
                   eps = 1.5, mindist = 0.1, fovbuffer = 0.2, # FOV condensing args
                   tissueorder = NULL, tissuebuffer = 0.2, widthheightratio = .85)
plot(xy, cex = 0.2, asp = 1, col = as.numeric(as.factor(metadata$slide_ID_numeric)))
for (tiss in unique(metadata$TMA_coreID)) {
  use <- metadata$TMA_coreID == tiss
  rect(min(xy[use,1]), min(xy[use,2]), max(xy[use,1]), max(xy[use,2]))
  text(median(range(xy[use,1])), max(xy[use,2])+0.2, tiss, cex = 0.5)
}

#### load individual condcors --------------------------------

# load individual core's condcor matrices:
loc <- "~//NAS_data/lwu/testRun/SpatialTest/giotto_test/SMITAP_proj/MGH_WillHwang_PDAC_2023Jan/rnaRes/perCore_spatialCor"
filenames <- dir(loc)[grepl("insituCor_res.rds", dir(loc))]
# remove low-count cores:
filenames <- setdiff(filenames, c("2B-2__R50__insituCor_res.rds",
                                  "2B-3__R50__insituCor_res.rds",
                                  "2B-8__R50__insituCor_res.rds"))
temp <- readRDS(paste0(loc, "/", filenames[1]))

condcors<- array(NA, dim = c(length(filenames), nrow(temp$condcor), ncol(temp$condcor)),
                 dimnames = list(filenames, rownames(temp$condcor), colnames(temp$condcor)))
for (name in filenames) {
  print(name)
  temp <- readRDS(paste0(loc, "/", name))
  condcors[name, ,] <- temp$condcor[dimnames(condcors)[[2]], dimnames(condcors)[[3]]]
  #saveRDS(temp$condcor, file = paste0("condcors/", substr(name, 1, gregexpr("_", name)[[1]][1] - 1), ".Rds"))
}

# summary stats:
maxcor <- apply(condcors, 2:3, max, na.rm = T)
meancor <- apply(condcors, 2:3, mean, na.rm = T)
qcor <- apply(condcors, 2:3, quantile, 17/19, na.rm = T)
lowqcor <- apply(condcors, 2:3, quantile, 3/19, na.rm = T)

diag(maxcor) = diag(meancor) = diag(qcor) = 0
diag(lowqcor) = 0
#sdcor <- apply(condcors, 2:3, sd)
#diag(maxcor) = diag(meancor) = diag(sdcor) = 0

#### define consensus network: high lowqcors -------------------------------------

consensus <- lowqcor > 0.2
consensus <- consensus[rowSums(consensus) != 0, colSums(consensus) != 0] * 1
# get a good layout:
xum <- uwot::umap(as.matrix(consensus), 
                  spread = 30, 
                  min_dist = 0.1,
                  n_neighbors = 15)
plot(xum, pch = 16, cex = 1)

# make an igraph from the adjacency matrix:
gr0 = igraph::graph_from_adjacency_matrix(consensus, weighted = TRUE,
                                          mode = "undirected")
V(gr0)$color = "black"
svg("results/consensus network.svg", width = 8, height = 8)
igraph::plot.igraph(gr0, vertex.size = 0.1)
dev.off()

# define a module from one promising cluster:
modgenes <- c("S100A6", "JUP", "ANXA2", "EZR", "CEACAM6", "KRT19", "KRT8", 
              "TMSB10", "TMSB4X", "LMNA", "PKM", "TACSTD2", "CLDN4", "CRIP1", "SPINT2")
score <- Matrix::colMeans(rna_counts[modgenes, ]) / colMeans(rna_counts)

# png("results/SP100A6 module spatial plot.png", width = 9, height = 9, units = "in", res = 300)
# par(mar = c(0,0,0,0))
# plot(xy, cex = 0.05, asp = 1, 
#      xlab = "", ylab = "", xaxt = "n", yaxt = "n",
#      col = viridis_pal(option="B")(101)[round(pmin(1 + 100 * score / quantile(score, 0.99)), 101)],
#      ylim = c(-0.5,9))
# for (tiss in unique(metadata$TMA_coreID)) {
#   use <- metadata$TMA_coreID == tiss
#   rect(min(xy[use,1]), min(xy[use,2]), max(xy[use,1]), max(xy[use,2]))
#   text(median(range(xy[use,1])), max(xy[use,2])+0.1, tiss, cex = 0.5)
# }
# lines(c(7.5,8.5), rep(-0.25,2))
# text(7.25, -0.25, "1 mm", cex = 0.8)
# dev.off()

# show only cancer cells:
png("results/SP100A6 module spatial plot - cancer only.png", width = 7.5, height = 9, units = "in", res = 400)
show1 = metadata$RNA_nbclust_clusters == "Malignancy"
par(mar = c(0,0,0,0))
plot(xy, cex = 0.01, asp = 1, col = "grey60",
     xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     ylim = c(-0.5,9.5))
points(xy[show1, ], cex = 0.05, 
       col = viridis_pal(option="B")(101)[round(pmin(1 + 100 * score[show1] / quantile(score, 0.99)), 101)])
for (tiss in unique(metadata$TMA_coreID)) {
  use <- metadata$TMA_coreID == tiss
  rect(min(xy[use,1]), min(xy[use,2]), max(xy[use,1]), max(xy[use,2]))
  text(median(range(xy[use,1])), max(xy[use,2])+0.1, tiss, cex = 0.5)
}
lines(c(6.5,7.5), rep(-0.25,2))
text(6.25, -0.25, "1 mm", cex = 0.8)
dev.off()

# show 10B-7
png("results/SP100A6 module spatial plot - 10B7 cancer only.png", width = 5, height = 5, units = "in", res = 500)
show1 = (metadata$TMA_coreID == "10B-7")
show = (metadata$RNA_nbclust_clusters == "Malignancy") & (metadata$TMA_coreID == "10B-7")
par(mar = c(0,0,0,0))
plot(xy[show1, ], cex = 0.05, asp = 1, 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     col = "grey60")
points(xy[show, ], pch = 16, cex = 0.2,
       col = viridis_pal(option="B")(101)[round(pmin(1 + 100 * score[show] / quantile(score, 0.99)), 101)])
lines(c(5,5.5), rep(3.2,2))
text(5.25, 3.3, "0.5 mm", cex = 0.8)
dev.off()


#### find highly variable gene pairs: -------------------------------
# identify pairs with strong positives and strong negatives. Say a q3/19?
bigrange <- (qcor > 0.4) & (lowqcor < 0.05)
sum(bigrange) / 2

# heatmap of highly variable pairs:
keepvec <- bigrange[upper.tri(bigrange)]
gene1 <- matrix(rep(rownames(maxcor), ncol(maxcor)), nrow(maxcor))[upper.tri(maxcor)][keepvec]
gene2 <- matrix(rep(colnames(maxcor), each = nrow(maxcor)), nrow(maxcor))[upper.tri(maxcor)][keepvec]
paircors.range <- c()
for (name in dimnames(condcors)[[1]]) {
  print(name)
  paircors.range <- cbind(paircors.range, condcors[name,,][upper.tri(maxcor)][keepvec])
}
rownames(paircors.range) <- paste0(gene1, "_", gene2)
colnames(paircors.range) <- dimnames(condcors)[[1]]
colnames(paircors.range) <- gsub("__R50__insituCor_res.rds", "", colnames(paircors.range))


# classify tumor morphologies based on gland size:   
small <- c("2B-8", "10B-11", "10B-4", "2B-2", "2B-3", "10B-7", "10B-8")
large <- c("10B-1", "10B-10", "10B-2", "10B-9", "2B-9", "2B-4", "2B-6", "2B-10") 
med <- c("10B-12", "10B-3", "10B-6", "2B-1", "2B-11", "2B-12", "2B-7") 

coreannot <- data.frame(gland_size = c("small", "med", "large")[
  1 * is.element(colnames(paircors.range), small) + 
    2 * is.element(colnames(paircors.range), med) + 
    3 * is.element(colnames(paircors.range), large)])
rownames(coreannot) = colnames(paircors.range)

# tall heatmap
pdf("results/discordant pairs heatmap.pdf", height = 20)
p1 = pheatmap(paircors.range, 
              col = colorRampPalette(c("darkblue","blue","white","red","darkred"))(100),
              breaks = seq(-0.5,0.5,length.out=101), 
              fontsize_row = 5,
              annotation_col = coreannot,
              annotation_colors = list("gland_size" = c("small" = "cornflowerblue", "med" = "grey", "large" = "orange")))
dev.off()

# show whole heatmap for figure:
svg("results/discordant pairs heatmap - whole.svg", height = 6, width = 4.3)
pheatmap(paircors.range[p1$tree_row$order, p1$tree_col$order], 
         cluster_rows = F, cluster_cols = F,
         show_rownames = F, show_colnames = F,
         col = colorRampPalette(c("darkblue","blue","white","red","darkred"))(100),
         breaks = seq(-0.5,0.5,length.out=101), fontsize_row = 5,
         annotation_col = coreannot,
         annotation_colors = list("gland_size" = c("small" = "cornflowerblue", "med" = "grey", "large" = "orange")),
         annotation_names_col = FALSE)
dev.off()
# show selected subsets for figure:
pairs = c("COX1_SPINT2","COX2_SPINT2","COX1_S100A6","COX2_S100A6","COX1_S100A4","COX2_S100A4","COX1_S100A14","COX2_S100A14")
p1 <- pheatmap(paircors.range[pairs, ], 
               cluster_rows = T, cluster_cols = T,
               #show_rownames = F, show_colnames = F,
               col = colorRampPalette(c("darkblue","blue","white","red","darkred"))(100),
               breaks = seq(-0.75,0.75,length.out=101),
               annotation_col = coreannot,
               annotation_colors = list("gland_size" = c("small" = "cornflowerblue", "med" = "grey", "large" = "orange")),
               annotation_names_col = FALSE)
dev.off()
svg("results/discordant pairs heatmap - selected.svg", height = 2.5, width = 6)
pheatmap(paircors.range[pairs[p1$tree_row$order], p1$tree_col$order], 
         cluster_rows = F, cluster_cols = F,
         #show_rownames = F, show_colnames = F,
         col = colorRampPalette(c("darkblue","blue","white","red","darkred"))(100),
         breaks = seq(-0.75,0.75,length.out=101),
         annotation_col = coreannot,
         annotation_colors = list("gland_size" = c("small" = "cornflowerblue", "med" = "grey", "large" = "orange")),
         annotation_names_col = FALSE)
dev.off()

#### spatial plots of modules with discordant results ------------------------------
# sets of similar genes:
sets <- list(
  set1 = c("COX1", "SPINT2", "COX2", "S100A6", "S100A4", "S100A14"),
  set2 = c("LMNA", "S100A4", "CCND1", "KRT19", "NQO1", "S100P", "CEACAM6", "MLPH", "CPB1", "SPINK1", "PERP"),
  set.small = c("KRT16", "KRT17", "MUC4", "KRT19", "S100A10", "ITGA2", "KRT7", "YWHAZ", "HMGA1", "LAMC2"), #selected for high in small glands
  set.large = c("TFF1", "TFF2", "S100P", "AGR2", "IER3", "LYZ", "CEACAM6"))


for (i in 1:length(sets)) {
  png(paste0("results/set", i, " module spatial plot - cancer only.png"), width = 7.5, height = 9, units = "in", res = 400)
  
  show = (metadata$RNA_nbclust_clusters == "Malignancy")
  setscore <- Matrix::colMeans(rna_counts[sets[[i]], ]) / Matrix::colMeans(rna_counts)
  par(mar = c(0,0,0,0))
  plot(xy, cex = 0.01, asp = 1, 
       col = "grey60",
       xlab = "", ylab = "", xaxt = "n", yaxt = "n",
       ylim = c(-0.5,9.5))
  points(xy[show, ], cex = 0.05, 
         col = viridis_pal(option="B")(101)[round(pmin(1 + 100 * setscore[show] / quantile(setscore[show], 0.99)), 101)])
  for (tiss in unique(metadata$TMA_coreID)) {
    use <- metadata$TMA_coreID == tiss
    rect(min(xy[use,1]), min(xy[use,2]), max(xy[use,1]), max(xy[use,2]))
    text(median(range(xy[use,1])), max(xy[use,2])+0.1, tiss, cex = 0.5)
  }
  lines(c(6.5,7.5), rep(-0.25,2))
  text(6.25, -0.25, "1 mm", cex = 0.8)
  dev.off()
}

# zoom in on selected tumors:
i = 1
tissues = c("2B-7", "2B-9", "2B-12")
tempuse <- is.element(metadata$TMA_coreID, tissues)
tempxy <- condenseXY(xy[tempuse, ], fov = paste0(metadata$fov, metadata$slide_ID_numeric)[tempuse], 
                      tissue = metadata$TMA_coreID[tempuse], 
                      condensefovs = TRUE, condensetissues = TRUE, # which operations to perform
                      eps = 1.5, mindist = 0.1, fovbuffer = 0.2, # FOV condensing args
                      tissueorder = NULL, tissuebuffer = 0.2, widthheightratio = 3)
png(paste0("results/set1 module spatial plot - focused.png"), width = 8, height = 5, units = "in", res = 500)
show = (metadata$RNA_nbclust_clusters[tempuse] == "Malignancy")
setscore <- Matrix::colMeans(rna_counts[sets[[i]], ])[tempuse] / Matrix::colMeans(rna_counts)[tempuse]
par(mar = c(0,0,0,0))
plot(tempxy, cex = 0.01, asp = 1, 
     col = "grey60",
     xlab = "", ylab = "", xaxt = "n", yaxt = "n")
points(tempxy[show, ], cex = 0.2, 
       col = viridis_pal(option="B")(101)[round(pmin(1 + 100 * setscore[show] / quantile(setscore[show], 0.99)), 101)])
for (tiss in tissues) {
  use <- metadata$TMA_coreID[tempuse] == tiss
  rect(min(tempxy[use,1]), min(tempxy[use,2]), max(tempxy[use,1]), max(tempxy[use,2]))
  text(median(range(tempxy[use,1])), max(tempxy[use,2])+0.1, tiss, cex = 0.8)
}
lines(c(1.8,2.8), rep(0.15,2))
text(2.3, .2, "1 mm", cex = 0.8)
dev.off()



