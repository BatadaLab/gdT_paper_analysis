library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(scID)
library(ggpubr)
library(plyr)
library(forcats)
library(easyGgplot2)
source("~/Google Drive/bin/batadalab_scrnaseq_utils.R")

# -------------------------------------------------------
# Load data 
# -------------------------------------------------------

# BC1
gem1 <- get_scrnaseq_data("p01e01_custom") %>%
  GetAssayData(slot = "counts") %>%
  as.data.frame() %>%
  set_colnames(paste("BC1", colnames(.), sep = "_"))

# BC2
gem2 <- get_scrnaseq_data("p01e02_bedtools") %>%
  GetAssayData(slot = "counts") %>%
  as.data.frame() %>%
  set_colnames(paste("BC2", colnames(.), sep = "_"))

# -------------------------------------------------------
# Set up objects
# -------------------------------------------------------
sobj1 <- CreateSeuratObject(counts = gem1, project = "BC1", min.cells = 5)
sobj1$stim <- "BC1"
sobj1 <- subset(sobj1, subset = nFeature_RNA > 500)
sobj1 <- NormalizeData(sobj1, verbose = FALSE)
sobj1 <- FindVariableFeatures(sobj1, selection.method = "vst", nfeatures = 2000)

sobj2 <- CreateSeuratObject(counts = gem2, project = "BC2", min.cells = 5)
sobj2$stim <- "BC2"
sobj2 <- subset(sobj2, subset = nFeature_RNA > 500)
sobj2 <- NormalizeData(sobj2, verbose = FALSE)
sobj2 <- FindVariableFeatures(sobj2, selection.method = "vst", nfeatures = 2000)

# -------------------------------------------------------
# Integrate
# -------------------------------------------------------
anchors <- FindIntegrationAnchors(object.list = list(sobj1, sobj2), dims = 1:20)
sobj.combined <- IntegrateData(anchorset = anchors, dims = 1:20)

# -------------------------------------------------------
# Perform integrated analysis
# -------------------------------------------------------
DefaultAssay(sobj.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
sobj.combined <- ScaleData(sobj.combined, verbose = FALSE, features = rownames(sobj.combined))
sobj.combined <- RunPCA(sobj.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
sobj.combined <- RunUMAP(sobj.combined, reduction = "pca", dims = 1:20)
sobj.combined <- FindNeighbors(sobj.combined, reduction = "pca", dims = 1:20)
sobj.combined <- FindClusters(sobj.combined, resolution = 0.5)

# Remove TRDC- population from cluster 4
filtered_cells <- read.delim(
  "data/raw/BC_filtered_cells.txt", 
  stringsAsFactors = F, 
  header = F
)

sobj.combined <- subset(
  sobj.combined, 
  cells = filtered_cells$V1
)





pdf("Figures/Figure3/Fig3A.pdf", width = 2, height = 2)
DimPlot(
  sobj.combined, 
  reduction = "umap", 
  label = TRUE, 
  pt.size = 0.01
) + 
  NoAxes() + 
  NoLegend()
dev.off()

pdf("Figures/Figure3/Fig3B.pdf", width = 4, height = 2.2)
DimPlot(
  sobj.combined, 
  reduction = "umap", 
  group.by = "stim", 
  pt.size = 0.01, 
  split.by = 'stim'
) + 
  NoAxes() + 
  NoLegend()
dev.off()


# -------------------------------------------------------
# Delta1/Delta2 genetics
# -------------------------------------------------------
# BC1
vdj_BC1 <- read.delim(
  "data/raw/VDJ_data/BC1_vdj_list.txt", 
  stringsAsFactors = F, 
  header = F
) %>%
  mutate(V1 = unlist(lapply(V1, function(x) strsplit(x, "#|_")[[1]][2]))) %>%
  mutate(V1 = paste("BC1", V1, sep = "_")) %>%
  filter(V1 %in% colnames(gem1))

TRDV2_BC1 <- vdj_BC1$V1[grep("TRDV2", vdj_BC1$V2)]
TRDC_BC1 <- vdj_BC1$V1[grep("TRDC", vdj_BC1$V2)]
TRGV9_BC1 <- vdj_BC1$V1[grep("TRGV9", vdj_BC1$V2)]

# BC2
vdj_BC2 <- read.delim(
  "data/raw/VDJ_data/BC2_vdj_list.txt", 
  stringsAsFactors = F,
  header = F
) %>%
  mutate(V1 = lapply(V1, function(x) strsplit(x, "#|_")[[1]][2])) %>%
  mutate(V1 = paste("BC2", V1, sep = "_")) %>%
  filter(V1 %in% colnames(gem2))
TRDV2_BC2 <- vdj_BC2$V1[grep("TRDV2", vdj_BC2$V2)]
TRDC_BC2 <- vdj_BC2$V1[grep("TRDC", vdj_BC2$V2)]
TRGV9_BC2 <- vdj_BC2$V1[grep("TRGV9", vdj_BC2$V2)]

# TRDC - Figure 3C(left)
vdj_labels <- rep("TRDC-", ncol(sobj.combined))
names(vdj_labels) <- colnames(sobj.combined)
vdj_labels[c(TRDC_BC1, TRDC_BC2)] <- "TRDC+"
df <- data.frame(sobj.combined@reductions$umap@cell.embeddings)
df$label <- vdj_labels[rownames(df)]
table(df$label)
dropouts <- rownames(df)[which(df$label == "TRDC-")]
tiff("Figures/Figure3/Fig3C_TRDC.tiff")
ggplot(
  df[dropouts, ], 
  aes(x=UMAP_1, y=UMAP_2, color=factor(label), fill=factor(label))
) + 
  geom_point(size=1) +
  theme_void() + 
  theme(legend.position="none") + 
  scale_color_manual(values=c("lightgrey","black")) +
  scale_fill_manual(values=c("lightgrey","black")) + 
  geom_point(data=df[c(TRDC_BC1, TRDC_BC2), ], colour="black", size=1)
dev.off()
# # alternative TRDC plot
# p <- FeaturePlot(
#   sobj.combined, 
#   "TRDC.NT-026437.11.3", 
#   pt.size = 0.01, 
#   combine = FALSE, 
#   min.cutoff = 0, 
#   max.cutoff = 3
# )
# pdf("Figures/Figure3/Fig3C_TRDC.pdf", width = 2, height = 2.3)
# plot(p[[1]] + NoAxes() + NoLegend())
# dev.off()

# TRDV2 - Figure 3C(middle)
vdj_labels <- rep("nothing", ncol(sobj.combined))
names(vdj_labels) <- colnames(sobj.combined)
vdj_labels[c(TRDV2_BC1, TRDV2_BC2)] <- "TRDV2+"
df <- data.frame(sobj.combined@reductions$umap@cell.embeddings)
df$label <- vdj_labels[rownames(df)]
table(df$label)
dropouts <- rownames(df)[which(df$label == "nothing")]
tiff("Figures/Figure3/Fig3C_TRDV2.tiff")
ggplot(
  df[dropouts, ], 
  aes(x=UMAP_1, y=UMAP_2, color=factor(label), fill=factor(label))
) + 
  geom_point(size=1) +
  theme_void() + 
  theme(legend.position="none") + 
  scale_color_manual(values=c("lightgrey", "#0073C2FF")) +
  scale_fill_manual(values=c("lightgrey", "#0073C2FF")) + 
  geom_point(
    data = df[c(TRDV2_BC1, TRDV2_BC2), ], 
    colour = "#0073C2FF", 
    size = 1
  )
dev.off()

# TRGV9 - Figure 3C(right)
vdj_labels <- rep("nothing", ncol(sobj.combined))
names(vdj_labels) <- colnames(sobj.combined)
vdj_labels[c(TRGV9_BC1, TRGV9_BC2)] <- "TRGV9+"
df <- data.frame(sobj.combined@reductions$umap@cell.embeddings)
df$label <- vdj_labels[rownames(df)]
table(df$label)
dropouts <- rownames(df)[which(df$label == "nothing")]
tiff("Figures/Figure3/Fig3C_TRGV9.tiff")
ggplot(
  df[dropouts, ], 
  aes(x=UMAP_1, y=UMAP_2, color=factor(label), fill=factor(label))
) + 
  geom_point(size=1) +
  theme_void() + 
  theme(legend.position="none") + 
  scale_color_manual(values=c("lightgrey", "#A73030FF")) +
  scale_fill_manual(values=c("lightgrey", "#A73030FF")) + 
  geom_point(
    data = df[c(TRGV9_BC1, TRGV9_BC2), ], 
    colour = "#A73030FF", 
    size = 1
  )
dev.off()

# -------------------------------------------------------
# Figure 3D
# -------------------------------------------------------
TRDC_cells <- c(TRDC_BC1, TRDC_BC2)
TRDV2_cells <- c(TRDV2_BC1, TRDV2_BC2)
TRGV9_cells <- c(TRGV9_BC1, TRGV9_BC2)

stats <- data.frame(
  matrix(NA, nrow = 3, ncol = length(unique(Idents(sobj.combined)))), 
  row.names = c("TRDC", "TRDV2", "TRGV9")
)
colnames(stats) <- unique(Idents(sobj.combined))
for (ID in colnames(stats)) {
  cells = WhichCells(sobj.combined, idents = ID)
  stats["TRDC", ID] <- length(intersect(cells, TRDC_cells))*100/length(cells)
  stats["TRDV2", ID] <- length(intersect(cells, TRDV2_cells))*100/length(cells)
  stats["TRGV9", ID] <- length(intersect(cells, TRGV9_cells))*100/length(cells)
}
df <- reshape2::melt(t(stats[, c("4", "6", "9")]))
colnames(df) <- c("ClusterID", "gene", "pct")
df$ClusterID <- factor(df$ClusterID, levels = c("9", "6", "4"))
pdf("Figures/Figure3/Fig3D.pdf", width = 4, height = 2)
ggplot(df, aes(x=ClusterID, y=pct, fill=factor(gene))) + geom_bar(stat = "identity", position=position_dodge()) +
  labs(title = "", y="", x="") + theme_minimal() + 
  theme(
    axis.text = element_text(size=7), 
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")
  )
dev.off()

# Print counts for plot
counts <- data.frame(
  matrix(NA, nrow = 3, ncol = length(unique(Idents(sobj.combined)))), 
  row.names = c("TRDC", "TRDV2", "TRGV9")
)
colnames(counts) <- unique(Idents(sobj.combined))
for (ID in colnames(counts)) {
  cells <- WhichCells(sobj.combined, idents = ID)
  counts["TRDC", ID] <- length(intersect(cells, TRDC_cells))
  counts["TRDV2", ID] <- length(intersect(cells, TRDV2_cells))
  counts["TRGV9", ID] <- length(intersect(cells, TRGV9_cells))
}
counts[, c("4", "6", "9")]


# -------------------------------------------------------
# Figure 3E
# -------------------------------------------------------
cells <- WhichCells(sobj.combined, idents = c(4, 6, 9))
p <- FeaturePlot(
  sobj.combined, 
  features = "KLRK1", 
  pt.size = 0.01, 
  combine = FALSE, 
  min.cutoff = 0, 
  max.cutoff = 4, 
  cells = cells
)
pdf("Figures/Figure3/FeaturePlots/KLRK1.pdf", width = 2, height = 2)
plot(p[[1]] + NoAxes() + NoLegend())
dev.off()

genes <- c("IFNG", "CCR6")
p <- FeaturePlot(
  sobj.combined, 
  features = genes, 
  pt.size = 0.01, 
  combine = FALSE, 
  min.cutoff = 0, 
  max.cutoff = 3, 
  cells = cells
)
for(i in 1:length(p)) {
  pdf(paste("Figures/Figure3/FeaturePlots/", genes[i], ".pdf", sep = ""), width = 2, height = 2)
  plot(p[[i]] + NoAxes() + NoLegend())
  dev.off()
}



