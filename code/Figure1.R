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
# Read gene expression data
# -------------------------------------------------------
# P03E06
gem1 <- get_scrnaseq_data("p03e06/gdt_pbmc") %>%
  GetAssayData(slot = "counts") %>%
  as.data.frame() %>%
  set_colnames(paste("P03E06", colnames(.), sep = "_"))
filtered_cells <- read.delim("data/raw/HD6_cells_filtered.txt", stringsAsFactors = F, header = F)$V1
gem1 <- gem1[, filtered_cells]

# P03E04_P03E05
gem2 <- get_scrnaseq_data("p03e04_plus_p03e05/gdt_pbmc") %>%
  GetAssayData(slot = "counts") %>%
  as.data.frame() %>%
  set_colnames(paste("P03E04E05", colnames(.), sep = "_"))
filtered_cells <- read.delim("data/raw/HD45_cells_filtered.txt", stringsAsFactors = F, header = F)$V1
gem2 <- gem2[, filtered_cells]

# -------------------------------------------------------
# Set up objects
# -------------------------------------------------------
sobj1 <- CreateSeuratObject(counts = gem1, project = "P03E06", min.cells = 5)
sobj1$stim <- "P03E06"
sobj1 <- subset(sobj1, subset = nFeature_RNA > 500)
sobj1 <- NormalizeData(sobj1, verbose = FALSE)
sobj1 <- FindVariableFeatures(sobj1, selection.method = "vst", nfeatures = 2000)

sobj2 <- CreateSeuratObject(counts = gem2, project = "P03E04E05", min.cells = 5)
sobj2$stim <- "P03E04E05"
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
# Clustering and UMAP 
sobj.combined <- RunUMAP(sobj.combined, reduction = "pca", dims = 1:20)
sobj.combined <- FindNeighbors(sobj.combined, reduction = "pca", dims = 1:20)
sobj.combined <- FindClusters(sobj.combined, resolution = 0.5)

p1 <- DimPlot(sobj.combined, label = T) + NoLegend() + NoAxes()

# -------------------------------------------------------
# Check expression of TRDC for cluster filtering
# -------------------------------------------------------
vdj_HD6 <- read.delim(
  "data/raw/VDJ_data/HD6_vdj_list.txt", 
  stringsAsFactors = F, 
  header = F
) %>%
  mutate(V1 = unlist(lapply(V1, function(x) strsplit(x, "_")[[1]][3]))) %>%
  mutate(V1 = paste("P03E06", V1, sep = "_")) %>%
  filter(V1 %in% colnames(gem1))
TRDC_E06 <- vdj_HD6$V1[grep("TRDC", vdj_HD6$V2)]

vdj_HD45 <- read.delim(
  "data/raw/VDJ_data/HD45_vdj_list.txt", 
  stringsAsFactors = F,
  header = F
) %>%
  mutate(V1 = lapply(V1, function(x) strsplit(x, "_")[[1]][3])) %>%
  mutate(V1 = paste("P03E04E05", V1, sep = "_")) %>%
  filter(V1 %in% colnames(gem2))
TRDC_E045 <- vdj_HD45$V1[grep("TRDC", vdj_HD45$V2)]

TRDC_cells <- c(TRDC_E045, TRDC_E06)
vdj_labels <- rep("nothing", ncol(sobj.combined))
names(vdj_labels) <- colnames(sobj.combined)
vdj_labels[TRDC_cells] <- "TRDC"
df <- data.frame(sobj.combined@reductions$umap@cell.embeddings)
df$label <- vdj_labels[rownames(df)]

table(df$label)
dropouts <- rownames(df)[which(df$label == "nothing")]
p2 <- ggplot(
  df[dropouts, ], 
  aes(x=UMAP_1, y=UMAP_2, color=factor(label), 
      fill=factor(label))
  ) + 
  geom_point(size=1) +
  theme_void() + 
  theme(legend.position="none") + 
  scale_color_manual(values=c("lightgrey","#0073C2FF")) +
  scale_fill_manual(values=c("lightgrey","#0073C2FF")) + 
  geom_point(data=df[TRDC_cells, ],colour="#A73030FF",size=1) 

grid.arrange(p1, p2, ncol=2)

coord <- data.frame(sobj.combined@reductions$umap@cell.embeddings)
plot(coord$UMAP_1, coord$UMAP_2, pch = 16)
abline(h = -11)
sobj.combined <- subset(
  sobj.combined, 
  cells = rownames(coord)[which(coord$UMAP_2 >= -11)]
)

# -------------------------------------------------------
# Rename clusters
# -------------------------------------------------------
new.cluster.ids <- c("c.gd3","c.gd4","c.gd1","c.gd5","c.gd2")
names(new.cluster.ids) <- levels(sobj.combined)
sobj.combined <- RenameIdents(sobj.combined, new.cluster.ids)
Idents(sobj.combined) <- factor(
  Idents(sobj.combined), 
  levels = c("c.gd1", "c.gd2", "c.gd3","c.gd4", "c.gd5")
)

# -------------------------------------------------------
# Visualization 
# -------------------------------------------------------
# Figure 1A
pdf("Figures/Figure1/Fig1A.pdf", width = 2, height = 2)
DimPlot(sobj.combined, reduction = "umap", label = TRUE, pt.size = 0.01) + NoAxes() + NoLegend()
dev.off()
# Figure 1B
pdf("Figures/Figure1/Fig1B.pdf", width = 4, height = 2.2)
DimPlot(sobj.combined, reduction = "umap", group.by = "stim", pt.size = 0.01, split.by = 'stim') + NoAxes() + NoLegend()
dev.off()

# -------------------------------------------------------
# Delta1/Delta2 genetics
# -------------------------------------------------------
# P03E06
TRDV2_E06 <- vdj_HD6$V1[grep("TRDV2", vdj_HD6$V2)]
TRDC_E06 <- vdj_HD6$V1[grep("TRDC", vdj_HD6$V2)]
TRGV9_E06 <- vdj_HD6$V1[grep("TRGV9", vdj_HD6$V2)]

# P03E04_E05
TRDV2_E045 <- vdj_HD45$V1[grep("TRDV2", vdj_HD45$V2)]
TRDC_E045 <- vdj_HD45$V1[grep("TRDC", vdj_HD45$V2)]
TRGV9_E045 <- vdj_HD45$V1[grep("TRGV9", vdj_HD45$V2)]

# TRDC - Figure 1C(left)
vdj_labels <- rep("TRDC-", ncol(sobj.combined))
names(vdj_labels) <- colnames(sobj.combined)
vdj_labels[c(TRDC_E045, TRDC_E06)] <- "TRDC+"
df <- data.frame(sobj.combined@reductions$umap@cell.embeddings)
df$label <- vdj_labels[rownames(df)]
table(df$label)
dropouts <- rownames(df)[which(df$label == "TRDC-")]
tiff("Figures/Figure1/Fig1C_TRDC.tiff")
ggplot(
  df[dropouts, ], 
  aes(x=UMAP_1, y=UMAP_2, color=factor(label), fill=factor(label))
) + 
  geom_point(size=1) +
  theme_void() + 
  theme(legend.position="none") + 
  scale_color_manual(values=c("lightgrey","black")) +
  scale_fill_manual(values=c("lightgrey","black")) + 
  geom_point(data=df[c(TRDC_E045, TRDC_E06), ], colour="black", size=1)
dev.off()

# TRDV2 - Figure 1C(middle)
vdj_labels <- rep("nothing", ncol(sobj.combined))
names(vdj_labels) <- colnames(sobj.combined)
vdj_labels[c(TRDV2_E045, TRDV2_E06)] <- "TRDV2+"
df <- data.frame(sobj.combined@reductions$umap@cell.embeddings)
df$label <- vdj_labels[rownames(df)]
table(df$label)
dropouts <- rownames(df)[which(df$label == "nothing")]
tiff("Figures/Figure1/Fig1C_TRDV2.tiff")
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
    data = df[c(TRDV2_E045, TRDV2_E06), ], 
    colour = "#0073C2FF", 
    size = 1
  )
dev.off()

# TRGV9 - Figure 1C(right)
vdj_labels <- rep("nothing", ncol(sobj.combined))
names(vdj_labels) <- colnames(sobj.combined)
vdj_labels[c(TRGV9_E045, TRGV9_E06)] <- "TRGV9+"
df <- data.frame(sobj.combined@reductions$umap@cell.embeddings)
df$label <- vdj_labels[rownames(df)]
table(df$label)
dropouts <- rownames(df)[which(df$label == "nothing")]
tiff("Figures/Figure1/Fig1C_TRGV9.tiff")
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
    data = df[c(TRGV9_E045, TRGV9_E06), ], 
    colour = "#A73030FF", 
    size = 1
  )
dev.off()

# -------------------------------------------------------
# Figure 1D
# -------------------------------------------------------
TRDC_cells <- c(TRDC_E045, TRDC_E06)
TRDV2_cells <- c(TRDV2_E045, TRDV2_E06)
TRGV9_cells <- c(TRGV9_E045, TRGV9_E06)

gdT_gpc <- c(apply(gem1, 2, function(x) sum(x > 0)), 
            apply(gem2, 2, function(x) sum(x > 0)))

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
df <- reshape2::melt(t(stats))
colnames(df) <- c("ClusterID", "gene", "pct")
df$ClusterID <- factor(df$ClusterID, levels = c("c.gd1", "c.gd2", "c.gd3","c.gd4", "c.gd5"))
pdf("Figures/Figure1/Fig1D.pdf", width = 4, height = 2)
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
counts

# # -------------------------------------------------------
# # Chi-square test of independence for TCR delta genes
# # -------------------------------------------------------
# # prepare data for chi square test
# df <- t(counts) %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column("clusterID") %>%
#   mutate(macroCluster = ifelse(clusterID %in% c("c.gd1", "c.gd2"), 0, 1))
# 
# macro1 <- length(WhichCells(sobj.combined, idents = c("c.gd1", "c.gd2")))
# macro2 <- length(WhichCells(sobj.combined, idents = c("c.gd3", "c.gd4", "c.gd5")))
# 
# # chi-square test of independence for TRDC and TRGV9
# #res <- chisq.test(aggregate(df$TRDC, by=list(Category=df$macroCluster), FUN=sum))
# #res$expected
# #res$observed
# #res$p.value
# #res$statistic
# 
# observed <- c(
#   sum(df$TRDC[which(df$macroCluster == 0)]), 
#   sum(df$TRDC[which(df$macroCluster == 1)]) 
# )
# p <- sum(observed) / ncol(sobj.combined)
# expected <- c(p * macro1, p*macro2)
# chisq <- sum((observed - expected) ^ 2 / expected)
# pchisq(chisq, df=1)
# 
# 
# # chi-square test of independence for TRDV2 and macroCluster
# #chisq.test(aggregate(df$TRDV2, by=list(Category=df$macroCluster), FUN=sum))
# observed <- c(
#   sum(df$TRDV2[which(df$macroCluster == 0)]), 
#   sum(df$TRDV2[which(df$macroCluster == 1)]) 
# )
# p <- sum(observed) / ncol(sobj.combined)
# expected <- c(p * macro1, p*macro2)
# chisq <- sum((observed - expected) ^ 2 / expected)
# pchisq(chisq, df=1)
# 
# 
# # chi-square test of independence for TRGV9 and macroCluster
# #chisq.test(aggregate(df$TRGV9, by=list(Category=df$macroCluster), FUN=sum))
# observed <- c(
#   sum(df$TRGV9[which(df$macroCluster == 0)]), 
#   sum(df$TRGV9[which(df$macroCluster == 1)]) 
# )
# p <- sum(observed) / ncol(sobj.combined)
# expected <- c(p * macro1, p*macro2)
# chisq <- sum((observed - expected) ^ 2 / expected)
# pchisq(chisq, df=1)

# -------------------------------------------------------
# Show IL17 and IFNG pathway genes -Feature plots 1E
# -------------------------------------------------------
genes <- c("SOX4", "CCR7", "LEF1")
p <- FeaturePlot(
  sobj.combined, 
  features = genes, 
  pt.size = 0.01, 
  combine = FALSE, 
  min.cutoff = 0, 
  max.cutoff = 3
)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + NoLegend()
  pdf(paste("Figures/Figure1/FeaturePlots/", genes[i], ".pdf", sep = ""), width = 2, height = 2)
  plot(p[[i]])
  dev.off()
}

genes <- c("IFNG", "RORC", "CCR6")
p <- FeaturePlot(
  sobj.combined, 
  features = genes, 
  pt.size = 0.01, 
  combine = FALSE, 
  min.cutoff = 0, 
  max.cutoff = 0.125
)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + NoLegend()
  pdf(paste("Figures/Figure1/FeaturePlots/", genes[i], ".pdf", sep = ""), width = 2, height = 2)
  plot(p[[i]])
  dev.off()
}

# -------------------------------------------------------
# Supervised pathway enrichment with scID
# -------------------------------------------------------

AP_genes <- c("HLA.DQPB1", "HLA.DRA", "HLA.DPA1")
cytotoxic_genes <- c("GZMA", "GZMB", "GZMK", "GZMM", "GZMH", "PRF1", "TRAIL", "FAS", "IL12")
IFNG <- read.delim("data/curated_genesets/IFNG_hallmark_geneset.txt", stringsAsFactors = F)
IFNG <- unique(c("TBX21", "EOMES", "STAT1", "STAT4", "IL12RB", "IFNG", IFNG$GeneSymbol))
markers <- data.frame(
  gene = c(
    cytotoxic_genes,
    "RORC", "IL23R", "CCR6", "IL1R1", "RORA", "BLK", "IL17A", 
    IFNG, 
    AP_genes, 
    "KLRK1", "HCST"
  ), 
  cluster = c(
    rep("Cytotoxic", length(cytotoxic_genes)), 
    rep("Il17", 7), 
    rep("IFNG", length(IFNG)),
    rep("AntigenPresentation", length(AP_genes)), 
    rep("innate", 2)), 
  avg_logFC = 1
)
weights = list()
for (ct in unique(markers$cluster)) {
  g <- markers$gene[which(markers$cluster == ct)]
  s <- rep(1, length(g))
  names(s) <- g
  weights[[ct]] <- s
}

# HD6
HD6_labels <- Idents(sobj.combined)[grep("E06", names(Idents(sobj.combined)))]
# Keep only gdT cells
HD6_gem <- counts_to_cpm(gem1[, names(HD6_labels)])
scID_res <- scid_multiclass(
  target_gem = HD6_gem, 
  weights = weights, 
  markers = markers
)
table(scID_res$labels)
scores <- data.frame(t(scID_res$scores))

scores$labels <- HD6_labels[rownames(scores)]
p1 <- ggplot2.histogram(
  data = scores, 
  xName = 'Cytotoxic', 
  groupName = 'labels', 
  legendPosition = "top",
  alpha = 0.5, 
  addDensity = TRUE, 
  addMeanLine = TRUE, 
  meanLineColor = "white", 
  meanLineSize = 1.5
)

p2 <- ggplot2.histogram(
  data = scores, 
  xName = 'Il17', 
  groupName = 'labels', 
  legendPosition = "top",
  alpha = 0.5, 
  addDensity = TRUE, 
  addMeanLine = TRUE, 
  meanLineColor ="white", 
  meanLineSize = 1.5
)

p3 <- ggplot2.histogram(
  data = scores, 
  xName ='IFNG', 
  groupName = 'labels', 
  legendPosition = "top",
  alpha = 0.5, 
  addDensity = TRUE, 
  addMeanLine = TRUE, 
  meanLineColor = "white", 
  meanLineSize = 1.5
)

p4 <- ggplot2.histogram(
  data = scores, 
  xName = 'AntigenPresentation', 
  groupName = 'labels', 
  legendPosition = "top",
  alpha = 0.5, 
  addDensity = TRUE, 
  addMeanLine = TRUE, 
  meanLineColor = "white", 
  meanLineSize = 1.5
)

p5 <- ggplot2.histogram(
  data = scores, 
  xName = 'innate', 
  groupName = 'labels', 
  legendPosition = "top",
  alpha = 0.5, 
  addDensity = TRUE, 
  addMeanLine = TRUE, 
  meanLineColor = "white", 
  meanLineSize = 1.5
)

grid.arrange(p1, p2, p3, p4, p5, ncol=2)

scores_scaled <- data.frame(scale(scores[, 1:5], scale = TRUE))
scores_scaled$labels <- scores$labels
df.m <- reshape2::melt(scores_scaled)
pdf("Figures/Figure1/Fig1F.pdf", width = 3.5, height = 3)
ggplot(df.m, aes(x = variable, y = value, fill = labels)) + 
  geom_boxplot(notch = TRUE, width = 0.5, outlier.size = 0, lwd = 0.5) +
  labs(title = "", y = "", x = "") + 
  theme_minimal() + 
  theme(
    axis.text = element_text(size=7), 
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"), 
    legend.position = "top"
  )
dev.off()

# -------------------------------------------------------
# Find subtype markers
# -------------------------------------------------------
markers <- FindConservedMarkers(
  sobj.combined, 
  ident.1 = "c.gd1", 
  test.use="MAST", 
  only.pos=TRUE, 
  grouping.var = "stim", 
  logfc.threshold = 0.2
)
markers$gene <- rownames(markers)
markers$cluster <- "c.gd1"
for (ID in c("c.gd2", "c.gd3","c.gd4", "c.gd5")) {
  m <- FindConservedMarkers(
    sobj.combined, 
    ident.1 = ID, 
    test.use="MAST", 
    only.pos=TRUE, 
    grouping.var = "stim", 
    logfc.threshold = 0.2
  )
  m <- m %>%
    mutate(
      gene = rownames(m),
      cluster = ID
    )
  markers <- rbind(markers, m)
}
table(markers$cluster)
write.table(markers, file="data/processed/PBMC_subtype_markers.txt", quote = F, sep = "\t")

# Select top markers per cluster by p_val_adj
top60 <- markers[which(markers$gene %in% rownames(sobj.combined)), ] %>% 
    group_by(cluster) %>% 
    top_n(n = -60, wt = P03E06_p_val_adj)
table(top60$cluster)

# Heatmap top markers per cluster
tiff("Figures/Figure1/Fig1G.tiff")
DoHeatmap(sobj.combined, features = top60$gene, size = 3, raster = F) + FontSize(y.text = 5)
dev.off()

# -------------------------------------------------------
# Save sobj.combined 
# -------------------------------------------------------
saveRDS(sobj.combined, file="data/processed/PBMC_sobj_combined.rds")


