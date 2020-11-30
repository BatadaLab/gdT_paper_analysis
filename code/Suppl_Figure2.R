library(Seurat)
library(scID)
library(dplyr)
library(ggplot2)
library(ggpubr)
source("~/Google Drive/bin/batadalab_scrnaseq_utils.R")

sobj <- get_scrnaseq_data("p01e01_custom")
sobj.combined_bc <- readRDS("~/Google Drive/Data/Analysis/gdT_human_pbmc/data/BC_combined_22Nov.rds")
BC1_labels <- Idents(sobj.combined_bc)[grep("BC1", names(Idents(sobj.combined_bc)))]
names(BC1_labels) <- unlist(lapply(names(BC1_labels), function(x) strsplit(x, "_")[[1]][2]))
# Keep only gdT cells
cells <- intersect(colnames(sobj), names(BC1_labels))
BC1_gem <- as.data.frame(GetAssayData(sobj))[, cells]

sobj_new <- subset(sobj, cells = names(BC1_labels))
Idents(sobj_new) <- factor(BC1_labels[colnames(sobj_new)])
DimPlot(sobj_new, label = T)

gdT1 <- FindMarkers(
  sobj_new, 
  ident.1 = 9, 
  test.use = "MAST", 
  only.pos = T, 
  logfc.threshold = 0.7
) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(cluster = "gdT.1")

gdT2 <- FindMarkers(
  sobj_new, 
  ident.1 = 6, 
  test.use = "MAST", 
  only.pos = T, 
  logfc.threshold = 0.7
) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(cluster = "gdT.2")

gdT3 <- FindMarkers(
  sobj_new, 
  ident.1 = 4, 
  test.use = "MAST", 
  only.pos = T, 
  logfc.threshold = 0.7
) %>%
  tibble::rownames_to_column("gene") %>%
  mutate(cluster = "gdT.3")

markers <- rbind(gdT1, gdT2, gdT3)

target_gem <- as.data.frame(as.matrix(GetAssayData(sobj_new)))
scID_res <- scid_multiclass(
  target_gem = target_gem, 
  markers = markers, 
  estimate_weights_from_target = TRUE
)

scores <- t(scID_res$scores)
df <- reshape2::melt(scores)
colnames(df) <- c("cellID", "signature", "score")
df$clusterID <- Idents(sobj_new)[df$cellID]
df$cluster <- "CD4T"
df$cluster[which(df$clusterID == 0)] <- "CD4T"
df$cluster[which(df$clusterID == 1)] <- "CD8T"
df$cluster[which(df$clusterID == 2)] <- "CD8T"
df$cluster[which(df$clusterID == 3)] <- "B"
df$cluster[which(df$clusterID == 4)] <- "gdT.3"
df$cluster[which(df$clusterID == 5)] <- "Treg"
df$cluster[which(df$clusterID == 6)] <- "gdT.2"
df$cluster[which(df$clusterID == 7)] <- "CD4T"
df$cluster[which(df$clusterID == 8)] <- "Mph"
df$cluster[which(df$clusterID == 9)] <- "gdT.1"
df$cluster[which(df$clusterID == 10)] <- "unclassified"
df$cluster[which(df$clusterID == 11)] <- "unclassified"
df$cluster[which(df$clusterID == 12)] <- "unclassified"


df1 <- df %>%
  filter(!is.na(clusterID)) %>%
  filter(signature == "gdT.1")

my_comparisons <- list( 
  c("gdT.1", "B"), 
  c("gdT.1", "CD4T"), 
  c("gdT.1", "CD8T"), 
  c("gdT.1", "gdT.2"), 
  c("gdT.1", "gdT.3"),
  c("gdT.1", "Mph"),
  c("gdT.1", "Treg"), 
  c("gdT.1", "unclassified") 
)
pdf("gdT_paper_analysis/Figures/Supplementary/SF2_gdT1.pdf")
ggplot(df1, aes(x=cluster, y=score, color=cluster)) + 
  geom_boxplot(notch=TRUE, width = 0.5, outlier.size = 0, lwd = 0.5) +
  labs(title = "", y="scID score", x="signature") + 
  theme_minimal() + 
  stat_compare_means(comparisons = my_comparisons) +
  theme(
    axis.text=element_text(size=10), 
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")
  )
dev.off()


df2 <- df %>%
  filter(!is.na(clusterID)) %>%
  filter(signature == "gdT.2")

my_comparisons <- list( 
  c("gdT.2", "B"), 
  c("gdT.2", "CD4T"), 
  c("gdT.2", "CD8T"), 
  c("gdT.2", "gdT.1"), 
  c("gdT.2", "gdT.3"),
  c("gdT.2", "Mph"),
  c("gdT.2", "Treg"), 
  c("gdT.2", "unclassified") 
)
pdf("gdT_paper_analysis/Figures/Supplementary/SF2_gdT2.pdf")
ggplot(df2, aes(x=cluster, y=score, color=cluster)) + 
  geom_boxplot(notch=TRUE, width = 0.5, outlier.size = 0, lwd = 0.5) +
  labs(title = "", y="scID score", x="signature") + 
  theme_minimal() + 
  stat_compare_means(comparisons = my_comparisons) +
  theme(
    axis.text=element_text(size=10), 
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")
  )
dev.off()


df3 <- df %>%
  filter(!is.na(clusterID)) %>%
  filter(signature == "gdT.3")

my_comparisons <- list( 
  c("gdT.3", "B"), 
  c("gdT.3", "CD4T"), 
  c("gdT.3", "CD8T"), 
  c("gdT.3", "gdT.1"), 
  c("gdT.3", "gdT.2"),
  c("gdT.3", "Mph"),
  c("gdT.3", "Treg"), 
  c("gdT.3", "unclassified") 
)
pdf("gdT_paper_analysis/Figures/Supplementary/SF2_gdT3.pdf")
ggplot(df3, aes(x=cluster, y=score, color=cluster)) + 
  geom_boxplot(notch=TRUE, width = 0.5, outlier.size = 0, lwd = 0.5) +
  labs(title = "", y="scID score", x="signature") + 
  theme_minimal() + 
  stat_compare_means(comparisons = my_comparisons) +
  theme(
    axis.text=element_text(size=10), 
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black")
  )
dev.off()



