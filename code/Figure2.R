library(Seurat)
library(dplyr)
library(ggplot2)
library(forcats)
library(scID)
library(ggpubr)
source("~/Google Drive/bin/batadalab_scrnaseq_utils.R")

setwd("~/gdT_paper_analysis/")
# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------
sobj.combined = readRDS("data/processed/PBMC_sobj_combined.rds")

# -------------------------------------------------------
# Figure 2A,C: Show expression of known and new subtype markers
# -------------------------------------------------------
# FCGR3A is CD16
genes <- c("CD28", "FCGR3A", "GPR56", "CXCR6")
p <- FeaturePlot(sobj.combined, features = genes, pt.size = 0.01, combine = FALSE, min.cutoff = 0, max.cutoff = 3)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + NoLegend()
  pdf(paste("Figures/Figure2/FeaturePlots/", genes[i], ".pdf", sep = ""), width = 2, height = 2)
  plot(p[[i]])
  dev.off()
}

# -------------------------------------------------------
# Figure 2B: calculate scID scores of d2 clusters' cells
# for CD16+ and CD28+ gene sets from Ryan et al., PNAS 2016
# -------------------------------------------------------
labels <- Idents(sobj.combined)
HD6_labels <- labels[grep("E06", names(labels))]
HD45_labels <- labels[grep("E04", names(labels))]

# P03E06
sobj_e06 <- get_scrnaseq_data("p03e06/gdt_pbmc")
gem_e06 <- data.frame(GetAssayData(sobj_e06, slot = "counts"))
colnames(gem_e06) <- paste("P03E06", colnames(gem_e06), sep = "_")
gem_cpm_e06 <- counts_to_cpm(gem_e06)

# P03E04_P03E05
sobj_e0405 <- get_scrnaseq_data("p03e04_plus_p03e05/gdt_pbmc")
gem_e0405 <- data.frame(GetAssayData(sobj_e0405, slot = "counts"))
colnames(gem_e0405) <- paste("P03E04E05", colnames(gem_e0405), sep = "_")
gem_cpm_e0405 <- counts_to_cpm(gem_e0405)

CD28pos <- c("GZMK", "CD27", "LTB", "CCR7", "MYC", "CCR6", "CD160", "CD7", "IL12RB2", "CD28", "CXCR6", "RORC", "SIPA1", "IL7R", "IL23R", "IRF1", "IL18R1", "CCR2")
CD16pos <- c("CX3CR1", "GZMB", "KIR3DL1", "KIR2DL4", "GNLY", "KIR2DL3", "KIR2DL1", "KIR3DL2", "GZMH", "KIR2DS5", "LAIR2", "KIR3DL3", "ITGB1", "PRF1", "SMAD7", "BCL9L", "LY6E", "IL18RB", "FCAR", "KIR2DL5A", "ITGAM", "CCL4L2")

markers <- data.frame(
  gene = c(CD28pos, CD16pos), 
  cluster = c(
    rep("CD28+", length(CD28pos)), 
    rep("CD16+", length(CD16pos))
  ), 
  avg_logFC = 1
)

res_HD6 <- scid_multiclass(
  target_gem = gem_cpm_e06[, names(HD6_labels)], 
  markers = markers, 
  estimate_weights_from_target = T
)
df_HD6 <- data.frame(t(res_HD6$score[names(HD6_labels)]))
df_HD6$label <- HD6_labels
df.m <- reshape2::melt(df_HD6)
p1 <- ggplot(df.m, aes(x = label, y = value, color = variable)) + 
  geom_boxplot(notch = TRUE, width = 0.5, outlier.size = 0, lwd = 0.5) +
  labs(title = "", y = "", x = "") + 
  theme_minimal() + 
  stat_compare_means(method = "wilcox.test", paired = TRUE) +
  theme(
    axis.text = element_text(size = 7), 
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line.x = element_line(colour = "black"), 
    axis.line.y = element_line(colour = "black")
  )

res_HD45 <- scid_multiclass(
  target_gem = gem_cpm_e0405[, names(HD45_labels)], 
  markers = markers, 
  estimate_weights_from_target = T
)
df_HD45 <- data.frame(t(res_HD45$score[names(HD45_labels)]))
df_HD45$label <- HD45_labels
df.m <- reshape2::melt(df_HD45)
p2 <- ggplot(df.m, aes(x = label, y = value, color = variable)) + 
  geom_boxplot(notch = TRUE, width = 0.5, outlier.size = 0, lwd = 0.5) +
  labs(title = "", y = "", x = "") + 
  theme_minimal() + 
  stat_compare_means(method = "wilcox.test", paired = TRUE) +
  theme(
    axis.text = element_text(size = 7), 
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line.x = element_line(colour = "black"), 
    axis.line.y = element_line(colour = "black")
  )

pdf("Figures/Figure2/Fig2B.pdf", width = 4, height = 4)
gridExtra::grid.arrange(p1, p2, ncol=1)
dev.off()
