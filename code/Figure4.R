library(Seurat)

sobj.combined_pbmc <- readRDS("gdT_paper_analysis/data/processed/PBMC_sobj_combined.rds")

genes <- c("CCR6","NCAM1","FCGR3A","SOX4")

p <- FeaturePlot(
  sobj.combined_pbmc, 
  features = genes, 
  pt.size = 0.01,
  combine = FALSE, 
  min.cutoff = 0, 
  max.cutoff = 3
)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + NoLegend()
  pdf(paste0("gdT_paper_analysis/Figures/Figure4/PBMC_features/", genes[i], ".pdf"), width = 1.5, height = 1.6)
  plot(p[[i]])
  dev.off()
}

sobj.combined_bc <- readRDS("gdT_paper_analysis/data/processed/BRCA_sobj_combined.rds")
sobj <- subset(
  x = sobj.combined_bc,
  cells = WhichCells(object = sobj.combined_bc, idents = c(4,6,9))
)

p <- FeaturePlot(
  sobj, 
  features = genes, 
  pt.size = 0.01, 
  combine = FALSE, 
  min.cutoff = 0,
  max.cutoff = 3
)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + NoLegend()
  pdf(paste0("gdT_paper_analysis/Figures/Figure4/BC_features/", genes[i], ".pdf"), width = 1.5, height = 1.6)
  plot(p[[i]])
  dev.off()
}
