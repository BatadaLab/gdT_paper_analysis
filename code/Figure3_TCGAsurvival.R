# TCGA
library(Seurat)
library(dplyr)
library(tibble)

custom_ggsurvplot <- function(fit, data) {
  gg <- ggsurvplot(
    fit, 
    data = data, 
    risk.table = TRUE, 
    pval = TRUE, 
    conf.int = TRUE, 
    xlim = c(0, 2000), 
    break.time.by = 500, 
    ggtheme = theme_minimal(), 
    risk.table.y.text.col = TRUE, 
    risk.table.y.text = FALSE
  )
  p <- gg$plot + 
    theme(
      axis.text=element_text(size = 7), 
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black")
    )
  
  p
}

# ------------------------------------------------------------
# Get gdT subtype signatures from BC1 data
# No need to repeat this - load pre-extracted signatures
# ------------------------------------------------------------

# source("~/Google Drive/bin/batadalab_scrnaseq_utils.R")
# sobj <- get_scrnaseq_data("p01e01_custom")
# sobj.combined_bc <- readRDS("~/Google Drive/Data/Analysis/gdT_human_pbmc/data/BC_combined_22Nov.rds")
# BC1_labels <- Idents(sobj.combined_bc)[grep("BC1", names(Idents(sobj.combined_bc)))]
# names(BC1_labels) <- unlist(lapply(names(BC1_labels), function(x) strsplit(x, "_")[[1]][2]))
# # Keep only gdT cells
# cells <- intersect(colnames(sobj), names(BC1_labels))
# BC1_gem <- as.data.frame(GetAssayData(sobj))[, cells]
# 
# sobj_new <- sobj
# Idents(sobj_new) <- factor(BC1_labels[colnames(sobj_new)])
# 
# gdT1 <- FindMarkers(sobj_new, ident.1 = 4, test.use = "MAST", only.pos = T, logfc.threshold = 0.7)
# gdT2 <- FindMarkers(sobj_new, ident.1 = 6, test.use = "MAST", only.pos = T, logfc.threshold = 0.7)
# gdT3 <- FindMarkers(sobj_new, ident.1 = 9, test.use = "MAST", only.pos = T, logfc.threshold = 0.7)
# 
# write.table(gdT1, file = "../data/BC_gdT_signatures/gdT1.txt", quote = F, sep = "\t")
# write.table(gdT2, file = "../data/BC_gdT_signatures/gdT2.txt", quote = F, sep = "\t")
# write.table(gdT3, file = "../data/BC_gdT_signatures/gdT3.txt", quote = F, sep = "\t")

gdT1 <- read.delim("../data/BC_gdT_signatures/gdT1.txt", stringsAsFactors = FALSE)
gdT2 <- read.delim("../data/BC_gdT_signatures/gdT2.txt", stringsAsFactors = FALSE)
gdT3 <- read.delim("../data/BC_gdT_signatures/gdT3.txt", stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Read TCGA data
# ------------------------------------------------------------

tcga_data <- read.delim("~/Google Drive/Data/tcga/rnaseq/BRCA_rnaseq.tumour", stringsAsFactors = F, row.names = 1)
# Change sample names: Keep only first 12 characters
colnames(tcga_data) <- unlist(lapply(colnames(tcga_data), function(x) substr(x, 1, 12)))

gdT1 <- gdT1[which(rownames(gdT1) %in% rownames(tcga_data)), ]
gdT2 <- gdT2[which(rownames(gdT2) %in% rownames(tcga_data)), ]
gdT3 <- gdT3[which(rownames(gdT3) %in% rownames(tcga_data)), ]

gd1Score <- colSums(tcga_data[rownames(gdT1), ] * gdT1$avg_logFC) / sum(gdT1$avg_logFC)
gd2Score <- colSums(tcga_data[rownames(gdT2), ] * gdT2$avg_logFC) / sum(gdT2$avg_logFC)
gd3Score <- colSums(tcga_data[rownames(gdT3), ] * gdT3$avg_logFC) / sum(gdT3$avg_logFC)

# Correct for purity
purity_data <- read.delim("~/Google Drive/Data/tcga/purity/aran2015_tumourpurity_estimate_tcga.txt", stringsAsFactors = F)
purity_data$Sample.ID <- unlist(lapply(purity_data$Sample.ID, function(x) gsub("-", ".", x)))
purity_data$Sample.ID <- unlist(lapply(purity_data$Sample.ID, function(x) substr(x, 1, 12)))

# Keep only common samples
purity_data <- purity_data[which(purity_data$Sample.ID %in% colnames(tcga_data)), ]
# remove duplicated samples
keepSamples <- setdiff(purity_data$Sample.ID, purity_data$Sample.ID[which(purity_data$Sample.ID %in% purity_data$Sample.ID[which(duplicated(purity_data$Sample.ID))])])
purity_data <- purity_data[which(purity_data$Sample.ID %in% keepSamples), ]
rownames(purity_data) <- purity_data$Sample.ID
tcga_data <- tcga_data[, keepSamples]

hist(purity_data$ESTIMATE, breaks = 50)
abline(v=c(0.6, 0.7))

# Select samples by purity: between 0.6 and 0.7
samples <- purity_data$Sample.ID[intersect(which(purity_data$ESTIMATE >= 0.6), which(purity_data$ESTIMATE <= 0.7))]

s1 <- gd1Score[samples]
G1_top <- names(which(s1 >= quantile(s1, 2/3)))
G1_bottom <- names(which(s1 <= quantile(s1, 1/3)))

s2 <- gd2Score[samples]
G2_top <- names(which(s2 >= quantile(s2, 2/3)))
G2_bottom <- names(which(s2 <= quantile(s2, 1/3)))

s3 <- gd3Score[samples]
G3_top <- names(which(s3 >= quantile(s3, 2/3)))
G3_bottom <- names(which(s3 <= quantile(s3, 1/3)))

# -------------------------------------------------------------
# Survival Analysis
# -------------------------------------------------------------
library(survminer)
library(survival)
library(RTCGA.clinical)

BRCAOV.survInfo <- survivalTCGA(
  BRCA.clinical, 
  OV.clinical,extract.cols = "admin.disease_code"
) %>%
  mutate(bcr_patient_barcode = unlist(lapply(bcr_patient_barcode, function(x) gsub("-", ".", x))))

# G1
BRCAOV.survInfo_G1 <- BRCAOV.survInfo %>%
  filter(bcr_patient_barcode %in% c(G1_bottom, G1_top)) %>%
  #column_to_rownames("bcr_patient_barcode") %>%
  mutate(gdStatus = case_when(
    bcr_patient_barcode %in% G1_bottom ~ "low", 
    bcr_patient_barcode %in% G1_top ~ "high", 
  ))
fit_G1 <- survfit(Surv(times, patient.vital_status) ~ gdStatus,data = BRCAOV.survInfo_G1)
p1 <- custom_ggsurvplot(fit = fit_G1, data = BRCAOV.survInfo_G1)

# G2
BRCAOV.survInfo_G2 <- BRCAOV.survInfo %>%
  filter(bcr_patient_barcode %in% c(G2_bottom, G2_top)) %>%
  mutate(gdStatus = case_when(
    bcr_patient_barcode %in% G2_bottom ~ "low", 
    bcr_patient_barcode %in% G2_top ~ "high", 
  ))

fit_G2 <- survfit(Surv(times, patient.vital_status) ~ gdStatus,data = BRCAOV.survInfo_G2)
p2 <- custom_ggsurvplot(fit = fit_G2, data = BRCAOV.survInfo_G2)

# G3
BRCAOV.survInfo_G3 <- BRCAOV.survInfo %>%
  filter(bcr_patient_barcode %in% c(G3_bottom, G3_top)) %>%
  mutate(gdStatus = case_when(
    bcr_patient_barcode %in% G3_bottom ~ "low", 
    bcr_patient_barcode %in% G3_top ~ "high", 
  ))
fit_G3 <- survfit(Surv(times, patient.vital_status) ~ gdStatus,data = BRCAOV.survInfo_G3)
p3 <- custom_ggsurvplot(fit = fit_G3, data = BRCAOV.survInfo_G3)

gridExtra::grid.arrange(p1, p2, p3, ncol=1)
