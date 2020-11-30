library(pheatmap)
library(ggplot2)
library(ggpubr)

geometric.mean <- function(x, na.rm=TRUE) { 
  exp(mean(log(x), na.rm = na.rm)) 
}

# --------------------------------------------------------------------------
# Load data
# --------------------------------------------------------------------------
gem <- read.delim(
  "~/Google Drive/Data/tcga/rnaseq/BRCA_rnaseq.tumour", 
  stringsAsFactors = F, 
  row.names = 1
)
colnames(gem) <- unlist(lapply(colnames(gem), function(x) substr(x, 1, 12)))

samples_low <- read.delim(
  "~/Google Drive/Data/Analysis/gdT_human_pbmc/Manuscript/Figures/Figure3/TCGA_further/gdT2_low_samples.txt", 
  stringsAsFactors = F, 
  header = F
)$V1
samples_high <- read.delim(
  "~/Google Drive/Data/Analysis/gdT_human_pbmc/Manuscript/Figures/Figure3/TCGA_further/gdT2_high_samples.txt", 
  stringsAsFactors = F, 
  header = F
)$V1

metadata <- read.delim(
  "~/Google Drive/Data/tcga/brca_pid2subtype.tab", 
  stringsAsFactors = F
)
metadata$pid <- unlist(lapply(metadata$pid, function(x) gsub("-", ".", x)))

metadata <- metadata %>%
  filter(pid %in% c(samples_low, samples_high))
rownames(metadata) <- metadata$pid
metadata <- metadata[, c("mut", "er")]
metadata$label <- "high"
metadata[samples_low, "label"] <- "low"
metadata <- metadata[-which(is.na(metadata$er)), ]

# --------------------------------------------------------------------------
# Expression of markers
# --------------------------------------------------------------------------
df <- data.frame(t(gem[c("CD3D", "CD4", "CD8A"), c(samples_low, samples_high)]))
df$labels <- rep("high")
df[samples_low, "labels"] <- "low"

df.m <- reshape2::melt(df)
pdf("~/Google Drive/Data/Analysis/gdT_human_pbmc/Manuscript/Figures/Figure3/TCGA_further/Tcellmarkers.pdf", width = 3, height = 3)
ggplot(df.m, aes(x=variable, y=value, fill=labels)) + geom_boxplot(notch=TRUE, width = 0.5, outlier.size = 0, lwd = 0.5) +
  labs(title = "", y="", x="") + theme_minimal() +  stat_compare_means() +
  theme(axis.text=element_text(size=7), plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))
dev.off()

# --------------------------------------------------------------------------
# Cytolytic score
# --------------------------------------------------------------------------
samples <- c(samples_high, samples_low)
score <- apply(
  gem[c("PRF1", "GZMA"), samples], 
  2, function(x) geometric.mean(x)
)
cyto_score <- data.frame(
  score = score[samples], 
  sampleID = samples, 
  label = c(rep("high", length(samples_high)), rep("low", length(samples_low)))
)

pdf("~/gdT_paper_analysis/Figures/Supplementary/SF3_cytolytic_score.pdf", width = 3, height = 3)
ggplot(cyto_score, aes(x=label, y=score)) + 
    geom_boxplot(notch=TRUE, width = 0.5, outlier.size = 0, lwd = 0.5) +
    labs(title = "", y="", x="") + 
    theme_minimal() +  
    stat_compare_means() +
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
# --------------------------------------------------------------------------
# Types of cancer
# --------------------------------------------------------------------------

df = data.frame(ER = rep(unique(metadata$er), 2), 
                label = c(rep("high", length(unique(metadata$er))), 
                          rep("low", length(unique(metadata$er)))),
                count = NA)
for (i in 1:nrow(df)) {
  df[i, "count"] = length(intersect(which(metadata$er == df[i, "ER"]), 
                                    which(metadata$label == df[i, "label"])))
}
df$count[which(df$label == "high")] = df$count[which(df$label == "high")] *100 / length(samples_high)
df$count[which(df$label == "low")] = df$count[which(df$label == "low")] *100 / length(samples_low)

ggplot(df, aes(x=label, y=count, fill=ER)) + 
  geom_bar(stat="identity", position = "dodge") + 
  labs(title = "", y="", x="") + theme_minimal() + #geom_jitter(width = 0.1, size = 0.5)+
  theme(axis.text=element_text(size=7), plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"), legend.position = "top")


ggplot(df, aes(x=label, y=count, fill="ER")) + geom_bar(stat="identity")

ggplot(df, aes(fill=ER, y=count, x=label)) + 
  geom_bar(stat="identity") +
  labs(title = "", y="Er", x="") 

# --------------------------------------------------------------------------
# Mutational load
# --------------------------------------------------------------------------
pdf("~/Google Drive/Data/Analysis/gdT_human_pbmc/Manuscript/Figures/Figure3/TCGA_further/mutation_load.pdf", width = 3, height = 3)
ggplot(metadata, aes(x=label, y=mut)) + geom_boxplot(notch=TRUE, width = 0.5, outlier.size = 0, lwd = 0.5) +
  labs(title = "", y="", x="") + theme_minimal() +  stat_compare_means() +
  theme(axis.text=element_text(size=7), plot.margin = unit(c(0,0,0,0), "cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"))
dev.off()


