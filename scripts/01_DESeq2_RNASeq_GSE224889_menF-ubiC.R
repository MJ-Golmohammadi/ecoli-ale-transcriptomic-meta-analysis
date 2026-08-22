# Clear the workspace
rm(list = ls())

library(DESeq2)
library(readxl)
library(dplyr)
library(readr)

#################################################
# Read raw counts sheet
#################################################

counts <- read_csv(
  "data/GSE224889/GSE224889_quinone_counts.csv.gz"
)

colnames(counts)
############################################################
# Rename gene ID column
############################################################

colnames(counts)[1] <- "Gene ID"

############################################################
# Read annotation
############################################################

library(rtracklayer)

gff <- import("data/Ecoli_K12_MG1655_annotation.gff")

annot <- data.frame(
  locus_tag = gff$locus_tag,
  gene = gff$gene
)

annot_map <- annot %>%
  filter(!is.na(locus_tag)) %>%
  distinct(locus_tag, gene)


############################################################
# Replace Gene IDs with Gene Names when available
############################################################

counts_new <- counts %>%
  left_join(
    annot_map,
    by = c(
      "Gene ID" = "locus_tag"
    )
  ) %>%
  mutate(
    `Gene ID` = ifelse(
      !is.na(gene) &
      gene != "",
      gene,
      `Gene ID`
    )
  ) %>%
  select(-gene)


############################################################
# Summary
############################################################

cat(
  "Mapped genes:",
  sum(!is.na(
    left_join(
      counts,
      annot_map,
      by = c("Gene ID" = "locus_tag")
    )$gene
  )),
  "\n"
)

############################################################
# Save updated file
############################################################

write.csv(
  counts_new,
  "data/GSE224889/GSE224889_count_data_GeneSymbol.csv",
  row.names = FALSE
)


################################################
# convert counts to counts_new
################################################
counts <- counts_new
sum(duplicated(counts$`Gene ID`))

############################################################
# Collapse duplicated gene symbols
############################################################

library(dplyr)

counts <- counts_new

counts <- counts %>%
  group_by(`Gene ID`) %>%
  summarise(
    across(
      where(is.numeric),
      sum,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

############################################################
# Set gene IDs as row names
############################################################

rownames(counts) <- counts$`Gene ID`

############################################################
# Remove Gene ID column for DESeq2
############################################################

counts_matrix <- counts %>%
  select(-`Gene ID`)

############################################################
# Quality control
############################################################

cat(
  "Duplicated genes:",
  sum(duplicated(rownames(counts_matrix))),
  "\n"
)

#################################################
# Define Samples
#################################################

# to solve the problem of converting genes to numbers
counts <- as.data.frame(counts)


counts_sub <- counts[,c(
"Gene ID",
"menF/ubiC KO, rep 1",
"menF/ubiC KO, rep 2",
"menF/ubiC KO adapted 1, rep 1",
"menF/ubiC KO adapted 1, rep 2",
"menF/ubiC KO adapted 2, rep 1",
"menF/ubiC KO adapted 2, rep 2",
"menF/ubiC KO adapted 3, rep 1",
"menF/ubiC KO adapted 3, rep 2",
"menF/ubiC KO adapted 4, rep 1",
"menF/ubiC KO adapted 4, rep 2"
)]

#################################################
# Convert to matrix
#################################################

counts_sub <- as.data.frame(counts_sub)

rownames(counts_sub) <- counts_sub$`Gene ID`
counts_sub$`Gene ID` <- NULL

count_matrix <- as.matrix(counts_sub)
rownames(count_matrix) <- rownames(counts_sub)
storage.mode(count_matrix) <- "integer"

summary(count_matrix)
#################################################
# Filter low expression genes
#################################################
cat(
  "Genes before filtering:",
  nrow(count_matrix),
  "\n"
)

keep <- rowSums(count_matrix >= 10) >= 3
count_matrix <- count_matrix[keep, ]

cat(
  "Genes after filtering:",
  nrow(count_matrix),
  "\n"
)


#################################################
# Metadata
#################################################

condition <- factor(c(
rep("KO",2),
rep("ALE",8)
))

coldata <- data.frame(
  row.names = colnames(count_matrix),
  condition
)

coldata$condition <- factor(
    coldata$condition,
    levels = c("KO","ALE")
)

levels(coldata$condition)

#################################################
# DESeq2
#################################################
condition <- factor(condition, levels = c("KO", "ALE"))

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)

resultsNames(dds)
#################################################
# VST Boxplot - Publication Quality
#################################################

library(reshape2)
library(ggplot2)

#################################################
# Extract VST matrix
#################################################

vst_mat <- assay(vsd)

#################################################
# Convert to data.frame
#################################################

vst_df <- as.data.frame(vst_mat)

vst_df$Gene <- rownames(vst_df)

#################################################
# Long format
#################################################

vst_long <- melt(
  vst_df,
  id.vars = "Gene",
  variable.name = "Sample",
  value.name = "Expression"
)

#################################################
# Sample annotation
#################################################

vst_long$Condition <- ifelse(
  grepl("^menF/ubiC KO, rep", vst_long$Sample),
  "KO",
  "ALE"
)

vst_long$Condition <- factor(
  vst_long$Condition,
  levels = c(
    "KO",
    "ALE"
  )
)

#################################################
# Boxplot
#################################################

p_box <- ggplot(
  vst_long,
  aes(
    x = Sample,
    y = Expression,
    fill = Condition
  )
) +

  geom_boxplot(
    linewidth = 0.4,
    outlier.size = 0.5
  ) +

  scale_fill_manual(
    values = c(
      "KO" = "#56B4E9",
      "ALE" = "#E69F00"
    )
  ) +

  labs(
    title = "VST-normalized expression distribution-menF-ubiC(GSE224889)",
    subtitle = "DESeq2 variance stabilizing transformation",
    x = NULL,
    y = "VST expression"
  ) +

  theme_bw(base_size = 13) +

  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),

    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10
    ),

    axis.text.y = element_text(
      size = 10
    ),

    legend.position = "top",

    legend.title = element_blank(),

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      linewidth = 0.8
    )
  )

#################################################
# Show plot
#################################################

p_box

#################################################
# Save figure
#################################################

ggsave(
  "figures/GSE224889/VST_boxplot_GSE224889-menF-ubiC.png",
  plot = p_box,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  "figures/GSE224889/VST_boxplot_GSE224889-menF-ubiC.pdf",
  plot = p_box,
  width = 7,
  height = 5
)

#################################################
# PCA
#################################################

# create new folder for results
dir.create("results/GSE224889", showWarnings = FALSE)

# create new folder for figures
dir.create("figures/GSE224889", showWarnings = FALSE)

# Save the normalized expression data to a CSV file
write.csv(assay(vsd), "results/GSE224889/VST_expression-menF-ubiC.csv")

p <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA of VST-normalized counts-menF-ubiC (GSE224889)") +
  theme(
    plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

p

#save PCA plot to figures folder
ggsave(
  "figures/GSE224889/PCA_plot_GSE224889-menF-ubiC.png",
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/GSE224889/PCA_plot_GSE224889-menF-ubiC.pdf",
  width = 8,
  height = 6
)

#################################################
## Sample distance heatmap
#################################################

library(pheatmap)

# Sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Sample annotation
annotation_col <- data.frame(
  Condition = colData(vsd)$condition
)

rownames(annotation_col) <- colnames(vsd)

#################################################
# Save PDF
#################################################

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col = annotation_col,
  annotation_names_col = TRUE,
  main = "Sample-to-Sample Distance after VST normalization - menF-ubiC (GSE224889)",
  filename = "figures/GSE224889/Sample_distance_heatmap_GSE224889-menF-ubiC.pdf",
  width = 6,
  height = 5,
  fontsize = 7,
  fontsize_row = 6,
  fontsize_col = 6
)

#################################################
# Save PNG
#################################################

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col = annotation_col,
  annotation_names_col = TRUE,
  main = "Sample-to-Sample Distance after VST normalization - menF-ubiC (GSE224889)",
  filename = "figures/GSE224889/Sample_distance_heatmap_GSE224889-menF-ubiC.png",
  width = 6,
  height = 5,
  fontsize = 7,
  fontsize_row = 6,
  fontsize_col = 6
)

#################################################
# Differential expression
#################################################

res <- results(
  dds,
  contrast = c(
    "condition",
    "ALE",
    "KO"
  )
)

#################################################
# Summary
#################################################

resultsNames(dds)

summary(res)

#################################################
# Convert results to data frame
#################################################

deg <- as.data.frame(res)

deg$Gene <- rownames(deg)

#################################################
# Remove NA adjusted p-values
#################################################

deg <- deg[
  !is.na(deg$padj),
]

#################################################
# Sort by adjusted p-value
#################################################

deg <- deg[
  order(deg$padj),
]

#################################################
# remove duplicates (keep best one)
#################################################

deg <- deg[!duplicated(deg$Gene), ]


#################################################
# Sort by adjusted p-value
#################################################

deg <- deg[
  order(deg$padj),
]

#################################################
# Save full DEG table
#################################################

write.csv(
  deg,
  "results/GSE224889/GSE224889_DEGs-menF-ubiC.csv",
  row.names = FALSE
)


#################################################
# Create meta-analysis input file
#################################################

meta_deg <- deg[, c(
  "Gene",
  "log2FoldChange",
  "padj"
)]

colnames(meta_deg) <- c(
  "Gene",
  "logFC",
  "adj.P.Val"
)

#################################################
# Check for duplicates
#################################################

sum(duplicated(rownames(counts)))
sum(duplicated(meta_deg$Gene))

#################################################
# Sort by adjusted p-value and log fold change
#################################################

meta_deg <- meta_deg[
  order(meta_deg$adj.P.Val,
        -abs(meta_deg$logFC)),
]

#################################################
# Save meta-analysis file
#################################################

write.csv(
  meta_deg,
  "results/GSE224889/GSE224889_meta_input-menF-ubiC.csv",
  row.names = FALSE
)

#################################################
# Significant DEGs
#################################################

deg_sig <- subset(
  meta_deg,
  adj.P.Val < 0.05 &
    abs(logFC) > 1
)

cat(
  "Significant DEGs:",
  nrow(deg_sig),
  "\n"
)
cat(
  "Upregulated:",
  nrow(subset(deg_sig, logFC > 1)),
  "\n"
)
cat(
  "Downregulated:",
  nrow(subset(deg_sig, logFC < -1)),
  "\n"
)

write.csv(
  deg_sig,
  "results/GSE224889/GSE224889_significant-menF-ubiC.csv",
  row.names = FALSE
)
#################################################
# Summary statistics
#################################################

cat(
  "Total genes:",
  nrow(meta_deg),
  "\n"
)

cat(
  "Significant genes:",
  sum(
    meta_deg$adj.P.Val < 0.05 &
      abs(meta_deg$logFC) > 1,
    na.rm = TRUE
  ),
  "\n"
)

sum(
  meta_deg$adj.P.Val < 0.05,
  na.rm = TRUE
)
resultsNames(dds)

summary(res)

head(
  deg[,c(
    "Gene",
    "log2FoldChange",
    "pvalue",
    "padj"
  )],
  20
)
#################################################
# Upregulated genes
#################################################

up_genes <- subset(
  meta_deg,
  adj.P.Val < 0.05 &
    logFC > 1
)

#################################################
# Downregulated genes
#################################################

down_genes <- subset(
  meta_deg,
  adj.P.Val < 0.05 &
    logFC < -1
)

cat(
  "Upregulated:",
  nrow(up_genes),
  "\n"
)

cat(
  "Downregulated:",
  nrow(down_genes),
  "\n"
)

#################################################
# Save upregulated genes
#################################################

write.csv(
  up_genes,
  "results/GSE224889/GSE224889_upregulated-menF-ubiC.csv",
  row.names = FALSE
)

#################################################
# Save downregulated genes
#################################################

write.csv(
  down_genes,
  "results/GSE224889/GSE224889_downregulated-menF-ubiC.csv",
  row.names = FALSE
)

#################################################
# Show top genes
#################################################

head(meta_deg, 20)

#################################################
# DESeq2 summary
#################################################

summary(res)

#################################################
# LFC Shrinkage
#################################################
resultsNames(dds)


library(apeglm)

res_shrink <- lfcShrink(
  dds,
  coef = "condition_ALE_vs_KO",
  type = "apeglm"
)

#################################################
# Convert to data frame
#################################################

deg_shrink <- as.data.frame(res_shrink)

deg_shrink$Gene <- rownames(deg_shrink)

deg_shrink <- deg_shrink[
  !is.na(deg_shrink$padj),
]

deg_shrink <- deg_shrink[
  order(deg_shrink$padj),
]

#################################################
# Top genes for labels
#################################################

rownames(deg_shrink) <- deg_shrink$Gene

topGenes <- head(rownames(deg_shrink), 10)

#################################################
# Volcano Plot
#################################################

library(EnhancedVolcano)
library(ggplot2)


(
  EnhancedVolcano(
    deg_shrink,
    lab = rownames(deg_shrink),
    selectLab = topGenes, 
    x = 'log2FoldChange',
    y = 'padj',
    xlab = 'log2 Fold Change',
    ylab = '-log10 Adjusted P-value',
    title = 'menF-ubiC KO  adapted vs menF-ubiC KO (GSE224889)',
    subtitle = '',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.5,
    labSize = 6,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = TRUE,
    colAlpha = 4/5,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1,
    colConnectors = 'black',
    max.overlaps = 20,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    border = 'full',
    borderWidth = 1.2,
    borderColour = 'black'
  )
) +
    theme(
    plot.title = element_text(hjust = 0.5),     
    plot.subtitle = element_text(hjust = 0.5) 
    )

#################################################
# Save PNG
#################################################

ggsave(
  "figures/GSE224889/Volcano_plot_shrink_GSE224889-menF-ubiC.png",
  width = 9.5,
  height = 9.5,
  dpi = 600
)

ggsave(
  "figures/GSE224889/Volcano_plot_shrink_GSE224889-menF-ubiC.pdf",
  width = 9.5,
  height = 9.5
)

#################################################
# MA Plot - Publication Quality
#################################################

library(ggplot2)
library(ggrepel)

# Convert to data frame
ma_df <- as.data.frame(res_shrink)

ma_df$Gene <- rownames(ma_df)

# Remove NAs
ma_df <- ma_df[
  !is.na(ma_df$padj) &
  !is.na(ma_df$baseMean),
]

#################################################
# Classification
#################################################

ma_df$Status <- "Not Significant"

ma_df$Status[
  ma_df$padj < 0.05 &
    ma_df$log2FoldChange > 1
] <- "Upregulated"

ma_df$Status[
  ma_df$padj < 0.05 &
    ma_df$log2FoldChange < -1
] <- "Downregulated"

#################################################
# Top genes for labeling
#################################################

top_genes <- ma_df[
  order(ma_df$padj),
][1:15, ]

#################################################
# Plot
#################################################

p_ma <- ggplot(
  ma_df,
  aes(
    x = log10(baseMean + 1),
    y = log2FoldChange
  )
) +

  geom_point(
    aes(color = Status),
    alpha = 0.75,
    size = 1.8
  ) +

  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_hline(
    yintercept = 0,
    color = "black",
    linewidth = 0.4
  ) +

  geom_text_repel(
    data = top_genes,
    aes(label = Gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2
  ) +

  scale_color_manual(
    values = c(
      "Upregulated" = "#D73027",
      "Downregulated" = "#4575B4",
      "Not Significant" = "grey80"
    )
  ) +

  labs(
    title = "menF-ubiC KO  adapted vs menF-ubiC KO (GSE224889)",
    x = expression(Log[10]~"(Mean normalized counts + 1)"),
    y = expression(Log[2]~Fold~Change),
    color = NULL
  ) +

  theme_bw(base_size = 13) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),

    legend.position = "top",

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      linewidth = 0.8
    )
  )

#################################################
# Show plot
#################################################

p_ma

#################################################
# Save figure
#################################################

ggsave(
  "figures/GSE224889/MA_plot_GSE224889-menF-ubiC.png",
  plot = p_ma,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/GSE224889/MA_plot_GSE224889-menF-ubiC.pdf",
  plot = p_ma,
  width = 8,
  height = 6
)


#################################################
# Signed Z-score
#################################################

meta_deg$P <- meta_deg$adj.P.Val

meta_deg$P[meta_deg$P < 1e-300] <- 1e-300

meta_deg$SignedZ <-
    sign(meta_deg$logFC) *
    qnorm(meta_deg$P/2, lower.tail = FALSE)

meta_deg <- meta_deg[
    order(-abs(meta_deg$SignedZ)),
]

write.csv(
    meta_deg,
    "results/GSE224889/GSE224889_menF-ubiC_SignedZ.csv",
    row.names = FALSE
)
