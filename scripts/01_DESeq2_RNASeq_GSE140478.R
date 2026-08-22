# Clear the workspace
rm(list = ls())

library(dplyr)
library(readxl)
library(DESeq2)
library(ggplot2)

#################################################
# Load packages
#################################################

library(data.table)

#################################################
# Extract GSE140478 RAW files
#################################################

raw_dir <- "data/GSE140478/GSE140478_RAW"

dir.create(
  raw_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

untar(
  "data/GSE140478/GSE140478_RAW.tar",
  exdir = raw_dir
)

#################################################
# Read GSE140478 expression files
#################################################

files <- list.files(
  raw_dir,
  pattern = "\\.csv\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

cat(
  "Number of expression files:",
  length(files),
  "\n"
)

#################################################
# Read individual expression files
#################################################

expr_list <- lapply(files, function(f){
  x <- fread(f)
  x <- x[, 2:3]
  colnames(x)[1] <- "Gene ID"
  sample_name <- basename(f)
  sample_name <- sub("\\.csv\\.gz$", "", sample_name)
  sample_name <- sub("GSM[0-9]+_", "", sample_name)
  colnames(x)[2] <- sample_name
  x
})

#################################################
# Merge all samples into expression matrix
#################################################

expr <- Reduce(
  function(x, y)
    merge(
      x,
      y,
      by = "Gene ID",
      all = TRUE
    ),
  expr_list
)

expr[is.na(expr)] <- 0

#################################################
# Merge duplicated Gene IDs
#################################################

expr$`Gene ID` <- sub(
  "_[0-9]+$",
  "",
  expr$`Gene ID`
)

expr <- as.data.table(expr)

expr <- expr[
  ,
  lapply(.SD, sum),
  by = `Gene ID`
]

#################################################
# Summary
#################################################

cat(
  "Genes:",
  nrow(expr),
  "\n"
)

cat(
  "Samples:",
  ncol(expr) - 1,
  "\n"
)

head(expr)
colnames(expr)

#################################################
# Keep only 44C samples
#################################################

keep_samples <- c(
  "Gene ID",
  grep(
    "_44C$",
    colnames(expr),
    value = TRUE
  )
)

expr <- expr[, ..keep_samples]

expr[is.na(expr)] <- 0
#################################################
# Read annotation
#################################################

library(rtracklayer)
library(dplyr)

gff <- import(
  "data/Ecoli_K12_MG1655_annotation.gff"
)

annot <- data.frame(
  locus_tag = gff$locus_tag,
  gene = gff$gene
)

annot_map <- annot %>%
  filter(!is.na(locus_tag)) %>%
  distinct(locus_tag, gene)

#################################################
# Replace Gene IDs with Gene Names when available
#################################################

expr_annot <- expr %>%
  left_join(
    annot_map,
    by = c("Gene ID" = "locus_tag")
  ) %>%
  mutate(
    `Gene ID` = ifelse(
      !is.na(gene) & gene != "",
      gene,
      `Gene ID`
    )
  ) %>%
  select(-gene)

#################################################
# Summary
#################################################

cat(
  "Mapped genes:",
  sum(expr$`Gene ID` %in% annot_map$locus_tag),
  "\n"
)

#################################################
# Collapse duplicated gene symbols
#################################################

cat(
  "Duplicated genes before collapse:",
  sum(duplicated(expr_annot$`Gene ID`)),
  "\n"
)

expr_annot <- expr_annot %>%
  group_by(`Gene ID`) %>%
  summarise(
    across(
      where(is.numeric),
      sum,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

#################################################
# Set gene IDs as row names
#################################################

expr_annot <- as.data.frame(expr_annot)

gene_names <- expr_annot$`Gene ID`

expr_annot$`Gene ID` <- NULL

expr_matrix <- as.matrix(expr_annot)

rownames(expr_matrix) <- gene_names

storage.mode(expr_matrix) <- "numeric"

head(rownames(expr_matrix))

#################################################
# Quality control
#################################################

cat(
  "Duplicated genes:",
  sum(
    duplicated(
      rownames(expr_matrix)
    )
  ),
  "\n"
)

summary(expr_matrix)

dim(expr_matrix)
#################################################
# Filter low expression genes
#################################################
cat(
  "Genes before filtering:",
  nrow(expr_matrix),
  "\n"
)

keep <- rowSums(expr_matrix >= 10) >= 3
expr_matrix <- expr_matrix[keep, ]

cat(
  "Genes after filtering:",
  nrow(expr_matrix),
  "\n"
)

#################################################
# Metadata
#################################################

condition <- factor(c(
  rep("Ancestor",2),
  rep("Evolved",16)
))

coldata <- data.frame(
  row.names = colnames(expr_matrix),
  condition
)

coldata$condition <- factor(
    coldata$condition,
    levels = c("Ancestor","Evolved")
)

levels(coldata$condition)

data.frame(
  Sample = colnames(expr_matrix),
  Condition = condition
)

#################################################
# DESeq2
#################################################
condition <- factor(condition, levels = c("Ancestor", "Evolved"))

dds <- DESeqDataSetFromMatrix(
  countData = expr_matrix,
  colData = coldata,
  design = ~ condition
)

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)

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

vst_long <- reshape2::melt(
  vst_df,
  id.vars = "Gene",
  variable.name = "Sample",
  value.name = "Expression"
)

#################################################
# Sample annotation
#################################################

vst_long$Condition <- ifelse(
  grepl("^42C_A3", vst_long$Sample),
  "Ancestor",
  "Evolved"
)

vst_long$Condition <- factor(
  vst_long$Condition,
  levels = c(
    "Ancestor",
    "Evolved"
  )
)


levels(vst_long$Condition)
resultsNames(dds)

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
  labs(
    title = "VST-normalized expression distribution(GSE140478)",
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
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

p_box

#################################################
# Save boxplot
#################################################

dir.create(
  "figures/GSE140478",
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  "figures/GSE140478/Expression_boxplot_GSE140478.png",
  plot = p_box,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  "figures/GSE140478/Expression_boxplot_GSE140478.pdf",
  plot = p_box,
  width = 7,
  height = 5
)

#################################################
# PCA
#################################################

# create new folder for results
dir.create("results/GSE140478", showWarnings = FALSE)

# Save the normalized expression data to a CSV file
write.csv(assay(vsd), "results/GSE140478/VST_expression.csv")

p <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA of VST-normalized counts (GSE140478)") +
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
  "figures/GSE140478/PCA_plot_GSE140478.png",
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/GSE140478/PCA_plot_GSE140478.pdf",
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
  main = "Sample-to-Sample Distance after VST normalization (GSE140478)",
  filename = "figures/GSE140478/Sample_distance_heatmap_GSE140478.pdf",
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
  main = "Sample-to-Sample Distance after VST normalization (GSE140478)",
  filename = "figures/GSE140478/Sample_distance_heatmap_GSE140478.png",
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
    "Evolved",
    "Ancestor"
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
  "results/GSE140478/GSE140478_DEGs.csv",
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
  "results/GSE140478/GSE140478_meta_input.csv",
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
  "results/GSE140478/GSE140478_significant.csv",
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
  "results/GSE140478/GSE140478_upregulated.csv",
  row.names = FALSE
)

#################################################
# Save downregulated genes
#################################################

write.csv(
  down_genes,
  "results/GSE140478/GSE140478_downregulated.csv",
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
  coef = "condition_Evolved_vs_Ancestor",
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
    title = 'Evolved vs Ancestor (GSE140478)',
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
  "figures/GSE140478/Volcano_plot_shrink_GSE140478.png",
  width = 9.5,
  height = 9.5,
  dpi = 600
)

ggsave(
  "figures/GSE140478/Volcano_plot_shrink_GSE140478.pdf",
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
    title = "Evolved vs Ancestor (GSE140478)",
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
  "figures/GSE140478/MA_plot_GSE140478.png",
  plot = p_ma,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/GSE140478/MA_plot_GSE140478.pdf",
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
    "results/GSE140478/GSE140478_SignedZ.csv",
    row.names = FALSE
)
