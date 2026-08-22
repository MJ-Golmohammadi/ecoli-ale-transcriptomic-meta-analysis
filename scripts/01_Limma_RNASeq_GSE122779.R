# Clear the workspace
rm(list = ls())

library(data.table)
library(dplyr)
library(readr)
library(limma)
library(ggplot2)

library(readxl)

#################################################
# Read TPM matrix
#################################################

counts <- read_csv(
  "data/GSE122779/GSE122779_log_tpm_GEO.csv.gz",
  show_col_types = FALSE
)

colnames(counts)
dim(counts)

str(counts)

#################################################
# Collapse duplicated gene IDs
# (remove suffixes such as _1, _2, _3)
#################################################

counts$`Gene ID` <- sub(
  "_[0-9]+$",
  "",
  counts$`Gene ID`
)

counts <- counts %>%
  group_by(`Gene ID`) %>%
  summarise(
    across(
      where(is.numeric),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

#################################################
# Summary
#################################################

dim(counts)
head(counts)
sapply(counts, class)
#################################################
# Select samples for GSE122779 ALE analysis
#################################################

samples <- c(
  "QN20A",
  "QN20B",
  "QN21A",
  "QN21B",
  "QN22A",
  "QN22B",
  "QN23A",
  "QN23B"
)

#################################################
# Create selected expression matrix
#################################################

expr_sub <- counts[, c("Gene ID", samples)]

cat(
"Genes:",
nrow(expr_sub),
"\n"
)

cat(
"Samples:",
ncol(expr_sub)-1,
"\n"
)

#################################################
# Read annotation
#################################################

library(rtracklayer)

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

expr_annot <- expr_sub %>%
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

#################################################
# Summary
#################################################

cat(
  "Mapped genes:",
  sum(
    expr_sub$`Gene ID` %in%
      annot_map$locus_tag
  ),
  "\n"
)

#################################################
# Collapse duplicated gene symbols
#################################################

cat(
  "Duplicated genes before collapse:",
  sum(
    duplicated(
      expr_annot$`Gene ID`
    )
  ),
  "\n"
)

expr_annot <- expr_annot %>%
  group_by(`Gene ID`) %>%
  summarise(
    across(
      where(is.numeric),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

head(expr_annot)

cat(
  "Duplicated genes before collapse:",
  sum(
    duplicated(
      expr_annot$`Gene ID`
    )
  ),
  "\n"
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


#################################################
# Filter low expression genes
#################################################

keep <- rowSums(expr_matrix > 1) >= 3

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
  "Ancestor", "Ancestor",
  rep("Evolved", 6)
))

coldata <- data.frame(
  row.names = colnames(expr_matrix),
  condition
)

coldata$condition <- factor(
  coldata$condition,
  levels = c(
    "Ancestor",
    "Evolved"
  )
)

levels(coldata$condition)

#################################################
# Quality control - Expression distribution
#################################################

library(reshape2)

expr_df <- as.data.frame(expr_matrix)

expr_df$Gene <- rownames(expr_df)

expr_long <- melt(
  expr_df,
  id.vars = "Gene",
  variable.name = "Sample",
  value.name = "Expression"
)

expr_long$Condition <- ifelse(
  expr_long$Sample %in% c(
    "QN20A",
    "QN20B"
  ),
  "Ancestor",
  "Evolved"
)

expr_long$Condition <- factor(
  expr_long$Condition,
  levels = c(
    "Ancestor",
    "Evolved"
  )
)

expr_long$Expression_plot <- expr_long$Expression
#################################################
# Boxplot
#################################################

p_box <- ggplot(
  expr_long,
  aes(
    x = Sample,
    y = Expression_plot,
    fill = Condition
  )
) +
  geom_boxplot(
    linewidth = 0.4,
    outlier.size = 0.5
  ) +
  labs(
	title = "Expression distribution (Evolved vs Ancestor, GSE122779)",
	subtitle = "log-transformed TPM values",
    x = NULL,
    y = "log TPM expression"
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

# Save boxplot
#################################################

dir.create(
  "figures/GSE122779",
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  "figures/GSE122779/Expression_boxplot_GSE122779.png",
  plot = p_box,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  "figures/GSE122779/Expression_boxplot_GSE122779.pdf",
  plot = p_box,
  width = 7,
  height = 5
)

#################################################
# PCA analysis
#################################################

library(ggplot2)

expr_matrix_log <- expr_matrix

pca <- prcomp(
  t(expr_matrix_log),
  scale. = TRUE
)

#################################################
# PCA dataframe
#################################################

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Sample = rownames(pca$x)
)

#################################################
# Define conditions
#################################################

pca_df$Condition <- coldata$condition[
  match(
    pca_df$Sample,
    rownames(coldata)
  )
]

#################################################
# PCA plot
#################################################

p <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Condition,
    label = Sample
  )
) +

  geom_point(
    size = 3
  ) +

  labs(
    title = "PCA of log expression values (Evolved vs Ancestor, GSE122779)",
    x = paste0(
      "PC1 (",
      round(summary(pca)$importance[2,1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(pca)$importance[2,2] * 100, 1),
      "%)"
    )
  ) +

  theme_bw(base_size = 13) +

  theme(
    plot.title = element_text(
      size = 11,
      hjust = 0.5,
      face = "bold"
    ),

    axis.title = element_text(
      size = 10
    ),

    axis.text = element_text(
      size = 10
    ),

    legend.title = element_text(
      size = 9
    ),

    legend.text = element_text(
      size = 8
    )
  )

p

#################################################
# Save PCA
#################################################

ggsave(
  "figures/GSE122779/PCA_plot_GSE122779.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/GSE122779/PCA_plot_GSE122779.pdf",
  plot = p,
  width = 8,
  height = 6
)


#################################################
## Sample distance heatmap
#################################################

library(pheatmap)

#################################################
# Sample-to-sample distances
#################################################

sampleDists <- dist(
  t(expr_matrix)
)

sampleDistMatrix <- as.matrix(
  sampleDists
)

#################################################
# Sample annotation
#################################################

annotation_col <- data.frame(
  Condition = coldata$condition
)

rownames(annotation_col) <- rownames(coldata)

#################################################
# Save PDF
#################################################

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col = annotation_col,
  annotation_names_col = TRUE,
  main = "Sample-to-Sample Distance of log expression values (Evolved vs Ancestor, GSE122779)",
  filename = "figures/GSE122779/Sample_distance_heatmap_GSE122779.pdf",
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
  main = "Sample-to-Sample Distance of log expression values (Evolved vs Ancestor, GSE122779)",
  filename = "figures/GSE122779/Sample_distance_heatmap_GSE122779.png",
  width = 6,
  height = 5,
  fontsize = 7,
  fontsize_row = 6,
  fontsize_col = 6
)

#################################################
# Differential expression analysis
#################################################

library(limma)

#################################################
# Design matrix
#################################################

design <- model.matrix(
  ~ condition,
  data = coldata
)

colnames(design)

#################################################
# Fit linear model
#################################################

fit <- lmFit(
  expr_matrix,
  design
)

fit <- eBayes(
  fit
)

#################################################
# Differential expression
#
# Evolved vs Ancestor
#################################################

res <- topTable(
  fit,
  coef = "conditionEvolved",
  number = Inf,
  adjust.method = "BH"
)

res$Gene <- rownames(res)

res <- res[, c(
  "Gene",
  "logFC",
  "AveExpr",
  "t",
  "P.Value",
  "adj.P.Val",
  "B"
)]

#################################################
# Summary
#################################################

summary(res)

head(res)

#################################################
# Convert results to data frame
#################################################

deg <- as.data.frame(res)

deg$Gene <- rownames(deg)

#################################################
# Remove NA adjusted p-values
#################################################

deg <- deg[
  !is.na(deg$adj.P.Val),
]

#################################################
# Sort by adjusted p-value
#################################################

deg <- deg[
  order(deg$adj.P.Val),
]

#################################################
# Remove duplicate genes
#################################################

deg <- deg[
  !duplicated(deg$Gene),
]

#################################################
# Save full DEG table
#################################################
# create results directory if it doesn't exist
dir.create(
  "results/GSE122779",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  deg,
  "results/GSE122779/GSE122779_DEGs.csv",
  row.names = FALSE
)

#################################################
# Create meta-analysis input file
#################################################

meta_deg <- deg[, c(
  "Gene",
  "logFC",
  "adj.P.Val"
)]

colnames(meta_deg) <- c(
  "Gene",
  "logFC",
  "adj.P.Val"
)

#################################################
# Sort by significance
#################################################

meta_deg <- meta_deg[
  order(
    meta_deg$adj.P.Val,
    -abs(meta_deg$logFC)
  ),
]

#################################################
# Save meta-analysis file
#################################################

write.csv(
  meta_deg,
  "results/GSE122779/GSE122779_meta_input.csv",
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
  nrow(subset(
    deg_sig,
    logFC > 1
  )),
  "\n"
)

cat(
  "Downregulated:",
  nrow(subset(
    deg_sig,
    logFC < -1
  )),
  "\n"
)

#################################################
# Save significant genes
#################################################

write.csv(
  deg_sig,
  "results/GSE122779/GSE122779_significantEvolved.csv",
  row.names = FALSE
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
  "results/GSE122779/GSE122779_upregulated.csv",
  row.names = FALSE
)

#################################################
# Save downregulated genes
#################################################

write.csv(
  down_genes,
  "results/GSE122779/GSE122779_downregulated.csv",
  row.names = FALSE
)

#################################################
# Signed Z-score
#################################################

meta_deg$P <- meta_deg$adj.P.Val

meta_deg$P[
  meta_deg$P < 1e-300
] <- 1e-300

meta_deg$SignedZ <-
  sign(meta_deg$logFC) *
  qnorm(
    meta_deg$P / 2,
    lower.tail = FALSE
  )

meta_deg <- meta_deg[
  order(
    -abs(meta_deg$SignedZ)
  ),
]

#################################################
# Save Signed Z file
#################################################

write.csv(
  meta_deg,
  "results/GSE122779/GSE122779_SignedZ.csv",
  row.names = FALSE
)

#################################################
# Volcano Plot
#################################################

library(EnhancedVolcano)

#################################################
# Prepare volcano data
#################################################

volcano_df <- deg

rownames(volcano_df) <- volcano_df$Gene

topGenes <- head(
  volcano_df$Gene[
    order(volcano_df$adj.P.Val)
  ],
  10
)

#################################################
# Volcano plot
#################################################
(
  EnhancedVolcano(
    volcano_df,
    lab = rownames(volcano_df),
    selectLab = topGenes,
    x = "logFC",
    y = "adj.P.Val",
    xlab = "log2 Fold Change",
    ylab = "-log10 Adjusted P-value",
    title = "Evolved vs Ancestor (GSE122779)",
    subtitle = "",
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.5,
    labSize = 6,
    labCol = "black",
    labFace = "bold",
    boxedLabels = TRUE,
    colAlpha = 4/5,
    legendPosition = "right",
    drawConnectors = TRUE,
    widthConnectors = 1,
    max.overlaps = 20,
    gridlines.major = FALSE,
    gridlines.minor = FALSE
  )
) +
  coord_cartesian(
    ylim = c(0, 6)
  )

#################################################
# Save Volcano plot
#################################################

ggsave(
  "figures/GSE122779/Volcano_plot_GSE122779.png",
  width = 9.5,
  height = 9.5,
  dpi = 600
)

ggsave(
  "figures/GSE122779/Volcano_plot_GSE122779.pdf",
  width = 9.5,
  height = 9.5
)


#################################################
# MA Plot
#################################################

library(ggrepel)

ma_df <- deg

ma_df$Gene <- rownames(deg)

#################################################
# Classification
#################################################

ma_df$Status <- "Not Significant"

ma_df$Status[
  ma_df$adj.P.Val < 0.05 &
    ma_df$logFC > 1
] <- "Upregulated"


ma_df$Status[
  ma_df$adj.P.Val < 0.05 &
    ma_df$logFC < -1
] <- "Downregulated"


#################################################
# Top genes for labeling
#################################################

top_genes <- ma_df[
  order(ma_df$adj.P.Val),
][1:15, ]


#################################################
# MA plot
#################################################
p_ma <- ggplot(
  ma_df,
  aes(
    x = AveExpr,
    y = logFC
  )
) +

  geom_point(
    aes(color = Status),
    alpha = 0.75,
    size = 2
  ) +

  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.5
  ) +

  geom_hline(
    yintercept = 0,
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
    title = "MA Plot: Evolved vs Ancestor (GSE122779)",
    subtitle = "",
    x = "Average log TPM expression",
    y = expression(log[2]~Fold~Change),
    color = NULL
  ) +

  theme_bw(base_size = 13) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      hjust = 0.5
    ),

    legend.position = "top",

    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      linewidth = 0.8
    )
  )
#################################################
# Show MA plot
#################################################

p_ma


#################################################
# Save MA plot
#################################################

ggsave(
  "figures/GSE122779/MA_plot_GSE122779.png",
  plot = p_ma,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/GSE122779/MA_plot_GSE122779.pdf",
  plot = p_ma,
  width = 8,
  height = 6
)