# Clear the workspace
rm(list = ls())

library(data.table)
library(dplyr)
library(readr)
library(readxl)
library(limma)
library(ggplot2)

#################################################
# Read GSE17276 log2 ratio matrix
# Values represent log2(Evolved / Ancestor)
#################################################

expr_data <- as.data.frame(
  read_excel(
    "data/GSE17276/GSE17276.xlsx"
  )
)

dim(expr_data)
head(expr_data[,1:10])

str(expr_data)

#################################################
# Prepare gene identifiers
# Keep ORF (locus tag) as gene ID
#################################################

expr_data$ORF <- as.character(expr_data$ORF)

#################################################
# Remove duplicated ORF IDs if present
#################################################

sum(duplicated(expr_data$ORF))

expr_data <- expr_data %>%
  group_by(ORF) %>%
  summarise(
    across(
      where(is.numeric),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

#################################################
# Create expression matrix (remove ID_REF and ORF)
#################################################

expr_matrix <- as.matrix(
  expr_data[, -(1:2)]
)

rownames(expr_matrix) <- expr_data$ORF

#################################################
# Convert to numeric
#################################################

expr_matrix <- apply(
  expr_matrix,
  2,
  as.numeric
)

rownames(expr_matrix) <- expr_data$ORF

#################################################
# Check matrix
#################################################

dim(expr_matrix)

head(expr_matrix[,1:5])

summary(expr_matrix)
#################################################
# Convert NaN values to NA
#################################################

expr_matrix[is.nan(expr_matrix)] <- NA

#################################################
# Check missing values
#################################################

sum(is.na(expr_matrix))
#################################################
# Check number of available measurements per gene
#################################################

gene_coverage <- rowSums(!is.na(expr_matrix))

summary(gene_coverage)
table(gene_coverage)

#################################################
# Select evolved lineage samples (exclude consortium)
#################################################

samples <- c(
  "GSM432640",
  "GSM432641",
  "GSM432642",
  "GSM432643",
  "GSM432644",
  "GSM432645",
  "GSM432646",
  "GSM432647",
  "GSM432648",
  "GSM432649",
  "GSM432650",
  "GSM432651",
  "GSM432652",
  "GSM432653",
  "GSM432654",
  "GSM432655",
  "GSM432656",
  "GSM432657",
  "GSM432658",
  "GSM432659",
  "GSM432660",
  "GSM432661",
  "GSM432662",
  "GSM432663",
  "GSM432664",
  "GSM432665",
  "GSM432666",
  "GSM432667",
  "GSM432668",
  "GSM432669",
  "GSM432670",
  "GSM432671",
  "GSM432672",
  "GSM432673",
  "GSM432674"
)

expr_matrix <- expr_matrix[, samples]

cat(
  "Selected samples:",
  ncol(expr_matrix),
  "\n"
)
#################################################
# Filter genes with sufficient measurements
#################################################

keep <- rowSums(!is.na(expr_matrix)) >= 10

expr_matrix_filtered <- expr_matrix[keep, ]

dim(expr_matrix_filtered)
#################################################
# Check remaining missing values
#################################################

sum(is.na(expr_matrix_filtered))
#################################################
# Check sample-wise missing values
#################################################

colSums(is.na(expr_matrix_filtered))

#################################################
# Prepare GSE17276 expression matrix
#################################################
expr_sub <- as.data.frame(expr_matrix_filtered)
expr_sub$ORF <- rownames(expr_sub)
expr_sub <- expr_sub[, c("ORF", setdiff(colnames(expr_sub), "ORF"))]
cat("Genes:", nrow(expr_sub), "\n")
cat("Samples:", ncol(expr_sub)-1, "\n")

#################################################
# Read annotation
#################################################
library(rtracklayer)
gff <- import("data/Ecoli_K12_MG1655_annotation.gff")
annot <- data.frame(
  locus_tag = gff$locus_tag,
  gene = gff$gene
)
annot_map <- annot %>%
  filter(!is.na(locus_tag)) %>%
  distinct(locus_tag, gene)
#################################################
# Replace ORF with gene symbols when available
#################################################
expr_annot <- expr_sub %>%
  left_join(
    annot_map,
    by = c("ORF" = "locus_tag")
  ) %>%
  mutate(
    Gene_ID = ifelse(
      !is.na(gene) & gene != "",
      gene,
      ORF
    )
  ) %>%
  select(
    Gene_ID,
    everything(),
    -ORF,
    -gene
  )
#################################################
# Mapping summary
#################################################
cat(
  "Mapped genes:",
  sum(expr_sub$ORF %in% annot_map$locus_tag),
  "\n"
)
#################################################
# Collapse duplicated gene symbols
#################################################
cat(
  "Duplicated genes before collapse:",
  sum(duplicated(expr_annot$Gene_ID)),
  "\n"
)
expr_annot <- expr_annot %>%
  group_by(Gene_ID) %>%
  summarise(
    across(
      where(is.numeric),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )
#################################################
# Final summary
#################################################
cat(
  "Genes after collapse:",
  nrow(expr_annot),
  "\n"
)
cat(
  "Duplicated genes after collapse:",
  sum(duplicated(expr_annot$Gene_ID)),
  "\n"
)
head(expr_annot)

#################################################
# Set gene IDs as row names
#################################################
expr_annot <- as.data.frame(expr_annot)
gene_names <- expr_annot$Gene_ID
expr_annot$Gene_ID <- NULL
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
# Check missing values
#################################################
cat(
  "Missing values:",
  sum(is.na(expr_matrix)),
  "\n"
)
#################################################
# Metadata
#################################################
samples <- colnames(expr_matrix)

cat(
  "Samples:",
  length(samples),
  "\n"
)

gene_coverage <- rowSums(!is.na(expr_matrix))
summary(gene_coverage)

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
expr_long <- expr_long[
  !is.na(expr_long$Expression),
]
#################################################
# Boxplot
#################################################
p_box <- ggplot(
  expr_long,
  aes(
    x = Sample,
    y = Expression
  )
) +
  geom_boxplot(
    linewidth = 0.4,
    outlier.size = 0.5
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Distribution of log2(Evolved/Ancestor) values (GSE17276)",
    subtitle = "Processed expression ratios",
    x = NULL,
    y = "log2(Evolved/Ancestor)"
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
    panel.grid.minor = element_blank()
  )
p_box
#################################################
# Save boxplot
#################################################
dir.create(
  "figures/GSE17276",
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  "figures/GSE17276/Expression_boxplot_GSE17276.png",
  plot = p_box,
  width = 9,
  height = 5,
  dpi = 600
)

ggsave(
  "figures/GSE17276/Expression_boxplot_GSE17276.pdf",
  plot = p_box,
  width = 9,
  height = 5
)

#################################################
# Prepare matrix for PCA (remove missing values)
#################################################

expr_matrix_pca <- expr_matrix_filtered

expr_matrix_pca <- apply(
  expr_matrix_pca,
  1,
  function(x) {
    x[is.na(x)] <- median(
      x,
      na.rm = TRUE
    )
    x
  }
)

expr_matrix_pca <- t(expr_matrix_pca)

dim(expr_matrix_pca)

sum(is.na(expr_matrix_pca))
#################################################
# PCA analysis
#################################################
pca <- prcomp(
  t(expr_matrix_pca),
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
# Define lineage information
#################################################

pca_df$Lineage <- c(
  rep("CV101", 3),
  rep("CV103", 3),
  rep("CV115", 3),
  rep("CV116", 2),
  rep("CV101", 3),
  rep("CV103", 3),
  rep("CV115", 3),
  rep("CV116", 3),
  rep("CV101", 3),
  rep("CV103", 3),
  rep("CV115", 3),
  rep("CV116", 3)
)

pca_df$Lineage <- factor(
  pca_df$Lineage,
  levels = c(
    "CV101",
    "CV103",
    "CV115",
    "CV116"
  )
)

#################################################
# PCA plot
#################################################

p <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Lineage,
    label = Sample
  )
) +
  geom_point(
    size = 3
  ) +
  labs(
    title = "PCA of log2(Evolved/Ancestor) values (GSE17276)",
    x = paste0(
      "PC1 (",
      round(summary(pca)$importance[2,1] * 100, 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(summary(pca)$importance[2,2] * 100, 1),
      "%)"
    ),
    color = "Lineage"
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
    legend.position = "top"
  )

p
#################################################
# Save PCA
#################################################
ggsave(
  "figures/GSE17276/PCA_plot_GSE17276.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 600
)
ggsave(
  "figures/GSE17276/PCA_plot_GSE17276.pdf",
  plot = p,
  width = 8,
  height = 6
)

#################################################
# Sample distance heatmap
#################################################

library(pheatmap)

#################################################
# Prepare matrix for heatmap
#################################################

expr_heatmap <- expr_matrix_filtered

expr_heatmap[is.na(expr_heatmap)] <- 0

#################################################
# Sample-to-sample distances
#################################################

sampleDists <- dist(
  t(expr_heatmap)
)

sampleDistMatrix <- as.matrix(
  sampleDists
)

#################################################
# Sample annotation
#################################################

annotation_col <- data.frame(
  Lineage = c(
    rep("CV101", 3),
    rep("CV103", 3),
    rep("CV115", 3),
    rep("CV116", 2),
    rep("CV101", 3),
    rep("CV103", 3),
    rep("CV115", 3),
    rep("CV116", 3),
    rep("CV101", 3),
    rep("CV103", 3),
    rep("CV115", 3),
    rep("CV116", 3)
  )
)

rownames(annotation_col) <- colnames(expr_heatmap)

annotation_col$Lineage <- factor(
  annotation_col$Lineage,
  levels = c(
    "CV101",
    "CV103",
    "CV115",
    "CV116"
  )
)

#################################################
# Save PDF
#################################################

dir.create(
  "figures/GSE17276",
  showWarnings = FALSE,
  recursive = TRUE
)

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col = annotation_col,
  annotation_names_col = TRUE,
  main = "Sample-to-Sample Distance of log2(Evolved/Ancestor) values (GSE17276)",
  filename = "figures/GSE17276/Sample_distance_heatmap_GSE17276_lineage.pdf",
  width = 7,
  height = 6,
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
  main = "Sample-to-Sample Distance of log2(Evolved/Ancestor) values (GSE17276)",
  filename = "figures/GSE17276/Sample_distance_heatmap_GSE17276_lineage.png",
  width = 7,
  height = 6,
  fontsize = 7,
  fontsize_row = 6,
  fontsize_col = 6
)

#################################################
# Differential expression analysis
#
# Important note:
# The GSE17276 matrix does NOT contain raw expression values.
# The GEO processed data represent pre-calculated log2 fold-change
# values for evolved populations relative to the ancestor:
#
# log2FC = log2(Evolved / Ancestor)
#
# Therefore, each sample column already represents an
# evolved-versus-ancestor comparison.
#
# Positive values indicate genes upregulated in evolved populations
# compared with the ancestor, while negative values indicate
# genes downregulated in evolved populations.
#
# Because an ancestor expression matrix is not available, a direct
# two-group comparison (Ancestor vs Evolved) cannot be performed.
#
# Instead, a one-sample limma model is applied to test whether the
# average log2FC across evolved populations significantly differs
# from zero:
#
# Mean(log2(Evolved / Ancestor)) != 0
#
# This identifies genes showing consistent evolutionary changes
# across independent evolved populations.
#################################################

library(limma)

#################################################
# Design matrix
#################################################
design <- matrix(
  1,
  ncol(expr_matrix),
  1
)
colnames(design) <- "Mean_log2FC"
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
# Mean log2(Evolved/Ancestor)
#################################################
res <- topTable(
  fit,
  coef = "Mean_log2FC",
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
dir.create(
  "results/GSE17276",
  showWarnings = FALSE,
  recursive = TRUE
)
write.csv(
  deg,
  "results/GSE17276/GSE17276_DEGs_Evolved_vs_Ancestor.csv",
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
  "results/GSE17276/GSE17276_meta_input.csv",
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
  "results/GSE17276/GSE17276_significant_DEGs.csv",
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
  "results/GSE17276/GSE17276_upregulated.csv",
  row.names = FALSE
)
#################################################
# Save downregulated genes
#################################################
write.csv(
  down_genes,
  "results/GSE17276/GSE17276_downregulated.csv",
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
  "results/GSE17276/GSE17276_SignedZ.csv",
  row.names = FALSE
)
