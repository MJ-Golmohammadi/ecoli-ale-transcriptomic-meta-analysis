# Clear the workspace
rm(list = ls())

library(data.table)
library(dplyr)
library(readxl)
library(limma)
library(ggplot2)

#################################################
# Read normalized expression matrix
#################################################

counts <- as.data.frame(
  read_excel(
    "data/GSE206196/GSE206196_dataset_HMS_normalized.xlsx"
  )
)

colnames(counts)
head(counts)
dim(counts)
############################################################
# Rename gene ID column and remove column2
############################################################

colnames(counts)[1] <- "Gene ID"
counts <- counts[, -2]
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
  "data/GSE206196/GSE206196_dataset_HMS_normalized_annotated.csv",
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
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

############################################################
# Set gene IDs as row names
############################################################

rownames(counts) <- counts$`Gene ID`

#################################################
# Keep selected samples
#################################################
# to solve the problem of converting genes to numbers
counts <- as.data.frame(counts)

samples <- c(
  "WT_2_1",
  "WT_3_1",
  "WT_1_1",
  "4E1_1",
  "4E_3_1",
  "4E_2_1"
)

expr <- counts[, c(
  "Gene ID",
  samples
)]

dim(expr)
#################################################
# Remove genes without IDs
#################################################

expr <- expr %>%
  filter(
    !is.na(`Gene ID`),
    `Gene ID` != ""
  )

dim(expr)
#################################################
# Collapse duplicated gene IDs
#################################################

expr <- expr %>%
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
# Create expression matrix
#################################################

expr_matrix <- as.matrix(
  expr[, -1]
)

rownames(expr_matrix) <- expr$`Gene ID`

mode(expr_matrix) <- "numeric"

#################################################
# Quality check
#################################################

dim(expr_matrix)

head(expr_matrix)

summary(expr_matrix)

#################################################
# Expression distribution
#################################################

summary(expr_matrix)

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

#################################################
# Define conditions
#################################################

expr_long$Condition <- ifelse(
  expr_long$Sample %in% c(
    "WT_2_1",
    "WT_3_1",
    "WT_1_1"
  ),
  "WT",
  "4E"
)

expr_long$Condition <- factor(
  expr_long$Condition,
  levels = c(
    "WT",
    "4E"
  )
)

#################################################
# Boxplot
#################################################

p_box <- ggplot(
  expr_long,
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
    title = "Normalized expression distribution (GSE206196)",
    subtitle = "Log-transformed expression values",
    x = NULL,
    y = "Expression"
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
  "figures/GSE206196",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results/GSE206196",
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  "figures/GSE206196/Expression_boxplot_GSE206196.png",
  plot = p_box,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  "figures/GSE206196/Expression_boxplot_GSE206196.pdf",
  plot = p_box,
  width = 7,
  height = 5
)

library(ggplot2)

pca <- prcomp(
  t(expr_matrix),
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

pca_df$Condition <- ifelse(
  pca_df$Sample %in% c(
    "WT_2_1",
    "WT_3_1",
    "WT_1_1"
  ),
  "WT",
  "4E"
)

pca_df$Condition <- factor(
  pca_df$Condition,
  levels = c(
    "WT",
    "4E"
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
    color = Condition,
    label = Sample
  )
) +

  geom_point(
    size = 3
  ) +

  labs(
    title = "PCA of normalized log-expression values (GSE206196)",
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
# Save PCA plot
#################################################

ggsave(
  "figures/GSE206196/PCA_GSE206196.png",
  plot = p,
  width = 6,
  height = 5,
  dpi = 600
)

ggsave(
  "figures/GSE206196/PCA_GSE206196.pdf",
  plot = p,
  width = 6,
  height = 5
)

library(pheatmap)

#################################################
# Sample distance matrix
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
  Condition = ifelse(
    colnames(expr_matrix) %in% c(
      "WT_2_1",
      "WT_3_1",
      "WT_1_1"
    ),
    "WT",
    "4E"
  )
)

rownames(annotation_col) <- colnames(expr_matrix)

#################################################
# Sample distance heatmap
#################################################

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col = annotation_col,
  annotation_names_col = TRUE,
  main = "Sample-to-sample distance of normalized log-expression values (GSE206196)",
  filename = "figures/GSE206196/Sample_distance_heatmap_GSE206196.pdf",
  width = 6,
  height = 5,
  fontsize = 7,
  fontsize_row = 6,
  fontsize_col = 6
)

#################################################
# Save PNG
#################################################

png(
  "figures/GSE206196/Sample_distance_heatmap_GSE206196.png",
  width = 1800,
  height = 1500,
  res = 300
)

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col = annotation_col,
  annotation_names_col = TRUE,
  main = "Sample-to-sample distance of normalized log-expression values (GSE206196)",
  fontsize = 7,
  fontsize_row = 6,
  fontsize_col = 6
)

dev.off()

#################################################
# Sample information
#################################################

coldata <- data.frame(
  row.names = colnames(expr_matrix),
  condition = c(
    rep("WT", 3),
    rep("4E", 3)
  )
)

coldata$condition <- factor(
  coldata$condition,
  levels = c(
    "WT",
    "4E"
  )
)

levels(coldata$condition)
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
# Differentially expressed genes
#################################################

res <- topTable(
  fit,
  coef = "condition4E",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

res <- res[order(res$adj.P.Val), ]
#################################################
# Add gene names
#################################################

res$Gene <- rownames(res)

res <- res %>%
  relocate(
    Gene
  )

#################################################
# Summary
#################################################

summary(res)

head(res)

#################################################
# Significant genes
#################################################

deg <- res %>%
  filter(
    adj.P.Val < 0.05 &
      abs(logFC) >= 1
  )

nrow(deg)

table(
  ifelse(
    deg$logFC > 0,
    "Upregulated",
    "Downregulated"
  )
)

#################################################
# Save results
#################################################

write.csv(
  res,
  "results/GSE206196/All_DEGs_GSE206196.csv",
  row.names = FALSE
)

write.csv(
  deg,
  "results/GSE206196/Significant_DEGs_GSE206196.csv",
  row.names = FALSE
)

#################################################
# Upregulated genes
#################################################

up_genes <- subset(
  res,
  adj.P.Val < 0.05 &
    logFC > 1
)

#################################################
# Downregulated genes
#################################################

down_genes <- subset(
  res,
  adj.P.Val < 0.05 &
    logFC < -1
)

#################################################
# Save upregulated genes
#################################################

write.csv(
  up_genes,
  "results/GSE206196/GSE206196_upregulated.csv",
  row.names = FALSE
)

#################################################
# Save downregulated genes
#################################################

write.csv(
  down_genes,
  "results/GSE206196/GSE206196_downregulated.csv",
  row.names = FALSE
)

#################################################
# Signed Z-score
#################################################

res$P <- res$adj.P.Val

res$P[
  res$P < 1e-300
] <- 1e-300

res$SignedZ <-
  sign(res$logFC) *
  qnorm(
    res$P / 2,
    lower.tail = FALSE
  )

res <- res[
  order(
    -abs(res$SignedZ)
  ),
]

#################################################
# Save Signed Z file
#################################################

write.csv(
  res,
  "results/GSE206196/GSE206196_SignedZ.csv",
  row.names = FALSE
)

#################################################
# Select top genes
#################################################

topGenes <- res %>%
  filter(
    adj.P.Val < 0.05
  ) %>%
  arrange(
    adj.P.Val
  ) %>%
  slice_head(
    n = 10
  ) %>%
  pull(
    Gene
  )

#################################################
# Volcano plot
#################################################
library(EnhancedVolcano)

volcano_df <- res

rownames(volcano_df) <- volcano_df$Gene

(
p_volcano <- EnhancedVolcano(
  volcano_df,
  lab = rownames(volcano_df),
  selectLab = topGenes,
  x = "logFC",
  y = "adj.P.Val",
  xlab = "log2 Fold Change",
  ylab = "-log10 Adjusted P-value",
  title = "4E vs WT (GSE206196)",
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

p_volcano

#################################################
# Save volcano plot
#################################################

ggsave(
  "figures/GSE206196/Volcano_GSE206196.png",
  width = 8,
  height = 7,
  dpi = 600
)

ggsave(
  "figures/GSE206196/Volcano_GSE206196.pdf",
  width = 8,
  height = 7
)


#################################################
# MA Plot
#################################################

library(ggrepel)

ma_df <- res

ma_df$Gene <- rownames(res)

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
    title =
      "MA plot - ALE vs Parent (GSE206196)",

    subtitle =
      "",

    x = "Average Expression",

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
  "figures/GSE206196/MA_plot_GSE206196.png",
  plot = p_ma,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/GSE206196/MA_plot_GSE206196.pdf",
  plot = p_ma,
  width = 8,
  height = 6
)

