# Clear workspace
rm(list = ls())

library(affy)
library(limma)
library(GEOquery)
library(dplyr)
library(ggplot2)

#################################################
# Read CEL files
#################################################
# Unpack the raw data
untar(
  "data/GSE33147/GSE33147_RAW.tar",
  exdir = "data/GSE33147"
)

# Load raw microarray data
affy_raw <- ReadAffy(
  celfile.path = "data/GSE33147"
)

#################################################
# RMA normalization
#################################################

eset <- rma(affy_raw)
expr_matrix <- exprs(eset)

dim(expr_matrix)
head(expr_matrix)

annotation(eset)
colnames(expr_matrix)

#################################################
# Keep selected Glycerol samples
#################################################

samples <- c(

  # Ancestor (WT glycerol)
  "GSM820912.CEL.gz",
  "GSM820913.CEL.gz",
  "GSM820914.CEL.gz",
  "GSM820915.CEL.gz",
  "GSM820916.CEL.gz",


  # Evolved (day 44)

  # glycerol1
  "GSM820875.CEL.gz",
  "GSM820876.CEL.gz",
  "GSM820877.CEL.gz",

  # glycerol2
  "GSM820881.CEL.gz",
  "GSM820882.CEL.gz",
  "GSM820883.CEL.gz",

  # glycerolA
  "GSM820887.CEL.gz",
  "GSM820888.CEL.gz",
  "GSM820889.CEL.gz",

  # glycerolB
  "GSM820893.CEL.gz",
  "GSM820894.CEL.gz",
  "GSM820895.CEL.gz",

  # glycerolC
  "GSM820897.CEL.gz",
  "GSM820898.CEL.gz",
  "GSM820899.CEL.gz",

  # glycerolD
  "GSM820903.CEL.gz",
  "GSM820904.CEL.gz",
  "GSM820905.CEL.gz",

  # glycerolE
  "GSM820909.CEL.gz",
  "GSM820910.CEL.gz",
  "GSM820911.CEL.gz"
)


expr_matrix <- expr_matrix[, samples]
#################################################
# Convert Probe IDs to ORF using GPL199
#################################################

expr_df <- data.frame(
  ID = rownames(expr_matrix),
  expr_matrix,
  check.names = FALSE
)


#################################################
# Check matching probe IDs
#################################################

cat(
  "Matched probes:",
  sum(expr_df$ID %in% annot$ID),
  "\n"
)

#################################################
# Merge probe annotation
#################################################

expr_df <- expr_df %>%
  left_join(
    annot,
    by = "ID"
  )


#################################################
# Remove probes without ORF
#################################################

expr_df <- expr_df %>%
  filter(
    !is.na(ORF),
    ORF != ""
  )


cat(
  "Probes with ORF:",
  nrow(expr_df),
  "\n"
)


#################################################
# Collapse duplicated probes per ORF
#################################################

expr_df <- expr_df %>%
  group_by(ORF) %>%
  summarise(
    across(
      where(is.numeric),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


cat(
  "Genes after ORF collapse:",
  nrow(expr_df),
  "\n"
)


#################################################
# Create ORF expression matrix
#################################################

expr_matrix <- as.matrix(
  expr_df[, -1]
)

rownames(expr_matrix) <- expr_df$ORF

mode(expr_matrix) <- "numeric"

dim(expr_matrix)
head(expr_matrix)
#################################################
# Convert ORF to Gene Symbol using GFF
#################################################

library(rtracklayer)

gff <- import(
  "data/Ecoli_K12_MG1655_annotation.gff"
)

annot_map <- data.frame(
  ORF = gff$locus_tag,
  Gene = gff$gene
) %>%
  filter(
    !is.na(ORF)
  ) %>%
  distinct(ORF, .keep_all = TRUE)

expr_annot <- expr_df %>%
  left_join(
    annot_map,
    by = "ORF"
  )

#################################################
# Replace ORF with Gene Symbol when available
#################################################

expr_annot <- expr_annot %>%
  mutate(
    `Gene ID` = ifelse(
      !is.na(Gene) &
        Gene != "",
      Gene,
      ORF
    )
  ) %>%
  select(
    `Gene ID`,
    everything(),
    -ORF,
    -Gene
  )

#################################################
# Remove genes without symbols
#################################################

expr_annot <- expr_annot %>%
  filter(
    !is.na(`Gene ID`),
    `Gene ID` != ""
  )

#################################################
# Collapse duplicated Gene Symbols
#################################################

cat(
  "Duplicated genes before collapse:",
  sum(
    duplicated(expr_annot$`Gene ID`)
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

cat(
  "Genes after Gene Symbol collapse:",
  nrow(expr_annot),
  "\n"
)

#################################################
# Final expression matrix
#################################################
colnames(expr_annot)
head(expr_annot)
summary(expr_annot)

expr_matrix <- as.matrix(
  expr_annot[, -1]
)

rownames(expr_matrix) <- expr_annot$`Gene ID`
mode(expr_matrix) <- "numeric"

dim(expr_matrix)
head(expr_matrix)

#################################################
# Rename samples for readability
#################################################
colnames(expr_matrix) <- c(

  "WT_Glycerol_1",
  "WT_Glycerol_2",
  "WT_Glycerol_3",
  "WT_Glycerol_4",
  "WT_Glycerol_5",

  "GlycerolEvolved_1",
  "GlycerolEvolved_2",
  "GlycerolEvolved_3",
  "GlycerolEvolved_4",
  "GlycerolEvolved_5",
  "GlycerolEvolved_6",
  "GlycerolEvolved_7",
  "GlycerolEvolved_8",
  "GlycerolEvolved_9",
  "GlycerolEvolved_10",
  "GlycerolEvolved_11",
  "GlycerolEvolved_12",
  "GlycerolEvolved_13",
  "GlycerolEvolved_14",
  "GlycerolEvolved_15",
  "GlycerolEvolved_16",
  "GlycerolEvolved_17",
  "GlycerolEvolved_18",
  "GlycerolEvolved_19",
  "GlycerolEvolved_20",
  "GlycerolEvolved_21"
)

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
  grepl(
    "^WT_",
    expr_long$Sample
  ),
  "WT",
  "Evolved"
)

expr_long$Condition <- factor(
  expr_long$Condition,
  levels = c(
    "WT",
    "Evolved"
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
    title = "Normalized expression distribution (GSE33147)",
    subtitle = "RMA-normalized log2 expression values (Glycerol WT vs Evolved)",
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
  "figures/GSE33147",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results/GSE33147",
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  "figures/GSE33147/Expression_boxplot_GSE33147_Glycerol.png",
  plot = p_box,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  "figures/GSE33147/Expression_boxplot_GSE33147_Glycerol.pdf",
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
  grepl(
    "^WT_",
    pca_df$Sample
  ),
  "WT",
  "Evolved"
)

pca_df$Condition <- factor(
  pca_df$Condition,
  levels = c(
    "WT",
    "Evolved"
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
    title = "PCA of RMA-normalized log2 expression values (GSE33147_Glycerol)",
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
  "figures/GSE33147/PCA_GSE33147_Glycerol.png",
  plot = p,
  width = 6,
  height = 5,
  dpi = 600
)

ggsave(
  "figures/GSE33147/PCA_GSE33147_Glycerol.pdf",
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
    grepl(
      "^WT_",
      colnames(expr_matrix)
    ),
    "WT",
    "Evolved"
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
  main = "Sample-to-sample distance of RMA-normalized log2 expression values (GSE33147_Glycerol)",
  filename = "figures/GSE33147/Sample_distance_heatmap_GSE33147_Glycerol.pdf",
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
  "figures/GSE33147/Sample_distance_heatmap_GSE33147_Glycerol.png",
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
  main = "Sample-to-sample distance of RMA-normalized log2 expression values (GSE33147_Glycerol)",
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
  condition = ifelse(
    grepl(
      "^WT_",
      colnames(expr_matrix)
    ),
    "WT",
    "Evolved"
  )
)

coldata$condition <- factor(
  coldata$condition,
  levels = c(
    "WT",
    "Evolved"
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
  coef = "conditionEvolved",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

res <- res[
  order(res$adj.P.Val),
]

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
  "results/GSE33147/All_DEGs_GSE33147_Glycerol.csv",
  row.names = FALSE
)

write.csv(
  deg,
  "results/GSE33147/Significant_DEGs_GSE33147_Glycerol.csv",
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
  "results/GSE33147/GSE33147_Glycerol_upregulated.csv",
  row.names = FALSE
)


#################################################
# Save downregulated genes
#################################################

write.csv(
  down_genes,
  "results/GSE33147/GSE33147_Glycerol_downregulated.csv",
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
  "results/GSE33147/GSE33147_Glycerol_SignedZ.csv",
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
  title = "Glycerol Evolved vs WT (GSE33147)",
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
    ylim = c(0, 15)
  )


p_volcano


#################################################
# Save volcano plot
#################################################

ggsave(
  "figures/GSE33147/Volcano_GSE33147_Glycerol.png",
  width = 8,
  height = 7,
  dpi = 600
)

ggsave(
  "figures/GSE33147/Volcano_GSE33147_Glycerol.pdf",
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
      "MA plot - Glycerol Evolved vs WT (GSE33147)",

    subtitle =
      "RMA-normalized Affymetrix expression data",

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
  "figures/GSE33147/MA_plot_GSE33147_Glycerol.png",
  plot = p_ma,
  width = 8,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/GSE33147/MA_plot_GSE33147_Glycerol.pdf",
  plot = p_ma,
  width = 8,
  height = 6
)

