############################################################
# Figure 4: Hub Gene Network
# Publication-quality Top20 Degree Hub Gene Network
# Pseudomonas aeruginosa meta-analysis
############################################################

library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)

theme_set(theme_bw())

############################################################
# Import STRING network
############################################################

network <- read.delim(
  "results/String Files/network.tsv",
  header = TRUE,
  check.names = FALSE
)

############################################################
# Import Degree hub ranking
############################################################

hub <- read.csv(
  "results/String Files/Hubba_Degree_Top20_Ranking.csv",
  skip = 1,
  header = TRUE,
  check.names = FALSE
)

hub <- hub %>%
  rename(
    Gene = Name,
    Degree = Score
  )

hub_genes <- hub$Gene

############################################################
# Import RRA results
############################################################

rra <- read.csv(
  "results/RRA/Core_genes.csv",
  stringsAsFactors = FALSE
)

rra <- rra %>%
  select(
    Name,
    Score,
    FDR,
    Direction,
    Studies
  ) %>%
  rename(
    Gene = Name,
    RRA_score = Score
  )

############################################################
# Prepare STRING network
############################################################

network <- network %>%
  rename(
    protein1 = `#node1`,
    protein2 = node2
  )
colnames(network)
############################################################
# Keep only top hub network
############################################################

network20 <- network %>%
  filter(
    protein1 %in% hub_genes &
    protein2 %in% hub_genes
  )

############################################################
# Create functional annotation
############################################################

annotation <- tibble(
  Gene = hub_genes
)

annotation$Class <- "Other"
############################################################
# Flagellar assembly
############################################################

annotation$Class[
  annotation$Gene %in%
  c(
    "fliC",
    "fliD",
    "fliA",
    "fliF",
    "fliK",
    "fliS",
    "flgB",
    "flgD",
    "flgE",
    "flgG",
    "flgH",
    "flgJ",
    "flgK",
    "flgI"
  )
] <- "Flagellar assembly"


############################################################
# Flagellar motor
############################################################

annotation$Class[
  annotation$Gene %in%
  c(
    "motA",
    "motB"
  )
] <- "Flagellar motor"


############################################################
# Chemotaxis
############################################################

annotation$Class[
  annotation$Gene %in%
  c(
    "cheA",
    "cheB",
    "cheW",
    "cheY",
    "cheZ",
    "tar",
    "tsr"
  )
] <- "Chemotaxis"


############################################################
# Curli / Biofilm
############################################################

annotation$Class[
  annotation$Gene %in%
  c(
    "csgA",
    "csgD",
    "csgE",
    "csgF",
    "csgG",
    "flu"
  )
] <- "Biofilm / Curli"

############################################################
# Regulatory
############################################################
annotation$Class[
  annotation$Gene == "evgS"
] <- "Regulatory"

############################################################
# Merge Degree and RRA information
############################################################

annotation <- annotation %>%
  left_join(
    hub,
    by = "Gene"
  ) %>%
  left_join(
    rra,
    by = "Gene"
  )

############################################################
# Create node size variable
############################################################

annotation$Degree_plot <- annotation$Degree

############################################################
# Build igraph object
############################################################

g <- graph_from_data_frame(
  d = network20,
  vertices = annotation,
  directed = FALSE
)

############################################################
# Convert to tidygraph
############################################################

tg <- as_tbl_graph(g)

############################################################
# Colors
############################################################

my_colors <- c(
  "Flagellar assembly" = "#D73027",
  "Flagellar motor"    = "#E6AB02",
  "Chemotaxis"         = "#4575B4",
  "Biofilm / Curli"    = "#e5bbec",
  "Regulatory"         = "#b950b6",
  "Other"              = "#938d79"
)

############################################################
# Plot hub network
############################################################

fig4 <- ggraph(
  tg,
  layout = "kk"
) +

  geom_edge_link(
    aes(
      width = combined_score,
      alpha = combined_score
    ),
    colour = "grey70"
  ) +

  scale_edge_width(
    range = c(
      0.05,
      0.7
    )
  ) +

  scale_edge_alpha(
    range = c(
      0.2,
      0.8
    )
  ) +

  geom_node_point(
    aes(
      size = Degree_plot,
      fill = Class
    ),
    shape = 21,
    colour = "black",
    stroke = 0.3
  ) +

  scale_size_continuous(
    range = c(7,24),
    guide = guide_legend(
      override.aes = list(
        shape = 21,
        fill = "grey70",
        colour = "black"
      )
    )
  ) +

  guides(

    fill = guide_legend(

      override.aes = list(

        shape = 21,
        size = 14,
        colour = "black"

      )

    )

  ) +

  scale_fill_manual(
    values = my_colors
  ) +

  geom_node_text(
    aes(label = name),
    size = 5,
    fontface = "bold",
    repel = TRUE,
    max.overlaps = Inf,
    box.padding = 1,
    point.padding = 0.5,
    force = 10
) +

  theme_void() +

  theme(
    legend.position = "right",

    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5
    ),

    legend.title = element_text(
      face = "bold",
      size = 16
    ),

    legend.text = element_text(
      size = 14
    )
  ) +

  labs(
    title = "Top Degree Hub Gene Network",
    fill = "Functional category",
    size = "Degree"
  )
############################################################
# Display figure
############################################################

fig4

############################################################
# Save publication-quality figure
############################################################

ggsave(
  "figures/string/Figure4_HubNetwork.pdf",
  fig4,
  width = 12,
  height = 10
)

ggsave(
  "figures/string/Figure4_HubNetwork.svg",
  fig4,
  width = 12,
  height = 10
)







############################################################
# Hub gene heatmap across ALE datasets
# E. coli K-12 transcriptomic meta-analysis
############################################################

library(tidyverse)
library(pheatmap)
library(svglite)

############################################################
# Read dataset files
############################################################

files <- c(
  GSE224889_menF = "results/GSE224889/GSE224889_menF_SignedZ.csv",
  GSE224889_menF_ubiC = "results/GSE224889/GSE224889_menF-ubiC_SignedZ.csv",
  GSE122779 = "results/GSE122779/GSE122779_SignedZ.csv",
  GSE140478 = "results/GSE140478/GSE140478_SignedZ.csv",
  GSE158959 = "results/GSE158959/GSE158959_SignedZ.csv",
  GSE17276 = "results/GSE17276/GSE17276_SignedZ.csv",
  GSE206196 = "results/GSE206196/GSE206196_SignedZ.csv",
  GSE33147 = "results/GSE33147/GSE33147_Glycerol_SignedZ.csv",
  GSE316857 = "results/GSE316857/GSE316857_SignedZ_HPE.csv"
)

############################################################
# Extract hub genes
############################################################

hub_genes <- hub$Gene

############################################################
# Read logFC values from each dataset
############################################################

expr_list <- list()

for(i in seq_along(files)){
  tmp <- read.csv(files[i])
  tmp <- tmp %>%
    select(
      Gene,
      logFC
    ) %>%
    rename(
      !!names(files)[i] := logFC
    )
  expr_list[[i]] <- tmp
}

############################################################
# Merge logFC values
############################################################

expr <- reduce(
  expr_list,
  full_join,
  by = "Gene"
)

############################################################
# Keep only hub genes
############################################################

expr <- expr %>%
  filter(
    Gene %in% hub_genes
  )

############################################################
# Convert to matrix
############################################################

expr_mat <- expr %>%
  column_to_rownames("Gene") %>%
  as.matrix()

############################################################
# Replace missing values
############################################################

expr_mat[is.na(expr_mat)] <- 0

############################################################
# Row Z-score normalization for visualization
############################################################

expr_z <- t(
  scale(
    t(expr_mat)
  )
)

############################################################
# Create row annotation
############################################################

annotation_heat <- annotation %>%
  select(
    Gene,
    Class
  )

annotation_heat <- as.data.frame(annotation_heat)

rownames(annotation_heat) <- annotation_heat$Gene

annotation_heat$Gene <- NULL

annotation_heat <- annotation_heat[
  rownames(expr_z),
  ,
  drop = FALSE
]

############################################################
# Functional category colors
############################################################
ann_colors <- list(
  Class = c(
    "Flagellar assembly" = "#0072B2",
    "Chemotaxis"         = "#E69F00",
    "Flagellar motor"    = "#009E73",
    "Biofilm / Curli"    = "#f3b4d6",
    "Regulatory"         = "#b85db2",
    "Other"              = "#7F7F7F"
  )
)

############################################################
# Draw heatmap
############################################################

heatmap_plot <- pheatmap(
  expr_z,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 12,
  fontsize_col = 12,
  border_color = NA,
  annotation_row = annotation_heat,
  annotation_colors = ann_colors,
  color = colorRampPalette(
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )(100),
  main = "Hub Gene Expression Patterns Across ALE Datasets"
)

############################################################
# Save SVG
############################################################

svglite(
  "figures/string/Figure6_HubGene_Heatmap.svg",
  width = 9,
  height = 12
)

pheatmap(
  expr_z,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 12,
  fontsize_col = 12,
  border_color = NA,
  annotation_row = annotation_heat,
  annotation_colors = ann_colors,
  color = colorRampPalette(
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )(100),
  main = "Hub Gene Expression Patterns Across ALE Datasets"
)

dev.off()

############################################################
# Save PDF
############################################################

pdf(
  "figures/string/Figure6_HubGene_Heatmap.pdf",
  width = 9,
  height = 12
)

pheatmap(
  expr_z,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 12,
  fontsize_col = 12,
  border_color = NA,
  annotation_row = annotation_heat,
  annotation_colors = ann_colors,
  color = colorRampPalette(
    c(
      "#2166AC",
      "white",
      "#B2182B"
    )
  )(100),
  main = "Hub Gene Expression Patterns Across ALE Datasets"
)

dev.off()

