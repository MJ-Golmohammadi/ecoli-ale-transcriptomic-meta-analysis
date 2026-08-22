############################################################
# Robust Rank Aggregation (RRA) Meta-analysis
# Project:
# Cross-study transcriptomic meta-analysis of
# Adaptive Laboratory Evolution (ALE) in
# Escherichia coli K-12
############################################################

# Clear workspace
rm(list = ls())

library(dplyr)
library(readr)
library(RobustRankAggreg)

############################################################
# Load meta-analysis input files
############################################################

g1 <- read_csv("results/menF_contained_GSE224889-GSE122779_Combine_SignedZ.csv")
g2 <- read_csv("results/GSE140478/GSE140478_SignedZ.csv")
g3 <- read_csv("results/GSE158959/GSE158959_SignedZ.csv")
g4 <- read_csv("results/GSE17276/GSE17276_SignedZ.csv")
g5 <- read_csv("results/GSE206196/GSE206196_SignedZ.csv")
g6 <- read_csv("results/GSE33147/GSE33147_Glycerol_SignedZ.csv")
g7 <- read_csv("results/GSE316857/GSE316857_SignedZ_HPE.csv")
############################################################
# Remove genes with missing identifiers
############################################################

g1 <- g1 %>% filter(!is.na(Gene))
g2 <- g2 %>% filter(!is.na(Gene))
g3 <- g3 %>% filter(!is.na(Gene))
g4 <- g4 %>% filter(!is.na(Gene))
g5 <- g5 %>% filter(!is.na(Gene))
g6 <- g6 %>% filter(!is.na(Gene))
g7 <- g7 %>% filter(!is.na(Gene))
############################################################
# Generate ranked gene lists from Signed Z-scores
#
# Positive SignedZ:
# Conserved upregulated response
#
# Negative SignedZ:
# Conserved downregulated response
############################################################

rank_up <- function(df){
  df %>%
    filter(SignedZ > 0) %>%
    arrange(desc(SignedZ)) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    pull(Gene)
}

rank_down <- function(df){
  df %>%
    filter(SignedZ < 0) %>%
    arrange(SignedZ) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    pull(Gene)
}
############################################################
# Build ranked gene lists
############################################################

up_lists <- list(
  Combined = rank_up(g1),
  GSE140478 = rank_up(g2),
  GSE158959 = rank_up(g3),
  GSE17276 = rank_up(g4),
  GSE206196 = rank_up(g5),
  GSE33147 = rank_up(g6),
  GSE316857 = rank_up(g7)
)

down_lists <- list(
  Combined = rank_down(g1),
  GSE140478 = rank_down(g2),
  GSE158959 = rank_down(g3),
  GSE17276 = rank_down(g4),
  GSE206196 = rank_down(g5),
  GSE33147 = rank_down(g6),
  GSE316857 = rank_down(g7)
)
############################################################
# Remove invalid gene identifiers
############################################################

clean_rank_list <- function(x){
  x <- trimws(x)
  x <- x[!is.na(x)]
  x <- x[x != ""]
  x <- unique(x)
  return(x)
}

up_lists <- lapply(up_lists, clean_rank_list)
down_lists <- lapply(down_lists, clean_rank_list)
############################################################
# Define gene universe
############################################################

all_genes <- unique(
  unlist(
    c(
      up_lists,
      down_lists
    )
  )
)

Ngenes <- length(all_genes)

cat(
  "Gene universe size:",
  Ngenes,
  "\n"
)
############################################################
# Construct rank matrices
############################################################

rmat_up <- rankMatrix(
  up_lists,
  N = Ngenes,
  full = TRUE
)

rmat_down <- rankMatrix(
  down_lists,
  N = Ngenes,
  full = TRUE
)
############################################################
# Run Robust Rank Aggregation
############################################################

rra_up <- aggregateRanks(
  rmat = rmat_up,
  method = "RRA"
)

rra_down <- aggregateRanks(
  rmat = rmat_down,
  method = "RRA"
)
############################################################
# Adjust RRA scores using Benjamini-Hochberg correction
#
# FDR is reported but not used as the primary filtering
# criterion because RRA is performed across thousands
# of genes and BH correction becomes highly conservative.
############################################################

rra_up$FDR <- p.adjust(
  rra_up$Score,
  method = "BH"
)

rra_down$FDR <- p.adjust(
  rra_down$Score,
  method = "BH"
)
############################################################
# Inspect RRA results
############################################################

head(rra_up)

head(rra_down)
############################################################
# Select conserved meta-genes
#
# RRA Score < 0.01 was used to identify genes showing
# consistent ranking across independent datasets.
############################################################

core_up <- rra_up %>%
  filter(
    Score < 0.01
  )

core_down <- rra_down %>%
  filter(
    Score < 0.01
  )

cat(
  "RRA conserved upregulated genes:",
  nrow(core_up),
  "\n"
)

cat(
  "RRA conserved downregulated genes:",
  nrow(core_down),
  "\n"
)
############################################################
# Save RRA results
############################################################

dir.create(
  "results/RRA",
  showWarnings = FALSE,
  recursive = TRUE
)

write_csv(
  rra_up,
  "results/RRA/RRA_upregulated.csv"
)

write_csv(
  rra_down,
  "results/RRA/RRA_downregulated.csv"
)

############################################################
# Remove genes with inconsistent regulation direction
#
# Genes detected in both RRA directions are removed because
# they do not represent a consistent adaptive response.
############################################################

conflict_genes <- intersect(
  core_up$Name,
  core_down$Name
)

cat(
  "Number of conflicting genes:",
  length(conflict_genes),
  "\n"
)

print(conflict_genes)
############################################################
# Remove conflicting genes
############################################################

core_up <- core_up %>%
  filter(
    !Name %in% conflict_genes
  )

core_down <- core_down %>%
  filter(
    !Name %in% conflict_genes
  )
############################################################
# Merge upregulated and downregulated genes
############################################################

core_genes <- bind_rows(
  core_up %>%
    mutate(
      Direction = "Up"
    ),
  core_down %>%
    mutate(
      Direction = "Down"
    )
)
############################################################
# Verify duplicated genes
############################################################

cat(
  "Remaining duplicated genes:",
  sum(
    duplicated(core_genes$Name)
  ),
  "\n"
)
############################################################
# Count supporting studies
#
# A gene is counted once per independent dataset.
# The merged GSE122779-GSE224889 dataset is considered
# one independent comparison.
############################################################

gene_presence <- data.frame(
  Gene = core_genes$Name
)

study_lists <- list(
  Combined = c(
    g1$Gene
  ),
  GSE140478 = c(
    g2$Gene
  ),
  GSE158959 = c(
    g3$Gene
  ),
  GSE17276 = c(
    g4$Gene
  ),
  GSE206196 = c(
    g5$Gene
  ),
  GSE33147 = c(
    g6$Gene
  ),
  GSE316857 = c(
    g7$Gene
  )
)

gene_presence$Studies <- sapply(
  gene_presence$Gene,
  function(g){
    sum(
      sapply(
        study_lists,
        function(x){
          g %in% x
        }
      )
    )
  }
)
############################################################
# Keep genes supported by at least four independent studies
############################################################

core_genes <- core_genes %>%
  left_join(
    gene_presence,
    by = c(
      "Name" = "Gene"
    )
  ) %>%
  filter(
    Studies >= 6
  )
############################################################
# Update final upregulated and downregulated gene sets
############################################################

core_up <- core_genes %>%
  filter(
    Direction == "Up"
  )

core_down <- core_genes %>%
  filter(
    Direction == "Down"
  )


write_csv(
  core_up,
  "results/RRA/Core_upregulated_genes.csv"
)

write_csv(
  core_down,
  "results/RRA/Core_downregulated_genes.csv"
)

write_csv(
  core_genes,
  "results/RRA/Core_genes.csv"
)
############################################################
# Generate final gene list for downstream analyses
############################################################

genes_for_enrichment <- unique(
  core_genes$Name
)

cat(
  "Final number of core genes:",
  length(genes_for_enrichment),
  "\n"
)

table(core_genes$Direction)

table(core_genes$Studies)

nrow(core_up)

nrow(core_down)
sum(duplicated(core_genes$Name))

############################################################
# Export gene list for STRING
############################################################
write.table(
  genes_for_enrichment,
  "results/RRA/core_genes_STRING.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
############################################################
# Summary
############################################################

cat(
  "Core upregulated genes:",
  nrow(core_up),
  "\n"
)

cat(
  "Core downregulated genes:",
  nrow(core_down),
  "\n"
)

cat(
  "Total core genes:",
  nrow(core_genes),
  "\n"
)
############################################################
# Statistics
############################################################

head(core_genes, 20)

table(gene_presence$Studies)

nrow(core_up)

nrow(core_down)

nrow(core_genes)

head(rra_up)

head(rra_down)

length(genes_for_enrichment)

head(
  genes_for_enrichment,
  20
)

############################################################
## Figures
## E. coli K-12 ALE Meta-analysis
############################################################
############################################################
# UpSet plot showing gene occurrence across datasets
############################################################

library(dplyr)

library(tidyr)

library(UpSetR)
############################################################
# Create gene universe from all studies
############################################################

all_genes <- unique(
  c(
    g1$Gene,
    g2$Gene,
    g3$Gene,
    g4$Gene,
    g5$Gene,
    g6$Gene,
    g7$Gene
  )
)
############################################################
# Build presence/absence matrix
############################################################

upset_df <- data.frame(
  Gene = all_genes
)

upset_df$Combined <- ifelse(
  all_genes %in% g1$Gene,
  1,
  0
)

upset_df$GSE140478 <- ifelse(
  all_genes %in% g2$Gene,
  1,
  0
)

upset_df$GSE158959 <- ifelse(
  all_genes %in% g3$Gene,
  1,
  0
)

upset_df$GSE17276 <- ifelse(
  all_genes %in% g4$Gene,
  1,
  0
)

upset_df$GSE206196 <- ifelse(
  all_genes %in% g5$Gene,
  1,
  0
)

upset_df$GSE33147 <- ifelse(
  all_genes %in% g6$Gene,
  1,
  0
)

upset_df$GSE316857 <- ifelse(
  all_genes %in% g7$Gene,
  1,
  0
)
############################################################
# Save UpSet plot
############################################################

dir.create(
  "figures/meta_analysis",
  showWarnings = FALSE,
  recursive = TRUE
)

svg(
  "figures/meta_analysis/Figure1_UpSet_StudyOverlap.svg",
  width = 11,
  height = 7.5
)

upset(
  upset_df,
  sets = c(
    "Combined",
    "GSE140478",
    "GSE158959",
    "GSE17276",
    "GSE206196",
    "GSE33147",
    "GSE316857"
  ),
  keep.order = TRUE,
  order.by = "freq",
  nintersects = 20,
  mb.ratio = c(0.65,0.35),
  text.scale = 1.5,
  point.size = 4,
  line.size = 1.2
)

dev.off()
############################################################
# Save PDF version
############################################################

pdf(
  "figures/meta_analysis/Figure1_UpSet_StudyOverlap.pdf",
  width = 11,
  height = 7.5,
  onefile = FALSE
)

upset(
  upset_df,
  sets = c(
    "Combined",
    "GSE140478",
    "GSE158959",
    "GSE17276",
    "GSE206196",
    "GSE33147",
    "GSE316857"
  ),
  keep.order = TRUE,
  order.by = "freq",
  mb.ratio = c(0.65,0.35),
  text.scale = 1.5,
  nintersects = 20,
  point.size = 4,
  line.size = 1.2
)

dev.off()
############################################################
# Top 20 Upregulated Genes Identified by RRA
############################################################

library(ggplot2)

library(dplyr)

top20_up <- core_up %>%
  arrange(Score) %>%
  head(20) %>%
  mutate(
    MinusLog10Score = -log10(Score)
  )

p_up <- ggplot(
  top20_up,
  aes(
    x = reorder(Name, MinusLog10Score),
    y = MinusLog10Score
  )
) +
  geom_col(
    fill = "#B2182B",
    width = 0.8
  ) +
  coord_flip() +
  theme_bw(
    base_size = 14
  ) +
  labs(
    title = "Top 20 Conserved Upregulated Genes",
    x = "",
    y = expression(-log[10](RRA~Score))
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/meta_analysis/Figure2A_Top20_Upregulated_RRA.svg",
  p_up,
  width = 8,
  height = 6
)
ggsave(
  "figures/meta_analysis/Figure2A_Top20_Upregulated_RRA.pdf",
  p_up,
  width = 8,
  height = 6
)
############################################################
# Top 20 Downregulated Genes Identified by RRA
############################################################

top20_down <- core_down %>%
  arrange(Score) %>%
  head(20) %>%
  mutate(
    MinusLog10Score = -log10(Score)
  )

p_down <- ggplot(
  top20_down,
  aes(
    x = reorder(Name, MinusLog10Score),
    y = MinusLog10Score
  )
) +
  geom_col(
    fill = "#2166AC",
    width = 0.8
  ) +
  coord_flip() +
  theme_bw(
    base_size = 14
  ) +
  labs(
    title = "Top 20 Conserved Downregulated Genes",
    x = "",
    y = expression(-log[10](RRA~Score))
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/meta_analysis/Figure2B_Top20_Downregulated_RRA.svg",
  p_down,
  width = 8,
  height = 6
)
ggsave(
  "figures/meta_analysis/Figure2B_Top20_Downregulated_RRA.pdf",
  p_down,
  width = 8,
  height = 6
)


##############################################
# test results
#############################################

core_genes %>%
  select(Name, Direction, Studies, Score) %>%
  arrange(Direction, Score)



datasets <- list(
  Combined = g1,
  GSE140478 = g2,
  GSE158959 = g3,
  GSE17276 = g4,
  GSE206196 = g5,
  GSE33147 = g6,
  GSE316857 = g7
)



# analyse the rank of a specific gene
flu_rank <- lapply(
  names(datasets),
  function(x){

    datasets[[x]] %>%
      filter(SignedZ < 0) %>%
      arrange(SignedZ) %>%
      mutate(Rank = row_number()) %>%
      filter(Gene == "flu") %>%
      mutate(Dataset = x) %>%
      select(Dataset, Gene, SignedZ, Rank)

  }
) %>%
  bind_rows()

rra_rank <- rra_down %>%
  arrange(Score) %>%
  mutate(
    RRA_Rank = row_number()
  ) %>%
  filter(Name == "flu") %>%
  select(
    Gene = Name,
    RRA_Rank,
    Score
  )

flu_rank <- flu_rank %>%
  left_join(
    rra_rank,
    by = "Gene"
  )

flu_rank