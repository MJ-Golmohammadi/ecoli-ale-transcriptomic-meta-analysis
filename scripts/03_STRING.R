############################################################
# Load packages
############################################################

# Clear workspace
rm(list = ls())

library(readr)
library(dplyr)
library(rtracklayer)

############################################################
# Read input files
############################################################

file1 <- read.csv(
  "results/RRA/Core_genes.csv",
  stringsAsFactors = FALSE
)

############################################################
# Import E. coli K-12 MG1655 GFF annotation
############################################################

gff <- import(
  "data/Ecoli_K12_MG1655_annotation.gff"
)

############################################################
# Extract gene annotation table
############################################################

annotation <- data.frame(
  Gene_Name = mcols(gff)$gene,
  Locus_Tag = mcols(gff)$locus_tag,
  stringsAsFactors = FALSE
)

############################################################
# Remove invalid annotations
############################################################

annotation <- annotation %>%
  filter(
    !is.na(Gene_Name),
    !is.na(Locus_Tag)
  ) %>%
  distinct()

############################################################
# Convert gene names to locus tags
############################################################

file1$string_format <- sapply(
  file1$Name,
  function(gene){

    hits <- annotation %>%
      filter(
        Gene_Name == gene
      )

    # No match found
    if(nrow(hits) == 0){

      return(gene)

    }

    # Multiple matches found
    if(nrow(hits) > 1){

      return("double")

    }

    # Single match found
    return(hits$Locus_Tag[1])

  }
)

############################################################
# Save output
############################################################

write.csv(
  file1,
  "results/string_format_Ecoli_locus_tags.csv",
  row.names = FALSE
)

############################################################
# Summary
############################################################

cat(
  "No match:",
  sum(file1$string_format == file1$Name),
  "\n"
)

cat(
  "Multiple matches:",
  sum(file1$string_format == "double"),
  "\n"
)

cat(
  "Converted to locus tags:",
  sum(
    file1$string_format != file1$Name &
    file1$string_format != "double"
  ),
  "\n"
)









############################################################
# Load packages
############################################################

library(readr)
library(dplyr)

############################################################
# Read STRING enrichment files
############################################################

go_bp <- read_tsv(
  "results/String Files/GO_BP.tsv",
  show_col_types = FALSE
)

kegg <- read_tsv(
  "results/String Files/KEGG.tsv",
  show_col_types = FALSE
)


head(go_bp)
head(kegg)


head(go_bp,20)

head(kegg,10)

head(hub,50)


############################################################
# Load packages
############################################################
library(dplyr)
library(ggplot2)
############################################################
# Prepare GO data
############################################################
go_plot <- dplyr::slice(
  arrange(
    go_bp,
    `false discovery rate`
  ),
  1:10
) %>%
  mutate(
    Term = factor(
      `term description`,
      levels = rev(`term description`)
    ),
    FDR = -log10(`false discovery rate`)
  )
############################################################
# Plot GO enrichment
############################################################

ggplot(
  go_plot,
  aes(
    x = FDR,
    y = Term,
    size = `observed gene count`,
    color = strength
  )
) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(
    low = "#74add1",
    high = "#d73027"
  ) +
  labs(
    title = "GO Biological Process Enrichment",
    x = expression(-log[10](FDR)),
    y = NULL,
    color = "Strength",
    size = "Gene count"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )


ggsave(
  "figures/string/Figure3A_GO_BP.svg",
  width = 8,
  height = 6
)

ggsave(
  "figures/string/Figure3A_GO_BP.pdf",
  width = 8,
  height = 6
)


#############################################################
# Prepare KEGG data
############################################################

kegg_plot <- kegg %>%
  mutate(
    Term = factor(
      `term description`,
      levels = rev(`term description`)
    ),
    FDR = -log10(`false discovery rate`)
  )

############################################################
# Plot KEGG enrichment
############################################################

ggplot(
  kegg_plot,
  aes(
    x = FDR,
    y = Term,
    size = `observed gene count`,
    color = strength
  )
) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(
    low = "#4575b4",
    high = "#d73027"
  ) +
  labs(
    title = "KEGG Pathway Enrichment",
    x = expression(-log[10](FDR)),
    y = NULL,
    color = "Strength",
    size = "Gene count"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )


ggsave(
  "figures/string/Figure3B_KEGG.svg",
  width = 8,
  height = 6
)

ggsave(
  "figures/string/Figure3B_KEGG.pdf",
  width = 8,
  height = 6
)
