#################################################
# Combine two Signed Z-score files
# Using Stouffer's method
#################################################
# Clear the workspace
rm(list = ls())

library(readxl)
library(readr)
#################################################
# Load files
#################################################

# File 1
deg1 <- read.csv(
  "results/GSE224889/GSE224889_menF_SignedZ.csv",
  stringsAsFactors = FALSE
)

# File 2
deg2 <- read.csv(
  "results/GSE224889/GSE224889_menF-ubiC_SignedZ.csv",
  stringsAsFactors = FALSE
)

# File 3
deg3 <- read.csv(
  "results/GSE122779/GSE122779_SignedZ.csv",
  stringsAsFactors = FALSE
)

#################################################
# Select required columns
#################################################

# We only need:
# Gene ID
# Signed Z-score

z1 <- deg1[, c(
  "Gene",
  "SignedZ"
)]

z2 <- deg2[, c(
  "Gene",
  "SignedZ"
)]

z3 <- deg3[, c(
  "Gene",
  "SignedZ"
)]
#################################################
# Rename columns
#################################################

# To distinguish the two signatures

colnames(z1)[2] <- "SignedZ_1"
colnames(z2)[2] <- "SignedZ_2"
colnames(z3)[2] <- "SignedZ_3"

#################################################
# Number of common genes among three comparisons
#################################################

cat("z1:", length(unique(z1$Gene)), "\n")
cat("z2:", length(unique(z2$Gene)), "\n")
cat("z3:", length(unique(z3$Gene)), "\n")

common_genes <- Reduce(intersect, list(z1$Gene, z2$Gene, z3$Gene))

cat("Common:", length(common_genes), "\n")
#################################################
# Merge three signatures by gene
#################################################

# Only genes present in 3 comparisons
# will be combined

merged <- merge(
  z1,
  z2,
  by = "Gene",
  all = FALSE
)

merged <- merge(
  merged,
  z3,
  by = "Gene",
  all = FALSE
)

#################################################
# Weighted Stouffer method
# Weight based on number of genes in each file
#################################################


#################################################
# Count number of genes in each DEG file
#################################################

n1 <- nrow(deg1)
n2 <- nrow(deg2)
n3 <- nrow(deg3)

#################################################
# Calculate weights
# Weight = sqrt(number of genes)
#################################################

w1 <- sqrt(n1)
w2 <- sqrt(n2)
w3 <- sqrt(n3)

#################################################
# Combine Signed Z scores
#################################################

merged$SignedZ <- (
  w1 * merged$SignedZ_1 +
  w2 * merged$SignedZ_2 +
  w3 * merged$SignedZ_3
) / sqrt(
  w1^2 + w2^2 +w3^2 
)


#################################################
# Check
#################################################

cat("Number of genes in file 1:", n1, "\n")
cat("Number of genes in file 2:", n2, "\n")
cat("Number of genes in file 2:", n3, "\n")

cat("Weights:", w1, w2, w3, "\n")

#################################################
# Check results
#################################################

head(merged, 20)

#################################################
# Prepare final file for RRA
#################################################

# Keep only:
# Gene
# Combined Signed Z

final_z <- merged[, c(
  "Gene",
  "SignedZ"
)]

final_z <- final_z[
    order(-abs(final_z$SignedZ)),
]
#################################################
# Save output
#################################################

write.csv(
  final_z,
  "results/menF_contained_GSE224889-GSE122779_Combine_SignedZ.csv",
  row.names = FALSE
)
