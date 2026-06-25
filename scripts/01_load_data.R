# Script 1: Load and Explore GEO Data
# Disease Pathway Mapping - PCOD Analysis
# Purpose: Download GSE114419 and explore expression data
# Author: projectanishka
# Date: 2026-06-25

# ============================================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================================

cat("Loading packages...\n")

library(GEOquery)
library(dplyr)
library(readr)

cat("✓ Packages loaded\n\n")

# ============================================================================
# 2. DOWNLOAD GEO DATASET
# ============================================================================

cat("Downloading GSE114419 from GEO...\n")
cat("(This may take a few minutes on first run)\n\n")

# Increase timeout for large downloads
options(timeout = 300)

# Download the dataset with matrix files
gse <- getGEO("GSE114419", GSEMatrix = TRUE)

cat(sprintf("✓ Successfully downloaded %d dataset(s)\n\n", length(gse)))

# ============================================================================
# 3. EXTRACT EXPRESSION DATA
# ============================================================================

cat("Extracting expression matrix...\n")

# Get expression matrix (rows = genes, columns = samples)
expr <- exprs(gse[[1]])

cat(sprintf("✓ Expression matrix dimensions: %d genes × %d samples\n", 
            nrow(expr), ncol(expr)))
cat(sprintf("  - Genes (probes): %d\n", nrow(expr))
cat(sprintf("  - Samples: %d\n\n", ncol(expr)))

# Display first few rows and columns
cat("First 5 genes, first 5 samples:\n")
print(head(expr[, 1:5]))

# ============================================================================
# 4. EXPLORE SAMPLE METADATA
# ============================================================================

cat("\n\nExploring sample information...\n")

# Get phenotype data
pheno <- pData(gse[[1]])

# Display available columns
cat("Available sample metadata:\n")
print(colnames(pheno))

cat("\nFirst few samples:\n")
print(head(pheno))

# ============================================================================
# 5. IDENTIFY SAMPLE GROUPS
# ============================================================================

cat("\n\nAnalyzing sample characteristics...\n")

# Display unique values in key columns
for (col in colnames(pheno)) {
  unique_vals <- unique(pheno[[col]])
  if (length(unique_vals) <= 10) {
    cat(sprintf("\n%s (%d unique values):\n", col, length(unique_vals)))
    print(unique_vals)
  }
}

# ============================================================================
# 6. DATA QUALITY CHECKS
# ============================================================================

cat("\n\nData Quality Summary:\n")
cat("========================================\n")

# Check for missing values
missing_per_gene <- rowSums(is.na(expr))
missing_per_sample <- colSums(is.na(expr))

cat(sprintf("Missing values per gene - Min: %d, Max: %d, Mean: %.2f\n",
            min(missing_per_gene), max(missing_per_gene), mean(missing_per_gene)))
cat(sprintf("Missing values per sample - Min: %d, Max: %d, Mean: %.2f\n",
            min(missing_per_sample), max(missing_per_sample), mean(missing_per_sample)))

# Basic expression statistics
cat("\nExpression value statistics:\n")
cat(sprintf("  Min: %.2f\n", min(expr, na.rm = TRUE)))
cat(sprintf("  Max: %.2f\n", max(expr, na.rm = TRUE)))
cat(sprintf("  Mean: %.2f\n", mean(expr, na.rm = TRUE)))
cat(sprintf("  Median: %.2f\n", median(expr, na.rm = TRUE)))

# ============================================================================
# 7. SAVE DATA FOR DOWNSTREAM ANALYSIS
# ============================================================================

cat("\n\nSaving processed data...\n")

# Create data directory if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Save expression matrix
saveRDS(expr, "data/expression_matrix.rds")
write.csv(expr, "data/expression_matrix.csv")

# Save phenotype data
saveRDS(pheno, "data/phenotype_data.rds")
write.csv(pheno, "data/phenotype_data.csv")

cat("✓ Expression matrix saved to data/expression_matrix.rds/.csv\n")
cat("✓ Phenotype data saved to data/phenotype_data.rds/.csv\n\n")

# ============================================================================
# 8. SUMMARY
# ============================================================================

cat("========================================\n")
cat("LOADING COMPLETE\n")
cat("========================================\n\n")

cat("Next Steps:\n")
cat("1. Review the sample characteristics above\n")
cat("2. Run 02_deg_analysis.R to perform differential expression\n")
cat("3. Use View(pheno) to explore metadata in RStudio\n\n")

cat("Data saved in:\n")
cat("  - data/expression_matrix.rds (for R)\n")
cat("  - data/expression_matrix.csv (for spreadsheet/Python)\n")
cat("  - data/phenotype_data.rds (for R)\n")
cat("  - data/phenotype_data.csv (for spreadsheet/Python)\n")
