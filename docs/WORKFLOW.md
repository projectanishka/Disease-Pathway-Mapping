# PCOD Pathway Mapping - Detailed Workflow

## Phase 1: Data Acquisition & Setup

### Step 1.1: Identify the Dataset
- **GEO Accession:** GSE114419
- **Title:** Coding and Non-coding gene expression signatures of granulosa cells from Polycystic Ovary Syndrome Patients
- **Link:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114419

### Step 1.2: Choose the Right Data Format

⚠️ **Important Lesson Learned:**
- **CEL files** = Raw, unnormalized data. These are individual cell intensity files that require normalization preprocessing. Not suitable for quick analysis.
- **Processed Expression Matrix (TXT)** = Pre-normalized, ready-to-use expression values. **This is what we use.**

✅ Use the **processed expression matrix** (normalized TXT format) available on GEO2R.

### Step 1.3: R Environment Configuration

```r
# Verify R version
R.version.string

# Check for Rtools (Windows users)
Sys.which("make")

# Install required packages
install.packages(c("dplyr", "ggplot2", "readr", "tidyr", "pheatmap", "jsonlite"))

# Install Bioconductor packages
install.packages("BiocManager")
BiocManager::install("GEOquery")
```

**Verification Test:**
```r
library(jsonlite)
library(dplyr)

test <- data.frame(
  gene = c("INSR", "AKT1", "IRS1"),
  logFC = c(2.1, 1.8, -1.2)
)
print(test)
```

## Phase 2: Data Loading & Exploration

### Step 2.1: Load GEO Dataset

```r
library(GEOquery)

# Download the dataset (first time only, takes time)
gse <- getGEO("GSE114419", GSEMatrix = TRUE)

# Check how many datasets were returned
length(gse)
```

### Step 2.2: Extract Expression Matrix

```r
# Get the expression data
expr <- exprs(gse[[1]])

# Check dimensions: genes × samples
dim(expr)

# View first few genes and samples
head(expr[,1:5])
```

**Expected Output:**
- **Rows:** Gene probes/IDs (typically thousands)
- **Columns:** Samples (healthy and PCOD)
- **Values:** Log2-transformed normalized expression

### Step 2.3: Explore Sample Metadata

```r
# Extract phenotype (sample) information
pheno <- pData(gse[[1]])

# View the metadata
View(pheno)

# Key information to check:
# - Which samples are "healthy" vs "PCOD"
# - Sample identifiers
# - Any other treatment/grouping variables
```

## Phase 3: Differential Expression Analysis

### Step 3.1: Prepare Groups

```r
library(dplyr)

# Create a grouping variable based on sample characteristics
# Example: if "disease status" is in a column called "characteristics_ch1"

pheno <- pheno %>%
  mutate(group = case_when(
    grepl("PCOD|polycystic|diseased", characteristics_ch1, ignore.case = TRUE) ~ "PCOD",
    grepl("control|healthy|normal", characteristics_ch1, ignore.case = TRUE) ~ "Healthy",
    TRUE ~ "Unknown"
  ))

# Verify groups
table(pheno$group)
```

### Step 3.2: Perform Differential Expression Analysis

```r
# Option A: Using limma (recommended for microarray)
BiocManager::install("limma")
library(limma)

# Create design matrix
design <- model.matrix(~ 0 + pheno$group)
colnames(design) <- c("Healthy", "PCOD")

# Fit linear model
fit <- lmFit(expr, design)

# Define contrast: PCOD vs Healthy
contrast.matrix <- makeContrasts(PCOD - Healthy, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get results
deg_results <- topTable(fit2, adjust.method = "BH", number = Inf)
```

### Step 3.3: Select Top 100 DEGs

```r
# Sort by p-value and logFC
deg_top100 <- deg_results %>%
  arrange(adj.P.Val) %>%
  head(100)

# Add gene names (if available)
deg_top100$gene <- rownames(deg_top100)

# Display top 10
head(deg_top100, 10)

# Save results
write.csv(deg_top100, "results/top_100_degs.csv", row.names = FALSE)
```

## Phase 4: Pathway Enrichment Analysis

### Step 4.1: Prepare Gene List

```r
# Get gene IDs for pathway tools
gene_list <- deg_top100$gene

# Remove any duplicates or NA values
gene_list <- na.omit(unique(gene_list))

# Print for manual queries or save
writeLines(gene_list, "results/top_100_genes.txt")
```

### Step 4.2: Query g:Profiler

```r
# Install gprofiler2
install.packages("devtools")
devtools::install_github("cran/gProfileR")

library(gProfileR)

# Functional enrichment analysis
gp_results <- gprofiler(
  query = gene_list,
  organism = "hsapiens",  # Human
  src_filter = c("GO", "KEGG", "REAC"),  # Include GO, KEGG, Reactome
  ordered_query = TRUE
)

# View results
View(gp_results)

# Save results
write.csv(gp_results, "results/gprofiler_enrichment.csv", row.names = FALSE)
```

### Step 4.3: KEGG Pathway Analysis

```r
# For direct KEGG queries (requires additional setup)
# Option: Use online KEGG tool or dedicated R packages

# Install clusterProfiler (more comprehensive)
BiocManager::install("clusterProfiler")
library(clusterProfiler)

# Note: May require KEGG API key for large queries
```

### Step 4.4: Reactome Pathway Analysis

```r
# Similar to KEGG, Reactome can be queried via:
# 1. Web interface: https://reactome.org/
# 2. R packages: ReactomePA

BiocManager::install("ReactomePA")
library(ReactomePA)

# Direct API queries or web-based analysis
```

## Phase 5: Visualization

### Step 5.1: Volcano Plot

```r
library(ggplot2)

ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < 0.05 & abs(logFC) > 1)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  labs(title = "Volcano Plot: PCOD vs Healthy",
       x = "log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme_minimal()

ggsave("results/volcano_plot.png", width = 8, height = 6)
```

### Step 5.2: Heatmap of Top DEGs

```r
library(pheatmap)

# Extract expression for top 100 genes
heatmap_data <- expr[rownames(deg_top100), ]

# Standardize for better visualization
heatmap_data <- scale(t(heatmap_data))

# Create heatmap
pheatmap(
  t(heatmap_data),
  main = "Top 100 DEGs: PCOD vs Healthy",
  annotation_col = pheno[, c("group")],
  show_colnames = FALSE,
  filename = "results/heatmap_top100.png",
  width = 12,
  height = 10
)
```

### Step 5.3: Pathway Enrichment Visualization

```r
# Barplot of top pathways
top_pathways <- gp_results %>%
  top_n(15, -p.value)

ggplot(top_pathways, aes(x = reorder(term.name, -p.value), y = -log10(p.value))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 15 Enriched Pathways",
       x = "Pathway",
       y = "-log10 P-value") +
  theme_minimal()

ggsave("results/pathway_enrichment.png", width = 10, height = 8)
```

## Key Output Files

After completing this workflow, you should have:

```
results/
├── top_100_degs.csv              # Differentially expressed genes
├── top_100_genes.txt             # Gene list for pathway tools
├── gprofiler_enrichment.csv      # g:Profiler results
├── kegg_enrichment.csv           # KEGG pathway results
├── reactome_enrichment.csv       # Reactome results
├── volcano_plot.png              # DE visualization
├── heatmap_top100.png            # Expression heatmap
└── pathway_enrichment.png        # Pathway bar plot
```

## Interpretation Guide

### What the results mean:

1. **Top 100 DEGs:** Genes most different between healthy and PCOD granulosa cells
2. **logFC > 0:** Gene is upregulated in PCOD (more expression)
3. **logFC < 0:** Gene is downregulated in PCOD (less expression)
4. **Enriched Pathways:** Biological processes disrupted in PCOD

### Example Interpretation:
If "Insulin Signaling" is an enriched pathway with multiple upregulated genes, it suggests PCOD involves dysregulated glucose/insulin metabolism in ovarian cells—consistent with PCOD's known insulin resistance.

---

**Status:** In Progress - Phase 3 (DEG Analysis) currently being implemented.
