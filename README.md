# Disease Pathway Mapping in PCOD

## Overview

This project identifies disease pathways in Polycystic Ovary Syndrome (PCOD) using gene expression data from the Gene Expression Omnibus (GEO). By analyzing differentially expressed genes (DEGs) from healthy vs. PCOD-affected granulosa cells, we map the biological pathways dysregulated in this condition.

**Key Insight:** PCOD is a polygenic disorder—not caused by a single gene, but by the dysregulation of multiple genes working together. This analysis identifies the top 100 most significant DEGs and maps them to known biological pathways.

## Dataset

- **Source:** [GEO2R - GSE114419](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114419)
- **Title:** Coding and Non-coding gene expression signatures of granulosa cells from Polycystic Ovary Syndrome Patients
- **Comparison:** Healthy vs. PCOD-affected granulosa cells
- **Format:** Processed expression matrix (normalized values, TXT format)

## Workflow

```
1. Obtain Dataset (GEO2R)
   ↓
2. Load & Explore Expression Data in R
   ↓
3. Perform Differential Expression Analysis
   ↓
4. Select Top 100 DEGs (by logFC and p-value)
   ↓
5. Pathway Enrichment Analysis
   ├─ g:Profiler
   ├─ KEGG
   └─ Reactome
   ↓
6. Visualize Results
```

## R Environment Setup

### Prerequisites
- **R Version:** 4.6.0+
- **Rtools:** Installed (for Windows users)
- **RStudio:** Recommended

### Installed Packages

```r
# Base packages
library(GEOquery)  # Access GEO datasets
library(dplyr)     # Data manipulation
library(readr)     # Data import
library(tidyr)     # Data tidying

# Visualization
library(ggplot2)   # Plots
library(pheatmap)  # Heatmaps

# Pathway analysis (to be installed)
library(jsonlite)  # JSON parsing for APIs
```

### Installation Commands

```r
# Install CRAN packages
install.packages(c("dplyr", "ggplot2", "readr", "tidyr", "pheatmap", "jsonlite"))

# Install Bioconductor packages
install.packages("BiocManager")
BiocManager::install("GEOquery")
```

## Quick Start

### 1. Load the Dataset

```r
library(GEOquery)

# Download GEO dataset
gse <- getGEO("GSE114419", GSEMatrix = TRUE)

# Extract expression matrix
expr <- exprs(gse[[1]])

# Check dimensions
dim(expr)

# View first few genes
head(expr[,1:5])
```

### 2. Explore Sample Metadata

```r
# Get phenotype data
pheno <- pData(gse[[1]])

# View sample information
View(pheno)
```

### 3. Next Steps

- [ ] Perform differential expression analysis (limma or edgeR)
- [ ] Select top 100 DEGs
- [ ] Query pathway databases (g:Profiler, KEGG, Reactome)
- [ ] Generate visualizations (volcano plots, heatmaps, pathway diagrams)

## Project Structure

```
Disease-Pathway-Mapping/
├── README.md                 # Project overview
├── docs/
│   ├── WORKFLOW.md          # Detailed workflow
│   ├── DATA_SOURCES.md      # GEO dataset info
│   └── SETUP.md             # R environment setup
├── scripts/
│   ├── 01_load_data.R       # Load & explore GEO data
│   ├── 02_deg_analysis.R    # Differential expression
│   ├── 03_pathway_analysis.R # Enrichment analysis
│   └── 04_visualize.R       # Create plots & heatmaps
├── data/
│   └── (datasets will be downloaded via R)
├── results/
│   ├── deg_results.csv      # Top 100 DEGs
│   ├── pathway_results.csv  # Pathway enrichment
│   └── plots/               # Visualizations
└── requirements.R           # R packages needed
```

## Key Concepts

### DEG (Differentially Expressed Genes)
Genes that show statistically significant differences in expression between healthy and PCOD samples.

**Metrics:**
- **logFC (log Fold Change):** How much a gene's expression changed (positive = upregulated, negative = downregulated)
- **p-value:** Statistical significance of the change

### Pathway Enrichment
Once we have the top 100 DEGs, we'll use:
- **g:Profiler** - Gene ontology and pathway annotation
- **KEGG** - Kyoto Encyclopedia of Genes and Genomes
- **Reactome** - Detailed molecular interaction database

This tells us *what biological processes* are disrupted in PCOD.

## Status

🚧 **In Progress**
- [x] R environment setup
- [x] GEO dataset identified
- [x] Basic data loading working
- [ ] DEG analysis pipeline
- [ ] Pathway enrichment integration
- [ ] Visualization suite

## References

- GEO Database: https://www.ncbi.nlm.nih.gov/geo/
- g:Profiler: https://biit.cs.ut.ee/gprofiler/
- KEGG Pathway Database: https://www.kegg.jp/
- Reactome: https://reactome.org/

## Author

projectanishka

---

**Note:** This is an educational bioinformatics project exploring PCOD disease mechanisms through gene expression analysis.
