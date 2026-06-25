# Data Sources & Background

## GEO Dataset: GSE114419

### Basic Information
- **Accession:** GSE114419
- **Title:** Coding and Non-coding gene expression signatures of granulosa cells from Polycystic Ovary Syndrome Patients
- **Link:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114419
- **Organism:** Homo sapiens (Human)
- **Platform:** Gene expression microarray

### Dataset Description

This dataset compares gene expression profiles between:
- **Case Group:** Granulosa cells from PCOD/PCOS patients
- **Control Group:** Granulosa cells from healthy individuals

**Granulosa cells** are important because they are steroid-producing cells in the ovary that are directly affected in PCOD. Their dysregulation contributes to the hormonal imbalances seen in the condition.

### Available Data Formats

#### CEL Files (Raw Data) ❌
- Individual probe intensity files
- Represent raw, unnormalized microarray data
- Require background correction, normalization, and summarization
- **Issue encountered:** Downloaded CEL files showed 0KB size, indicating extraction problems
- **Not recommended** for this quick analysis

#### Processed Expression Matrix (Normalized) ✅
- Already normalized and background-corrected
- Contains log2-transformed expression values
- One gene per row, one sample per column
- Ready for immediate differential expression analysis
- **Recommended:** Use this format

### How to Access the Data

1. **Web Interface:** Visit https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114419
2. **Download directly:** Look for "Series Matrix File(s)" section
3. **Via GEOquery in R:**
   ```r
   library(GEOquery)
   gse <- getGEO("GSE114419", GSEMatrix = TRUE)
   expr <- exprs(gse[[1]])
   ```

### Dataset Specifications

| Feature | Value |
|---------|-------|
| Number of Samples | ~50 (varies by matrix version) |
| Number of Genes | ~20,000+ |
| Expression Values | log2-transformed, normalized |
| Sample Groups | PCOD vs. Control |
| Cell Type | Granulosa cells |
| Technology | Microarray |

### Key Metadata

The phenotype data includes:
- **Characteristics:** Disease status (control/PCOD)
- **Title:** Sample description
- **Source Name:** Tissue type and condition

## Understanding PCOD/PCOS

### What is PCOD?
Polycystic Ovary Syndrome (PCOS/PCOD) is an endocrine disorder characterized by:
- Irregular menstrual cycles
- Elevated androgen levels (male hormones)
- Ovarian cysts
- **Insulin resistance** (key mechanism)
- Fertility issues

### Why Gene Expression Matters

In PCOD, the dysregulation of multiple genes affects:
1. **Steroidogenesis** - Hormone production
2. **Insulin Signaling** - Glucose metabolism
3. **Inflammation** - Chronic low-grade inflammation
4. **Cell proliferation** - Abnormal granulosa cell growth
5. **Cell differentiation** - Impaired cell development

By comparing healthy vs. PCOD granulosa cells, we can identify which genes and pathways are dysregulated.

### Polygenic Nature of PCOD

❌ **NOT monogenic:** One gene mutation doesn't cause PCOD
✅ **Polygenic:** Multiple genes contribute to the phenotype

This is why pathway analysis is crucial—we look for patterns across many genes rather than single-gene explanations.

## Related Research

### Gene Sets of Interest
- **Insulin Signaling Pathway:** INSR, IRS1, AKT1, GLUT4
- **Steroid Metabolism:** CYP11A1, CYP17A1, 3β-HSD
- **Inflammation:** TNF, IL-6, IL-8
- **Androgen Pathway:** AR, CYP19A1
- **Cell Cycle Regulation:** TP53, RB1

### Complementary Databases
- **NCBI Gene:** https://www.ncbi.nlm.nih.gov/gene/
- **UniProt:** https://www.uniprot.org/
- **PubMed:** https://pubmed.ncbi.nlm.nih.gov/ (for literature context)

## Pathway Databases Used

### 1. g:Profiler
- **URL:** https://biit.cs.ut.ee/gprofiler/
- **Strengths:** Comprehensive gene ontology and pathway annotation
- **Output:** GO terms, KEGG pathways, disease associations
- **Best for:** Quick functional enrichment

### 2. KEGG (Kyoto Encyclopedia of Genes and Genomes)
- **URL:** https://www.kegg.jp/
- **Strengths:** Detailed molecular interaction networks
- **Coverage:** ~500 pathways with disease associations
- **Best for:** Understanding pathway topology and cross-talk

### 3. Reactome
- **URL:** https://reactome.org/
- **Strengths:** Manually curated molecular reactions
- **Coverage:** Most comprehensive for signaling and metabolism
- **Best for:** Detailed mechanistic understanding

## Data Quality Considerations

### Normalization
The expression matrix is already normalized using standard microarray normalization procedures (likely RMA - Robust Multichip Average or similar).

### Batch Effects
- Already accounted for in the processed matrix
- Important to note when comparing multiple datasets

### Technical Variability
- Minimized through preprocessing
- Biological replicates (different individuals) represented

---

**Last Updated:** 2026-06-25  
**Status:** Initial data assessment complete
