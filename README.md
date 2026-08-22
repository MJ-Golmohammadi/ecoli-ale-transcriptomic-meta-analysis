# Cross-Study Transcriptomic Meta-Analysis Reveals Conserved Adaptive Programs in *Escherichia coli* K-12

[![R](https://img.shields.io/badge/R-4.3.3-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)

## Overview

This repository contains the complete computational workflow, analysis scripts, results, and figures generated for a cross-study transcriptomic meta-analysis of adaptive laboratory evolution (ALE) experiments in *Escherichia coli* K-12.

The study investigates whether transcriptional responses observed during independent evolutionary experiments are reproducible across diverse genetic and environmental evolutionary contexts.

Publicly available transcriptomic datasets were integrated using dataset-specific differential-expression analyses, directional Signed Z-score harmonization, weighted Stouffer integration of related *menF*-associated comparisons, and Robust Rank Aggregation (RRA).

The resulting conserved transcriptional signature was further characterized using functional enrichment and protein–protein interaction (PPI) network analyses.

---

## Research Question

**Which transcriptional responses are reproducibly associated with adaptation across independent *E. coli* K-12 ALE experiments?**

The study specifically aims to distinguish recurrent transcriptional responses that extend across evolutionary contexts from responses that are specific to individual experimental conditions.

---

## Study Design

The cross-study meta-analysis consisted of **seven study-level inputs**:

1. A combined *menF*-associated signature derived from three related comparisons.
2. GSE140478
3. GSE158959
4. GSE17276
5. GSE206196
6. GSE33147
7. GSE316857

The three *menF*-associated comparisons were integrated before the principal meta-analysis and treated as a single independent study-level input.

### Transcriptomic datasets

| GEO accession | Evolutionary context | Analysis |
|---|---|---|
| GSE224889 | *menF*-associated evolution | DESeq2 |
| GSE122779 | Δ*menF*Δ*entC* evolution | Processed expression analysis |
| GSE140478 | Thermal adaptation | DESeq2 |
| GSE158959 | Isoprenol-associated adaptation | DESeq2 |
| GSE17276 | Long-term evolution | Microarray / limma |
| GSE206196 | Metabolic evolution | Microarray / limma |
| GSE33147 | Glycerol-associated evolution | Affymetrix / RMA + limma |
| GSE316857 | Homoserine-associated evolution | Processed expression / limma |

All original datasets were obtained from the **NCBI Gene Expression Omnibus (GEO)** and should be cited according to the original study publications.

---

## Analytical Workflow

```text
Publicly available GEO datasets
              │
              ▼
     Dataset-specific preprocessing
              │
              ▼
    Differential-expression analysis
              │
       ┌──────┴──────┐
       │             │
     DESeq2         limma
       │             │
       └──────┬──────┘
              ▼
       log2 fold change
       + adjusted P value
              │
              ▼
      Signed Z-score
      harmonization
              │
              ▼
 Integration of related
 menF-associated comparisons
              │
       Weighted Stouffer
              │
              ▼
     Seven study-level
       ranked inputs
              │
              ▼
 Robust Rank Aggregation
            (RRA)
              │
       ┌──────┴──────┐
       ▼             ▼
  Upregulated    Downregulated
    genes           genes
       └──────┬──────┘
              ▼
 Directional consistency
       + recurrence
        filtering
              │
              ▼
   Conserved core genes
              │
       ┌──────┴──────┐
       ▼             ▼
 Functional       PPI network
 enrichment       analysis
                       │
                ┌──────┴──────┐
                ▼             ▼
            Hub genes       MCODE
