# Bulk RNA-seq Analysis

## Overview

This repository contains the complete bulk RNA-seq analysis workflow for HEK293T cells transfected with lipid nanoparticle (LNP)-encapsulated mRNA formulations.

The workflow includes:

- Quality control (FastQC)
- Adapter trimming (FastP)
- Transcript quantification (Salmon)
- Gene-level summarization (tximport)
- Differential expression analysis (DESeq2)
- Functional enrichment (GO and KEGG)
- Visualization (PCA, Volcano plots, Heatmaps)

---

## Workflow

FASTQ
↓
FastQC
↓
FastP
↓
Salmon
↓
tximport
↓
DESeq2
↓
Differentially Expressed Genes
↓
GO / KEGG Enrichment
↓
Visualization

---

## Repository Structure

```
scripts/
metadata/
results/
figures/
docs/
environment/
```

---

## Author

Pratheeba Pandiaraj

PhD Student

University of Würzburg
