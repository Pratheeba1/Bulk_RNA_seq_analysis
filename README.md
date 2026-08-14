Bulk RNA-seq Analysis

Overview

This repository contains the scripts, results, and visualizations generated for bulk RNA-seq analysis of HEK293T cells transfected with mRNA-loaded lipid nanoparticle (LNP) formulations.
The analysis was performed to investigate transcriptomic changes induced by different LNP/mRNA formulations compared with untreated control cells.


Experimental Conditions

The following experimental groups were analyzed:

•	Untreated control
•	BNT formulation
•	TriB35/mRNA 1:1
•	TriB35/mRNA 3:1
•	TriB35/mRNA 10:1

Biological replicates were included for each experimental condition.

RNA-seq Analysis Workflow

The overall workflow was:

FASTQ → FastQC → fastp → Salmon → tximport → DESeq2 → Differential Expression Analysis → Functional Enrichment → Visualization

1. Quality Control
Raw sequencing reads were assessed using FastQC.
2. Read Processing
Adapter trimming and quality filtering were performed using fastp.
3. Transcript Quantification
Transcript abundance was quantified using Salmon.
4. Gene-level Quantification
Salmon transcript-level abundance estimates were imported into R using tximport and summarized at the gene level.
5. Differential Expression Analysis
Differential gene expression analysis was performed using DESeq2.
Individual LNP/mRNA formulations were compared with the untreated control group.
6. Functional Enrichment Analysis
Differentially expressed genes were further investigated using functional enrichment analyses, including:
•	Gene Ontology (GO)
•	KEGG pathway enrichment
7. Visualization
The analysis includes several visualization approaches:
•	Principal Component Analysis (PCA)
•	Volcano plots
•	Differential gene expression heatmaps
•	Pathway and functional enrichment visualizations

Repository Structure
Bulk_RNA_seq_analysis/
│
├── scripts/
│   └── Analysis scripts used for RNA-seq processing,
│       differential expression and visualization
│
├── figures/
│   └── PCA, volcano plots, heatmaps and enrichment figures
│
├── results/
│   └── Differential expression and enrichment results
│
├── metadata/
│   └── Sample information and analysis metadata
│
├── docs/
│   └── Analysis documentation and methodological information
│
├── environment/
│   └── Software and R package information required
│       for reproducibility
│
└── README.md
Main Tools

The workflow uses:
•	FastQC
•	fastp
•	Salmon
•	R
•	tximport
•	DESeq2
•	ggplot2
•	EnhancedVolcano
•	pheatmap
•	clusterProfiler

Figures
Generated figures are available in the figures/ directory.

These include visualizations for comparisons between the individual LNP/mRNA formulations and untreated control cells.

Results
Differential expression and pathway enrichment results are available in the results/ directory.

Reproducibility
Information about the computational environment, software versions, and R packages used for the analysis is provided in the environment/ directory.

Author
Pratheeba Pandiaraj
PhD Student
University of Würzburg

PhD Student
University of Würzburg
