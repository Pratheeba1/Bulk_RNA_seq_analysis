# Computational Methods

## Transcript quantification

RNA-sequencing reads were quantified at the transcript level using
Salmon against a human transcriptome reference.

## Gene-level quantification

Transcript-level abundance estimates were imported into R using
tximport and summarized for downstream gene-level analysis.

## Differential expression

Differential expression analysis was performed using DESeq2.

DESeq2 performs library-size normalization using size factors, estimates
gene-wise dispersion, and fits negative binomial generalized linear
models to evaluate differential expression between experimental
conditions.

P-values were corrected for multiple testing using the
Benjamini-Hochberg procedure. Genes with an adjusted p-value
(padj) < 0.05 were considered statistically significant.

The following comparisons were evaluated:

- BNT vs untreated control
- TriB35/mRNA 1:1 vs untreated control
- TriB35/mRNA 3:1 vs untreated control
- TriB35/mRNA 10:1 vs untreated control

## Functional enrichment

Differentially expressed genes were subjected to Gene Ontology (GO)
and KEGG pathway enrichment analyses.

GO enrichment was evaluated for:

- Biological Process (BP)
- Cellular Component (CC)
- Molecular Function (MF)

## Visualization

Differential expression and pathway-level responses were visualized
using volcano plots, heatmaps, and enrichment plots.

Detailed R package and software version information is provided in
`environment/R_sessionInfo.txt`.
