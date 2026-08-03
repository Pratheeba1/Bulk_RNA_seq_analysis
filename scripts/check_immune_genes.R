library(openxlsx)

analysis_dir <- path.expand("~/bulk_RNA_seq/deseq2_analysis")
table_dir <- file.path(analysis_dir, "tables")
out_dir <- file.path(analysis_dir, "immune_gene_check")

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

immune_genes <- c(
  "IFNB1","IFNA1","IFNA2","IFNG",
  "IFIT1","IFIT2","IFIT3","IFI6","IFI27","IFI44","IFI44L","IFIH1",
  "ISG15","ISG20","OAS1","OAS2","OAS3","OASL","MX1","MX2",
  "RSAD2","BST2","DDX58","IRF1","IRF3","IRF7","IRF9",
  "STAT1","STAT2","STAT3","STAT6",
  "JAK1","JAK2","TYK2",
  "IL1A","IL1B","IL6","IL8","IL10","IL12A","IL12B","IL18",
  "TNF","TNFAIP3","TNFRSF1A","TNFRSF11A",
  "CXCL1","CXCL2","CXCL8","CXCL9","CXCL10","CXCL11","CXCL16",
  "CCL2","CCL3","CCL4","CCL5","CCL20",
  "TLR1","TLR2","TLR3","TLR4","TLR7","TLR8","TLR9",
  "NFKB1","NFKB2","NFKBIA","RELA","RELB",
  "SOCS1","SOCS2","SOCS3",
  "HLA-A","HLA-B","HLA-C","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1",
  "B2M","TAP1","TAP2","CD74",
  "C1QA","C1QB","C1QC","C3","C4A","C4B",
  "CD14","CD40","CD80","CD86","ICAM1","VCAM1"
)

files <- list.files(
  table_dir,
  pattern="^(BNTvsCtrl|ONK_).*_DEG_padj005_annotated.csv$",
  full.names=TRUE
)

immune_summary <- data.frame()

for (f in files) {
  comp <- sub("_DEG_padj005_annotated.csv", "", basename(f))
  deg <- read.csv(f, check.names=FALSE)

  hits <- deg[deg$gene_name %in% immune_genes, ]

  if (nrow(hits) > 0) {
    hits$comparison <- comp
    hits$immune_category <- "curated_immune_gene"

    hits <- hits[, c(
      "comparison",
      "gene_name",
      "gene_id",
      "log2FoldChange",
      "padj",
      "regulation",
      "gene_chr",
      "gene_start",
      "gene_end",
      "gene_strand",
      "gene_type",
      "immune_category"
    )]

    immune_summary <- rbind(immune_summary, hits)
  }
}

immune_summary <- immune_summary[order(immune_summary$comparison, immune_summary$padj), ]

write.csv(
  immune_summary,
  file.path(out_dir, "immune_DEG_summary.csv"),
  row.names=FALSE
)

write.xlsx(
  immune_summary,
  file.path(out_dir, "immune_DEG_summary.xlsx"),
  overwrite=TRUE
)

count_summary <- aggregate(
  gene_name ~ comparison + regulation,
  data=immune_summary,
  FUN=length
)

colnames(count_summary)[3] <- "immune_gene_count"

write.csv(
  count_summary,
  file.path(out_dir, "immune_DEG_count_summary.csv"),
  row.names=FALSE
)

write.xlsx(
  count_summary,
  file.path(out_dir, "immune_DEG_count_summary.xlsx"),
  overwrite=TRUE
)

print(immune_summary)
print(count_summary)
print("Immune gene check completed.")
