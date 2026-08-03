library(tximport)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

base_dir <- path.expand("~/bulk_RNA_seq")
salmon_dir <- file.path(base_dir, "salmon_quant")
analysis_dir <- file.path(base_dir, "deseq2_analysis")

dir.create(file.path(analysis_dir, "tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(analysis_dir, "plots"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(analysis_dir, "results"), recursive=TRUE, showWarnings=FALSE)

samples <- read.csv(file.path(analysis_dir, "sample_sheet.csv"), stringsAsFactors=FALSE)
samples$condition <- factor(samples$condition)
samples$condition <- relevel(samples$condition, ref="Untreated")

files <- file.path(salmon_dir, samples$sample, "quant.sf")
names(files) <- samples$sample

if (!all(file.exists(files))) {
  print(files[!file.exists(files)])
  stop("Some quant.sf files are missing.")
}

tx2gene <- read.csv(file.path(analysis_dir, "tx2gene.csv"), stringsAsFactors=FALSE)

gene_anno <- read.csv(file.path(analysis_dir, "gene_annotation.csv"), stringsAsFactors=FALSE)
gene_anno$gene_id_clean <- sub("[.].*$", "", gene_anno$gene_id)

txi <- tximport(
  files,
  type="salmon",
  tx2gene=tx2gene,
  ignoreAfterBar=FALSE,
  ignoreTxVersion=FALSE,
  dropInfReps=TRUE
)

dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData=samples,
  design=~condition
)

keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]

dds <- DESeq(dds)

saveRDS(dds, file.path(analysis_dir, "results", "dds_object.rds"))

raw_counts <- counts(dds, normalized=FALSE)
norm_counts <- counts(dds, normalized=TRUE)
write.csv(norm_counts, file.path(analysis_dir, "tables", "normalized_counts.csv"))

vsd <- vst(dds, blind=FALSE)

pdf(file.path(analysis_dir, "plots", "PCA_plot.pdf"))
print(plotPCA(vsd, intgroup="condition"))
dev.off()

pdf(file.path(analysis_dir, "plots", "Sample_distance_heatmap.pdf"))
pheatmap(
  as.matrix(dist(t(assay(vsd)))),
  main="Sample distance heatmap"
)
dev.off()

comparisons <- list(
  BNTvsCtrl = c("condition", "BioNTech_formulation", "Untreated"),
  ONK_1to1vsCtrl = c("condition", "Oligo_1to1", "Untreated"),
  ONK_3to1vsCtrl = c("condition", "Oligo_3to1", "Untreated"),
  ONK_10to1vsCtrl = c("condition", "Oligo_10to1", "Untreated")
)

make_table <- function(res, comp_name) {
  group1 <- comparisons[[comp_name]][2]
  group2 <- comparisons[[comp_name]][3]

  group1_samples <- samples$sample[samples$condition == group1]
  group2_samples <- samples$sample[samples$condition == group2]

  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_id_clean <- sub("[.].*$", "", res_df$gene_id)

  count_df <- as.data.frame(raw_counts)
  count_df$gene_id <- rownames(count_df)
  count_df$gene_id_clean <- sub("[.].*$", "", count_df$gene_id)

  keep_count_cols <- c("gene_id_clean", group1_samples, group2_samples)
  count_df <- count_df[, keep_count_cols, drop=FALSE]
  colnames(count_df) <- c("gene_id_clean", paste0(colnames(count_df)[-1], "_count"))

  final <- merge(res_df, gene_anno, by="gene_id_clean", all.x=TRUE)
  final <- merge(final, count_df, by="gene_id_clean", all.x=TRUE)

  final$gene_id <- final$gene_id.x

  group1_count_cols <- paste0(group1_samples, "_count")
  group2_count_cols <- paste0(group2_samples, "_count")

  final[[group1]] <- rowMeans(final[, group1_count_cols, drop=FALSE], na.rm=TRUE)
  final[[group2]] <- rowMeans(final[, group2_count_cols, drop=FALSE], na.rm=TRUE)

  final$regulation <- "NOT_SIGNIFICANT"
  final$regulation[final$padj < 0.05 & final$log2FoldChange > 0] <- "UP"
  final$regulation[final$padj < 0.05 & final$log2FoldChange < 0] <- "DOWN"

  final <- final[, c(
  "gene_id",
  group1_count_cols,
  group2_count_cols,
  group1,
  group2,
  "log2FoldChange",
  "pvalue",
  "padj",
  "regulation",
  "gene_name",
  "gene_chr",
  "gene_start",
  "gene_end",
  "gene_strand",
  "gene_length",
  "gene_type",
  "gene_description"
)]

  final <- final[order(final$padj), ]
  return(final)
}

all_deg_genes <- c()

for (comp_name in names(comparisons)) {
  res <- results(dds, contrast=comparisons[[comp_name]])
  final <- make_table(res, comp_name)

  write.csv(
    final,
    file.path(analysis_dir, "tables", paste0(comp_name, "_all_results_annotated.csv")),
    row.names=FALSE
  )

  deg <- final[!is.na(final$padj) & final$padj < 0.05, ]
  up <- deg[deg$regulation == "UP", ]
  down <- deg[deg$regulation == "DOWN", ]

  write.csv(deg, file.path(analysis_dir, "tables", paste0(comp_name, "_DEG_padj005_annotated.csv")), row.names=FALSE)
  write.csv(up, file.path(analysis_dir, "tables", paste0(comp_name, "_UP_padj005_annotated.csv")), row.names=FALSE)
  write.csv(down, file.path(analysis_dir, "tables", paste0(comp_name, "_DOWN_padj005_annotated.csv")), row.names=FALSE)

  all_deg_genes <- unique(c(all_deg_genes, deg$gene_id))

  pdf(file.path(analysis_dir, "plots", paste0(comp_name, "_volcano.pdf")), width=8, height=7)
  print(
    EnhancedVolcano(
      final,
      lab=final$gene_name,
      x="log2FoldChange",
      y="padj",
      title=comp_name,
      pCutoff=0.05,
      FCcutoff=1,
      selectLab=head(final$gene_name[!is.na(final$padj) & final$padj < 0.05], 20)
    )
  )
  dev.off()

  topgenes <- head(final$gene_id[order(final$padj)], 50)
  topgenes <- topgenes[topgenes %in% rownames(vsd)]

  if (length(topgenes) > 2) {
    heatmat <- assay(vsd)[topgenes, ]
    gene_labels <- final$gene_name[match(topgenes, final$gene_id)]
    gene_labels[is.na(gene_labels) | gene_labels==""] <- topgenes[is.na(gene_labels) | gene_labels==""]
    rownames(heatmat) <- gene_labels

    pdf(file.path(analysis_dir, "plots", paste0(comp_name, "_heatmap_top50_gene_names.pdf")), width=10, height=12)
    pheatmap(
      heatmat,
      cluster_rows=TRUE,
      cluster_cols=TRUE,
      show_rownames=TRUE,
      annotation_col=data.frame(condition=samples$condition, row.names=samples$sample),
      main=paste0(comp_name, " top 50 DEGs")
    )
    dev.off()
  }
}

all_deg_genes <- unique(all_deg_genes)
all_deg_genes <- all_deg_genes[all_deg_genes %in% rownames(vsd)]

if (length(all_deg_genes) > 2) {
  combined_mat <- assay(vsd)[all_deg_genes, ]
  combined_anno <- gene_anno[
    match(sub("[.].*$", "", all_deg_genes), gene_anno$gene_id_clean),
  ]

  gene_labels <- combined_anno$gene_name
  gene_labels[is.na(gene_labels) | gene_labels==""] <- all_deg_genes[is.na(gene_labels) | gene_labels==""]
  rownames(combined_mat) <- gene_labels

  pdf(file.path(analysis_dir, "plots", "Combined_all_DEGs_heatmap.pdf"), width=11, height=16)
  pheatmap(
    combined_mat,
    cluster_rows=TRUE,
    cluster_cols=TRUE,
    show_rownames=FALSE,
    annotation_col=data.frame(condition=samples$condition, row.names=samples$sample),
    main="Combined all significant DEGs"
  )
  dev.off()
}

summary_table <- data.frame(
  comparison=character(),
  total_DEG=integer(),
  UP=integer(),
  DOWN=integer()
)

for (comp_name in names(comparisons)) {
  deg <- read.csv(
    file.path(analysis_dir, "tables", paste0(comp_name, "_DEG_padj005_annotated.csv")),
    check.names=FALSE
  )

  summary_table <- rbind(
    summary_table,
    data.frame(
      comparison=comp_name,
      total_DEG=nrow(deg),
      UP=sum(deg$log2FoldChange > 0, na.rm=TRUE),
      DOWN=sum(deg$log2FoldChange < 0, na.rm=TRUE)
    )
  )
}

write.csv(
  summary_table,
  file.path(analysis_dir, "tables", "DEG_summary_statistics.csv"),
  row.names=FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(analysis_dir, "results", "sessionInfo.txt")
)

print("DESeq2 analysis completed successfully.")
