library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(openxlsx)
library(ggplot2)

analysis_dir <- path.expand("~/bulk_RNA_seq/deseq2_analysis")
table_dir <- file.path(analysis_dir, "tables")
out_dir <- file.path(analysis_dir, "enrichment")

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

files <- list.files(
  table_dir,
  pattern="^(BNTvsCtrl|ONK_).*_DEG_padj005_annotated.csv$",
  full.names=TRUE
)

for (f in files) {
  comp <- sub("_DEG_padj005_annotated.csv", "", basename(f))
  message("Running enrichment for: ", comp)

  deg <- read.csv(f, check.names=FALSE)

  genes <- unique(na.omit(deg$gene_id))
  genes <- sub("[.].*$", "", genes)

  entrez <- bitr(
    genes,
    fromType="ENSEMBL",
    toType="ENTREZID",
    OrgDb=org.Hs.eg.db
  )

  if (nrow(entrez) < 5) {
    message("Too few mapped genes for ", comp)
    next
  }

  for (ont in c("BP", "CC", "MF")) {
    ego <- enrichGO(
      gene=entrez$ENTREZID,
      OrgDb=org.Hs.eg.db,
      keyType="ENTREZID",
      ont=ont,
      pAdjustMethod="BH",
      qvalueCutoff=0.05,
      readable=TRUE
    )

    write.xlsx(
      as.data.frame(ego),
      file.path(out_dir, paste0(comp, "_GO_", ont, ".xlsx")),
      overwrite=TRUE
    )

    if (nrow(as.data.frame(ego)) > 0) {
      pdf(file.path(out_dir, paste0(comp, "_GO_", ont, "_dotplot.pdf")),
          width=9, height=7)
      print(dotplot(ego, showCategory=20) + ggtitle(paste(comp, "GO", ont)))
      dev.off()
    }
  }

  kegg <- enrichKEGG(
    gene=entrez$ENTREZID,
    organism="hsa",
    pAdjustMethod="BH",
    qvalueCutoff=0.05
  )

  write.xlsx(
    as.data.frame(kegg),
    file.path(out_dir, paste0(comp, "_KEGG.xlsx")),
    overwrite=TRUE
  )

  if (nrow(as.data.frame(kegg)) > 0) {
    pdf(file.path(out_dir, paste0(comp, "_KEGG_dotplot.pdf")),
        width=9, height=7)
    print(dotplot(kegg, showCategory=20) + ggtitle(paste(comp, "KEGG")))
    dev.off()
  }
}

print("GO and KEGG enrichment completed successfully.")
