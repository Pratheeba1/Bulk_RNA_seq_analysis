ann <- read.csv("~/bulk_RNA_seq/deseq2_analysis/gene_annotation.csv",
                stringsAsFactors = FALSE)

ann$gene_id_clean <- sub("[.].*$", "", ann$gene_id)

files <- list.files("~/bulk_RNA_seq/deseq2_analysis/tables",
                    pattern = "_DEG_padj005.csv$",
                    full.names = TRUE)

for (f in files) {
  deg <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)

  deg$gene_id_clean <- sub("[.].*$", "", deg[[1]])

  deg2 <- merge(deg, ann, by = "gene_id_clean", all.x = TRUE)

  out <- sub(".csv$", "_annotated.csv", f)

  write.csv(deg2, out, row.names = FALSE)

  cat("Written:", out, "\n")
}
