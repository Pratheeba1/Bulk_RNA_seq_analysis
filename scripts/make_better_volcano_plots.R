library(ggplot2)
library(ggrepel)

analysis_dir <- path.expand("~/bulk_RNA_seq/deseq2_analysis")
table_dir <- file.path(analysis_dir, "tables")
plot_dir <- file.path(analysis_dir, "plots", "better_volcano")

dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)

files <- list.files(
  table_dir,
  pattern="^(BNTvsCtrl|ONK_).*_all_results_annotated.csv$",
  full.names=TRUE
)

for (f in files) {
  comp <- sub("_all_results_annotated.csv", "", basename(f))
  message("Making volcano for: ", comp)

  df <- read.csv(f, check.names=FALSE)

  df$gene_label <- df$gene_name
  df$gene_label[
    is.na(df$gene_label) |
    df$gene_label == "" |
    grepl("^ENSG", df$gene_label)
  ] <- NA

  df$neglog10padj <- -log10(df$padj)
  df$neglog10padj[is.infinite(df$neglog10padj)] <- NA

  df$category <- "Not significant"
  df$category[
    !is.na(df$padj) &
    df$padj < 0.05 &
    df$log2FoldChange > 0
  ] <- "Up"

  df$category[
    !is.na(df$padj) &
    df$padj < 0.05 &
    df$log2FoldChange < 0
  ] <- "Down"

  up_n <- sum(df$category == "Up", na.rm=TRUE)
  down_n <- sum(df$category == "Down", na.rm=TRUE)

  label_df <- df[
    !is.na(df$padj) &
    df$padj < 0.05 &
    !is.na(df$gene_label),
  ]

  label_df <- label_df[order(label_df$padj), ]
  label_df <- head(label_df, 20)

  p <- ggplot(df, aes(x=log2FoldChange, y=neglog10padj)) +
    geom_point(aes(color=category), alpha=0.65, size=1.5) +
    scale_color_manual(
      values=c(
        "Down"="#2B6CB0",
        "Not significant"="grey75",
        "Up"="#D62728"
      )
    ) +
    geom_vline(xintercept=c(-1, 1), linetype="dashed", linewidth=0.4) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth=0.4) +
    geom_text_repel(
      data=label_df,
      aes(label=gene_label),
      size=3.2,
      max.overlaps=Inf,
      box.padding=0.35,
      point.padding=0.25,
      segment.size=0.25
    ) +
    scale_x_continuous(
      limits=c(-5, 5),
      breaks=c(-5, -3, -1, 0, 1, 3, 5)
    ) +
    labs(
      title=comp,
      subtitle=paste0("Significant genes: Up = ", up_n, " | Down = ", down_n),
      x=expression(Log[2]~fold~change),
      y=expression(-Log[10]~adjusted~italic(P)),
      color=""
    ) +
    theme_bw(base_size=13) +
    theme(
      plot.title=element_text(face="bold", size=18),
      plot.subtitle=element_text(size=13),
      legend.position="top",
      panel.grid.minor=element_blank(),
      panel.grid.major=element_line(linewidth=0.25),
      axis.title=element_text(size=14),
      axis.text=element_text(size=12)
    )

  ggsave(
    filename=file.path(plot_dir, paste0(comp, "_better_volcano.pdf")),
    plot=p,
    width=9,
    height=7
  )

  ggsave(
    filename=file.path(plot_dir, paste0(comp, "_better_volcano.png")),
    plot=p,
    width=9,
    height=7,
    dpi=300
  )
}

print("Better volcano plots completed.")
