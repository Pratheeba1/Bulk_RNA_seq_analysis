library(openxlsx)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)

analysis_dir <- path.expand("~/bulk_RNA_seq/deseq2_analysis")
table_dir <- file.path(analysis_dir, "tables")
out_dir <- file.path(analysis_dir, "excel_broad_functional_categories")

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

categorize_go <- function(go_text) {
  if (is.na(go_text) || go_text == "") return("Other/uncategorized")

  x <- tolower(go_text)
  cats <- c()

  if (grepl("immune|immunity|cytokine|chemokine|interferon|inflammatory|antigen|mhc|t cell|b cell|leukocyte|lymphocyte|macrophage|nf-kappa|toll-like", x))
    cats <- c(cats, "Immune/inflammatory")

  if (grepl("signal|signaling|receptor|kinase|phosphorylation|jak|stat|mapk|nf-kappa|growth factor", x))
    cats <- c(cats, "Signalling")

  if (grepl("cholesterol|sterol|lipid|fatty acid|isoprenoid|terpenoid", x))
    cats <- c(cats, "Lipid/cholesterol metabolism")

  if (grepl("amino acid|glutamine|glutamate|serine|cysteine|methionine|transporter|transport", x))
    cats <- c(cats, "Amino acid metabolism/transport")

  if (grepl("endoplasmic reticulum|unfolded protein|protein folding|heat shock|chaperone|stress|oxidative stress", x))
    cats <- c(cats, "Stress/ER stress/protein folding")

  if (grepl("ribosome|translation|trna|rrna|protein synthesis|peptide biosynthetic", x))
    cats <- c(cats, "Translation/ribosome")

  if (grepl("rna processing|mrna|splicing|rna binding|transcription|transcription factor", x))
    cats <- c(cats, "RNA/transcription regulation")

  if (grepl("cell cycle|mitosis|dna replication|dna repair|chromosome|chromatin", x))
    cats <- c(cats, "Cell cycle/DNA/chromatin")

  if (grepl("mitochondri|oxidative phosphorylation|respiratory chain|atp synthesis|electron transport", x))
    cats <- c(cats, "Mitochondria/energy metabolism")

  if (grepl("apoptosis|cell death|autophagy|lysosome|proteasome|ubiquitin", x))
    cats <- c(cats, "Cell death/protein degradation")

  if (grepl("extracellular matrix|collagen|adhesion|cytoskeleton|actin|microtubule", x))
    cats <- c(cats, "ECM/cytoskeleton/adhesion")

  if (length(cats) == 0) return("Other/uncategorized")
  paste(unique(cats), collapse="; ")
}

csv_files <- list.files(table_dir, pattern="[.]csv$", full.names=TRUE)

for (f in csv_files) {
  message("Processing: ", basename(f))

  df <- read.csv(f, check.names=FALSE)

  if (!"gene_name" %in% colnames(df)) {
    message("Skipping file without gene_name: ", basename(f))
    next
  }

  symbols <- unique(na.omit(df$gene_name))
  symbols <- symbols[!grepl("^ENSG", symbols)]

  anno <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys=symbols,
    columns=c("SYMBOL", "GENENAME", "ENTREZID", "GO", "ONTOLOGY"),
    keytype="SYMBOL"
  )

  gene_desc <- anno[!duplicated(anno$SYMBOL), c("SYMBOL", "GENENAME", "ENTREZID")]

  go_terms <- AnnotationDbi::select(
    GO.db,
    keys=unique(na.omit(anno$GO)),
    columns=c("GOID", "TERM"),
    keytype="GOID"
  )

  anno$GO_term <- go_terms$TERM[match(anno$GO, go_terms$GOID)]

  go_collapsed <- aggregate(
    GO_term ~ SYMBOL,
    data=anno[!is.na(anno$GO_term), ],
    FUN=function(x) paste(unique(x), collapse="; ")
  )

  df$gene_description <- gene_desc$GENENAME[match(df$gene_name, gene_desc$SYMBOL)]
  df$entrez_id <- gene_desc$ENTREZID[match(df$gene_name, gene_desc$SYMBOL)]
  df$GO_terms <- go_collapsed$GO_term[match(df$gene_name, go_collapsed$SYMBOL)]
  df$pathway_or_category <- sapply(df$GO_terms, categorize_go)
  df$immune_gene <- ifelse(grepl("Immune", df$pathway_or_category), "YES", "NO")

  first_cols <- c(
    "gene_id",
    "gene_name",
    "gene_description",
    "entrez_id",
    "immune_gene",
    "pathway_or_category",
    "GO_terms"
  )

  existing_first <- first_cols[first_cols %in% colnames(df)]
  remaining <- setdiff(colnames(df), existing_first)
  df <- df[, c(existing_first, remaining), drop=FALSE]

  out_xlsx <- file.path(
    out_dir,
    paste0(tools::file_path_sans_ext(basename(f)), "_GO_functional_categories.xlsx")
  )

  write.xlsx(df, out_xlsx, overwrite=TRUE)
}

print("Broad GO-based functional categories added to all eligible Excel files.")
