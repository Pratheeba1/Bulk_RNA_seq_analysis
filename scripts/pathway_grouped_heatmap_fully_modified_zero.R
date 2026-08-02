# ============================================================
# Pathway-grouped log2FoldChange heatmap
# Comirnaty and TriB35/mRNA formulations versus untreated control
# ============================================================

# ------------------------------------------------------------
# 1. Install and load required packages
# ------------------------------------------------------------

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "stringr",
  "pheatmap",
  "openxlsx"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(pheatmap)
  library(openxlsx)
})

# ------------------------------------------------------------
# 2. Input and output paths
# ------------------------------------------------------------

input_file <- "C:/Users/ENVY/Downloads/integrated_results.xlsx"

output_dir <- "C:/Users/ENVY/Downloads/Pathway_grouped_heatmap"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    paste0(
      "The Excel file was not found:\n",
      input_file,
      "\n\nCheck that the file is named integrated_results.xlsx ",
      "and is located in the Downloads folder."
    )
  )
}

message("Input file found: ", input_file)
message("Output folder: ", output_dir)

# ------------------------------------------------------------
# 3. Excel sheet names and labels
# ------------------------------------------------------------

sheet_labels <- c(
  ComirnatyvsCntrl = "Comirnaty vs Control",
  ONK1_1vsCntrl    = "TriB35/mRNA 1:1 vs Control",
  ONK3_1vsCntrl    = "TriB35/mRNA 3:1 vs Control",
  ONK10_1vsCntrl   = "TriB35/mRNA 10:1 vs Control"
)

available_sheets <- excel_sheets(input_file)

missing_sheets <- setdiff(
  names(sheet_labels),
  available_sheets
)

if (length(missing_sheets) > 0) {
  stop(
    paste0(
      "The following expected sheets are missing:\n",
      paste(missing_sheets, collapse = ", "),
      "\n\nAvailable sheets are:\n",
      paste(available_sheets, collapse = ", ")
    )
  )
}

# ------------------------------------------------------------
# 4. Read and combine all four sheets
# ------------------------------------------------------------

required_columns <- c(
  "gene_name",
  "gene_description",
  "immune_gene",
  "pathway_or_category",
  "log2FoldChange",
  "padj"
)

all_results_list <- list()

for (sheet_name in names(sheet_labels)) {

  message("Reading sheet: ", sheet_name)

  df <- read_excel(
    input_file,
    sheet = sheet_name
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(df)
  )

  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Missing columns in sheet ", sheet_name, ":\n",
        paste(missing_columns, collapse = ", ")
      )
    )
  }

  df <- df %>%
    mutate(
      comparison = sheet_labels[[sheet_name]],

      gene_name = toupper(
        trimws(
          as.character(gene_name)
        )
      ),

      gene_description = ifelse(
        is.na(gene_description),
        "",
        as.character(gene_description)
      ),

      immune_gene = ifelse(
        is.na(immune_gene),
        "",
        as.character(immune_gene)
      ),

      pathway_or_category = ifelse(
        is.na(pathway_or_category),
        "",
        as.character(pathway_or_category)
      ),

      log2FoldChange = suppressWarnings(
        as.numeric(log2FoldChange)
      ),

      padj = suppressWarnings(
        as.numeric(padj)
      ),

      regulation = if ("regulation" %in% colnames(df)) {
        as.character(regulation)
      } else {
        ifelse(log2FoldChange > 0, "UP", "DOWN")
      }
    ) %>%
    filter(
      !is.na(gene_name),
      gene_name != "",
      !is.na(log2FoldChange)
    )

  all_results_list[[sheet_name]] <- df
}

all_results <- bind_rows(all_results_list)

if (nrow(all_results) == 0) {
  stop("No usable gene rows were found in the Excel workbook.")
}

# ------------------------------------------------------------
# 5. Curated pathway gene lists
# ------------------------------------------------------------

lipid_genes <- c(
  "HMGCR", "HMGCS1", "SQLE", "MSMO1",
  "INSIG1", "INSIG2", "LDLR", "FDFT1",
  "FDPS", "IDI1", "MVK", "MVD", "LSS",
  "DHCR7", "DHCR24", "SC5D",
  "SREBF1", "SREBF2",
  "ABCG1", "ABCA1", "ACAT2",
  "MLYCD", "VLDLR", "SLCO2A1"
)

stress_genes <- c(
  "CHAC1", "TRIB3", "SESN2", "GDF15",
  "STC1", "HSPA5", "HSPA8", "HSPH1",
  "HSPB8", "DDIT3", "ATF3", "ATF4",
  "XBP1", "ERN1", "EIF2AK3",
  "DNAJB9", "HERPUD1", "LAMP3",
  "ASNS", "PPP1R15A"
)

amino_carbon_genes <- c(
  "CTH", "GPT2", "PSAT1", "PHGDH",
  "SHMT2", "ASNS", "SLC7A11",
  "SLC1A5", "SLC7A5", "SLC38A2",
  "GOT1", "GOT2", "GLS", "GLUL",
  "BCAT1", "BCAT2", "MTHFD2",
  "ALDOA", "ENO1", "ENO2",
  "GAPDH", "PKM", "LDHA", "PDK1",
  "SLC2A1", "NDUFA4L2"
)

immune_genes <- c(
  "SOCS1", "SOCS2", "SOCS3",
  "STAT1", "STAT2", "STAT3", "STAT6",
  "NFKB1", "NFKB2", "NFKBIA",
  "RELA", "RELB",
  "TAP1", "TAP2", "B2M",
  "TNFRSF1A", "TNFRSF11A",
  "CXCL16", "CXCL8", "CXCL10",
  "CCL2", "CCL5",
  "IFIT1", "IFIT2", "IFIT3",
  "ISG15", "OAS1", "OAS2", "OAS3",
  "MX1", "MX2",
  "IFIH1", "DDX58",
  "IRF1", "IRF3", "IRF7",
  "HLA-A", "HLA-B", "HLA-C",
  "CD74", "TRAF1"
)

hif_genes <- c(
  "HIF1A", "ARNT", "NDUFA4L2",
  "PDK1", "LDHA", "SLC2A1",
  "VEGFA", "BNIP3",
  "ADM", "ADM2",
  "CA9", "EGLN1", "EGLN3",
  "ENO1", "ENO2", "ALDOA",
  "PGK1", "HK2", "GDF15"
)

trafficking_genes <- c(
  "LAMP1", "LAMP2", "LAMP3",
  "RIN1", "RAB5A", "RAB7A", "RAB11A",
  "EEA1", "VPS4A", "VPS4B",
  "CHMP2A", "CHMP3", "CHMP4B",
  "TSG101", "VCP",
  "SQSTM1", "ATG5", "ATG7",
  "BECN1", "TFEB"
)

# ------------------------------------------------------------
# 6. Assign genes to pathway groups
# ------------------------------------------------------------

assign_pathway <- function(
    gene,
    description,
    category,
    immune_annotation
) {

  gene <- toupper(
    trimws(
      as.character(gene)
    )
  )

  category_lower <- tolower(
    as.character(category)
  )

  text_to_search <- tolower(
    paste(
      gene,
      description,
      category,
      immune_annotation,
      sep = " "
    )
  )

  # Use workbook annotations first
  if (
    str_detect(
      category_lower,
      "cholesterol|sterol|lipid"
    )
  ) {
    return("Lipid / cholesterol metabolism")
  }

  if (
    str_detect(
      category_lower,
      "stress response|er stress|protein folding"
    )
  ) {
    return("Cell stress response / protein folding")
  }

  if (
    str_detect(
      category_lower,
      "amino acid|carbon metabolism"
    )
  ) {
    return("Amino acid / carbon metabolism")
  }

  if (
    str_detect(
      category_lower,
      "immune"
    )
  ) {
    return("Immune signalling")
  }

  # Supplement using curated genes and keywords
  if (
    gene %in% lipid_genes ||
    str_detect(
      text_to_search,
      paste0(
        "cholesterol|sterol|steroid|lipid|",
        "fatty acid|terpenoid|isoprenoid"
      )
    )
  ) {
    return("Lipid / cholesterol metabolism")
  }

  if (
    gene %in% stress_genes ||
    str_detect(
      text_to_search,
      paste0(
        "endoplasmic reticulum|er stress|",
        "unfolded protein|protein folding|",
        "chaperone|integrated stress"
      )
    )
  ) {
    return("Cell stress response / protein folding")
  }

  if (
    gene %in% amino_carbon_genes ||
    str_detect(
      text_to_search,
      paste0(
        "amino acid|carbon metabolism|glycolysis|",
        "glucose|glutamate|serine|cysteine|",
        "transaminase|tricarboxylic"
      )
    )
  ) {
    return("Amino acid / carbon metabolism")
  }

  if (
    gene %in% immune_genes ||
    str_detect(
      text_to_search,
      paste0(
        "immune|cytokine|chemokine|interferon|",
        "jak.?stat|nf.?kappa|nfkb|antigen|inflamm"
      )
    )
  ) {
    return("Immune signalling")
  }

  if (
    gene %in% hif_genes ||
    str_detect(
      text_to_search,
      "hif|hypoxia|oxygen response"
    )
  ) {
    return("HIF-1 / hypoxia response")
  }

  if (
    gene %in% trafficking_genes ||
    str_detect(
      text_to_search,
      paste0(
        "endosome|lysosome|vesicle|",
        "autophagy|intracellular trafficking|endocyt"
      )
    )
  ) {
    return("Endosomal / lysosomal trafficking")
  }

  return(NA_character_)
}

all_results$pathway_group <- mapply(
  FUN = assign_pathway,
  gene = all_results$gene_name,
  description = all_results$gene_description,
  category = all_results$pathway_or_category,
  immune_annotation = all_results$immune_gene,
  USE.NAMES = FALSE
)

grouped_results <- all_results %>%
  filter(
    !is.na(pathway_group)
  )

if (nrow(grouped_results) == 0) {
  stop(
    "No genes could be assigned to the selected pathway groups."
  )
}

# ------------------------------------------------------------
# 7. Keep genes significant in at least one comparison
# ------------------------------------------------------------

grouped_results <- grouped_results %>%
  group_by(gene_name) %>%
  filter(
    any(
      !is.na(padj) &
      padj < 0.05
    )
  ) %>%
  ungroup()

if (nrow(grouped_results) == 0) {
  stop(
    "No pathway-assigned genes passed padj < 0.05."
  )
}

# ------------------------------------------------------------
# 8. Rank genes and select top genes per pathway
# ------------------------------------------------------------

pathway_order <- c(
  "Lipid / cholesterol metabolism",
  "Cell stress response / protein folding",
  "Amino acid / carbon metabolism",
  "Endosomal / lysosomal trafficking",
  "HIF-1 / hypoxia response",
  "Immune signalling"
)

gene_ranking <- grouped_results %>%
  group_by(
    gene_name,
    pathway_group
  ) %>%
  summarise(
    maximum_absolute_log2FC = max(
      abs(log2FoldChange),
      na.rm = TRUE
    ),

    minimum_padj = if (
      all(is.na(padj))
    ) {
      NA_real_
    } else {
      min(
        padj,
        na.rm = TRUE
      )
    },

    .groups = "drop"
  ) %>%
  mutate(
    pathway_group = factor(
      pathway_group,
      levels = pathway_order
    )
  ) %>%
  filter(
    !is.na(pathway_group)
  )

# Change this if you want more/fewer genes per pathway
genes_per_pathway <- 8

selected_genes <- gene_ranking %>%
  group_by(pathway_group) %>%
  arrange(
    desc(maximum_absolute_log2FC),
    minimum_padj,
    .by_group = TRUE
  ) %>%
  slice_head(
    n = genes_per_pathway
  ) %>%
  ungroup()

if (nrow(selected_genes) == 0) {
  stop(
    "No genes remained after pathway ranking."
  )
}

message("")
message("Selected genes per pathway:")
print(
  table(
    selected_genes$pathway_group,
    useNA = "ifany"
  )
)

# ------------------------------------------------------------
# 9. Build heatmap matrix
# ------------------------------------------------------------

comparison_order <- unname(
  sheet_labels
)

plot_data <- grouped_results %>%
  filter(
    gene_name %in% selected_genes$gene_name
  ) %>%
  select(
    gene_name,
    pathway_group,
    comparison,
    log2FoldChange,
    padj
  ) %>%
  arrange(
    gene_name,
    comparison,
    padj
  ) %>%
  distinct(
    gene_name,
    comparison,
    .keep_all = TRUE
  )

heatmap_data <- plot_data %>%
  select(
    gene_name,
    comparison,
    log2FoldChange
  ) %>%
  pivot_wider(
    names_from = comparison,
    values_from = log2FoldChange
  )

heatmap_matrix <- as.data.frame(
  heatmap_data
)

rownames(heatmap_matrix) <- heatmap_matrix$gene_name
heatmap_matrix$gene_name <- NULL

for (comparison_name in comparison_order) {

  if (
    !comparison_name %in%
    colnames(heatmap_matrix)
  ) {
    heatmap_matrix[[comparison_name]] <- NA_real_
  }
}

heatmap_matrix <- heatmap_matrix[
  ,
  comparison_order,
  drop = FALSE
]

heatmap_matrix <- as.matrix(
  heatmap_matrix
)

storage.mode(heatmap_matrix) <- "numeric"

if (all(is.na(heatmap_matrix))) {
  stop(
    "The heatmap matrix contains only missing values."
  )
}

# ------------------------------------------------------------
# 10. Order genes by pathway and effect size
# ------------------------------------------------------------

row_annotation_order <- selected_genes %>%
  select(
    gene_name,
    pathway_group,
    maximum_absolute_log2FC
  ) %>%
  distinct() %>%
  arrange(
    pathway_group,
    desc(maximum_absolute_log2FC),
    gene_name
  )

row_order <- row_annotation_order$gene_name[
  row_annotation_order$gene_name %in%
    rownames(heatmap_matrix)
]

heatmap_matrix <- heatmap_matrix[
  row_order,
  ,
  drop = FALSE
]

annotation_row <- data.frame(
  `Biological pathway` = factor(
    row_annotation_order$pathway_group[
      match(
        rownames(heatmap_matrix),
        row_annotation_order$gene_name
      )
    ],
    levels = pathway_order
  ),
  check.names = FALSE
)

rownames(annotation_row) <- rownames(
  heatmap_matrix
)

valid_rows <- !is.na(
  annotation_row$`Biological pathway`
)

heatmap_matrix <- heatmap_matrix[
  valid_rows,
  ,
  drop = FALSE
]

annotation_row <- annotation_row[
  valid_rows,
  ,
  drop = FALSE
]

if (nrow(heatmap_matrix) == 0) {
  stop(
    "No valid genes remained for plotting."
  )
}

# ------------------------------------------------------------
# 11. Add gaps between pathway groups
# ------------------------------------------------------------

pathway_vector <- as.character(
  annotation_row$`Biological pathway`
)

pathway_runs <- rle(
  pathway_vector
)

if (
  length(pathway_runs$lengths) > 1
) {

  gaps_row <- cumsum(
    pathway_runs$lengths
  )

  gaps_row <- gaps_row[
    -length(gaps_row)
  ]

} else {

  gaps_row <- NULL
}

# ------------------------------------------------------------
# 12. Define colours
# ------------------------------------------------------------

annotation_colours <- list(
  `Biological pathway` = c(
    "Lipid / cholesterol metabolism" = "#E69F00",
    "Cell stress response / protein folding" = "#CC79A7",
    "Amino acid / carbon metabolism" = "#009E73",
    "Endosomal / lysosomal trafficking" = "#56B4E9",
    "HIF-1 / hypoxia response" = "#0072B2",
    "Immune signalling" = "#D55E00"
  )
)

maximum_scale <- max(
  abs(heatmap_matrix),
  na.rm = TRUE
)

if (
  !is.finite(maximum_scale) ||
  maximum_scale == 0
) {
  maximum_scale <- 1
}

# Limit scale so that one extreme gene does not dominate
colour_limit <- min(
  maximum_scale,
  3
)

heatmap_breaks <- seq(
  -colour_limit,
  colour_limit,
  length.out = 101
)

heatmap_colours <- colorRampPalette(
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)(100)

# Replace missing values by 0.00 for visualization.
# This means genes absent from a comparison are displayed as zero.
heatmap_matrix[is.na(heatmap_matrix)] <- 0

number_matrix <- matrix(
  sprintf("%.2f", heatmap_matrix),
  nrow = nrow(heatmap_matrix),
  ncol = ncol(heatmap_matrix),
  dimnames = dimnames(heatmap_matrix)
)

heatmap_title <-
  "Pathway-grouped transcriptional responses to mRNA-LNP formulations"

# ------------------------------------------------------------
# 13. Heatmap drawing function
# ------------------------------------------------------------

draw_heatmap <- function() {

  pheatmap(
    heatmap_matrix,

    color = heatmap_colours,
    breaks = heatmap_breaks,

    cluster_rows = FALSE,
    cluster_cols = FALSE,

    gaps_row = gaps_row,

    annotation_row = annotation_row,
    annotation_colors = annotation_colours,
    drop_levels = TRUE,

    display_numbers = number_matrix,
    number_color = "black",

    fontsize_number = 6,
    fontsize_row = 8,
    fontsize_col = 9,

    angle_col = 45,

    border_color = "grey85",

    treeheight_row = 0,
    treeheight_col = 0,

    main = heatmap_title
  )
}

# ------------------------------------------------------------
# 14. Save PDF
# ------------------------------------------------------------

pdf_height <- max(
  10,
  nrow(heatmap_matrix) * 0.31
)

pdf(
  file.path(
    output_dir,
    "Pathway_grouped_log2FC_heatmap.pdf"
  ),
  width = 12,
  height = pdf_height,
  useDingbats = FALSE
)

draw_heatmap()

dev.off()

# ------------------------------------------------------------
# 15. Save PNG
# ------------------------------------------------------------

png_height <- max(
  3000,
  nrow(heatmap_matrix) * 120
)

png(
  file.path(
    output_dir,
    "Pathway_grouped_log2FC_heatmap.png"
  ),
  width = 4000,
  height = png_height,
  res = 300
)

draw_heatmap()

dev.off()

# ------------------------------------------------------------
# 16. Export Excel values and diagnostic tables
# ------------------------------------------------------------

heatmap_values <- as.data.frame(
  heatmap_matrix,
  check.names = FALSE
)

heatmap_values <- cbind(
  gene_name = rownames(heatmap_matrix),

  pathway_group = as.character(
    annotation_row$`Biological pathway`
  ),

  heatmap_values
)

selected_gene_details <- all_results %>%
  filter(
    gene_name %in%
      rownames(heatmap_matrix)
  ) %>%
  select(
    comparison,
    gene_name,
    gene_description,
    pathway_group,
    pathway_or_category,
    immune_gene,
    log2FoldChange,
    padj,
    regulation
  ) %>%
  arrange(
    factor(
      pathway_group,
      levels = pathway_order
    ),
    gene_name,
    factor(
      comparison,
      levels = comparison_order
    )
  )

assignment_summary <- all_results %>%
  count(
    pathway_or_category,
    pathway_group,
    name = "number_of_rows",
    sort = TRUE
  )

genes_per_pathway_summary <- grouped_results %>%
  distinct(
    gene_name,
    pathway_group
  ) %>%
  count(
    pathway_group,
    name = "number_of_unique_genes"
  )

write.xlsx(
  list(
    Heatmap_values = heatmap_values,
    Selected_gene_details = selected_gene_details,
    Selected_gene_ranking = as.data.frame(selected_genes),
    Assignment_summary = assignment_summary,
    Genes_per_pathway = genes_per_pathway_summary
  ),

  file.path(
    output_dir,
    "Pathway_grouped_heatmap_values.xlsx"
  ),

  overwrite = TRUE
)

# ------------------------------------------------------------
# 17. Completion message
# ------------------------------------------------------------

message("")
message("Heatmap analysis completed successfully.")
message("Number of plotted genes: ", nrow(heatmap_matrix))
message("Results saved in:")
message(normalizePath(output_dir))

message("")
message("Generated files:")
print(
  list.files(
    output_dir,
    full.names = TRUE
  )
)

# Automatically open the result folder on Windows.
if (.Platform$OS.type == "windows") {
  shell.exec(normalizePath(output_dir))
}
