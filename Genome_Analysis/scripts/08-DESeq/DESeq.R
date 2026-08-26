#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# 1. USER SETTINGS
# ------------------------------------------------------------------------------
setwd("E:/Genome analysis/")

count_files <- c(
  "trim_paired_ERR1797969_counts.txt",
  "trim_paired_ERR1797970_counts.txt",
  "trim_paired_ERR1797971_counts.txt",
  "trim_paired_ERR1797972_counts.txt",
  "trim_paired_ERR1797973_counts.txt",
  "trim_paired_ERR1797974_counts.txt"
)

# The order must match count_files exactly.
sample_names <- c(
  "serum_1",
  "serum_2",
  "serum_3",
  "BHI_1",
  "BHI_2",
  "BHI_3"
)

# The order must match count_files and sample_names exactly.
sample_conditions <- c(
  "serum",
  "serum",
  "serum",
  "BHI",
  "BHI",
  "BHI"
)


reference_condition <- "BHI"
comparison_condition <- "serum"

# Path to the Prokka TSV annotation file.
annotation_file <- "E.faecium-annotation.tsv"

# Statistical settings.
adjusted_p_value_threshold <- 0.001
absolute_log2_fold_change_threshold <- 2

# Low-count filtering.
minimum_count <- 10
minimum_samples <- 2

# Plot settings.
number_heatmap_genes <- 30
number_volcano_labels <- 30

# Output directory.
output_dir <- "results/deseq2_annotated"


# ------------------------------------------------------------------------------
# 2. LOAD PACKAGES
# ------------------------------------------------------------------------------

required_packages <- c(
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "ggrepel"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall DESeq2 using BiocManager and the remaining packages ",
    "using install.packages()."
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
})


# ------------------------------------------------------------------------------
# 3. VALIDATE INPUT SETTINGS
# ------------------------------------------------------------------------------

if (length(count_files) != length(sample_names)) {
  stop(
    "The number of count files must equal the number of sample names."
  )
}

if (length(count_files) != length(sample_conditions)) {
  stop(
    "The number of count files must equal the number of conditions."
  )
}

if (anyDuplicated(sample_names)) {
  stop("Sample names must be unique.")
}

missing_count_files <- count_files[!file.exists(count_files)]

if (length(missing_count_files) > 0) {
  stop(
    "These count files were not found:\n",
    paste(missing_count_files, collapse = "\n")
  )
}

if (!file.exists(annotation_file)) {
  stop(
    "The annotation TSV file was not found: ",
    annotation_file
  )
}

valid_conditions <- c(
  reference_condition,
  comparison_condition
)

invalid_conditions <- setdiff(
  unique(sample_conditions),
  valid_conditions
)

if (length(invalid_conditions) > 0) {
  stop(
    "Unrecognized condition labels: ",
    paste(invalid_conditions, collapse = ", "),
    "\nExpected only: ",
    paste(valid_conditions, collapse = ", ")
  )
}


# ------------------------------------------------------------------------------
# 4. CREATE OUTPUT DIRECTORIES
# ------------------------------------------------------------------------------

plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")
object_dir <- file.path(output_dir, "objects")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# 5. FUNCTION TO READ ONE HTSEQ FILE
# ------------------------------------------------------------------------------

read_htseq_file <- function(file_path, sample_name) {
  
  count_data <- read.delim(
    file_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    col.names = c("gene_id", "count")
  )
  
  if (nrow(count_data) == 0) {
    stop("Count file is empty: ", file_path)
  }
  
  if (ncol(count_data) != 2) {
    stop(
      "Expected two columns in HTSeq file: ",
      file_path
    )
  }
  
  count_data$count <- suppressWarnings(
    as.numeric(count_data$count)
  )
  
  if (anyNA(count_data$count)) {
    stop(
      "Nonnumeric count values were found in: ",
      file_path
    )
  }
  
  if (any(count_data$count < 0)) {
    stop(
      "Negative counts were found in: ",
      file_path
    )
  }
  
  if (any(count_data$count %% 1 != 0)) {
    stop(
      "Noninteger counts were found in: ",
      file_path
    )
  }
  
  if (anyDuplicated(count_data$gene_id)) {
    stop(
      "Duplicate gene IDs were found in: ",
      file_path
    )
  }
  
  colnames(count_data)[2] <- sample_name
  
  count_data
}


# ------------------------------------------------------------------------------
# 6. IMPORT AND MERGE THE HTSEQ FILES
# ------------------------------------------------------------------------------

count_tables <- Map(
  read_htseq_file,
  file_path = count_files,
  sample_name = sample_names
)

# Merge by gene ID rather than assuming identical row order.
merged_counts <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "gene_id",
      all = TRUE,
      sort = FALSE
    )
  },
  count_tables
)

if (anyNA(merged_counts)) {
  
  incomplete_rows <- merged_counts[
    !complete.cases(merged_counts),
    ,
    drop = FALSE
  ]
  
  write.table(
    incomplete_rows,
    file = file.path(
      table_dir,
      "genes_missing_from_some_files.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  warning(
    "Some IDs were absent from one or more count files. ",
    "Missing values were replaced with zero."
  )
  
  merged_counts[is.na(merged_counts)] <- 0
}


# ------------------------------------------------------------------------------
# 7. REMOVE HTSEQ SUMMARY ROWS
# ------------------------------------------------------------------------------

htseq_summary_rows <- grepl(
  "^__",
  merged_counts$gene_id
)

if (any(htseq_summary_rows)) {
  
  htseq_summary <- merged_counts[
    htseq_summary_rows,
    ,
    drop = FALSE
  ]
  
  write.table(
    htseq_summary,
    file = file.path(
      table_dir,
      "htseq_summary_counts.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

gene_counts <- merged_counts[
  !htseq_summary_rows,
  ,
  drop = FALSE
]

if (nrow(gene_counts) == 0) {
  stop(
    "No gene counts remain after removing HTSeq summary rows."
  )
}


# ------------------------------------------------------------------------------
# 8. CONSTRUCT THE COUNT MATRIX
# ------------------------------------------------------------------------------

count_matrix <- as.matrix(
  gene_counts[, sample_names, drop = FALSE]
)

rownames(count_matrix) <- gene_counts$gene_id
storage.mode(count_matrix) <- "integer"

if (anyDuplicated(rownames(count_matrix))) {
  stop("The count matrix contains duplicate gene IDs.")
}

write.table(
  data.frame(
    gene_id = rownames(count_matrix),
    count_matrix,
    check.names = FALSE
  ),
  file = file.path(
    table_dir,
    "combined_count_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 9. CHECK THE NUMBER OF ASSIGNED COUNTS
# ------------------------------------------------------------------------------

assigned_counts <- colSums(count_matrix)

counting_summary <- data.frame(
  sample = sample_names,
  condition = sample_conditions,
  assigned_gene_counts = assigned_counts,
  stringsAsFactors = FALSE
)

write.table(
  counting_summary,
  file = file.path(
    table_dir,
    "assigned_count_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(counting_summary)

if (all(assigned_counts == 0)) {
  stop(
    "All assigned gene counts are zero. ",
    "The HTSeq reference and annotation mismatch must be fixed first."
  )
}

if (any(assigned_counts == 0)) {
  stop(
    "At least one sample has zero assigned gene counts."
  )
}


# ------------------------------------------------------------------------------
# 10. CREATE SAMPLE INFORMATION
# ------------------------------------------------------------------------------

sample_information <- data.frame(
  condition = factor(
    sample_conditions,
    levels = c(
      reference_condition,
      comparison_condition
    )
  ),
  row.names = sample_names,
  stringsAsFactors = FALSE
)

if (anyNA(sample_information$condition)) {
  stop(
    "At least one condition became NA. Check the spelling and ",
    "capitalization of sample_conditions."
  )
}

if (!identical(
  rownames(sample_information),
  colnames(count_matrix)
)) {
  stop(
    "The order of sample information does not match the count matrix."
  )
}

write.table(
  data.frame(
    sample = rownames(sample_information),
    sample_information,
    row.names = NULL
  ),
  file = file.path(
    table_dir,
    "sample_information.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("\nSample information:")
print(sample_information)

message("\nSamples per condition:")
print(
  table(
    sample_information$condition,
    useNA = "always"
  )
)


# ------------------------------------------------------------------------------
# 11. READ THE PROKKA ANNOTATION TABLE
# ------------------------------------------------------------------------------

annotation <- read.delim(
  annotation_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

message("\nAnnotation columns:")
print(colnames(annotation))

if (!"locus_tag" %in% colnames(annotation)) {
  stop(
    "The annotation table does not have a 'locus_tag' column.\n",
    "Available columns: ",
    paste(colnames(annotation), collapse = ", ")
  )
}

if (anyDuplicated(annotation$locus_tag)) {
  warning(
    "Duplicate locus tags were found in the annotation. ",
    "Only the first annotation for each locus tag will be used."
  )
  
  annotation <- annotation[
    !duplicated(annotation$locus_tag),
    ,
    drop = FALSE
  ]
}

# Determine how well the annotation matches the HTSeq gene IDs.
annotation_match_count <- sum(
  rownames(count_matrix) %in% annotation$locus_tag
)

annotation_match_percentage <- round(
  100 * annotation_match_count / nrow(count_matrix),
  2
)

message(
  "\nAnnotation matched ",
  annotation_match_count,
  " of ",
  nrow(count_matrix),
  " count-matrix gene IDs (",
  annotation_match_percentage,
  "%)."
)

if (annotation_match_count == 0) {
  stop(
    "None of the count-matrix IDs matched annotation$locus_tag.\n",
    "Inspect head(rownames(count_matrix)) and head(annotation$locus_tag)."
  )
}

if (annotation_match_percentage < 90) {
  warning(
    "Fewer than 90% of count-matrix IDs matched the annotation. ",
    "Check that the TSV and HTSeq annotation came from the same Prokka run."
  )
}


# ------------------------------------------------------------------------------
# 12. PREPARE AN ANNOTATION LOOKUP TABLE
# ------------------------------------------------------------------------------

# Create missing annotation columns when Prokka did not include them.
optional_annotation_columns <- c(
  "ftype",
  "gene",
  "product",
  "EC_number",
  "COG"
)

for (column_name in optional_annotation_columns) {
  if (!column_name %in% colnames(annotation)) {
    annotation[[column_name]] <- NA_character_
  }
}

annotation_lookup <- annotation[
  ,
  c(
    "locus_tag",
    "ftype",
    "gene",
    "product",
    "EC_number",
    "COG"
  ),
  drop = FALSE
]

# Retain only annotations relevant to counted genes.
annotation_lookup <- annotation_lookup[
  annotation_lookup$locus_tag %in% rownames(count_matrix),
  ,
  drop = FALSE
]

write.table(
  annotation_lookup,
  file = file.path(
    table_dir,
    "annotation_used.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 13. FILTER GENES WITH LOW COUNTS
# ------------------------------------------------------------------------------

keep_gene <- rowSums(
  count_matrix >= minimum_count
) >= minimum_samples

filtered_count_matrix <- count_matrix[
  keep_gene,
  ,
  drop = FALSE
]

filtering_summary <- data.frame(
  total_genes = nrow(count_matrix),
  retained_genes = nrow(filtered_count_matrix),
  removed_genes = sum(!keep_gene),
  minimum_count = minimum_count,
  minimum_samples = minimum_samples
)

write.table(
  filtering_summary,
  file = file.path(
    table_dir,
    "filtering_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (nrow(filtered_count_matrix) == 0) {
  stop(
    "No genes remained after low-count filtering."
  )
}

message(
  "\nRetained ",
  nrow(filtered_count_matrix),
  " of ",
  nrow(count_matrix),
  " genes after filtering."
)


# ------------------------------------------------------------------------------
# 14. CONSTRUCT AND RUN THE DESEQ2 MODEL
# ------------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = filtered_count_matrix,
  colData = sample_information,
  design = ~ condition
)

dds <- DESeq(dds)

deseq_results <- results(
  dds,
  contrast = c(
    "condition",
    comparison_condition,
    reference_condition
  ),
  alpha = adjusted_p_value_threshold
)

results_table <- as.data.frame(deseq_results)
results_table$gene_id <- rownames(results_table)


# ------------------------------------------------------------------------------
# 15. ADD ANNOTATION WITHOUT CHANGING RESULT ORDER
# ------------------------------------------------------------------------------

# match() is preferable here because it preserves the exact DESeq2 row order.
annotation_positions <- match(
  results_table$gene_id,
  annotation_lookup$locus_tag
)

results_table$locus_tag <- annotation_lookup$locus_tag[
  annotation_positions
]

results_table$feature_type <- annotation_lookup$ftype[
  annotation_positions
]

results_table$gene_symbol <- annotation_lookup$gene[
  annotation_positions
]

results_table$product <- annotation_lookup$product[
  annotation_positions
]

results_table$EC_number <- annotation_lookup$EC_number[
  annotation_positions
]

results_table$COG <- annotation_lookup$COG[
  annotation_positions
]


# ------------------------------------------------------------------------------
# 16. CLEAN GENE SYMBOLS AND PRODUCT DESCRIPTIONS
# ------------------------------------------------------------------------------

clean_annotation_value <- function(x) {
  
  x <- as.character(x)
  x <- trimws(x)
  
  x[
    x %in% c(
      "",
      "-",
      "NA",
      "na",
      "N/A",
      "unknown"
    )
  ] <- NA_character_
  
  x
}

results_table$gene_symbol <- clean_annotation_value(
  results_table$gene_symbol
)

results_table$product <- clean_annotation_value(
  results_table$product
)

# Use the following display priority:
# 1. Gene symbol
# 2. Product description
# 3. Locus tag
results_table$display_name <- results_table$gene_id

has_product <- !is.na(results_table$product)

results_table$display_name[has_product] <-
  results_table$product[has_product]

has_gene_symbol <- !is.na(results_table$gene_symbol)

results_table$display_name[has_gene_symbol] <-
  results_table$gene_symbol[has_gene_symbol]

# Add the locus tag to short gene symbols so labels remain unambiguous.
results_table$plot_label <- results_table$display_name

results_table$plot_label[has_gene_symbol] <- paste0(
  results_table$gene_symbol[has_gene_symbol],
  "\n(",
  results_table$gene_id[has_gene_symbol],
  ")"
)


# ------------------------------------------------------------------------------
# 17. CLASSIFY DIFFERENTIALLY EXPRESSED GENES
# ------------------------------------------------------------------------------

results_table$classification <- "Not significant"

results_table$classification[
  !is.na(results_table$padj) &
    results_table$padj < adjusted_p_value_threshold &
    results_table$log2FoldChange >=
    absolute_log2_fold_change_threshold
] <- paste("Up in", comparison_condition)

results_table$classification[
  !is.na(results_table$padj) &
    results_table$padj < adjusted_p_value_threshold &
    results_table$log2FoldChange <=
    -absolute_log2_fold_change_threshold
] <- paste("Down in", comparison_condition)

results_table$classification <- factor(
  results_table$classification,
  levels = c(
    paste("Down in", comparison_condition),
    "Not significant",
    paste("Up in", comparison_condition)
  )
)

results_table <- results_table[
  order(results_table$padj, na.last = TRUE),
  ,
  drop = FALSE
]

significant_results <- results_table[
  !is.na(results_table$padj) &
    results_table$padj < adjusted_p_value_threshold &
    abs(results_table$log2FoldChange) >=
    absolute_log2_fold_change_threshold,
  ,
  drop = FALSE
]

upregulated_results <- significant_results[
  significant_results$log2FoldChange >=
    absolute_log2_fold_change_threshold,
  ,
  drop = FALSE
]

downregulated_results <- significant_results[
  significant_results$log2FoldChange <=
    -absolute_log2_fold_change_threshold,
  ,
  drop = FALSE
]


# ------------------------------------------------------------------------------
# 18. EXPORT ANNOTATED RESULT TABLES
# ------------------------------------------------------------------------------

preferred_column_order <- c(
  "gene_id",
  "gene_symbol",
  "product",
  "display_name",
  "feature_type",
  "EC_number",
  "COG",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj",
  "classification",
  "plot_label",
  "locus_tag"
)

preferred_column_order <- intersect(
  preferred_column_order,
  colnames(results_table)
)

results_table <- results_table[
  ,
  preferred_column_order,
  drop = FALSE
]

write.csv(
  results_table,
  file = file.path(
    table_dir,
    "all_deseq2_results_annotated.csv"
  ),
  row.names = FALSE
)

write.csv(
  significant_results,
  file = file.path(
    table_dir,
    "significant_genes_annotated.csv"
  ),
  row.names = FALSE
)

write.csv(
  upregulated_results,
  file = file.path(
    table_dir,
    paste0(
      "upregulated_in_",
      comparison_condition,
      "_annotated.csv"
    )
  ),
  row.names = FALSE
)

write.csv(
  downregulated_results,
  file = file.path(
    table_dir,
    paste0(
      "downregulated_in_",
      comparison_condition,
      "_annotated.csv"
    )
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 19. EXPORT NORMALIZED COUNTS WITH ANNOTATION
# ------------------------------------------------------------------------------

normalized_counts <- counts(
  dds,
  normalized = TRUE
)

normalized_annotation_positions <- match(
  rownames(normalized_counts),
  annotation_lookup$locus_tag
)

normalized_gene_symbols <- clean_annotation_value(
  annotation_lookup$gene[
    normalized_annotation_positions
  ]
)

normalized_products <- clean_annotation_value(
  annotation_lookup$product[
    normalized_annotation_positions
  ]
)

normalized_counts_table <- data.frame(
  gene_id = rownames(normalized_counts),
  gene_symbol = normalized_gene_symbols,
  product = normalized_products,
  normalized_counts,
  check.names = FALSE
)

write.csv(
  normalized_counts_table,
  file = file.path(
    table_dir,
    "normalized_counts_annotated.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 20. VARIANCE-STABILIZING TRANSFORMATION
# ------------------------------------------------------------------------------

vsd <- vst(
  dds,
  blind = FALSE
)

saveRDS(
  dds,
  file = file.path(
    object_dir,
    "deseq2_dataset.rds"
  )
)

saveRDS(
  vsd,
  file = file.path(
    object_dir,
    "variance_stabilized_data.rds"
  )
)


# ------------------------------------------------------------------------------
# 21. PCA PLOT
# ------------------------------------------------------------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = "condition",
  returnData = TRUE
)

percent_variance <- round(
  100 * attr(pca_data, "percentVar")
)

pca_data$sample <- rownames(pca_data)

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = condition,
    label = sample
  )
) +
  geom_point(size = 4) +
  geom_text_repel(
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  labs(
    title = "PCA of RNA-seq samples",
    x = paste0(
      "PC1: ",
      percent_variance[1],
      "% variance"
    ),
    y = paste0(
      "PC2: ",
      percent_variance[2],
      "% variance"
    ),
    color = "Condition"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(
    plot_dir,
    "pca_plot.png"
  ),
  plot = pca_plot,
  width = 7,
  height = 6,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 22. SAMPLE CORRELATION HEATMAP
# ------------------------------------------------------------------------------

sample_annotation <- data.frame(
  Condition = sample_information$condition,
  row.names = rownames(sample_information)
)

sample_correlations <- cor(
  assay(vsd),
  method = "pearson"
)

png(
  filename = file.path(
    plot_dir,
    "sample_correlation_heatmap.png"
  ),
  width = 2200,
  height = 2000,
  res = 300
)

pheatmap(
  sample_correlations,
  annotation_col = sample_annotation,
  annotation_row = sample_annotation,
  main = "RNA-seq sample correlation",
  border_color = NA
)

dev.off()


# ------------------------------------------------------------------------------
# 23. DISPERSION PLOT
# ------------------------------------------------------------------------------

png(
  filename = file.path(
    plot_dir,
    "dispersion_plot.png"
  ),
  width = 2200,
  height = 1800,
  res = 300
)

plotDispEsts(
  dds,
  main = "DESeq2 dispersion estimates"
)

dev.off()


# ------------------------------------------------------------------------------
# 24. ANNOTATED MA PLOT
# ------------------------------------------------------------------------------

ma_data <- results_table

ma_labels <- ma_data[
  !is.na(ma_data$padj) &
    ma_data$classification != "Not significant",
  ,
  drop = FALSE
]

ma_labels <- head(
  ma_labels[
    order(ma_labels$padj),
    ,
    drop = FALSE
  ],
  number_volcano_labels
)

ma_plot <- ggplot(
  ma_data,
  aes(
    x = baseMean,
    y = log2FoldChange,
    color = classification
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.4,
    na.rm = TRUE
  ) +
  geom_hline(
    yintercept = c(
      -absolute_log2_fold_change_threshold,
      absolute_log2_fold_change_threshold
    ),
    linetype = "dashed"
  ) +
  geom_text_repel(
    data = ma_labels,
    aes(label = plot_label),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_x_log10() +
  labs(
    title = paste(
      "MA plot:",
      comparison_condition,
      "versus",
      reference_condition
    ),
    x = "Mean normalized count",
    y = paste0(
      "log2 fold change (",
      comparison_condition,
      " / ",
      reference_condition,
      ")"
    ),
    color = "Classification"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    plot_dir,
    "annotated_ma_plot.png"
  ),
  plot = ma_plot,
  width = 9,
  height = 7,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 25. ANNOTATED VOLCANO PLOT
# ------------------------------------------------------------------------------

volcano_data <- results_table

positive_adjusted_p_values <- volcano_data$padj[
  !is.na(volcano_data$padj) &
    volcano_data$padj > 0
]

if (length(positive_adjusted_p_values) > 0) {
  smallest_adjusted_p_value <- min(
    positive_adjusted_p_values
  )
} else {
  smallest_adjusted_p_value <- .Machine$double.xmin
}

volcano_data$plot_padj <- volcano_data$padj

volcano_data$plot_padj[
  !is.na(volcano_data$plot_padj) &
    volcano_data$plot_padj == 0
] <- smallest_adjusted_p_value / 10

volcano_data$negative_log10_padj <- -log10(
  volcano_data$plot_padj
)

volcano_labels <- volcano_data[
  !is.na(volcano_data$padj) &
    volcano_data$classification != "Not significant",
  ,
  drop = FALSE
]

volcano_labels <- head(
  volcano_labels[
    order(volcano_labels$padj),
    ,
    drop = FALSE
  ],
  number_volcano_labels
)

classification_colors <- c(
  "#377EB8",
  "grey70",
  "#E41A1C"
)

names(classification_colors) <- c(
  paste("Down in", comparison_condition),
  "Not significant",
  paste("Up in", comparison_condition)
)

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = log2FoldChange,
    y = negative_log10_padj,
    color = classification
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.6,
    na.rm = TRUE
  ) +
  geom_vline(
    xintercept = c(
      -absolute_log2_fold_change_threshold,
      absolute_log2_fold_change_threshold
    ),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(
      adjusted_p_value_threshold
    ),
    linetype = "dashed"
  ) +
  geom_text_repel(
    data = volcano_labels,
    aes(label = plot_label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = classification_colors,
    drop = FALSE
  ) +
  labs(
    title = paste(
      "Differential expression:",
      comparison_condition,
      "versus",
      reference_condition
    ),
    subtitle = paste0(
      "Adjusted p-value < ",
      adjusted_p_value_threshold,
      " and |log2 fold change| >= ",
      absolute_log2_fold_change_threshold
    ),
    x = paste0(
      "log2 fold change (",
      comparison_condition,
      " / ",
      reference_condition,
      ")"
    ),
    y = "-log10 adjusted p-value",
    color = "Classification"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(
    plot_dir,
    "annotated_volcano_plot.png"
  ),
  plot = volcano_plot,
  width = 10,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------------------------
# 26. FUNCTION FOR CREATING UNIQUE HEATMAP LABELS
# ------------------------------------------------------------------------------

get_heatmap_labels <- function(gene_ids, annotated_results) {
  
  result_positions <- match(
    gene_ids,
    annotated_results$gene_id
  )
  
  labels <- annotated_results$display_name[
    result_positions
  ]
  
  labels[
    is.na(labels) |
      trimws(labels) == ""
  ] <- gene_ids[
    is.na(labels) |
      trimws(labels) == ""
  ]
  
  # Different genes can share the same symbol.
  # make.unique prevents duplicate heatmap row names.
  make.unique(
    labels,
    sep = " [duplicate "
  )
}


# ------------------------------------------------------------------------------
# 27. HEATMAP OF TOP VARIABLE GENES
# ------------------------------------------------------------------------------

gene_variances <- apply(
  assay(vsd),
  1,
  var
)

number_to_plot <- min(
  number_heatmap_genes,
  length(gene_variances)
)

top_variable_gene_ids <- names(
  sort(
    gene_variances,
    decreasing = TRUE
  )
)[seq_len(number_to_plot)]

variable_heatmap_matrix <- assay(vsd)[
  top_variable_gene_ids,
  ,
  drop = FALSE
]

variable_heatmap_matrix <- t(
  scale(
    t(variable_heatmap_matrix)
  )
)

valid_rows <- apply(
  variable_heatmap_matrix,
  1,
  function(x) all(is.finite(x))
)

variable_heatmap_matrix <- variable_heatmap_matrix[
  valid_rows,
  ,
  drop = FALSE
]

if (nrow(variable_heatmap_matrix) >= 2) {
  
  rownames(variable_heatmap_matrix) <- get_heatmap_labels(
    rownames(variable_heatmap_matrix),
    results_table
  )
  
  png(
    filename = file.path(
      plot_dir,
      "annotated_top_variable_genes_heatmap.png"
    ),
    width = 2600,
    height = 3000,
    res = 300
  )
  
  pheatmap(
    variable_heatmap_matrix,
    annotation_col = sample_annotation,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    fontsize_row = 7,
    main = paste(
      "Top",
      nrow(variable_heatmap_matrix),
      "variable genes"
    )
  )
  
  dev.off()
}


# ------------------------------------------------------------------------------
# 28. HEATMAP OF TOP SIGNIFICANT GENES
# ------------------------------------------------------------------------------

if (nrow(significant_results) >= 2) {
  
  significant_gene_number <- min(
    number_heatmap_genes,
    nrow(significant_results)
  )
  
  top_significant_gene_ids <- head(
    significant_results$gene_id,
    significant_gene_number
  )
  
  significant_heatmap_matrix <- assay(vsd)[
    top_significant_gene_ids,
    ,
    drop = FALSE
  ]
  
  significant_heatmap_matrix <- t(
    scale(
      t(significant_heatmap_matrix)
    )
  )
  
  valid_rows <- apply(
    significant_heatmap_matrix,
    1,
    function(x) all(is.finite(x))
  )
  
  significant_heatmap_matrix <- significant_heatmap_matrix[
    valid_rows,
    ,
    drop = FALSE
  ]
  
  if (nrow(significant_heatmap_matrix) >= 2) {
    
    rownames(significant_heatmap_matrix) <- get_heatmap_labels(
      rownames(significant_heatmap_matrix),
      results_table
    )
    
    png(
      filename = file.path(
        plot_dir,
        "annotated_significant_genes_heatmap.png"
      ),
      width = 2600,
      height = 3000,
      res = 300
    )
    
    pheatmap(
      significant_heatmap_matrix,
      annotation_col = sample_annotation,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      border_color = NA,
      fontsize_row = 7,
      main = "Top differentially expressed genes"
    )
    
    dev.off()
  }
} else {
  message(
    "Fewer than two significant genes were found. ",
    "The significant-gene heatmap was not generated."
  )
}


# ------------------------------------------------------------------------------
# 29. NORMALIZED COUNTS FOR TOP SIGNIFICANT GENES
# ------------------------------------------------------------------------------

if (nrow(significant_results) > 0) {
  
  number_count_plot_genes <- min(
    12,
    nrow(significant_results)
  )
  
  selected_gene_ids <- head(
    significant_results$gene_id,
    number_count_plot_genes
  )
  
  count_plot_tables <- lapply(
    selected_gene_ids,
    function(gene_id) {
      
      gene_data <- plotCounts(
        dds,
        gene = gene_id,
        intgroup = "condition",
        returnData = TRUE
      )
      
      result_position <- match(
        gene_id,
        results_table$gene_id
      )
      
      gene_data$gene_id <- gene_id
      gene_data$display_name <-
        results_table$display_name[result_position]
      gene_data$sample <- rownames(gene_data)
      
      gene_data
    }
  )
  
  count_plot_data <- do.call(
    rbind,
    count_plot_tables
  )
  
  # Prevent duplicate facet labels.
  label_lookup <- unique(
    count_plot_data[
      ,
      c("gene_id", "display_name")
    ]
  )
  
  label_lookup$facet_label <- make.unique(
    label_lookup$display_name,
    sep = " [duplicate "
  )
  
  facet_positions <- match(
    count_plot_data$gene_id,
    label_lookup$gene_id
  )
  
  count_plot_data$facet_label <-
    label_lookup$facet_label[facet_positions]
  
  normalized_count_plot <- ggplot(
    count_plot_data,
    aes(
      x = condition,
      y = count,
      color = condition
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.2
    ) +
    geom_jitter(
      width = 0.12,
      size = 2.3
    ) +
    facet_wrap(
      ~ facet_label,
      scales = "free_y"
    ) +
    scale_y_log10() +
    labs(
      title = "Normalized counts for top significant genes",
      x = "Condition",
      y = "Normalized count, log10 scale",
      color = "Condition"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 8)
    )
  
  ggsave(
    filename = file.path(
      plot_dir,
      "annotated_normalized_count_plots.png"
    ),
    plot = normalized_count_plot,
    width = 12,
    height = 9,
    dpi = 300
  )
}


# ------------------------------------------------------------------------------
# 30. ANNOTATION AND ANALYSIS SUMMARY
# ------------------------------------------------------------------------------

genes_with_symbol <- sum(
  !is.na(results_table$gene_symbol)
)

genes_with_product <- sum(
  !is.na(results_table$product)
)

analysis_summary <- data.frame(
  comparison = paste(
    comparison_condition,
    "versus",
    reference_condition
  ),
  total_counted_genes = nrow(count_matrix),
  genes_after_filtering = nrow(filtered_count_matrix),
  annotation_matches = annotation_match_count,
  annotation_match_percentage = annotation_match_percentage,
  genes_with_gene_symbol = genes_with_symbol,
  genes_with_product_description = genes_with_product,
  significant_genes = nrow(significant_results),
  upregulated_in_comparison = nrow(upregulated_results),
  downregulated_in_comparison = nrow(downregulated_results),
  adjusted_p_value_threshold =
    adjusted_p_value_threshold,
  absolute_log2_fold_change_threshold =
    absolute_log2_fold_change_threshold
)

write.table(
  analysis_summary,
  file = file.path(
    table_dir,
    "analysis_summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("\nAnalysis summary:")
print(analysis_summary)


