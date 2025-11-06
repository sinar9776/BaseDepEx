## ============================================================
## 1_DEAnalysis.R — Modular Differential Expression  (gene/isoform)
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ggrepel)
})

set.seed(1)

## ---------------- Parameters ----------------
padj_cutoff <- 0.01
lfc_cutoff  <- 0
shrink_type <- "apeglm"  

## ---------------- Paths ----------------
base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
res_dir  <- file.path(base_dir, "results/1_DEAnalysis")
tbl_dir  <- file.path(res_dir, "tables/raw")
fig_dir  <- file.path(res_dir, "figures/raw")
shr_tbl  <- file.path(res_dir, "tables/shrunk")
shr_fig  <- file.path(res_dir, "figures/shrunk")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(shr_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(shr_fig, recursive = TRUE, showWarnings = FALSE)

## ============================================================
## Load data automatically based on analysis level
## ============================================================

# Choose: "gene" or "isoform"
analysis_level <- "gene"

metadata <- read.delim(file.path(data_dir, "metadata.txt"), check.names = FALSE)
stopifnot(all(c("SampleID", "group") %in% colnames(metadata)))

annotation <- read.delim(file.path(data_dir, "annotation.txt"), check.names = FALSE)

if (analysis_level == "isoform") {
  df <- readr::read_tsv(file.path(data_dir, "isoform_counts.tsv"), show_col_types = FALSE) %>% as.data.frame()
  id_col <- "isoform_id"
} else if (analysis_level == "gene") {
  df <- readr::read_tsv(file.path(data_dir, "gene_counts.tsv"), show_col_types = FALSE) %>% as.data.frame()
  id_col <- "gene_id"
} else {
  stop("analysis_level must be either 'gene' or 'isoform'")
}

rownames(df) <- df[[1]]
count_mat <- df[, intersect(colnames(df), metadata$SampleID), drop = FALSE]
metadata  <- metadata[match(colnames(count_mat), metadata$SampleID), ]
stopifnot(all(colnames(count_mat) == metadata$SampleID))

# Annotation
sym_col <- if ("symbol" %in% names(annotation)) "symbol" else if ("gene_name" %in% names(annotation)) "gene_name" else NA
ann_keep_cols <- unique(c("gene_id", id_col, sym_col))
ann_keep_cols <- ann_keep_cols[ann_keep_cols %in% names(annotation)]
annotation2 <- annotation[, ann_keep_cols, drop = FALSE]
if (!is.na(sym_col)) names(annotation2)[names(annotation2) == sym_col] <- "SYMBOL"

if (analysis_level == "gene") {
  # Deduplicate per gene_id for gene-level
  annotation2 <- annotation2 %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(SYMBOL = dplyr::first(na.omit(SYMBOL)), .groups = "drop")
}

if (!"gene_id" %in% colnames(annotation2)) {
  stop("Annotation file must include a 'gene_id' column for enrichment compatibility.")
}

## ---------------- Filter low counts ----------------
keep <- rowSums(count_mat >= 10) >= min(table(metadata$group))
count_mat <- count_mat[keep, , drop = FALSE]

dds1 <- DESeqDataSetFromMatrix(round(as.matrix(count_mat)), metadata, design = ~ group)
dds1$group <- relevel(dds1$group, ref = "DMSO")
dds2 <- dds1
if ("DMSO_OHT" %in% levels(dds2$group)) dds2$group <- relevel(dds2$group, ref = "DMSO_OHT")

dds1 <- DESeq(dds1, parallel = TRUE)
dds2 <- DESeq(dds2, parallel = TRUE)

## ============================================================
## Part 1 — DE calculation 
## ============================================================

counts_summary_raw <- data.frame()
for (g in levels(dds1$group)) {
  # determine reference
  ref <- if (grepl("_OHT$", g)) "DMSO_OHT" else "DMSO"
  
  # special baseline comparison
  if (g == "DMSO_OHT") ref <- "DMSO"
  
  if (g == ref || !(ref %in% levels(dds1$group))) next
  nm <- if (g == "DMSO_OHT") "DMSO_OHT_vs_DMSO" else g
  message(" [RAW] ", nm)
  
  res_raw <- results(dds1, contrast = c("group", g, ref))  
  res_df <- as.data.frame(res_raw)
  res_df[[id_col]] <- sub("\\..*$", "", rownames(res_df))  # strip version
  res_df <- merge(res_df, annotation2, by = id_col, all.x = TRUE, sort = FALSE)
  
  nc1 <- counts(dds1, normalized = TRUE)
  idx_t <- which(metadata$group == g)
  idx_c <- which(metadata$group == ref)
  res_df$avg_treated <- rowMeans(nc1[, idx_t, drop = FALSE], na.rm = TRUE)
  res_df$avg_control <- rowMeans(nc1[, idx_c, drop = FALSE], na.rm = TRUE)
  
  # Significance for raw
  res_df$sig_raw <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff, "up",
                           ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))
  
  # Save
  write.table(res_df, file.path(tbl_dir, paste0("tT_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_raw == "up"),   file.path(tbl_dir, paste0("up_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_raw == "down"), file.path(tbl_dir, paste0("down_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  counts_summary_raw <- rbind(counts_summary_raw, data.frame(
    contrast = nm,
    n_up     = sum(res_df$sig_raw == "up"),
    n_down   = sum(res_df$sig_raw == "down")
  ))
}
write.table(counts_summary_raw, file.path(tbl_dir, "DE_counts_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ============================================================
## Part 2 — Plots 
## ============================================================

# Dispersion plots
png(file.path(fig_dir, "Dispersion_dds1.png"), width = 1600, height = 1200, res = 150)
plotDispEsts(dds1)
dev.off()

if (exists("dds2")) {
  png(file.path(fig_dir, "Dispersion_dds2.png"), width = 1600, height = 1200, res = 150)
  plotDispEsts(dds2)
  dev.off()
}


## PCA plots per comparison
# Define color scheme for PCA (treated = red, control = blue)
pca_colors <- c("treated" = "#D62728", "control" = "#1F77B4")

# Loop through same contrasts as DESeq2 comparisons
for (g in levels(dds1$group)) {
  ref <- if (grepl("_OHT$", g)) "DMSO_OHT" else "DMSO"
  if (g == ref || !(ref %in% levels(dds1$group))) next
  nm <- g
  message(" [PCA] ", nm, " vs ", ref)
  
  # Subset metadata and counts for the two groups
  keep_samples <- metadata$group %in% c(ref, g)
  sub_metadata <- metadata[keep_samples, ]
  sub_counts <- count_mat[, sub_metadata$SampleID, drop = FALSE]
  
  # Make DESeq object and VST transform
  dds_sub <- DESeqDataSetFromMatrix(round(sub_counts), sub_metadata, design = ~ group)
  dds_sub$group <- droplevels(dds_sub$group)
  vsd_sub <- vst(dds_sub, blind = TRUE)
  
  # PCA data
  pca_data <- plotPCA(vsd_sub, intgroup = "group", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  # Rename groups for legend clarity
  pca_data$group <- factor(pca_data$group, levels = c(ref, g))
  
  # Assign colors consistently (blue = control, red = treated)
  color_map <- setNames(c(pca_colors["control"], pca_colors["treated"]), c(ref, g))
  
  # Plot PCA
  p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = color_map,
                       name = "Condition",
                       labels = c(ref, g)) +
    labs(x = paste0("PC1 (", percentVar[1], "%)"),
         y = paste0("PC2 (", percentVar[2], "%)"),
         title = paste("PCA:", g, "vs", ref)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "right")
  
  # Save PCA plot
  ggsave(file.path(fig_dir, paste0("PCA_", g, "_vs_", ref, ".png")),
         p_pca, width = 6, height = 5, dpi = 300)
}

## ---------- MA & Volcano  ----------
col_scale <- c(up = "#D62728", down = "#1F77B4", ns = "gray80")
files <- list.files(tbl_dir, pattern = "^tT_.*\\.tsv$", full.names = TRUE)

for (fp in files) {
  nm_raw <- tools::file_path_sans_ext(basename(fp))
  nm <- sub("^tT_", "", nm_raw)
  
  res_df <- tryCatch(read.delim(fp, check.names = FALSE), error = function(e) NULL)
  if (is.null(res_df) || nrow(res_df) == 0) {
    message(" [skip] ", nm, " — file empty or unreadable.")
    next
  }
  
  # Ensure required columns exist
  needed <- c("baseMean", "log2FoldChange", "padj")
  if (!all(needed %in% names(res_df))) {
    message(" [skip] ", nm, " — missing required columns: ", paste(setdiff(needed, names(res_df)), collapse = ", "))
    next
  }
  
  # If sig_raw is missing, compute it (fix for your crash)
  if (!"sig_raw" %in% names(res_df)) {
    res_df$sig_raw <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff, "up",
                             ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))
  }
  
  # Build plotting category and order so red/blue sit on top of gray
  res_df$cat <- factor(res_df$sig_raw, levels = c("ns", "down", "up"))
  plot_df <- res_df %>% arrange(cat)  # ns first, then down, then up
  
  ## --- MA plot (with legend) ---
  p_ma <- ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = cat)) +
    geom_point(size = 1.2, alpha = 0.85) +
    scale_color_manual(values = col_scale,
                       breaks = c("up","down","ns"),
                       labels = c("Up","Down","NS"),
                       name = "DE category") +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    labs(title = paste("MA:", nm), x = "log10(baseMean + 1)", y = "log2(Fold Change)") +
    theme_bw(10)
  
  ## --- Volcano plot (with legend) ---
  p_vol <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj), color = cat)) +
    geom_point(size = 1.2, alpha = 0.85) +
    scale_color_manual(values = col_scale,
                       breaks = c("up","down","ns"),
                       labels = c("Up","Down","NS"),
                       name = "DE category") +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(title = paste("Volcano:", nm), x = "log2(Fold Change)", y = "-log10(padj)") +
    theme_bw(10)
  
  ggsave(file.path(fig_dir, paste0("MA_", nm, ".png")), p_ma,  width = 6,  height = 5, dpi = 300)
  ggsave(file.path(fig_dir, paste0("Volcano_", nm, ".png")), p_vol, width = 6.5, height = 5, dpi = 300)
}

## ============================================================
## Part 3 — Shrinkage-based DE calculation 
## ============================================================

counts_summary_shr <- data.frame()
for (g in levels(dds1$group)) {
  # define reference
  ref <- if (grepl("_OHT$", g)) "DMSO_OHT" else "DMSO"
  
  # special baseline comparison
  if (g == "DMSO_OHT") ref <- "DMSO"
  
  # skip invalid contrasts
  if (g == ref || !(ref %in% levels(dds1$group))) next
  
  # set readable contrast name
  nm <- if (g == "DMSO_OHT") "DMSO_OHT_vs_DMSO" else g  
  message(" [SHRUNKEN] ", nm)
  
  dds_use <- if (ref == "DMSO_OHT") dds2 else dds1
  coef_name <- paste0("group_", g, "_vs_", ref)
  rn <- resultsNames(dds_use)
  if (!coef_name %in% rn) next
  
  res_raw <- results(dds_use, contrast = c("group", g, ref))
  res_shr <- tryCatch(
    lfcShrink(dds_use, coef = coef_name, res = res_raw, type = shrink_type),
    error = function(e) {
      message("apeglm failed, using 'normal'.")
      lfcShrink(dds_use, coef = coef_name, res = res_raw, type = "normal")
    }
  )
  
  res_df <- as.data.frame(res_shr)
  res_df[[id_col]] <- sub("\\..*$", "", rownames(res_df))
  res_df <- merge(res_df, annotation2, by = id_col, all.x = TRUE, sort = FALSE)
  
  nc <- counts(dds_use, normalized = TRUE)
  idx_t <- which(metadata$group == g)
  idx_c <- which(metadata$group == ref)
  res_df$avg_treated <- rowMeans(nc[, idx_t, drop = FALSE], na.rm = TRUE)
  res_df$avg_control <- rowMeans(nc[, idx_c, drop = FALSE], na.rm = TRUE)
  
  res_df$sig_shr <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff &
                             res_df$log2FoldChange >  lfc_cutoff, "up",
                           ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff &
                                    res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))
  
  write.table(res_df, file.path(shr_tbl, paste0("tT_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_shr == "up"),   file.path(shr_tbl, paste0("up_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_shr == "down"), file.path(shr_tbl, paste0("down_", nm, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  counts_summary_shr <- rbind(counts_summary_shr, data.frame(
    contrast = nm,
    n_up_shr   = sum(res_df$sig_shr == "up"),
    n_down_shr = sum(res_df$sig_shr == "down")
  ))
}
write.table(counts_summary_shr, file.path(shr_tbl, "DE_counts_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


## ============================================================
## Part 4 — Shrunken plots
## ============================================================

col_scale <- c(up = "#D62728", down = "#1F77B4", ns = "gray80")

files_shr <- list.files(shr_tbl, pattern = "^tT_.*\\.tsv$", full.names = TRUE)

for (fp in files_shr) {
  nm_raw <- tools::file_path_sans_ext(basename(fp))
  nm <- sub("^tT_", "", nm_raw)
  
  res_df <- tryCatch(read.delim(fp, check.names = FALSE), error = function(e) NULL)
  if (is.null(res_df) || nrow(res_df) == 0) {
    message(" [skip] ", nm, " — empty or unreadable file.")
    next
  }
  
  if (!"sig_shr" %in% colnames(res_df)) next
  res_df$cat <- factor(res_df$sig_shr, levels = c("ns", "down", "up"))
  plot_df <- res_df %>% arrange(cat)
  
  ## --- MA plot ---
  p_ma <- ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = cat)) +
    geom_point(size = 1.2, alpha = 0.85) +
    scale_color_manual(values = col_scale,
                       breaks = c("up","down","ns"),
                       labels = c("Up","Down","NS"),
                       name = "DE category") +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    labs(title = paste("MA (shrunk):", nm),
         x = "log10(baseMean + 1)", y = "log2(Fold Change)") +
    theme_bw(10)
  
  ## --- Volcano plot ---
  p_vol <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj), color = cat)) +
    geom_point(size = 1.2, alpha = 0.85) +
    scale_color_manual(values = col_scale,
                       breaks = c("up","down","ns"),
                       labels = c("Up","Down","NS"),
                       name = "DE category") +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(title = paste("Volcano (shrunk):", nm),
         x = "log2(Fold Change)", y = "-log10(adjusted p-value)") +
    theme_bw(10)
  
  ## --- Save all ---
  ggsave(file.path(shr_fig, paste0("MA_", nm, ".png")), p_ma,  width = 6, height = 5, dpi = 300)
  ggsave(file.path(shr_fig, paste0("Volcano_", nm, ".png")), p_vol, width = 6.5, height = 5, dpi = 300)
}

