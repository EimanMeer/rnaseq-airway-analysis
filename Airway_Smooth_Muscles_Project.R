################################################################################
# Project: Differential Expression in Airway Smooth Muscle Cells
# Description: Identifying genes regulated by Dexamethasone treatment
# Data Source: Himes et al. (2014) via 'airway' Bioconductor package
################################################################################

# 1. SETUP & REQUISITES --------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

packages <- c("DESeq2", "airway", "org.Hs.eg.db", "AnnotationDbi", 
              "EnhancedVolcano", "ggplot2", "pheatmap", "RColorBrewer", "dplyr")

# Install missing packages
for (pkg in packages) {
  if (!require(pkg, character_only = TRUE)) {
    if (pkg %in% c("ggplot2", "pheatmap", "RColorBrewer", "dplyr")) {
      install.packages(pkg)
    } else {
      BiocManager::install(pkg)
    }
  }
}

# 2. DIRECTORY STRUCTURE -------------------------------------------------------
project_path <- "E:/Airway_Smooth_Muscles_Project"
if (!dir.exists(project_path)) dir.create(project_path, recursive = TRUE)
setwd(project_path)

folders <- c("data", "results", "figures", "scripts", "report")
invisible(lapply(folders, function(f) if(!dir.exists(f)) dir.create(f)))

# 3. DATA LOADING & PRE-FILTERING ----------------------------------------------
data("airway")
# Pre-filter: keep genes with total counts >= 10
keep <- rowSums(assay(airway)) >= 10
airway_filtered <- airway[keep, ]

# 4. EXPLORATORY QC (Library Sizes) --------------------------------------------
lib_sizes <- colSums(assay(airway_filtered))
lib_df <- data.frame(sample = names(lib_sizes), total = lib_sizes, group = colData(airway_filtered)$dex)

ggplot(lib_df, aes(x = sample, y = total / 1e6, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("trt" = "#e07b54", "untrt" = "#5b8db8")) +
  labs(title = "Library sizes per sample", x = "Sample", y = "Total reads (millions)") +
  theme_minimal()
ggsave("figures/01_library_sizes.png", width = 8, height = 5)

# 5. DESeq2 PIPELINE -----------------------------------------------------------
dds <- DESeqDataSet(airway_filtered, design = ~ dex)
dds$dex <- relevel(dds$dex, ref = "untrt") # untreated as baseline
dds <- DESeq(dds)

# VST Transformation for visualization
vsd <- vst(dds, blind = FALSE)

# 6. SAMPLE CORRELATION & PCA --------------------------------------------------
# Correlation Heatmap
sample_cor <- cor(assay(vsd))
anno_df <- data.frame(Treatment = colData(vsd)$dex, row.names = colnames(vsd))
pheatmap(sample_cor, annotation_col = anno_df, 
         color = colorRampPalette(c("white", "#2166ac"))(100),
         filename = "figures/02_sample_correlation.png", width = 7, height = 6)

# PCA Plot
pca_data <- plotPCA(vsd, intgroup = c("dex", "cell"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x=PC1, y=PC2, colour=dex, shape=cell)) +
  geom_point(size = 5, alpha = 0.8) +
  scale_colour_manual(values = c("trt" = "#e07b54", "untrt" = "#5b8db8")) +
  labs(title = "PCA (VST normalised)", x = paste0("PC1: ", percent_var[1], "%"), 
       y = paste0("PC2: ", percent_var[2], "%")) +
  theme_minimal()
ggsave("figures/03_PCA.png", width=8, height=6)

# 7. RESULTS & ANNOTATION ------------------------------------------------------
res <- results(dds, contrast = c("dex", "trt", "untrt"), alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = "dex_trt_vs_untrt", type = "apeglm")

res_df <- as.data.frame(res_shrunk)
res_df$ensembl_id <- rownames(res_df)
res_df$gene_name <- mapIds(org.Hs.eg.db, keys = res_df$ensembl_id,
                           column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

sig_genes <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)

write.csv(res_df, "results/all_genes_results.csv", row.names = FALSE)
write.csv(sig_genes, "results/significant_DEGs.csv", row.names = FALSE)

# 8. FINAL VISUALIZATIONS ------------------------------------------------------
# Volcano Plot
res_df$label <- ifelse(is.na(res_df$gene_name), res_df$ensembl_id, res_df$gene_name)
volcano <- EnhancedVolcano(res_df, lab = res_df$label, x = "log2FoldChange", y = "padj",
                           pCutoff = 0.05, FCcutoff = 1.0, title = "Volcano Plot")
ggsave("figures/04_volcano.png", volcano, width=10, height=8)

# Top 50 Heatmap
top50_ids <- sig_genes %>% arrange(padj) %>% head(50) %>% pull(ensembl_id)
heatmap_mat <- assay(vsd)[top50_ids, ]
rownames(heatmap_mat) <- sig_genes$gene_name[match(top50_ids, sig_genes$ensembl_id)]
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))

pheatmap(heatmap_mat_scaled, annotation_col = anno_df,
         color = colorRampPalette(c("#2166ac","white","#d6604d"))(100),
         filename = "figures/05_top50_heatmap.png", width = 7, height = 9)

# 9. BACKUP WORKSPACE ----------------------------------------------------------
save.image("data/Project_Backup.RData")
message("--- Analysis Complete ---")