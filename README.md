# Transcriptomic Analysis: Airway Smooth Muscle Cells

This project analyzes RNA-seq data from human airway smooth muscle cells to identify genes differentially expressed following **Dexamethasone** (steroid) treatment.

## 📊 Summary of Findings
- **Organism:** *Homo sapiens*
- **Experimental Design:** 8 samples (4 Treated vs 4 Untreated Controls)
- **Primary Tool:** DESeq2
- **Key Results:** PCA shows strong separation of treatment groups on PC1. Significant genes (padj < 0.05, |LFC| > 1) were identified and categorized into up- and down-regulated pathways.

## 📂 Project Structure
- `/data`: Contains raw count previews and RData workspace backups.
- `/scripts`: The full R pipeline (`analysis.R`).
- `/figures`: QC plots (Library sizes), Sample Correlations, PCA, Volcano Plot, and Expression Heatmaps.
- `/results`: Full CSV tables of differential expression results.
- `/report`: RMarkdown template for automated HTML reports.

## 🛠️ Requirements
This project requires **R 4.0+** and the following Bioconductor packages:
- `DESeq2`, `airway`, `org.Hs.eg.db`, `EnhancedVolcano`, `pheatmap`.

## 🚀 How to Run
1. Clone this repository.
2. Open `scripts/analysis.R` in RStudio.
3. Ensure your working directory  and file names are added to the script correctly.
4. Run the script to regenerate all figures and results.
