# Illumina 16S ASVs Microbiome Analysis pipeline
### This is the latest backup for Midgut - FTA comparisson analysis (August 2024)
From raw reads (paired end Illumina amplicon) to microbiome analysis, including the following publication plots:
   - Rarefaction Curves
   - Alpha diversity (Observed, Shannon, Simpson, Pielou) and boxplots with paired wised Wilcoxon test (Benjamini Huchberg p-value adjustment)
   - Beta diversity (Bray-Curtis dissimilarity distance) and NMDS (stress score) and PCoA (Beta disperssion and PERMANOVA)
   - Bar Plots of top 50 most abundant taxa, with LCBD value included (beta diversity showing grouping influence of each sample within categories
   - Differential abundance analysis (Log2Fold change, DeSeq2)
