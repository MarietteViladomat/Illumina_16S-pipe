############################################################################
##############                                                   ##############
##############                                                   ##############
##############             MARIETTE VILADOMAT JASSO              ##############
##############                                                   ##############
##############                                                   ##############
##############   16S rRNA amplicon Microbiome Analysis  with R   ##############
##############                                                   ##############
##############               Mosquito Microbiome                 ##############
##############                  Midgut & FTA                     ##############
##############                                                   ##############
##############                                                   ##############
############################################################################

# Personal R Script for analysing 16S Illumina Raw reads
# for microbiome comparisson between 2 or more categorical groups

# Analysis includes
# a) Adapter trimming and min length of 100pb filtering, using Rfastp
# b) Quality trimming of min Q 20
# c) DADA2 pre-processing: learn error rates of F & R reads, derreplication,
  # sample error inference, merge paired end reads, sequence table construction, 
  # removing chimeras & making a CHECKPOINT OF reads without chimeras and generating DADA 2 REPORT
# d) Taxonomic assignation using SILVA train set and species database. (with DADA2)
# e) Generating phyloseq object (ps)
# f) Data normalization through 2 possible methods: DeSeq2 and TMM.
# g) Generating new normalized phyloseq object (ps_norm_deseq) or (ps_norm_tmm)
# h) Statistical analysis and figure making of (using vegan):
      # Rarefaction curves
      # Alpha Diversity (Observed, Shannon, Simpson, Pielou) and boxplots (Willcoxon with Benjamini Huchberg p-value correction)
      # Beta Diversity (Bray-Curtis distance) and NMDS (stress score) and PCoA (beta dispertion and PERMANOVA) plots
      # Bar Plots that include de LCBD value (beta: grouping influence of each sample)
      # Differential Abundance Plot


###-------STEP 0: LOAD LIBRARIES -----
library(dplyr)
library(Rfastp)
library(ggplot2)
library(dada2)
library(stats)

###-------STEP 1: IMPORT DATA-------------

# Folder with input data
file.path="RealData/"
filepath = "RealData/unzipped_fastqs/" 

# Extract sample IDs
metadata <- read.table(paste0(file.path,"metadata.txt"), header = T)
sample.names <- sapply(strsplit(basename(metadata$sample_id), "-"), `[`, 1)

# Identify forward and reverse reads
fnFs <- sort(list.files(filepath, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(filepath, pattern="_R2.fastq", full.names = TRUE))



# Check to ensure the order matches the sample names
if(length(fnFs) != length(fnRs)) {
  stop("Mismatched number of forward and reverse reads")
}else{
  print("All samples have paired reads.")
}


###-------STEP 2.1: DATA PREPROCESSING: barcode cut and filter with fastp -----------

# Define output directory
output_dir <- "RealData/processed"

# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Iterate over each file pair and run rfastp
for (i in seq_along(fnFs)) {
  in_filtF <- fnFs[i]
  in_filtR <- fnRs[i]
  
  # Use sub to remove the suffix "_F_filt" or "_R_filt" to get the core sample name
  core_name <- sub("_R[12]\\.fastq", "", basename(fnFs[i]))
  
  # Construct the output prefix by combining output directory and core name
  output_prefix <- file.path(output_dir, paste0(core_name, "_processed"))
  
  # Run rfastp
  rfastp(
    read1 = in_filtF,
    read2 = in_filtR,
    outputFastq = output_prefix,
    adapterTrimming = TRUE,
    qualityFiltering = TRUE,
    verbose = TRUE,
    minReadLength = 100,
    averageQualFilter = 30
  )
  
  # Optionally, print progress
  cat("Processed:", in_filtF, "and", in_filtR, "\n")
}




###-------STEP 2.2: DATA PREPROCESSING: filters, q reports and calculate error rates -----------


# Identify forward and reverse reads
fnFs_processed <- sort(list.files("RealData/processed/", pattern="_R1.fastq", full.names = TRUE))
fnRs_processed <- sort(list.files("RealData/processed/", pattern="_R2.fastq", full.names = TRUE))


#set output directory for trimmed data
outdir_trimmed <- file.path(file.path, "filtered")
dir.create(outdir_trimmed, showWarnings = FALSE)

filtFs <- file.path(outdir_trimmed, paste0(unique(sample.names), "_F_filt.fastq"))
filtRs <- file.path(outdir_trimmed, paste0(unique(sample.names), "_R_filt.fastq"))
names(filtFs) <- unique(sample.names)
names(filtRs) <- unique(sample.names)

# check raw fastq quality
raw_F_QC <- plotQualityProfile(fnFs_processed)
ggsave(paste0(outdir_trimmed,"/fastp_F_QC.pdf"), raw_F_QC, dpi = 300, height = 10, width = 10)

raw_R_QC <- plotQualityProfile(fnRs_processed)
ggsave(paste0(outdir_trimmed,"/fastp_R_QC.pdf"), raw_F_QC, dpi = 300, height = 10, width = 10)


# filter/trim
filtered_report <- filterAndTrim(fnFs_processed, filtFs, fnRs_processed, filtRs, truncLen = 0,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimRight = 0, trimLeft = 0, 
                     compress=TRUE, multithread=FALSE, verbose = TRUE,) #WHEN IN LINUX, multithread = TRUE!!!

#format and output report
filtered_report <- as.data.frame(filtered_report)
filtered_report$sample_id <- rownames(filtered_report)
rownames(filtered_report) <- NULL
filtered_report$sample_id <- sub("-.*", "", filtered_report$sample_id)
filtered_report <- filtered_report[,c("sample_id", "reads.in", "reads.out")]
filtered_report$percentage_reads_kept <-  round((filtered_report$reads.out/filtered_report$reads.in* 100), 2)


# check filtered fastq quality
filt_F_QC <- plotQualityProfile(filtFs)
ggsave(paste0(outdir_trimmed,"/filt_F_QC.pdf"), filt_F_QC, dpi = 300, height = 10, width = 10)

filt_R_QC <- plotQualityProfile(filtRs)
ggsave(paste0(outdir_trimmed,"/filt_R_QC.pdf"), filt_F_QC, dpi = 300, height = 10, width = 10)


#learn error rates for F and R
#(paths to filtered reads: filtFs and filtRs)
set.seed(420)

err_F = learnErrors(filtFs, randomize=TRUE) #WHEN IN LINUX, multithread = TRUE!!!
err_R = learnErrors(filtRs, randomize=TRUE) #WHEN IN LINUX, multithread = TRUE!!!

plot_err_F <- plotErrors(err_F, nominalQ=TRUE)
ggsave(paste0(outdir_trimmed,"/err_rates_F.pdf"), plot_err_F, dpi = 300, height = 10, width = 10)

plot_err_R <-plotErrors(err_R, nominalQ=TRUE)
ggsave(paste0(outdir_trimmed,"/err_rates_R.pdf"), plot_err_R, dpi = 300, height = 10, width = 10)


#PENDIENTES EN LA SECCIÓN
# tamaño de plots
# mean quality of raw and trimmed reads en reporte



###-------STEP 3: DADA2 PROCESSING AND REPORT //  CHECKPOINT ---------

# dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# sample inference
dadaFs <- dada(derepFs, err=err_F) #WHEN IN LINUX, multithread = TRUE!!!
dadaRs <- dada(derepRs, err=err_R) #WHEN IN LINUX, multithread = TRUE!!!

# merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# construct sequence table
seqtab <- makeSequenceTable(mergers)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", verbose=TRUE)  #WHEN IN LINUX, multithread = TRUE!!!

# checkpoint
saveRDS(seqtab.nochim, file="seqtab_nochim.rds")


##   DADA2 REPORT

getN <- function(x) sum(getUniques(x))
track <- cbind(filtered_report, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("sample_id", "input", "filtered","percentage_reads_kept", "denoisedF", "denoisedR", "merged", "nonchim")

track

write.csv(track, "dada2_report.csv", row.names = F)



###-------STEP 4: ASSIGN TAXONOMY  // CHECKPOINT no chimera file ---------

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", verbose = TRUE) #WHEN IN LINUX, multithread = TRUE!!! 
taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa.gz", verbose = TRUE)

# checkpoint
saveRDS(taxa, file="taxa.rds")

###-------STEP 5: CREATE PHYLOSEQ OBJECT  // CHECKPOINT Equivalences of ASVs and seqs -----

library(phyloseq)
library(vegan)
library(reshape2)
library(colorspace)
library(adespatial)
library(DESeq2)
library(microbiomeMarker)
library(edgeR)
library(RColorBrewer)

# Import metadata
file.path = "RealData/" 
metadata <- read.csv(paste0(file.path,"metadata.txt"), sep = "\t")

#edit metadata
#metadata$sample_id <- gsub("_S.*_L.*_R.*", "", metadata$sample_id)
metadata <- read.table(paste0(file.path,"metadata.txt"), header = T)
metadata$sample_id <- sapply(strsplit(basename(metadata$sample_id), "-"), `[`, 1)
metadata <- metadata[!duplicated(metadata),]
samplenames <- metadata$sample_id
#sample_metadata<- as.data.frame(metadata[,-1])
sample_metadata<- as.data.frame(metadata)
rownames(sample_metadata) <- samplenames
#colnames(sample_metadata) <- "category"

sample_metadata<- sample_data(sample_metadata)

# Import dada2 results
seqtab.nochim <- readRDS("seqtab_nochim.rds")
taxa <- readRDS("taxa.rds")

# Create taxonomy table object
tax.table <- tax_table(taxa)

# Create OTU table object
otu.table <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
sample_names(otu.table) <- samplenames

# Assign new taxa names to rows of otu.table and tax.table
new_taxa_names <- paste0("ASV_", seq_len(nrow(tax.table)))

### CHECKPOINT  --- NO FUNCIONA AQUÍ. NO SE GUARDA CON SEQUENCES
# Save a csv with the equivalences between ASV names and sequences from DADA2
# equivalences <- as.data.frame(cbind(ASV_id = new_taxa_names, seq_id = colnames(otu.table)))
# write.csv(equivalences, "equivalences.csv", row.names = F)

rownames(tax.table) <- new_taxa_names
colnames(otu.table) <- new_taxa_names

# Combine into a phyloseq object
ps <- phyloseq(otu.table, tax.table, sample_metadata)



###-------STEP 6: NORMALIZATION and FILTERING with DeSeq2-----

# Convert phyloseq object to DESeq2 object
dds <- phyloseq_to_deseq2(ps, ~category)

# Pre-filter low-count taxa (e.g., keep taxa with at least 4 counts in at least 20% of the samples)
dds <- dds[rowSums(counts(dds) >= 5) >= round(length(dds$sample_id)*0.1, 0), ] # counts(dds) >= 4) >= round(length(dds$sample_id)*0.2, 0)

#estimate normalization factors
dds <- estimateSizeFactors(dds, type = 'poscounts')

# Perform DESeq normalization
dds <- DESeq(dds, test = "Wald", fitType = "local")

# Extract normalized counts
norm_counts <- counts(dds, normalized = TRUE)

# Update the OTU table with normalized counts
otu_table_norm <- otu_table(norm_counts, taxa_are_rows = TRUE)

# Create a new phyloseq object with the normalized OTU table
ps_norm_deseq <- phyloseq(otu_table_norm, tax_table(ps), sample_data(ps))

# Verify the normalization
ps_norm_deseq





###-------STEP 7: RAREFACTION CURVES-----

# Extract OTU table as matrix
otu_matrix <- as(otu_table(ps_norm_deseq), "matrix")

otu_matrix <- round(otu_matrix, 0)
otu_matrix<- t(otu_matrix)

# Inspect the count distribution and check for singletons
otu_counts <- colSums(otu_matrix)
summary(otu_counts)
sum(otu_counts == 1)

# Plot rarefaction curves
png(filename = "rarefaction_curves_deseq2.png",width = 12,height = 10, units = "in", res = 300)
rarecurve(otu_matrix, cex=0.6,step = 20,lwd = 2,ylab= "Number of ASVs",xlab = "Sequences",label = T,col = palette("r3"))
dev.off()



###-------STEP 8: ALPHA DIVERSITY------

#round normalized counts to integers
ps_norm_deseq <- transform_sample_counts(ps_norm_deseq,  function(x) round(x,0) )

# Get metrics
alpha_diversity <- estimate_richness(ps_norm_deseq, measures = c("Observed", "Shannon", "Simpson"))
alpha_diversity$Pielou <- alpha_diversity$Shannon / log(alpha_diversity$Observed) # Calculate Pielou's J manually

# Merge alpha diversity with metadata
alpha_diversity$sample_id <- rownames(alpha_diversity)
merged_data <- merge(alpha_diversity, metadata, by = "sample_id")

# Define function to perform pairwise Wilcoxon tests
perform_pairwise_comparisons <- function(data, metric, group_col) {
  
  # Perform pairwise Wilcoxon test
  pairwise_results <- pairwise.wilcox.test(data[[metric]], data[[group_col]], p.adjust.method = "BH")
  
  #melt data
  results_df <- melt(pairwise_results$p.value)

  return(results_df)
}


# Loop through each metric and perform pairwise Wilcoxon tests
results_list <- list()

metrics <- c("Observed", "Shannon", "Simpson", "Pielou")
for (metric in metrics) {
  results_list[[metric]] <- perform_pairwise_comparisons(merged_data, metric, "category")
}

# Combine all results into a single data frame
results_df <- do.call(rbind, results_list)
results_df$metric <- rownames(results_df)
rownames(results_df) <- NULL

write.csv(results_df, "alpha_results.csv", row.names = F)


# Exclude sample_id from merged_data
data_for_boxplot <- merged_data %>%
  select(-sample_id)

# Melt the data frame for easier plotting
data_for_boxplot_long <- tidyr::gather(data_for_boxplot, key = "Metric", value = "Value", -category)

#order metrics
data_for_boxplot_long$Metric <- factor(data_for_boxplot_long$Metric, levels = metrics)

# Define the Set3 color palette
set3_palette <- brewer.pal(12, "Set3")

# Plot using Set3 palette
plot_alpha <- ggplot(data_for_boxplot_long, aes(x = category, y = Value, fill = category)) +
  geom_boxplot(width = 0.5) +
  #geom_violin(width = 0.8, alpha = 0.5) +
  facet_wrap(~ Metric, scales = "free_y", ncol = length(unique(data_for_boxplot_long$Metric))) +
  labs(x = "Category", y = "Value", fill = "Category") +
  scale_fill_manual(values = set3_palette, name = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),  
        panel.grid.major = element_blank(), 
        strip.background = element_blank(),  
        strip.text = element_text(face = "bold", size = 13),
        axis.title.x = element_blank())

plot_alpha

ggsave("plot_alpha_deseq2.png", plot_alpha, dpi = 300, height = 5, width = 10, bg = "white")

#PENDIENTES
# estilo plot
# significativos plot
# export pairwise comparisons
# abundancias relativas??



###-------STEP 9: BETA DIVERSITY------

# Transform counts to relative abundances
ps_rel <- transform_sample_counts(ps_norm_deseq, function(x) x / sum(x))

# Convert phyloseq data to dataframe
otu_table_norm <- as.data.frame(otu_table(ps_rel))
otu_table_norm <- t(otu_table_norm)

# Calculate Bray-Curtis distance matrix
bray_dist_rel <- vegdist(otu_table_norm, method = "bray")

#permanova 
permanova_results_rel <- adonis2(bray_dist_rel ~ metadata$category, permutations = 9999)

# beta dispersion
betadisper_results_rel <- betadisper(bray_dist_rel, metadata$category)
betadisper_results_rel_final<- anova(betadisper_results_rel)

# Summary table for PERMANOVA and betadisper
permanova_summary_rel <- data.frame(
  Method = "PERMANOVA",
  stat = permanova_results_rel$R2[1],
  P.value = permanova_results_rel$`Pr(>F)`[1]
)

betadisper_summary_rel <- data.frame(
  Method = "Betadisper",
  stat = anova(betadisper_results_rel)$`F value`[1],
  P.value = anova(betadisper_results_rel)$`Pr(>F)`[1]
)

beta_diversity_summary_rel <- rbind(permanova_summary_rel, betadisper_summary_rel)
print(beta_diversity_summary_rel)

write.csv(beta_diversity_summary_rel, "beta_results.csv", row.names = F)


#NMDS
nmds_results_rel <- metaMDS(bray_dist_rel, k = 2, trymax = 1000)
nmds_plot_rel <- ggplot(as.data.frame(nmds_results_rel$points), aes(MDS1, MDS2, color = metadata$category)) +
  geom_point(size = 5, alpha = 0.75) +
  geom_text(aes(label = metadata$sample_id), vjust = -1.5, hjust = 0.5, size = 3) +
  labs(x = "NMDS1", y = "NMDS2", title = "NMDS Plot (Relative Abundances)") +
  theme_minimal()

stress_score <- nmds_results_rel$stress
stress_score

ggsave("nmds_plot_rel.png", nmds_plot_rel, dpi = 300, height = 6, width = 7, bg = "white")


# PCoA
pcoa_results_rel <- cmdscale(bray_dist_rel, eig = TRUE, k = 2)
pcoa_df_rel <- as.data.frame(pcoa_results_rel$points)
pcoa_df_rel$category <- metadata$category
pcoa_plot_rel <- ggplot(pcoa_df_rel, aes(x = V1, y = V2, color = category)) +
  geom_point(size = 5, alpha = 0.75) +
  geom_text(aes(label = metadata$sample_id), vjust = -1.5, hjust = 0.5, size = 3) +
  labs(x = "PCoA Axis 1", y = "PCoA Axis 2", title = "PCoA Plot (Relative Abundances)") +
  theme_minimal()

ggsave("pcoa_plot_rel.png", pcoa_plot_rel, dpi = 300, height = 6, width = 7, bg = "white")



# Calculate LCBD values using the beta.div function
beta_div_results <- beta.div(bray_dist_rel, method = "hellinger", nperm = 9999)
lcbd_values <- beta_div_results$LCBD

lcbd_values <- as.data.frame(lcbd_values)
lcbd_values$Sample <- rownames(lcbd_values)


# Melt the data to long format
ps_melt <- psmelt(ps_rel)

# GENUS STACKED BARPLOTS
ps_agg <- ps_melt %>%
  group_by(Sample, category, Genus) %>%
  summarize(Abundance = sum(Abundance)) %>%
  ungroup()

ps_agg$Sample <- factor(ps_agg$Sample, levels = unique(ps_agg$Sample))

#add category to lcbds
lcbd_values <- dplyr::left_join(lcbd_values, ps_agg %>% dplyr::select(Sample, category) %>% distinct(), by = "Sample")



# Plot stacked barplots TOP 50 ONLY
mean_abundances <- ps_agg %>%
  group_by(Genus) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance))

mean_abundances_top <- mean_abundances[1:50,]

ps_agg_top <- ps_agg[ps_agg$Genus %in% mean_abundances_top$Genus,]


set.seed(420)
n_taxa <- length(unique(ps_agg_top$Genus))
color_palette <- brewer.pal(min(n_taxa, 12), "Set3")
if (n_taxa > 12) {
  color_palette_2 <- rainbow_hcl(n_taxa-12)
  color_palette_2 <- sample(color_palette_2)
}

color_palette <- c(color_palette, color_palette_2)
color_palette <- sample(color_palette, n_taxa)


top_lcbd <- ggplot(ps_agg_top, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_point(data = lcbd_values, aes(x = Sample, y = -0.05, size = lcbd_values), shape = 21, fill = "black") +
  facet_wrap(~ category, scales = "free_x", ncol = 1) +
  labs(x = "Sample", y = "Relative Abundance", fill = "Genus") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "black", fill = NA, size = 1),  
        panel.grid.major = element_blank(), 
        strip.background = element_blank(),  
        strip.text = element_text(face = "bold", size = 13),
        axis.title.x = element_blank())  +
  scale_fill_manual(values = color_palette, guide = guide_legend(ncol = 2))

top_lcbd

ggsave("top_lcbd_deseq2_top50.png", top_lcbd, dpi = 300, height = 10, width = 12, bg = "white")

#PENDIENTES DE LA SECCIÓN
# añadir permanova, betadisper y stress score a ordination plots
# ordenar taxa en stacked barplot
# Y SACAR TOP 50 generos mas abundantes antes de plotear
# lcbd desordena el stacked bar. fix.


###-------STEP 10: DIFFERENTIAL ABUNDANCE-------

# extract significant results
alpha <- 0.01
res <- results(dds, alpha = alpha, contrast = c("category", "Midgut", "FTA")) # pairwise has to be done here. LOOP.
res <- res[order(res$padj), ]
res <-  res[!is.na(res$padj),]

res_sig <- res[res$padj < alpha, ]
res_sig <- cbind(as(res_sig, "data.frame"), as(tax_table(ps)[rownames(res_sig), ], "matrix"))
res_sig

res <- cbind(as(res, "data.frame"), as(tax_table(ps)[rownames(res), ], "matrix"))


# choose colors
set.seed(424) 
unique_genera <- unique(res$Genus)
num_to_select <- length(unique_genera) # Or specify a number
selected_genera <- sample(unique_genera, num_to_select)
color_palette <- sample(colors(), num_to_select)
color_mapping <- setNames(color_palette, selected_genera)

# Add a column in the dataframe for colors
res$Color <- ifelse(res$Genus %in% selected_genera, color_mapping[res$Genus], "gray")
res_sig <- res[res$padj < alpha, ]

diffabun_plot <- ggplot(res, aes(x = baseMean, y = log2FoldChange, color = Color)) + # use x = Genus and color = Phylum of res_sig for clean plot
  geom_point(size = 3, alpha = 0.6) +
  scale_color_identity()+
  geom_text(data = res_sig, aes(x = baseMean, y = log2FoldChange, label = Genus), vjust = 0, hjust = 0, size = 3) + # Text labels for significant points
  theme_minimal() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))+
  guides(color = guide_legend(ncol = 3))

diffabun_plot

ggsave("diffabun_deseq2_FTA_is_base.png", diffabun_plot, dpi = 300, height = 6, width = 6, bg = "white")

