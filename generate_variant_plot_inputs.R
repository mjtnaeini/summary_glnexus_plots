#!/bin/bash
#PBS -P ox63
#PBS -q hugemem
#PBS -N generate_variant_plot
#PBS -l walltime=48:00:00
#PBS -l mem=1400GB
#PBS -l ncpus=1
#PBS -l wd
#PBS -l storage=gdata/ox63+scratch/ox63
#PBS -M m.naeini@garvan.org.au
#PBS -m abe

module load R/4.4.2

Rscript /scratch/ox63/mn4616/generate_variant_plot_inputs.R \
  --mode "SNV" \
  --vcf_rdata /scratch/ox63/mn4616/plots_snv_indel/All_autosomes_combined_VCF_df.RData \
  --output_dir /scratch/ox63/mn4616/plots_snv_indel/

(base) [mn4616@gadi-login-07 summary_plots_glnexus]$ cat generate_variant_plot_inputs.R
#!/usr/bin/env Rscript

# =====================================================
# Argument parsing
# =====================================================
library(argparser)
p <- arg_parser("Generate summary tables and plots from pre-processed VCF data")
p <- add_argument(p, "--mode", help = "SNV, Indel, or All", default = "All")
p <- add_argument(p, "--vcf_rdata", help = "Pre-processed VCF .RData path", default = "All_autosomes_combined_VCF_df.RData")
p <- add_argument(p, "--output_dir", help = "Directory to save outputs", default = "./output")
args <- parse_args(p)

dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# =====================================================
# Load libraries
# =====================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(data.table)

# =====================================================
# Load preprocessed data
# =====================================================
load(args$vcf_rdata)  # expects combined_vcf_df already processed
vcf_df <- combined_vcf_df
sample_info <- read.delim("/g/data/ox63/marjan/projects/indo_genome/sample_info/indo_info_sheet_with_coverage.tsv")

# =====================================================
# Subset by mode (SNV, Indel, All)
# =====================================================
if (args$mode == "SNV") {
  vcf_df <- vcf_df[nchar(vcf_df$REF) == 1 & nchar(vcf_df$ALT) == 1, ]
} else if (args$mode == "Indel") {
  vcf_df <- vcf_df[!(nchar(vcf_df$REF) == 1 & nchar(vcf_df$ALT) == 1), ]
}

variant_label <- args$mode
output_prefix <- file.path(args$output_dir, variant_label)
sample_cols <- intersect(sample_info$ID, colnames(vcf_df))

# =====================================================
# Plot per sample
# =====================================================
plot_df <- data.frame(ID = vcf_df$ID, vcf_df[, sample_cols], Shared_Status = vcf_df$Shared_Status_Sample)
plot_long <- reshape(plot_df,
                     varying = sample_cols,
                     v.names = "Genotype",
                     times = sample_cols,
                     timevar = "Sample",
                     idvar = c("ID", "Shared_Status"),
                     direction = "long")

plot_long$Present <- !sapply(strsplit(as.character(plot_long$Genotype), ":"), `[`, 1) %in% c("./.", "0/0")
plot_long$Population <- sample_info$Population.Groups.1[match(plot_long$Sample, sample_info$ID)]
plot_long <- subset(plot_long, Present)

agg_df <- plot_long %>%
  group_by(Sample, Shared_Status, Population) %>%
  summarise(SNV_Indel_Count = n(), .groups = "drop")

agg_df$Sample <- factor(agg_df$Sample, levels = agg_df %>%
                          summarise(total = sum(SNV_Indel_Count), .by = Sample) %>%
                          arrange(desc(total)) %>%
                          pull(Sample))

agg_df$Population <- factor(agg_df$Population, levels = unique(agg_df[order(agg_df$SNV_Indel_Count, decreasing = TRUE), "Population"])[[1]])
agg_df$Shared_Status <- factor(agg_df$Shared_Status, levels = c("Common", "Shared", "Private"))

plot1 <- ggplot(agg_df, aes(x = Sample, y = SNV_Indel_Count, fill = Shared_Status)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.9) +
  geom_text(aes(label = SNV_Indel_Count), position = position_stack(vjust = 0.5), size = 2.5, angle = 90) +
  scale_fill_manual(values = c("Common" = "lightgreen", "Shared" = "steelblue", "Private" = "tomato")) +
  facet_grid(~ Population, scales = "free_x", space = "free") +
  theme_bw() +
  labs(title = "SNV/Indel Burden per Sample", x = "Sample", y = "Variant Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(output_prefix, "_sample_burden.pdf"), plot1, width = 12, height = 5)

# =====================================================
# Burden Summary
# =====================================================
dt <- plot_long %>%
  count(Sample, name = "SV_Count") %>%
  arrange(desc(SV_Count)) %>%
  mutate(Sample = factor(Sample, levels = Sample))

plot2 <- ggplot(dt, aes(x = Sample, y = SV_Count)) +
  geom_bar(stat = "identity", fill = "lightpink", alpha = 0.9) +
  geom_text(aes(label = SV_Count), position = position_stack(vjust = 0.5), size = 2.5, angle = 90) +
  theme_bw() +
  labs(title = "Burden per Sample", x = "Sample", y = "SV Count (n)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(paste0(output_prefix, "_sample_burden.pdf"), plot2, width = 9, height = 4)

# =====================================================
# Population group-based burden (Group 1)
# =====================================================
extract_gt <- function(x) sapply(strsplit(x, ":"), `[`, 1)
genotype_matrix <- sapply(vcf_df[, sample_cols], extract_gt)
if (is.vector(genotype_matrix)) {
  genotype_matrix <- matrix(genotype_matrix, ncol = 1)
  colnames(genotype_matrix) <- sample_cols
}

group_names <- unique(sample_info$Population.Groups.1)

presence_matrix_group1 <- sapply(group_names, function(group) {
      group_samples <- intersect(sample_info$ID[sample_info$Population.Groups.1 == group], sample_cols)
      apply(genotype_matrix[, group_samples, drop = FALSE], 1, function(gts) any(!(gts %in% c("./.", "0/0"))))
})
rownames(presence_matrix_group1) <- vcf_df$ID

group1_support <- rowSums(presence_matrix_group1)
vcf_df$Group1_Support_Count <- group1_support
vcf_df$Shared_Status_group1 <- ifelse(group1_support == 1, "Population-specific",
                                       ifelse(group1_support == length(group_names), "Pan-population", "Multi-population"))

plot_df <- data.frame(ID = rownames(presence_matrix_group1),
                      presence_matrix_group1,
                      Shared_Status = vcf_df$Shared_Status_group1)

plot_long <- reshape(plot_df,
                     varying = setdiff(names(plot_df), c("ID", "Shared_Status")),
                     v.names = "Present",
                     times = setdiff(names(plot_df), c("ID", "Shared_Status")),
                     timevar = "Population",
                     idvar = c("ID", "Shared_Status"),
                     direction = "long")

plot_long <- subset(plot_long, Present)

agg_df <- aggregate(Present ~ Population + Shared_Status, data = plot_long, FUN = length)
colnames(agg_df)[3] <- "SNV_Indel_Count"

agg_df$Population <- factor(agg_df$Population, levels = agg_df %>%
                                  group_by(Population) %>%
                                  summarise(total = sum(SNV_Indel_Count)) %>%
                                  arrange(desc(total)) %>%
                                  pull(Population))

agg_df$Shared_Status <- factor(agg_df$Shared_Status, levels = c("Pan-population", "Multi-population", "Population-specific"))

plot3 <- ggplot(agg_df, aes(x = Population, y = SNV_Indel_Count, fill = Shared_Status)) +
      geom_bar(stat = "identity", position = "stack", alpha = 0.9) +
      geom_text(aes(label = SNV_Indel_Count), position = position_stack(vjust = 0.5), size = 3, angle = 90) +
      scale_fill_brewer(palette = "Accent") +
      theme_bw() +
      labs(title = "SNV/Indel Burden per Population", x = "Population", y = "Variant Count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(output_prefix, "_population_burden.pdf"), plot3, width = 9, height = 6)


# =====================================================
# PCoA plot - Group 1
# =====================================================

snv_matrix <- apply(genotype_matrix, 2, function(gt) as.numeric(!gt %in% c("./.", "0/0")))
rownames(snv_matrix) <- rownames(genotype_matrix)
snv_binary_t <- t(snv_matrix)

snv_dist <- vegan::vegdist(snv_binary_t, method = "jaccard", binary = TRUE)
pcoa_res <- cmdscale(snv_dist, eig = TRUE, k = 2)

# PCoA grouped by Population.Groups.1
pcoa_df <- data.frame(SampleID = rownames(pcoa_res$points),
                      PC1 = pcoa_res$points[, 1],
                      PC2 = pcoa_res$points[, 2])
pcoa_df <- merge(pcoa_df, sample_info[, c("ID", "Population.Groups.1")], by.x = "SampleID", by.y = "ID")

plot4 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Population.Groups.1)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_brewer(palette = "Set3") +
  theme_bw() +
  labs(title = "PCoA of SNVs/Indels (Group 1)", x = "PCoA 1", y = "PCoA 2") +
  theme(legend.title = element_blank())

ggsave(paste0(output_prefix, "_pcoa_group1.pdf"), plot4, width = 7, height = 6)

# Labeled plot
plot5 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Population.Groups.1)) +
  geom_text(aes(label = SampleID), size = 3) +
  scale_color_brewer(palette = "Set3") +
  theme_bw() +
  labs(title = "PCoA of SNVs/Indels (labeled, Group 1)", x = "PCoA 1", y = "PCoA 2") +
  theme(legend.title = element_blank())

ggsave(paste0(output_prefix, "_pcoa_group1_labelled.pdf"), plot5, width = 8, height = 6)

# =====================================================
# PCoA plot - Continental Groups
# =====================================================
pcoa_df <- data.frame(SampleID = rownames(pcoa_res$points),
                      PC1 = pcoa_res$points[, 1],
                      PC2 = pcoa_res$points[, 2])
pcoa_df <- merge(pcoa_df, sample_info[, c("ID", "Continental.Groups")], by.x = "SampleID", by.y = "ID")

plot6 <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Continental.Groups)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  labs(title = "PCoA of SNVs/Indels (by Continental Group)", x = "PCoA 1", y = "PCoA 2") +
  theme(legend.title = element_blank())

ggsave(paste0(output_prefix, "_pcoa_continental.pdf"), plot6, width = 7, height = 6)

message(paste0("Finished summary and plots for mode: ", args$mode))
