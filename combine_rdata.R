#!/usr/bin/env Rscript

# =====================================================
# Load libraries
# =====================================================
library(argparser)
library(dplyr)

# =====================================================
# Parse arguments
# =====================================================
p <- arg_parser("Combine autosomal vcf_chr RData files")
p <- add_argument(p, "--input_dir", help = "Directory containing RData files")
p <- add_argument(p, "--output_file", help = "Output file to save combined RData")

argv <- parse_args(p)

# =====================================================
# Define chromosomes to include (chr1 to chr22)
# =====================================================
chromosomes <- paste0("chr", 1:22)

vcf_list <- list()

# =====================================================
# Load and combine RData files containing `vcf_chr`
# =====================================================
for (chr in chromosomes) {
  file_path <- file.path(argv$input_dir, paste0("All_", chr, "_VCF_df.RData"))
  if (file.exists(file_path)) {
    # Load into environment and retrieve vcf_chr
    e <- new.env()
    load(file_path, envir = e)
    vcf_list[[chr]] <- e$vcf_chr
    rm(e)
  } else {
    warning(paste("Missing file for", chr))
  }
}

# =====================================================
# Combine and save
# =====================================================
combined_vcf_df <- do.call(rbind, vcf_list)
save(combined_vcf_df, file = argv$output_file)
