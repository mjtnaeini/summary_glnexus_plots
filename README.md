# Joint SNV/INDEL Variants Analysis

This repository contains a PBS- and R-based pipeline to perform per-chromosome variant analysis, burden calculation, and visualization for SNVs and indels. It includes chromosome-wise analysis, data merging, and population-stratified plotting.

---

## 🧬 Overview

This pipeline processes jointly-called VCF data (e.g., from GLnexus), extracts SNVs and indels, and generates summary statistics and plots:

- Chromosome-wise variant processing
- Population-based shared/private/common classification
- Sample and population burden plots
- Principal Coordinates Analysis (PCoA)
- Final merged results across autosomes

---

## 📁 Project Structure

summary_plots_glnexus/
├── analyse_variants_chr.pbs.sh # Submit one job per chromosome
├── analyse_variants_chr.R # R script for SNV/indel processing
├── combine_rdata.pbs # PBS job to merge chromosome outputs
├── combine_rdata.R # R script to merge all RData files
├── generate_variant_plot_inputs.pbs # Final summary plotting PBS job
├── generate_variant_plot_inputs.R # R script for final population plots
└── README.md # This file

---

## ⚙️ Workflow

### 1. Run Per-Chromosome Analysis

Submit per-chromosome jobs using:

```bash
bash analyse_variants_chr.pbs.sh
Each job:

Loads .RData input (joint-called VCF data)

Filters by variant type (SNV/Indel/All)

Assigns Private, Shared, or Common status

Saves plots and summary for each chromosome

2. Combine RData Files
After all chromosomes are processed, merge them:

qsub combine_rdata.pbs
This step creates:
All_autosomes_combined_VCF_df.RData
3. Generate Final Summary Plots
Create summary plots from the combined dataset:

qsub generate_variant_plot_inputs.pbs
This includes:

Population burden barplots

Optional per-sample barplots

Optional PCoA (Jaccard distance) plots

📊 Output Files
For each mode (SNV, Indel, or All), the following plots may be generated:

*_sample_burden.pdf: Sample-wise burden by shared status

*_population_burden.pdf: Population-wise burden by shared status

*_pcoa_group1.pdf: PCoA colored by population group

*_pcoa_continental.pdf: PCoA colored by continental group

*_VCF_df.RData: Processed data per chromosome or full genome

🔧 Requirements
Load modules on Gadi:

module load R/4.4.2
Required R Packages:
argparser

dplyr

tidyr

ggplot2

data.table

vegan

📝 Notes
The main input .RData file should contain a vcf_df object.

Sample metadata file must be:

tsv file
Use --mode to specify "SNV", "Indel", or "All" in each script.

Job scripts are PBS-compatible for execution on HPC environments.

