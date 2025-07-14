# Joint SNV/INDEL Variants Summary Plots

This repository contains a PBS- and R-based pipeline to perform per-chromosome variant analysis, burden calculation, and visualization for SNVs and indels. It is designed for HPC systems like **Gadi (NCI Australia)** and includes chromosome-wise analysis, data merging, and population-stratified plotting.

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

```
summary_plots_glnexus/
├── analyse_variants_chr.pbs.sh         # Submit one job per chromosome
├── analyse_variants_chr.R              # R script for SNV/indel processing
├── combine_rdata.pbs                   # PBS job to merge chromosome outputs
├── combine_rdata.R                     # R script to merge all RData files
├── generate_variant_plot_inputs.pbs    # Final summary plotting PBS job
├── generate_variant_plot_inputs.R      # R script for final population plots
└── README.md                           # This file
```

---

## ⚙️ Workflow

### 1. Run Per-Chromosome Analysis

Submit per-chromosome jobs:

```bash
bash analyse_variants_chr.pbs.sh
```

Each job will:
- Load `.RData` input (joint-called VCF)
- Filter by variant type (SNV / Indel / All)
- Classify variants as `Private`, `Shared`, or `Common`
- Save plots and intermediate data per chromosome

### 2. Combine All Chromosome Outputs

After all jobs complete, combine per-chromosome outputs:

```bash
qsub combine_rdata.pbs
```

This creates a merged file:
```
All_autosomes_combined_VCF_df.RData
```

### 3. Generate Summary Plots Across Genome

Generate genome-wide burden plots:

```bash
qsub generate_variant_plot_inputs.pbs
```

This creates publication-ready plots grouped by:
- Sample
- Population group
- Variant type (SNV/Indel)

---

## 📊 Output Files

Output files per mode (SNV / Indel / All) include:

- `*_sample_burden.pdf`: SNV/Indel burden per sample  
- `*_population_burden.pdf`: Burden stratified by population  
- `*_pcoa_group1.pdf`: PCoA (colored by population group)  
- `*_pcoa_continental.pdf`: PCoA (colored by continental group)  
- `*_VCF_df.RData`: Chromosome- or genome-level RData output  

---

## 🔧 Requirements

### Gadi Modules

```bash
module load R/4.4.2
```

### Required R Packages

- `argparser`  
- `dplyr`  
- `tidyr`  
- `ggplot2`  
- `data.table`  
- `vegan`  

---

## 📝 Notes

- The pipeline assumes preprocessed joint VCF converted to `.RData` format with `vcf_df` loaded.  
- The metadata file used for sample population info is:

  ```
  tsv file
  ```

- Each script accepts `--mode` argument with options:
  - `SNV`  
  - `Indel`  
  - `All`  

- Final outputs are saved under the specified `--output_dir`.

