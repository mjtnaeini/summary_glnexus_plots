#!/bin/bash

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
        chr20 chr21 chrX chrY)

for CHR in "${CHROMS[@]}"; do
  qsub -v CHR="$CHR" <<EOF
#!/bin/bash
#PBS -P ox63
#PBS -q normal
#PBS -N analyse_variant_$CHR
#PBS -l walltime=48:00:00
#PBS -l mem=128GB
#PBS -l jobfs=256GB
#PBS -l ncpus=1
#PBS -l wd
#PBS -l storage=gdata/ox63+scratch/ox63
#PBS -M m.naeini@garvan.org.au

module load R/4.4.2

Rscript /scratch/ox63/mn4616/analyse_variants_chr.R \
  --chr "\$CHR" \
  --mode "All" \
  --vcf_rdata /scratch/ox63/mn4616/snv_indel.RData \
  --output_dir /scratch/ox63/mn4616/plots_indel
EOF
done

