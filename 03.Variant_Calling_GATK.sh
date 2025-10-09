#!/bin/bash
# ============================================================
# Script: 03.Variant_Calling_GATK.sh
# Author: Doris Repusic
# Description:
#   PBS job submission script for variant calling and base quality 
#   score recalibration (BQSR) using GATK. 
#
#   Workflow per sample:
#     1. Initial variant calling (HaplotypeCaller)
#     2. Extract SNPs and apply filtering
#     3. Generate BQSR recalibration table
#     4. Apply BQSR to BAM file
#     5. Second round of variant calling on recalibrated BAM
#     6. SNP selection and filtering
#     7. Merge VCFs across all samples
#     8. Filter merged VCF to retain only complete chromosomes
#
# Usage:
#   qsub -J 0-N 03.Variant_Calling_GATK.sh
#   (where N = number of BAM samples - 1)
#
# Requirements:
#   - GATK 4.x
#   - bcftools
#   - Reference genome fasta (indexed)
# ============================================================

# ---------------- PBS Directives ----------------
# Queue selection
#PBS -q q2
# Job name
#PBS -N gatk_varcalling
# Resources
#PBS -l select=ncpus=16:mem=400g
# Job array (adjust to number of input BAMs)
#PBS -J 0-200
# Merge stdout and stderr
#PBS -j oe

# Move to working directory
cd $PBS_O_WORKDIR

# ---------------- Configurable Parameters ----------------
REF=/hg38.fa              # Reference genome fasta
INPUT_DIR=.               # Input directory containing BAMs

# ---------------- File Selection ----------------
# Select input BAMs with marked duplicates and read groups
IN_SEQ=($INPUT_DIR/*_sorted_dedupe_groups.bam)
FILE=${IN_SEQ[$PBS_ARRAY_INDEX]}
BASE=${FILE##*/}          # Remove path
BASE=${BASE%%_*}          # Remove suffix for sample name

echo $BASE

# ---------------- Variant Calling: Round 1 ----------------
gatk HaplotypeCaller \
    -R $REF \
    -I $FILE \
    -O ${BASE}_variants.vcf

# Extract SNPs
gatk SelectVariants \
    -R $REF \
    -V ${BASE}_variants.vcf \
    -select-type SNP \
    -O ${BASE}_raw_snps.vcf

# Apply hard filters
gatk VariantFiltration \
    -V ${BASE}_raw_snps.vcf \
    -filter "QD < 2.0"        --filter-name "QD2" \
    -filter "QUAL < 30.0"     --filter-name "QUAL30" \
    -filter "SOR > 3.0"       --filter-name "SOR3" \
    -filter "FS > 60.0"       --filter-name "FS60" \
    -filter "MQ < 40.0"       --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${BASE}_filtered_snps.vcf

# Keep only unfiltered SNPs
gatk SelectVariants \
    --exclude-filtered \
    -V ${BASE}_filtered_snps.vcf \
    -O ${BASE}_bqsr_snps.vcf

# ---------------- Base Quality Score Recalibration ----------------
gatk BaseRecalibrator \
    -I ${BASE}_sorted_dedupe_groups.bam \
    -R $REF \
    --known-sites ${BASE}_bqsr_snps.vcf \
    -O ${BASE}_recal_data.table

gatk ApplyBQSR \
    -R $REF \
    -I ${BASE}_sorted_dedupe_groups.bam \
    --bqsr-recal-file ${BASE}_recal_data.table \
    -O ${BASE}_sorted_dedupe_recal.bam

# ---------------- Variant Calling: Round 2 ----------------
gatk HaplotypeCaller \
    -R $REF \
    -I ${BASE}_sorted_dedupe_recal.bam \
    -O ${BASE}_recal_variants.vcf

gatk SelectVariants \
    -R $REF \
    -V ${BASE}_recal_variants.vcf \
    -select-type SNP \
    -O ${BASE}_recal_snps.vcf

gatk VariantFiltration \
    -V ${BASE}_recal_snps.vcf \
    -filter "QD < 2.0"        --filter-name "QD2" \
    -filter "QUAL < 30.0"     --filter-name "QUAL30" \
    -filter "SOR > 3.0"       --filter-name "SOR3" \
    -filter "FS > 60.0"       --filter-name "FS60" \
    -filter "MQ < 40.0"       --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${BASE}_recal_filtered_snps.vcf

# ---------------- Joint Genotyping / Merging ----------------
# Merge all filtered SNPs across samples into one VCF
bcftools merge /*_recal_filtered_snps.vcf.gz \
    -m snps -0 -O z -o /snps_recal_all.vcf.gz

# Index merged VCF
bcftools index /snps_recal_all.vcf.gz

# ---------------- Chromosome Filtering ----------------
# Keep only canonical chromosomes (1–22, X, Y, M)
bcftools view /snps_recal_all.vcf.gz | \
    grep -w '^#\|^#CHROM\|chr[1-9]\|chr[1-2][0-9]\|chr[X]\|chr[Y]\|chr[M]' \
    > /snps_recal_all_xym.vcf.gz

# ---------------- Notes ----------------
# Per-sample outputs:
#   <sample>_variants.vcf
#   <sample>_filtered_snps.vcf
#   <sample>_sorted_dedupe_recal.bam
#   <sample>_recal_filtered_snps.vcf
#
# Final outputs:
#   snps_recal_all.vcf.gz       (merged across samples)
#   snps_recal_all_xym.vcf.gz   (filtered to canonical chromosomes)
# ============================================================
